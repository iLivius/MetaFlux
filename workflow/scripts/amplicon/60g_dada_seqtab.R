#!/usr/bin/env Rscript
# DADA2's core denoising chain: turns filtered reads into a table of exact
# amplicon sequence variants (ASVs) — DADA2's alternative to OTU clustering,
# where every ASV is one exact sequence rather than a cluster of similar ones.
#
# Input: every sample's filtered FASTQ pair from 04.filtered/ (one pair per
# sample, produced by dada_filter.R) plus that sample's filter_stats JSON
# (reads in/out from the filtering step), collected here across the whole
# sample set for the final read-tracking table; the two error models that
# learn_errors (60d_learn_errors.R) learned from the same filtered reads,
# 05.dada2/error_model_R1.rds and _R2.rds; and the verdict of check_error_model,
# stats/dada2/error_model_check.txt, which only exists when the check let the run
# go on (both models converged, or the user accepted a model by its code). The
# verdict is copied into the log, so the log says which models the tables were
# made with.
#
# The error model says how often, say, a true A is read as a G at each quality
# score, separately for R1 and R2. Up to v2.6.1 this script learned it itself
# (learnErrors) before denoising; learning is now a step of its own, so that a
# model that did not converge can be checked, and accepted, before any denoising.
# The models are DADA2's own objects, saved unchanged, so dada() below uses them
# exactly as it used the objects in memory before.
#
# What it does, in order:
#   1. readRDS      — reads the two saved error models.
#   2. derepFastq   — dereplication: collapses each sample's reads down to
#      its unique sequences, each carrying a count (abundance). DADA2 works
#      on these unique sequences rather than every individual read, which is
#      what keeps it fast on large samples.
#   3. dada         — the actual denoising step: using the error models of
#      step 1, decides for each unique sequence whether it looks like a real,
#      distinct biological sequence (an ASV) or just an error-derived variant
#      of a more abundant one, and folds the latter into their true parent.
#   4. mergePairs   — stitches each sample's denoised R1 with its denoised R2
#      where they overlap, keeping only pairs whose overlap is long enough
#      and close enough a match to trust (see minOverlap/maxMismatch below).
#   5. makeSequenceTable — assembles every sample's merged ASVs into one
#      table: rows = samples, columns = ASV sequences, values = read counts.
#   6. removeBimeraDenovo — drops chimeras: PCR artefacts where an
#      incompletely-extended copy of one template re-primes on a different
#      template mid-reaction, so the resulting "sequence" reads as two
#      organisms stitched together. DADA2 specifically detects "bimeras" —
#      chimeras built from exactly two parents — purely from the data itself
#      ("de novo"), by testing whether a candidate ASV could be explained as
#      one more-abundant ASV's first part plus another's second part.
#
# Output: the chimera-free ASV table in two flavours (raw sequences as column
# names, and short ASV_N ids), a FASTA of the ASV sequences, and a per-sample
# read-tracking table (stripped -> filtered -> denoised -> merged ->
# non_chimeric). seqs.fasta and the ASV_N seqtab feed dada_length_filter (60h)
# next, which applies the post-hoc ASV length window; assign_taxonomy (80a)
# consumes the length-filtered versions of both further downstream. The plots
# of the error models are drawn by error_model_report (60e), before the check.

# ── Load packages and redirect output to the log file ───────────
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, jsonlite))

# sink() below redirects everything this script prints from here on —
# messages, warnings, DADA2's own progress output — into the rule's log file
# instead of the console, so the log: file Snakemake tracks has the full record.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

#' The R2 file of each R1 file, found by name.
#'
#' dada_filter writes both files of a sample under one name: {sample}_R1_filt.fastq.gz
#' and {sample}_R2_filt.fastq.gz, in the same folder. So the R2 partner of an R1 path is
#' that path with "_R1_filt.fastq.gz" at its end replaced by "_R2_filt.fastq.gz", looked
#' up among the R2 paths Snakemake handed over, whatever their order.
#'
#' @param r1_files The R1 paths, in the order the tables will have.
#' @param r2_inputs The R2 paths, in any order.
#' @return The R2 paths, one per R1 file and in the same order; NA where the partner is
#'   not among r2_inputs.
r2_partners_of <- function(r1_files, r2_inputs) {
  partners <- sub("_R1_filt\\.fastq\\.gz$", "_R2_filt.fastq.gz", r1_files)
  partners[!(partners %in% r2_inputs)] <- NA
  partners
}

# ── Read input file paths from Snakemake ─────────────────────────
# The R1 files are sorted by full path with R's default sort(), as up to v2.6.1. Their
# order becomes the row order of the two sequence tables and of the read tracking.
# sort() follows the alphabetical order of the shell's language setting (locale), so
# sample names that mix capitals and small letters, or "-", "_" and ".", can come out
# in one row order under LANG/LC_ALL=C and in another under a UTF-8 setting. That
# changes the order of the rows and, among ASVs with the same total read count, which
# one gets the lower ASV_N number (makeSequenceTable keeps such ties in the order the
# sequences first appear, which follows the rows). The counts stay the same: with the
# default pool: false every sample is denoised on its own, and with pool: true or
# "pseudo" the 16S test reads, renamed to sort differently, also gave the same counts
# in both settings. The sort is kept as it was, so a folder's tables keep the row order
# and the ASV numbers they had up to v2.6.1. (learn_errors sorts byte by byte instead,
# because there the order decides which samples the error models are learned from.)
filt_r1            <- sort(unlist(snakemake@input[["r1"]]))
filter_stats_files <- unlist(snakemake@input[["filter_stats"]])

# Each sample's name, from its filtered-R1 file name ("{sample}_R1_filt.fastq.gz").
# dada() and mergePairs() use these names to label every result by sample.
sample_names <- sub("_R1_filt\\.fastq\\.gz$", "", basename(filt_r1))

# The R2 file of each sample, found by its name (r2_partners_of above), in the order of
# the R1 files. Up to v2.6.1 the R2 list was sorted on its own and paired with the R1
# list by position, which assumed that both lists always sort alike. They do not: for
# samples named "Soil" and "Soil_R2", LANG/LC_ALL=C sorts the R1 files "Soil, Soil_R2"
# ("Soil_R1_filt" before "Soil_R2_R1_filt") but the R2 files "Soil_R2, Soil"
# ("Soil_R2_R2_filt" before "Soil_R2_filt", because a capital R sorts before a small f),
# so mergePairs would have been handed each sample's R1 reads with the other sample's
# R2 reads. When the two samples have different numbers of filtered reads (the usual
# case), mergePairs then stops the job with "Non-corresponding derep-class and
# dada-class objects." (it compares the number of reads in R1 and R2); when the numbers
# happen to be equal, it merges the wrong reads without any error. Pairing by name
# cannot mix up two samples, in any language setting.
filt_r2 <- r2_partners_of(filt_r1, unlist(snakemake@input[["r2"]]))

# A sample whose R2 file is not among the inputs stops the job with the reason. This
# should not happen (the rule hands over both files of every sample, as dada_filter
# names them); it would mean that the two input lists no longer follow one naming. The
# text goes into the log, then, after the log is closed, onto the terminal, as in
# learn_errors; the short stop() line makes the job fail.
if (anyNA(filt_r2)) {
  missing_samples <- sample_names[is.na(filt_r2)]
  shown_samples <- paste(head(missing_samples, 10), collapse = ", ")
  if (length(missing_samples) > 10) {
    shown_samples <- paste0(shown_samples, " and ", length(missing_samples) - 10, " more")
  }
  text <- sprintf(paste0(
    "[dada_seqtab] STOP: the filtered R2 file of sample(s) %s is not among the inputs of ",
    "this step. dada_seqtab pairs each <sample>_R1_filt.fastq.gz with the ",
    "<sample>_R2_filt.fastq.gz of the same name, so that the R1 reads of one sample are ",
    "never merged with the R2 reads of another, and it cannot go on without the partner. ",
    "This should not happen; please report it with %s."),
    shown_samples, snakemake@log[[1]])
  message(text)                                       # full text into the log
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  message(text)                                       # full text onto the terminal
  stop("dada_seqtab stopped: reason above and in its log", call. = FALSE)
}

names(filt_r1) <- sample_names
names(filt_r2) <- sample_names

message("[dada_seqtab] ", length(sample_names), " sample(s): ", paste(sample_names, collapse = ", "))

# ── Read parameters from Snakemake ───────────────────────────────
# No seed here: nothing in this script draws random numbers (dada() is
# deterministic). The one random step of DADA2, the draw of the samples
# learnErrors learns from, is in learn_errors, which sets amplicon.seed itself.
pool          <- snakemake@params[["pool"]]
omega_a       <- as.numeric(snakemake@params[["omega_a"]])
min_overlap   <- as.integer(snakemake@params[["merge_min_overlap"]])
max_mismatch  <- as.integer(snakemake@params[["merge_max_mismatch"]])
just_concat   <- as.logical(snakemake@params[["merge_just_concatenate"]])
trim_overhang <- as.logical(snakemake@params[["merge_trim_overhang"]])
chim_method   <- snakemake@params[["chimera_method"]]
min_fold      <- as.numeric(snakemake@params[["chimera_min_fold_parent"]])
allow_one_off <- as.logical(snakemake@params[["chimera_allow_one_off"]])
threads       <- as.integer(snakemake@threads)

# YAML false/true → Python False/True → R FALSE/TRUE (logical); "pseudo" → character.
# The is.character guard converts string "false"/"true" if YAML used quotes.
if (is.character(pool)) {
  pool <- switch(tolower(trimws(pool)), "false" = FALSE, "true" = TRUE, "pseudo" = "pseudo", FALSE)
}

# ── Create output directories ────────────────────────────────────
dir.create(dirname(snakemake@output[["seqtab_seqs"]]),  recursive = TRUE, showWarnings = FALSE)

# ── Read the error models learned by learn_errors ────────────────
# R1 and R2 have separate error models — on Illumina, reverse reads are
# typically noisier than forward reads (signal quality degrades further into
# a longer sequencing run), so a single shared model would fit neither well.
# Each .rds file holds what learnErrors returned, unchanged; dada() reads the
# final table from it. The verdict line says whether each model converged or was
# accepted by its code (see check_error_model).
message("[dada_seqtab] error models: ", readLines(snakemake@input[["verdict"]]))
err_r1 <- readRDS(snakemake@input[["model_r1"]])
err_r2 <- readRDS(snakemake@input[["model_r2"]])

# ── Dereplicate reads ────────────────────────────────────────────
message("[dada_seqtab] Dereplicating...")
derep_r1 <- derepFastq(filt_r1, verbose = TRUE)
derep_r2 <- derepFastq(filt_r2, verbose = TRUE)
# A run with ONE sample. derepFastq() here, and dada() and mergePairs() below, give
# back a list with one result per sample, named after the samples. With a single
# sample, each of them returns that sample's result on its own instead, not inside a
# list (each function does this itself; dada() and mergePairs() do it even when they
# are handed a one-element list). Left like that, the lone results of dada() and
# mergePairs() break two steps further down:
#   - the read tracking at the end runs sapply() over the samples in dada_r1 and in
#     mergers; over a lone result it runs over that result's parts instead, and the
#     job stops with an error in getUniques();
#   - makeSequenceTable() takes the name of each table row from its list element;
#     the lone data frame of mergePairs() has no name, so the sample would be called
#     "1" in the two sequence tables written below, and in every table built from
#     them later.
# So right after each of the three calls, a lone result is put back into a
# one-element list named after the sample, and every later step works as in a run
# with many samples. Each guard tests what DADA2 returned (one derep object, one
# dada object, one data frame) rather than counting the samples, so it would stay
# right if a later DADA2 release returned a list for one sample too. With two or
# more samples DADA2 returns lists, and the three guards change nothing.
#
# The lone result of derepFastq() does no harm by itself: dada() and mergePairs()
# accept one and put it into a list themselves, and the read tracking never reads
# it. But DADA2 refuses to rename one, so the two names() lines below would only
# print the warning "derep-class objects cannot be renamed" and leave it without
# the sample's name. Putting it into a list first lets them name it, as in a run
# with many samples.
if (inherits(derep_r1, "derep")) derep_r1 <- list(derep_r1)
if (inherits(derep_r2, "derep")) derep_r2 <- list(derep_r2)
names(derep_r1) <- sample_names
names(derep_r2) <- sample_names

# ── Denoise with DADA2 ───────────────────────────────────────────
# pool controls how much samples inform each other during denoising:
#   FALSE  (config default) — each sample denoised on its own. A variant only
#          seen in one sample can still be called if it's abundant enough
#          within that sample, but a variant that's rare in every sample
#          individually — even if it recurs across many of them — can be missed.
#   TRUE   — all samples pooled before denoising, so a variant rare in each
#            sample but recurring across many gets enough combined evidence to
#            be called. More sensitive to that kind of variant, but slower and
#            more memory-hungry as sample count grows.
#   "pseudo" — a middle ground: samples are still denoised individually, but
#            using a prior built from a first, pooled pass, at a fraction of
#            full pooling's cost.
# With a single sample there is nothing to pool: dada() then denoises it on its
# own whatever pool says (the log line below still prints the configured value).
message("[dada_seqtab] Denoising (pool=", pool, ", OMEGA_A=", omega_a, ")...")
dada_r1 <- dada(derep_r1, err = err_r1, multithread = threads,
                 pool = pool, OMEGA_A = omega_a, verbose = TRUE)
dada_r2 <- dada(derep_r2, err = err_r2, multithread = threads,
                 pool = pool, OMEGA_A = omega_a, verbose = TRUE)
# One sample: dada() returns one dada object, not a list; put it back into a
# one-element list named after the sample (see "A run with ONE sample" above).
if (inherits(dada_r1, "dada")) dada_r1 <- setNames(list(dada_r1), sample_names)
if (inherits(dada_r2, "dada")) dada_r2 <- setNames(list(dada_r2), sample_names)

# ── Merge paired reads ───────────────────────────────────────────
# trimOverhang handles reads that run past the far end of the amplicon into
# the OTHER primer's binding site — when that happens R1 and R2 extend past
# each other in the alignment ("overhang"); trimOverhang cuts that excess off
# rather than keeping it (which could leave primer leftovers in the merged
# sequence) or rejecting the pair outright.
message("[dada_seqtab] Merging pairs (minOverlap=", min_overlap,
        ", maxMismatch=", max_mismatch, ", justConcatenate=", just_concat, ")...")
mergers <- mergePairs(
  dada_r1, derep_r1, dada_r2, derep_r2,
  minOverlap      = min_overlap,
  maxMismatch     = max_mismatch,
  justConcatenate = just_concat,
  trimOverhang    = trim_overhang,
  verbose         = TRUE
)
# One sample: mergePairs() returns one data frame, not a list of them; put it back
# into a one-element list named after the sample, which makeSequenceTable() below
# uses as the row name (see "A run with ONE sample" above).
if (is.data.frame(mergers)) mergers <- setNames(list(mergers), sample_names)

# ── Build the sequence table ─────────────────────────────────────
# Columns (ASVs) come back ordered by decreasing total abundance across
# samples — makeSequenceTable's own default ordering; this call doesn't set
# it explicitly.
seqtab <- makeSequenceTable(mergers)
message("[dada_seqtab] Sequence table: ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs")
message("[dada_seqtab] Amplicon length distribution:")
print(table(nchar(getSequences(seqtab))))

# ── Remove chimeras ──────────────────────────────────────────────
# method="consensus" (config default): each sample is screened for chimeras
# on its own, then an ASV is only called chimeric overall if enough samples
# independently flagged it — a vote across samples, more conservative than
# pooling every sample together first and screening once ("pooled" method).
# min_fold_parent sets how much more abundant a candidate "parent" sequence
# must be than the ASV being tested before it's even considered as a possible
# source of that ASV. This run's config value (4) is stricter than DADA2's
# own default for the consensus method (1.5) — fewer ASVs get flagged as
# chimeric on a weak abundance signal alone.
# allow_one_off toggles whether ASVs one mismatch/indel away from an exact
# two-parent match are ALSO caught, not just exact matches; this run leaves it
# at DADA2's own default — see config.yaml for the config-level rationale.
message("[dada_seqtab] Removing chimeras (method=", chim_method,
        ", minFoldParent=", min_fold, ", allowOneOff=", allow_one_off, ")...")
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method                     = chim_method,
  minFoldParentOverAbundance = min_fold,
  allowOneOff                = allow_one_off,
  multithread                = threads,
  verbose                    = TRUE
)
chimera_frac <- 1 - sum(seqtab_nochim) / sum(seqtab)
message("[dada_seqtab] Chimeric read fraction: ", round(chimera_frac * 100, 2), "%")
message("[dada_seqtab] Final: ", ncol(seqtab_nochim), " non-chimeric ASVs")

# ── Build per-sample read tracking table ─────────────────────────
# Helper: total reads surviving a DADA2 step for one sample
getN <- function(x) sum(getUniques(x))

filter_stats <- do.call(rbind, lapply(filter_stats_files, function(f) {
  j <- fromJSON(f)
  data.frame(sample = j$sample, reads_in = j$reads_in, reads_out = j$reads_out,
             stringsAsFactors = FALSE)
}))

# denoised is expected to sit close to filtered: dada() reassigns each read to
# its most likely true ASV, it does not discard reads on its own, so a big gap
# here would point to a problem at the denoising step rather than filtering.
# merged can drop below denoised, though: mergePairs rejects any pair whose R1
# and R2 don't overlap by minOverlap bases with at most maxMismatch
# differences, and (with returnRejects left at its FALSE default here) those
# rejected pairs are simply left out of mergers rather than kept and flagged.
track <- data.frame(
  sample       = sample_names,
  stripped     = filter_stats$reads_in[match(sample_names,  filter_stats$sample)],
  filtered     = filter_stats$reads_out[match(sample_names, filter_stats$sample)],
  denoised     = as.integer(sapply(dada_r1, getN)),
  merged       = as.integer(sapply(mergers, getN)),
  non_chimeric = as.integer(rowSums(seqtab_nochim)),
  stringsAsFactors = FALSE
)
rownames(track) <- track$sample
# Drop the sample column before writing: it is already encoded as row names
write.table(track[, -1L], file = snakemake@output[["read_counts"]], sep = "\t", col.names = NA)
message("[dada_seqtab] Read counts written.")

# ── Write outputs ────────────────────────────────────────────────
# IDs are assigned in the sequence table's existing column order, which (per
# the makeSequenceTable note above) is decreasing total abundance — so ASV_1
# is this run's single most abundant ASV overall, ASV_2 the next, and so on.
asv_ids <- paste0("ASV_", seq_len(ncol(seqtab_nochim)))

# FASTA of unique ASV sequences
uniquesToFasta(getUniques(seqtab_nochim),
               fout = snakemake@output[["seqs_fasta"]], ids = asv_ids)

# Sequence table with raw sequences as column names
write.table(seqtab_nochim, file = snakemake@output[["seqtab_seqs"]],
            col.names = NA, sep = "\t")

# Sequence table with ASV IDs (ASV_1, ASV_2, ...) as column names
seqtab_named <- seqtab_nochim
colnames(seqtab_named) <- asv_ids
write.table(seqtab_named, file = snakemake@output[["seqtab_names"]],
            col.names = NA, sep = "\t")

message("[dada_seqtab] DONE. seqs.fasta, seqtab_head_seqs.txt, seqtab_head_names.txt written.")

# ── Close log ────────────────────────────────────────────────────
sink(type = "message")
sink(type = "output")
close(log_con)
