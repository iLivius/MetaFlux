#!/usr/bin/env Rscript
# DADA2's core denoising chain: turns filtered reads into a table of exact
# amplicon sequence variants (ASVs) — DADA2's alternative to OTU clustering,
# where every ASV is one exact sequence rather than a cluster of similar ones.
#
# Input: every sample's filtered FASTQ pair from 4.filtered/ (one pair per
# sample, produced by dada_filter.R) plus that sample's filter_stats JSON
# (reads in/out from the filtering step), collected here across the whole
# sample set for the final read-tracking table.
#
# What it does, in order:
#   1. learnErrors  — learns this run's actual sequencing-error profile (how
#      often, say, a true A gets read as a G, as a function of quality score)
#      from the filtered reads themselves, separately for R1 and R2, by
#      alternating between calling sequence variants and re-estimating error
#      rates from them until the estimate stops changing much. DADA2 needs
#      this error model in the next step to tell a real biological sequence
#      apart from one that is just a sequencing mistake away from another.
#   2. derepFastq   — dereplication: collapses each sample's reads down to
#      its unique sequences, each carrying a count (abundance). DADA2 works
#      on these unique sequences rather than every individual read, which is
#      what keeps it fast on large samples.
#   3. dada         — the actual denoising step: using the error model from
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
# names, and short ASV_N ids), a FASTA of the ASV sequences, a per-sample
# read-tracking table (stripped -> filtered -> denoised -> merged ->
# non_chimeric), and diagnostic error-model plots. seqs.fasta and the ASV_N
# seqtab feed dada_length_filter (60d) next, which applies the post-hoc ASV
# length window; assign_taxonomy (80a) consumes the length-filtered versions
# of both further downstream.

# ── Load packages and redirect output to the log file ───────────
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, jsonlite, ggplot2))

# sink() below redirects everything this script prints from here on —
# messages, warnings, DADA2's own progress output — into the rule's log file
# instead of the console, so the log: file Snakemake tracks has the full record.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# ── Read input file paths from Snakemake ─────────────────────────
# Both lists are sorted independently, by full path. Because an R1 and R2 path
# for the same sample differ only in that "R1"/"R2" substring (same position,
# same length), sorting each list alphabetically still lines up corresponding
# samples at the same index in both — which matters here because sample_names
# below is derived once from filt_r1's filenames and then reused for filt_r2
# by POSITION, not by matching filename.
filt_r1          <- sort(unlist(snakemake@input[["r1"]]))
filt_r2          <- sort(unlist(snakemake@input[["r2"]]))
filter_stats_files <- unlist(snakemake@input[["filter_stats"]])

# Recover each sample's name from its filtered-R1 filename (dada_filter.R's own
# naming convention: "{sample}_R1_filt.fastq.gz") and attach it to both file
# vectors — dada()/mergePairs() use these names to label every result by sample.
sample_names <- sub("_R1_filt\\.fastq\\.gz$", "", basename(filt_r1))
names(filt_r1) <- sample_names
names(filt_r2) <- sample_names

message("[dada_seqtab] ", length(sample_names), " sample(s): ", paste(sample_names, collapse = ", "))

# ── Read parameters from Snakemake ───────────────────────────────
pool          <- snakemake@params[["pool"]]
nbases        <- as.numeric(snakemake@params[["learn_errors_nbases"]])
max_consist   <- as.integer(snakemake@params[["learn_errors_max_consist"]])
omega_a       <- as.numeric(snakemake@params[["omega_a"]])
min_overlap   <- as.integer(snakemake@params[["merge_min_overlap"]])
max_mismatch  <- as.integer(snakemake@params[["merge_max_mismatch"]])
just_concat   <- as.logical(snakemake@params[["merge_just_concatenate"]])
trim_overhang <- as.logical(snakemake@params[["merge_trim_overhang"]])
chim_method   <- snakemake@params[["chimera_method"]]
min_fold      <- as.numeric(snakemake@params[["chimera_min_fold_parent"]])
allow_one_off <- as.logical(snakemake@params[["chimera_allow_one_off"]])
threads       <- as.integer(snakemake@threads)
seed          <- as.integer(snakemake@params[["seed"]])

# learnErrors(randomize=TRUE) draws its nbases sample of reads via R's RNG;
# fix it so the error model (and everything downstream) is reproducible.
set.seed(seed)

# YAML false/true → Python False/True → R FALSE/TRUE (logical); "pseudo" → character.
# The is.character guard converts string "false"/"true" if YAML used quotes.
if (is.character(pool)) {
  pool <- switch(tolower(trimws(pool)), "false" = FALSE, "true" = TRUE, "pseudo" = "pseudo", FALSE)
}

# ── Create output directories ────────────────────────────────────
dir.create(dirname(snakemake@output[["seqtab_seqs"]]),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake@output[["err_plot_r1"]]), recursive = TRUE, showWarnings = FALSE)

# ── Learn the error model from filtered reads ────────────────────
# R1 and R2 get separate error models — on Illumina, reverse reads are
# typically noisier than forward reads (signal quality degrades further into
# a longer sequencing run), so a single shared model would fit neither well.
message("[dada_seqtab] Learning error rates (R1, nbases=", nbases, ")...")
err_r1 <- learnErrors(filt_r1, nbases = nbases, multithread = threads,
                       randomize = TRUE, MAX_CONSIST = max_consist, verbose = TRUE)

message("[dada_seqtab] Learning error rates (R2, nbases=", nbases, ")...")
err_r2 <- learnErrors(filt_r2, nbases = nbases, multithread = threads,
                       randomize = TRUE, MAX_CONSIST = max_consist, verbose = TRUE)

# ── Save error-model diagnostic plots ───────────────────────────
# plotErrors shows, for each possible base transition (e.g. A misread as G),
# the observed error frequency by quality score (points) against the learned
# model (line). nominalQ=TRUE adds a reference line for what the error rate
# WOULD be if quality scores were taken at face value — comparing the learned
# line to that reference is a quick check that error learning actually
# converged to something sensible for this run, rather than something odd.
# log10 warnings from zero-count quality bins are cosmetic — suppressed.
p_err_r1 <- suppressWarnings(plotErrors(err_r1, nominalQ = TRUE))
ggsave(snakemake@output[["err_plot_r1"]],     plot = p_err_r1, device = "png", width = 30, height = 20, units = "cm", dpi = 150)
ggsave(snakemake@output[["err_plot_r1_pdf"]], plot = p_err_r1, device = "pdf", width = 30, height = 20, units = "cm")
p_err_r2 <- suppressWarnings(plotErrors(err_r2, nominalQ = TRUE))
ggsave(snakemake@output[["err_plot_r2"]],     plot = p_err_r2, device = "png", width = 30, height = 20, units = "cm", dpi = 150)
ggsave(snakemake@output[["err_plot_r2_pdf"]], plot = p_err_r2, device = "pdf", width = 30, height = 20, units = "cm")
message("[dada_seqtab] Error plots written.")

# ── Dereplicate reads ────────────────────────────────────────────
message("[dada_seqtab] Dereplicating...")
derep_r1 <- derepFastq(filt_r1, verbose = TRUE)
derep_r2 <- derepFastq(filt_r2, verbose = TRUE)
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
message("[dada_seqtab] Denoising (pool=", pool, ", OMEGA_A=", omega_a, ")...")
dada_r1 <- dada(derep_r1, err = err_r1, multithread = threads,
                 pool = pool, OMEGA_A = omega_a, verbose = TRUE)
dada_r2 <- dada(derep_r2, err = err_r2, multithread = threads,
                 pool = pool, OMEGA_A = omega_a, verbose = TRUE)

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
