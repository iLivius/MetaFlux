#!/usr/bin/env Rscript
# DADA2 core: learnErrors → dada → mergePairs → makeSequenceTable → removeBimeraDenovo.
# Reads all per-sample filtered FASTQs produced by dada_filter.R.
# Emits: seqtab (two flavours), seqs.fasta, read.counts, error-model diagnostic plots.

# ── Load packages and redirect output to the log file ───────────
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, jsonlite, ggplot2))

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# ── Read input file paths from Snakemake ─────────────────────────
filt_r1          <- sort(unlist(snakemake@input[["r1"]]))
filt_r2          <- sort(unlist(snakemake@input[["r2"]]))
filter_stats_files <- unlist(snakemake@input[["filter_stats"]])

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
message("[dada_seqtab] Learning error rates (R1, nbases=", nbases, ")...")
err_r1 <- learnErrors(filt_r1, nbases = nbases, multithread = threads,
                       randomize = TRUE, MAX_CONSIST = max_consist, verbose = TRUE)

message("[dada_seqtab] Learning error rates (R2, nbases=", nbases, ")...")
err_r2 <- learnErrors(filt_r2, nbases = nbases, multithread = threads,
                       randomize = TRUE, MAX_CONSIST = max_consist, verbose = TRUE)

# ── Save error-model diagnostic plots ───────────────────────────
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
message("[dada_seqtab] Denoising (pool=", pool, ", OMEGA_A=", omega_a, ")...")
dada_r1 <- dada(derep_r1, err = err_r1, multithread = threads,
                 pool = pool, OMEGA_A = omega_a, verbose = TRUE)
dada_r2 <- dada(derep_r2, err = err_r2, multithread = threads,
                 pool = pool, OMEGA_A = omega_a, verbose = TRUE)

# ── Merge paired reads ───────────────────────────────────────────
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
seqtab <- makeSequenceTable(mergers)
message("[dada_seqtab] Sequence table: ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs")
message("[dada_seqtab] Amplicon length distribution:")
print(table(nchar(getSequences(seqtab))))

# ── Remove chimeras ──────────────────────────────────────────────
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
