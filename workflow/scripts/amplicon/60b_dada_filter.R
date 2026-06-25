#!/usr/bin/env Rscript
# Per-sample DADA2 filterAndTrim.
# TruncLen is read from the pipeline's trunclen.json (produced by pick_trunclen.py).
# Outputs the filtered FASTQ pair and a per-sample JSON with read-in/out counts
# that dada_seqtab.R collects for the full pipeline read-tracking table.

# ── Load packages and redirect output to the log file ───────────
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, jsonlite))

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# ── Read inputs and parameters from Snakemake ────────────────────
sample <- snakemake@wildcards[["sample"]]
message("[dada_filter] sample=", sample)

tl      <- fromJSON(snakemake@input[["trunclen_json"]])
trunc_r1 <- as.integer(tl[["r1"]])
trunc_r2 <- as.integer(tl[["r2"]])
message("[dada_filter] truncLen: R1=", trunc_r1, " R2=", trunc_r2, " (mode=", tl[["mode"]], ")")

fn_r1   <- snakemake@input[["r1"]]
fn_r2   <- snakemake@input[["r2"]]
filt_r1 <- snakemake@output[["r1"]]
filt_r2 <- snakemake@output[["r2"]]

max_ee       <- unlist(snakemake@params[["max_ee"]])
trunc_q      <- as.integer(snakemake@params[["trunc_q"]])
amp_type     <- snakemake@params[["amp_type"]]
min_len_stat <- snakemake@params[["min_len_stat"]]

# min_len: filterAndTrim's min_len is a PER-READ length floor.
#   ITS: the probe distribution (extracted ITS2 lengths, ~140-260 bp) is
#        below per-read length, so auto-derive min_len from probe[min_len_stat].
#   16S: the probe distribution is the full in-silico PCR amplicon length
#        (~400 bp for V3-V4 or V5-V7), which exceeds per-read length and would
#        drop every read. Fall back to the config's manual min_len.
probe_json_path <- snakemake@input[["probe_json"]]
if (amp_type == "ITS" && length(probe_json_path) > 0L && file.exists(probe_json_path[[1]])) {
  probe   <- fromJSON(probe_json_path[[1]])
  min_len <- as.integer(probe[[min_len_stat]])
  message("[dada_filter] min_len auto-set from probe ", min_len_stat, "=", min_len, " bp (ITS)")
} else {
  min_len <- as.integer(snakemake@params[["min_len"]])
  message("[dada_filter] min_len from config: ", min_len, " bp")
}

# ITS reads are not hard-truncated: variable amplicon length means a fixed cut
# would discard genuinely short sequences. Quality is controlled by maxEE and
# minLen only. For 16S, truncLen stays as picked by pick_trunclen.
if (amp_type == "ITS") {
  trunc_r1 <- 0L
  trunc_r2 <- 0L
  message("[dada_filter] ITS mode: truncLen overridden to c(0, 0)")
}

# max_len is optional; if not set in config it arrives as NULL or NA → use Inf (no upper limit)
max_len_raw <- snakemake@params[["max_len"]]
max_len <- if (is.null(max_len_raw) || (length(max_len_raw) == 1L && is.na(max_len_raw))) Inf else as.numeric(max_len_raw)

dir.create(dirname(filt_r1), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake@output[["stats"]]), recursive = TRUE, showWarnings = FALSE)

# ── Filter and trim reads ────────────────────────────────────────
out <- filterAndTrim(
  fn_r1, filt_r1,
  fn_r2, filt_r2,
  truncLen   = c(trunc_r1, trunc_r2),
  maxN       = 0L,
  maxEE      = max_ee,
  truncQ     = trunc_q,
  minLen     = min_len,
  maxLen     = max_len,
  rm.phix    = FALSE,
  compress   = TRUE,
  multithread = as.integer(snakemake@threads),
  verbose    = TRUE
)

# ── Write per-sample read counts to JSON ─────────────────────────
reads_in  <- as.integer(out[1L, "reads.in"])
reads_out <- as.integer(out[1L, "reads.out"])
message("[dada_filter] reads in=", reads_in, " out=", reads_out,
        " (", round(100 * reads_out / max(reads_in, 1L), 1), "% retained)")

write(
  toJSON(list(sample = sample, reads_in = reads_in, reads_out = reads_out),
         auto_unbox = TRUE, pretty = TRUE),
  snakemake@output[["stats"]]
)

# ── Close log ────────────────────────────────────────────────────
sink(type = "message")
sink(type = "output")
close(log_con)
