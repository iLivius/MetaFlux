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

# min_len: filterAndTrim's min_len is a PER-READ length floor. It normally comes
# straight from the config (DADA2's own default is 20).
#
# ITS only, and only when min_len_stat is set: the floor can instead be derived from
# the probe distribution, which for ITS is the extracted ITS1/ITS2 subregion length.
# That is OPT-IN because it is a trade-off, not a free improvement: the UNITE
# reference bottoms out at 17 bp (ITS1) / 38 bp (ITS2), so a stat like q1 (~175 bp)
# sits above roughly a quarter of known fungal ITS diversity and discards reads from
# genuinely short amplicons — the same bias that setting truncLen to c(0,0) for ITS
# exists to avoid. Leave min_len_stat null to keep the config floor for every marker.
#
# The other markers never take this path: their probe measures the full amplicon
# (~400 bp for 16S V3-V4), which exceeds per-read length and would drop every read.
probe_json_path <- snakemake@input[["probe_json"]]
use_probe_min_len <- (
  amp_type == "ITS" &&
  !is.null(min_len_stat) && !is.na(min_len_stat) && nzchar(as.character(min_len_stat)) &&
  length(probe_json_path) > 0L && file.exists(probe_json_path[[1]])
)
if (use_probe_min_len) {
  probe   <- fromJSON(probe_json_path[[1]])
  min_len <- as.integer(probe[[min_len_stat]])
  message("[dada_filter] min_len auto-set from probe ", min_len_stat, "=", min_len, " bp (ITS)")
} else {
  min_len <- as.integer(snakemake@params[["min_len"]])
  message("[dada_filter] min_len from config: ", min_len, " bp")
}

# ITS: disable the fixed-length truncation. In DADA2, truncLen does two things —
# it cuts every read to a set length AND discards any read shorter than that.
# ITS amplicons vary in length, so a fixed cut would throw away genuinely short
# sequences. The 3' low-quality tail is NOT left untrimmed, though: truncQ in the
# filterAndTrim call below still trims each read adaptively at its first low-quality
# base, and maxEE drops reads with too many expected errors.
# This override is unconditional, so it also wins over manual mode: manual_r1 /
# manual_r2 have no effect for ITS. (16S keeps truncLen as picked by
# pick_trunclen, auto or manual.)
if (amp_type == "ITS") {
  trunc_r1 <- 0L
  trunc_r2 <- 0L
  message("[dada_filter] ITS mode: truncLen overridden to c(0, 0); 3' quality handled by truncQ/maxEE")
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
