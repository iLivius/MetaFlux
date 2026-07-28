#!/usr/bin/env Rscript
# Quality-filters and length-trims one sample's read pair before DADA2 denoising.
#
# What this is for: the cleanup step between primer trimming and DADA2's own
# denoising (rule dada_seqtab, next in the chain). It removes reads that are
# too low-quality or the wrong length to trust, and cuts every surviving read
# to a fixed length so all reads entering the denoiser are directly comparable
# position-by-position (DADA2's error model is learned per-position, so mixed
# read lengths would misalign it).
#
# Input: this sample's primer-stripped FASTQ pair from 3.stripped/ (rule
# trim_primers, Cutadapt, see 30_preprocess.smk), plus trunclen.json — the
# truncLen decision made once for the whole run by rule pick_trunclen
# (50a_pick_trunclen.py). This script runs once per sample (Snakemake
# wildcard {sample}), so it only ever sees one R1/R2 pair at a time.
#
# What it does: dada2::filterAndTrim() (the call near the bottom of this
# file) is DADA2's read-filtering function. In one pass over each read it:
# cuts the read to the truncLen decided above (a fixed-position 3' cut —
# anything shorter than that is dropped outright, not padded); trims further,
# adaptively, at the first base whose quality score falls to trunc_q or below;
# then drops the read if it now has more than max_ee "expected errors" (the
# sum of the per-base error probabilities implied by the quality scores — a
# read-wide error budget, not a single bad base) or falls outside the
# [min_len, max_len] length window.
#
# Output: the filtered FASTQ pair (4.filtered/) that dada_seqtab.R denoises
# next, plus a small JSON with this sample's reads-in/reads-out counts, which
# dada_seqtab.R reads back in to build the pipeline-wide read-tracking table.

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

# ── Read inputs and parameters from Snakemake ────────────────────
sample <- snakemake@wildcards[["sample"]]
message("[dada_filter] sample=", sample)

# trunclen.json holds ONE truncLen decision shared by the whole run (both
# directions), not a per-sample value — every sample calls filterAndTrim with
# the same truncLen so all ASVs later line up position-for-position.
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
# maxN is hardcoded to 0 below, not read from config: DADA2's denoiser cannot
# handle any ambiguous (N) bases at all, so this is a hard requirement of the
# method, not a tunable setting.
# rm.phix is hardcoded FALSE below, independent of amplicon.decontamination.
# remove_phix. When that flag is on, PhiX removal already happened upstream
# via bowtie2 (rule rm_phix, 30_preprocess.smk) before primer trimming even
# ran — an alignment against the actual PhiX genome, more sensitive than
# filterAndTrim's own built-in PhiX check — so repeating a weaker check here
# would add nothing. When remove_phix is off, no PhiX screening happens at any
# stage of the amplicon path; that is what turning the flag off means.
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
# out is filterAndTrim's own return value: a matrix with one row per input
# file pair and reads.in/reads.out columns. This script always passes exactly
# one sample's pair, so out always has exactly one row — [1L, ...] below.
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
