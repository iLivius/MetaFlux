#!/usr/bin/env Rscript
# Quality-filters and length-trims one sample's read pair before DADA2 denoising.
#
# What this is for: the cleanup step between primer trimming and DADA2's own steps:
# rule learn_errors learns the run's error models from the filtered reads of all
# samples, and rule dada_seqtab then denoises them. It removes reads that are
# too low-quality or the wrong length to trust, and cuts every surviving read
# to a fixed length so all reads entering the denoiser are directly comparable
# position-by-position (DADA2's error model is learned per-position, so mixed
# read lengths would misalign it).
#
# Input: this sample's primer-stripped FASTQ pair from 03.stripped/ (rule
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
# When no read pair of the sample passes, filterAndTrim writes no files at all;
# the script then stops with a message that names the settings that remove reads.
# Then it counts which quality values the filtered reads carry (the first 100,000
# reads of each file; count_quality_values below). A MiSeq writes some 30 distinct
# values, a NovaSeq or NextSeq only 3 to 6 ("binned" quality scores), and that decides
# which error function DADA2 should learn the error model with.
#
# Output: the filtered FASTQ pair (04.filtered/), which learn_errors learns the
# error models from before dada_seqtab denoises it (rule falco also profiles each
# file for the MultiQC report, stage "filtered"), plus two small JSON files:
#   - stats/dada2/{sample}.filter_stats.json: this sample's reads-in/reads-out counts,
#     which dada_seqtab reads back in to build the pipeline-wide read-tracking table
#     (the same file, field for field, as up to v2.6.1);
#   - stats/dada2/{sample}.quality_values.json: the quality values found, which
#     choose_error_model (60c_choose_error_model.py) pools over all samples to pick the
#     error function. A file of its own rather than a field of filter_stats.json: an
#     output folder made before this file existed lacks it, and Snakemake then runs
#     this step again there by itself, whatever --rerun-triggers says, because
#     choose_error_model, which every run asks for, cannot run without it.

# ── Load packages and redirect output to the log file ───────────
# ShortRead (in the dada2 environment, which DADA2 itself needs) reads the filtered
# FASTQ files in chunks for the quality count.
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, jsonlite, ShortRead))

# sink() below redirects everything this script prints from here on —
# messages, warnings, DADA2's own progress output — into the rule's log file
# instead of the console, so the log: file Snakemake tracks has the full record.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")


#' Stop the job with a text the user can read in full, in the log and on the terminal.
#'
#' R's own stop() prints at most about 1,000 characters of its message, so the full
#' text goes out with message() instead: once into the log (the sinks are still on),
#' then, after the sinks are closed, once more onto the terminal (on a cluster, into
#' the scheduler's log of the job). The short stop() line at the end makes the job
#' fail, which is what stops Snakemake. (60d_learn_errors.R carries the same function;
#' MetaFlux scripts share no code.)
#'
#' @param text The full message, starting with "[dada_filter] STOP:".
#' @param log_name The log of this job, for the short line, e.g. "logs/dada_filter/S1.log".
#' @return Never returns.
stop_with_message <- function(text, log_name) {
  message(text)                                       # full text into the log
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  message(text)                                       # full text onto the terminal
  stop("dada_filter stopped: reason above and in ", log_name, call. = FALSE)
}

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

# ── Per-sample read counts ───────────────────────────────────────
# out is filterAndTrim's own return value: a matrix with one row per input
# file pair and reads.in/reads.out columns. This script always passes exactly
# one sample's pair, so out always has exactly one row — [1L, ...] below.
reads_in  <- as.integer(out[1L, "reads.in"])
reads_out <- as.integer(out[1L, "reads.out"])
message("[dada_filter] reads in=", reads_in, " out=", reads_out,
        " (", round(100 * reads_out / max(reads_in, 1L), 1), "% retained)")

# ── Stop when no read pair passed ────────────────────────────────
# When no read pair passes, filterAndTrim writes neither filtered file ("The filter
# removed all reads: ... not written") and removes any it had started. Without this
# check the job would end with Snakemake's bare "missing output files" error, which
# names no cause; this names the settings that remove reads, so the user can either
# leave the sample out or change them.
if (!file.exists(filt_r1) || !file.exists(filt_r2)) {
  max_ee_text  <- paste0("[", paste(max_ee, collapse = ", "), "]")
  max_len_text <- if (is.infinite(max_len)) "null (no limit)" else as.character(max_len)
  log_name <- file.path("logs", "dada_filter", paste0(sample, ".log"))
  stop_with_message(sprintf(paste0(
    "[dada_filter] STOP: filterAndTrim kept none of the %d read pairs of sample %s, so ",
    "it wrote no filtered files. The settings that remove reads: truncLen R1=%d, R2=%d ",
    "(a read shorter than truncLen is dropped), max_ee %s, trunc_q %d, min_len %d, ",
    "max_len %s. Leave the sample out of the run, or check these settings; ",
    "filterAndTrim's own lines are in %s."),
    reads_in, sample, trunc_r1, trunc_r2, max_ee_text, trunc_q, min_len, max_len_text,
    log_name), log_name)
}

# ── Count the quality values the filtered reads carry ────────────
# Why: DADA2 learns its error model from the quality scores of the reads, and the
# function that fits it should suit the kind of quality values they carry. A MiSeq
# writes continuous values (filtered MiSeq reads carry 28-32 distinct values, roughly
# 7 to 38), for which DADA2's usual loess fit is made; a NovaSeq, a NextSeq 1000/2000 or
# another recent instrument writes "binned" values (only 3 to 6, for example 9, 20 and
# 34), for which DADA2's author recommends a binned fit. The values differ between
# instruments, their control software and their kits, so they are read from the reads
# themselves, never from a table, and from the filtered reads, which are exactly what
# learn_errors learns from (whatever trunc_q and max_ee did has already happened).
# Rule choose_error_model pools the counts of all samples and decides.
# How many reads: the first 100,000 of each file. Every distinct value is visible
# after about 1,000 reads, and 100,000 reads of 2 x 300 bp take a few seconds (3 to 8 s
# measured on binned NextSeq 2000 files); a file with fewer reads is read whole.
QUALITY_SCAN_READS <- 100000L


#' Count how many bases carry each quality value in the first reads of a FASTQ file.
#'
#' Reads the file in chunks with ShortRead's FastqStreamer, so that a large file is
#' never held in memory whole. qualityType = "Auto" works out how the quality letters
#' are encoded (Phred+33 in every current FASTQ) exactly as DADA2's derepFastq does,
#' and as(quality(chunk), "matrix") turns each read's letters into Phred values, one
#' row per read; reads shorter than the longest in the chunk are padded with NA, which
#' is dropped. Draws no random numbers (checked: R's random-number state is the same
#' before and after). learnErrors runs in another job (learn_errors), which sets the
#' seed itself, so nothing here could change the samples it draws anyway.
#'
#' @param fastq_path A FASTQ file, gzipped or not.
#' @param max_reads How many reads to look at, from the start of the file.
#' @param chunk_size How many reads to read at once.
#' @return list(reads_scanned = the number of reads looked at; values = the distinct
#'   quality values found, ascending; bases_per_value = the number of bases at each of
#'   those values, in the same order).
count_quality_values <- function(fastq_path, max_reads = QUALITY_SCAN_READS, chunk_size = 10000L) {
  # One counter per Phred value 0..93 (Phred+33 letters reach value 93 at most);
  # tabulate() counts value q in slot q + 1.
  bases_per_phred <- numeric(94)
  reads_scanned <- 0L
  streamer <- FastqStreamer(fastq_path, n = chunk_size)
  on.exit(close(streamer))
  while (reads_scanned < max_reads) {
    chunk <- yield(streamer, qualityType = "Auto")
    if (length(chunk) == 0) {
      break                                            # end of the file
    }
    reads_left <- max_reads - reads_scanned
    if (length(chunk) > reads_left) {
      chunk <- chunk[seq_len(reads_left)]
    }
    phred <- as(quality(chunk), "matrix")
    phred <- phred[!is.na(phred)]
    bases_per_phred <- bases_per_phred + tabulate(phred + 1L, nbins = 94L)
    reads_scanned <- reads_scanned + length(chunk)
  }
  values <- which(bases_per_phred > 0) - 1L
  list(reads_scanned = reads_scanned,
       values = as.integer(values),
       bases_per_value = bases_per_phred[values + 1L])
}


#' The quality values of one file, as a short text for the log: "3 values (9, 20,
#' 34)", or "32 values (7-38)" when there are more than ten.
#'
#' @param counted What count_quality_values returned.
#' @return Character string.
quality_values_log_text <- function(counted) {
  values <- counted$values
  if (length(values) == 0) {
    return("no values")
  }
  if (length(values) > 10) {
    return(sprintf("%d values (%d-%d)", length(values), min(values), max(values)))
  }
  sprintf("%d values (%s)", length(values), paste(values, collapse = ", "))
}


#' One read direction's quality values, as they go into quality_values.json below.
#'
#' Vectors are wrapped in I() so that a single value stays a JSON list; the bases are
#' written per value, as {"9": 2617, "20": 5559, ...}.
#'
#' @param counted What count_quality_values returned.
#' @return A named list.
quality_values_record <- function(counted) {
  bases <- as.list(counted$bases_per_value)
  names(bases) <- as.character(counted$values)
  list(reads_scanned = counted$reads_scanned,
       values = I(counted$values),
       bases_per_value = bases)
}


quality_r1 <- count_quality_values(filt_r1)
quality_r2 <- count_quality_values(filt_r2)
message("[dada_filter] quality values in the first ", format(QUALITY_SCAN_READS, big.mark = ","),
        " filtered reads: R1 ", quality_values_log_text(quality_r1),
        ", R2 ", quality_values_log_text(quality_r2))

# ── Write the two per-sample JSON files ──────────────────────────
# filter_stats.json: reads_in and reads_out, for the read-tracking table (dada_seqtab).
# Written with exactly the call of v2.6.1, so the file keeps its layout byte for byte.
write(
  toJSON(list(sample = sample, reads_in = reads_in, reads_out = reads_out),
         auto_unbox = TRUE, pretty = TRUE),
  snakemake@output[["stats"]]
)
# quality_values.json: the quality values of the filtered reads, for the choice of
# error function (choose_error_model). scan_limit records how many reads per file were
# looked at, so the counts can be read correctly later. digits = NA: numbers are
# written at full precision, never rounded to a few decimals (the counts are whole
# numbers anyway).
write(
  toJSON(list(sample = sample,
              scan_limit = QUALITY_SCAN_READS,
              R1 = quality_values_record(quality_r1),
              R2 = quality_values_record(quality_r2)),
         auto_unbox = TRUE, pretty = TRUE, digits = NA),
  snakemake@output[["quality_values"]]
)

# ── Close log ────────────────────────────────────────────────────
sink(type = "message")
sink(type = "output")
close(log_con)
