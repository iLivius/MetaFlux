#!/usr/bin/env Rscript
# Quality-profile diagnostic plot for primer-trimmed reads, R1 and R2 separately.
#
# What this is for: a visual check of how base-call quality changes along the
# read, computed AFTER primer removal but BEFORE any DADA2 filtering. It is
# the human-readable counterpart of the falco quality tables that rule
# pick_trunclen (script 50a_pick_trunclen.py) uses to choose truncLen
# automatically — this script does not feed any decision back into the
# pipeline, it is only for a human to eyeball a run.
#
# Input: the primer-stripped FASTQ pair for every sample, from 3.stripped/
# (produced by rule trim_primers, Cutadapt, see 30_preprocess.smk). All
# samples' R1 files go into one plot, all R2 files into another.
#
# What it does: dada2::plotQualityProfile() reads the FASTQ files and, for
# each position along the read, plots the distribution of quality scores as a
# grey heat map (darker = more reads carry that score at that position),
# overlaid with a green mean line and an orange median line (dashed orange =
# 25th/75th percentile). If reads in the set vary in length, a red line shows
# what percentage of reads are still that long at each position — handy for
# spotting how much of the library actually reaches a given base.
# aggregate=TRUE pools every sample into one profile per direction instead of
# drawing a separate panel per file, which would be unreadable at real sample
# counts.
#
# Output: one PNG + one PDF per direction (R1, R2) under stats/dada2/. Nothing
# downstream reads these back in, and they are not even picked up by the
# MultiQC report (MultiQC only scans stats/falco and stats/cutadapt in
# amplicon mode — see 90_report.smk): open them by hand when sanity-checking a
# run or deciding whether trunc_len needs adjusting.

# ── Load packages and redirect output to the log file ───────────
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, ggplot2))

# sink() below redirects everything this script prints from here on —
# messages, warnings, DADA2's own progress output — into the rule's log file
# instead of the console, so the log: file Snakemake tracks has the full record.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# ── Read input file paths from Snakemake ─────────────────────────
# Sorted purely for a stable, reproducible file order in the log lines below.
# With aggregate=TRUE the plot itself pools every file together, so plotting
# order has no effect on the result.
fn_r1 <- sort(unlist(snakemake@input[["fnFs"]]))
fn_r2 <- sort(unlist(snakemake@input[["fnRs"]]))

message("[dada_quality_plots] R1 files: ", length(fn_r1))
message("[dada_quality_plots] R2 files: ", length(fn_r2))

dir.create(dirname(snakemake@output[["plot_r1"]]), recursive = TRUE, showWarnings = FALSE)

# ── Generate and save quality plots ─────────────────────────────
# n = 1e6 caps how many reads plotQualityProfile samples per file before
# building the plot (a memory guard on huge files) — well above any realistic
# per-sample amplicon read count, so in practice every read is used. Each
# profile is saved twice: PNG for quick viewing, PDF as a scalable copy for a
# report.
p_r1 <- plotQualityProfile(fn_r1, n = 1e6, aggregate = TRUE)
ggsave(snakemake@output[["plot_r1"]],     plot = p_r1, device = "png", width = 30, height = 20, units = "cm", dpi = 150)
ggsave(snakemake@output[["plot_r1_pdf"]], plot = p_r1, device = "pdf", width = 30, height = 20, units = "cm")
message("[dada_quality_plots] Saved R1 plot: ", snakemake@output[["plot_r1"]])

# Same as above, R2 direction.
p_r2 <- plotQualityProfile(fn_r2, n = 1e6, aggregate = TRUE)
ggsave(snakemake@output[["plot_r2"]],     plot = p_r2, device = "png", width = 30, height = 20, units = "cm", dpi = 150)
ggsave(snakemake@output[["plot_r2_pdf"]], plot = p_r2, device = "pdf", width = 30, height = 20, units = "cm")
message("[dada_quality_plots] Saved R2 plot: ", snakemake@output[["plot_r2"]])

# ── Close log ────────────────────────────────────────────────────
sink(type = "message")
sink(type = "output")
close(log_con)
