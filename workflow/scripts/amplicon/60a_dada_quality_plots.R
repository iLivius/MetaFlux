#!/usr/bin/env Rscript
# DADA2 per-base quality profile plots for primer-trimmed (stripped) reads.
# Aggregates across all samples; one PNG and one PDF per direction.
# Diagnostic only — does not influence truncLen (that comes from pick_trunclen.py).

# ── Load packages and redirect output to the log file ───────────
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, ggplot2))

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# ── Read input file paths from Snakemake ─────────────────────────
fn_r1 <- sort(unlist(snakemake@input[["fnFs"]]))
fn_r2 <- sort(unlist(snakemake@input[["fnRs"]]))

message("[dada_quality_plots] R1 files: ", length(fn_r1))
message("[dada_quality_plots] R2 files: ", length(fn_r2))

dir.create(dirname(snakemake@output[["plot_r1"]]), recursive = TRUE, showWarnings = FALSE)

# ── Generate and save quality plots ─────────────────────────────
p_r1 <- plotQualityProfile(fn_r1, n = 1e6, aggregate = TRUE)
ggsave(snakemake@output[["plot_r1"]],     plot = p_r1, device = "png", width = 30, height = 20, units = "cm", dpi = 150)
ggsave(snakemake@output[["plot_r1_pdf"]], plot = p_r1, device = "pdf", width = 30, height = 20, units = "cm")
message("[dada_quality_plots] Saved R1 plot: ", snakemake@output[["plot_r1"]])

p_r2 <- plotQualityProfile(fn_r2, n = 1e6, aggregate = TRUE)
ggsave(snakemake@output[["plot_r2"]],     plot = p_r2, device = "png", width = 30, height = 20, units = "cm", dpi = 150)
ggsave(snakemake@output[["plot_r2_pdf"]], plot = p_r2, device = "pdf", width = 30, height = 20, units = "cm")
message("[dada_quality_plots] Saved R2 plot: ", snakemake@output[["plot_r2"]])

# ── Close log ────────────────────────────────────────────────────
sink(type = "message")
sink(type = "output")
close(log_con)
