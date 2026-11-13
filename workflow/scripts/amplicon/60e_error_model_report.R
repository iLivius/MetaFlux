#!/usr/bin/env Rscript
# The report on the two DADA2 error models of the run: how closely each model follows
# the error rates it was learned from, how its learning ended, and which quality values
# the choice of its error function rested on. Made in every amplicon run, as tables,
# figures and a section of the MultiQC report.
#
# Why it exists: dada_seqtab denoises every sample with one error model per read
# direction (rule learn_errors, 60d_learn_errors.R), so a model that describes the reads
# badly affects every ASV table of the run. When a model did not converge, rule
# check_error_model stops the run before denoising and sends the user here: these files
# are the evidence for accepting the model or not. The check takes this report's
# error_model_report.json as an input, so the report is on disk before the check can
# stop the run, and it is drawn in every run, so that a converged model can be judged
# the same way.
#
# Input (rule error_model_report):
#   - record:   stats/dada2/error_model.json, written by learn_errors. Per read
#               direction: the error function used and why, the quality values of the
#               voting samples, the anchors of a binned fit, the rounds, how the
#               learning ended (stopped_by), the change in every round, the distance of
#               every round's starting table to the final one, the model code, and the
#               transition counts the final table was fitted to (counts_behind_model).
#   - model_r1, model_r2: 05.dada2/error_model_R1.rds and _R2.rds, what learnErrors
#               returned (the final table err_out, the starting table of every round
#               err_in, the counts of the last round trans).
#   - choice:   stats/dada2/error_model_choice.json, written by choose_error_model: per
#               direction the bases at each quality value in the first 100,000 filtered
#               reads of each voting sample (the quality-values figure).
#   - error_functions, vendored_binned: the registry of the error functions users can
#               choose (60_error_functions.R) and DADA2's binned fit that the registry
#               calls (vendor/dada2_makeBinnedQualErrfun.R). Both are sourced below.
#   - custom_function_file, with the params custom_function_name and
#               custom_function_file_as_configured: the R file of the user's own error
#               function (EXPERIMENTAL), whenever amplicon.dada2.learn_errors.
#               custom_error_function names one (empty otherwise). It is refitted and
#               drawn like the registry's functions, labelled "custom (<name>)": as the
#               model used under error_model: custom, as a comparison under any other
#               setting.
#
# What it does, for each read direction:
#   1. Rebuilds the transition counts the model was fitted to: at each quality score,
#      how many bases of each true nucleotide were read as A, C, G or T in the reads
#      learnErrors drew, after its last round of denoising (for a model ended by the
#      flat-rounds rule, the counts of the round before the last, which the kept table
#      was fitted to). When they should equal DADA2's own counts (trans), it checks so.
#   2. Refits every function of the registry on those counts, once: only the smoothing
#      step, which turns counts into rates, not a new learning (which would take hours
#      and give other counts). The refit of the function used must give back the model
#      exactly; the refits of the others show what they would draw through the same
#      points. The statement REPORT_STATEMENT below says what this comparison can and
#      cannot show; it is repeated in the report JSON, in the MultiQC section and in the
#      documentation.
#   3. Measures how closely each fitted table follows the counted rates, separately at
#      the quality values the reads carry ("read values") and at the other columns that
#      hold bases ("averaged columns": DADA2 gives a sequence seen in several reads the
#      mean quality of those reads at each position, rounded, so binned reads with only
#      9, 20 and 34 produce counts at many values in between; most of their bases come
#      from reads identical to their ASV, and how many of the misreads they hold depends
#      on the data, so the notes of the fit table give both shares). The measures are
#      defined at goodness_of_fit() below: the deviance of the counts under the fitted
#      rates, the typical fold difference between counted and fitted rates, and how
#      often a fitted rate rises with quality. No measure ranks the functions: each one
#      wins on what it is built for.
#   4. Writes the tables, draws the figures, and writes the files of the MultiQC section.
#
# Output (22 files under stats/dada2/; the list is _error_model_report_files() in
# 00_common.smk, which is also the rule's output list):
#   - error_model_fit.tsv:    one row per direction and error function, the measures of
#                             step 3 (columns described at fit_table_columns below);
#   - error_model_rates.tsv:  one row per direction, kind of misread and quality score:
#                             bases, errors, counted rate, nominal rate of the quality
#                             score, the model's rate and each refit's rate;
#   - error_model_rounds.tsv: one row per direction and round: the change in the round,
#                             the lowest change so far, the distance to the final model,
#                             and where the learning stopped;
#   - error_model_report.json: the statement, per direction how the learning ended, the
#                             model code, the shares of bases and errors at the read values
#                             and the averaged columns, the rows of error_model_fit.tsv, and
#                             the warnings of this report;
#   - error_model_fit_R1.png/.pdf and _R2: counted rates, the model and the refits, per
#                             kind of misread, against the quality score;
#   - error_model_convergence.png/.pdf: the change per round and the distance to the
#                             final model, with DADA2's stop, the flat-rounds rule and the
#                             cap marked;
#   - error_model_quality_values.png/.pdf: the bases at each quality value of the
#                             filtered reads, behind the choice of error function;
#   - filtered_read_R1_error_plot.png/.pdf and _R2: DADA2's own plotErrors picture of each
#                             model, as MetaFlux has drawn it since v1.0.0;
#   - multiqc/error_model_summary_mqc.tsv, multiqc/error_model_goodness_of_fit_mqc.tsv
#                             and PNG copies of the four figures (multiqc/*_mqc.png): the
#                             "DADA2 error model" section of the MultiQC report (MultiQC
#                             needs files of its own in a folder it scans; the section's
#                             texts and order are in workflow/multiqc_config.yaml).
# Who reads them: check_error_model takes error_model_report.json, only to wait for the
# report and to repeat its warnings; rule multiqc reads the six files in multiqc/;
# nothing else reads any of them.
#
# Robustness: every figure and every table is made inside its own tryCatch. A figure
# that cannot be drawn is replaced by a picture that states the error, and every such
# failure becomes a warning in error_model_report.json, in the Warnings column of the
# MultiQC summary and in the log (check_error_model repeats it as a note). The job fails
# only when the record or a model cannot be read. The reason: the check waits for this
# report, so a report that failed would block the denoising with no setting to get past
# it.
#
# Determinism: the TSV, JSON and PNG files hold no date, time or run time, so the same
# models give byte-identical files and a rerun does not count as a change for a later
# run. A PDF embeds its creation date, so no rule reads a PDF.

# ── Load packages and redirect output to the log file ───────────
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, jsonlite, ggplot2))

# sink() below redirects everything this script prints from here on into the rule's
# log file instead of the console, so the log: file Snakemake tracks has the record.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
# A warning goes into the log where it happens, next to the step that raised it,
# instead of being collected until the end of the script.
options(warn = 1)


#' Stop the job with a text the user can read in full, in the log and on the terminal.
#'
#' Used only when the record or a model cannot be read: everything else that goes wrong
#' becomes a warning and a placeholder (see "Robustness" in the header). R's own stop()
#' prints at most about 1,000 characters of its message, so the full text goes out with
#' message(), once into the log and, after the sinks are closed, once onto the terminal.
#'
#' @param text The full message, starting with "[error_model_report] STOP:".
#' @return Never returns.
stop_with_message <- function(text) {
  message(text)                                       # full text into the log
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  message(text)                                       # full text onto the terminal
  stop("error_model_report stopped: reason above and in logs/error_model_report.log", call. = FALSE)
}


# ── The error functions MetaFlux offers ───────────────────────────
# The registry (60_error_functions.R) lists every function users can choose, with a
# label and a fit function; the learning step learned each model with one of them. The
# report refits all of them, in their order in the registry, so a function added there
# appears here without a change to this script. The binned entry calls the vendored copy
# of DADA2's makeBinnedQualErrfun, which is therefore sourced first. Sourcing only
# defines functions; it reads no data and draws no random numbers.
source(snakemake@input[["vendored_binned"]])
source(snakemake@input[["error_functions"]])


# ── The user's own error function (EXPERIMENTAL) ──────────────────
# When amplicon.dada2.learn_errors.custom_error_function names an R file and a function
# (checked at start-up in 00_common.smk), the rule hands the file over as an input and
# the name and the path as written in the config as params; otherwise the input is
# empty (NULL here) and the params are "". The function is read with the registry's
# load_custom_error_function() and added to the registry as the entry "custom", so that
# everything below that goes through the registry (the refits, the tables, the figures)
# takes it in like loess and binned, whatever error_model says: under error_model:
# custom it is the function the models were learned with, and its refit must give back
# the model like any other; under any other setting it is one more comparison. It is
# done here, before the constants below, because some of them are built from the
# registry's names (the rate_<function> columns).
# A file that cannot be read now does not fail the report (a failed report would block
# the denoising): the entry then stops with the reason whenever it is refitted, so the
# function has no curve, and the reason becomes a warning of the report (added below,
# once the warnings exist).
custom_function_file_input         <- snakemake@input[["custom_function_file"]]
custom_function_name               <- snakemake@params[["custom_function_name"]]
custom_function_file_as_configured <- snakemake@params[["custom_function_file_as_configured"]]
custom_function      <- NULL     # what load_custom_error_function() returned, when it worked
custom_load_problem  <- NULL     # why it did not, as a sentence
if (length(custom_function_file_input) > 0) {
  loaded <- tryCatch(load_custom_error_function(custom_function_file_input[[1]], custom_function_name,
                                                custom_function_file_as_configured),
                     error = function(e) e)
  if (inherits(loaded, "error")) {
    custom_load_problem <- conditionMessage(loaded)
    error_function_registry[["custom"]] <- list(
      label = custom_function_label(custom_function_name),
      suits = NA_character_,
      fit   = function(trans, quality_values, extend_range = FALSE) stop(custom_load_problem, call. = FALSE))
  } else {
    custom_function <- loaded
    error_function_registry[["custom"]] <- custom_registry_entry(custom_function)
  }
}


# ── Constants ─────────────────────────────────────────────────────

# The 16 rows of DADA2's tables, in its order: a true A read as A, C, G or T, then a
# true C, and so on. The four rows of bases read correctly carry no error rate of their
# own (1 minus the three misreads of that nucleotide); the other 12 are the kinds of
# misread ("substitutions") the figures and tables show.
TRANSITION_ROWS <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
CORRECT_ROWS    <- c("A2A", "C2C", "G2G", "T2T")
SUBSTITUTIONS   <- setdiff(TRANSITION_ROWS, CORRECT_ROWS)

# A fitted rate "rises with quality" where it is more than 10 % above the rate at the
# next lower quality score (rising_steps below). A model's rates should fall as the
# quality score rises; small rises are noise of the fit, hence the 10 % margin.
RISE_FACTOR <- 1.10

# Colours. The error functions take blue, orange and aqua in the order of the registry
# (a palette checked for colour-blind readers, all three pairs); everything that is not
# a function is drawn in shades of ink, so a colour always means a function. Aqua has
# low contrast on white, which is why the figures always carry a legend and give every
# function its own line type as well. A 4th function would get no new colour but the
# secondary ink and a line type of its own.
FUNCTION_COLOURS <- c("#2a78d6", "#eb6834", "#1baf7a")
INK_PRIMARY      <- "#0b0b0b"     # text, the nominal line, the change per round
INK_SECONDARY    <- "#52514e"     # counted rates, the distance to the final model
INK_MUTED        <- "#898781"     # cells without errors, values at or below Q2
BAND_LIGHT       <- "#e1e0d9"     # bands that mark DADA2's repeat or the flat window
GRID_LIGHT       <- "#ecebe6"     # gridlines, one step lighter than the bands

# Line types and widths of the curves. The model used is solid and thicker; a refit is
# dashed and thinner; from a 4th function on, each refit gets another line type.
MODEL_LINE_WIDTH   <- 1.0
REFIT_LINE_WIDTH   <- 0.6
NOMINAL_LINE_WIDTH <- 0.5
REFIT_LINE_TYPES   <- c("dashed", "dashed", "dashed", "dotdash", "longdash", "twodash")

# The legend texts of the points and of the nominal line (figures are plain ASCII: the
# PDF device's standard fonts have no other characters).
POINT_READ_VALUE <- "counted rate at a quality value the reads carry"
POINT_AVERAGED   <- "counted rate at an averaged quality of identical reads"
POINT_NO_ERROR   <- "no error counted"
NOMINAL_LABEL    <- "nominal rate of the quality score"

# Statement S: what the comparison of the error functions shows and what it cannot. The
# same text is the description of the MultiQC section (parent_description in
# workflow/multiqc_config.yaml) and is quoted in the documentation;
# tests/test_error_model_settings.py keeps the copies equal (it reads the text between
# the quotes below, with every run of spaces and line ends made one space, as the
# script does).
REPORT_STATEMENT_TEXT <- "
  The comparison refits only the smoothing step: each error function MetaFlux offers is
  applied to the transition counts that the model used was fitted to (for a model ended
  by the flat-rounds rule, those of the round before the last). It shows how closely
  each function describes those counts, not which function would have given the better
  ASVs: learning with another function would have produced other counts, and a more
  flexible function always comes closer to the counts it is fitted to. The numbers
  compare functions within this run only; they cannot judge a run, nor tell faulty reads
  from good ones."
REPORT_STATEMENT <- trimws(gsub("\\s+", " ", REPORT_STATEMENT_TEXT))

# The descriptions of the two MultiQC tables, written into the headers of their files
# (the four figures get theirs from workflow/multiqc_config.yaml). MultiQC reads these
# headers as YAML, in double quotes, so the texts hold no double quote or backslash.
# The ranges quoted in GOODNESS_OF_FIT_DESCRIPTION ("on the models MetaFlux was tested
# on") are those of the table "Reference values" in docs/amplicon/error-model.md: the
# MiSeq models of the eight test sets (loess) and the NextSeq 2000 16S models (one
# 16S run with binned values and the Zymo D6331 mock). A model added to that table that
# falls outside them means changing both.
SUMMARY_DESCRIPTION <- paste(
  "One row per read direction: the error function and how it was chosen, how the",
  "learning ended, and the model code. Converged no means the model did not settle. If",
  "stats/dada2/error_model_check.txt exists and says accepted, the run went on with this",
  "model; if that file is missing, the run stopped before denoising (MultiQC does not",
  "wait for that check).")
GOODNESS_OF_FIT_DESCRIPTION <- paste(
  "Each error function MetaFlux offers, refitted on the counts behind the model used",
  "(see the module description). Rising steps: places where the curve gives a higher",
  "quality score more errors (a rise of more than 10 % between neighbouring scores), out",
  "of those compared. Typical difference: how many times the curve is above or below the",
  "counted rates, on average weighted by bases (1.00x = on the counted rates). Averaged",
  "columns are consensus qualities of identical reads that no read carries itself; they",
  "exist only with binned quality values, and the notes in stats/dada2/error_model_fit.tsv",
  "give the share of the bases and of the misreads they hold. On the models MetaFlux was",
  "tested on, loess on continuous values (MiSeq) rose in 2 to 40 of 348 to 372 steps and",
  "sat 1.11x to 1.78x from the counted rates; loess on binned values (NextSeq 2000) rose in",
  "52 to 77 of 300 steps; the binned fit on binned values never rose and sits on the",
  "counted rates at the reads' values by construction. The function used is not always",
  "the closest: each function wins on what it is built for, and the choice follows the",
  "kind of quality values (DADA2's author recommends the binned fit for binned values) and",
  "whether the rates fall as quality rises, not closeness to these counts. The full",
  "numbers are in stats/dada2/error_model_fit.tsv.")

# What each column of the two MultiQC tables shows, for the tooltip of its title.
SUMMARY_COLUMN_DESCRIPTIONS <- c(
  "Error function"                 = "The function the model was learned with; a binned fit names the quality values it is anchored on.",
  "Chosen by"                      = "auto: chosen from the quality values of the filtered reads; setting: amplicon.dada2.learn_errors.error_model.",
  "Quality values above Q2"        = "How many distinct quality values above Q2 the filtered reads of the voting samples carry, and which.",
  "Rounds"                         = "Rounds of learning.",
  "Cap"                            = "amplicon.dada2.learn_errors.max_consist: the most rounds allowed.",
  "How learning ended"             = "DADA2's own stop (a round gave back the table of an earlier round: a repeat or a cycle), the flat-rounds rule, or the cap.",
  "Converged"                      = "yes only for DADA2's own stop.",
  "Change in rounds 1, 2 and last" = "How much the model changed in round 1, round 2 and the last round, summed over all its rates. Close to convergence the last is several orders of magnitude below the first two.",
  "Bases learned from"             = "What learnErrors read: whole samples in random order until amplicon.dada2.learn_errors.nbases bases were in.",
  "Model code"                     = "The code that accepts this model in amplicon.dada2.learn_errors.accept_unconverged when it did not converge.",
  "Warnings"                       = "Warnings of the learning and of this report.",
  "Custom error function"          = "Experimental, only in a run that names one: your own error function (amplicon.dada2.learn_errors.custom_error_function), the first 8 characters of the md5 sum of its file (they change with any edit), and whether this model was learned with it (error_model: custom) or it was only refitted on the same counts for comparison.")
GOODNESS_OF_FIT_COLUMN_DESCRIPTIONS <- c(
  "Model used"                                      = "used: the function the model was learned with (its refit gives back the model); refit only: a function refitted on the same counts for comparison.",
  "Anchors"                                         = "The quality values a binned fit is anchored on.",
  "Rising steps"                                    = "Neighbouring quality scores where the curve rises by more than 10 %, out of those compared (12 kinds of misread).",
  "Substitution types with a rise"                  = "Of the 12 kinds of misread, how many have at least one rising step.",
  "Typical difference at the reads' quality values" = "How many times the curve typically lies above or below the counted rates at the quality values the reads carry, weighted by bases (1.00x = on them).",
  "Typical difference at averaged columns"          = "The same at the averaged qualities of identical reads (binned data only).",
  "Note"                                            = "Where the numbers of this row mislead.")
# Shorter titles for the long column names of the goodness-of-fit table, so that the
# Note column keeps some width on the page; the descriptions above say the rest.
GOODNESS_OF_FIT_COLUMN_TITLES <- c(
  "Substitution types with a rise"                  = "Types with a rise",
  "Typical difference at the reads' quality values" = "Typical difference: read values",
  "Typical difference at averaged columns"          = "Typical difference: averaged columns")

# The columns of error_model_fit.tsv (and of the fit objects in the report JSON), in
# order. What each holds is described at fit_rows_of_direction() and goodness_of_fit().
FIT_TABLE_COLUMNS <- c(
  "read_direction", "error_function", "model_used", "anchors", "refit_status",
  "refit_reproduces_model", "bases", "share_bases_at_read_values", "loglik_per_base",
  "deviance_per_million_bases", "deviance_per_million_bases_at_read_values",
  "deviance_per_million_bases_at_averaged_columns", "typical_fold_difference",
  "typical_fold_difference_at_read_values", "typical_fold_difference_at_averaged_columns",
  "cells_with_errors", "rising_steps", "steps_compared", "substitutions_with_rise", "note")
# The rows built for those files carry one more column, the short note of the MultiQC
# table (fit_note()).
FIT_ROW_COLUMNS <- c(FIT_TABLE_COLUMNS, "note_short")
RATES_TABLE_COLUMNS <- c(
  "read_direction", "substitution", "quality", "column_kind", "bases", "errors",
  "observed_rate", "nominal_rate", "model_rate",
  paste0("rate_", names(error_function_registry)))
ROUNDS_TABLE_COLUMNS <- c(
  "read_direction", "round", "change", "lowest_change_so_far", "distance_to_final",
  "depth_limit", "in_flat_window", "event")


# ── Warnings of this report ───────────────────────────────────────
# Everything that goes wrong without stopping the job (a figure that could not be
# drawn, counts that differ from DADA2's, a refit of the function used that does not
# give back the model) is collected here, with the read direction it concerns ("both"
# for a figure of both directions). The texts go into error_model_report.json, which
# check_error_model reads to print a note, and into the Warnings column of the MultiQC
# summary; each is also written to the log when it happens.
report_state <- new.env()
report_state$warning_texts      <- character(0)
report_state$warning_directions <- character(0)

#' Record a warning of this report.
#'
#' @param text What went wrong, in one sentence.
#' @param direction "R1", "R2" or "both".
#' @return Nothing.
add_warning <- function(text, direction = "both") {
  full_text <- if (direction == "both") text else paste0(direction, ": ", text)
  report_state$warning_texts      <- c(report_state$warning_texts, full_text)
  report_state$warning_directions <- c(report_state$warning_directions, direction)
  message("[error_model_report] WARNING: ", full_text)
}


# ── Small helpers ─────────────────────────────────────────────────

#' A value read from the JSON as a plain vector: NULL (a JSON null) and list() (an
#' empty JSON list) become a vector of length 0.
int_values  <- function(x) if (is.null(x)) integer(0) else as.integer(unlist(x))
num_values  <- function(x) if (is.null(x)) numeric(0) else as.numeric(unlist(x))
text_values <- function(x) if (is.null(x)) character(0) else as.character(unlist(x))

#' A number with three significant digits, for texts (27.8, 0.0824, 1.13e-05).
format_number <- function(x) trimws(formatC(x, digits = 3, format = "g"))

#' A count with thousands separators (52,975,572).
format_count <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

#' A fold difference as MultiQC shows it ("1.15x"), or "-" when there is none.
format_fold <- function(x) if (is.na(x)) "-" else sprintf("%.2fx", x)

#' "1 round" or "200 rounds".
rounds_text <- function(count) if (count == 1) "1 round" else paste(count, "rounds")

#' How the figures, the MultiQC tables and the log name an error function: "loess fit"
#' and "binned fit" for MetaFlux's own, "custom (<name>)" for the user's own function
#' (the name as set in amplicon.dada2.learn_errors.custom_error_function).
#'
#' @param function_name The registry name ("loess", "binned" or "custom").
#' @return The text.
function_display_name <- function(function_name) {
  if (identical(function_name, "custom")) {
    return(sprintf("custom (%s)", custom_function_name))
  }
  paste(function_name, "fit")
}

#' Quality values as a short text: "9, 20, 34", or "7-38" when there are more than ten.
values_text <- function(values) {
  values <- sort(unique(values))
  if (length(values) == 0) return("none")
  if (length(values) > 10) return(sprintf("%d-%d", min(values), max(values)))
  paste(values, collapse = ", ")
}

#' The anchors of a binned fit for a figure title or a table cell: "9, 20, 34", or "16
#' values (7-38)" when there are more than ten.
#'
#' Why the short form: ggplot does not wrap a title, it cuts it at the right edge of the
#' picture, after about 133 characters in the fit figure and 138 in DADA2's plot (bold
#' and plain titles at the sizes drawn here). A list of more than ten anchors, which
#' error_model: binned set by hand on 11 to 19 quality values can give, would push the end
#' of the title, the outcome and the model code, off the picture. The exact anchors stay
#' in the record (anchors_used), in the learning log and in error_model_fit.tsv.
#'
#' @param anchors Integer vector of the anchors.
#' @return The text.
anchor_list_text <- function(anchors) {
  anchors <- sort(unique(anchors))
  if (length(anchors) > 10) return(sprintf("%d values (%s)", length(anchors), values_text(anchors)))
  values_text(anchors)
}

#' Text made safe for a figure: plain ASCII (R's own error messages can hold curly
#' quotes, which the PDF device's standard fonts cannot draw).
ascii_text <- function(text) iconv(text, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "?")

#' Text made safe for one cell of a MultiQC table: no tab and no line end, and ", "
#' written as "; ". MultiQC 1.35 shows every comma in a table cell as a thin space (its
#' way of showing thousands separators, which is right for 1,234,567 but turns
#' "9, 20, 34" into "9 20 34" and takes the commas out of a sentence).
cell_text <- function(text) {
  text <- gsub("[\t\r\n]+", " ", text)
  gsub(", ", "; ", text, fixed = TRUE)
}

#' Fill in the output folder and the read direction in a sentence of the registry,
#' which knows neither and writes the placeholders {out} and {dir}.
fill_placeholders <- function(text, direction) {
  text <- gsub("{out}", output_folder, text, fixed = TRUE)
  gsub("{dir}", direction, text, fixed = TRUE)
}

#' A table with the given columns and no rows, for a table that could not be built.
empty_table <- function(column_names) {
  table <- data.frame(matrix(nrow = 0, ncol = length(column_names)))
  names(table) <- column_names
  table
}

#' Write a table as a tab-separated file with a header line and no row names.
write_tsv <- function(table, path) {
  write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}


# ── Titles and texts about how the learning ended ─────────────────

#' The title of one model's figures, from its part of the record.
#'
#' Examples:
#'   "R1 error model: loess fit - converged (round 4 repeated its starting table)"
#'   "R2 error model: loess fit - converged (cycle of 11 rounds, rounds 22-32)"
#'   "R2 error model: loess fit - NOT converged: flat-rounds rule, round 42 (code R2:7ebe3bed)"
#'   "R1 error model: loess fit - NOT converged: cap of 200 rounds (code R1:26dab1e9)"
#'   "R2 error model: binned fit on 9, 20, 34 - converged (cycle of 3 rounds, rounds 5-7)"
#'   "R1 error model: binned fit on 16 values (7-38) - NOT converged: cap of 3 rounds (code R1:4be06c1d)"
#'   "R1 error model: custom (gam_errfun) - converged (round 5 repeated its starting table)"
#' A fit anchored on quality values (the binned fit) names them, because they are
#' where its straight lines bend; more than ten are shortened to their number and range
#' (anchor_list_text), so that the outcome and the model code at the end of the title
#' stay on the picture. Used for the fit figures and DADA2's own plots.
#'
#' @param direction "R1" or "R2".
#' @param part The record's entry of this direction (a list read from the JSON).
#' @param max_consist The cap on rounds (top level of the record).
#' @return The title text.
plot_title <- function(direction, part, max_consist) {
  function_text <- function_display_name(part$error_function)
  # anchors_used is null in the record for a fit without anchors (loess).
  anchors <- int_values(part$anchors_used)
  if (length(anchors) > 0 && !all(is.na(anchors))) {
    function_text <- paste0(function_text, " on ", anchor_list_text(anchors))
  }
  if (part$stopped_by == "dada2_repeat") {
    if (part$cycle_length == 1) {
      outcome <- sprintf("converged (round %d repeated its starting table)", part$rounds)
    } else {
      outcome <- sprintf("converged (cycle of %d rounds, rounds %d-%d)",
                         part$cycle_length, part$matched_round, part$rounds)
    }
  } else if (part$stopped_by == "flat_rounds") {
    outcome <- sprintf("NOT converged: flat-rounds rule, round %d (code %s)",
                       part$flat_stop_round, part$model_code)
  } else {
    outcome <- sprintf("NOT converged: cap of %s (code %s)", rounds_text(max_consist), part$model_code)
  }
  sprintf("%s error model: %s - %s", direction, function_text, outcome)
}


#' How the learning of one direction ended, in a sentence for the report JSON.
#'
#' It describes the learning only. Whether the run went on with a model that did not
#' converge is decided later, by check_error_model, which this report runs before, so
#' the text never says "accepted".
#'
#' @param part The record's entry of the direction.
#' @param max_consist The cap on rounds.
#' @param max_flat_rounds N of the flat-rounds rule.
#' @return The sentence.
learning_outcome_text <- function(part, max_consist, max_flat_rounds) {
  if (part$stopped_by == "dada2_repeat") {
    if (part$cycle_length == 1) {
      return(sprintf("converged: round %d gave back the table that went into it", part$rounds))
    }
    return(sprintf("converged: round %d gave back the table that went into round %d (a cycle of %d rounds)",
                   part$rounds, part$matched_round, part$cycle_length))
  }
  if (part$stopped_by == "flat_rounds") {
    return(sprintf("not converged: ended by the flat-rounds rule in round %d (max_flat_rounds %d)",
                   part$flat_stop_round, max_flat_rounds))
  }
  sprintf("not converged: stopped at the cap of %s (max_consist)", rounds_text(max_consist))
}


#' How the learning of one direction ended, in a few words for the MultiQC summary.
#'
#' @param part The record's entry of the direction.
#' @param max_consist The cap on rounds.
#' @return Text such as "round 5 repeated its starting table".
short_outcome_text <- function(part, max_consist) {
  if (part$stopped_by == "dada2_repeat") {
    if (part$cycle_length == 1) {
      return(sprintf("round %d repeated its starting table", part$rounds))
    }
    return(sprintf("cycle of %d rounds (rounds %d-%d)", part$cycle_length, part$matched_round, part$rounds))
  }
  if (part$stopped_by == "flat_rounds") {
    return(sprintf("flat-rounds rule in round %d", part$flat_stop_round))
  }
  sprintf("cap of %s", rounds_text(max_consist))
}


# ── Step 1: the counts behind each model ──────────────────────────

#' The transition counts the final table of one direction was fitted to, rebuilt from
#' the record as the matrix DADA2 hands an error function.
#'
#' learn_errors writes them as counts_behind_model: the quality scores (0, 1, ..., one
#' per column) and one list of counts per row, A2A to T2T. They are DADA2's counts of the
#' last round (err$trans), or, for a model ended by the flat-rounds rule, those of the
#' round before, because the kept table was fitted to them; DADA2 keeps no copy of those.
#'
#' @param part The record's entry of the direction.
#' @return A 16 x Q matrix with row names A2A ... T2T and column names "0", "1", ...
#'   (the counts as the JSON gave them: whole numbers, as DADA2 1.34.0 counts).
counts_behind_model <- function(part) {
  behind <- part$counts_behind_model
  if (is.null(behind) || is.null(behind$rows) || !identical(names(behind$rows), TRANSITION_ROWS)) {
    stop("the record holds no counts_behind_model with the 16 rows A2A to T2T")
  }
  counts <- do.call(rbind, behind$rows)
  quality_scores <- int_values(behind$quality_scores)
  if (ncol(counts) != length(quality_scores)) {
    stop("counts_behind_model has a different number of counts and quality scores")
  }
  dimnames(counts) <- list(TRANSITION_ROWS, as.character(quality_scores))
  counts
}


#' Whether the counts rebuilt from the record equal DADA2's own counts.
#'
#' Compared as numbers, with their row and column names: dada2 1.34.0 stores its counts
#' as whole numbers (integers) and dada2 1.38.0 as decimal numbers (doubles), which are
#' the same counts but not identical() objects, so both are made doubles first.
#'
#' @param counts From counts_behind_model().
#' @param trans DADA2's counts of the last round (err$trans).
#' @return TRUE or FALSE.
same_counts <- function(counts, trans) {
  counts_as_numbers <- counts
  storage.mode(counts_as_numbers) <- "double"
  trans_as_numbers <- trans
  storage.mode(trans_as_numbers) <- "double"
  identical(counts_as_numbers, trans_as_numbers)
}


#' What kind of quality value each column of the counts is.
#'
#' "read_value": a quality value the reads carry (the voting samples' values in the
#' record); "averaged": a column that holds bases at a value no read carries, a
#' consensus quality of identical reads (binned data only); "no_data": no bases.
#'
#' @param counts The counts (columns "0", "1", ...).
#' @param read_values The quality values the reads carry.
#' @return Character vector, one entry per column.
column_kinds <- function(counts, read_values) {
  quality <- as.integer(colnames(counts))
  has_bases <- colSums(counts) > 0
  kinds <- rep("no_data", length(quality))
  kinds[has_bases & quality %in% read_values] <- "read_value"
  kinds[has_bases & !(quality %in% read_values)] <- "averaged"
  kinds
}


# ── Step 2: refit every error function on those counts ────────────

#' Why a refitted table cannot be compared, or NULL when it can.
#'
#' @param table What a fit returned.
#' @param counts The counts it was fitted to.
#' @return NULL, or a few words about the problem.
table_problem <- function(table, counts) {
  if (!is.matrix(table) || !is.numeric(table) || !identical(dim(table), dim(counts))) {
    return("its table does not have the shape of the counts")
  }
  if (!identical(dimnames(table), dimnames(counts))) {
    return("its table has other row or column names than the counts")
  }
  if (!all(is.finite(table)) || any(table < 0) || any(table > 1)) {
    return("its table holds a missing, infinite or out-of-range rate")
  }
  NULL
}


#' Refit every function of the registry once on the counts behind a model.
#'
#' Only the smoothing step is repeated: each function turns the same counts into
#' rates, as it would in a round of learning. The function the model was learned with
#' is refitted exactly as in learning (extend_range = FALSE), so its table must equal
#' the model. The other functions are refitted with extend_range = TRUE: when the counts
#' hold bases outside the quality values a binned fit is anchored on, the end columns
#' become anchors and the note says so, instead of no curve at all. Warnings of a fit
#' are not shown (DADA2 hides them in learning too). A fit that fails gives a status
#' "failed: <reason>" and no table; it never fails the job.
#'
#' @param counts The counts behind the model.
#' @param read_values The quality values of the voting samples (the binned fit's anchors).
#' @param function_used The registry name of the function the model was learned with.
#' @param direction "R1" or "R2", to fill in the placeholders of a reason.
#' @return A list with one entry per registry name, each list(status = "ok" or
#'   "failed: ...", table = the fitted table or NULL, anchors, note).
refit_every_function <- function(counts, read_values, function_used, direction) {
  refits <- list()
  for (function_name in names(error_function_registry)) {
    fit_function <- error_function_registry[[function_name]]$fit
    comparison_only <- function_name != function_used
    result <- tryCatch(
      suppressWarnings(fit_function(counts, read_values, extend_range = comparison_only)),
      error = function(e) e)
    if (inherits(result, "error")) {
      refits[[function_name]] <- list(
        status = paste("failed:", fill_placeholders(conditionMessage(result), direction)),
        table = NULL, anchors = NA_integer_, note = "")
      next
    }
    problem <- table_problem(result$table, counts)
    if (!is.null(problem)) {
      refits[[function_name]] <- list(status = paste("failed:", problem), table = NULL,
                                      anchors = NA_integer_, note = "")
      next
    }
    refits[[function_name]] <- list(status = "ok", table = result$table,
                                    anchors = result$anchors,
                                    note = fill_placeholders(result$note, direction))
  }
  refits
}


# ── Step 3: how closely each fitted table follows the counts ──────

#' How often the fitted rates rise with quality.
#'
#' Over the 12 kinds of misread, the number of neighbouring quality scores q-1, q
#' between the lowest and the highest column that holds bases where the fitted rate at
#' q is more than 10 % above the rate at q-1. A model's rates should fall as quality
#' rises, so a rise is a place where the curve gives better bases more errors.
#'
#' @param fitted A fitted table.
#' @param counts The counts (to find the columns that hold bases).
#' @return list(rises, steps = 12 x (highest - lowest column), rows_with_rise = how many
#'   of the 12 kinds of misread have at least one rise).
rising_steps <- function(fitted, counts) {
  columns_with_bases <- which(colSums(counts) > 0)
  lowest <- min(columns_with_bases)
  highest <- max(columns_with_bases)
  rises <- 0
  rows_with_rise <- 0
  for (substitution in SUBSTITUTIONS) {
    rates <- fitted[substitution, lowest:highest]
    rises_here <- sum(rates[-1] > rates[-length(rates)] * RISE_FACTOR)
    rises <- rises + rises_here
    if (rises_here > 0) rows_with_rise <- rows_with_rise + 1
  }
  list(rises = rises, steps = 12 * (highest - lowest), rows_with_rise = rows_with_rise)
}


#' How closely one fitted table follows the counts it was fitted to.
#'
#' Notation: N is a count (bases of true nucleotide i read as j at quality q), E the
#' fitted rate of the same cell, and n the bases of true nucleotide i at quality q (the
#' sum of its four counts), so N/n is the counted rate.
#'   - loglik_per_base: the sum of N x log(E) over all cells with N > 0 (the 16 rows),
#'     divided by the bases: the mean log-probability of a counted base under the fitted
#'     rates. Higher (closer to 0) is closer.
#'   - deviance: 2 x (the sum of N x log(N/n) minus the sum of N x log(E)), over the same
#'     cells: how much less probable the counts are under the fitted rates than under
#'     their own counted rates (a multinomial deviance, computed per column). 0 means
#'     the fitted rates equal the counted ones. It is given per million bases (of all
#'     the bases), in total and split into the read values and the averaged columns;
#'     the two parts add up to the total.
#'   - typical fold difference: over the misread cells with at least one error, 10 to
#'     the power of the mean of |log10(N/n) - log10(E)|, weighted by n: how many times
#'     the curve typically lies above or below the counted rates (1.00 = on them). In
#'     total and per kind of column; missing when that kind has no such cell.
#'   - cells_with_errors: misread cells with at least one error, over all columns.
#'   - rising steps: see rising_steps().
#' No measure ranks the functions (see the header).
#'
#' @param fitted A fitted table (16 x Q).
#' @param counts The counts it was fitted to.
#' @param kinds column_kinds() of the counts.
#' @return A named list of the measures.
goodness_of_fit <- function(fitted, counts, kinds) {
  # As decimal numbers: R gives up (NA) on a sum of whole numbers above about 2.1
  # billion, and the bases of a run can come close to that.
  storage.mode(counts) <- "double"
  column_count <- ncol(counts)
  deviance_by_column <- numeric(column_count)
  difference_by_column <- numeric(column_count)   # sum of n x |log10 difference|
  weight_by_column <- numeric(column_count)       # sum of n over the same cells
  loglik_sum <- 0
  cells_with_errors <- 0
  for (nucleotide in c("A", "C", "G", "T")) {
    rows_of_nucleotide <- paste0(nucleotide, "2", c("A", "C", "G", "T"))
    bases_of_nucleotide <- colSums(counts[rows_of_nucleotide, , drop = FALSE])
    for (row in rows_of_nucleotide) {
      count <- counts[row, ]
      counted <- count > 0
      log_fitted <- log(fitted[row, counted])
      loglik_sum <- loglik_sum + sum(count[counted] * log_fitted)
      deviance_by_column[counted] <- deviance_by_column[counted] +
        2 * count[counted] * (log(count[counted] / bases_of_nucleotide[counted]) - log_fitted)
      if (row %in% SUBSTITUTIONS) {
        log10_difference <- abs(log10(count[counted] / bases_of_nucleotide[counted]) -
                                log10(fitted[row, counted]))
        difference_by_column[counted] <- difference_by_column[counted] +
          bases_of_nucleotide[counted] * log10_difference
        weight_by_column[counted] <- weight_by_column[counted] + bases_of_nucleotide[counted]
        cells_with_errors <- cells_with_errors + sum(counted)
      }
    }
  }
  total_bases <- sum(counts)
  every_column <- rep(TRUE, column_count)
  at_read_values <- kinds == "read_value"
  at_averaged_columns <- kinds == "averaged"
  deviance_per_million <- function(selected) sum(deviance_by_column[selected]) / total_bases * 1e6
  typical_fold <- function(selected) {
    weight <- sum(weight_by_column[selected])
    if (weight == 0) return(NA_real_)
    10^(sum(difference_by_column[selected]) / weight)
  }
  rise <- rising_steps(fitted, counts)
  list(
    bases                                          = total_bases,
    share_bases_at_read_values                     = sum(counts[, at_read_values]) / total_bases,
    loglik_per_base                                = loglik_sum / total_bases,
    deviance_per_million_bases                     = deviance_per_million(every_column),
    deviance_per_million_bases_at_read_values      = deviance_per_million(at_read_values),
    deviance_per_million_bases_at_averaged_columns = deviance_per_million(at_averaged_columns),
    typical_fold_difference                        = typical_fold(every_column),
    typical_fold_difference_at_read_values         = typical_fold(at_read_values),
    typical_fold_difference_at_averaged_columns    = typical_fold(at_averaged_columns),
    cells_with_errors                              = cells_with_errors,
    rising_steps                                   = rise$rises,
    steps_compared                                 = rise$steps,
    substitutions_with_rise                        = rise$rows_with_rise)
}


#' What the averaged columns of one direction hold, in words, for the notes of the fit
#' table: "7.6 % of the bases and 20.8 % of the counted misreads".
#'
#' Why the shares and not a fixed sentence: how many of the misreads the averaged
#' columns hold depends on the data. On the NextSeq 2000 models MetaFlux was first
#' checked on, they held 8 % and 32 % of the bases but only 1.4 % and 2.8 % of the
#' misreads; one sample of the same run, learned alone, gave 7.6 % and 35 % of the bases
#' and 21 % and 28 % of the misreads. So the note gives the numbers of this run.
#'
#' @param analysis analyse_direction() of the direction (its two shares).
#' @return The text, or "" when the counts have no averaged column or the shares are
#'   missing.
averaged_columns_text <- function(analysis) {
  share_at_read_values <- analysis$share_bases_at_read_values
  share_of_misreads <- analysis$share_substitutions_at_averaged_columns
  # Both shares are single numbers, NA when they could not be computed.
  if (length(share_at_read_values) != 1 || is.na(share_at_read_values)) return("")
  share_of_bases <- 1 - share_at_read_values      # the columns without bases hold none
  if (share_of_bases <= 0) return("")
  bases_part <- sprintf("%.1f %% of the bases", 100 * share_of_bases)
  if (length(share_of_misreads) != 1 || is.na(share_of_misreads)) return(bases_part)
  sprintf("%s and %.1f %% of the counted misreads", bases_part, 100 * share_of_misreads)
}


#' The note of the fit table where the numbers of a function mislead, in two lengths.
#'
#' Each function wins on what it is built for, so a number read on its own can point the
#' wrong way: the binned fit sits on the counted rates at the reads' values by
#' construction (its deviance there is about 0, wherever it is used), and loess on
#' binned values also passes near the averaged columns, which can lower its total
#' deviance below the binned fit's. The full sentence goes into error_model_fit.tsv and
#' the report JSON, with the shares of the bases and the misreads the averaged columns
#' hold in this run (averaged_columns_text()); the short form, a few words, into the
#' MultiQC table, where a long text would be squeezed into a narrow column (the table's
#' description says the rest).
#'
#' @param function_name The registry name of the function.
#' @param quality_class The kind of quality values of the reads ("binned",
#'   "continuous" or "undetermined"; from the record).
#' @param averaged_text What the averaged columns hold (averaged_columns_text()), or ""
#'   when the counts have none.
#' @param refit_note The fit's own note (extra anchors of a comparison refit), or "".
#' @return list(full = the sentence(s) joined with "; ", short = the short form); both
#'   "" when there is nothing to say. The user's own function (custom) always gets a
#'   note that it is experimental, with its name and the first 8 characters of its
#'   file's md5 sum, which tell two versions of the file apart.
fit_note <- function(function_name, quality_class, averaged_text, refit_note) {
  full <- character(0)
  short <- character(0)
  suits <- error_function_registry[[function_name]]$suits
  if (function_name == "custom") {
    md5_text <- if (is.null(custom_function)) "file not read" else paste("md5", substr(custom_function$md5, 1, 8))
    full <- c(full, sprintf("experimental: your own error function %s (%s, %s)",
                            custom_function_name, custom_function_file_as_configured, md5_text))
    short <- c(short, sprintf("experimental: your own function (%s)", md5_text))
  } else if (function_name == "binned" && identical(quality_class, "continuous")) {
    full <- c(full, paste("not the function for continuous quality values: the binned",
                          "fit joins the counted rates with straight lines"))
    short <- c(short, "made for binned values: joins the counted rates point to point")
  } else if (function_name == "binned" && identical(quality_class, "binned")) {
    sentence <- "sits on the counted rates at the reads' values by construction"
    if (nzchar(averaged_text)) {
      sentence <- paste0(sentence, "; the averaged columns hold ", averaged_text)
    }
    full <- c(full, sentence)
    short <- c(short, "on the counted rates at the reads' values by construction")
  } else if (function_name == "loess" && identical(quality_class, "binned")) {
    sentence <- "made for continuous quality values (DADA2's author recommends the binned fit for binned ones)"
    if (nzchar(averaged_text)) {
      sentence <- paste0(sentence, "; its total deviance also covers the averaged columns (",
                         averaged_text, ") and can come out lower than the binned fit's")
    }
    full <- c(full, sentence)
    short <- c(short, "made for continuous values")
  } else if (quality_class %in% c("binned", "continuous") && !identical(suits, quality_class)) {
    full <- c(full, sprintf("made for %s quality values and not for %s ones", suits, quality_class))
    short <- c(short, sprintf("made for %s values", suits))
  }
  if (nzchar(refit_note)) {
    full <- c(full, refit_note)
    short <- c(short, "end columns added as anchors")
  }
  list(full = paste(full, collapse = "; "), short = paste(short, collapse = "; "))
}


# ── Step 4a: the tables of one direction ──────────────────────────

#' The rows of error_model_fit.tsv for one direction: one per registered function.
#'
#' Columns (FIT_TABLE_COLUMNS): the direction; the function; model_used (TRUE for the
#' function the model was learned with); anchors (the quality values a binned fit is
#' anchored on, NA for loess); refit_status ("ok" or "failed: <reason>");
#' refit_reproduces_model (for the function used: whether its refit equals the model
#' exactly; NA for the others); the measures of goodness_of_fit(); note (fit_note()).
#' One more column, note_short, carries the short form of the note for the MultiQC
#' table; it is left out of error_model_fit.tsv and the report JSON.
#'
#' @param direction "R1" or "R2".
#' @param analysis What analyse_direction() has built so far (counts, kinds, refits).
#' @param part The record's entry of the direction.
#' @return A data frame with the columns FIT_ROW_COLUMNS, unrounded.
fit_rows_of_direction <- function(direction, analysis, part) {
  # The shares of the bases and of the misreads at the averaged columns, in words, for
  # the notes ("" when the counts have no averaged column).
  averaged_text <- if (any(analysis$kinds == "averaged")) averaged_columns_text(analysis) else ""
  rows <- list()
  for (function_name in names(error_function_registry)) {
    refit <- analysis$refits[[function_name]]
    is_used <- function_name == part$error_function
    measures <- NULL
    note <- list(full = "", short = "")
    if (is.null(refit$table)) {
      note$short <- "refit failed (see error_model_fit.tsv)"
    } else {
      measures <- tryCatch(goodness_of_fit(refit$table, analysis$counts, analysis$kinds),
                           error = function(e) e)
      if (inherits(measures, "error")) {
        note <- list(full = paste("the measures could not be computed:", conditionMessage(measures)),
                     short = "measures failed (see error_model_fit.tsv)")
        measures <- NULL
      } else {
        note <- fit_note(function_name, part$quality_class, averaged_text, refit$note)
      }
    }
    anchors <- int_values(refit$anchors)
    row <- data.frame(
      read_direction         = direction,
      error_function         = function_name,
      model_used             = is_used,
      anchors                = if (length(anchors) > 0 && !all(is.na(anchors))) paste(anchors, collapse = ", ") else NA_character_,
      refit_status           = refit$status,
      refit_reproduces_model = if (is_used) analysis$refit_reproduces_model else NA,
      stringsAsFactors = FALSE)
    measure_names <- setdiff(FIT_TABLE_COLUMNS, c(names(row), "note"))
    # A rate of 0 where an error was counted (no function MetaFlux ships gives one)
    # would make a measure infinite; it is written as missing, with the reason.
    has_infinite_measure <- FALSE
    for (measure in measure_names) {
      value <- if (is.null(measures)) NA_real_ else measures[[measure]]
      if (!is.na(value) && !is.finite(value)) {
        value <- NA_real_
        has_infinite_measure <- TRUE
      }
      row[[measure]] <- value
    }
    if (has_infinite_measure) {
      reason <- "the table gives a rate of 0 where errors were counted"
      note$full <- paste(c(note$full[nzchar(note$full)], reason), collapse = "; ")
      note$short <- paste(c(note$short[nzchar(note$short)], reason), collapse = "; ")
    }
    row$note <- note$full
    row$note_short <- note$short
    rows[[function_name]] <- row
  }
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  table[, FIT_ROW_COLUMNS]
}


#' The rows of error_model_rates.tsv for one direction: one per kind of misread and
#' quality score.
#'
#' Columns (RATES_TABLE_COLUMNS): the direction; the substitution (A2C: a true A read as
#' C); the quality score; column_kind (column_kinds()); bases (of the true nucleotide at
#' that score); errors (of this kind); observed_rate (errors / bases, NA without bases);
#' nominal_rate (one third of 10^(-Q/10), the rate of this one misread the quality
#' score promises, as DADA2's plotErrors draws it); model_rate (the model used, err_out);
#' rate_<function> for each registered function (its refit; NA when the refit failed).
#'
#' @param direction "R1" or "R2".
#' @param analysis What analyse_direction() has built so far.
#' @param err_out The model used (16 x Q).
#' @return A data frame, unrounded (the figure of this direction is drawn from it).
rates_of_direction <- function(direction, analysis, err_out) {
  counts <- analysis$counts
  storage.mode(counts) <- "double"
  quality <- as.integer(colnames(counts))
  pieces <- list()
  for (substitution in SUBSTITUTIONS) {
    true_nucleotide <- substr(substitution, 1, 1)
    bases <- colSums(counts[paste0(true_nucleotide, "2", c("A", "C", "G", "T")), , drop = FALSE])
    errors <- counts[substitution, ]
    piece <- data.frame(
      read_direction = direction,
      substitution   = substitution,
      quality        = quality,
      column_kind    = analysis$kinds,
      bases          = unname(bases),
      errors         = unname(errors),
      observed_rate  = ifelse(bases > 0, unname(errors / bases), NA_real_),
      nominal_rate   = (1 / 3) * 10^(-quality / 10),
      # The model's columns are matched by name, so a model with other columns than its
      # counts fails here (and becomes a warning) instead of being drawn shifted.
      model_rate     = unname(err_out[substitution, colnames(counts)]),
      stringsAsFactors = FALSE)
    for (function_name in names(error_function_registry)) {
      refit_table <- analysis$refits[[function_name]]$table
      piece[[paste0("rate_", function_name)]] <-
        if (is.null(refit_table)) NA_real_ else unname(refit_table[substitution, ])
    }
    pieces[[substitution]] <- piece
  }
  table <- do.call(rbind, pieces)
  rownames(table) <- NULL
  table
}


#' The rows of error_model_rounds.tsv for one direction: one per round of learning.
#'
#' Columns (ROUNDS_TABLE_COLUMNS): the direction; the round; change (how much the table
#' moved in that round, summed over all its rates); lowest_change_so_far; distance_to_final
#' (how far the table that went into the round is from the final model: DADA2's
#' checkConvergence); depth_limit (1/1000 of the change in round 1, the depth the
#' flat-rounds rule needs); in_flat_window (TRUE for the last max_flat_rounds rounds of
#' a model the flat-rounds rule ended); event ("" or, on its round: "matched_round", the
#' round whose starting table DADA2's last round gave back in a cycle; "dada2_repeat",
#' the last round of a converged model; "flat_stop"; "cap").
#'
#' @param direction "R1" or "R2".
#' @param part The record's entry of the direction.
#' @param max_flat_rounds N of the flat-rounds rule.
#' @return A data frame, unrounded (the convergence figure is drawn from it).
rounds_of_direction <- function(direction, part, max_flat_rounds) {
  rounds <- as.integer(part$rounds)
  change <- num_values(part$change_per_round)
  distance <- num_values(part$distance_to_final)
  if (length(change) != rounds || length(distance) != rounds) {
    stop(sprintf("the record lists %d rounds but %d changes and %d distances",
                 rounds, length(change), length(distance)))
  }
  event <- rep("", rounds)
  in_flat_window <- rep(FALSE, rounds)
  if (part$stopped_by == "dada2_repeat") {
    if (part$matched_round < rounds) event[part$matched_round] <- "matched_round"
    event[rounds] <- "dada2_repeat"
  } else if (part$stopped_by == "flat_rounds") {
    event[rounds] <- "flat_stop"
    in_flat_window[max(1, rounds - max_flat_rounds + 1):rounds] <- TRUE
  } else {
    event[rounds] <- "cap"
  }
  data.frame(read_direction = direction, round = seq_len(rounds), change = change,
             lowest_change_so_far = cummin(change), distance_to_final = distance,
             depth_limit = change[1] / 1000, in_flat_window = in_flat_window,
             event = event, stringsAsFactors = FALSE)
}


#' Everything the report says about one direction: counts, refits, measures and tables.
#'
#' Each part is made inside its own tryCatch: what cannot be made leaves an empty table
#' and a warning, never a failed job.
#'
#' @param direction "R1" or "R2".
#' @param part The record's entry of the direction.
#' @param err The model as learnErrors returned it.
#' @param max_flat_rounds N of the flat-rounds rule.
#' @return list(counts, kinds, read_values, refits, refit_reproduces_model, shares, and
#'   the three tables fit_rows, rates and rounds).
analyse_direction <- function(direction, part, err, max_flat_rounds) {
  analysis <- list(counts = NULL, kinds = NULL, read_values = integer(0), refits = list(),
                   refit_reproduces_model = NA,
                   share_bases_at_read_values = NA_real_,
                   share_substitutions_at_averaged_columns = NA_real_,
                   fit_rows = empty_table(FIT_ROW_COLUMNS),
                   rates = empty_table(RATES_TABLE_COLUMNS),
                   rounds = empty_table(ROUNDS_TABLE_COLUMNS))

  # The rounds table needs only the record, so it is made first and on its own.
  analysis$rounds <- tryCatch(
    rounds_of_direction(direction, part, max_flat_rounds),
    error = function(e) {
      add_warning(paste("the table of the rounds could not be made:", conditionMessage(e)), direction)
      empty_table(ROUNDS_TABLE_COLUMNS)
    })

  # Step 1: the counts. Without them nothing below can be computed.
  counts <- tryCatch(counts_behind_model(part), error = function(e) {
    add_warning(paste("the counts behind the model could not be read:", conditionMessage(e)), direction)
    NULL
  })
  if (is.null(counts)) {
    return(analysis)
  }
  analysis$counts <- counts
  # For a model DADA2's loop ended itself, or the cap, the record's counts are DADA2's
  # counts of the last round, and must equal them. After a flat-rounds stop they are the
  # counts of the round before, which DADA2 does not keep, so there is nothing to compare.
  if (part$stopped_by != "flat_rounds" && !same_counts(counts, err$trans)) {
    add_warning(paste("the counts in error_model.json differ from those in the model's .rds",
                      "file; the comparison below uses those of error_model.json"), direction)
  }

  # The quality values the reads carry: the voting samples' values, which
  # choose_error_model found and learn_errors copied into the record.
  analysis$read_values <- int_values(part$quality_values)
  if (length(analysis$read_values) == 0) {
    add_warning(paste("the record lists no quality values of the reads, so every column",
                      "with bases counts as a read value"), direction)
    analysis$read_values <- as.integer(colnames(counts))[colSums(counts) > 0]
  }
  analysis$kinds <- column_kinds(counts, analysis$read_values)

  # The shares the report JSON gives: of the bases at the read values, and of the counted
  # misreads at the averaged columns.
  counts_as_numbers <- counts
  storage.mode(counts_as_numbers) <- "double"
  bases_per_column <- colSums(counts_as_numbers)
  errors_per_column <- colSums(counts_as_numbers[SUBSTITUTIONS, , drop = FALSE])
  analysis$share_bases_at_read_values <-
    sum(bases_per_column[analysis$kinds == "read_value"]) / sum(bases_per_column)
  analysis$share_substitutions_at_averaged_columns <-
    if (sum(errors_per_column) > 0) sum(errors_per_column[analysis$kinds == "averaged"]) / sum(errors_per_column) else NA_real_

  # Step 2: the refits, and whether the function used gives back the model exactly.
  analysis$refits <- refit_every_function(counts, analysis$read_values, part$error_function, direction)
  used_table <- analysis$refits[[part$error_function]]$table
  analysis$refit_reproduces_model <- !is.null(used_table) && identical(used_table, err$err_out)
  if (!analysis$refit_reproduces_model) {
    add_warning(sprintf(paste("the refit of the %s fit on the counts behind the model does not",
                              "give back the model exactly (%s), so the comparison may not",
                              "describe the model used"),
                        part$error_function, analysis$refits[[part$error_function]]$status), direction)
  }

  # Steps 3 and 4: the measures and the two tables built from the counts.
  analysis$fit_rows <- tryCatch(fit_rows_of_direction(direction, analysis, part), error = function(e) {
    add_warning(paste("the goodness-of-fit table could not be made:", conditionMessage(e)), direction)
    empty_table(FIT_ROW_COLUMNS)
  })
  analysis$rates <- tryCatch(rates_of_direction(direction, analysis, err$err_out), error = function(e) {
    add_warning(paste("the table of the error rates could not be made:", conditionMessage(e)), direction)
    empty_table(RATES_TABLE_COLUMNS)
  })
  analysis
}


# ── Step 4b: the figures ──────────────────────────────────────────

#' The look shared by the report's figures: white background, hairline gridlines,
#' text in ink, legend at the bottom.
report_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID_LIGHT, linewidth = 0.3),
          axis.text        = element_text(colour = INK_SECONDARY),
          axis.title       = element_text(colour = INK_PRIMARY),
          strip.text       = element_text(colour = INK_PRIMARY, face = "bold", hjust = 0),
          plot.title       = element_text(colour = INK_PRIMARY, face = "bold"),
          plot.subtitle    = element_text(colour = INK_SECONDARY),
          plot.caption     = element_text(colour = INK_SECONDARY, hjust = 0),
          legend.position  = "bottom",
          legend.box       = "vertical",
          legend.title     = element_text(colour = INK_PRIMARY),
          legend.text      = element_text(colour = INK_PRIMARY))
}


#' A long text wrapped into lines of about `width` characters, for a subtitle or caption.
wrap_text <- function(text, width = 110) paste(strwrap(text, width = width), collapse = "\n")


#' The fit figure of one direction: counted error rates and the fitted curves.
#'
#' One panel per kind of misread (A2C ... T2G), quality score along the bottom, error
#' rate up the side on a log scale. Points: the rates counted in the reads learnErrors
#' drew, in grey ink, sized by the bases at that score; filled at the quality values the
#' reads carry, open at averaged qualities of identical reads. A cell with bases but no
#' error has a rate of 0, which a log scale cannot show; leaving it out would make the
#' visible points look higher than the counts are, so it is drawn as a small open
#' triangle at the bottom of the panel. Lines: the model used (solid, in its function's
#' colour), the other functions refitted on the same counts (dashed), and the nominal
#' rate of the quality score (dotted, in ink). Ticks along the bottom: the quality
#' values of the filtered reads (muted where the binned fit left a value out). The x axis
#' spans the columns that hold bases, one score more on each side; the y axis spans the
#' counted rates, the model and the nominal line, and a refit that runs outside it is
#' cut off at the panel's edge.
#'
#' @param direction "R1" or "R2".
#' @param analysis analyse_direction() of this direction.
#' @param part The record's entry of the direction.
#' @param max_consist The cap on rounds (for the title).
#' @return A ggplot object.
fit_figure <- function(direction, analysis, part, max_consist) {
  rates <- analysis$rates
  if (nrow(rates) == 0) {
    stop("there are no error rates to draw (see the warnings of this direction)")
  }
  function_used <- part$error_function
  function_names <- names(error_function_registry)

  # The x range: the columns that hold bases, one score more on each side.
  scores_with_bases <- rates$quality[rates$bases > 0]
  x_low  <- max(min(rates$quality), min(scores_with_bases) - 1)
  x_high <- min(max(rates$quality), max(scores_with_bases) + 1)
  in_range <- rates$quality >= x_low & rates$quality <= x_high

  # The points: cells with at least one error (filled or open), and cells without one.
  counted <- rates[rates$bases > 0 & rates$errors > 0, ]
  counted$point_kind <- ifelse(counted$column_kind == "read_value", POINT_READ_VALUE, POINT_AVERAGED)
  no_error <- rates[rates$bases > 0 & rates$errors == 0, ]

  # The curves, each a label with its colour, line type and width. The model used is
  # err_out itself; its refit equals it (refit_reproduces_model), so it is not drawn twice.
  curve_pieces <- list()
  curve_colours <- character(0)
  curve_types <- character(0)
  curve_widths <- numeric(0)
  model_label <- paste0(function_display_name(function_used), ": model used")
  curve_pieces[[model_label]] <- data.frame(substitution = rates$substitution, quality = rates$quality,
                                            rate = rates$model_rate, curve = model_label)[in_range, ]
  used_index <- match(function_used, function_names)
  curve_colours[model_label] <- if (!is.na(used_index) && used_index <= length(FUNCTION_COLOURS)) FUNCTION_COLOURS[used_index] else INK_SECONDARY
  curve_types[model_label] <- "solid"
  curve_widths[model_label] <- MODEL_LINE_WIDTH
  for (index in seq_along(function_names)) {
    function_name <- function_names[index]
    if (function_name == function_used || is.null(analysis$refits[[function_name]]$table)) next
    label <- paste0(function_display_name(function_name), ": refitted on the same counts")
    curve_pieces[[label]] <- data.frame(substitution = rates$substitution, quality = rates$quality,
                                        rate = rates[[paste0("rate_", function_name)]], curve = label)[in_range, ]
    curve_colours[label] <- if (index <= length(FUNCTION_COLOURS)) FUNCTION_COLOURS[index] else INK_SECONDARY
    curve_types[label] <- REFIT_LINE_TYPES[min(index, length(REFIT_LINE_TYPES))]
    curve_widths[label] <- REFIT_LINE_WIDTH
  }
  curve_pieces[[NOMINAL_LABEL]] <- data.frame(substitution = rates$substitution, quality = rates$quality,
                                              rate = rates$nominal_rate, curve = NOMINAL_LABEL)[in_range, ]
  curve_colours[NOMINAL_LABEL] <- INK_PRIMARY
  curve_types[NOMINAL_LABEL] <- "dotted"
  curve_widths[NOMINAL_LABEL] <- NOMINAL_LINE_WIDTH
  curves <- do.call(rbind, curve_pieces)
  curve_labels <- names(curve_colours)
  curves$curve <- factor(curves$curve, levels = curve_labels)

  # The y range on the log scale, shared by the 12 panels: the counted rates, the model
  # used and the nominal line, with the triangles of the cells without errors just
  # below the lowest of them. The refits of the other functions do not set it: where
  # the counts are thin, a refit can drop to the floor of its function (1e-7 for the
  # binned fit at a single cell without errors on continuous data, 1e-6 dips of loess
  # on binned data), and one such dip would stretch every panel down by orders of
  # magnitude and squeeze the points and the model into the top of the panels. A refit
  # that runs outside the range is cut off at the panel's edge (coord_cartesian below
  # zooms in without removing data, so the line is cut, not broken into pieces), and
  # the caption says so.
  is_core_curve <- curves$curve %in% c(model_label, NOMINAL_LABEL)
  shown_rates <- c(counted$observed_rate, curves$rate[is_core_curve])
  shown_rates <- shown_rates[is.finite(shown_rates) & shown_rates > 0]
  floor_rate <- min(shown_rates) / 2
  no_error$observed_rate <- rep(floor_rate, nrow(no_error))
  no_error$point_kind <- rep(POINT_NO_ERROR, nrow(no_error))
  y_limits <- c(floor_rate / 1.5, max(shown_rates) * 1.5)
  refit_rates <- curves$rate[!is_core_curve]
  refits_cut <- any(!is.na(refit_rates) & (refit_rates < y_limits[1] | refit_rates > y_limits[2]))

  # The ticks of the quality values of the filtered reads, inside the x range; a value
  # the binned fit left out (no bases in the reads drawn) gets a muted tick.
  read_values <- analysis$read_values
  read_values <- read_values[read_values >= x_low & read_values <= x_high]
  dropped <- int_values(part$anchors_dropped)
  ticks_used <- data.frame(quality = read_values[!(read_values %in% dropped)])
  ticks_dropped <- data.frame(quality = read_values[read_values %in% dropped])

  # Every panel gets the same order and a label that says what the misread is.
  panel_labels <- setNames(paste0(SUBSTITUTIONS, ": a true ", substr(SUBSTITUTIONS, 1, 1),
                                  " read as ", substr(SUBSTITUTIONS, 3, 3)), SUBSTITUTIONS)
  counted$substitution  <- factor(counted$substitution, levels = SUBSTITUTIONS)
  no_error$substitution <- factor(no_error$substitution, levels = SUBSTITUTIONS)
  curves$substitution   <- factor(curves$substitution, levels = SUBSTITUTIONS)

  # The point legend shows only the kinds this figure has (no averaged qualities in
  # continuous data, for example), each in the ink it is drawn in.
  point_kinds <- c(POINT_READ_VALUE, POINT_AVERAGED, POINT_NO_ERROR)
  point_shapes <- setNames(c(16, 1, 6), point_kinds)
  point_inks <- setNames(c(INK_SECONDARY, INK_SECONDARY, INK_MUTED), point_kinds)
  kinds_present <- point_kinds[point_kinds %in% c(counted$point_kind, no_error$point_kind)]

  subtitle <- wrap_text(paste(
    "Points: rates counted in the reads learnErrors drew (size = bases; filled = a quality",
    "value the reads carry; open = an averaged quality of identical reads; triangles at the",
    "bottom = no error counted). Solid line: model used. Dashed: the other functions refitted",
    "on the same counts (smoothing step only). Black dotted: nominal rate of the quality score.",
    "Ticks: quality values in the filtered reads."))
  # The caption lists what a reader of the figure must know about the refits: a failed
  # refit (no curve), extra anchors, a refit cut off at the edge of a panel, and a refit
  # of the function used that does not give back the model.
  caption_lines <- character(0)
  for (function_name in function_names) {
    refit <- analysis$refits[[function_name]]
    if (refit$status != "ok") {
      caption_lines <- c(caption_lines, paste0(function_display_name(function_name), ": refit ", refit$status, " (no curve)"))
    } else if (nzchar(refit$note)) {
      caption_lines <- c(caption_lines, paste0(function_display_name(function_name), ": ", refit$note))
    }
  }
  if (refits_cut) {
    caption_lines <- c(caption_lines, "Dashed curves that run outside a panel are cut off at its edge.")
  }
  if (!isTRUE(analysis$refit_reproduces_model)) {
    caption_lines <- c(caption_lines, "Warning: the refit of the function used does not give back the model exactly.")
  }
  # Each line ends with a full stop, so that the lines read as sentences once joined.
  needs_stop <- !grepl("[.]$", caption_lines)
  caption_lines[needs_stop] <- paste0(caption_lines[needs_stop], ".")
  caption <- if (length(caption_lines) > 0) wrap_text(ascii_text(paste(caption_lines, collapse = " "))) else NULL

  plot <- ggplot() +
    geom_line(data = curves, aes(x = quality, y = rate, colour = curve, linetype = curve, linewidth = curve)) +
    geom_point(data = counted, aes(x = quality, y = observed_rate, size = bases, shape = point_kind),
               colour = INK_SECONDARY, alpha = 0.7) +
    geom_point(data = no_error, aes(x = quality, y = observed_rate, shape = point_kind),
               colour = INK_MUTED, size = 1.2) +
    facet_wrap(~ substitution, ncol = 3, labeller = as_labeller(panel_labels)) +
    scale_colour_manual(name = NULL, values = curve_colours, breaks = curve_labels) +
    scale_linetype_manual(name = NULL, values = curve_types, breaks = curve_labels) +
    scale_linewidth_manual(name = NULL, values = curve_widths, breaks = curve_labels) +
    scale_shape_manual(name = NULL, values = point_shapes[kinds_present], breaks = kinds_present) +
    # Point area grows with the bases at the score, from a smallest dot that stays wider
    # than the solid line of the model: with sizes strictly proportional to the bases
    # (scale_size_area), or a smaller first size, the counted rates at quality scores
    # that hold few bases would vanish under the line. On binned data those are the
    # filled points at 9 and 20 (0.2 to 6 % of the bases), the very points the line
    # must be compared with.
    scale_size(name = "bases at this quality score", range = c(1.6, 4), labels = format_count) +
    scale_x_continuous(limits = c(x_low, x_high), expand = expansion(add = 0.5)) +
    # The log scale without limits, and the range set by coord_cartesian: it zooms in
    # on y_limits and cuts a refit that runs outside at the panel's edge, where scale
    # limits would remove those stretches of the line.
    scale_y_log10() +
    coord_cartesian(ylim = y_limits) +
    guides(colour    = guide_legend(order = 1, nrow = 1, override.aes = list(alpha = 1)),
           linetype  = guide_legend(order = 1, nrow = 1),
           linewidth = guide_legend(order = 1, nrow = 1),
           shape     = guide_legend(order = 2, nrow = 1,
                                    override.aes = list(colour = unname(point_inks[kinds_present]),
                                                        size = 2, alpha = 1)),
           size      = guide_legend(order = 3, nrow = 1)) +
    labs(title = plot_title(direction, part, max_consist), subtitle = subtitle, caption = caption,
         x = "Consensus quality score (mean over identical reads)", y = "Error rate (log scale)") +
    report_theme() +
    theme(legend.key.width = unit(1.4, "cm"))
  if (nrow(ticks_used) > 0) {
    plot <- plot + geom_rug(data = ticks_used, aes(x = quality), sides = "b", colour = INK_PRIMARY,
                            linewidth = 0.3, length = unit(0.03, "npc"))
  }
  if (nrow(ticks_dropped) > 0) {
    plot <- plot + geom_rug(data = ticks_dropped, aes(x = quality), sides = "b", colour = INK_MUTED,
                            linewidth = 0.3, length = unit(0.03, "npc"))
  }
  plot
}


#' The label of one direction's panel in the convergence figure: how the learning
#' ended and whether the cap was reached, on two lines.
#'
#' @param direction "R1" or "R2".
#' @param part The record's entry of the direction.
#' @param max_consist The cap on rounds.
#' @return The label text.
convergence_panel_label <- function(direction, part, max_consist) {
  if (part$stopped_by == "dada2_repeat") {
    outcome <- sprintf("converged in round %d", part$rounds)
  } else if (part$stopped_by == "flat_rounds") {
    outcome <- sprintf("NOT converged: flat-rounds rule, round %d", part$flat_stop_round)
  } else {
    outcome <- "NOT converged: stopped at the cap"
  }
  cap_text <- if (part$rounds >= max_consist) {
    sprintf("cap %d (reached)", max_consist)
  } else {
    sprintf("cap %d (not reached)", max_consist)
  }
  sprintf("%s: %s, %s\n%s", direction, function_display_name(part$error_function), outcome, cap_text)
}


#' The convergence figure: how much the model changed in each round, for both directions.
#'
#' One panel per direction, the round along the bottom, a log scale up the side. Black:
#' the change in the round (the sum, over all rates, of how much the table moved). Grey,
#' dashed: the distance of the table that went into each round to the final model
#' (DADA2's checkConvergence). Both are sums of rate differences, so they share one axis.
#' A muted dashed line marks 1/1000 of the change in round 1, the depth the flat-rounds
#' rule needs. A value of 0 cannot be drawn on a log scale, so zeros are markers at the
#' bottom of the panel: open triangles where a table repeated (DADA2's stop), an open
#' square where the distance is 0 only because MetaFlux handed back the starting table
#' of the round (the flat-rounds rule). A light band marks DADA2's repeat (from the round
#' whose table came back to the last round) or the window of the flat-rounds rule; a
#' dashed line marks the flat-rounds stop and a solid line the cap, when it was reached.
#'
#' @param analyses analyse_direction() of R1 and R2 (their rounds tables).
#' @param record The record.
#' @return A ggplot object.
convergence_figure <- function(analyses, record) {
  max_consist <- record$max_consist
  max_flat_rounds <- record$max_flat_rounds
  series_pieces <- list()
  zero_pieces <- list()
  depth_pieces <- list()
  band_pieces <- list()
  stop_pieces <- list()
  cap_pieces <- list()
  any_flat_stop <- FALSE
  for (direction in c("R1", "R2")) {
    rounds <- analyses[[direction]]$rounds
    if (nrow(rounds) == 0) {
      stop(sprintf("the rounds of %s could not be read (see the warnings)", direction))
    }
    part <- record[[direction]]
    panel <- convergence_panel_label(direction, part, max_consist)
    series_pieces[[paste(direction, "change")]] <- data.frame(
      panel = panel, round = rounds$round, value = rounds$change, series = "change in the round")
    series_pieces[[paste(direction, "distance")]] <- data.frame(
      panel = panel, round = rounds$round, value = rounds$distance_to_final,
      series = "distance to the final model")
    # Zeros: a repeated table gives a change of 0 (a plain repeat) and a distance of 0 at
    # the round whose table came back. After a flat-rounds stop, the last distance is 0
    # because MetaFlux handed back that round's starting table.
    zero_change <- rounds$round[rounds$change == 0]
    zero_distance <- rounds$round[rounds$distance_to_final == 0]
    zero_markers <- data.frame(panel = character(0), round = numeric(0), marker = character(0))
    for (round in zero_change) {
      zero_markers <- rbind(zero_markers, data.frame(panel = panel, round = round, marker = "0 = repeat"))
    }
    for (round in zero_distance) {
      handed_back <- part$stopped_by == "flat_rounds" && round == part$flat_stop_round
      marker <- if (handed_back) "handed back (flat-rounds rule)" else "0 = repeat"
      zero_markers <- rbind(zero_markers, data.frame(panel = panel, round = round, marker = marker))
    }
    zero_pieces[[direction]] <- zero_markers
    depth_pieces[[direction]] <- data.frame(panel = panel, depth = rounds$depth_limit[1])
    if (part$stopped_by == "dada2_repeat") {
      band_label <- if (part$cycle_length == 1) {
        "repeat"
      } else {
        sprintf("repeat of round %d (cycle of %d rounds)", part$matched_round, part$cycle_length)
      }
      band_pieces[[direction]] <- data.frame(panel = panel, xmin = part$matched_round - 0.5,
                                             xmax = part$rounds + 0.5, label = band_label)
    } else if (part$stopped_by == "flat_rounds") {
      any_flat_stop <- TRUE
      first_round_of_window <- max(1, part$flat_stop_round - max_flat_rounds + 1)
      band_pieces[[direction]] <- data.frame(panel = panel, xmin = first_round_of_window - 0.5,
                                             xmax = part$flat_stop_round + 0.5,
                                             label = "flat-rounds rule")
      stop_pieces[[direction]] <- data.frame(panel = panel, round = part$flat_stop_round)
    }
    if (part$rounds >= max_consist) {
      cap_pieces[[direction]] <- data.frame(panel = panel, round = max_consist)
    }
  }
  series <- do.call(rbind, series_pieces)
  zeros <- do.call(rbind, zero_pieces)
  depths <- do.call(rbind, depth_pieces)

  # The y range: every value above 0, the depth lines, and room for the zero markers.
  shown_values <- c(series$value, depths$depth)
  shown_values <- shown_values[is.finite(shown_values) & shown_values > 0]
  floor_value <- min(shown_values) / 3
  zeros$value <- rep(floor_value, nrow(zeros))
  y_limits <- c(floor_value / 1.5, max(shown_values) * 3)
  # A value of 0 becomes a gap in its line (NA), not a dropped row: geom_line joins the
  # rows it is given, so a dropped zero would draw a straight line from the round before
  # to the round after and hide that the value was 0 there (the matched round of a
  # cycle, for example). The marker at the bottom of the panel shows the 0.
  series$value[which(series$value <= 0)] <- NA

  series_colours <- c("change in the round" = INK_PRIMARY, "distance to the final model" = INK_SECONDARY)
  series_types <- c("change in the round" = "solid", "distance to the final model" = "dashed")
  marker_shapes <- c("0 = repeat" = 2, "handed back (flat-rounds rule)" = 0)
  markers_present <- names(marker_shapes)[names(marker_shapes) %in% zeros$marker]

  plot <- ggplot()
  if (length(band_pieces) > 0) {
    bands <- do.call(rbind, band_pieces)
    plot <- plot +
      geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = y_limits[1], ymax = y_limits[2]),
                fill = BAND_LIGHT, alpha = 0.7) +
      geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = y_limits[2], label = label),
                vjust = 1.4, size = 3, colour = INK_SECONDARY)
  }
  plot <- plot +
    geom_hline(data = depths, aes(yintercept = depth), colour = INK_MUTED, linetype = "dashed", linewidth = 0.4) +
    geom_text(data = depths, aes(x = 1, y = depth, label = "1/1000 of round 1"),
              hjust = 0, vjust = -0.5, size = 3, colour = INK_MUTED)
  if (length(stop_pieces) > 0) {
    plot <- plot + geom_vline(data = do.call(rbind, stop_pieces), aes(xintercept = round),
                              colour = INK_SECONDARY, linetype = "dashed", linewidth = 0.5)
  }
  if (length(cap_pieces) > 0) {
    plot <- plot + geom_vline(data = do.call(rbind, cap_pieces), aes(xintercept = round),
                              colour = INK_PRIMARY, linewidth = 0.6)
  }
  # na.rm = TRUE: the gaps of the zeros are intended, so ggplot2's note about removed
  # missing values is left out of the log.
  plot <- plot +
    geom_line(data = series, aes(x = round, y = value, colour = series, linetype = series),
              linewidth = 0.6, na.rm = TRUE) +
    geom_point(data = series, aes(x = round, y = value, colour = series), size = 1.3, na.rm = TRUE)
  if (nrow(zeros) > 0) {
    plot <- plot +
      geom_point(data = zeros, aes(x = round, y = value, shape = marker), colour = INK_PRIMARY, size = 2.2) +
      scale_shape_manual(name = NULL, values = marker_shapes[markers_present], breaks = markers_present)
  }
  caption <- if (any_flat_stop) {
    "The last distance is 0 because MetaFlux handed back the table of that round."
  } else {
    NULL
  }
  plot +
    facet_wrap(~ panel, ncol = 2, scales = "free_x") +
    scale_colour_manual(name = NULL, values = series_colours) +
    scale_linetype_manual(name = NULL, values = series_types) +
    scale_y_log10(limits = y_limits) +
    scale_x_continuous(breaks = function(limits) unique(pmax(1, round(pretty(limits))))) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    labs(title = "How the learning of the error models ended",
         subtitle = wrap_text(paste(
           "Black: how much the model changed in each round (summed over all its rates). Grey,",
           "dashed: how far the table that went into each round is from the final model (DADA2's",
           "checkConvergence). Both on a log scale; a value of 0 breaks its line and is drawn as a",
           "marker at the bottom.",
           "Muted dashed line: 1/1000 of the change in round 1, the depth the flat-rounds rule needs.",
           "A grey band marks DADA2's own stop (a repeat of an earlier table, or a cycle) or the window",
           "of the flat-rounds rule; a vertical line marks the flat-rounds stop (dashed) or the cap.")),
         caption = caption,
         x = "Round of learning", y = "Sum of rate differences (log scale)") +
    report_theme() +
    theme(legend.key.width = unit(1.4, "cm"))
}


#' The quality-values figure: the bases at each quality value of the filtered reads.
#'
#' One panel per direction. A point on a stem for each Phred value, its height the
#' bases at that value in the first 100,000 filtered reads of each voting sample, pooled
#' (choose_error_model's bases_per_value). The axis is logarithmic, so a bar's length
#' would mean nothing: the point carries the value, the stem only leads the eye. Values
#' at or below Q2 (no-call values) are muted. With ten values or fewer, each is labelled
#' with its share of the bases. The panel label gives the number of values above Q2, the
#' kind of quality values they make, and the error function chosen, for example
#' "R1: 3 values above Q2 (9, 20, 34) -> binned (auto)".
#'
#' @param choice error_model_choice.json, read.
#' @return A ggplot object.
quality_figure <- function(choice) {
  pieces <- list()
  for (direction in c("R1", "R2")) {
    chosen <- choice[[direction]]
    bases_per_value <- chosen$bases_per_value
    if (is.null(bases_per_value) || length(bases_per_value) == 0) {
      stop(sprintf("error_model_choice.json holds no bases per quality value for %s", direction))
    }
    values <- as.integer(names(bases_per_value))
    bases <- num_values(bases_per_value)
    above_q2 <- sort(values[values > 2])
    panel <- sprintf("%s: %d values above Q2 (%s) -> %s (%s)", direction, length(above_q2),
                     values_text(above_q2), chosen$error_function, chosen$decided_by)
    label <- if (length(values) <= 10) sprintf("%.1f %%", 100 * bases / sum(bases)) else ""
    pieces[[direction]] <- data.frame(panel = panel, quality = values, bases = bases, label = label,
                                      group = ifelse(values <= 2, "Q2 and below (no-call values)", "above Q2"))
  }
  values <- do.call(rbind, pieces)
  values <- values[values$bases > 0, ]
  # The stems start at the power of ten below the smallest count (one decade lower when
  # the smallest count is itself a power of ten), so every stem has a length.
  stem_floor <- 10^(ceiling(log10(min(values$bases))) - 1)
  # The two inks get a legend only when both kinds of value are there (one kind needs no
  # legend: the subtitle says what is drawn).
  group_inks <- c("above Q2" = INK_PRIMARY, "Q2 and below (no-call values)" = INK_MUTED)
  groups_present <- names(group_inks)[names(group_inks) %in% values$group]
  legend_position <- if (length(groups_present) > 1) "bottom" else "none"
  ggplot(values, aes(x = quality, y = bases, colour = group)) +
    geom_segment(aes(xend = quality, y = stem_floor, yend = bases), linewidth = 0.5) +
    geom_point(size = 2.2) +
    geom_text(aes(label = label), vjust = -0.9, size = 3, colour = INK_SECONDARY) +
    facet_wrap(~ panel, ncol = 2) +
    scale_colour_manual(name = NULL, values = group_inks[groups_present], breaks = groups_present) +
    scale_y_log10(limits = c(stem_floor, max(values$bases) * 4)) +
    scale_x_continuous(limits = c(0, max(values$quality) + 1), breaks = seq(0, 45, by = 5)) +
    labs(title = "Quality values in the filtered reads",
         subtitle = wrap_text(paste(
           "Bases per quality value in the first 100,000 filtered reads of each voting sample, pooled.",
           "3 to 10 values above Q2 with gaps between them are read as binned quality scores,",
           "20 or more as continuous ones;",
           "error_model: auto chooses the error function from this.")),
         x = "Phred quality value", y = "Bases (log scale)") +
    report_theme() +
    theme(legend.position = legend_position)
}


#' DADA2's own picture of one model (plotErrors), with the report's title.
#'
#' Points: the counted rates; black line: the model; red line: the rate the quality
#' score promises (nominalQ = TRUE). MetaFlux has drawn it since v1.0.0, and it is kept
#' for users who know it. Its log10 warnings, raised when the picture is saved, come from
#' cells without errors (a rate of 0 has no logarithm); they are cosmetic, so draw_figure
#' leaves them out of the log for this picture (hide_warnings).
#'
#' @param err The model as learnErrors returned it (after a flat-rounds stop, a copy whose
#'   trans holds the counts the kept table was fitted to; see where it is called).
#' @param title The title (plot_title()).
#' @return A ggplot object.
dada2_plot <- function(err, title) {
  suppressWarnings(plotErrors(err, nominalQ = TRUE)) + ggtitle(title)
}


#' A picture that stands in for a figure that could not be drawn, and says why.
#'
#' @param what The figure, as the warning names it.
#' @param problem The error message.
#' @return A ggplot object with the text only.
placeholder_figure <- function(what, problem) {
  first_line <- paste0(toupper(substr(what, 1, 1)), substr(what, 2, nchar(what)), " could not be drawn.")
  text <- paste0(
    wrap_text(first_line, 80), "\n\n",
    wrap_text(ascii_text(problem), 80), "\n\n",
    wrap_text(paste("The numbers are in error_model_fit.tsv, error_model_rates.tsv and",
                    "error_model_rounds.tsv in stats/dada2/, and the reason is also in",
                    "logs/error_model_report.log."), 80))
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = text, colour = INK_PRIMARY, size = 4.5) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA))
}


#' Save one figure as PNG (150 dpi) and PDF.
save_figure <- function(plot, png_path, pdf_path, width_cm, height_cm) {
  ggsave(png_path, plot = plot, device = "png", width = width_cm, height = height_cm,
         units = "cm", dpi = 150, bg = "white")
  ggsave(pdf_path, plot = plot, device = "pdf", width = width_cm, height = height_cm, units = "cm")
}


#' Draw one figure, or a placeholder that states why it could not be drawn.
#'
#' ggplot2 builds a figure only when it is saved, so building and saving are tried
#' together. When either fails, the placeholder is saved under the same names and the
#' failure becomes a warning of the report.
#'
#' @param what The figure, for the warning ("the fit figure of the R1 reads").
#' @param build A function without arguments that returns the figure.
#' @param png_path,pdf_path The two files.
#' @param width_cm,height_cm The size.
#' @param direction "R1", "R2" or "both", for the warning.
#' @param hide_warnings TRUE to keep R's warnings of building and saving out of the log
#'   (DADA2's own plots only, whose warnings are known to be cosmetic).
#' @return Nothing.
draw_figure <- function(what, build, png_path, pdf_path, width_cm, height_cm, direction = "both",
                        hide_warnings = FALSE) {
  draw <- function() save_figure(build(), png_path, pdf_path, width_cm, height_cm)
  problem <- tryCatch({
    if (hide_warnings) suppressWarnings(draw()) else draw()
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(problem)) {
    add_warning(sprintf("%s could not be drawn: %s", what, problem), direction)
    save_figure(placeholder_figure(what, problem), png_path, pdf_path, width_cm, height_cm)
  }
}


# ── Step 4c: the files of the MultiQC section ─────────────────────

#' Write a MultiQC custom-content table: a commented YAML header, then the table.
#'
#' MultiQC reads the lines starting with "#" as YAML: the section id, the module it
#' belongs to (parent_id; the module's name and description come from the image
#' entries in workflow/multiqc_config.yaml), its name, its description, that it is a
#' table (kept in the order of the file: sort_rows False), and per column a
#' description (MultiQC shows it in its column settings) and scale: False (MultiQC would
#' otherwise draw a coloured bar behind every number, which means nothing in a table of
#' two or four rows). The first column of the table is what MultiQC shows as the row
#' name.
#'
#' @param table A data frame whose column names are the headers MultiQC shows.
#' @param path The file (in stats/dada2/multiqc/, which MultiQC scans).
#' @param section_id,section_name,description The section's id, name and description.
#' @param column_descriptions Named character vector: a description per column (all
#'   columns but the first).
#' @param column_titles Named character vector: a shorter title to show for some
#'   columns (the column name stays the column's id in MultiQC's data files).
#' @return Nothing.
write_multiqc_table <- function(table, path, section_id, section_name, description,
                                column_descriptions, column_titles = character(0)) {
  header <- c(
    sprintf("# id: '%s'", section_id),
    "# parent_id: 'metaflux_error_model'",
    sprintf("# section_name: '%s'", section_name),
    sprintf("# description: \"%s\"", description),
    "# plot_type: 'table'",
    "# pconfig:",
    sprintf("#     id: '%s_table'", section_id),
    "#     namespace: 'DADA2 error model'",
    "#     sort_rows: False",
    "# headers:")
  for (column in names(table)[-1]) {
    header <- c(header,
                sprintf("#     \"%s\":", column),
                sprintf("#         description: \"%s\"", column_descriptions[[column]]),
                "#         scale: False")
    if (column %in% names(column_titles)) {
      header <- c(header, sprintf("#         title: \"%s\"", column_titles[[column]]))
    }
  }
  lines <- c(header, paste(names(table), collapse = "\t"))
  for (row_index in seq_len(nrow(table))) {
    cells <- vapply(table[row_index, ], function(value) cell_text(as.character(value)), character(1))
    lines <- c(lines, paste(cells, collapse = "\t"))
  }
  writeLines(lines, path)
}


#' The rows of the MultiQC summary: one per read direction.
#'
#' @param record The record.
#' @param choice error_model_choice.json, read (NULL when it could not be read).
#' @return A data frame with the columns MultiQC shows.
summary_rows <- function(record, choice) {
  rows <- list()
  for (direction in c("R1", "R2")) {
    part <- record[[direction]]
    anchors <- int_values(part$anchors_used)
    # MetaFlux's own functions by their registry name, the user's own with its name. The
    # anchors of a binned fit as in the figure titles (more than ten: their number and
    # range, anchor_list_text).
    function_text <- if (identical(part$error_function, "custom")) function_display_name("custom") else part$error_function
    if (length(anchors) > 0 && !all(is.na(anchors))) {
      function_text <- paste0(function_text, " on ", anchor_list_text(anchors))
    }
    values_above <- if (is.null(choice)) integer(0) else int_values(choice[[direction]]$values_above_q2)
    above_text <- if (length(values_above) == 0) "-" else sprintf("%d (%s)", length(values_above), values_text(values_above))
    changes <- num_values(part$change_per_round)
    shown_rounds <- unique(c(1, 2, length(changes)))
    shown_rounds <- shown_rounds[shown_rounds <= length(changes)]
    # learn_errors reads the three numbers from DADA2's own line in its log; any of
    # them can be missing (null) when that line could not be read.
    learned_parts <- character(0)
    if (length(num_values(part$bases_learned)) == 1 && !is.na(part$bases_learned)) {
      learned_parts <- c(learned_parts, paste(format_count(part$bases_learned), "bases"))
    }
    if (length(num_values(part$reads_learned)) == 1 && !is.na(part$reads_learned)) {
      learned_parts <- c(learned_parts, paste("in", format_count(part$reads_learned), "reads"))
    }
    if (length(num_values(part$samples_learned)) == 1 && !is.na(part$samples_learned)) {
      sample_word <- if (part$samples_learned == 1) "sample" else "samples"
      learned_parts <- c(learned_parts, paste("from", format_count(part$samples_learned), sample_word))
    }
    learned <- if (length(learned_parts) == 0) "not recorded" else paste(learned_parts, collapse = " ")
    # The warnings of the learning (in the record) and those of this report that concern
    # this direction or both.
    report_warnings <- report_state$warning_texts[report_state$warning_directions %in% c(direction, "both")]
    warnings <- c(text_values(part$warnings), report_warnings)
    row <- data.frame(
      "Read direction"                 = direction,
      "Error function"                 = function_text,
      "Chosen by"                      = if (is.null(part$decided_by)) "-" else part$decided_by,
      "Quality values above Q2"        = above_text,
      "Rounds"                         = part$rounds,
      "Cap"                            = record$max_consist,
      "How learning ended"             = short_outcome_text(part, record$max_consist),
      "Converged"                      = if (isTRUE(part$converged)) "yes" else "no",
      "Change in rounds 1, 2 and last" = paste(format_number(changes[shown_rounds]), collapse = ", "),
      "Bases learned from"             = learned,
      "Model code"                     = part$model_code,
      "Warnings"                       = if (length(warnings) == 0) "none" else paste(warnings, collapse = " | "),
      check.names = FALSE, stringsAsFactors = FALSE)
    # One more column only in a run that names a function of the user's own
    # (experimental), so that the summary of every other run stays as it was.
    if (length(custom_function_file_input) > 0) {
      row[["Custom error function"]] <- custom_summary_text(part)
    }
    rows[[direction]] <- row
  }
  do.call(rbind, rows)
}


#' The MultiQC summary's cell about the user's own error function (experimental): its
#' name, the first 8 characters of its file's md5 sum (they change with any edit of the
#' file), and whether this model was learned with it, such as
#' "gam_errfun (md5 1a2b3c4d): used to learn this model" or
#' "gam_errfun (md5 1a2b3c4d): refitted for comparison only".
#'
#' @param part The record's entry of the direction.
#' @return The text.
custom_summary_text <- function(part) {
  md5_text <- if (is.null(custom_function)) "file not read" else paste("md5", substr(custom_function$md5, 1, 8))
  role <- if (identical(part$error_function, "custom")) "used to learn this model" else "refitted for comparison only"
  sprintf("%s (%s): %s", custom_function_name, md5_text, role)
}


#' The anchors of a binned fit as the MultiQC table shows them: "9, 20, 34", or
#' "32 values (7-38)" when there are more than ten; "-" for a fit without anchors.
#'
#' @param anchors The anchors column of the fit rows ("9, 20, 34" or NA).
#' @return The text.
anchors_text <- function(anchors) {
  if (is.na(anchors)) return("-")
  anchor_list_text(as.integer(strsplit(anchors, ", ")[[1]]))
}


#' The rows of the MultiQC goodness-of-fit table: one per direction and function.
#'
#' @param fit_table The fit rows of both directions (fit_rows_of_direction(): unrounded,
#'   with the short notes).
#' @return A data frame with the columns MultiQC shows.
goodness_of_fit_rows <- function(fit_table) {
  rows <- list()
  for (row_index in seq_len(nrow(fit_table))) {
    fit <- fit_table[row_index, ]
    rising <- if (is.na(fit$rising_steps)) "-" else sprintf("%d of %d", as.integer(fit$rising_steps), as.integer(fit$steps_compared))
    # MetaFlux's own functions by their registry name ("R1 | loess"), the user's own with
    # its name ("R1 | custom (gam_errfun)").
    function_label <- if (fit$error_function == "custom") function_display_name("custom") else fit$error_function
    rows[[row_index]] <- data.frame(
      "Read direction and function"                     = paste(fit$read_direction, "|", function_label),
      # Not yes/no: MultiQC paints every "no" in a table as a red badge, which would
      # mark the functions refitted only for comparison as failures.
      "Model used"                                      = if (isTRUE(fit$model_used)) "used" else "refit only",
      "Anchors"                                         = anchors_text(fit$anchors),
      "Rising steps"                                    = rising,
      "Substitution types with a rise"                  = if (is.na(fit$substitutions_with_rise)) "-" else sprintf("%d of 12", as.integer(fit$substitutions_with_rise)),
      "Typical difference at the reads' quality values" = format_fold(fit$typical_fold_difference_at_read_values),
      "Typical difference at averaged columns"          = format_fold(fit$typical_fold_difference_at_averaged_columns),
      "Note"                                            = if (nzchar(fit$note_short)) fit$note_short else "-",
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) {
    return(data.frame("Read direction and function" = "none", "Note" = "the goodness of fit could not be computed",
                      check.names = FALSE, stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}


#' A copy of a table with its decimal numbers rounded to six significant digits, as the
#' TSV and JSON files carry them. Counts (whole numbers) are left exact.
round_for_files <- function(table, columns) {
  for (column in intersect(columns, names(table))) {
    table[[column]] <- signif(as.numeric(table[[column]]), 6)
  }
  # A deviance that is 0 by construction (the binned fit at the reads' own values) comes
  # out of the arithmetic as a tiny remainder such as 4.2e-14; anything below one
  # millionth per million bases is written as 0.
  for (column in grep("^deviance", intersect(columns, names(table)), value = TRUE)) {
    tiny <- !is.na(table[[column]]) & abs(table[[column]]) < 1e-6
    table[[column]][tiny] <- 0
  }
  table
}


#' A deviance for the log: 0 when it is below one millionth (the rule round_for_files()
#' applies to the tables), so that the log and error_model_fit.tsv agree.
zero_if_tiny <- function(x) ifelse(!is.na(x) & abs(x) < 1e-6, 0, x)


# ══ Main ══════════════════════════════════════════════════════════

# ── Read the record, the models and the choice ────────────────────
# The record and the two models are what the report is about: when one cannot be read,
# the job stops (the only reason it ever does).
read_or_stop <- function(read, path, what) {
  tryCatch(read(path), error = function(e) {
    stop_with_message(sprintf("[error_model_report] STOP: %s (%s) could not be read: %s",
                              what, path, conditionMessage(e)))
  })
}
record <- read_or_stop(function(path) fromJSON(path, simplifyVector = TRUE),
                       snakemake@input[["record"]], "the record of learn_errors")
models <- list(
  R1 = read_or_stop(readRDS, snakemake@input[["model_r1"]], "the R1 error model"),
  R2 = read_or_stop(readRDS, snakemake@input[["model_r2"]], "the R2 error model"))
# The choice is needed only for the quality-values figure and one column of the summary.
choice <- tryCatch(fromJSON(snakemake@input[["choice"]], simplifyVector = TRUE), error = function(e) {
  add_warning(paste("error_model_choice.json could not be read:", conditionMessage(e)))
  NULL
})
# The user's own error function (experimental), when its file could not be read at the
# top of this script: the report goes on without its curve, and says why.
if (!is.null(custom_load_problem)) {
  add_warning(paste("your own error function (amplicon.dada2.learn_errors.custom_error_function)",
                    "could not be read, so it has no curve:", custom_load_problem))
}
if (!is.null(custom_function)) {
  message(sprintf("[error_model_report] %s from %s (md5 %s), refitted with the others%s",
                  custom_function_label(custom_function$name), custom_function$file,
                  substr(custom_function$md5, 1, 8),
                  if (identical(record$requested, "custom")) "; the models were learned with it" else ", for comparison only"))
}
# The run's output folder (the report JSON is <out>/stats/dada2/error_model_report.json;
# OUT is an absolute path), to fill in the paths of a refit's reason.
output_folder <- dirname(dirname(dirname(snakemake@output[["report_json"]])))
max_consist <- record$max_consist
max_flat_rounds <- if (is.null(record$max_flat_rounds)) 0 else record$max_flat_rounds

dir.create(dirname(snakemake@output[["report_json"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake@output[["mqc_summary"]]), recursive = TRUE, showWarnings = FALSE)

# ── Steps 1-3: counts, refits and measures, per direction ─────────
analyses <- list()
for (direction in c("R1", "R2")) {
  part <- record[[direction]]
  message("[error_model_report] ", plot_title(direction, part, max_consist))
  analyses[[direction]] <- analyse_direction(direction, part, models[[direction]], max_flat_rounds)
  # One line per function in the log, for a quick look without opening the tables. The
  # deviances go through zero_if_tiny() (above), so the log shows the 0 the tables show.
  fit_rows <- analyses[[direction]]$fit_rows
  for (row_index in seq_len(nrow(fit_rows))) {
    fit <- fit_rows[row_index, ]
    message(sprintf(paste("[error_model_report] %s %s%s: %s; deviance per million bases %s",
                          "(read values %s, averaged columns %s); typical difference %s",
                          "(read values %s, averaged columns %s); rising steps %s of %s"),
                    direction, function_display_name(fit$error_function), if (isTRUE(fit$model_used)) " (model used)" else "",
                    fit$refit_status, format_number(zero_if_tiny(fit$deviance_per_million_bases)),
                    format_number(zero_if_tiny(fit$deviance_per_million_bases_at_read_values)),
                    format_number(zero_if_tiny(fit$deviance_per_million_bases_at_averaged_columns)),
                    format_fold(fit$typical_fold_difference),
                    format_fold(fit$typical_fold_difference_at_read_values),
                    format_fold(fit$typical_fold_difference_at_averaged_columns),
                    fit$rising_steps, fit$steps_compared))
  }
}

# ── Step 4b: the figures ──────────────────────────────────────────
# Each in its own tryCatch (draw_figure), so a figure that fails leaves a placeholder
# and a warning, never a failed job. Sizes in cm; PNG at 150 dpi.
output <- snakemake@output
for (direction in c("R1", "R2")) {
  key <- tolower(direction)
  draw_figure(sprintf("the fit figure of the %s reads", direction),
              function() fit_figure(direction, analyses[[direction]], record[[direction]], max_consist),
              output[[paste0("fit_png_", key)]], output[[paste0("fit_pdf_", key)]], 30, 25, direction)
}
draw_figure("the convergence figure", function() convergence_figure(analyses, record),
            output[["convergence_png"]], output[["convergence_pdf"]], 30, 20)
draw_figure("the quality-values figure", function() {
              if (is.null(choice)) stop("error_model_choice.json could not be read")
              quality_figure(choice)
            },
            output[["quality_png"]], output[["quality_pdf"]], 30, 15)
# DADA2's own plots, under the names they have had since v1.0.0. plotErrors draws the
# points from the model's trans (the counts of the last round) and the line from its
# err_out. After a flat-rounds stop these belong to different rounds: the kept table was
# fitted to the counts of the round before the last (counts_behind_model), which DADA2
# does not keep. So for such a model the plot is drawn from a copy whose trans is
# replaced by those counts, as the fit figures use them; the saved .rds stays as
# learnErrors returned it.
for (direction in c("R1", "R2")) {
  key <- tolower(direction)
  draw_figure(sprintf("DADA2's error plot of the %s reads", direction),
              function() {
                shown <- models[[direction]]
                flat_stop <- identical(record[[direction]]$stopped_by, "flat_rounds")
                if (flat_stop && !is.null(analyses[[direction]]$counts)) {
                  shown$trans <- analyses[[direction]]$counts
                }
                dada2_plot(shown, plot_title(direction, record[[direction]], max_consist))
              },
              output[[paste0("dada2_plot_", key)]], output[[paste0("dada2_plot_", key, "_pdf")]], 30, 20, direction,
              hide_warnings = TRUE)
}

# ── Step 4a: the tables ───────────────────────────────────────────
fit_table <- rbind(analyses$R1$fit_rows, analyses$R2$fit_rows)
rates_table <- rbind(analyses$R1$rates, analyses$R2$rates)
rounds_table <- rbind(analyses$R1$rounds, analyses$R2$rounds)
fit_numbers <- c("share_bases_at_read_values", "loglik_per_base", "deviance_per_million_bases",
                 "deviance_per_million_bases_at_read_values",
                 "deviance_per_million_bases_at_averaged_columns", "typical_fold_difference",
                 "typical_fold_difference_at_read_values",
                 "typical_fold_difference_at_averaged_columns")
rate_numbers <- c("observed_rate", "nominal_rate", "model_rate", paste0("rate_", names(error_function_registry)))
round_numbers <- c("change", "lowest_change_so_far", "distance_to_final", "depth_limit")
# The short notes are for the MultiQC table only (fit_note()).
fit_table_rounded <- round_for_files(fit_table[, FIT_TABLE_COLUMNS, drop = FALSE], fit_numbers)
write_tsv(fit_table_rounded, output[["fit_tsv"]])
write_tsv(round_for_files(rates_table, rate_numbers), output[["rates_tsv"]])
write_tsv(round_for_files(rounds_table, round_numbers), output[["rounds_tsv"]])

# The report JSON: the statement, per direction how the learning ended and the fit rows
# (the same fields and values as error_model_fit.tsv), and every warning of this report.
# Written after the figures, so that the warnings of a failed figure are in it.
report <- list(statement = REPORT_STATEMENT)
for (direction in c("R1", "R2")) {
  part <- record[[direction]]
  direction_rows <- fit_table_rounded[fit_table_rounded$read_direction == direction, , drop = FALSE]
  rownames(direction_rows) <- NULL
  report[[direction]] <- list(
    learning_outcome                        = learning_outcome_text(part, max_consist, max_flat_rounds),
    model_code                              = part$model_code,
    error_function                          = part$error_function,
    share_bases_at_read_values              = signif(analyses[[direction]]$share_bases_at_read_values, 6),
    share_substitutions_at_averaged_columns = signif(analyses[[direction]]$share_substitutions_at_averaged_columns, 6),
    fit                                     = direction_rows)
}
report$warnings <- I(report_state$warning_texts)
write_json(report, output[["report_json"]], auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")

# ── Step 4c: the MultiQC section ──────────────────────────────────
# Two tables and copies of the four figures, in stats/dada2/multiqc/, a folder MultiQC
# scans (90_report.smk). The copies get readable names ending in _mqc.png, which
# MultiQC shows as image sections; workflow/multiqc_config.yaml gives them their names,
# descriptions and order, and puts all six sections in one module, "DADA2 error model".
write_multiqc_table(summary_rows(record, choice), output[["mqc_summary"]],
                    "metaflux_error_model_summary", "Summary", SUMMARY_DESCRIPTION,
                    SUMMARY_COLUMN_DESCRIPTIONS)
write_multiqc_table(goodness_of_fit_rows(fit_table), output[["mqc_fit_table"]],
                    "metaflux_error_model_goodness_of_fit", "Goodness of fit",
                    GOODNESS_OF_FIT_DESCRIPTION, GOODNESS_OF_FIT_COLUMN_DESCRIPTIONS,
                    GOODNESS_OF_FIT_COLUMN_TITLES)
copies <- c(mqc_fit_r1 = "fit_png_r1", mqc_fit_r2 = "fit_png_r2",
            mqc_convergence = "convergence_png", mqc_quality = "quality_png")
for (copy_key in names(copies)) {
  copied <- file.copy(output[[copies[[copy_key]]]], output[[copy_key]], overwrite = TRUE)
  if (!copied) {
    stop_with_message(sprintf("[error_model_report] STOP: could not copy %s to %s.",
                              output[[copies[[copy_key]]]], output[[copy_key]]))
  }
}

message(sprintf("[error_model_report] DONE. Report written with %d warning(s).",
                length(report_state$warning_texts)))

# ── Close log ────────────────────────────────────────────────────
sink(type = "message")
sink(type = "output")
close(log_con)
