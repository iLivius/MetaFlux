# Checks of the R scripts of DADA2's error learning that need R but no real reads: the
# order in which the learning step sorts the filtered files, how dada_seqtab pairs each
# R1 file with its R2 file, three texts (the flat-rounds line of the learning log, the
# reason the binned fit stops on quality scores outside its anchors, and the figure
# titles and MultiQC summary of the report), and the two files dada_filter writes, on a
# made-up read pair. A change that breaks one of these fails a test instead of surfacing
# in a real run.
#
# Who runs it: tests/test_error_model_r.py, only when the environment variable
# METAFLUX_RSCRIPT names an Rscript (the CI runs no R, so the test is skipped there).
# Called with two arguments: the folder workflow/scripts/amplicon of the repository, and
# an empty folder where the dada_filter case writes its files (the test reads them
# there; without it, a temporary folder is used).
#
# What it needs: base R for every case but two. The case that runs dada_seqtab's script
# up to its stop for a missing R2 file loads the script's packages (pacman and dada2),
# and the case that runs dada_filter's script whole loads its own (pacman, dada2,
# jsonlite and ShortRead), all in MetaFlux's DADA2 environment; without them a case
# prints SKIP.
#
# What it loads. Sourcing 60_error_functions.R only defines functions. The other three
# scripts read the snakemake object of a running job, so they are parsed (read as R code
# without running it) and only some of their top-level statements are run, as
# .test/custom_error_functions/custom_function_checks.R does for 60d_learn_errors.R:
#   - from 60d_learn_errors.R, 60e_error_model_report.R and 60g_dada_seqtab.R: every
#     function definition ("name <- function(...) ...");
#   - from 60d_learn_errors.R: also its two constants of the flat-rounds rule, and the
#     statements that sort and name the filtered files, run on made-up paths.
# Two scripts are also run whole, with a stand-in snakemake object: 60g_dada_seqtab.R up
# to its stop for a missing R2 file (case 2), and 60b_dada_filter.R on a made-up read
# pair (case 6).
#
# Output: one line per case on stdout, "<case><TAB><OK, SKIP: reason, or what was found
# instead>", which the test reads.

args <- commandArgs(trailingOnly = TRUE)
scripts_folder <- args[1]
source(file.path(scripts_folder, "60_error_functions.R"))

# One line of output.
report <- function(case, text) {
  cat(case, "\t", gsub("\n", " ", text), "\n", sep = "")
}

# The top-level statements of an R file whose deparsed text starts with one of `starts`
# (for example "filt_r1 <- "), in file order. parse() reads the file as code without
# running it; deparse() turns a statement back into text.
statements_starting_with <- function(path, starts) {
  found <- list()
  for (statement in parse(path, keep.source = FALSE)) {
    text <- paste(deparse(statement, width.cutoff = 500L), collapse = " ")
    if (any(startsWith(text, starts))) found[[length(found) + 1]] <- statement
  }
  found
}

# Run the function definitions of an R file (and the constants named in `constants`) in
# the global environment, and nothing else of it.
load_definitions <- function(path, constants = character(0)) {
  for (statement in parse(path, keep.source = FALSE)) {
    is_assignment <- is.call(statement) && identical(statement[[1]], as.name("<-")) &&
      is.name(statement[[2]])
    if (!is_assignment) next
    assigned_value <- statement[[3]]
    is_function <- is.call(assigned_value) && identical(assigned_value[[1]], as.name("function"))
    if (is_function || as.character(statement[[2]]) %in% constants) eval(statement, globalenv())
  }
}

load_definitions(file.path(scripts_folder, "60d_learn_errors.R"),
                 constants = c("FLAT_DEPTH_FACTOR", "FLAT_MIN_GAIN"))
load_definitions(file.path(scripts_folder, "60e_error_model_report.R"))
load_definitions(file.path(scripts_folder, "60g_dada_seqtab.R"))


# ── 1. learn_errors sorts the filtered files byte by byte, in every locale ──
# The statements of 60d that sort and name the files are run on made-up paths, once
# with the collation of a C locale and once with that of a UTF-8 locale. The sample
# names mix capitals, small letters, "-" and "_", so R's default sort orders them
# differently in the two. That is checked below: where the default sort gives one order
# in both, the case is skipped, so it cannot pass by accident.
SnakemakeStandIn <- setClass("SnakemakeStandIn", representation(input = "list"))
sort_statements <- statements_starting_with(
  file.path(scripts_folder, "60d_learn_errors.R"),
  c("filt_r1 <- ", "filt_r2 <- ", "sample_names <- ", "names(filt_r1) <- ", "names(filt_r2) <- "))
folder <- "/run/out/04.filtered/"
run_sort_statements <- function(samples) {
  snakemake <<- SnakemakeStandIn(input = list(
    r1 = as.list(paste0(folder, samples, "_R1_filt.fastq.gz")),
    r2 = as.list(paste0(folder, rev(samples), "_R2_filt.fastq.gz"))))
  for (statement in sort_statements) eval(statement, globalenv())
  list(r1 = filt_r1, r2 = filt_r2)
}
utf8_locale <- ""
for (candidate in c("C.UTF-8", "C.utf8", "en_US.UTF-8", "en_US.utf8")) {
  utf8_locale <- suppressWarnings(Sys.setlocale("LC_COLLATE", candidate))
  if (nzchar(utf8_locale)) break
}
if (!nzchar(utf8_locale)) {
  report("sort_same_in_every_locale", "SKIP: no UTF-8 locale can be set")
} else {
  samples <- c("soil-2", "Soil_10", "Soil1")
  default_utf8 <- sort(samples)
  sorted_utf8 <- run_sort_statements(samples)
  invisible(Sys.setlocale("LC_COLLATE", "C"))
  default_c <- sort(samples)
  sorted_c <- run_sort_statements(samples)
  expected_names <- c("Soil1", "Soil_10", "soil-2")       # byte order (C)
  if (identical(default_c, default_utf8)) {
    # The default sort gives one order in both locales, so 60d's sort could not be
    # told apart from it here: skip rather than fail, because the code under test is
    # not at fault. This happens when LC_ALL or LC_COLLATE is C in R's environment: R
    # then stops using its own (ICU) alphabetical order, even after Sys.setlocale()
    # above, and sorts these names alike in both locales (tests/test_error_model_r.py
    # runs Rscript without the two variables for that reason).
    report("sort_same_in_every_locale",
           paste0("SKIP: R's default sort orders the sample names alike under C and ",
                  utf8_locale, " in this R session (as when LC_ALL or LC_COLLATE is C ",
                  "in its environment), so this case can show nothing"))
  } else {
    problems <- character(0)
    if (!identical(unname(sorted_c$r1), unname(sorted_utf8$r1))) problems <- c(problems, "R1 order depends on the locale")
    if (!identical(unname(sorted_c$r2), unname(sorted_utf8$r2))) problems <- c(problems, "R2 order depends on the locale")
    if (!identical(names(sorted_c$r1), expected_names)) problems <- c(problems, paste("R1 names", paste(names(sorted_c$r1), collapse = ", ")))
    if (!identical(names(sorted_c$r2), expected_names)) problems <- c(problems, paste("R2 names", paste(names(sorted_c$r2), collapse = ", ")))
    report("sort_same_in_every_locale", if (length(problems) == 0) "OK" else paste(problems, collapse = "; "))
  }
}

# Each R2 file carries its own sample's name, also where byte order puts the R2 list in
# another order than the R1 list ("Soil" and "Soil_R2": "Soil, Soil_R2" for R1, but
# "Soil_R2, Soil" for R2, because a capital R comes before a small f).
invisible(Sys.setlocale("LC_COLLATE", "C"))
sorted <- run_sort_statements(c("Soil", "Soil_R2"))
r2_own_names <- sub("_R2_filt\\.fastq\\.gz$", "", basename(sorted$r2))
if (identical(unname(names(sorted$r2)), r2_own_names) &&
    identical(r2_own_names, c("Soil_R2", "Soil")) && identical(names(sorted$r1), c("Soil", "Soil_R2"))) {
  report("r2_files_named_after_their_own_sample", "OK")
} else {
  report("r2_files_named_after_their_own_sample",
         paste("R1:", paste(names(sorted$r1), collapse = ", "), "| R2 files:", paste(r2_own_names, collapse = ", "),
               "| R2 names:", paste(names(sorted$r2), collapse = ", ")))
}


# ── 2. dada_seqtab pairs each R1 file with its R2 file by name ──
# Under a C locale R's default sort puts the R1 files of "Soil" and "Soil_R2" as "Soil,
# Soil_R2" and their R2 files as "Soil_R2, Soil"; pairing by position, as up to v2.6.1,
# hands mergePairs each sample's R1 reads with the other's R2 reads. DADA2 then stops
# with an error when the two samples have different read counts, and merges the wrong
# reads without one when their read counts are equal.
r1_sorted <- sort(paste0(folder, c("Soil_R2", "Soil"), "_R1_filt.fastq.gz"))
r2_inputs <- paste0(folder, c("Soil", "Soil_R2"), "_R2_filt.fastq.gz")
partners <- r2_partners_of(r1_sorted, r2_inputs)
by_position <- sort(r2_inputs)
expected_partners <- paste0(folder, c("Soil", "Soil_R2"), "_R2_filt.fastq.gz")
if (!identical(partners, expected_partners)) {
  report("seqtab_pairs_by_name", paste("got", paste(basename(partners), collapse = ", ")))
} else if (identical(by_position, expected_partners)) {
  report("seqtab_pairs_by_name", "the sort by position pairs these names right as well, so this case proves nothing")
} else {
  report("seqtab_pairs_by_name", "OK")
}
partners <- r2_partners_of(paste0(folder, c("A", "B"), "_R1_filt.fastq.gz"), paste0(folder, "A_R2_filt.fastq.gz"))
report("seqtab_missing_partner_is_na",
       if (identical(partners, c(paste0(folder, "A_R2_filt.fastq.gz"), NA))) "OK" else paste(partners, collapse = ", "))

# The stop itself: 60g run whole, as Snakemake runs it, with the R2 file of sample B
# missing from its inputs. It must stop before reading any file, with its text in the
# log. Needs the packages 60g loads; it would try to install pacman otherwise.
if (!requireNamespace("pacman", quietly = TRUE) || !requireNamespace("dada2", quietly = TRUE)) {
  report("seqtab_stops_without_partner", "SKIP: the R packages pacman and dada2 are needed")
} else {
  Snakemake <- setClass("Snakemake", representation(input = "list", output = "list", params = "list",
                                                     log = "list", threads = "numeric"))
  work <- tempfile("seqtab_stop_")
  dir.create(work)
  log_file <- file.path(work, "dada_seqtab.log")
  snakemake <- Snakemake(
    input = list(r1 = as.list(paste0(work, "/", c("A", "B"), "_R1_filt.fastq.gz")),
                 r2 = list(paste0(work, "/A_R2_filt.fastq.gz")),
                 filter_stats = list(), verdict = file.path(work, "error_model_check.txt")),
    output = list(), params = list(), log = list(log_file), threads = 1)
  # Not wrapped in suppressMessages(): 60g writes its log with message(), through a sink.
  # Its stop text also goes to stderr, which the test does not read.
  outcome <- tryCatch(source(file.path(scripts_folder, "60g_dada_seqtab.R"), local = new.env()),
                      error = function(e) e)
  # Put the console back, should the script have failed before closing its sinks.
  while (sink.number() > 0) sink()
  if (sink.number(type = "message") != 2) sink(type = "message")
  log_text <- if (file.exists(log_file)) paste(readLines(log_file, warn = FALSE), collapse = " ") else ""
  expected <- "[dada_seqtab] STOP: the filtered R2 file of sample(s) B is not among the inputs of this step."
  if (inherits(outcome, "error") && grepl(expected, log_text, fixed = TRUE) &&
      grepl(log_file, log_text, fixed = TRUE)) {
    report("seqtab_stops_without_partner", "OK")
  } else {
    report("seqtab_stops_without_partner",
           paste("outcome:", if (inherits(outcome, "error")) conditionMessage(outcome) else "no error",
                 "| log:", log_text))
  }
}


# ── 3. The flat-rounds line of the learning log (L2b) ──
# The round keeper of 60d with a made-up fit whose tables are set by hand: round 1
# changes the table by 4.92 (the four "read correctly" rows go from 1 to 0.97), every
# later round by 1e-4, which is below 1/1000 of round 1's. The rule then ends the
# learning in round 3 with a window of one round, in round 4 with a window of two.
rows <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
trans <- matrix(1000L, nrow = 16, ncol = 41, dimnames = list(rows, as.character(0:40)))
base_table <- matrix(0.01, nrow = 16, ncol = 41, dimnames = dimnames(trans))
base_table[c("A2A", "C2C", "G2G", "T2T"), ] <- 0.97
tables <- list(base_table, base_table)                      # rounds 0 and 1
for (round in 2:4) {
  next_table <- tables[[round]]
  next_table[round, 1] <- next_table[round, 1] + 1e-4       # a change of 1e-4
  tables[[round + 1]] <- next_table
}
flat_line <- function(max_flat_rounds) {
  calls <- new.env()
  calls$count <- 0
  error_function_registry[["handmade"]] <<- list(
    label = "handmade fit", suits = NA_character_,
    fit = function(trans, quality_values, extend_range = FALSE) {
      calls$count <- calls$count + 1
      list(table = tables[[calls$count]], anchors = NA_integer_, dropped = integer(0), note = "")
    })
  keeper <- make_round_keeper("R1", "handmade", integer(0), max_flat_rounds, max_consist = 10)
  lines <- character(0)
  withCallingHandlers(
    for (call in seq_len(max_flat_rounds + 3)) keeper$errfun(trans),
    message = function(m) {
      lines <<- c(lines, trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  list(line = grep("the flat-rounds rule ends the learning", lines, value = TRUE),
       round = keeper$state$flat_stop_round)
}
for (window in 1:2) {
  found <- flat_line(window)
  stop_round <- window + 2
  window_text <- if (window == 1) "over the last round (it is" else "over the last 2 rounds (it is"
  good <- length(found$line) == 1 && identical(as.integer(found$round), as.integer(stop_round)) &&
    grepl(window_text, found$line, fixed = TRUE) &&
    grepl(sprintf("so DADA2 ends the loop and reports convergence after %d rounds;", stop_round), found$line, fixed = TRUE) &&
    !grepl("Convergence after", found$line, fixed = TRUE) && !grepl("last 1 rounds", found$line, fixed = TRUE)
  report(sprintf("flat_line_window_of_%d", window),
         if (good) "OK" else paste("round", found$round, "|", paste(found$line, collapse = " / ")))
}


# ── 4. Quality scores outside the anchors of the binned fit ──
# A case that compares a text (or a few texts) with the one expected.
check_text <- function(case, got, expected) {
  report(case, if (identical(got, expected)) "OK" else paste("got:", paste(got, collapse = " | ")))
}
check_text("out_of_range_both_sides", out_of_range_phrase(c(7:19, 31:38), c(20, 30)),
           "21 quality scores, 13 below 20 (7-19) and 8 above 30 (31-38)")
check_text("out_of_range_below_only", out_of_range_phrase(7:19, c(20, 25, 30)),
           "13 quality scores below 20 (7-19)")
check_text("out_of_range_above_only", out_of_range_phrase(31:42, c(9, 20, 30)),
           "12 quality scores above 30 (31-42)")
check_text("out_of_range_single_score_on_one_side", out_of_range_phrase(c(7:16, 31), c(20, 30)),
           "11 quality scores, 10 below 20 (7-16) and 1 above 30 (31)")
check_text("out_of_range_ten_or_fewer_listed", out_of_range_phrase(c(7, 8, 35, 36, 37, 38), c(9, 20, 34)),
           "quality 7, 8, 35, 36, 37, 38")
# The stop of the binned fit itself, on counts with bases at quality scores 7 to 38 and
# the declared values 20 and 30 (both hold bases, so both are anchors).
counts <- trans
counts[, as.character(c(0:6, 39:40))] <- 0L
fit_stop <- tryCatch(error_function_registry$binned$fit(counts, c(20L, 30L)), error = conditionMessage)
report("binned_fit_stop_names_both_sides",
       if (grepl("have bases at 21 quality scores, 13 below 20 (7-19) and 8 above 30 (31-38), outside the quality values of the voting samples (20, 30)",
                 fit_stop, fixed = TRUE)) "OK" else paste("got:", fit_stop))
check_text("extended_note_names_both_sides", binned_extended_note(c(7:19, 31:38), c(20, 30), c(7, 38)),
           paste0("the counts hold bases at 21 quality scores, 13 below 20 (7-19) and 8 above 30 (31-38), ",
                  "outside the quality values of the voting samples; for this comparison the end column(s) ",
                  "7, 38 were added as anchors"))


# ── 5. The report's titles and MultiQC summary with many anchors ──
# A binned fit on 16 anchors: the list is shortened to "16 values (7-22)", so the title
# keeps the outcome and the model code, and stays within the about 133 characters the
# fit figure shows. Three or ten anchors are listed in full, as before.
binned_part <- function(anchors) {
  list(error_function = "binned", anchors_used = as.list(anchors), stopped_by = "cap",
       model_code = "R1:4be06c1d", rounds = 3, converged = FALSE, cycle_length = NA,
       matched_round = NA, flat_stop_round = NA, decided_by = "setting",
       change_per_round = list(30.1, 0.52, 0.04), bases_learned = 1e6, reads_learned = 4000,
       samples_learned = 1, warnings = list())
}
title_16 <- plot_title("R1", binned_part(7:22), 3)
expected_16 <- "R1 error model: binned fit on 16 values (7-22) - NOT converged: cap of 3 rounds (code R1:4be06c1d)"
report("title_shortens_more_than_ten_anchors",
       if (identical(title_16, expected_16) && nchar(title_16) <= 133) "OK" else paste("got:", title_16))
check_text("title_lists_three_anchors", plot_title("R1", binned_part(c(9, 20, 34)), 3),
           "R1 error model: binned fit on 9, 20, 34 - NOT converged: cap of 3 rounds (code R1:4be06c1d)")
check_text("title_lists_ten_anchors", plot_title("R1", binned_part(11:20), 3),
           "R1 error model: binned fit on 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 - NOT converged: cap of 3 rounds (code R1:4be06c1d)")
report_state <- new.env()
report_state$warning_texts <- character(0)
report_state$warning_directions <- character(0)
custom_function_file_input <- character(0)
loess_part <- binned_part(integer(0))
loess_part$error_function <- "loess"
loess_part$anchors_used <- NULL
summary_table <- summary_rows(list(max_consist = 3, R1 = binned_part(7:22), R2 = loess_part), NULL)
check_text("summary_shortens_more_than_ten_anchors", summary_table[["Error function"]],
           c("binned on 16 values (7-22)", "loess"))
check_text("goodness_of_fit_anchors_unchanged", c(anchors_text("9, 20, 34"), anchors_text(paste(7:38, collapse = ", "))),
           c("9, 20, 34", "32 values (7-38)"))


# ── 6. dada_filter writes the read counts and the quality values to two files ──
# 60b run whole, as Snakemake runs it, on a made-up pair of 60 reads of 120 bases whose
# quality values are binned: 100 bases at Q34, 15 at Q20 and 5 at Q9, the values of a
# NextSeq 1000/2000 with the 600-cycle kit. Every read passes the filter (about 0.8
# expected errors per read). filter_stats.json must keep the layout it had up to v2.6.1,
# byte for byte, because dada_seqtab reads it for the read tracking; quality_values.json
# must hold what choose_error_model reads: the sample, scan_limit, and per read direction
# the reads scanned, the values and the bases at each. The files stay in the folder given
# as the second argument, where tests/test_error_model_r.py runs choose_error_model on
# them. Needs the packages 60b loads (pacman, dada2, jsonlite and ShortRead).
filter_folder <- if (length(args) >= 2) args[2] else tempfile("dada_filter_")
if (!requireNamespace("pacman", quietly = TRUE) || !requireNamespace("dada2", quietly = TRUE) ||
    !requireNamespace("ShortRead", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
  report("dada_filter_writes_two_files", "SKIP: the R packages pacman, dada2, jsonlite and ShortRead are needed")
} else {
  # A gzipped FASTQ file of n_reads reads: the bases cycle through A, C, G and T from a
  # different start in each read; the quality letters are C (Q34), 5 (Q20) and * (Q9).
  write_fastq <- function(path, n_reads) {
    connection <- gzfile(path, "w")
    for (read in seq_len(n_reads)) {
      bases <- c("A", "C", "G", "T")[(seq_len(120) + read) %% 4 + 1]
      cat("@read", read, "\n", paste(bases, collapse = ""), "\n+\n",
          strrep("C", 100), strrep("5", 15), strrep("*", 5), "\n", sep = "", file = connection)
    }
    close(connection)
  }
  dir.create(file.path(filter_folder, "03.stripped"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(filter_folder, "logs", "dada_filter"), recursive = TRUE, showWarnings = FALSE)
  r1_in <- file.path(filter_folder, "03.stripped", "S1_R1.fastq.gz")
  r2_in <- file.path(filter_folder, "03.stripped", "S1_R2.fastq.gz")
  write_fastq(r1_in, 60)
  write_fastq(r2_in, 60)
  # truncLen 0 for both reads, as pick_trunclen writes it in manual mode.
  trunclen_file <- file.path(filter_folder, "trunclen.json")
  writeLines('{"r1": 0, "r2": 0, "mode": "manual"}', trunclen_file)
  stats_file <- file.path(filter_folder, "stats", "dada2", "S1.filter_stats.json")
  quality_file <- file.path(filter_folder, "stats", "dada2", "S1.quality_values.json")
  log_file <- file.path(filter_folder, "logs", "dada_filter", "S1.log")
  SnakemakeForFilter <- setClass("SnakemakeForFilter", representation(
    input = "list", output = "list", params = "list", log = "list", wildcards = "list",
    threads = "numeric"))
  # The rule's inputs, outputs and params as 60_dada2.smk declares them, with the
  # template's filter settings (max_ee 2 and 5, trunc_q 2, min_len 20, no max_len).
  snakemake <- SnakemakeForFilter(
    input = list(r1 = r1_in, r2 = r2_in, trunclen_json = trunclen_file, probe_json = list()),
    output = list(r1 = file.path(filter_folder, "04.filtered", "S1_R1_filt.fastq.gz"),
                  r2 = file.path(filter_folder, "04.filtered", "S1_R2_filt.fastq.gz"),
                  stats = stats_file, quality_values = quality_file),
    params = list(max_ee = list(2, 5), trunc_q = 2L, min_len = 20L, max_len = NULL,
                  amp_type = "16S", min_len_stat = NULL),
    log = list(log_file), wildcards = list(sample = "S1"), threads = 1)
  outcome <- tryCatch(source(file.path(scripts_folder, "60b_dada_filter.R"), local = new.env()),
                      error = function(e) e)
  # Put the console back, should the script have failed before closing its sinks.
  while (sink.number() > 0) sink()
  if (sink.number(type = "message") != 2) sink(type = "message")
  problems <- character(0)
  if (inherits(outcome, "error")) {
    problems <- c(problems, paste("the script failed:", conditionMessage(outcome)))
  }
  # filter_stats.json, byte for byte as v2.6.1 wrote it (write() adds the last line end).
  expected_stats <- '{\n  "sample": "S1",\n  "reads_in": 60,\n  "reads_out": 60\n}\n'
  stats_text <- if (file.exists(stats_file)) readChar(stats_file, file.size(stats_file), useBytes = TRUE) else "(missing)"
  if (!identical(stats_text, expected_stats)) {
    problems <- c(problems, paste("filter_stats.json is", stats_text))
  }
  if (!file.exists(quality_file)) {
    problems <- c(problems, "quality_values.json is missing")
  } else {
    quality <- jsonlite::fromJSON(quality_file, simplifyVector = FALSE)
    if (!identical(names(quality), c("sample", "scan_limit", "R1", "R2"))) {
      problems <- c(problems, paste("quality_values.json has the fields", paste(names(quality), collapse = ", ")))
    } else {
      for (direction in c("R1", "R2")) {
        part <- quality[[direction]]
        values <- as.numeric(unlist(part$values))
        bases <- unlist(part$bases_per_value)
        bases_found <- as.numeric(bases[c("9", "20", "34")])
        if (!identical(quality$sample, "S1") || !identical(as.numeric(part$reads_scanned), 60) ||
            !identical(values, c(9, 20, 34)) || !identical(bases_found, c(300, 900, 6000)) ||
            length(bases) != 3) {
          problems <- c(problems, paste(direction, "of quality_values.json:",
                                        jsonlite::toJSON(part, auto_unbox = TRUE)))
        }
      }
    }
  }
  report("dada_filter_writes_two_files", if (length(problems) == 0) "OK" else paste(problems, collapse = " | "))
}
