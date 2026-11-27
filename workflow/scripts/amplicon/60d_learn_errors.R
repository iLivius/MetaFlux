#!/usr/bin/env Rscript
# Learns DADA2's error model of the run, once for the R1 reads and once for the R2
# reads, and records how the learning ended.
#
# What the model is: for every quality score, how often each base is read as each
# other base (a true A read as C, as G, as T, and so on: 12 kinds of misread). DADA2
# needs it to tell a real sequence variant from a sequencing error of a more abundant
# sequence, so dada_seqtab denoises every sample with it. R1 and R2 get separate
# models: on Illumina the reverse read is usually noisier (its signal degrades further
# into the run), so one shared model would fit neither well.
#
# Input:
#   - every sample's filtered FASTQ pair from 04.filtered/ (rule dada_filter);
#   - stats/dada2/error_model_choice.json (rule choose_error_model): per read direction
#     the error function to learn with ("loess" for continuous quality values, "binned"
#     for binned ones, or the one set in amplicon.dada2.learn_errors.error_model) and
#     the quality values the voting samples' filtered reads carry, which are the
#     anchors the binned fit may use;
#   - 60_error_functions.R, the registry of the error functions, shared with the report
#     (60e_error_model_report.R), and vendor/dada2_makeBinnedQualErrfun.R, DADA2's
#     binned fit, which the registry calls. Both are sourced below; the rule lists them
#     as inputs, so a change to either makes this step run again.
#   - only under error_model: custom (EXPERIMENTAL): the R file of the user's own error
#     function (input custom_function_file; params custom_function_name and
#     custom_function_file_as_configured, from amplicon.dada2.learn_errors.
#     custom_error_function). It is read with the registry's
#     load_custom_error_function() and learned with like loess or binned, with every
#     table it returns checked first (check_custom_table() there). As an input, an edit
#     of the file makes this step run again.
#
# What it does, in order:
#   1. Sorts the R1 and the R2 files byte by byte (the same order whatever the
#      language setting of the shell; see the comment at the sort below), names them
#      after the samples, and reads the choice of error function for each direction.
#   2. Sets R's random-number seed (amplicon.seed), then calls learnErrors for R1 and
#      then for R2, in this one R session. The order matters for reproducibility:
#      learnErrors (randomize = TRUE) shuffles the list of files with R's random
#      numbers and reads whole samples in that order until it has nbases bases, and
#      R2's shuffle continues from where R1's left the random numbers. Nothing between
#      set.seed() and the R1 call draws a random number, so the same sorted files and
#      the same seed give the same samples, and the same models. They are the samples
#      dada_seqtab drew up to v2.6.1 wherever v2.6.1 sorted the files the same way:
#      under LANG/LC_ALL=C, or for sample names that sort alike in every language
#      setting (as those of the test sets do).
#   3. Inside learnErrors, DADA2 repeats two steps on those reads: denoise them with
#      the current error rates, then fit new rates to the misreads it counted (the
#      "transition counts": bases of each true nucleotide read as each nucleotide, at
#      each quality score), with the error function chosen for that direction. A
#      "round keeper" (make_round_keeper below) sits between DADA2 and the fit: it hands
#      the fit's table to DADA2 unchanged, notes how much the table changed in the
#      round, and applies MetaFlux's flat-rounds rule. When the fit cannot fit (for
#      example the binned fit finds bases at a quality value outside the anchors), the
#      job stops with the fit's own reason, which DADA2 would hide behind "Error rates
#      could not be estimated". The loop ends in one of three ways:
#        - dada2_repeat: a round gives back a table that already went into an earlier
#          round. This is DADA2's own stop; MetaFlux calls it "converged". The table
#          can repeat the one that went into the same round, or one of a few rounds
#          before (a cycle: a few reads keep switching between two ASVs, so the loop
#          goes round in a circle).
#        - flat_rounds: the model has stopped improving (the rule is described above
#          make_round_keeper). The keeper then hands DADA2 back the table that went
#          into the round, which DADA2 sees as a repeat, so it ends the loop normally
#          and keeps that table. MetaFlux counts such a model as NOT converged.
#        - cap: the last round allowed (amplicon.dada2.learn_errors.max_consist) gave
#          a table not seen before. NOT converged.
#   4. Checks that its own record of the rounds matches the tables DADA2 kept, works
#      out each model's code (see model_code_of below) and writes the outcome to the
#      log.
#
# Output:
#   - 05.dada2/error_model_R1.rds and _R2.rds: what learnErrors returned, unchanged
#     (the final table err_out, the table that went into each round err_in, and the
#     counts trans), read by error_model_report (the plots) and by dada_seqtab, which
#     denoises with them.
#   - stats/dada2/error_model.json, the record: under error_model: custom, which
#     function of the user's (its name, its file as configured and as resolved, the
#     md5 sum of the file, the R packages it loaded); per read direction the function
#     used and why (the kind of quality values, auto or the setting), the anchors of a
#     binned fit, how many bases, reads and samples learnErrors learned from, the
#     rounds, how the learning ended, the change in every round, the distance of every
#     round's starting table to the final one (DADA2's checkConvergence), the model
#     code, and the counts the final table was fitted to. Read by error_model_report,
#     which takes from it the counts behind each model (for its refits, the rates table
#     and the fit figures), the change and the distance of every round (the rounds
#     table and the convergence figure), the quality values and the anchors, the
#     warnings, and the bases, reads and samples learned from (the MultiQC summary);
#     and by check_error_model, which stops the run before denoising when a model did
#     not converge and the user has not accepted its code. The record holds no time
#     stamp, so the same model gives the same file.
#
# The script always finishes when DADA2 does, whatever the outcome, except when its own
# check of the rounds fails (step 4: the record of the rounds does not match DADA2's
# tables, stop text L7); it then stops without writing the models or the record, and the
# reason is in the log. Otherwise the check is a later step, so the models and the record
# stay on disk for the user to look at.

# ── Load packages and redirect output to the log file ───────────
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, jsonlite))

# sink() below redirects everything this script prints from here on — messages,
# warnings, DADA2's own progress lines — into the rule's log file instead of the
# console, so the log: file Snakemake tracks has the full record.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")


#' Stop the job with a text the user can read in full, in the log and on the terminal.
#'
#' R's own stop() prints at most about 1,000 characters of its message, so the full
#' text goes out with message() instead: once into the log (the sinks are still on),
#' then, after the sinks are closed, once more onto the terminal (on a cluster, into
#' the scheduler's log of the job). The short stop() line at the end makes the job
#' fail, which is what stops Snakemake.
#'
#' @param text The full message, starting with "[learn_errors] STOP:".
#' @return Never returns.
stop_with_message <- function(text) {
  message(text)                                       # full text into the log
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  message(text)                                       # full text onto the terminal
  stop("learn_errors stopped: reason above and in logs/learn_errors.log", call. = FALSE)
}


# ── The error functions MetaFlux can learn with ───────────────────
# DADA2 calls an error function once per round with the transition counts of that
# round (a 16 x Q matrix: rows A2A, A2C, ... T2T, columns the quality scores 0..Q-1)
# and expects back a table of the same shape with the probability of each transition
# at each quality score. The functions MetaFlux offers are listed in one file shared
# with the report, 60_error_functions.R (the "registry": loess, DADA2's default smooth
# curve, and binned, DADA2's makeBinnedQualErrfun, which joins the rates counted at the
# few values of binned quality scores with straight lines), so that the learning and
# the report always mean the same thing by a function's name. Each entry gives a label
# for the logs and a fit function (see that file's header). The binned entry calls the
# vendored copy of DADA2's function, which is therefore sourced first. Sourcing only
# defines functions: it reads no data and draws no random numbers.
source(snakemake@input[["vendored_binned"]])
source(snakemake@input[["error_functions"]])


# ── The flat-rounds rule: its two constants ───────────────────────
# The rule ends the learning in round r when BOTH hold:
#   (a) deep enough: the change in round r (the sum, over all 16 x Q cells, of how
#       much the table moved in that round) is at most 1/FLAT_DEPTH_FACTOR of the
#       change in round 1, so the model has already improved by about three orders
#       of magnitude (DADA2's author asks for checkConvergence values "several orders
#       of magnitude" below the first);
#   (b) flat: the lowest change seen so far fell by less than FLAT_MIN_GAIN (10 %)
#       over the last N rounds (N = amplicon.dada2.learn_errors.max_flat_rounds),
#       which also needs more than N rounds to judge.
# Why it is needed: on some libraries a few reads keep switching between two ASVs, so
# the table never repeats exactly but barely changes, for a hundred rounds or more;
# DADA2's author has said that ending the loop where the checkConvergence values stop
# decreasing, and using that table, makes no biological difference.
# Calibration: the rule was calibrated by replaying four learning runs round by round,
# all with the loess fit: R1 (two samples, 38 rounds) and R2 (one sample, 101 rounds)
# of one binned NextSeq 2000 16S run, a MiSeq COI library (32 rounds), and a COI
# library with reads of both strands in one file, which never converged. With N = 20
# the rule stopped the 101-round run at round 42, where denoising that sample with the
# stopped model changed 0.02 % of reads, and the mixed library at round 32. Over every
# N tested the largest change was 0.04 % of reads. The other two runs ended by DADA2's
# own repeat before the rule could act. The same replays with the binned fit, on six
# learning runs: N = 20 stopped one 95-round run at round 49, where denoising with the
# stopped model changed 0.007 % of reads (at most 0.02 % over every N tested), and no
# run reached the cap of 200.
# The rule alone cannot tell a usable model from one fitted to faulty reads (the mixed
# library looked just like a slow healthy one), so a model it ends counts as not
# converged and goes through check_error_model like a model stopped at the cap.
FLAT_DEPTH_FACTOR <- 1000
FLAT_MIN_GAIN     <- 0.10


# ── Helpers ───────────────────────────────────────────────────────

#' A number with three significant digits, for the log (27.7, 0.0824, 1.13e-05).
#'
#' @param x Numeric vector.
#' @return Character vector of the same length, without padding (formatC pads a zero
#'   to the width of three digits, hence trimws).
format_number <- function(x) {
  trimws(formatC(x, digits = 3, format = "g"))
}


#' Write one line of the round keeper into the log, next to DADA2's progress lines.
#'
#' DADA2 prints its "selfConsist step N" lines without a line end (it starts each with
#' one instead), so a line written while the loop runs must start with its own line
#' end, and must not end with one, or the log would show a blank line between rounds.
#'
#' @param text The line, without line ends.
#' @return Nothing.
log_round_line <- function(text) {
  message("\n", text, appendLF = FALSE)
}


#' The log line for quality values the fit left out (the binned fit only), such as
#' "[learn_errors] R2: quality value 2 was seen in the filtered reads but has no bases
#' in the reads learnErrors drew; the binned fit uses 9, 20, 34."
#'
#' A value is left out when it has no bases in the counts of the round, or, between the
#' lowest and the highest anchor, no bases of one of the four nucleotides (see
#' binned_anchors_with_data in 60_error_functions.R); the line says which, for each.
#' The keeper writes it when the set of left-out values changes, so it appears once in
#' a normal learning.
#'
#' @param direction "R1" or "R2".
#' @param fit What the fit returned (its dropped and anchors).
#' @param trans The counts of the round.
#' @param function_name The registry name of the fit, for the line.
#' @return The line.
dropped_values_text <- function(direction, fit, trans, function_name) {
  uses <- sprintf("the %s fit uses %s.", function_name, paste(fit$anchors, collapse = ", "))
  if (length(fit$dropped) == 0) {
    return(sprintf("[learn_errors] %s: every quality value seen in the filtered reads now has bases in the reads learnErrors drew; %s",
                   direction, uses))
  }
  quality_of_column <- as.integer(colnames(trans))
  columns_with_bases <- quality_of_column[colSums(trans) > 0]
  without_bases <- fit$dropped[!(fit$dropped %in% columns_with_bases)]
  without_one_nucleotide <- setdiff(fit$dropped, without_bases)
  reasons <- character(0)
  if (length(without_bases) == 1) {
    reasons <- c(reasons, sprintf("quality value %d was seen in the filtered reads but has no bases in the reads learnErrors drew",
                                  without_bases))
  } else if (length(without_bases) > 1) {
    reasons <- c(reasons, sprintf("quality values %s were seen in the filtered reads but have no bases in the reads learnErrors drew",
                                  paste(without_bases, collapse = ", ")))
  }
  if (length(without_one_nucleotide) == 1) {
    reasons <- c(reasons, sprintf("quality value %d has no bases of one of the four nucleotides there",
                                  without_one_nucleotide))
  } else if (length(without_one_nucleotide) > 1) {
    reasons <- c(reasons, sprintf("quality values %s have no bases of one of the four nucleotides there",
                                  paste(without_one_nucleotide, collapse = ", ")))
  }
  sprintf("[learn_errors] %s: %s; %s", direction, paste(reasons, collapse = ", and "), uses)
}


#' Check that a fitted table is one DADA2 can use, before handing it over.
#'
#' DADA2 itself only checks that the table has 16 rows and values from 0 to 1, and a
#' table with a missing value makes it fail later with "missing value where
#' TRUE/FALSE needed", which names no cause. This check names it.
#'
#' @param table The table the error function returned.
#' @param trans The counts it was fitted to.
#' @param function_name The registry name of the function, for the message.
#' @return Nothing; stops with "the <function> fit returned a table DADA2 cannot use
#'   (...)" when the table is not a numeric matrix of 16 rows with the row and column
#'   names of the counts and every value finite and between 0 and 1.
check_table <- function(table, trans, function_name) {
  problem <- NULL
  if (!is.matrix(table) || !is.numeric(table) || nrow(table) != 16) {
    problem <- "not 16 rows"
  } else if (!identical(rownames(table), rownames(trans)) ||
             !identical(colnames(table), colnames(trans))) {
    problem <- "columns differ from the counts"
  } else if (!all(is.finite(table)) || any(table < 0) || any(table > 1)) {
    problem <- "a missing, infinite or out-of-range value"
  }
  if (!is.null(problem)) {
    stop(sprintf("the %s fit returned a table DADA2 cannot use (%s)", function_name, problem),
         call. = FALSE)
  }
}


#' Make the round keeper of one read direction: the error function handed to DADA2.
#'
#' DADA2's learnErrors calls its error function once in a start-up pass (round 0, on
#' the counts of a first denoising with the highest possible error rates) and then
#' once per round, with that round's transition counts. The keeper wraps the chosen
#' fit and, each round:
#'   1. fits the table (the registry's fit function) and checks it (check_table);
#'   2. remembers the table that went into the round (the previous round's table,
#'      or, after round 0, the start-up table with the four "read correctly" rows set
#'      to 1, exactly as DADA2 sets them) and the change: the sum over all cells of
#'      |new table - table that went in|;
#'   3. if the new table repeats a table that already went into a round, DADA2 will
#'      stop by itself (converged): it hands the table over and does nothing else;
#'   4. otherwise, if the flat-rounds rule holds (see FLAT_DEPTH_FACTOR above), it
#'      hands back the table that went into this round instead. DADA2 finds it
#'      identical to that round's starting table and ends the loop, keeping it.
#'   5. otherwise it hands the new table over unchanged.
#' The keeper never changes a table it hands over (only the hand-back of step 4
#' differs from the fit), and draws no random numbers, so without a flat-rounds stop
#' the model is exactly the one learnErrors gives with the fit on its own.
#' Warnings of the fit are written to the log and the record once per distinct text,
#' with the round in which the text first appeared (DADA2 would hide them), and the
#' reason of an error is kept, because DADA2 replaces it with "Error rates could not be
#' estimated".
#'
#' @param direction "R1" or "R2", for the log.
#' @param function_name Registry name of the fit ("loess", "binned", or "custom" under
#'   error_model: custom).
#' @param quality_values The quality values the voting samples' filtered reads carry
#'   (error_model_choice.json), handed to the fit: the anchors the binned fit may use;
#'   the loess fit ignores them.
#' @param max_flat_rounds N of the flat-rounds rule; 0 switches the rule off.
#' @param max_consist The cap on rounds, to say in the log when a stop falls in the
#'   last allowed round.
#' @return list(errfun = the function to hand to learnErrors, state = an environment
#'   with what the keeper saw: tables_in (the table that went into each round),
#'   change and lowest (the change of each round and the lowest so far),
#'   flat_stop_round, last_new (the last table handed over unchanged),
#'   output_counts, output_anchors and output_dropped (the counts, the anchors and the
#'   left-out quality values behind the table DADA2 ends with), failure and
#'   record_mismatch (why the fit failed), fit_warnings).
make_round_keeper <- function(direction, function_name, quality_values,
                              max_flat_rounds, max_consist) {
  state <- new.env()
  state$calls           <- 0L
  state$tables_in       <- list()
  state$change          <- numeric(0)
  state$lowest          <- numeric(0)
  state$next_in         <- NULL
  state$last_new        <- NULL
  state$output_counts   <- NULL
  state$output_anchors  <- NA_integer_
  state$output_dropped  <- integer(0)
  state$flat_stop_round <- NULL
  state$failure         <- NULL
  state$record_mismatch <- FALSE
  state$last_dropped    <- integer(0)
  state$warning_texts   <- character(0)
  state$fit_warnings    <- character(0)
  state$last_call_at    <- Sys.time()

  fit_function <- error_function_registry[[function_name]]$fit

  # One round: what errfun below runs inside its error and warning handlers.
  keeper_step <- function(trans, round) {
    # After a flat-rounds stop DADA2 must not call again: the hand-back is identical to
    # a table it holds. A further call means MetaFlux's record of the tables differs
    # from DADA2's, so nothing the rule decided can be trusted.
    if (!is.null(state$flat_stop_round)) {
      state$record_mismatch <- TRUE
      stop(sprintf("DADA2 went on to round %d after the flat-rounds rule had ended the learning in round %d",
                   round, state$flat_stop_round), call. = FALSE)
    }

    fit <- fit_function(trans, quality_values)
    new_table <- fit$table
    check_table(new_table, trans, function_name)

    # Quality values the fit left out (the binned fit only; the loess fit leaves out
    # none). Said when the set of left-out values changes, so once in a normal run.
    if (!identical(fit$dropped, state$last_dropped)) {
      log_round_line(dropped_values_text(direction, fit, trans, function_name))
      state$last_dropped <- fit$dropped
    }

    seconds <- as.numeric(difftime(Sys.time(), state$last_call_at, units = "secs"))
    state$last_call_at <- Sys.time()

    # Round 0, DADA2's start-up pass: no table went in. DADA2 sets the four "read
    # correctly" rows (A2A, C2C, G2G, T2T) of this table to 1 before the next round;
    # the keeper does the same on its own copy, which is then the table that goes
    # into round 1. The table itself is handed over unchanged.
    if (round == 0) {
      start_table <- new_table
      start_table[c(1, 6, 11, 16), ] <- 1
      state$next_in        <- start_table
      state$output_counts  <- trans
      state$output_anchors <- fit$anchors
      state$output_dropped <- fit$dropped
      state$last_new       <- new_table
      return(new_table)
    }

    # Rounds 1, 2, ...: how much the table moved in this round, and the lowest change
    # so far.
    table_in <- state$next_in
    state$tables_in[[round]] <- table_in
    change <- sum(abs(new_table - table_in))
    state$change[round] <- change
    state$lowest[round] <- if (round == 1) change else min(state$lowest[round - 1], change)
    round_line <- if (round == 1) {
      sprintf("[learn_errors] %s round 1: change %s; %.0f s", direction, format_number(change), seconds)
    } else {
      sprintf("[learn_errors] %s round %d: change %s (round 1: %s; lowest so far %s); %.0f s",
              direction, round, format_number(change), format_number(state$change[1]),
              format_number(state$lowest[round]), seconds)
    }

    # The new table repeats one that went into a round (this one, or an earlier one
    # in a cycle): DADA2 stops by itself, converged. Hand it over unchanged.
    repeats_a_round <- any(vapply(state$tables_in, identical, logical(1), new_table))
    if (repeats_a_round) {
      state$output_counts  <- trans
      state$output_anchors <- fit$anchors
      state$output_dropped <- fit$dropped
      state$last_new       <- new_table
      log_round_line(round_line)
      return(new_table)
    }

    # The flat-rounds rule (see FLAT_DEPTH_FACTOR above). lowest[round - N] is the
    # lowest change of the rounds up to N rounds ago.
    if (max_flat_rounds > 0 && round > max_flat_rounds &&
        change <= state$change[1] / FLAT_DEPTH_FACTOR &&
        state$lowest[round] >= (1 - FLAT_MIN_GAIN) * state$lowest[round - max_flat_rounds]) {
      state$flat_stop_round <- round
      lowest_round <- which.min(state$change[seq_len(round)])
      # What DADA2 prints next. Its convergence line is described, not quoted: DADA2
      # 1.34.0 writes it with two spaces on each side of the number ("Convergence
      # after  8  rounds."), so a quote with single spaces would not be found by a
      # search of the log. "terminated before convergence" is part of DADA2's own line
      # ("Self-consistency loop terminated before convergence."), so it is quoted.
      dada2_next_line <- if (round >= max_consist) {
        "and DADA2's next line says \"terminated before convergence\", because this was also the last round allowed"
      } else {
        sprintf("so DADA2 ends the loop and reports convergence after %d rounds", round)
      }
      # The window of the rule: "the last round" for N = 1, "the last 20 rounds" otherwise.
      flat_window <- if (max_flat_rounds == 1) "the last round" else sprintf("the last %d rounds", max_flat_rounds)
      log_round_line(round_line)
      log_round_line(sprintf(paste0(
        "[learn_errors] %s round %d: the flat-rounds rule ends the learning: the change, %s, ",
        "is below 1/1000 of round 1's (%s), and the lowest change so far fell by less than ",
        "10 %% over %s (it is %s, from round %d). MetaFlux hands DADA2 back ",
        "the table that went into this round, %s; MetaFlux counts this model as not converged."),
        direction, round, format_number(change), format_number(state$change[1]),
        flat_window, format_number(state$lowest[round]), lowest_round, dada2_next_line))
      # output_counts, output_anchors and output_dropped stay those of the round
      # before: the counts (and anchors) the handed-back table was fitted to.
      return(table_in)
    }

    # A regular round: hand the new table over unchanged; it goes into the next round.
    state$next_in        <- new_table
    state$output_counts  <- trans
    state$output_anchors <- fit$anchors
    state$output_dropped <- fit$dropped
    state$last_new       <- new_table
    log_round_line(round_line)
    new_table
  }

  # What DADA2 calls. withCallingHandlers sees a warning of the fit before DADA2's
  # suppressWarnings() would hide it; the inner tryCatch keeps the reason of an error
  # (the first one only) and passes the error on, so DADA2 still fails as it would.
  errfun <- function(trans) {
    round <- state$calls                               # 0 = DADA2's start-up pass
    state$calls <- round + 1L
    withCallingHandlers(
      tryCatch(keeper_step(trans, round), error = function(e) {
        if (is.null(state$failure)) {
          state$failure <- sprintf("round %d: %s", round, conditionMessage(e))
        }
        stop(e)
      }),
      warning = function(w) {
        text <- conditionMessage(w)
        # Each distinct text once per read direction, with the round it first appeared
        # in. A function that gives several warnings in every round (one per kind of
        # misread, say, which a model fitted separately to each of the 12 kinds can
        # give) would otherwise add a copy of each text to the log and the record in
        # every round.
        if (!(text %in% state$warning_texts)) {
          log_round_line(sprintf("[learn_errors] %s round %d: the %s fit warned: %s",
                                 direction, round, function_name, text))
          state$warning_texts <- c(state$warning_texts, text)
          state$fit_warnings <- c(state$fit_warnings,
                                  sprintf("round %d: the %s fit warned: %s", round, function_name, text))
        }
        invokeRestart("muffleWarning")
      })
  }

  list(errfun = errfun, state = state)
}


#' Fill in the output folder and the read direction in a reason written by the
#' registry (60_error_functions.R), which knows neither and writes {out} and {dir}.
#'
#' @param text The text with the placeholders.
#' @param direction "R1" or "R2".
#' @return The text with {out} replaced by the run's output folder (absolute) and
#'   {dir} by the direction.
fill_placeholders <- function(text, direction) {
  text <- gsub("{out}", output_folder, text, fixed = TRUE)
  gsub("{dir}", direction, text, fixed = TRUE)
}


#' Run learnErrors on one read direction, and stop the job with the reason if it fails.
#'
#' @param files Sorted filtered FASTQ files of one direction, named after the samples.
#' @param keeper The round keeper of this direction (make_round_keeper).
#' @param direction "R1" or "R2".
#' @return What learnErrors returns: list(err_out, err_in, trans).
learn_one_direction <- function(files, keeper, direction) {
  message("[learn_errors] Learning error rates (", direction, ", nbases=", nbases, ")...")
  tryCatch(
    learnErrors(files, nbases = nbases, multithread = threads, randomize = TRUE,
                MAX_CONSIST = max_consist, errorEstimationFunction = keeper$errfun,
                verbose = TRUE),
    error = function(e) {
      if (isTRUE(keeper$state$record_mismatch)) {
        stop_with_message(record_mismatch_text(direction, keeper$state$flat_stop_round))
      }
      if (!is.null(keeper$state$failure)) {
        stop_with_message(fill_placeholders(sprintf(paste0(
          "[learn_errors] STOP: the %s error function failed in %s. DADA2 reports this only ",
          "as \"Error rates could not be estimated\", which is why MetaFlux repeats the ",
          "reason here."), direction, keeper$state$failure), direction))
      }
      stop_with_message(sprintf("[learn_errors] STOP: learnErrors failed on the %s reads: %s",
                                direction, conditionMessage(e)))
    })
}


#' The stop text (L7) for a record of the rounds that does not match DADA2's tables.
#'
#' @param direction "R1" or "R2".
#' @param round The first round that differs.
#' @return The text.
record_mismatch_text <- function(direction, round) {
  sprintf(paste0(
    "[learn_errors] STOP: MetaFlux's record of the %s rounds does not match DADA2's own ",
    "tables (round %d). The early-stop rule may have acted on wrong numbers, so the model ",
    "is not used. This should not happen with dada2 %s; please report it with ",
    "logs/learn_errors.log."), direction, round, as.character(packageVersion("dada2")))
}


#' How the learning of one direction ended, from DADA2's own object and the keeper.
#'
#' DADA2's own end line can be wrong: it says "terminated before convergence" whenever
#' the last allowed round was reached, also when the table did repeat in that round.
#' So the outcome is worked out here from the tables themselves:
#'   - the keeper ended it (flat_stop_round set): "flat_rounds", not converged;
#'   - the final table equals a table that went into a round: "dada2_repeat",
#'     converged; matched_round is that round, and the cycle runs from it to the last
#'     round (a cycle of 1 is a plain repeat);
#'   - otherwise "cap", not converged.
#' Before that, it checks that the keeper's record matches DADA2's tables, round by
#' round; the flat-rounds rule acted on the keeper's numbers, so a mismatch stops the
#' job (L7) rather than let a model through on wrong numbers.
#'
#' @param err What learnErrors returned.
#' @param state The keeper's state of this direction.
#' @param direction "R1" or "R2".
#' @return list(rounds, converged, stopped_by, matched_round, cycle_length,
#'   flat_stop_round, change_per_round, distance_to_final, counts, anchors, dropped).
describe_learning <- function(err, state, direction) {
  rounds <- length(err$err_in)

  # The keeper's record against DADA2's own tables.
  mismatch_round <- NA_integer_
  for (k in seq_len(rounds)) {
    if (k > length(state$tables_in) || !identical(state$tables_in[[k]], err$err_in[[k]])) {
      mismatch_round <- k
      break
    }
  }
  if (is.na(mismatch_round)) {
    final_table_expected <- if (!is.null(state$flat_stop_round)) state$tables_in[[rounds]] else state$last_new
    if (length(state$tables_in) != rounds || length(state$change) != rounds ||
        !identical(err$err_out, final_table_expected)) {
      mismatch_round <- rounds
    }
  }
  if (!is.na(mismatch_round)) {
    stop_with_message(record_mismatch_text(direction, mismatch_round))
  }

  matched <- which(vapply(err$err_in, identical, logical(1), err$err_out))
  matched_round   <- NA_integer_
  cycle_length    <- NA_integer_
  flat_stop_round <- NA_integer_
  if (!is.null(state$flat_stop_round)) {
    stopped_by      <- "flat_rounds"
    converged       <- FALSE
    flat_stop_round <- as.integer(state$flat_stop_round)
  } else if (length(matched) > 0) {
    stopped_by    <- "dada2_repeat"
    converged     <- TRUE
    matched_round <- as.integer(matched[1])
    cycle_length  <- as.integer(rounds - matched_round + 1)
  } else {
    stopped_by <- "cap"
    converged  <- FALSE
  }

  # The counts the final table was fitted to: those of the last round, or, after a
  # flat-rounds stop, those of the round before (the handed-back table was fitted to
  # them).
  counts <- if (stopped_by == "flat_rounds") state$output_counts else err$trans

  list(rounds            = rounds,
       converged         = converged,
       stopped_by        = stopped_by,
       matched_round     = matched_round,
       cycle_length      = cycle_length,
       flat_stop_round   = flat_stop_round,
       change_per_round  = state$change,
       # DADA2's checkConvergence: how far each round's starting table is from the
       # final one.
       distance_to_final = vapply(err$err_in, function(table_in) sum(abs(table_in - err$err_out)),
                                  numeric(1)),
       counts            = counts,
       # The anchors of the fit behind the final table, and the quality values it left
       # out (binned fit only; NA and none for the loess fit).
       anchors           = state$output_anchors,
       dropped           = state$output_dropped)
}


#' The code that names a model, such as "R2:0a55be94".
#'
#' The md5 sum of the final table's numbers, each written with 17 significant digits
#' (enough to tell any two different numbers apart), one per line, column by column.
#' Its first 8 characters, after the read direction, are the code the user lists in
#' accept_unconverged. It comes from the numbers, not from the bytes of the .rds file,
#' which change with the R version and the compression without the model changing.
#' R reads and writes numbers with a "." whatever the computer's language settings, so
#' the code does not depend on them.
#'
#' @param err_out The final table (16 x Q).
#' @param direction "R1" or "R2".
#' @return list(code = "R1:xxxxxxxx", md5 = the full md5 sum).
model_code_of <- function(err_out, direction) {
  numbers_file <- tempfile(fileext = ".txt")
  writeLines(sprintf("%.17g", err_out), numbers_file)
  md5 <- unname(tools::md5sum(numbers_file))
  unlink(numbers_file)
  list(code = paste0(direction, ":", substr(md5, 1, 8)), md5 = md5)
}


#' Write the outcome of one direction to the log (the L3, L4 and L5 lines).
#'
#' @param outcome From describe_learning.
#' @param code From model_code_of.
#' @param direction "R1" or "R2".
#' @return Nothing.
log_outcome <- function(outcome, code, direction) {
  last_round_allowed <- outcome$rounds >= max_consist
  if (outcome$stopped_by == "dada2_repeat") {
    if (outcome$cycle_length == 1) {
      text <- sprintf("[learn_errors] %s: converged: round %d gave back the table that went into it.",
                      direction, outcome$rounds)
    } else {
      text <- sprintf(paste0("[learn_errors] %s: converged: round %d gave back the table that went ",
                             "into round %d (a cycle of %d rounds)."),
                      direction, outcome$rounds, outcome$matched_round, outcome$cycle_length)
    }
    if (last_round_allowed) {
      text <- paste0(text, " DADA2's line above says \"terminated before convergence\" only ",
                     "because this was the last round allowed; the table did repeat.")
    }
  } else if (outcome$stopped_by == "flat_rounds") {
    text <- sprintf(paste0("[learn_errors] %s: NOT converged: stopped by the flat-rounds rule in ",
                           "round %d (amplicon.dada2.learn_errors.max_flat_rounds: %d)."),
                    direction, outcome$flat_stop_round, max_flat_rounds)
  } else {
    rounds_text <- if (max_consist == 1) "1 round" else paste(max_consist, "rounds")
    text <- sprintf(paste0("[learn_errors] %s: NOT converged: stopped at the cap of %s ",
                           "(amplicon.dada2.learn_errors.max_consist). DADA2's line above says ",
                           "\"terminated before convergence\"."),
                    direction, rounds_text)
  }
  message(text)
  message("[learn_errors] ", direction, ": change per round: ",
          paste(format_number(outcome$change_per_round), collapse = ", "))
  message("[learn_errors] ", direction, ": distance of each round's starting table to the final one: ",
          paste(format_number(outcome$distance_to_final), collapse = ", "))
  message("[learn_errors] ", direction, ": model code ", code$code)
}


#' How many bases, reads and samples learnErrors learned from in each direction, read
#' back from DADA2's own line in this job's log.
#'
#' learnErrors prints, before its first round, a line such as "33435402 total bases in
#' 217113 reads from 6 samples will be used for learning the error rates." (with cat(),
#' so a round number can come out as 1.2e+08). It returns these numbers nowhere else,
#' so the log is read back here: the first such line is R1's, the second R2's. A
#' direction without its line gets missing values (null in the record) and a warning.
#'
#' @param log_path This job's log file.
#' @return list(R1 = list(bases, reads, samples, warning), R2 = ...), the warning being
#'   NULL or a text for the record.
read_learned_amounts <- function(log_path) {
  flush(log_con)
  log_lines <- readLines(log_path, warn = FALSE)
  pattern <- "^([0-9.e+]+) total bases in ([0-9.e+]+) reads from ([0-9]+) samples"
  found <- regmatches(log_lines, regexec(pattern, log_lines))
  found <- found[lengths(found) == 4]
  amounts <- list()
  for (index in 1:2) {
    direction <- c("R1", "R2")[index]
    if (length(found) >= index) {
      amounts[[direction]] <- list(bases = as.numeric(found[[index]][2]),
                                   reads = as.numeric(found[[index]][3]),
                                   samples = as.numeric(found[[index]][4]),
                                   warning = NULL)
    } else {
      amounts[[direction]] <- list(bases = NA, reads = NA, samples = NA,
                                   warning = "the bases, reads and samples learnErrors learned from could not be read from its line in the log")
    }
  }
  amounts
}


#' The record of one direction, as it goes into stats/dada2/error_model.json.
#'
#' Vectors are wrapped in I() so that a vector of one value stays a JSON list; a
#' missing value is written as null.
#'
#' @param outcome From describe_learning.
#' @param code From model_code_of.
#' @param state The keeper's state (for the fit's warnings).
#' @param chosen This direction's part of error_model_choice.json (function, class,
#'   how it was decided, quality values).
#' @param learned This direction's part of read_learned_amounts().
#' @return A named list.
direction_record <- function(outcome, code, state, chosen, learned) {
  counts <- outcome$counts
  count_rows <- list()
  for (row_name in rownames(counts)) {
    count_rows[[row_name]] <- I(as.vector(counts[row_name, ]))
  }
  function_name <- chosen$error_function
  # The anchors exist for a fit that has them (binned); the loess fit has none, which
  # the record writes as null (anchors) and an empty list (dropped).
  has_anchors <- !all(is.na(outcome$anchors))
  anchors_used <- if (has_anchors) I(as.integer(outcome$anchors)) else NA
  anchors_dropped <- if (has_anchors) I(as.integer(outcome$dropped)) else I(integer(0))
  warnings <- state$fit_warnings
  if (!is.null(learned$warning)) {
    warnings <- c(warnings, learned$warning)
  }
  list(
    error_function       = function_name,
    error_function_label = error_function_registry[[function_name]]$label,
    quality_class        = chosen$quality_class,
    decided_by           = chosen$decided_by,
    quality_values       = I(as.integer(chosen$quality_values)),
    anchors_used         = anchors_used,
    anchors_dropped      = anchors_dropped,
    bases_learned        = learned$bases,
    reads_learned        = learned$reads,
    samples_learned      = learned$samples,
    rounds               = outcome$rounds,
    converged            = outcome$converged,
    stopped_by           = outcome$stopped_by,
    matched_round        = outcome$matched_round,
    cycle_length         = outcome$cycle_length,
    flat_stop_round      = outcome$flat_stop_round,
    change_per_round     = I(outcome$change_per_round),
    distance_to_final    = I(outcome$distance_to_final),
    depth_limit          = outcome$change_per_round[1] / FLAT_DEPTH_FACTOR,
    model_code           = code$code,
    model_md5            = code$md5,
    counts_behind_model  = list(quality_scores = I(as.integer(colnames(counts))),
                                rows = count_rows),
    warnings             = I(warnings)
  )
}


#' The record's description of the user's own error function (EXPERIMENTAL), or a
#' missing value (null in the record) when the models were learned with MetaFlux's own.
#'
#' What it holds, so that a model can be traced to the exact function behind it: the
#' function's name; the file as written in the config and as resolved at start-up (an
#' absolute path); the md5 sum of the file, which changes with any edit; and the R
#' packages that were loaded while the file was read or the function ran, with their
#' versions (packages MetaFlux had loaded before, dada2 and the packages it loads, are
#' not listed: the record's dada2_version and the environment cover them).
#'
#' @param custom What load_custom_error_function() returned, or NULL.
#' @return A named list, or NA.
custom_function_record <- function(custom) {
  if (is.null(custom)) {
    return(NA)
  }
  packages <- custom$state$packages_loaded
  versions <- list()
  for (package in packages) {
    versions[[package]] <- tryCatch(as.character(utils::packageVersion(package)),
                                    error = function(e) NA_character_)
  }
  # No package loaded: an empty JSON object ({}) rather than an empty list.
  if (length(versions) == 0) {
    versions <- structure(list(), names = character(0))
  }
  list(name               = custom$name,
       file_as_configured = custom$file_as_configured,
       file               = custom$file,
       md5                = custom$md5,
       packages_loaded    = versions)
}


#' The log line (L1) that says which error function a direction is learned with, and
#' why, such as "[learn_errors] R1: loess fit (loessErrfun, DADA2's default), chosen by
#' auto: 32 distinct quality values above Q2 (7-38)" or "[learn_errors] R2: binned fit
#' (makeBinnedQualErrfun from dada2 1.38.0) on quality values 9, 20, 34, chosen by auto".
#' The user's own function (experimental) is named with its file and the first 8
#' characters of the file's md5 sum, which tell two versions of the file apart:
#' "[learn_errors] R1: custom error function gam_errfun (experimental), as set in
#' amplicon.dada2.learn_errors.error_model, from /data/gam_errfun.R (md5 1a2b3c4d)".
#'
#' @param direction "R1" or "R2".
#' @param chosen This direction's part of error_model_choice.json.
#' @return The line.
plan_line <- function(direction, chosen) {
  label <- error_function_registry[[chosen$error_function]]$label
  how <- if (identical(chosen$decided_by, "auto")) {
    "chosen by auto"
  } else {
    "as set in amplicon.dada2.learn_errors.error_model"
  }
  if (identical(chosen$error_function, "custom")) {
    return(sprintf("[learn_errors] %s: %s, %s, from %s (md5 %s)", direction, label, how,
                   custom_function$file, substr(custom_function$md5, 1, 8)))
  }
  if (identical(chosen$error_function, "binned")) {
    return(sprintf("[learn_errors] %s: %s on quality values %s, %s", direction, label,
                   paste(chosen$quality_values, collapse = ", "), how))
  }
  above <- as.integer(chosen$values_above_q2)
  shown <- if (length(above) > 10) {
    sprintf("%d-%d", min(above), max(above))
  } else {
    paste(above, collapse = ", ")
  }
  sprintf("[learn_errors] %s: %s, %s: %d distinct quality values above Q2 (%s)", direction,
          label, how, length(above), shown)
}


# ── Read input file paths from Snakemake ─────────────────────────
# The order of the files decides which samples the models are learned from:
# learnErrors shuffles this list with the seeded random numbers and reads whole samples
# in the shuffled order until it has nbases bases. So the same list, in the same order,
# with the same seed, draws the same samples.
# method = "radix" sorts the paths byte by byte (by the codes of their characters, as
# in a C locale: every capital letter before every small letter, for example), the same
# in every language setting (locale) of the shell that runs MetaFlux. R's default sort()
# follows the alphabetical order of that setting instead, which differs for sample
# names that mix capitals and small letters, or "-", "_" and ".": under LANG/LC_ALL=C
# (common on cluster nodes, in containers and in cron jobs) three samples sort "Soil1,
# Soil_10, soil-2", under a UTF-8 setting "Soil_10, soil-2, Soil1". With the default
# sort, learning the same folder again from a shell with another setting could draw
# other samples, and so give another model, with another model code, and a second stop
# at check_error_model after the user had accepted the first code.
# Up to v2.6.1 dada_seqtab sorted with the default; the samples drawn are the same as
# then wherever v2.6.1 ran under LANG/LC_ALL=C, or wherever the names sort alike in
# every setting (as those of the test sets do).
# The two lists are sorted separately. Nothing in this script pairs an R1 file with an
# R2 file: each direction is learned on its own, from the samples its own shuffle
# draws, and the names below only label the files.
filt_r1 <- sort(unlist(snakemake@input[["r1"]]), method = "radix")
filt_r2 <- sort(unlist(snakemake@input[["r2"]]), method = "radix")

# Name every file after its sample, read from its own file name (dada_filter's naming:
# "{sample}_R1_filt.fastq.gz" and "{sample}_R2_filt.fastq.gz"). Each list is named from
# its own files because the two lists need not be in the same order: byte order puts
# samples called "Soil" and "Soil_R2" as "Soil, Soil_R2" in the R1 list but "Soil_R2,
# Soil" in the R2 list (a capital R comes before a small f). learnErrors does not use
# the names; the log line below lists the R1 names.
sample_names <- sub("_R1_filt\\.fastq\\.gz$", "", basename(filt_r1))
names(filt_r1) <- sample_names
names(filt_r2) <- sub("_R2_filt\\.fastq\\.gz$", "", basename(filt_r2))

message("[learn_errors] ", length(sample_names), " sample(s): ", paste(sample_names, collapse = ", "))

# ── Read parameters from Snakemake ───────────────────────────────
# nbases arrives as text when the config writes 1e8 (YAML reads a number without a
# decimal point as text), hence as.numeric(). max_consist and max_flat_rounds were
# checked at start-up (00_common.smk): whole numbers, at least 1 and at least 0.
nbases          <- as.numeric(snakemake@params[["nbases"]])
max_consist     <- as.integer(snakemake@params[["max_consist"]])
max_flat_rounds <- as.integer(snakemake@params[["max_flat_rounds"]])
seed            <- as.integer(snakemake@params[["seed"]])
threads         <- as.integer(snakemake@threads)

# The user's own error function (EXPERIMENTAL), handed over only under error_model:
# custom: its file (an input, as resolved at start-up: an absolute path), its name, and
# the path as written in the config. Under any other setting the input is empty (it
# arrives as NULL) and both params are "".
custom_function_file_input         <- snakemake@input[["custom_function_file"]]
custom_function_file               <- if (length(custom_function_file_input) > 0) custom_function_file_input[[1]] else ""
custom_function_name               <- snakemake@params[["custom_function_name"]]
custom_function_file_as_configured <- snakemake@params[["custom_function_file_as_configured"]]

# The run's output folder (the record is <out>/stats/dada2/error_model.json; OUT is an
# absolute path), for the paths in a stop message.
output_folder <- dirname(dirname(dirname(snakemake@output[["record"]])))

message("[learn_errors] at most ", max_consist, " rounds per read direction; flat-rounds rule: ",
        if (max_flat_rounds > 0) paste0("N = ", max_flat_rounds) else "off")

# ── Read the choice of error function ────────────────────────────
# Written by choose_error_model from the quality values of the filtered reads: per
# direction the function (a registry name), the kind of quality values, how it was
# decided (auto or the setting) and the quality values of the voting samples, which
# the binned fit takes as its anchors. Its warnings go into the record as well.
choice <- fromJSON(snakemake@input[["choice"]], simplifyVector = TRUE)

# ── The user's own error function (EXPERIMENTAL) ──────────────────
# Under error_model: custom the choice names "custom" for both directions. The function
# is then read from its file (load_custom_error_function in 60_error_functions.R: into
# an environment of its own, fetched by its name, the file's md5 sum noted) and added to
# the registry as the entry "custom", whose fit calls it and checks every table it
# returns before DADA2 sees it (check_custom_table there). This happens before the seed
# is set, and load_custom_error_function puts back R's kind of random numbers if the
# file switched it, so whatever reading the file does with random numbers cannot change
# the samples learnErrors draws. A file that cannot be read, a name it does not define,
# or a name that is not a function stops the job here, with the reason.
custom_function <- NULL
uses_custom <- identical(choice$R1$error_function, "custom") || identical(choice$R2$error_function, "custom")
if (uses_custom) {
  loaded <- tryCatch(load_custom_error_function(custom_function_file, custom_function_name,
                                                custom_function_file_as_configured),
                     error = function(e) e)
  if (inherits(loaded, "error")) {
    stop_with_message(sprintf(paste0(
      "[learn_errors] STOP: MetaFlux could not load your error function ",
      "(amplicon.dada2.learn_errors.custom_error_function, experimental): %s."),
      conditionMessage(loaded)))
  }
  custom_function <- loaded
  error_function_registry[["custom"]] <- custom_registry_entry(custom_function)
}

for (direction in c("R1", "R2")) {
  function_name <- choice[[direction]]$error_function
  if (is.null(error_function_registry[[function_name]])) {
    stop_with_message(sprintf(paste0(
      "[learn_errors] STOP: error_model_choice.json names the error function \"%s\" for %s, ",
      "which 60_error_functions.R does not define. The choice and the registry come from ",
      "different versions of MetaFlux; start MetaFlux again, so that choose_error_model ",
      "runs with this version."), function_name, direction))
  }
  message(plan_line(direction, choice[[direction]]))
}

# The two round keepers, made before the seed is set: nothing between set.seed() and
# the R1 call may draw a random number (the keepers draw none anyway).
keeper_r1 <- make_round_keeper("R1", choice$R1$error_function, as.integer(choice$R1$quality_values),
                               max_flat_rounds, max_consist)
keeper_r2 <- make_round_keeper("R2", choice$R2$error_function, as.integer(choice$R2$quality_values),
                               max_flat_rounds, max_consist)

# ── Learn the error models: seed, then R1, then R2 ───────────────
set.seed(seed)
err_r1 <- learn_one_direction(filt_r1, keeper_r1, "R1")
err_r2 <- learn_one_direction(filt_r2, keeper_r2, "R2")

# ── How each learning ended, the model codes, the log ────────────
outcome_r1 <- describe_learning(err_r1, keeper_r1$state, "R1")
outcome_r2 <- describe_learning(err_r2, keeper_r2$state, "R2")
code_r1 <- model_code_of(err_r1$err_out, "R1")
code_r2 <- model_code_of(err_r2$err_out, "R2")
log_outcome(outcome_r1, code_r1, "R1")
log_outcome(outcome_r2, code_r2, "R2")

# ── Save the models and the record ───────────────────────────────
# The models are saved exactly as learnErrors returned them, so dada() and
# plotErrors() take them as they took the objects in memory before.
dir.create(dirname(snakemake@output[["model_r1"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake@output[["record"]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(err_r1, snakemake@output[["model_r1"]])
saveRDS(err_r2, snakemake@output[["model_r2"]])

# The bases, reads and samples each learning used, from DADA2's own line in the log.
learned <- read_learned_amounts(snakemake@log[[1]])

# requested is the setting error_model (auto, loess, binned or custom), as
# choose_error_model received it; custom_error_function describes the user's own
# function when the models were learned with it (null otherwise); the run-level
# warnings are those of the choice (mixed or unusual quality values, a function set on
# the other kind of values, small samples).
record <- list(
  dada2_version         = as.character(packageVersion("dada2")),
  requested             = choice$requested,
  custom_error_function = custom_function_record(custom_function),
  nbases                = nbases,
  max_consist           = max_consist,
  seed                  = seed,
  max_flat_rounds       = max_flat_rounds,
  flat_rule             = list(depth_factor = FLAT_DEPTH_FACTOR, min_gain = FLAT_MIN_GAIN),
  R1                    = direction_record(outcome_r1, code_r1, keeper_r1$state, choice$R1, learned$R1),
  R2                    = direction_record(outcome_r2, code_r2, keeper_r2$state, choice$R2, learned$R2),
  warnings              = I(as.character(unlist(choice$warnings)))
)
write_json(record, snakemake@output[["record"]], auto_unbox = TRUE, pretty = TRUE,
           digits = NA, na = "null")
message("[learn_errors] DONE. Models and record written.")

# ── Close log ────────────────────────────────────────────────────
sink(type = "message")
sink(type = "output")
close(log_con)
