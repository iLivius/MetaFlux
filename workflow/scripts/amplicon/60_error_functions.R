# The error functions MetaFlux can learn DADA2's error model with: one list, shared by
# the learning step and the report.
#
# What an error function is: DADA2's learnErrors (rule learn_errors) learns the error
# model of a run in rounds. In each round it denoises the reads with the current error
# rates and counts, at each quality score, how often each true nucleotide was read as
# each nucleotide (the "transition counts": a 16 x Q matrix, rows A2A, A2C, ... T2T,
# columns the quality scores 0, 1, ..., Q-1). It then hands these counts to an error
# function, which turns them into the error rates of the next round: a table of the
# same shape with the probability of each transition at each quality score. How it
# does that matters most where the counts are thin or uneven, and the kind of quality
# values the reads carry decides which function fits:
#   - loess: DADA2's default (loessErrfun). For each of the 12 kinds of misread it fits
#     a smooth curve (a loess regression of the log10 error rate on the quality score,
#     weighted by the bases at each score). Made for continuous quality values (a MiSeq
#     writes every value from 2 to 38 or so). It ignores the quality values handed to it.
#   - binned: makeBinnedQualErrfun, written by DADA2's author for binned quality values
#     (a NovaSeq or NextSeq writes only a few, for example 9, 20 and 34). It joins the
#     error rates counted at those values with straight lines. MetaFlux carries a copy
#     of it (vendor/dada2_makeBinnedQualErrfun.R, unchanged, under DADA2's license),
#     because the function first appeared in a newer DADA2 than the one MetaFlux pins.
# Rule choose_error_model (60c_choose_error_model.py) picks one per read direction from
# the quality values it found in the filtered reads, or takes the one named by
# amplicon.dada2.learn_errors.error_model, and writes the choice to
# stats/dada2/error_model_choice.json.
#
# Who sources this file: 60d_learn_errors.R (rule learn_errors), which learns with the
# chosen function, and 60e_error_model_report.R (rule error_model_report), which draws
# the result. Both run in the dada2 environment (workflow/envs/dada2.yaml; or
# dada2_custom.yaml, the same builds plus the R packages mgcv and scam, when a run names
# an error function of the user's own), and both source
# vendor/dada2_makeBinnedQualErrfun.R first, because the binned entry calls it. MetaFlux
# scripts otherwise never share code (each carries its own copy of what it needs); this
# file is the one exception, so that the learning and the report can never disagree on
# what a function name means. Both rules list it as an input.
#
# The user's own error function (EXPERIMENTAL): a run config can name an R function of
# the user's in amplicon.dada2.learn_errors.custom_error_function (file and name). It is
# not an entry of the registry below, because it exists only in the runs that name it:
# the functions at the end of this file ("Your own error function") read it from its
# file, wrap it in an entry called "custom" and check every table it returns; the
# learning step adds that entry to its registry under error_model: custom, and the report
# whenever the function is named.
#
# What a change to this file costs: because both rules list it as an input, ANY change
# to it (a new function, a reworded label, a comment) makes learn_errors run again in
# every amplicon output folder after `git pull`, and with it the report, the check,
# dada_seqtab and every step after it, even when the function a run uses is unchanged.
# On a large data set that means hours of learning. So:
#   - A function that is only to be compared with the others in the report, not
#     learned with, does not belong here. It goes in a separate file sourced only by
#     60e_error_model_report.R and listed only as an input of that rule; adding it then
#     draws the report again without learning the models again. The check reads the
#     report, so the check, dada_seqtab and the steps after it still run again (see the
#     note on ancient() at the top of workflow/rules/amplicon/60_dada2.smk). (That file
#     is created when the first such function arrives.) A function of a user's own, for
#     their runs only, needs no such file: it is named in
#     amplicon.dada2.learn_errors.custom_error_function (see "Your own error function"
#     at the end of this file).
#   - Colours and line types of the functions in the report's figures are set in
#     60e_error_model_report.R, not here, so a colour change does not make the models be
#     learned again either (the same steps as above run again).
#
# How to add a function users can choose:
#   1. add an entry to error_function_registry below, between the two marker lines,
#      written like the two that are there (name = list(label, suits, fit));
#   2. add its name to ERROR_FUNCTION_NAMES in workflow/rules/shared/00_common.smk, in
#      the same order (the start-up check of error_model reads that list);
#   3. add it to SUITED_CLASS and ON_OTHER_CLASS in 60c_choose_error_model.py (which kind
#      of quality values it is made for, and what happens when it is set on the other
#      kind: "warn" and go on, or "stop");
#   4. write the tests and the documentation. tests/test_error_model_settings.py checks
#      that the three lists of names agree.
# A function that needs a package the dada2 environment does not have, such as the
# modified loess fits of the Ernakovich lab (loessErrfun_mod1 to mod4, which use dplyr),
# needs that package added to workflow/envs/dada2.yaml first, and to dada2_custom.yaml,
# which must stay dada2.yaml plus mgcv and scam (the report refits every function of
# this list there when a run names a function of its own). That is a new dependency
# (ask before adding one), and it changes the environment of every DADA2 step, so every
# amplicon step from the filtering on runs again in every folder. The alternative is a
# port to base R, shown to give identical tables.
#
# What every entry holds:
#   label : how the logs and the plot titles name the function;
#   suits : the kind of quality values it is made for, "continuous", "binned", or "both"
#           (60c_choose_error_model.py keeps the same fact in SUITED_CLASS; NA for the
#           user's own function, which is made for no particular kind);
#   fit   : function(trans, quality_values, extend_range = FALSE). trans is the matrix
#           of transition counts; quality_values are the quality values the voting
#           samples' filtered reads carry (stats/dada2/error_model_choice.json);
#           extend_range is FALSE when learning, and TRUE only when the report refits a
#           function that was not used, on the counts behind the model that was. It
#           returns list(table = the fitted 16 x Q table, anchors = the quality values
#           the fit is anchored on (NA when it has none), dropped = quality values it
#           left out, note = a sentence for the report, or ""). It stops with a
#           sentence that says why when it cannot fit; the learning step repeats that
#           sentence in its stop message.
# None of the code below draws random numbers, and no fit sets attributes on its
# table: DADA2 compares the tables of two rounds with identical(), so anything extra
# would change when the learning stops.


#' The quality values, as a short text for the messages: "9, 20, 34", or
#' "32 values, 7-38" when there are more than ten.
#'
#' @param values Integer vector of quality values.
#' @return Character string ("none" for an empty vector).
#' @seealso quality_scores_phrase, for the same values inside a sentence.
quality_values_text <- function(values) {
  values <- sort(unique(as.integer(values)))
  if (length(values) == 0) {
    return("none")
  }
  if (length(values) > 10) {
    return(sprintf("%d values, %d-%d", length(values), min(values), max(values)))
  }
  paste(values, collapse = ", ")
}


#' Quality scores inside a sentence: "quality 35, 36, 37", or "14 quality scores
#' (21-34)" when there are more than ten.
#'
#' @param values Integer vector of quality scores.
#' @return Character string.
quality_scores_phrase <- function(values) {
  values <- sort(unique(as.integer(values)))
  if (length(values) > 10) {
    return(sprintf("%d quality scores (%d-%d)", length(values), min(values), max(values)))
  }
  paste("quality", paste(values, collapse = ", "))
}


#' Which of the quality values the reads carry can anchor the binned fit, given the
#' counts of one round.
#'
#' makeBinnedQualErrfun draws a straight line between two neighbouring anchors only
#' when the counts hold bases at both of them. An anchor without bases makes it skip the
#' lines on both sides, and what the skipped stretch then gets depends on where it lies:
#' between two stretches that were fitted, missing rates, on which DADA2 fails with
#' "missing value where TRUE/FALSE needed"; at the low or high end of the fitted range,
#' a flat copy of the nearest fitted rate, which is a wrong model without any message
#' (for anchors 9, 20, 25, 34 with no bases at 25, every score above 20 gets the rate
#' of 20); and when every line of one nucleotide is skipped, the function itself
#' fails. learnErrors reads only some whole samples (until it has nbases bases), so a
#' quality value seen in the filtered reads of all samples can be missing from the reads
#' it drew. Hence, in order:
#'   1. sort the quality values; a value above the last column of the counts has no
#'      bases (the counts have one column per quality score from 0 up to the highest
#'      score in the reads learnErrors drew, at least 41 columns);
#'   2. keep a value whose column holds bases;
#'   3. drop a kept value between the lowest and the highest kept one at which one of the
#'      four true nucleotides (A, C, G or T) has no bases: its error rates would be 0/0,
#'      and the lines of that nucleotide on both sides would be skipped, as above. At the
#'      lowest or highest value this is tolerated (that nucleotide's rate is then copied
#'      flat from the next anchor); with values read from the reads it would need a whole
#'      nucleotide missing at that value, which no data checked showed;
#'   4. list the columns that hold bases below the lowest or above the highest kept
#'      value: there the function has no anchor, and it would stop with "Input data
#'      contains a lower (higher) quality score than the provided binned values".
#'
#' @param trans Transition counts of one round (16 rows, columns "0", "1", ...).
#' @param declared_values The quality values the voting samples' filtered reads carry,
#'   in any order.
#' @return list(anchors = the kept values, ascending; dropped = the values left out by
#'   steps 1-3, ascending; out_of_range = the columns of step 4, ascending; empty when
#'   no value was kept).
binned_anchors_with_data <- function(trans, declared_values) {
  declared <- sort(unique(as.integer(declared_values)))
  quality_of_column <- as.integer(colnames(trans))          # 0, 1, 2, ...
  bases_per_column <- colSums(trans)
  columns_with_bases <- quality_of_column[bases_per_column > 0]

  # Steps 1 and 2: a value keeps its place when its column holds bases (a value beyond
  # the last column is not in columns_with_bases, so it has none).
  has_bases <- declared %in% columns_with_bases
  anchors <- declared[has_bases]
  dropped <- declared[!has_bases]

  # Step 3: values between the lowest and the highest anchor need bases of every true
  # nucleotide. The rows of each nucleotide are its four transitions (A2A, A2C, A2G, A2T
  # for a true A, and so on).
  if (length(anchors) > 2) {
    inner_values <- anchors[-c(1, length(anchors))]
    thin_values <- integer(0)
    for (value in inner_values) {
      column <- which(quality_of_column == value)
      for (nucleotide in c("A", "C", "G", "T")) {
        rows_of_nucleotide <- paste0(nucleotide, "2", c("A", "C", "G", "T"))
        if (sum(trans[rows_of_nucleotide, column]) == 0) {
          thin_values <- c(thin_values, value)
          break
        }
      }
    }
    anchors <- setdiff(anchors, thin_values)
    dropped <- sort(c(dropped, thin_values))
  }

  # Step 4: bases outside the range of the anchors.
  if (length(anchors) == 0) {
    out_of_range <- integer(0)
  } else {
    outside <- columns_with_bases < min(anchors) | columns_with_bases > max(anchors)
    out_of_range <- columns_with_bases[outside]
  }

  list(anchors = as.integer(anchors), dropped = as.integer(dropped),
       out_of_range = as.integer(out_of_range))
}


#' Quality scores that lie outside the anchors of the binned fit, inside a sentence:
#' "quality 7, 8, 35, 36" when there are ten or fewer; when there are more, counted on
#' each side of the anchors, "21 quality scores, 13 below 20 (7-19) and 8 above 30
#' (31-38)", or "13 quality scores below 20 (7-19)" when all lie on one side.
#'
#' Why the two sides: quality_scores_phrase would write "21 quality scores (7-38)", a
#' range that spans the anchors themselves and reads as if every score from 7 to 38 lay
#' outside them. The scores outside the anchors are, by definition, below the lowest
#' anchor or above the highest (binned_anchors_with_data, step 4), so splitting at those
#' two anchors names every one of them. The split uses the anchors the fit actually uses,
#' not the declared quality values: a declared value without bases is no anchor.
#'
#' @param out_of_range The quality scores with bases outside the anchors.
#' @param anchors The anchors the fit uses (ascending, or in any order).
#' @return Character string.
out_of_range_phrase <- function(out_of_range, anchors) {
  values <- sort(unique(as.integer(out_of_range)))
  if (length(values) <= 10 || length(anchors) == 0) {
    return(quality_scores_phrase(values))
  }
  lowest_anchor  <- min(anchors)
  highest_anchor <- max(anchors)
  below <- values[values < lowest_anchor]
  above <- values[values > highest_anchor]
  # One side as "13 below 20 (7-19)", or "1 above 30 (31)" for a single score.
  side_text <- function(scores, where, anchor) {
    span <- if (length(scores) == 1) as.character(scores) else sprintf("%d-%d", min(scores), max(scores))
    sprintf("%d %s %d (%s)", length(scores), where, anchor, span)
  }
  if (length(below) > 0 && length(above) > 0) {
    return(sprintf("%d quality scores, %s and %s", length(values),
                   side_text(below, "below", lowest_anchor),
                   side_text(above, "above", highest_anchor)))
  }
  if (length(below) > 0) {
    return(sprintf("%d quality scores below %d (%d-%d)", length(below), lowest_anchor,
                   min(below), max(below)))
  }
  sprintf("%d quality scores above %d (%d-%d)", length(above), highest_anchor,
          min(above), max(above))
}


#' The reason the binned fit stops when the counts hold bases outside the quality values
#' of the voting samples (learning only; the learning step repeats it in its stop
#' message).
#'
#' The output folder and the read direction are not known here, so the text holds the
#' placeholders {out} and {dir}; 60d_learn_errors.R fills them in before printing.
#'
#' @param out_of_range The quality scores with bases outside the anchors.
#' @param declared_values The quality values of the voting samples.
#' @param anchors The anchors the fit uses (those of declared_values that hold bases),
#'   for out_of_range_phrase.
#' @return The sentence.
binned_out_of_range_text <- function(out_of_range, declared_values, anchors) {
  sprintf(paste0(
    "the counts of the reads learnErrors drew, which include mean qualities of ",
    "identical reads, have bases at %s, outside the quality values of the ",
    "voting samples (%s), so the binned fit has no anchor there. A sample with fewer ",
    "than 1,000 filtered reads (it does not vote) or with rewritten quality values is ",
    "the likely source: {out}/stats/dada2/quality_values.tsv lists the values of every ",
    "sample, and {out}/stats/falco/<sample>_{dir}_filtered/fastqc_report.html shows ",
    "each one. Leave that sample out, or set ",
    "amplicon.dada2.learn_errors.error_model: loess"),
    out_of_range_phrase(out_of_range, anchors), quality_values_text(declared_values))
}


#' The reason the binned fit stops when fewer than two quality values hold bases.
#'
#' The output folder is not known here, so the text holds the placeholder {out};
#' 60d_learn_errors.R fills it in before printing.
#'
#' @param declared_values The quality values of the voting samples.
#' @param anchors The values that hold bases (fewer than two).
#' @return The sentence.
binned_too_few_anchors_text <- function(declared_values, anchors) {
  kept_text <- if (length(anchors) == 0) {
    "none of them does"
  } else {
    sprintf("only %s does", quality_values_text(anchors))
  }
  sprintf(paste0(
    "the binned fit needs at least two of the quality values seen in the filtered reads ",
    "(%s) to carry bases in the reads learnErrors drew; %s. Set ",
    "amplicon.dada2.learn_errors.error_model: loess, or check ",
    "{out}/stats/dada2/quality_values.tsv"),
    quality_values_text(declared_values), kept_text)
}


#' The report's note for a comparison refit of the binned fit that needed extra anchors.
#'
#' @param out_of_range The quality scores with bases outside the anchors.
#' @param anchors The anchors before the end columns were added, for
#'   out_of_range_phrase.
#' @param added The end columns added as anchors.
#' @return The sentence.
binned_extended_note <- function(out_of_range, anchors, added) {
  sprintf(paste0(
    "the counts hold bases at %s, outside the quality values of the voting ",
    "samples; for this comparison the end column(s) %s were added as anchors"),
    out_of_range_phrase(out_of_range, anchors), paste(added, collapse = ", "))
}


# ── The registry ──────────────────────────────────────────────────
# One entry per function users can choose, looked up by its name (the value of
# amplicon.dada2.learn_errors.error_model, or the choice of auto). The names between the
# two marker lines are read by tests/test_error_model_settings.py, which checks them
# against ERROR_FUNCTION_NAMES in 00_common.smk and against SUITED_CLASS and
# ON_OTHER_CLASS in 60c_choose_error_model.py: keep one entry per name, written as
# `  name = list(` at the start of a line, in the order of ERROR_FUNCTION_NAMES (the
# report draws the functions in this order).
error_function_registry <- list(
  # ---- registry: start ----
  loess = list(
    label = "loess fit (loessErrfun, DADA2's default)",
    suits = "continuous",
    # DADA2's own function, unchanged: the model is exactly the one learnErrors gives
    # without an error function named. It needs no anchors and ignores the quality
    # values and extend_range.
    fit   = function(trans, quality_values, extend_range = FALSE) {
      list(table = dada2::loessErrfun(trans), anchors = NA_integer_,
           dropped = integer(0), note = "")
    }
  ),
  binned = list(
    label = "binned fit (makeBinnedQualErrfun from dada2 1.38.0)",
    suits = "binned",
    # The vendored makeBinnedQualErrfun, anchored on the quality values of the voting
    # samples that hold bases in these counts (binned_anchors_with_data above).
    fit   = function(trans, quality_values, extend_range = FALSE) {
      anchor_check <- binned_anchors_with_data(trans, quality_values)
      anchors <- anchor_check$anchors
      out_of_range <- anchor_check$out_of_range
      note <- ""
      # A refit for the report's comparison only (never when learning, and never for
      # the function the model was learned with, so that its refit reproduces the
      # model): instead of no curve at all, the lowest and the highest column with
      # bases become anchors where they lie outside the quality values, and the note
      # says so.
      if (length(out_of_range) > 0 && extend_range) {
        below <- out_of_range[out_of_range < min(anchors)]
        above <- out_of_range[out_of_range > max(anchors)]
        added <- c(if (length(below) > 0) min(below), if (length(above) > 0) max(above))
        # The note names the scores outside the anchors as they were before the end
        # columns were added, so it is written first.
        note <- binned_extended_note(out_of_range, anchors, added)
        anchors <- sort(unique(c(anchors, added)))
      }
      # A straight line needs two ends.
      if (length(anchors) < 2) {
        stop(binned_too_few_anchors_text(quality_values, anchors), call. = FALSE)
      }
      # When learning, bases outside the quality values of the voting samples stop the
      # fit with the reason (makeBinnedQualErrfun itself would stop with a sentence that
      # names neither the values nor the likely cause).
      if (length(out_of_range) > 0 && !extend_range) {
        stop(binned_out_of_range_text(out_of_range, quality_values, anchors), call. = FALSE)
      }
      list(table = makeBinnedQualErrfun(anchors)(trans), anchors = as.integer(anchors),
           dropped = anchor_check$dropped, note = note)
    }
  ),
  scam = list(
    label = "monotone decreasing P-spline (SCAM)",
    suits = "both",
    fit = function(trans, quality_values, extend_range = FALSE) {
      if (!requireNamespace("scam", quietly = TRUE)) {
        stop("The r-scam/SCAM package is not available in the DADA2 Conda environment",
             call. = FALSE)
      }
      suppressPackageStartupMessages(library(scam, character.only = TRUE))
      nts <- c("A", "C", "G", "T")
      qq <- as.numeric(colnames(trans))
      est <- matrix(0, nrow = 0, ncol = length(qq))

      for (nti in nts) {
        for (ntj in nts) {
          if (nti == ntj) next
          errs <- trans[paste0(nti, "2", ntj), ]
          tot <- colSums(trans[paste0(nti, "2", nts), ])
          rlogp <- log10((errs + 1) / tot)
          rlogp[is.infinite(rlogp)] <- NA
          df <- data.frame(q = qq, tot = tot, rlogp = rlogp)
          df$w <- log10(df$tot)
          df <- df[is.finite(df$rlogp) & is.finite(df$w) & df$w > 0, ]
          df$w <- df$w / mean(df$w)
          nq <- length(unique(df$q))
          if (nq < 4) {
            stop("Only ", nq, " Q values with data for ", nti, "2", ntj,
                 call. = FALSE)
          }
          mod <- scam::scam(
            rlogp ~ s(q, bs = "mpd", k = min(6, nq)),
            data = df,
            weights = w
          )
          qc <- pmin(pmax(qq, min(df$q)), max(df$q))
          est <- rbind(est, 10^as.numeric(predict(mod, newdata = data.frame(q = qc))))
        }
      }

      est[est > 0.25] <- 0.25
      est[est < 1e-7] <- 1e-7
      err <- rbind(
        1 - colSums(est[1:3, ]), est[1:3, ], est[4, ],
        1 - colSums(est[4:6, ]), est[5:6, ], est[7:8, ],
        1 - colSums(est[7:9, ]), est[9, ], est[10:12, ],
        1 - colSums(est[10:12, ])
      )
      rownames(err) <- paste0(rep(nts, each = 4), "2", nts)
      colnames(err) <- colnames(trans)
      list(table = err, anchors = NA_integer_, dropped = integer(0), note = "")
    }
  )
  # ---- registry: end ----
)


# ── Your own error function (EXPERIMENTAL) ─────────────────────────
# A run config can name an R function of the user's in
# amplicon.dada2.learn_errors.custom_error_function: the R file that defines it (file)
# and its name there (name). 00_common.smk checks at start-up that both are set, that
# the file can be read and that the name is one R can give a function; the rest can
# only be checked in R, here.
#
# The contract is DADA2's own errorEstimationFunction, so the same function works in
# plain DADA2 (learnErrors(files, errorEstimationFunction = the_function)): it is called
# with one argument, trans, DADA2's matrix of transition counts of one round (whole
# numbers; 16 rows named A2A, A2C, A2G, A2T, C2A, ... T2T, a true base and the base it
# was read as; one column per quality score, named "0", "1", ..., at least 41 of them
# and more when the reads carry scores above 40), and it returns a numeric matrix of the
# same shape: in each cell, the probability that the true base of the row is read as
# the second base of the row at that quality score.
#
# How it is used: 60d_learn_errors.R (under error_model: custom) and
# 60e_error_model_report.R (whenever the function is named) read it with
# load_custom_error_function() and add custom_registry_entry() to their registry under
# the name "custom", next to loess and binned. Its fit calls the function with
# call_custom_error_function() and checks the table with check_custom_table() before
# anything else sees it; in learning, a problem stops the job with the reason. In the
# report it leaves the function without a curve, with the reason in error_model_fit.tsv
# and in the figure caption; it is also a warning of the report when the models were
# learned with this function (its refit then does not give back the model), or when its
# file could not be read.


#' The label of the user's own function in the logs and the record, such as
#' "custom error function gam_errfun (experimental)".
#'
#' @param name The function's name.
#' @return The label.
custom_function_label <- function(name) {
  sprintf("custom error function %s (experimental)", name)
}


#' Read the user's own error function from its file.
#'
#' The file is read with source() into an R environment of its own (a separate set of
#' named objects in the R session; nothing to do with a conda environment), made for
#' this file. Its parent is R's global environment, where the MetaFlux script that reads
#' the file keeps its own variables and functions. What the file defines with an
#' ordinary assignment (<-), the function and any helpers or constants it needs, stays
#' in that environment, so an ordinary assignment cannot replace a function or a
#' variable of the MetaFlux script. Two things still reach the script, and the
#' documentation asks users to avoid both: <<- or assign(..., envir = globalenv()),
#' which write into the global environment, and variables the function uses without
#' defining them, which R then looks up in the global environment (the learning step and
#' the report define different ones there, so such a function can behave differently in
#' the two). The file sees R's own packages, those MetaFlux has loaded (dada2 among
#' them) and those it loads itself with library().
#' Reading the file draws no random number of MetaFlux's: the learning step reads it
#' before it sets R's seed, and R's kind of random numbers is put back after reading,
#' in case the file switched it (see the comment at the source() call).
#'
#' @param file The R file, as an absolute path (resolved at start-up).
#' @param name The function's name in that file.
#' @param file_as_configured The path as written in the config, for the record.
#' @return list(name, file, file_as_configured, md5 = the md5 sum of the file, fn = the
#'   function, state = an environment whose packages_loaded holds, sorted, the R
#'   packages loaded while the file was read or the function ran: those a function needs
#'   beyond what MetaFlux had loaded already). Stops with a sentence (the calling script
#'   puts its own "STOP:" before it) when no file or name was handed over, when the file
#'   cannot be read, when reading it fails with an R error, when it defines nothing of
#'   that name, and when what it defines under that name is not a function.
load_custom_error_function <- function(file, name, file_as_configured = file) {
  if (length(file) != 1 || length(name) != 1 || !nzchar(file) || !nzchar(name)) {
    stop("no file or no name was handed over for the custom error function ",
         "(amplicon.dada2.learn_errors.custom_error_function)", call. = FALSE)
  }
  md5 <- unname(tools::md5sum(file))
  if (is.na(md5)) {
    stop(sprintf("the file of the custom error function, %s, cannot be read", file), call. = FALSE)
  }

  # Read the file into its own environment, noting the packages that reading it loads.
  state <- new.env()
  namespaces_before <- loadedNamespaces()
  file_environment <- new.env(parent = globalenv())
  # R's kind of random numbers (RNGkind: which generator R uses) is noted before reading
  # and put back after it. A file may switch it: code written for parallel work often
  # starts with RNGkind("L'Ecuyer-CMRG"), and a package loaded with library() can do the
  # same. The seed the learning step sets afterwards would then start another series of
  # numbers, so learnErrors would draw other samples than with loess or binned, and the
  # models would no longer rest on the same reads. The kind is only touched when the
  # file changed it.
  rng_kind_before <- RNGkind()
  sourced <- tryCatch(source(file, local = file_environment), error = function(e) e)
  if (!identical(RNGkind(), rng_kind_before)) {
    RNGkind(kind = rng_kind_before[1], normal.kind = rng_kind_before[2],
            sample.kind = rng_kind_before[3])
  }
  if (inherits(sourced, "error")) {
    stop(sprintf("reading %s, the file of the custom error function, failed with this R error: %s",
                 file, conditionMessage(sourced)), call. = FALSE)
  }
  state$packages_loaded <- sort(setdiff(loadedNamespaces(), namespaces_before))

  # Fetch the function by its name, from what the file itself defines (inherits = FALSE:
  # a name the file does not define is not looked up in R's packages, where name: loess
  # would find a function that has nothing to do with the file).
  if (!exists(name, envir = file_environment, inherits = FALSE)) {
    # all.names = TRUE lists names that start with a dot too (.my_errfun), which the
    # start-up check accepts as a function name.
    defined <- ls(file_environment, all.names = TRUE)
    defined_text <- if (length(defined) == 0) "nothing" else paste(defined, collapse = ", ")
    stop(sprintf("%s does not define %s (it defines: %s). Check the name in the config against the file",
                 file, name, defined_text), call. = FALSE)
  }
  fn <- get(name, envir = file_environment, inherits = FALSE)
  if (!is.function(fn)) {
    stop(sprintf("%s in %s is not a function but an object of class %s", name, file, class(fn)[1]),
         call. = FALSE)
  }
  list(name = name, file = file, file_as_configured = file_as_configured, md5 = md5, fn = fn,
       state = state)
}


#' Call the user's own error function on one round's transition counts.
#'
#' Three things happen around the call:
#'   - R's random-number state is saved before the call and put back after it. The only
#'     random numbers of the learning step are learnErrors' own: it draws the order in
#'     which it reads the samples, and R2's draw continues from where R1's left the
#'     random numbers. A function that drew random numbers, or set a seed, would change
#'     that state, and with it which samples the R2 model is learned from. With the state
#'     put back, the samples drawn are the same whatever the error function, as with
#'     loess and binned, so models learned with different functions rest on the same
#'     reads. This does not replace a seed of the function's own: the report refits the
#'     function outside the learning, and a function that draws random numbers without
#'     setting its own seed then returns another table there, so its refit does not give
#'     back the model and the report warns.
#'   - An R error inside the function becomes an error that names the function, its file
#'     and R's own message. The learning step repeats it in its stop message, because
#'     DADA2 would only say "Error rates could not be estimated".
#'   - The R packages loaded during the call (a function may load mgcv, say, the first
#'     time it runs) are added to custom$state$packages_loaded, for the record.
#' Warnings of the function pass through unchanged: the learning step writes each
#' distinct text once to its log and to the record, as for the other functions.
#'
#' @param custom What load_custom_error_function() returned.
#' @param trans The transition counts of the round.
#' @return What the function returned, not yet checked (check_custom_table does that).
call_custom_error_function <- function(custom, trans) {
  namespaces_before <- loadedNamespaces()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) {
    saved_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  # Runs when the call ends, whether the function returned or failed.
  on.exit({
    if (had_seed) {
      assign(".Random.seed", saved_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())   # no state before the call: none after it
    }
    loaded_now <- setdiff(loadedNamespaces(), namespaces_before)
    custom$state$packages_loaded <- sort(unique(c(custom$state$packages_loaded, loaded_now)))
  })
  tryCatch(custom$fn(trans), error = function(e) {
    stop(sprintf("the custom error function %s (%s) stopped with this R error: %s",
                 custom$name, custom$file, conditionMessage(e)), call. = FALSE)
  })
}


#' Where the first marked cell of a table is, for a message: "A2G at quality 37".
#'
#' @param table A matrix with row names A2A ... T2T and column names "0", "1", ...
#' @param marked A logical matrix of the same shape, with at least one TRUE.
#' @return The text.
first_cell_text <- function(table, marked) {
  position <- which(marked, arr.ind = TRUE)[1, ]
  sprintf("%s at quality %s", rownames(table)[position[1]], colnames(table)[position[2]])
}


#' Check the table the user's own function returned, before DADA2 or the report uses it,
#' and give it the row and column names of the counts.
#'
#' DADA2 checks only that a table is numeric, has 16 rows and values from 0 to 1, and it
#' only warns about a rate of 0. A table that passes those checks can still change every
#' ASV of the run without a message, so MetaFlux asks more of the user's function, in
#' this order (the first problem found stops):
#'   1. a numeric matrix with as many rows and columns as the counts (16 x Q);
#'   2. row and column names: when the function gives none, those of the counts are
#'      copied; when it gives some, they must be those of the counts, in their order (a
#'      table whose columns are other quality scores would be read shifted). Any other
#'      attribute the function set on its table is dropped (see the comment in the code);
#'   3. no missing (NA, NaN) or infinite value;
#'   4. every value from 0 to 1: a rate is a probability;
#'   5. no value of exactly 0. DADA2 reads a rate of 0 as "this misread is impossible at
#'      this quality score", so sequences that carry it can no longer be corrected, and
#'      those seen in more than one read can become ASVs of their own: on one sample of
#'      the 16S test set, A2G set to 0 at every quality score turned 138 ASVs into 226
#'      (A2G at Q37 alone: 147). DADA2 only warns ("Zero in error matrix."), so the run
#'      would go on; its own error functions never return a 0, both keep every rate at
#'      or above 1e-7. Hence a stop here;
#'   6. for each true base, its four rows (A2A, A2C, A2G and A2T for A, and so on) sum to
#'      1 within 1e-6 at every quality score: a true base is read as one of the four.
#' Then, a warning and not a stop: rates below 1e-7, the lowest rate DADA2's own error
#' functions return (see the comment in the code).
#' The numbers are otherwise handed on as the function returned them.
#'
#' @param table What the function returned.
#' @param trans The counts it was given.
#' @param custom What load_custom_error_function() returned (its name and file, for the
#'   message).
#' @return The table, with the row and column names of trans. Stops otherwise, with a
#'   sentence that starts "the custom error function <name> (<file>) returned".
check_custom_table <- function(table, trans, custom) {
  returned <- sprintf("the custom error function %s (%s) returned", custom$name, custom$file)

  # 1. A numeric matrix of the shape of the counts.
  if (!is.matrix(table) || !is.numeric(table)) {
    what <- if (is.matrix(table)) paste("a matrix of type", typeof(table)) else paste("an object of class", class(table)[1])
    stop(sprintf("%s %s, not a numeric matrix", returned, what), call. = FALSE)
  }
  if (!identical(dim(table), dim(trans))) {
    stop(sprintf(paste0(
      "%s a matrix of %d rows and %d columns; it must have the shape of the counts it is ",
      "given: %d rows (A2A ... T2T) and %d columns (quality scores %s to %s)"),
      returned, nrow(table), ncol(table), nrow(trans), ncol(trans),
      colnames(trans)[1], colnames(trans)[ncol(trans)]), call. = FALSE)
  }

  # 2. The names: copied where missing, compared where given.
  if (is.null(rownames(table))) {
    rownames(table) <- rownames(trans)
  } else if (!identical(rownames(table), rownames(trans))) {
    stop(sprintf("%s a table whose rows are named %s ...; they must be those of the counts, in their order: %s",
                 returned, paste(head(rownames(table), 4), collapse = ", "),
                 paste(rownames(trans), collapse = ", ")), call. = FALSE)
  }
  if (is.null(colnames(table))) {
    colnames(table) <- colnames(trans)
  } else if (!identical(colnames(table), colnames(trans))) {
    stop(sprintf("%s a table whose columns are named %s ...; they must be those of the counts, the quality scores %s to %s in order",
                 returned, paste(head(colnames(table), 4), collapse = ", "),
                 colnames(trans)[1], colnames(trans)[ncol(trans)]), call. = FALSE)
  }
  # The names exactly as the counts carry them (also when the function gave the same
  # names in another form, such as a named list of names), and nothing else: any other
  # attribute the function set on its table (attributes are labels R can attach to an
  # object, such as a fitted model kept with the table for inspection, or a counter of
  # calls) is dropped. DADA2 compares the tables of two rounds with identical(), which
  # also compares attributes, so an attribute that changes from round to round would
  # keep it from ever seeing a repeat: the learning would run on to the flat-rounds
  # rule or the cap and count as not converged, with the numbers of a converged model.
  # The numbers themselves are not touched.
  attributes(table) <- list(dim = dim(trans), dimnames = dimnames(trans))

  # 3. Missing or infinite values.
  not_finite <- !is.finite(table)
  if (any(not_finite)) {
    stop(sprintf("%s missing or infinite values (NA, NaN or Inf) in %d cell(s), the first %s",
                 returned, sum(not_finite), first_cell_text(table, not_finite)), call. = FALSE)
  }

  # 4. Values outside 0 to 1.
  outside <- table < 0 | table > 1
  if (any(outside)) {
    stop(sprintf("%s values below 0 or above 1 in %d cell(s), the first %s (%s); a rate is a probability",
                 returned, sum(outside), first_cell_text(table, outside),
                 format(table[outside][1], digits = 6)), call. = FALSE)
  }

  # 5. Rates of exactly 0.
  zero <- table == 0
  if (any(zero)) {
    stop(sprintf(paste0(
      "%s a rate of exactly 0 in %d cell(s), the first %s. DADA2 reads a rate of 0 as a ",
      "misread that cannot happen, so sequences that carry it can no longer be corrected, ",
      "and those seen in more than one read can become ASVs of their own (DADA2 itself only ",
      "warns, \"Zero in error matrix.\"). ",
      "DADA2's own error functions keep every rate at or above 1e-7; keep yours there too"),
      returned, sum(zero), first_cell_text(table, zero)), call. = FALSE)
  }

  # 6. The four rows of each true base sum to 1 at every quality score.
  for (nucleotide in c("A", "C", "G", "T")) {
    rows <- paste0(nucleotide, "2", c("A", "C", "G", "T"))
    sums <- colSums(table[rows, , drop = FALSE])
    off <- abs(sums - 1) > 1e-6
    if (any(off)) {
      first <- which(off)[1]
      stop(sprintf(paste0(
        "%s rows %s that do not sum to 1 at %d quality score(s), the first at quality %s ",
        "(sum %s). They are the probabilities of reading a true %s as A, C, G or T, so they ",
        "must sum to 1 (within 1e-6) at every quality score"),
        returned, paste(rows, collapse = ", "), sum(off), colnames(table)[first],
        format(sums[[first]], digits = 10), nucleotide), call. = FALSE)
    }
  }

  # Then a warning, not a stop: rates below 1e-7, the lowest rate DADA2's own error
  # functions return (loessErrfun clamps its estimates to 1e-7 ... 0.25). A rate far
  # below the counted one acts like a weaker 0: on the R1 reads of one sample of the
  # 16S test set (28,648 reads), the loess model gave 138 ASVs, and the same model
  # with A2G at every quality score set to 1e-7 gave 146, to 1e-12 159, to 1e-30 186,
  # and to 1e-300 the 226 of a 0. It is a warning because the effect grows gradually
  # and a fitted curve can dip below 1e-7 at the highest quality scores, or at scores
  # without bases, where no read carries the rate. The text holds no round number, so
  # the learning step writes it once per read direction; the report does not show
  # warnings of its refits.
  if (any(table < 1e-7)) {
    warning(sprintf(paste0(
      "%s rates below 1e-7, the lowest rate DADA2's own error functions return. DADA2 ",
      "reads such a rate as a misread that almost never happens, so sequences that carry ",
      "it are rarely corrected, and those seen in more than one read can become ASVs of ",
      "their own"), returned), call. = FALSE)
  }
  table
}


#' The registry entry of the user's own function, in the form of the entries above.
#'
#' Its fit ignores quality_values and extend_range (a user's function gets the counts
#' only, as in plain DADA2), calls the function (call_custom_error_function) and checks
#' its table (check_custom_table). It has no anchors and leaves no value out.
#'
#' @param custom What load_custom_error_function() returned.
#' @return list(label, suits = NA, fit).
custom_registry_entry <- function(custom) {
  list(
    label = custom_function_label(custom$name),
    suits = NA_character_,
    fit   = function(trans, quality_values, extend_range = FALSE) {
      table <- call_custom_error_function(custom, trans)
      table <- check_custom_table(table, trans, custom)
      list(table = table, anchors = NA_integer_, dropped = integer(0), note = "")
    }
  )
}
