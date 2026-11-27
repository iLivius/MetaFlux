# An error function of your own for MetaFlux (EXPERIMENTAL: error_model: custom):
# a monotone GAM, an error rate that never rises as the quality score rises.
#
# Shown as an example under "Your own error function" in docs/amplicon/error-model.md.
# It shows the parts a fitted error function needs besides the fit itself: settings,
# fall-backs for thin counts, rates held flat where no read gives a rate, DADA2's
# bounds for a misread rate, and read-correctly rows that make each true base's rows
# sum to 1.
#
# How to use it: in the amplicon: block your run config already has,
#   amplicon:
#     dada2:
#       learn_errors:
#         error_model: custom
#         custom_error_function:
#           file: .test/custom_error_functions/monotone_gam.R
#           name: monotone_gam
# (a relative path is read from the folder snakemake is started in). With error_model
# auto, loess or binned and the same custom_error_function block, the models are learned
# as usual and the report draws this function next to them, refitted on the same counts.
# In plain DADA2 (with the R package scam installed):
#   source("monotone_gam.R")
#   err <- learnErrors(filtered_files, errorEstimationFunction = monotone_gam)
#
# What DADA2 hands the function, and what it expects back (DADA2's own contract for an
# errorEstimationFunction):
#   Input, trans: the transition counts of one round of learning. In each round DADA2
#     denoises the reads it drew, compares every read with the ASV it was assigned to,
#     and counts, at each quality score, how often each true base was read as each
#     base. 16 rows named A2A, A2C, A2G, A2T, C2A, ... T2T (the true base, "2", the base
#     it was read as), one column per quality score, named "0", "1", "2", ... (at least
#     41 columns, scores 0 to 40).
#   Output: a table of the same shape. Each cell is the probability that the true base
#     of the row is read as the second base of the row at that quality score; the four
#     rows of each true base sum to 1 at every score.
#
# How it fits, for each of the 12 kinds of misread (A read as C, A as G, ... T as G):
#   1. At each quality score that holds bases of the true base, the counted error rate
#      is the misreads of this kind divided by all bases of the true base there.
#   2. A curve of error rate against quality score is fitted to these rates. The curve
#      is a GAM (generalised additive model: a regression whose curve is not a straight
#      line chosen in advance but a smooth curve whose shape is learned from the data)
#      from the R package scam (shape constrained additive models), which can force the
#      curve to only fall: bs = "mpd" is a monotone decreasing smooth. A higher quality
#      score promises fewer errors, so the rate may not rise with it.
#   3. The fit is binomial: each quality score counts as so many bases (the bases there
#      are the weights), of which so many were misread. The curve is fitted on the
#      log-odds scale, so it never leaves 0 to 1, and a score with two million bases
#      counts far more than a score with twenty.
#   4. The curve is read off at every score from the lowest to the highest score with
#      bases, gaps included. Below and above that range it is held flat at its end
#      value, as DADA2's own loessErrfun does: no read gives a rate there.
#   5. Every misread rate is kept between 1e-7 and 0.25, the bounds of DADA2's own error
#      functions (loessErrfun and makeBinnedQualErrfun). A rate of 0 would tell DADA2
#      that the misread cannot happen, so reads that carry it could become ASVs of their
#      own; MetaFlux stops on a 0 for that reason.
#   Then, for each true base, the read-correctly row (A2A, C2C, G2G, T2T) is 1 minus
#   its three misread rates, so its four rows sum to 1 at every score.
#
# Fall-backs, when a kind of misread has too little data for the curve:
#   - Fewer than MIN_SCORES_FOR_CURVE scores hold bases of the true base, or scam fails:
#     a logit-linear fit (a binomial GLM, generalised linear model, whose log-odds are a
#     straight line in the quality score), used only when the line falls or stays flat
#     as quality rises. When scam fails the function warns, so a fall-back that was not
#     planned does not pass unnoticed (MetaFlux writes the warning to its log).
#   - The line would rise, or only one score holds bases: one pooled rate for every
#     score (all misreads of this kind divided by all bases of the true base).
#   - No misread of this kind was counted at all: the lowest allowed rate, 1e-7.
# On binned reads (a NextSeq or NovaSeq writes only a few quality values), the counts
# still hold bases at many scores, because DADA2 counts at consensus qualities, the
# mean quality of the identical reads of a sequence: every score from 9 to 34 on the
# NextSeq 2000 test set. The curve is then fitted as on MiSeq reads; the fall-backs are
# for counts much thinner than that.
#
# Reproducibility: nothing here draws random numbers, so the same counts always give the
# same table. DADA2 ends the learning only when a round gives back a table it has seen
# before, and MetaFlux's report refits the function on the counts behind the saved model
# and expects the same table back. For the same reason the function returns a plain
# table, with no attributes of its own.
#
# R packages: scam (it loads mgcv, on which it is built, and nlme), and base R's stats.
# MetaFlux runs your function in workflow/envs/dada2_custom.yaml, which holds them
# (scam 1.2-22, mgcv 1.9-4). conda-forge has no scam for Apple Silicon Macs or Linux on
# ARM, so this file cannot run there. library(scam) attaches scam alone to R's search
# path (for the whole R session); mgcv and nlme are loaded without being attached.
#
# License: written for MetaFlux; MIT license, like MetaFlux.

library(scam)


# ---- Settings --------------------------------------------------------------------

# DADA2's own bounds for a misread rate (loessErrfun and makeBinnedQualErrfun).
MIN_ERROR_RATE <- 1e-7
MAX_ERROR_RATE <- 0.25

# The curve needs enough quality scores with bases. scam's smallest monotone curve is
# built from 4 pieces (its basis size k; scam refuses k = 3 for bs = "mpd"), and at
# least 2 scores more than pieces are kept, so that the curve cannot simply pass through
# every point: 6 scores, k = 4.
MIN_SCORES_FOR_CURVE <- 6

# The largest k. 10 is the usual size for a curve of one variable in mgcv and scam;
# quality scores span about 40 values, and a curve that may only fall needs few pieces.
MAX_BASIS_SIZE <- 10


# ---- Helpers ---------------------------------------------------------------------

#' Hold a curve flat outside the quality scores that hold bases.
#'
#' Input: the rates read off the curve at the scores from the lowest to the highest
#'   score with bases (predicted_inside), and which of all scores those are (inside, a
#'   TRUE/FALSE value per column of trans).
#' Output: one rate per column of trans: the predicted rates inside, the first predicted
#'   rate below them and the last predicted rate above them.
#' Used by fit_one_misread() below, whichever fit it used.
hold_flat_outside <- function(predicted_inside, inside) {
  rates <- rep(NA_real_, length(inside))
  rates[inside] <- predicted_inside
  first_inside <- min(which(inside))
  last_inside <- max(which(inside))
  column_number <- seq_along(inside)
  rates[column_number < first_inside] <- predicted_inside[1]
  rates[column_number > last_inside] <- predicted_inside[length(predicted_inside)]
  rates
}


#' Fit the error rate of one kind of misread (say A read as G) against quality.
#'
#' Input (from monotone_gam below, for one misread row of trans):
#'   misreads       : the misreads of this kind at each quality score (one row of trans);
#'   bases          : all bases of the true base at each score (its four rows summed);
#'   quality_scores : the quality score of each column (from the column names of trans);
#'   label          : the row name, such as "A2G", for the warning text.
#' What it does: the curve of steps 1 to 5 at the top of the file, or a fall-back.
#' Output: list(rates = one misread rate per column, between MIN_ERROR_RATE and
#'   MAX_ERROR_RATE; method = which fit was used: "scam", "glm", "pooled" or "no
#'   misreads"). monotone_gam puts the rates into its table; the method is there for a
#'   person who tries the function by hand.
fit_one_misread <- function(misreads, bases, quality_scores, label) {
  has_bases <- bases > 0
  scores_with_bases <- sum(has_bases)

  # No base of the true base, or no misread of this kind anywhere: nothing to fit. The
  # lowest allowed rate says "rare", never "impossible".
  if (scores_with_bases == 0 || sum(misreads[has_bases]) == 0) {
    lowest <- rep(MIN_ERROR_RATE, length(quality_scores))
    return(list(rates = lowest, method = "no misreads"))
  }

  # One row per quality score that holds bases: the score, the counted rate, and the
  # bases there, which are the weight of that score in the fits below.
  counted <- data.frame(
    quality = quality_scores[has_bases],
    rate    = misreads[has_bases] / bases[has_bases],
    weight  = bases[has_bases]
  )

  # The scores the curve is read off at: every column from the lowest to the highest
  # score with bases, gaps included (the curve fills them).
  inside <- quality_scores >= min(counted$quality) & quality_scores <= max(counted$quality)
  scores_inside <- data.frame(quality = quality_scores[inside])

  predicted_inside <- NULL
  method <- NULL

  # 1. The monotone decreasing curve (scam, bs = "mpd"), binomial with the bases as
  #    weights. k grows with the scores available, up to MAX_BASIS_SIZE. An error of
  #    scam becomes a warning, and the fall-back below takes over.
  if (scores_with_bases >= MIN_SCORES_FOR_CURVE) {
    k <- min(MAX_BASIS_SIZE, scores_with_bases - 2)
    curve_model <- tryCatch(
      scam(rate ~ s(quality, bs = "mpd", k = k), family = binomial(link = "logit"),
           data = counted, weights = weight),
      error = function(e) {
        warning(sprintf("monotone_gam: scam failed for %s (%s); a logit-linear fit is used instead",
                        label, conditionMessage(e)), call. = FALSE)
        NULL
      })
    if (!is.null(curve_model)) {
      predicted_inside <- as.numeric(predict(curve_model, newdata = scores_inside, type = "response"))
      method <- "scam"
    }
  }

  # 2. Fall-back: a binomial GLM, log-odds a straight line in the quality score. Used
  #    only when the line falls or stays flat as quality rises.
  if (is.null(predicted_inside) && scores_with_bases >= 2) {
    line_model <- tryCatch(
      glm(rate ~ quality, family = binomial(link = "logit"), data = counted, weights = weight),
      error = function(e) NULL)
    if (!is.null(line_model)) {
      slope <- coef(line_model)[["quality"]]
      if (is.finite(slope) && slope <= 0) {
        predicted_inside <- as.numeric(predict(line_model, newdata = scores_inside, type = "response"))
        method <- "glm"
      }
    }
  }

  # 3. Fall-back: one pooled rate for every score.
  if (is.null(predicted_inside)) {
    pooled_rate <- sum(misreads[has_bases]) / sum(bases[has_bases])
    predicted_inside <- rep(pooled_rate, sum(inside))
    method <- "pooled"
  }

  # Flat beyond the scores with bases (step 4), then DADA2's bounds (step 5).
  rates <- hold_flat_outside(predicted_inside, inside)
  rates <- pmax(rates, MIN_ERROR_RATE)
  rates <- pmin(rates, MAX_ERROR_RATE)
  list(rates = rates, method = method)
}


# ---- The error function ----------------------------------------------------------

#' monotone_gam: one monotone decreasing curve per kind of misread.
#'
#' Input: trans, one round's transition counts (see the top of the file), handed over by
#'   DADA2's learnErrors (inside MetaFlux: rule learn_errors, through MetaFlux's checks),
#'   or by MetaFlux's report when it refits the function on the counts of a saved model.
#' What it does: fits the 12 misread rates with fit_one_misread(), then sets each
#'   read-correctly row to 1 minus the three misread rates of its true base.
#' Output: the 16 x Q table of error rates, rows and columns named as in trans. DADA2
#'   denoises the next round with it; MetaFlux first checks it (shape, names, no missing
#'   value, values from 0 to 1, no 0, the rows of each true base summing to 1).
monotone_gam <- function(trans) {
  # The quality score of each column. DADA2 names the columns "0", "1", ...; a table
  # without column names is read as starting at score 0.
  if (is.null(colnames(trans))) {
    quality <- seq_len(ncol(trans)) - 1
  } else {
    quality <- as.numeric(colnames(trans))
  }

  nucleotides <- c("A", "C", "G", "T")
  row_names <- paste0(rep(nucleotides, each = 4), "2", nucleotides)   # A2A, A2C, ... T2T
  rates <- matrix(NA_real_, nrow = 16, ncol = ncol(trans),
                  dimnames = list(row_names, colnames(trans)))

  for (true_base in nucleotides) {
    # All bases of this true base at each quality score: its four rows summed.
    rows_of_this_base <- paste0(true_base, "2", nucleotides)
    bases <- colSums(trans[rows_of_this_base, , drop = FALSE])

    # The three kinds of misread of this true base, one curve each.
    misread_rows <- paste0(true_base, "2", setdiff(nucleotides, true_base))
    for (row in misread_rows) {
      fit <- fit_one_misread(misreads = as.numeric(trans[row, ]), bases = bases,
                             quality_scores = quality, label = row)
      rates[row, ] <- fit$rates
    }

    # Read correctly: whatever probability the three misreads leave. With each misread
    # rate at most 0.25, this is at least 0.25.
    correct_row <- paste0(true_base, "2", true_base)
    rates[correct_row, ] <- 1 - colSums(rates[misread_rows, , drop = FALSE])
  }

  rates
}
