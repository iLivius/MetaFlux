# An error function of your own for MetaFlux (EXPERIMENTAL: error_model: custom):
# DADA2's binned fit, with the bins taken from the counts.
#
# Shown as an example under "Your own error function" in docs/amplicon/error-model.md.
# It shows two things a file of your own may need: carrying a function of a newer DADA2
# than the one MetaFlux runs, and taking from the counts the quality scores a fit is
# anchored on.
#
# How to use it: in the amplicon: block your run config already has,
#   amplicon:
#     dada2:
#       learn_errors:
#         error_model: custom
#         custom_error_function:
#           file: .test/custom_error_functions/binned_from_counts.R
#           name: binned_from_counts
# (a relative path is read from the folder snakemake is started in). In plain DADA2:
#   source("binned_from_counts.R")
#   err <- learnErrors(filtered_files, errorEstimationFunction = binned_from_counts)
#
# What DADA2 hands the function, and what it expects back (DADA2's own contract for an
# errorEstimationFunction):
#   Input, trans: the transition counts of one round of learning. In each round DADA2
#     denoises the reads it drew, compares every read with the ASV it was assigned to,
#     and counts, at each quality score, how often each true base was read as each
#     base. 16 rows named A2A, A2C, A2G, A2T, C2A, ... T2T (the true base, "2", the base
#     it was read as), one column per quality score, named "0", "1", "2", ... (at least
#     41 columns, scores 0 to 40). The quality score of a count is a consensus quality:
#     the mean quality, rounded, of the identical reads of a unique sequence at that
#     position.
#   Output: a table of the same shape. Each cell is the probability that the true base
#     of the row is read as the second base of the row at that quality score; the four
#     rows of each true base sum to 1 at every score.
#
# What it does, in each round:
#   1. It takes as bins every quality score at which each of the four true bases (A, C,
#      G and T) has bases in the counts, and adds the lowest and the highest score that
#      hold any bases (see the comment in the function for why).
#   2. It hands these bins to makeBinnedQualErrfun, DADA2's binned fit, copied below.
#      For each of the 12 kinds of misread (A read as C, A as G, and so on), that
#      function divides the misreads by all bases of the true base at each bin, which
#      gives the error rate counted there, and joins these rates with straight lines.
#      Below the lowest bin and above the highest it copies the rate of the nearest bin.
#      It keeps every misread rate between 1e-7 and 0.25 and sets each read-correctly row
#      (A2A, C2C, G2G, T2T) to 1 minus the three misread rates of its base.
#
# How it differs from MetaFlux's built-in binned fit (error_model: binned)
#   The built-in fit anchors on the quality values the filtered reads carry, which
#   MetaFlux reads from the reads themselves (9, 20 and 34 on a NextSeq 2000). This
#   function takes its bins from the counts, and the counts have bases at more scores
#   than the reads carry, because their scores are consensus qualities.
#   - On MiSeq reads the bins are the 30 or so scores the reads carry, and the fit joins
#     the counted rates without any smoothing.
#   - On binned reads the bins include the scores in between (on the NextSeq 2000 test
#     set, 25 or 26 bins from 9 to 34), at which few misreads were counted. The result
#     is not the built-in binned fit; use error_model: binned for binned reads.
#
# Why makeBinnedQualErrfun is copied into this file
#   It first appeared in DADA2 1.36.0, and MetaFlux runs DADA2 1.34.0, which does not
#   have it. MetaFlux carries its own copy (workflow/scripts/amplicon/vendor/), but
#   that copy is read by MetaFlux's scripts, not by a file of yours, and a file of your
#   own should not read another file with source(): MetaFlux tracks only this one file,
#   so an edit of the other file would not make the models be learned again. Hence the
#   copy here.
#
# Source and license
#   - The block between the two marker lines below is makeBinnedQualErrfun from DADA2
#     1.38.0 (Bioconductor 3.22, dada2_1.38.0.tar.gz, file R/errorModels.R, lines
#     96-156), copied unchanged, trailing spaces included. Its md5 (of the lines between
#     the markers, each ending with a line feed) is 9b320182dead53a31334364a026c957e, the
#     md5 of MetaFlux's own copy. Do not edit the block: an editor that strips trailing
#     spaces changes it.
#   - Authors, as DADA2 lists them: Benjamin Callahan, Paul McMurdie, Susan Holmes. The
#     block stays under DADA2's license: DADA2's DESCRIPTION file declares LGPL-2, and
#     the LICENSE file of the same tarball holds the text of the GNU Lesser General
#     Public License version 3, which MetaFlux ships as
#     workflow/scripts/amplicon/vendor/LICENSE.dada2.
#   - The rest of this file (these comments and binned_from_counts below the block) was
#     written for MetaFlux and is under MetaFlux's MIT license.
#
# R packages: none beyond base R. Nothing here draws random numbers, so the same counts
# always give the same table, as DADA2 needs to see a round repeat.


# ---- makeBinnedQualErrfun, DADA2 1.38.0 (LGPL; see "Source and license" above) ----
# Given the bins, it returns the function that does the fit (a function made by a
# function); binned_from_counts below calls it once per round with that round's bins.
# ---- vendored block: start ----
makeBinnedQualErrfun <- function(binnedQ) {
  if(is.null(binnedQ)) { stop("The quality scores used in your data must be provided.") } 
  function(trans, binnedQuals=binnedQ) {
    qq <- as.numeric(colnames(trans))
    # Get min and max observed quality scores
    qmax <- max(qq[colSums(trans)>0])
    qmin <- min(qq[colSums(trans)>0])
    # Check for data consistency with provided binned qualities
    if(qmax > max(binnedQuals)) stop("Input data contains a higher quality score than the provided binned values.")
    if(qmin < min(binnedQuals)) stop("Input data contains a lower quality score than the provided binned values.")
    if(!qmax %in% binnedQuals) warning("Maximum observed quality score is not in the provided binned values.")
    if(!qmin %in% binnedQuals) warning("Minimum observed quality score is not in the provided binned values.")
    
    est <- matrix(0, nrow=0, ncol=length(qq))
    for(nti in c("A","C","G","T")) {
      for(ntj in c("A","C","G","T")) {
        if(nti != ntj) {
          errs <- trans[paste0(nti,"2",ntj),]
          tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
          p <- errs/tot
          df <- data.frame(q=qq, errs=errs, tot=tot, p=p)
          # Check and enforce that this q scores start at zero
          if(!all(df$q == seq(nrow(df))-1)) stop("Unexpected Q score series.") ###!
          pred <- rep(NA, nrow(df))
          for(i in seq(length(binnedQuals)-1)) {
            loQ <- binnedQuals[i]
            hiQ <- binnedQuals[i+1]
            loP <- df$p[loQ+1]
            hiP <- df$p[hiQ+1]
            # Linear interpolation between the binned Q scores observed in the data
            if(!is.na(loP) && !is.na(hiP)) {
              pred[(loQ+1):(hiQ+1)] <- seq(loP, hiP, length.out=(hiQ-loQ+1))
            }
          }
          
          maxrli <- max(which(!is.na(pred)))
          minrli <- min(which(!is.na(pred)))
          pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
          pred[seq_along(pred)<minrli] <- pred[[minrli]]
          est <- rbind(est, pred)
        } # if(nti != ntj)
      } # for(ntj in c("A","C","G","T"))
    } # for(nti in c("A","C","G","T"))
    
    # HACKY
    MAX_ERROR_RATE <- 0.25
    MIN_ERROR_RATE <- 1e-7
    est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
    est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
    
    # Expand the err matrix with the self-transition probs
    err <- rbind(1-colSums(est[1:3,]), est[1:3,],
                 est[4,], 1-colSums(est[4:6,]), est[5:6,],
                 est[7:8,], 1-colSums(est[7:9,]), est[9,],
                 est[10:12,], 1-colSums(est[10:12,]))
    rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
    colnames(err) <- colnames(trans)
    # Return
    return(err)
  }
}
# ---- vendored block: end ----


# ---- The error function ----------------------------------------------------------

#' binned_from_counts: DADA2's binned fit, with the bins taken from the counts.
#'
#' Input: trans, the transition counts of one round (see the top of the file), handed
#'   over by DADA2's learnErrors (inside MetaFlux: rule learn_errors, through MetaFlux's
#'   checks), or by MetaFlux's report when it refits the function on the counts of a
#'   saved model.
#' What it does: chooses the bins (steps 1 and 2 below) and hands them, with the
#'   counts, to makeBinnedQualErrfun (step 3).
#' Output: the 16 x Q table of error rates, rows and columns named as in trans. DADA2
#'   denoises the next round with it; MetaFlux first checks it (shape, names, no missing
#'   value, values from 0 to 1, no 0, the rows of each true base summing to 1).
binned_from_counts <- function(trans) {
  quality_scores <- as.numeric(colnames(trans))
  nucleotides <- c("A", "C", "G", "T")

  # 1. The scores at which each of the four true bases has bases. At a score where one
  #    true base has none, its error rates would be 0 divided by 0 (missing), and
  #    makeBinnedQualErrfun would drop its straight lines on both sides of that score.
  #    Such a score is no bin; the line between its neighbouring bins covers it.
  has_all_four <- rep(TRUE, length(quality_scores))
  for (true_base in nucleotides) {
    rows_of_base <- paste0(true_base, "2", nucleotides)   # A2A, A2C, A2G, A2T for A
    bases_of_this_base <- colSums(trans[rows_of_base, , drop = FALSE])
    has_all_four <- has_all_four & bases_of_this_base > 0
  }

  # Two such scores are the least a straight line needs. Real counts have about 25 to
  # 32; fewer means that almost nothing was counted.
  if (sum(has_all_four) < 2) {
    stop("binned_from_counts needs two or more quality scores at which all four true ",
         "bases have bases; these counts have ", sum(has_all_four), call. = FALSE)
  }
  bins <- quality_scores[has_all_four]

  # 2. The lowest and the highest score that hold any bases are always bins.
  #    makeBinnedQualErrfun stops when the counts hold bases below its lowest or above
  #    its highest bin, and at these two scores one true base can be missing: on the
  #    COI fwhF2/fwhR2n test set, quality 7 of the R1 counts holds three bases, none of
  #    them a true A. At such a score the missing base gets the rate of the next bin
  #    (makeBinnedQualErrfun copies it flat); the other bases get their counted rate.
  scores_with_bases <- quality_scores[colSums(trans) > 0]
  bins <- sort(unique(c(min(scores_with_bases), bins, max(scores_with_bases))))

  # 3. The fit: makeBinnedQualErrfun(bins) makes the fitting function for these bins,
  #    which is then called on the counts.
  makeBinnedQualErrfun(bins)(trans)
}
