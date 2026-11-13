# makeBinnedQualErrfun, copied unchanged from DADA2 1.38.0.
#
# What it is: the error function DADA2's author wrote for binned quality scores. A
# NovaSeq, a NextSeq 1000/2000 or another recent Illumina instrument writes only a few
# quality values (for example 9, 20 and 34) instead of every value from 2 to about 40,
# and DADA2's usual fit (loessErrfun, a smooth curve over all quality scores) is not
# made for that. Given the quality values the reads carry, makeBinnedQualErrfun returns
# a function that DADA2's learnErrors calls once per round of learning, with that
# round's transition counts (the bases of each true nucleotide read as each nucleotide,
# at each quality score). For each of the 12 kinds of misread (A read as C, A as G, and
# so on) that function:
#   - divides the errors by the bases of the true nucleotide at each given quality
#     value, which gives the error rate counted there;
#   - joins these rates with straight lines between neighbouring given values, and
#     copies the rate of the lowest (highest) given value to every quality score below
#     (above) it;
#   - keeps every rate between 1e-7 and 0.25.
# It stops when the counts hold bases below the lowest or above the highest given
# value, and warns when the lowest or highest quality score that holds bases is not one
# of the given values.
#
# Who uses it: workflow/scripts/amplicon/60_error_functions.R, the list of the error
# functions MetaFlux can learn with, wraps it as the function called "binned". The
# wrapper (binned_anchors_with_data there) hands it only given values that hold bases in
# the counts, because a given value without bases makes it skip the straight lines on
# both sides of that value. The skipped stretch then gets missing rates (DADA2 fails on
# them) when fitted lines lie on both sides of it, and otherwise a flat copy of the
# nearest fitted rate, which is a wrong model without any message.
# 60d_learn_errors.R (rule learn_errors) and 60e_error_model_report.R (rule
# error_model_report) source this file first, then 60_error_functions.R.
#
# Why the file is here: MetaFlux pins dada2 1.34.0 (workflow/envs/dada2.yaml), and this
# function first appeared in dada2 1.36.0. It is self-contained base R (it calls no
# other DADA2 code), so it is copied here instead of moving MetaFlux to a newer DADA2:
# keeping 1.34.0 keeps both DADA2 environments (dada2.yaml and taxonomy.yaml) and every
# number of a run on continuous quality values (MiSeq) as it was.
#
# Source: dada2 1.38.0 from Bioconductor 3.22
#   https://bioconductor.org/packages/3.22/bioc/src/contrib/dada2_1.38.0.tar.gz
#   (md5 of the tarball: 8675dab95ed519215a2662e1f71f914d), file R/errorModels.R,
#   lines 96-156.
# Authors, as DADA2 lists them: Benjamin Callahan, Paul McMurdie, Susan Holmes.
# License: this file stays under DADA2's license, not under MetaFlux's MIT license.
#   DADA2's DESCRIPTION file declares "License: LGPL-2"; the LICENSE file shipped in
#   the same tarball holds the text of the GNU Lesser General Public License version 3,
#   which is copied unchanged next to this file as LICENSE.dada2.
# Changes: none. The lines between the two markers below are the original bytes,
#   trailing spaces included. tests/test_error_model_settings.py checks them against
#   the md5 on the next line, so a hand edit, or an editor that strips trailing spaces,
#   fails the tests.
# md5 of the vendored block: 9b320182dead53a31334364a026c957e
#   (of the lines between the markers, each ending with a line feed)
# Do not edit the block. A newer version of the function from DADA2 goes in as a new
# copy, with its own source line and md5, and a line in the CHANGELOG.
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
