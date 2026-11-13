# A minimal error function of your own for MetaFlux (EXPERIMENTAL: error_model: custom).
#
# Used by the dry-run config .test/config/16S_custom_error_function.yaml, and the example
# of "Your own error function" in docs/amplicon/error-model.md.
#
# What it does: nothing of its own. It hands DADA2's transition counts to DADA2's default
# error function, loessErrfun, and returns that table unchanged. A run with
# error_model: custom and this function therefore learns the same error models as a run
# with error_model: loess (same model codes, same ASV tables), which is how the hook was
# checked. A function of your own replaces the body with your own fit.
#
# The contract is DADA2's errorEstimationFunction, so the same function works in plain
# DADA2: learnErrors(files, errorEstimationFunction = loess_wrapper).
#   Input, trans: DADA2's matrix of transition counts of one round of learning, whole
#     numbers. 16 rows named A2A, A2C, A2G, A2T, C2A, ... T2T (the true base, "2", the base
#     it was read as) and one column per quality score, named "0", "1", ... (at least 41
#     columns, scores 0 to 40, and more when the reads carry higher scores).
#   Output: a numeric matrix of the same shape. Each cell is the probability that the
#     true base of the row is read as the second base of the row at that quality score:
#     above 0 and at most 1, and the four rows of each true base sum to 1 at every score.
loess_wrapper <- function(trans) {
  dada2::loessErrfun(trans)
}
