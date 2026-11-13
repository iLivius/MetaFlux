# The R side of the experimental error_model: custom, checked outside Snakemake on
# made-up tables, so that a change that breaks one of its promises fails a test.
#
# Who runs it: tests/test_custom_error_function_r.py, only when the environment variable
# METAFLUX_RSCRIPT names an Rscript (the CI runs no R, so the test is skipped there).
# Called with one argument, the folder workflow/scripts/amplicon of the repository.
#
# What it needs: base R only. No dada2, no reads and no conda environment: the
# functions under test never call DADA2 themselves, and the made-up error functions
# below return fixed tables.
#
# What it loads:
#   - 60_error_functions.R, whole: load_custom_error_function, call_custom_error_function,
#     check_custom_table and custom_registry_entry, the code that reads, calls and checks
#     a user's function (sourcing the file only defines functions);
#   - from 60d_learn_errors.R, only the function definitions (make_round_keeper, the
#     "round keeper" that sits between DADA2 and the error function, and its helpers)
#     and the two constants of the flat-rounds rule. The rest of that script reads the
#     snakemake object of a running job, so it is parsed but never run, as
#     tests/helpers.py does for the Python files with load_functions().
#
# What it checks, one line per case:
#   - each broken table stops with its reason, and a valid table passes (with or
#     without row and column names);
#   - an attribute the function attached to its table is dropped (DADA2 compares the
#     tables of two rounds with identical(), which also compares attributes);
#   - a rate below 1e-7 passes, with one warning;
#   - R's random-number state is the same after a call of the function as before it;
#   - R's kind of random numbers is the same after reading a file that switches it;
#   - a file that defines only a name starting with a dot lists that name in the stop;
#   - the round keeper records two warnings that come back in every round once each.
# Output: one line per case on stdout, "<case><TAB><the message, or OK>", which the test
# reads and compares with the fragments of docs/troubleshooting.md.

args <- commandArgs(trailingOnly = TRUE)
scripts_folder <- args[1]
source(file.path(scripts_folder, "60_error_functions.R"))

# The function definitions of 60d_learn_errors.R and its two constants, without the rest.
# parse() reads the file as R code without running it and returns its top-level
# statements one by one. A statement is kept when it is an assignment with <- to a
# plain name ("name <- ...") and what is assigned is either a function
# ("name <- function(...) ...") or one of the two constants; only those are run
# (eval), so the lines that read the snakemake object, open the log or call DADA2 are
# skipped.
for (statement in parse(file.path(scripts_folder, "60d_learn_errors.R"), keep.source = FALSE)) {
  is_assignment <- is.call(statement) && identical(statement[[1]], as.name("<-")) &&
    is.name(statement[[2]])
  if (!is_assignment) next
  assigned_name  <- as.character(statement[[2]])
  assigned_value <- statement[[3]]
  is_function <- is.call(assigned_value) && identical(assigned_value[[1]], as.name("function"))
  is_constant <- assigned_name %in% c("FLAT_DEPTH_FACTOR", "FLAT_MIN_GAIN")
  if (is_function || is_constant) eval(statement, globalenv())
}

# A made-up round of transition counts, shaped as DADA2 hands them to an error function
# (16 rows A2A ... T2T, one column per quality score 0 to 40, whole numbers), and a
# valid table of rates for it: every misread 0.001, every read-correctly row 0.997, so
# the four rows of each true base sum to 1.
rows <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
trans <- matrix(1000L, nrow = 16, ncol = 41, dimnames = list(rows, as.character(0:40)))
valid <- matrix(0.001, nrow = 16, ncol = 41, dimnames = dimnames(trans))
valid[c("A2A", "C2C", "G2G", "T2T"), ] <- 0.997

# What load_custom_error_function() would return, for the messages of the checks.
custom <- list(name = "my_errfun", file = "/data/my_errfun.R", state = new.env())

# One line of output.
report <- function(case, text) {
  cat(case, "\t", gsub("\n", " ", text), "\n", sep = "")
}

# A small R file in R's temporary folder, holding `lines`; returns its path.
write_function_file <- function(lines) {
  path <- tempfile(fileext = ".R")
  writeLines(lines, path)
  path
}


# ── 1. The six checks of check_custom_table, and two valid tables ──
# Each case is a table a function could return, made from the valid one; the comment
# says which check must stop it (the numbers of check_custom_table in
# 60_error_functions.R).

# check 1: the same numbers written as text
char_matrix <- matrix(as.character(valid), nrow = 16, dimnames = dimnames(valid))
# check 2: rows named in reverse order
other_rows <- valid
rownames(other_rows) <- rev(rownames(valid))
# check 2: columns named 1 to 41 instead of 0 to 40 (shifted by one)
other_cols <- valid
colnames(other_cols) <- as.character(1:41)
# check 3: a missing value
missing_value <- valid
missing_value["A2C", "0"] <- NA
# check 4: a rate above 1
above_one <- valid
above_one["A2C", "0"] <- 1.5
# check 5: A2G of exactly 0 at Q37 (A2A takes up the rest, so the rows still sum to 1)
zero_rate <- valid
zero_rate["A2G", "37"] <- 0
zero_rate["A2A", "37"] <- 0.998
# check 6: the four rows of a true A sum to 1.003
rows_not_one <- valid
rows_not_one["A2A", ] <- 1

table_cases <- list(
  valid        = valid,                  # passes
  no_names     = unname(valid),          # passes: the names of the counts are copied
  data_frame   = as.data.frame(valid),   # check 1: a data frame, not a matrix
  char_matrix  = char_matrix,
  wrong_shape  = valid[, 1:40],          # check 1: one quality column missing
  other_rows   = other_rows,
  other_cols   = other_cols,
  missing      = missing_value,
  above_one    = above_one,
  zero         = zero_rate,
  rows_not_one = rows_not_one
)
for (case in names(table_cases)) {
  outcome <- tryCatch({
    check_custom_table(table_cases[[case]], trans, custom)
    "OK"
  }, error = function(e) conditionMessage(e))
  report(case, outcome)
}


# ── 2. An attribute attached to the table is dropped ──
# A function that keeps its fitted model with the table would make each round's table
# differ from the last in that attribute, and DADA2 would never see a repeat.
with_attribute <- valid
attr(with_attribute, "fit") <- list(round = 3)
checked <- check_custom_table(with_attribute, trans, custom)
attributes_left <- paste(names(attributes(checked)), collapse = ", ")
report("attribute_dropped",
       if (identical(checked, valid)) "OK" else paste("attributes left:", attributes_left))


# ── 3. A rate below 1e-7 passes, with a warning ──
# withCallingHandlers catches each warning as it is raised, notes its text and lets the
# check go on (muffleWarning), as the round keeper of the learning step does.
tiny <- valid
tiny["A2G", "40"] <- 1e-12
tiny["A2A", "40"] <- 1 - sum(tiny[c("A2C", "A2G", "A2T"), "40"])
warnings_seen <- character(0)
outcome <- withCallingHandlers(
  tryCatch({
    check_custom_table(tiny, trans, custom)
    "passed"
  }, error = function(e) paste("stopped:", conditionMessage(e))),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
report("tiny_rate", sprintf("%s with %d warning(s): %s", outcome, length(warnings_seen),
                            paste(warnings_seen, collapse = " | ")))


# ── 4. The random-number state is put back after a call ──
# The function sets a seed of its own and draws numbers; learnErrors' draw of the
# samples must not notice.
set.seed(42)
state_before <- .Random.seed
custom$fn <- function(trans) {
  set.seed(1)
  runif(10)
  valid
}
invisible(call_custom_error_function(custom, trans))
report("rng_guard", if (identical(.Random.seed, state_before)) "OK" else "state changed")


# ── 5. The kind of random numbers is put back after reading the file ──
# A file written for parallel work often starts with RNGkind("L'Ecuyer-CMRG"); the
# learning step's seed would then start another series of numbers.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
kind_before <- RNGkind()
rngkind_file <- write_function_file(c(
  "RNGkind(\"L'Ecuyer-CMRG\")",
  "rngkind_errfun <- function(trans) trans"))
invisible(load_custom_error_function(rngkind_file, "rngkind_errfun"))
kind_after <- RNGkind()
report("rng_kind",
       if (identical(kind_after, kind_before)) "OK" else paste("kind left at", paste(kind_after, collapse = "/")))


# ── 6. A name that starts with a dot is listed in the stop ──
dot_file <- write_function_file(".my_errfun <- function(trans) trans")
outcome <- tryCatch({
  load_custom_error_function(dot_file, "my_errfun")
  "loaded"
}, error = function(e) conditionMessage(e))
report("dot_name", sub(".*does not define", "does not define", outcome))


# ── 7. The round keeper writes each distinct warning once ──
# A function that gives two different warnings at every call (as a model fitted to each
# kind of misread separately may), called three times: rounds 0, 1 and 2.
warning_file <- write_function_file(c(
  "two_warnings <- function(trans) {",
  "  warning(\"A2C: fit did not reach full convergence\")",
  "  warning(\"A2G: fit did not reach full convergence\")",
  "  table <- matrix(0.001, nrow = 16, ncol = ncol(trans), dimnames = dimnames(trans))",
  "  table[c(\"A2A\", \"C2C\", \"G2G\", \"T2T\"), ] <- 0.997",
  "  table",
  "}"))
loaded <- load_custom_error_function(warning_file, "two_warnings")
error_function_registry[["custom"]] <- custom_registry_entry(loaded)
keeper <- make_round_keeper("R1", "custom", integer(0), max_flat_rounds = 20L, max_consist = 200L)
suppressMessages(for (call in 1:3) keeper$errfun(trans))     # the keeper's log lines are messages
report("keeper_warnings", sprintf("%d call(s), %d warning(s) recorded: %s", keeper$state$calls,
                                  length(keeper$state$fit_warnings),
                                  paste(keeper$state$fit_warnings, collapse = " | ")))

unlink(c(rngkind_file, dot_file, warning_file))
