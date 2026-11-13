"""The R side of the user's own error function (EXPERIMENTAL, error_model: custom),
checked against what the documentation promises.

Why a test of its own: the rest of the suite checks the Python half of this setting
(the start-up check in 00_common.smk and the rule wiring), but nothing there would
fail if the R half broke: the six checks of every table the function returns, the
random-number guard around each call, the dropping of attributes, the warning on
rates below 1e-7, or the round keeper writing each warning once. Those live in
workflow/scripts/amplicon/60_error_functions.R and 60d_learn_errors.R.

How: the test runs .test/custom_error_functions/custom_function_checks.R once with
Rscript. That script feeds made-up tables to the R functions (base R only: no dada2,
no reads, no conda) and prints one line per case, "<case><TAB><message or OK>". The
tests below read those lines:
  - each broken table stops with its reason, and the fragment of that reason is a row
    of the table in docs/troubleshooting.md, so the page and the code cannot drift
    apart;
  - a valid table passes, also without row and column names, and an attribute the
    function attached is dropped;
  - a rate below 1e-7 passes with exactly one warning;
  - R's random-number state is put back after a call, and R's kind of random numbers
    after reading a file that switches it;
  - a file that defines only a dot name lists that name in the stop;
  - the round keeper records two warnings that come back in every round once each.

The CI runs no R, so these tests are skipped unless an Rscript is named, for example
the one of MetaFlux's DADA2 environment:
    METAFLUX_RSCRIPT=/path/to/the/dada2/env/bin/Rscript pytest tests
"""
import os
import re
import subprocess

import pytest

from helpers import REPO_ROOT, SCRIPTS, TEST_DIR

RSCRIPT = os.environ.get("METAFLUX_RSCRIPT", "")
pytestmark = pytest.mark.skipif(not RSCRIPT, reason="set METAFLUX_RSCRIPT to run the R checks")

# Each broken table, and the fragment of its row in the troubleshooting table ("..."
# stands for the parts that depend on the table, such as the counts of cells).
EXPECTED_STOPS = {
    "data_frame": "returned an object of class ..., not a numeric matrix",
    "char_matrix": "returned a matrix of type ..., not a numeric matrix",
    "wrong_shape": "returned a matrix of ... rows and ... columns; it must have the shape of the counts",
    "other_rows": "returned a table whose rows are named ...",
    "other_cols": "returned a table whose columns are named",
    "missing": "returned missing or infinite values",
    "above_one": "returned values below 0 or above 1",
    "zero": "returned a rate of exactly 0",
    "rows_not_one": "returned rows ... that do not sum to 1",
}


@pytest.fixture(scope="module")
def outcomes():
    """Run the R script once for all tests of this file; return {case: message or OK}."""
    script = TEST_DIR / "custom_error_functions" / "custom_function_checks.R"
    result = subprocess.run([RSCRIPT, "--vanilla", str(script), str(SCRIPTS / "amplicon")],
                            capture_output=True, text=True, check=True)
    lines = result.stdout.splitlines()
    return dict(line.split("\t", 1) for line in lines)


def fragment_pattern(fragment):
    """The fragment as a regular expression: each "..." may stand for any text."""
    pieces = [re.escape(piece.strip()) for piece in fragment.split("...")]
    return ".*".join(pieces)


def test_valid_tables_pass_and_attributes_are_dropped(outcomes):
    assert outcomes["valid"] == "OK"
    assert outcomes["no_names"] == "OK"
    assert outcomes["attribute_dropped"] == "OK"


def test_each_broken_table_stops_with_the_documented_reason(outcomes):
    troubleshooting = (REPO_ROOT / "docs" / "troubleshooting.md").read_text()
    for case, fragment in EXPECTED_STOPS.items():
        assert re.search(fragment_pattern(fragment), outcomes[case]), (case, outcomes[case])
        # The columns share the row of the rows ("returned a table whose rows are named
        # ... (or columns)"), so only the other fragments are looked up word for word.
        if case != "other_cols":
            assert f"| `{fragment}`" in troubleshooting, fragment


def test_a_rate_below_the_floor_of_dada2_passes_with_one_warning(outcomes):
    outcome = outcomes["tiny_rate"]
    assert outcome.startswith("passed with 1 warning(s): ")
    assert "rates below 1e-7, the lowest rate DADA2's own error functions return" in outcome


def test_the_random_number_state_and_kind_are_put_back(outcomes):
    assert outcomes["rng_guard"] == "OK"
    assert outcomes["rng_kind"] == "OK"


def test_a_dot_name_is_listed_when_the_name_is_not_found(outcomes):
    assert outcomes["dot_name"].startswith("does not define my_errfun (it defines: .my_errfun)")


def test_the_round_keeper_writes_each_distinct_warning_once(outcomes):
    # Three calls (rounds 0, 1, 2), two different texts at every call: two entries,
    # both from round 0, not six.
    outcome = outcomes["keeper_warnings"]
    assert outcome.startswith("3 call(s), 2 warning(s) recorded: "), outcome
    assert "round 0: the custom fit warned: A2C: fit did not reach full convergence" in outcome
    assert "round 0: the custom fit warned: A2G: fit did not reach full convergence" in outcome
