"""R checks of DADA2's error learning that need R but no real reads (the scripts
60b_dada_filter.R, 60d_learn_errors.R, 60g_dada_seqtab.R, 60e_error_model_report.R and
the registry 60_error_functions.R), against what the code comments and the
documentation promise.

What is tested, and why it matters:
  - learn_errors sorts the filtered files byte by byte, so the samples learnErrors draws,
    and with them the model and its code, do not depend on the language setting (locale)
    of the shell: a relearn from a cluster node or a container under LANG/LC_ALL=C must
    not give another model code than one from a desktop shell. Each R2 file is named
    after its own sample, also where the R2 list sorts in another order than the R1 list.
  - dada_seqtab pairs each R1 file with the R2 file of the same name. Paired by position,
    as up to v2.6.1, samples named "Soil" and "Soil_R2" would have had their R2 files
    swapped under LANG/LC_ALL=C: DADA2 stops with an error when the two samples have
    different read counts, and merges the wrong reads without one when the counts are
    equal. A missing R2 file stops the job with a text in its log (this case needs the R
    packages pacman and dada2, and is skipped without them).
  - The flat-rounds line of the learning log says "over the last round" for a window of
    one round, and describes DADA2's convergence line instead of quoting it (DADA2 writes
    it with two spaces around the number, so a quote would not be found in the log).
  - The reason the binned fit stops on quality scores outside its anchors, and the
    report's note for a refit with added anchors, count more than ten such scores on each
    side of the anchors ("13 below 20 (7-19) and 8 above 30 (31-38)"); one range over both
    sides would span the anchors and read as if every score in it lay outside them.
  - The report's figure titles and MultiQC summary shorten a list of more than ten anchors
    to "16 values (7-22)", so that the outcome and the model code at the end of a title are
    not cut off at the right edge of the picture; three or ten anchors are listed in full,
    as before, and the goodness-of-fit table is unchanged.
  - dada_filter, run whole on a made-up read pair with binned quality values (9, 20,
    34), writes two files: filter_stats.json with the read counts, byte for byte in the
    layout of v2.6.1 (dada_seqtab reads it for the read tracking), and
    quality_values.json with the quality scan; choose_error_model, run on that file as
    in a run, takes the binned fit on 9, 20 and 34. This joins the two sides of the
    file: the R script that writes it and the Python script that reads it (this case
    needs the R packages pacman, dada2, jsonlite and ShortRead).

How: the test runs .test/r_checks/error_model_checks.R once with Rscript. That script
loads only the function definitions of three scripts (and runs the few statements of
60d that sort the files on made-up paths), runs 60g up to a stop and 60b whole, feeds
them made-up input, and prints one line per case, "<case><TAB>OK", "SKIP: ..." or what
it found instead. 60b writes its files into a folder of this test (filter_folder).
Rscript runs without the variables LC_ALL and LC_COLLATE (see the fixture below for
why).

The CI runs no R, so these tests are skipped unless an Rscript is named, for example the
one of MetaFlux's DADA2 environment:
    METAFLUX_RSCRIPT=/path/to/the/dada2/env/bin/Rscript pytest tests
"""
import json
import os
import subprocess

import pytest

from helpers import SCRIPTS, TEST_DIR, fake_snakemake, run_script

RSCRIPT = os.environ.get("METAFLUX_RSCRIPT", "")
pytestmark = pytest.mark.skipif(not RSCRIPT, reason="set METAFLUX_RSCRIPT to run the R checks")


@pytest.fixture(scope="module")
def filter_folder(tmp_path_factory):
    """The empty folder where the dada_filter case of the R script writes its files."""
    return tmp_path_factory.mktemp("dada_filter")


@pytest.fixture(scope="module")
def outcomes(filter_folder):
    """Run the R script once for all tests of this file; return {case: OK or text}."""
    script = TEST_DIR / "r_checks" / "error_model_checks.R"
    # Rscript runs without LC_ALL and LC_COLLATE. When either of them is C in R's
    # environment (common on cluster nodes, in containers and in cron jobs), R stops
    # using its own (ICU) alphabetical order, even after the locale case switches to a
    # UTF-8 locale with Sys.setlocale(): that case's sample names then sort alike under
    # both locales, and the case could show nothing (it would report SKIP). R's help
    # page ?icuSetCollate gives this rule. LANG and the other variables are kept;
    # LANG=C does not have this effect.
    r_environment = dict(os.environ)
    r_environment.pop("LC_ALL", None)
    r_environment.pop("LC_COLLATE", None)
    result = subprocess.run([RSCRIPT, "--vanilla", str(script), str(SCRIPTS / "amplicon"),
                             str(filter_folder)],
                            capture_output=True, text=True, check=True, env=r_environment)
    lines = result.stdout.splitlines()
    return dict(line.split("\t", 1) for line in lines)


def passed(outcomes, case):
    """Assert that a case printed OK; skip the test when the case could not run here."""
    outcome = outcomes[case]
    if outcome.startswith("SKIP:"):
        pytest.skip(f"{case}: {outcome}")
    assert outcome == "OK", (case, outcome)


def test_learn_errors_sorts_the_files_the_same_way_in_every_locale(outcomes):
    passed(outcomes, "sort_same_in_every_locale")
    passed(outcomes, "r2_files_named_after_their_own_sample")


def test_dada_seqtab_pairs_r1_and_r2_files_by_name(outcomes):
    passed(outcomes, "seqtab_pairs_by_name")
    passed(outcomes, "seqtab_missing_partner_is_na")


def test_dada_seqtab_stops_when_an_r2_file_is_missing(outcomes):
    passed(outcomes, "seqtab_stops_without_partner")


@pytest.mark.parametrize("window", [1, 2])
def test_the_flat_rounds_line_names_its_window_and_describes_dada2s_line(outcomes, window):
    passed(outcomes, f"flat_line_window_of_{window}")


@pytest.mark.parametrize("case", [
    "out_of_range_both_sides", "out_of_range_below_only", "out_of_range_above_only",
    "out_of_range_single_score_on_one_side", "out_of_range_ten_or_fewer_listed",
    "binned_fit_stop_names_both_sides", "extended_note_names_both_sides",
])
def test_scores_outside_the_anchors_are_counted_on_each_side(outcomes, case):
    passed(outcomes, case)


@pytest.mark.parametrize("case", [
    "title_shortens_more_than_ten_anchors", "title_lists_three_anchors",
    "title_lists_ten_anchors", "summary_shortens_more_than_ten_anchors",
    "goodness_of_fit_anchors_unchanged",
])
def test_the_report_shortens_long_anchor_lists_only(outcomes, case):
    passed(outcomes, case)


def test_dada_filter_writes_the_read_counts_and_the_quality_values_to_two_files(outcomes):
    passed(outcomes, "dada_filter_writes_two_files")


def test_choose_error_model_reads_the_file_dada_filter_wrote(outcomes, filter_folder, tmp_path):
    # The file 60b wrote in the R case above, read by 60c as in a run: one sample whose
    # reads carry 9, 20 and 34 gets the binned fit on those values. 60 reads are fewer
    # than the 1,000 a sample needs to vote, so every sample votes, with a warning.
    passed(outcomes, "dada_filter_writes_two_files")
    quality_file = filter_folder / "stats" / "dada2" / "S1.quality_values.json"
    out = tmp_path / "out"
    snakemake = fake_snakemake(
        input={"quality_values": [str(quality_file)]},
        output={"choice": str(out / "stats" / "dada2" / "error_model_choice.json"),
                "quality_tsv": str(out / "stats" / "dada2" / "quality_values.tsv")},
        params={"error_model_requested": "auto", "trunc_q": 2},
        log=[str(out / "logs" / "choose_error_model.log")])
    run_script(SCRIPTS / "amplicon" / "60c_choose_error_model.py", snakemake)
    choice = json.loads((out / "stats" / "dada2" / "error_model_choice.json").read_text())
    for direction in ("R1", "R2"):
        assert choice[direction]["quality_class"] == "binned"
        assert choice[direction]["error_function"] == "binned"
        assert choice[direction]["quality_values"] == [9, 20, 34]
        assert choice[direction]["bases_per_value"] == {"9": 300, "20": 900, "34": 6000}
    assert choice["samples"] == [
        {"sample": "S1",
         "R1": {"reads_scanned": 60, "votes": True, "quality_values": [9, 20, 34]},
         "R2": {"reads_scanned": 60, "votes": True, "quality_values": [9, 20, 34]}}]
    assert any("no sample has 1,000 or more filtered R1 reads; every sample votes" in warning
               for warning in choice["warnings"])
