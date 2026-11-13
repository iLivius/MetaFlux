"""The choice of error function for DADA2's error learning (rule choose_error_model,
workflow/scripts/amplicon/60c_choose_error_model.py).

What is tested, and why it matters:
  - The class of the quality values decides the function under auto: continuous
    values (a MiSeq writes some 30) get DADA2's loess fit, binned values (a NovaSeq or
    NextSeq writes 3 to 6) get the binned fit, anchored on the values the reads carry,
    Q2 and below included when they survive the filtering. Picking the wrong one would
    change every ASV table of a run without a word.
  - Only samples with at least 1,000 filtered reads vote: a blank with a few dozen
    reads can lack values by chance, so its values are neither used for the choice nor
    handed to the learning as anchors, and it never stops a run (a note says what it
    carries). When no sample has 1,000 reads, every sample votes.
  - The stops, each with a text that says how to go on: no usable quality information
    (C2, in every mode; files re-encoded by ENA carry Q30 only), one such sample among
    good ones (C3, also with its own wording when no voter has quality information but
    their values together reach three), mixed instruments (C5) or two binned tables (C6)
    or an undetermined kind of values (C4) under auto, and binned set on continuous
    values (C7). loess set in the config turns C4-C6 into warnings and goes on, and
    loess set on binned values warns (W1) and goes on; binned does the same for C4 and
    C6, but on binned and continuous voters together their values make the run
    continuous, so it stops with C7, and before any "going on" warning. C5 offers loess
    only, says where the line goes, and, like C6, points to the Falco reports of the
    filtered reads; a stop text's table of samples, cut at 20 rows, names the log with
    its absolute path; the verbs agree with the number of samples named.
  - Neighbouring values: 3-10 values above Q2 that include two neighbours (29 and 30)
    count as undetermined, not binned. No instrument table has neighbouring values,
    but MiSeq reads under a trunc_q of 28 to 30 keep only the 8 to 10 values above it,
    all neighbours, which the count alone would read as binned: the binned fit would be
    taken without a word. Under auto they stop the run with C4 in its own form, which
    names trunc_q when it is above 2 (and the likely causes otherwise); with a function
    set in the config they give the matching warning. Every instrument table of the
    documentation stays binned, and neighbours that only two binned tables make
    together (a NovaSeq X writes 40 on some software versions and 41 on others) keep
    the run binned: C6 under auto, its warning with a function set, as for any two
    binned tables.
  - Both read directions are decided before anything stops; a stop text appears once
    (in the log and as the SystemExit text, not on the terminal before it), and the
    outputs are byte-identical between two runs, so a rerun that changes nothing makes
    no later step run again.
  - error_model: custom (the user's own function, experimental) gives "custom" for both
    directions, decided by the setting, whatever the quality values. The stops that
    describe the data still apply (C2, C3: no function can learn from made-up quality
    values); those that exist only because auto chooses become the warnings any
    function set in the config gets (C4, C5, C6); the rule for a function set on the
    other kind of values (C7, W1) and the note on small samples (W4, about the choice
    and the anchors of the binned fit) do not apply.

The script is run whole with a stand-in snakemake object (run_script in helpers.py) on
quality_values.json files, one per sample as dada_filter writes them, in a temporary
folder laid out like an output folder; classify() is also loaded on its own
(load_functions). That dada_filter itself writes such files is checked with R by
tests/test_error_model_r.py.
"""
import json

import pytest

from helpers import SCRIPTS, fake_snakemake, load_functions, run_script, write_text

CHOOSE_SCRIPT = SCRIPTS / "amplicon" / "60c_choose_error_model.py"

CONTINUOUS = list(range(7, 39))                     # 32 values: a filtered MiSeq file
BINNED = [9, 20, 34]                                # NextSeq 1000/2000, 600-cycle kit

# What a MiSeq file keeps under trunc_q 28, 29 and 30: only the values above it, up to
# 38 (filtered MiSeq reads carry every value from 7 to 38), so 10, 9 and 8 values, all
# neighbours, which the count alone would read as binned.
TRUNC_Q_CUT = {28: list(range(29, 39)), 29: list(range(30, 39)), 30: list(range(31, 39))}

# The quality values of the instrument tables in docs/amplicon/length-and-truncation.md
# (box "Binned quality scores"), Q2 included where the table lists it. None has two
# neighbouring values above Q2.
INSTRUMENT_TABLES = {
    "HiSeq 2000/2500": [6, 15, 22, 27, 33, 37, 40],
    "NextSeq 500/550": [2, 14, 21, 27, 32, 36],
    "NovaSeq 6000": [2, 11, 25, 37],
    "NovaSeq 6000, early software": [12, 23, 37],
    "NextSeq 1000/2000, 50-300-cycle kits": [12, 26, 34],
    "NextSeq 1000/2000, 600-cycle kit": [9, 20, 34],
    "NextSeq 1000/2000, XLEAP-SBS": [12, 24, 40],
    "NovaSeq X, before v1.2": [12, 20, 37],
    "NovaSeq X, v1.2": [2, 12, 24, 40],
    "NovaSeq X, later versions": [2, 9, 24, 40],
    "NovaSeq X, later versions with 41": [2, 9, 24, 41],
    "MiSeq i100": [2, 12, 24, 38],
    "MiSeq i100, later software": [2, 9, 23, 38],
}

# Two versions of one instrument's table in one run: each sample's own values are a
# binned table, but together they include neighbours (40 and 41, 23 and 24, 11 and 12).
# Three samples carry the first table, one the second.
TWO_TABLES_WITH_NEIGHBOURS_BETWEEN_THEM = {
    "NovaSeq X, 40 and 41": ([2, 9, 24, 40], [2, 9, 24, 41]),
    "MiSeq i100, two versions": ([2, 12, 24, 38], [2, 9, 23, 38]),
    "NovaSeq 6000, two versions": ([2, 11, 25, 37], [12, 23, 37]),
}


def quality_scan(folder, sample, r1, r2=None, reads=100000):
    """Write stats/dada2/<sample>.quality_values.json as dada_filter writes it (the
    quality scan of one sample), with the given quality values (R2 the same as R1
    unless given) and `reads` reads scanned; every value gets 100 bases, except 34,
    which gets 10,000. Returns the path."""
    data = {"sample": sample, "scan_limit": 100000}
    for direction, values in (("R1", r1), ("R2", r2 if r2 is not None else r1)):
        bases = {str(value): (10000 if value == 34 else 100) for value in values}
        data[direction] = {"reads_scanned": reads, "values": sorted(values),
                           "bases_per_value": bases}
    return write_text(folder / "out" / "stats" / "dada2" / f"{sample}.quality_values.json",
                      json.dumps(data))


def run_choose(folder, paths, requested="auto", trunc_q=2):
    """Run the script on the given quality_values.json files. Returns (choice, tsv
    text, log text) when it goes on; raises SystemExit when it stops."""
    out = folder / "out"
    sm = fake_snakemake(
        input={"quality_values": [str(path) for path in paths]},
        output={"choice": str(out / "stats" / "dada2" / "error_model_choice.json"),
                "quality_tsv": str(out / "stats" / "dada2" / "quality_values.tsv")},
        params={"error_model_requested": requested, "trunc_q": trunc_q},
        log=[str(out / "logs" / "choose_error_model.log")],
    )
    run_script(CHOOSE_SCRIPT, sm)
    choice = json.loads((out / "stats" / "dada2" / "error_model_choice.json").read_text())
    tsv = (out / "stats" / "dada2" / "quality_values.tsv").read_text()
    log = (out / "logs" / "choose_error_model.log").read_text()
    return choice, tsv, log


def stop_text(folder, paths, requested="auto", trunc_q=2):
    """Run the script where it must stop; return the stop text."""
    with pytest.raises(SystemExit) as stop:
        run_choose(folder, paths, requested, trunc_q)
    return str(stop.value)


# ── classify ──────────────────────────────────────────────────────────────────────

CLASSIFY = load_functions(CHOOSE_SCRIPT, [
    "NO_CALL_MAX", "MIN_QUALITY_VALUES", "MAX_BINNED_VALUES", "MIN_CONTINUOUS_VALUES",
    "first_neighbours", "classify"])


def test_classify_uses_the_measured_thresholds():
    classify = CLASSIFY.classify
    assert classify([30]) == "none"
    assert classify([3, 30]) == "none"
    assert classify(BINNED) == "binned"
    assert classify(range(11, 31, 2)) == "binned"         # 10 values, no two neighbours
    assert classify(range(11, 22)) == "undetermined"      # 11 values
    assert classify(range(11, 30)) == "undetermined"      # 19 values
    assert classify(range(11, 31)) == "continuous"        # 20 values
    assert classify(CONTINUOUS) == "continuous"


def test_classify_reads_neighbouring_values_as_undetermined():
    classify = CLASSIFY.classify
    # What MiSeq reads keep under trunc_q 28, 29 and 30 (values up to 38), and under 35.
    for values in TRUNC_Q_CUT.values():
        assert classify(values) == "undetermined"
    assert classify([36, 37, 38]) == "undetermined"       # three values, neighbours
    assert classify([37, 38]) == "none"                    # two values: no information
    # One pair of neighbours is enough; the first pair is the one named in messages.
    assert classify([9, 20, 21]) == "undetermined"
    assert CLASSIFY.first_neighbours([9, 20, 21, 30, 31]) == (20, 21)
    assert CLASSIFY.first_neighbours(BINNED) is None
    # Ten values with a gap between each pair are still binned.
    assert classify([3, 5, 7, 9, 11, 13, 15, 17, 19, 21]) == "binned"


@pytest.mark.parametrize("table", sorted(INSTRUMENT_TABLES))
def test_every_instrument_table_is_binned(table):
    values = INSTRUMENT_TABLES[table]
    above_q2 = [value for value in values if value > 2]
    assert CLASSIFY.first_neighbours(above_q2) is None
    assert CLASSIFY.classify(above_q2) == "binned"


# ── The function auto takes ───────────────────────────────────────────────────────

def test_1_continuous_values_get_loess(tmp_path):
    paths = [quality_scan(tmp_path, "S1", CONTINUOUS), quality_scan(tmp_path, "S2", CONTINUOUS)]
    choice, tsv, log = run_choose(tmp_path, paths)
    for direction in ("R1", "R2"):
        part = choice[direction]
        assert part["quality_class"] == "continuous"
        assert part["error_function"] == "loess"
        assert part["decided_by"] == "auto"
        assert part["n_values_above_q2"] == 32
        assert part["samples_voting"] == 2 and part["samples_not_voting"] == []
    assert choice["warnings"] == []
    assert choice["requested"] == "auto"


def test_2_binned_values_get_the_binned_fit_on_those_values(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED), quality_scan(tmp_path, "S2", BINNED)]
    choice, tsv, log = run_choose(tmp_path, paths)
    part = choice["R1"]
    assert part["quality_class"] == "binned"
    assert part["error_function"] == "binned"
    assert part["quality_values"] == [9, 20, 34]
    # Bases summed over the voters.
    assert part["bases_per_value"] == {"9": 200, "20": 200, "34": 20000}


def test_3_q2_that_survives_the_filtering_is_handed_to_learning(tmp_path):
    # trunc_q 0 keeps Q2 bases: Q2 is declared (an anchor) but does not count for the class.
    paths = [quality_scan(tmp_path, "S1", [2] + BINNED)]
    choice, tsv, log = run_choose(tmp_path, paths, trunc_q=0)
    part = choice["R1"]
    assert part["quality_values"] == [2, 9, 20, 34]
    assert part["values_above_q2"] == [9, 20, 34]
    assert part["n_values_above_q2"] == 3
    assert part["error_function"] == "binned"
    assert "\t2,9,20,34\tbinned\t100\n" in tsv            # 100 bases at or below Q2


# ── The stops ─────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("requested", ["auto", "loess", "binned"])
def test_4_one_quality_value_stops_in_every_mode(tmp_path, requested):
    paths = [quality_scan(tmp_path, "S1", [30]), quality_scan(tmp_path, "S2", [30])]
    text = stop_text(tmp_path, paths, requested)
    assert "fewer than three distinct quality values above Q2" in text
    assert "whatever error_model says" in text
    assert "sra-tools" in text and "Are the quality values real?" in text
    assert "trunc_q is" not in text
    # Both directions stop, in one text.
    assert "filtered R1 reads" in text and "filtered R2 reads" in text
    assert f"{tmp_path}/out/stats/falco/<sample>_R1_filtered/fastqc_report.html" in text


def test_4_a_raised_trunc_q_is_named_before_the_ena_sentence(tmp_path):
    paths = [quality_scan(tmp_path, "S1", [30])]
    text = stop_text(tmp_path, paths, trunc_q=20)
    assert "amplicon.dada2.trunc_q is 20" in text
    assert text.index("amplicon.dada2.trunc_q is 20") < text.index("Reads downloaded from ENA")


def test_5_fifteen_values_stop_under_auto_and_warn_when_set(tmp_path):
    fifteen = list(range(20, 35))
    paths = [quality_scan(tmp_path, "S1", fifteen)]
    text = stop_text(tmp_path, paths)
    assert "carry 15 distinct quality values above Q2" in text
    assert "does not guess" in text
    assert "some of them are neighbours" not in text      # the 11-19 form of C4
    assert "error_model: loess   # continuous quality values (MiSeq-like)" in text
    assert "do not add a second amplicon: block" in text
    choice, tsv, log = run_choose(tmp_path, paths, requested="loess")
    assert choice["R1"]["error_function"] == "loess"
    assert choice["R1"]["quality_class"] == "undetermined"
    assert choice["R1"]["decided_by"] == "setting"
    assert any("15 distinct quality values above Q2, neither binned nor continuous" in warning
               and "Going on with loess, as set." in warning for warning in choice["warnings"])


def test_6_binned_and_continuous_voters_stop_under_auto_and_warn_under_loess(tmp_path):
    paths = [quality_scan(tmp_path, "NOVA", BINNED), quality_scan(tmp_path, "MISEQ", CONTINUOUS)]
    text = stop_text(tmp_path, paths)
    assert "do not share one kind of quality values" in text
    # One sample of each kind: the verbs are singular.
    assert "1 sample carries binned values (9, 20, 34) and 1 carries continuous ones (32 values, 7-38)" in text
    assert "MISEQ" in text and "NOVA" in text           # the per-sample table
    # Only loess is offered (binned set on this mix stops with C7, see test_6b), with
    # where the line goes and the Falco reports of the filtered reads.
    assert "add error_model: loess under amplicon: -> dada2: -> learn_errors:" in text
    assert "error_model: binned" not in text
    assert "inside the amplicon: block it already has (do not add a second amplicon: block" in text
    assert f"{tmp_path}/out/stats/falco/<sample>_R1_filtered/fastqc_report.html" in text
    choice, tsv, log = run_choose(tmp_path, paths, requested="loess")
    assert choice["R1"]["error_function"] == "loess"
    assert any("different kinds of quality values; going on with one loess model" in warning
               for warning in choice["warnings"])


def test_6_with_several_samples_of_a_kind_the_verbs_are_plural(tmp_path):
    paths = [quality_scan(tmp_path, f"NOVA{number}", BINNED) for number in (1, 2, 3)]
    paths += [quality_scan(tmp_path, f"MISEQ{number}", CONTINUOUS) for number in (1, 2)]
    text = stop_text(tmp_path, paths)
    assert "3 samples carry binned values (9, 20, 34) and 2 carry continuous ones (32 values, 7-38)" in text


def test_6b_binned_set_on_binned_and_continuous_voters_stops_with_c7_only(tmp_path, capsys):
    # The union of the two samples' values has 32 values above Q2, so the run is
    # continuous and binned cannot pass. The stop (C7) comes before the warning that
    # would say "going on with one binned model ..., as set", which must not appear, and
    # it describes the mix instead of saying "as a MiSeq writes them".
    paths = [quality_scan(tmp_path, "NOVA", BINNED), quality_scan(tmp_path, "MISEQ", CONTINUOUS)]
    text = stop_text(tmp_path, paths, requested="binned")
    assert text.count("[choose_error_model] STOP:") == 2           # R1 and R2, nothing else
    assert "error_model is binned, but together the voting samples' filtered R1 reads carry" in text
    assert "although not every sample carries that many on its own" in text
    assert "as a MiSeq writes them" not in text
    assert "Run the groups in separate output folders, or set error_model to loess" in text
    assert "NOVA" in text and "MISEQ" in text                      # the per-sample table
    captured = capsys.readouterr()
    assert "going on with one binned model" not in captured.err
    log = (tmp_path / "out" / "logs" / "choose_error_model.log").read_text()
    assert "WARNING" not in log


def test_7_two_binned_tables_stop_under_auto_and_warn_under_binned(tmp_path):
    paths = [quality_scan(tmp_path, f"A{number}", [11, 25, 37]) for number in (1, 2, 3)]
    paths.append(quality_scan(tmp_path, "B1", [9, 25, 37]))
    text = stop_text(tmp_path, paths)
    assert "two different sets of binned quality values" in text
    # One sample named: the verb is singular.
    assert "most carry 11, 25, 37, but B1 also carries 9" in text
    assert "error_model: binned is not advised here" in text
    assert f"{tmp_path}/out/stats/falco/<sample>_R1_filtered/fastqc_report.html" in text
    choice, tsv, log = run_choose(tmp_path, paths, requested="binned")
    assert choice["R1"]["error_function"] == "binned"
    assert choice["R1"]["quality_values"] == [9, 11, 25, 37]          # the voters' union
    warning = next(warning for warning in choice["warnings"] if "two sets of binned values" in warning)
    assert "(11, 25, 37; 9, 25, 37)" in warning
    assert "the values of a group it does not draw are fitted on averaged qualities" in warning


def test_7_with_two_samples_named_the_verb_is_plural(tmp_path):
    paths = [quality_scan(tmp_path, f"A{number}", [11, 25, 37]) for number in (1, 2, 3)]
    paths += [quality_scan(tmp_path, f"B{number}", [9, 25, 37]) for number in (1, 2)]
    text = stop_text(tmp_path, paths)
    assert "most carry 11, 25, 37, but B1, B2 also carry 9" in text


def test_7b_loess_set_on_two_binned_tables_warns_twice_and_goes_on(tmp_path):
    # The general warning on two binned tables (W3, second variant) and the warning on
    # loess set on binned values (W1); the run goes on with loess.
    paths = [quality_scan(tmp_path, f"A{number}", [11, 25, 37]) for number in (1, 2, 3)]
    paths.append(quality_scan(tmp_path, "B1", [9, 25, 37]))
    choice, tsv, log = run_choose(tmp_path, paths, requested="loess")
    assert choice["R1"]["error_function"] == "loess"
    assert choice["R1"]["quality_class"] == "binned"
    r1_warnings = [warning for warning in choice["warnings"] if "WARNING: R1:" in warning]
    assert len(r1_warnings) == 2
    assert ("the voting samples carry two sets of binned values (11, 25, 37; 9, 25, 37); "
            "going on with one loess model for all of them, as set.") in r1_warnings[0]
    assert "error_model is loess, but the filtered reads carry binned quality values (9, 11, 25, 37)" in r1_warnings[1]


def test_8_a_voter_lacking_a_value_does_not_stop(tmp_path):
    # NextSeq 500/550 table (2, 14, 21, 27, 32, 36): a voter without 14 still carries 4
    # values above Q2, a subset of the common set, which is fine.
    nextseq_500 = [2, 14, 21, 27, 32, 36]
    paths = [quality_scan(tmp_path, "S1", nextseq_500), quality_scan(tmp_path, "S2", nextseq_500),
             quality_scan(tmp_path, "S3", [21, 27, 32, 36])]
    choice, tsv, log = run_choose(tmp_path, paths)
    assert choice["R1"]["error_function"] == "binned"
    assert choice["R1"]["quality_values"] == nextseq_500
    assert choice["warnings"] == []


def test_8b_on_a_three_value_table_a_voter_lacking_a_value_stops(tmp_path):
    # With only 9, 20, 34 a voter lacking 9 carries 2 values above Q2, fewer than three:
    # step 3 (C3) stops it, as it stops a rewritten file. A sample with 1,000 reads or
    # more carries every bin in practice, so this stop is not expected on real data.
    paths = [quality_scan(tmp_path, "S1", BINNED), quality_scan(tmp_path, "S2", BINNED),
             quality_scan(tmp_path, "S3", [20, 34])]
    text = stop_text(tmp_path, paths)
    assert "sample S3 (100,000 R1 reads scanned) carries fewer than three" in text


def test_9_one_voter_without_quality_information_stops(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED), quality_scan(tmp_path, "S2", BINNED),
             quality_scan(tmp_path, "ENA1", [3, 30])]
    text = stop_text(tmp_path, paths)
    assert "sample ENA1 (100,000 R1 reads scanned) carries fewer than three" in text
    assert "unlike the other voting samples (9, 20, 34)" in text
    assert f"{tmp_path}/out/stats/falco/ENA1_R1_filtered/fastqc_report.html" in text


@pytest.mark.parametrize("requested", ["auto", "loess", "binned", "custom"])
def test_9b_no_voter_with_quality_information_stops_without_naming_other_samples(tmp_path, requested):
    # Two files rewritten with different made-up values: each carries two values above
    # Q2, but together they carry three, so the run is classed binned and C2 does not
    # apply. The stop must not compare them with "the other voting samples (none)".
    paths = [quality_scan(tmp_path, "ENA_A", [3, 30]), quality_scan(tmp_path, "REWRITTEN_B", [3, 40])]
    text = stop_text(tmp_path, paths, requested)
    assert ("no voting sample carries three or more distinct quality values above Q2 in its R1 "
            "reads: ENA_A (100,000 reads scanned: 3, 30), REWRITTEN_B (100,000 reads scanned: 3, 40)") in text
    assert "(none)" not in text and "unlike the other voting samples" not in text
    assert "whatever error_model says" in text and "sra-tools" in text
    assert f"{tmp_path}/out/stats/falco/<sample>_R1_filtered/fastqc_report.html" in text


def test_9c_a_table_cut_short_names_the_log_with_its_absolute_path(tmp_path):
    # 25 samples with 15 values each stop with C4, whose table shows 20 rows and says
    # where the other 5 are: the absolute path of the log, as every path of the stop
    # texts (a relative one would be read against the MetaFlux folder).
    fifteen = list(range(20, 35))
    paths = [quality_scan(tmp_path, f"S{number:02d}", fifteen) for number in range(1, 26)]
    text = stop_text(tmp_path, paths)
    assert f"  ... and 5 more in {tmp_path}/out/logs/choose_error_model.log" in text
    assert " more in logs/choose_error_model.log" not in text


# ── Neighbouring values: continuous values cut by a raised trunc_q ────────────────

@pytest.mark.parametrize("trunc_q", [28, 29, 30])
def test_values_cut_by_a_raised_trunc_q_stop_under_auto_and_name_trunc_q(tmp_path, trunc_q):
    values = TRUNC_Q_CUT[trunc_q]
    paths = [quality_scan(tmp_path, "S1", values), quality_scan(tmp_path, "S2", values)]
    text = stop_text(tmp_path, paths, trunc_q=trunc_q)
    listed = ", ".join(str(value) for value in values)
    assert (f"the filtered R1 reads carry {len(values)} distinct quality values above Q2 "
            f"({listed}), and some of them are neighbours ({trunc_q + 1} and {trunc_q + 2})."
            in text)
    assert "only when there are gaps between them, as in every instrument table" in text
    assert (f"amplicon.dada2.trunc_q is {trunc_q}: filterAndTrim cut each read at its first "
            f"base at or below Q{trunc_q}, so only the values above Q{trunc_q} remain. Set it "
            "back to 2, the default") in text
    # The way on, where the line goes, the Falco reports and the per-sample table.
    assert "error_model: loess   # continuous quality values (MiSeq-like)" in text
    assert "do not add a second amplicon: block" in text
    assert f"{tmp_path}/out/stats/falco/<sample>_R1_filtered/fastqc_report.html" in text
    assert "  R1  S1  " in text and "  R1  S2  " in text
    # Both directions, once each; not the 11-19 form.
    assert text.count("[choose_error_model] STOP:") == 2
    assert "has seen nothing in between" not in text


def test_neighbouring_values_without_a_raised_trunc_q_name_other_causes(tmp_path):
    # Reads trimmed at Q30 before they reached MetaFlux, filtered with the default trunc_q.
    paths = [quality_scan(tmp_path, "S1", TRUNC_Q_CUT[30])]
    text = stop_text(tmp_path, paths, trunc_q=2)
    assert "and some of them are neighbours (31 and 32)" in text
    assert "Reads trimmed at a quality threshold before they reached MetaFlux" in text
    assert "amplicon.dada2.trunc_q is" not in text


@pytest.mark.parametrize("requested", ["loess", "binned", "custom"])
def test_values_cut_by_a_raised_trunc_q_warn_when_a_function_is_set(tmp_path, requested):
    paths = [quality_scan(tmp_path, "S1", TRUNC_Q_CUT[28])]
    choice, tsv, log = run_choose(tmp_path, paths, requested=requested, trunc_q=28)
    for direction in ("R1", "R2"):
        part = choice[direction]
        assert part["quality_class"] == "undetermined"
        assert part["error_function"] == requested
        assert part["decided_by"] == "setting"
        warnings = [warning for warning in choice["warnings"] if f"WARNING: {direction}:" in warning]
        # The matching warning (W2), and nothing else: no warning on loess set on binned
        # values (W1), no stop on binned set on continuous values (C7).
        assert warnings == [
            f"[choose_error_model] WARNING: {direction}: 10 distinct quality values above Q2 "
            "(29, 30, 31, 32, 33, 34, 35, 36, 37, 38), and some are neighbours (29 and 30), "
            "which binned quality values never are: neither binned nor continuous as MetaFlux "
            "knows them. amplicon.dada2.trunc_q is 28, so only the values above Q28 remain. "
            f"Going on with {requested}, as set."]
    # Each sample's own class in the table for people.
    assert "S1\tR1\t100000\tyes\t10\t29,30,31,32,33,34,35,36,37,38\tundetermined\t0\n" in tsv


@pytest.mark.parametrize("table", sorted(INSTRUMENT_TABLES))
def test_every_instrument_table_gets_the_binned_fit_without_a_warning(tmp_path, table):
    values = INSTRUMENT_TABLES[table]
    paths = [quality_scan(tmp_path, "S1", values), quality_scan(tmp_path, "S2", values)]
    choice, tsv, log = run_choose(tmp_path, paths)
    for direction in ("R1", "R2"):
        assert choice[direction]["quality_class"] == "binned"
        assert choice[direction]["error_function"] == "binned"
        assert choice[direction]["quality_values"] == values
    assert choice["warnings"] == []


@pytest.mark.parametrize("case", sorted(TWO_TABLES_WITH_NEIGHBOURS_BETWEEN_THEM))
def test_neighbours_between_two_binned_tables_keep_the_run_binned(tmp_path, case):
    # Each sample carries a binned table of its own; only together do their values
    # include neighbours. That is two binned tables (C6 under auto), not values cut by
    # trunc_q, so nothing of the neighbour rule appears.
    first_table, second_table = TWO_TABLES_WITH_NEIGHBOURS_BETWEEN_THEM[case]
    paths = [quality_scan(tmp_path, f"A{number}", first_table) for number in (1, 2, 3)]
    paths.append(quality_scan(tmp_path, "B1", second_table))
    text = stop_text(tmp_path, paths)
    assert "the voting samples carry two different sets of binned quality values" in text
    assert "some of them are neighbours" not in text
    # With a function set: the warning on two binned tables, and for loess the warning
    # on loess set on binned values; never the warning on undetermined values.
    for requested in ("loess", "binned", "custom"):
        folder = tmp_path / requested
        paths = [quality_scan(folder, f"A{number}", first_table) for number in (1, 2, 3)]
        paths.append(quality_scan(folder, "B1", second_table))
        choice, tsv, log = run_choose(folder, paths, requested=requested)
        assert choice["R1"]["quality_class"] == "binned"
        assert choice["R1"]["error_function"] == requested
        r1_warnings = [warning for warning in choice["warnings"] if "WARNING: R1:" in warning]
        assert "two sets of binned values" in r1_warnings[0]
        assert not any("neither binned nor continuous" in warning for warning in r1_warnings)
        if requested == "loess":
            assert len(r1_warnings) == 2
            assert "error_model is loess, but the filtered reads carry binned quality values" in r1_warnings[1]
        else:
            assert len(r1_warnings) == 1


# ── Small samples ─────────────────────────────────────────────────────────────────

def test_10_a_small_sample_with_an_extra_value_gets_a_note(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED), quality_scan(tmp_path, "S2", BINNED),
             quality_scan(tmp_path, "blank1", BINNED + [37], reads=500)]
    choice, tsv, log = run_choose(tmp_path, paths)
    assert choice["R1"]["quality_values"] == [9, 20, 34]
    assert choice["R1"]["samples_not_voting"] == ["blank1"]
    note = next(warning for warning in choice["warnings"] if "blank1" in warning)
    assert note.startswith("[choose_error_model] NOTE: blank1 (500 R1 reads scanned)")
    assert "carries quality values the voting samples do not (37)" in note
    assert "blank1\tR1\t500\tno\t" in tsv


def test_11_when_no_sample_has_1000_reads_every_sample_votes(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED, reads=400),
             quality_scan(tmp_path, "S2", BINNED, reads=600)]
    choice, tsv, log = run_choose(tmp_path, paths)
    assert choice["R1"]["samples_voting"] == 2
    assert any("no sample has 1,000 or more filtered R1 reads; every sample votes" in warning
               for warning in choice["warnings"])


def test_17_a_small_continuous_sample_among_binned_voters_is_ignored(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED), quality_scan(tmp_path, "S2", BINNED),
             quality_scan(tmp_path, "tiny", CONTINUOUS, reads=500)]
    choice, tsv, log = run_choose(tmp_path, paths)
    assert choice["R1"]["error_function"] == "binned"
    assert choice["R1"]["quality_values"] == [9, 20, 34]
    assert any("NOTE: tiny" in warning for warning in choice["warnings"])


def test_18_a_small_sample_with_rewritten_values_is_ignored(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED), quality_scan(tmp_path, "S2", BINNED),
             quality_scan(tmp_path, "ena_blank", [3, 30], reads=500)]
    choice, tsv, log = run_choose(tmp_path, paths)
    assert choice["R1"]["quality_values"] == [9, 20, 34]
    note = next(warning for warning in choice["warnings"] if "ena_blank" in warning)
    assert "(3, 30)" in note


def test_19_non_voters_never_stop_anything(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED), quality_scan(tmp_path, "S2", BINNED),
             quality_scan(tmp_path, "odd", list(range(20, 35)), reads=500)]
    choice, tsv, log = run_choose(tmp_path, paths)
    assert choice["R1"]["error_function"] == "binned"
    assert choice["R1"]["quality_class"] == "binned"


# ── A function set in the config ──────────────────────────────────────────────────

def test_12_binned_set_on_continuous_values_stops(tmp_path):
    paths = [quality_scan(tmp_path, "S1", CONTINUOUS)]
    text = stop_text(tmp_path, paths, requested="binned")
    assert "error_model is binned, but the filtered R1 reads carry continuous quality values" in text
    assert "(32 distinct values above Q2, 7-38)" in text
    assert "Set error_model to auto or loess" in text


def test_13_loess_set_on_binned_values_warns_and_goes_on(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED)]
    choice, tsv, log = run_choose(tmp_path, paths, requested="loess")
    assert choice["R1"]["error_function"] == "loess"
    assert choice["R1"]["decided_by"] == "setting"
    assert any("error_model is loess, but the filtered reads carry binned quality values (9, 20, 34)"
               in warning and "the fit MetaFlux used up to v2.6.1" in warning
               for warning in choice["warnings"])


# ── Mixed directions, reproducibility, output channels ────────────────────────

def test_15_r1_binned_and_r2_continuous_warn(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED, r2=CONTINUOUS)]
    choice, tsv, log = run_choose(tmp_path, paths)
    assert choice["R1"]["error_function"] == "binned"
    assert choice["R2"]["error_function"] == "loess"
    assert any("R1 carries binned and R2 continuous quality values" in warning
               for warning in choice["warnings"])


def test_16_two_runs_give_byte_identical_outputs(tmp_path):
    paths = [quality_scan(tmp_path, "S2", BINNED), quality_scan(tmp_path, "S1", BINNED),
             quality_scan(tmp_path, "blank", BINNED + [37], reads=500)]
    first = run_choose(tmp_path, paths)
    first_bytes = (tmp_path / "out" / "stats" / "dada2" / "error_model_choice.json").read_bytes()
    second = run_choose(tmp_path, list(reversed(paths)))
    second_bytes = (tmp_path / "out" / "stats" / "dada2" / "error_model_choice.json").read_bytes()
    assert first_bytes == second_bytes
    assert first[1] == second[1]                                  # the TSV
    # Samples sorted by name, values ascending, a final line end.
    assert [sample["sample"] for sample in first[0]["samples"]] == ["S1", "S2", "blank"]
    assert first_bytes.endswith(b"}\n")


def test_20_a_stop_text_appears_once(tmp_path, capsys):
    paths = [quality_scan(tmp_path, "S1", [30])]
    text = stop_text(tmp_path, paths)
    captured = capsys.readouterr()
    assert "STOP" not in captured.err and "STOP" not in captured.out
    log = (tmp_path / "out" / "logs" / "choose_error_model.log").read_text()
    assert text in log
    assert log.count("[choose_error_model] STOP:") == text.count("[choose_error_model] STOP:")


def test_the_quality_tsv_has_one_row_per_sample_and_direction(tmp_path):
    paths = [quality_scan(tmp_path, "S1", BINNED, r2=[2] + BINNED)]
    choice, tsv, log = run_choose(tmp_path, paths, trunc_q=0)
    lines = tsv.rstrip("\n").split("\n")
    assert lines[0] == ("sample\tread_direction\treads_scanned\tvotes\tn_values_above_q2\tvalues"
                        "\tquality_class\tbases_at_or_below_q2")
    assert lines[1] == "S1\tR1\t100000\tyes\t3\t9,20,34\tbinned\t0"
    assert lines[2] == "S1\tR2\t100000\tyes\t3\t2,9,20,34\tbinned\t100"


# ── error_model: custom, the user's own function (experimental) ──────────────────

@pytest.mark.parametrize("values", [CONTINUOUS, BINNED, list(range(20, 35))])
def test_custom_is_used_for_both_directions_whatever_the_values(tmp_path, values):
    # Continuous, binned and 11-19 values alike: the function is set in the config.
    paths = [quality_scan(tmp_path, "S1", values), quality_scan(tmp_path, "S2", values)]
    choice, tsv, log = run_choose(tmp_path, paths, requested="custom")
    assert choice["requested"] == "custom"
    for direction in ("R1", "R2"):
        part = choice[direction]
        assert part["error_function"] == "custom"
        assert part["decided_by"] == "setting"
        # The class and the values are still read and written down for the report.
        assert part["quality_values"] == sorted(values)
    assert ("R1: " in log and "-> your own error function (custom, experimental), "
            "as set in error_model" in log)


def test_custom_on_continuous_or_binned_values_gives_no_warning(tmp_path):
    # No C7 stop and no W1 warning: the user's function is made for no particular kind.
    for values in (CONTINUOUS, BINNED):
        paths = [quality_scan(tmp_path, "S1", values)]
        choice, tsv, log = run_choose(tmp_path, paths, requested="custom")
        assert choice["warnings"] == []


def test_custom_keeps_the_stop_for_reads_without_quality_information(tmp_path):
    # C2: no function can learn a rate per quality score from one made-up value.
    paths = [quality_scan(tmp_path, "S1", [30]), quality_scan(tmp_path, "S2", [30])]
    text = stop_text(tmp_path, paths, requested="custom")
    assert "fewer than three distinct quality values above Q2" in text
    assert "whatever error_model says" in text


def test_custom_keeps_the_stop_for_one_sample_without_quality_information(tmp_path):
    # C3: one sample with rewritten values would make the one model wrong for every sample.
    paths = [quality_scan(tmp_path, "S1", BINNED), quality_scan(tmp_path, "S2", BINNED),
             quality_scan(tmp_path, "ENA1", [3, 30])]
    text = stop_text(tmp_path, paths, requested="custom")
    assert "sample ENA1 (100,000 R1 reads scanned) carries fewer than three" in text


def test_custom_turns_the_stops_of_auto_into_warnings(tmp_path):
    # C5: binned and continuous voters.
    paths = [quality_scan(tmp_path, "NOVA", BINNED), quality_scan(tmp_path, "MISEQ", CONTINUOUS)]
    choice, tsv, log = run_choose(tmp_path, paths, requested="custom")
    assert choice["R1"]["error_function"] == "custom"
    assert any("different kinds of quality values; going on with one custom model for all of "
               "them, as set." in warning for warning in choice["warnings"])
    # C6: two binned tables (the general warning, not the binned fit's caveat).
    folder = tmp_path / "two_tables"
    paths = [quality_scan(folder, f"A{number}", [11, 25, 37]) for number in (1, 2, 3)]
    paths.append(quality_scan(folder, "B1", [9, 25, 37]))
    choice, tsv, log = run_choose(folder, paths, requested="custom")
    warning = next(warning for warning in choice["warnings"] if "two sets of binned values" in warning)
    assert "going on with one custom model for all of them, as set." in warning
    assert "averaged qualities" not in warning
    # C4: 11-19 values.
    folder = tmp_path / "fifteen"
    paths = [quality_scan(folder, "S1", list(range(20, 35)))]
    choice, tsv, log = run_choose(folder, paths, requested="custom")
    assert any("15 distinct quality values above Q2, neither binned nor continuous" in warning
               and "Going on with custom, as set." in warning for warning in choice["warnings"])


def test_custom_leaves_out_the_note_on_small_samples(tmp_path):
    # W4 says the small sample's values are used neither for the choice nor as anchors
    # of the binned fit; under custom there is neither, so the note would mislead.
    paths = [quality_scan(tmp_path, "S1", BINNED), quality_scan(tmp_path, "S2", BINNED),
             quality_scan(tmp_path, "blank1", BINNED + [37], reads=500)]
    choice, tsv, log = run_choose(tmp_path, paths, requested="custom")
    assert not any("blank1" in warning for warning in choice["warnings"])
    assert choice["R1"]["samples_not_voting"] == ["blank1"]


def test_custom_keeps_the_warning_on_r1_and_r2_of_different_kinds(tmp_path):
    # W6 describes the input files, whatever function learns from them.
    paths = [quality_scan(tmp_path, "S1", BINNED, r2=CONTINUOUS)]
    choice, tsv, log = run_choose(tmp_path, paths, requested="custom")
    assert choice["R1"]["error_function"] == choice["R2"]["error_function"] == "custom"
    assert any("R1 carries binned and R2 continuous quality values" in warning
               for warning in choice["warnings"])
