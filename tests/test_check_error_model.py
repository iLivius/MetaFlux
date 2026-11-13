"""The check between learning and denoising (rule check_error_model,
workflow/scripts/amplicon/60f_check_error_model.py).

What is tested, and why it matters:
  - Both models converged: the verdict is written and the run goes on, as before.
  - A model that did not converge (flat-rounds rule or cap) and whose code is not in
    amplicon.dada2.learn_errors.accept_unconverged stops the run before denoising,
    with a text that sends the user to the reads first and then to the model, gives
    absolute paths (the user's shell reads a relative one against the MetaFlux folder),
    prints the exact line that accepts this model, and warns against pasting a second
    amplicon: block (YAML keeps only the last of two). The text appears once on the
    terminal (as the text of SystemExit) and in full in the log; no verdict is written,
    so dada_seqtab cannot run.
  - A listed code accepts exactly that model: the run goes on with a warning, a code of
    the other read direction accepts nothing, and listed codes that name no model of
    the run get a note without changing the verdict. A code that reaches the check as
    written in the run config (one code without brackets, other capitals) still accepts
    its model, should the normalised list of the start-up check ever be missing.
  - A cap below the default of 200 gets its own first paragraph (raising it is then the
    likely fix); the orientation line appears only under orientation: mixed.
  - A cap of 1 or 2 rounds has its own wording of how the model was still changing
    ("in its only round", "in its last round, round 2, against ... in round 1"); the
    forced stop with max_consist: 1 meets the first. A flat-rounds rule with a window of
    one round (max_flat_rounds: 1) says "over the last round", not "the last 1 rounds".
  - A config without the key accepts nothing, and the verdict is byte-identical between
    two runs, so a rerun of the check does not count as a change for the next run.
  - The warnings of the error-model report (error_model_report.json, written by
    60e_error_model_report.R: a figure it could not draw, for example) are repeated as a
    note when the run goes on, without changing the verdict.

The script is run whole with a stand-in snakemake object (run_script in helpers.py) on
a learning record written to a temporary folder laid out like an output folder.
"""
import json

import pytest

from helpers import SCRIPTS, fake_snakemake, run_script, write_text

CHECK_SCRIPT = SCRIPTS / "amplicon" / "60f_check_error_model.py"

# The change per round of a real R2 learning run (a binned NextSeq 2000 16S library,
# loess fit) up to round 42, where the flat-rounds rule with N = 20 ends it.
FLAT_CHANGES = [
    41.5, 1.53, 0.209, 0.0905, 0.0345, 0.0157, 0.013, 0.0117, 0.0054, 0.0111, 0.0124,
    0.000877, 0.00449, 0.0127, 0.00859, 0.00498, 0.00629, 0.00566, 0.00434, 0.00236,
    0.00981, 0.000386, 0.00693, 0.00686, 0.0151, 0.00512, 0.00244, 0.00903, 0.00191,
    0.00886, 0.00987, 0.00885, 0.00479, 0.0068, 0.015, 0.0121, 0.0117, 0.00987, 0.00988,
    0.00727, 0.00514, 0.0097,
]


def converged_part(code):
    """The record entry of a model that converged after 4 rounds (a plain repeat)."""
    return {"error_function": "loess", "rounds": 4, "converged": True,
            "stopped_by": "dada2_repeat", "matched_round": 4, "cycle_length": 1,
            "flat_stop_round": None, "change_per_round": [27.7, 0.72, 0.0824, 0.0],
            "model_code": code}


def flat_part(code):
    """The record entry of a model ended by the flat-rounds rule in round 42."""
    return {"error_function": "loess", "rounds": 42, "converged": False,
            "stopped_by": "flat_rounds", "matched_round": None, "cycle_length": None,
            "flat_stop_round": 42, "change_per_round": FLAT_CHANGES, "model_code": code}


def cap_part(code, cap):
    """The record entry of a model stopped at a cap of `cap` rounds."""
    changes = [48.5, 3.05, 0.254, 0.0602, 0.0267, 0.0174, 0.0033, 0.00192, 0.00422, 0.00302]
    return {"error_function": "loess", "rounds": cap, "converged": False,
            "stopped_by": "cap", "matched_round": None, "cycle_length": None,
            "flat_stop_round": None, "change_per_round": changes[:cap], "model_code": code}


def run_check(tmp_path, r1, r2, accepted=None, orientation="fixed", max_consist=200,
              config=None, report_warnings=None, max_flat_rounds=20):
    """Write a learning record and the report JSON, and run the check on them, as rule
    check_error_model would.

    Returns (exit, out, verdict_path, log_path): `exit` is the SystemExit raised by a
    stop (None when the run goes on), `out` the output folder."""
    out = tmp_path / "out"
    record = {"dada2_version": "1.34.0", "requested": "loess", "nbases": 1e8,
              "max_consist": max_consist, "seed": 42, "max_flat_rounds": max_flat_rounds,
              "flat_rule": {"depth_factor": 1000, "min_gain": 0.1},
              "R1": r1, "R2": r2, "warnings": []}
    record_path = write_text(out / "stats" / "dada2" / "error_model.json", json.dumps(record))
    report = {"statement": "...", "warnings": report_warnings or []}
    report_path = write_text(out / "stats" / "dada2" / "error_model_report.json", json.dumps(report))
    verdict_path = out / "stats" / "dada2" / "error_model_check.txt"
    log_path = out / "logs" / "check_error_model.log"
    if config is None:
        config = {"amplicon": {"dada2": {"learn_errors": {"accept_unconverged": accepted or []}}}}
    sm = fake_snakemake(
        input={"record": record_path, "report": report_path},
        output={"verdict": str(verdict_path)},
        params={"orientation": orientation},
        log=[str(log_path)],
        config=config,
    )
    try:
        run_script(CHECK_SCRIPT, sm)
    except SystemExit as stop:
        return stop, out, verdict_path, log_path
    return None, out, verdict_path, log_path


def flat(text):
    """The text with every run of spaces and line ends made one space, so that a
    phrase can be found wherever the text breaks its lines."""
    return " ".join(text.split())


def test_both_converged_writes_the_verdict(tmp_path, capsys):
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        converged_part("R2:6544a5c3"))
    assert stop is None
    assert verdict.read_text() == "R1:0a55be94 converged; R2:6544a5c3 converged\n"
    assert capsys.readouterr().err == ""


def test_an_unaccepted_flat_stop_stops_the_run_with_the_whole_text(tmp_path, capsys):
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        flat_part("R2:7ebe3bed"))
    assert stop is not None
    text = str(stop.code)
    lines = text.splitlines()
    assert "STOP" in lines[0] and "R2:7ebe3bed" in lines[0]
    # The single line to paste comes before the nested blocks.
    single = text.index('    accept_unconverged: ["R2:7ebe3bed"]')
    nested = text.index("    amplicon:\n      dada2:\n        learn_errors:\n"
                        '          accept_unconverged: ["R2:7ebe3bed"]')
    assert single < nested
    words = flat(text)
    assert "do not paste a second amplicon: block" in words
    assert f"{out}/stats/dada2/error_model_fit_R2.png" in text
    assert f"{out}/stats/dada2/error_model_convergence.png" in text
    assert f"{out}/stats/dada2/error_model_quality_values.png" in text
    assert "the solid line should follow the filled points" in words
    # The page that explains the figures, and the troubleshooting entry (both must exist:
    # a strict docs build checks the page, this test the pointer to it).
    assert "docs/amplicon/error-model.md" in words
    assert "\"DADA2's error model never converges\" in docs/troubleshooting.md" in words
    assert "error_model_fit_R1.png" not in text            # R1 converged: nothing to look at
    assert "round 42" in words
    assert "the lowest change so far fell by less than 10 %" in words
    assert "41.5" in words and "(round 2: 1.53)" in words and "changed by 0.0097" in words
    assert "0.04 %" in words
    assert "R1:0a55be94" not in text.split("To go on")[1]   # the converged model needs no code
    # The first and the last line say where the text is kept.
    assert str(log) in lines[-1]
    # The log holds the same text; nothing went to the terminal before the exit, and
    # no verdict was written.
    assert log.read_text() == text + "\n"
    assert capsys.readouterr().err == ""
    assert not verdict.exists()


def test_a_listed_code_accepts_that_model(tmp_path, capsys):
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        flat_part("R2:7ebe3bed"), accepted=["R2:7ebe3bed"])
    assert stop is None
    assert verdict.read_text() == (
        "R1:0a55be94 converged; R2:7ebe3bed accepted (not converged: flat-rounds rule, round 42)\n")
    err = capsys.readouterr().err
    assert "WARNING: going on with R2:7ebe3bed" in err
    assert "flat-rounds rule, round 42" in err


def test_a_code_as_written_in_the_run_config_still_accepts_its_model(tmp_path, capsys):
    # 00_common.smk normalises the list at start-up and writes it back into the config
    # the check reads. Should that list ever arrive as written in the run config, the
    # check normalises it itself: a single code without brackets is one code (not eight
    # characters), and other capitals name the same model.
    single_text = {"amplicon": {"dada2": {"learn_errors": {"accept_unconverged": " r2:7EBE3BED "}}}}
    stop, out, verdict, log = run_check(tmp_path / "text", converged_part("R1:0a55be94"),
                                        flat_part("R2:7ebe3bed"), config=single_text)
    assert stop is None
    assert verdict.read_text().endswith("R2:7ebe3bed accepted (not converged: flat-rounds rule, round 42)\n")
    stop, out, verdict, log = run_check(tmp_path / "capitals", converged_part("R1:0a55be94"),
                                        flat_part("R2:7ebe3bed"), accepted=["r2:7EBE3BED"])
    assert stop is None
    assert "WARNING: going on with R2:7ebe3bed" in capsys.readouterr().err


def test_a_code_of_the_other_direction_accepts_nothing(tmp_path):
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        flat_part("R2:7ebe3bed"), accepted=["R1:7ebe3bed"])
    assert stop is not None
    assert not verdict.exists()
    # The current list is echoed at the end of the text.
    assert '(accept_unconverged in this run: ["R1:7ebe3bed"].)' in str(stop.code)


def test_a_stale_code_gets_a_note_and_changes_nothing(tmp_path, capsys):
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        converged_part("R2:6544a5c3"), accepted=["R2:12345678"])
    assert stop is None
    assert verdict.read_text() == "R1:0a55be94 converged; R2:6544a5c3 converged\n"
    assert "NOTE: accept_unconverged also lists R2:12345678" in capsys.readouterr().err


def test_a_low_cap_comes_first_and_has_no_flat_rounds_sentence(tmp_path):
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        cap_part("R2:26dab1e9", 10), max_consist=10)
    text = str(stop.code)
    paragraphs = text.split("\n\n")
    assert paragraphs[1].startswith("Your limit of 10 rounds")
    words = flat(text)
    assert "used all 10 rounds it is allowed" in words
    assert "by 0.00302 in its last round, against 48.5 and 3.05 in rounds 1 and 2" in words
    assert "0.04 %" not in words


def test_a_cap_of_one_round_has_its_own_wording(tmp_path):
    # max_consist: 1, the forced stop of the validation: one round, so "the one round it
    # is allowed" and "in its only round" (no round 2 to compare with).
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        cap_part("R2:1dd0b484", 1), max_consist=1)
    words = flat(str(stop.code))
    assert "learnErrors used the one round it is allowed (amplicon.dada2.learn_errors.max_consist)" in words
    assert "the model was still changing: by 48.5 in its only round." in words
    assert words.count("Your limit of 1 round (amplicon.dada2.learn_errors.max_consist)") == 1
    assert not verdict.exists()


def test_a_cap_of_two_rounds_names_both_rounds(tmp_path):
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        cap_part("R2:a5dd15a5", 2), max_consist=2)
    words = flat(str(stop.code))
    assert "learnErrors used all 2 rounds it is allowed" in words
    assert "by 3.05 in its last round, round 2, against 48.5 in round 1." in words
    assert "Your limit of 2 rounds" in words


def test_a_flat_stop_with_a_window_of_one_round_says_the_last_round(tmp_path):
    # max_flat_rounds: 1. The rule's window is one round: "over the last round", not
    # "over the last 1 rounds".
    part = flat_part("R2:7ebe3bed")
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"), part,
                                        max_flat_rounds=1)
    words = flat(str(stop.code))
    assert "(amplicon.dada2.learn_errors.max_flat_rounds: 1)" in words
    assert "fell by less than 10 % over the last round." in words
    assert "last 1 rounds" not in words
    # With the default window the number is given.
    stop = run_check(tmp_path / "default", converged_part("R1:0a55be94"), part)[0]
    assert "fell by less than 10 % over the last 20 rounds." in flat(str(stop.code))


def test_both_directions_unaccepted_need_both_codes(tmp_path):
    stop, out, verdict, log = run_check(tmp_path, cap_part("R1:26dab1e9", 10),
                                        flat_part("R2:7ebe3bed"), max_consist=10)
    text = str(stop.code)
    assert "R1 and R2 reads" in text.splitlines()[0]
    assert '    accept_unconverged: ["R1:26dab1e9", "R2:7ebe3bed"]' in text
    assert '          accept_unconverged: ["R1:26dab1e9", "R2:7ebe3bed"]' in text


def test_the_orientation_line_appears_under_mixed_only(tmp_path):
    fixed_stop = run_check(tmp_path / "fixed", converged_part("R1:0a55be94"),
                           flat_part("R2:7ebe3bed"), orientation="fixed")[0]
    mixed_stop = run_check(tmp_path / "mixed", converged_part("R1:0a55be94"),
                           flat_part("R2:7ebe3bed"), orientation="mixed")[0]
    assert "/stats/orientation/" not in str(fixed_stop.code)
    assert "/stats/orientation/ how many pairs each pass kept" in flat(str(mixed_stop.code))


def test_a_config_without_the_key_accepts_nothing(tmp_path):
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        flat_part("R2:7ebe3bed"), config={})
    assert stop is not None
    assert "(accept_unconverged in this run: [].)" in str(stop.code)


def test_the_verdict_is_byte_identical_between_two_runs(tmp_path):
    first = run_check(tmp_path / "a", converged_part("R1:0a55be94"),
                      flat_part("R2:7ebe3bed"), accepted=["R2:7ebe3bed"])[2].read_bytes()
    second = run_check(tmp_path / "b", converged_part("R1:0a55be94"),
                       flat_part("R2:7ebe3bed"), accepted=["R2:7ebe3bed"])[2].read_bytes()
    assert first == second


def test_the_warnings_of_the_report_become_a_note(tmp_path, capsys):
    warnings = ["R2: the fit figure of the R2 reads could not be drawn: no rates.",
                "the convergence figure could not be drawn: no rounds"]
    stop, out, verdict, log = run_check(tmp_path, converged_part("R1:0a55be94"),
                                        converged_part("R2:6544a5c3"),
                                        report_warnings=warnings)
    assert stop is None
    assert verdict.read_text() == "R1:0a55be94 converged; R2:6544a5c3 converged\n"
    err = capsys.readouterr().err
    assert ("NOTE: the error-model report has 2 warning(s), the first: R2: the fit figure "
            "of the R2 reads could not be drawn: no rates. All of them are in "
            f"{out}/stats/dada2/error_model_report.json") in err
    assert "warning(s)" in log.read_text()


def test_a_report_without_warnings_adds_no_note(tmp_path, capsys):
    run_check(tmp_path, converged_part("R1:0a55be94"), converged_part("R2:6544a5c3"))
    assert "error-model report" not in capsys.readouterr().err
