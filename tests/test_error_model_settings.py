"""The start-up checks of DADA2's error learning (amplicon.dada2.learn_errors) and the
files a run asks for because of it (workflow/rules/shared/00_common.smk).

What is tested, and why it matters:
  - max_consist, the cap on learning rounds. The script used to read it with
    as.integer(), so a wrong value went wrong without a word: 0 made DADA2 hand back its
    rough starting estimate as the model, 2.5 became 2, and Inf or text failed at the
    end of the first round. Now a whole number of at least 1 is required, and 999 (the
    value every template carried up to v2.6.1) gets a note, because with it a model that
    keeps changing runs for days before the check can stop the run.
  - max_flat_rounds, the N of the early stop: a whole number of at least 0 (0 = off),
    with a note when it can never act (N at or above max_consist).
  - accept_unconverged, the codes of models accepted although they did not converge.
    true or false must stop the run (an on/off switch would accept any model, however
    often the data change), a single code is wrapped into a list, codes are normalised
    (R1/R2 in capitals, the md5 characters in lower case) and sorted, and anything that
    is not a code stops the run with the item named.
  - error_model, the error function the models are learned with: auto, a function
    name or custom (the user's own function, experimental), in any capitals and with
    spaces around it; anything else (an instrument name such as novaseq, empty text,
    null, true, a number) stops the run with the allowed words named.
  - The learn_errors block itself (_learn_errors_block): missing blocks are created
    inside the config, and a dada2: or learn_errors: block left empty in the run config
    (its lines all commented out, which YAML reads as null and which replaces the
    template's block) stops the run with a message that names the block, instead of a
    bare Python error. The stop text of check_error_model asks users to write a
    learn_errors: block holding only accept_unconverged, so commenting that line out is
    a natural edit; .test/config/expected_stop/16S_learn_errors_block_empty.yaml checks
    the same end to end in a dry run.
  - custom_error_function (EXPERIMENTAL), the R file and the name of the user's own
    error function: every branch of its start-up check. Under error_model: custom both
    must be set; under any other setting an empty block means no function, a complete
    one is a comparison for the report, and half a block stops the run. The block must
    be a block, its values text, the file an existing, readable file (a relative path
    is read from the folder snakemake is started in, like every path of the config),
    and the name one R can give a function. A wrong value would otherwise surface as
    Snakemake's bare "missing input files" (a file) or, hours later, as an R error
    inside learn_errors (a name).
  - The file of the user's own function is an input of learn_errors only under
    error_model: custom (editing a function named for comparison must not make the
    models be learned again), and an input of the report whenever one is named.
  - The conda environment of those two steps: dada2_custom.yaml (dada2.yaml plus the R
    packages mgcv and scam) under the same two conditions, dada2.yaml otherwise, and
    dada2.yaml for every other DADA2 step, so that a run without a function of the
    user's own never builds dada2_custom.yaml (conda-forge has no r-scam for Apple
    Silicon Macs or Linux on ARM). dada2_custom.yaml stays dada2.yaml plus the two
    packages, and its pin holds every build of dada2's pin plus mgcv, nlme and scam.
  - The three lists of selectable error functions agree: the registry in
    workflow/scripts/amplicon/60_error_functions.R (read between its marker lines),
    ERROR_FUNCTION_NAMES in 00_common.smk (what the start-up check accepts), and
    SUITED_CLASS and ON_OTHER_CLASS in 60c_choose_error_model.py (what the choice step
    knows about each function). A function added to one list and not the others would
    be accepted at start-up and fail hours later, or never be offered.
  - The vendored copy of DADA2's makeBinnedQualErrfun is byte for byte the original:
    the md5 of the lines between its markers equals the md5 in its header. An editor
    that strips trailing spaces, or a hand edit, fails here.
  - The resources of learn_errors and error_model_report: learn_errors does the heavy
    half of what dada_seqtab did alone up to v2.6.1 and falls back to dada_seqtab's
    threads and memory; error_model_report draws figures in one thread and falls back
    to 1 thread and to the memory of dada_quality_plots, so that it neither holds eight
    cores of a workstation nor waits for an 8-core, 16-GB cluster slot. With the
    template that is 8 threads and 16000 MB for learning, 1 and 4000 for the report, as
    the configuration page says; a run config that names either rule still wins.
  - The quality values dada_filter counts have a file of their own per sample, which
    choose_error_model reads (filter_stats.json keeps its v2.6.1 layout), and which is
    not a target of the run.
  - _amplicon_targets() holds the choice of error function, the two models, the
    record, every file of the report and the verdict of the check. The verdict as a
    target is what makes sure an acceptance takes effect (see the comment in
    _amplicon_targets). The quality-values table is not a target: no step reads it,
    and as a target its deletion would make the models be learned again.
  - accept_unconverged reaches the check through the config only: 00_common.smk writes
    the normalised list back into the config (one line), and rule check_error_model
    does not take it as a param, so adding or removing a code makes no step run again.
  - The report's 22 files (_error_model_report_files): all under stats/dada2/, the six
    files of the MultiQC section in their own folder with names MultiQC recognises,
    and DADA2's own plots under the names they have had since v1.0.0.
  - The MultiQC section (workflow/multiqc_config.yaml): every figure of the section has
    its own search pattern (without it MultiQC 1.35 makes each figure a module of its
    own), every figure belongs to the module "DADA2 error model", the order lists the
    module and its six sections, and the two tables of 60e_error_model_report.R carry
    the section ids that order names.
  - The statement of what the comparison of the error functions can and cannot show is
    the same text in 60e_error_model_report.R (the report JSON), in the MultiQC config
    (the module description) and on the documentation page of the error model
    (docs/amplicon/error-model.md), so that the three never promise different things.
  - Falco on the filtered reads, the reads the error models are learned from:
    FALCO_STAGES in 00_common.smk ends with the stage "filtered", with and without PhiX
    removal, so a run asks for a Falco report of every sample's filtered R1 and R2
    reads; and the renames of the MultiQC config, applied in order as MultiQC applies
    them, give each stage its own label ("QC | quality-filtered reads | S1_R2"). The
    stop messages of choose_error_model and learn_errors send the user to these
    reports. The dry runs do not compare job counts, so these tests are what notices
    a stage left out.

The functions are loaded from 00_common.smk with load_functions() (helpers.py), which
runs their definitions only; a stop is a SystemExit whose text the tests read.
"""
import hashlib
import os
import re
from pathlib import Path

import pytest

import yaml

from helpers import COMMON_SMK, REPO_ROOT, SCRIPTS, TEST_DIR, WORKFLOW, load_functions

from test_orientation_table import amplicon_targets

checks = load_functions(COMMON_SMK, [
    "_whole_number", "_max_consist_setting", "_max_consist_note",
    "_max_flat_rounds_setting", "_max_flat_rounds_note", "_accept_unconverged_setting",
    "ERROR_FUNCTION_NAMES", "CUSTOM_ERROR_MODEL", "_error_model_setting",
    "R_RESERVED_WORDS", "_is_r_name", "_custom_error_function_text",
    "_custom_error_function_setting", "_learn_errors_block",
    "_error_model_report_files", "mem_mb_for_or",
])

REGISTRY_FILE = SCRIPTS / "amplicon" / "60_error_functions.R"
VENDORED_FILE = SCRIPTS / "amplicon" / "vendor" / "dada2_makeBinnedQualErrfun.R"
CHOOSE_SCRIPT = SCRIPTS / "amplicon" / "60c_choose_error_model.py"
REPORT_SCRIPT = SCRIPTS / "amplicon" / "60e_error_model_report.R"
MULTIQC_CONFIG = WORKFLOW / "multiqc_config.yaml"


def stop_text(function, value):
    """Call a check that must stop the run and return its message."""
    with pytest.raises(SystemExit) as stop:
        function(value)
    return str(stop.value)


# ── max_consist ───────────────────────────────────────────────────────────────

@pytest.mark.parametrize("value", [200, 200.0, "200", " 200 ", 1])
def test_max_consist_accepts_whole_numbers_of_at_least_one(value):
    expected = 1 if value == 1 else 200
    assert checks._max_consist_setting(value) == expected


@pytest.mark.parametrize("value", [0, -1, 2.5, "200.5", True, False, None, "Inf",
                                   float("inf"), float("nan"), [], "abc"])
def test_max_consist_refuses_everything_else(value):
    message = stop_text(checks._max_consist_setting, value)
    assert "max_consist must be a whole number of at least 1" in message
    assert repr(value) in message


def test_max_consist_note_only_for_999():
    note = checks._max_consist_note(999)
    assert "999" in note and "200" in note
    assert checks._max_consist_note(200) is None
    assert checks._max_consist_note(1000) is None


# ── max_flat_rounds ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("value, expected", [(0, 0), (20, 20), ("20", 20), (20.0, 20)])
def test_max_flat_rounds_accepts_whole_numbers_of_at_least_zero(value, expected):
    assert checks._max_flat_rounds_setting(value) == expected


@pytest.mark.parametrize("value", [-1, 2.5, True, None, "twenty", [20]])
def test_max_flat_rounds_refuses_everything_else(value):
    message = stop_text(checks._max_flat_rounds_setting, value)
    assert "max_flat_rounds must be a whole number >= 0" in message


def test_max_flat_rounds_note_when_the_rule_can_never_act():
    note = checks._max_flat_rounds_note(200, 200)
    assert "can never act" in note and "(200)" in note
    assert checks._max_flat_rounds_note(250, 200) is not None
    # One round below the cap, the rule can still act in the last round.
    assert checks._max_flat_rounds_note(199, 200) is None
    assert checks._max_flat_rounds_note(20, 200) is None
    # 0 means off: no note.
    assert checks._max_flat_rounds_note(0, 200) is None


# ── accept_unconverged ────────────────────────────────────────────────────────

def test_accept_unconverged_empty_values_give_an_empty_list():
    assert checks._accept_unconverged_setting(None) == []
    assert checks._accept_unconverged_setting([]) == []


def test_accept_unconverged_wraps_a_single_code():
    assert checks._accept_unconverged_setting("R2:0a55be94") == ["R2:0a55be94"]
    assert checks._accept_unconverged_setting("  r2:0A55BE94 ") == ["R2:0a55be94"]


def test_accept_unconverged_normalises_sorts_and_drops_duplicates():
    codes = ["r2:0A55BE94", "R1:26dab1e9", "R2:0a55be94"]
    assert checks._accept_unconverged_setting(codes) == ["R1:26dab1e9", "R2:0a55be94"]


@pytest.mark.parametrize("value", [True, False])
def test_accept_unconverged_refuses_true_and_false(value):
    message = stop_text(checks._accept_unconverged_setting, value)
    assert "accept_unconverged takes the codes" in message
    assert "not true or false" in message


@pytest.mark.parametrize("value, named", [
    (["R3:12345678"], "'R3:12345678'"),
    (["R2:0a55be9"], "'R2:0a55be9'"),       # 7 characters
    (["R2-0a55be94"], "'R2-0a55be94'"),     # a dash instead of the colon
    (["R2:0a55be9g"], "'R2:0a55be9g'"),     # g is not a hex character
    ([123], "123"),
    ({"R2": "x"}, "{'R2': 'x'}"),           # what "R2: 0a55be94" (with a space) becomes in YAML
    (12345678, "12345678"),
])
def test_accept_unconverged_refuses_what_is_not_a_code(value, named):
    message = stop_text(checks._accept_unconverged_setting, value)
    assert "accept_unconverged must be a list of model codes" in message
    assert f"got {named}" in message


# ── error_model ───────────────────────────────────────────────────────────────

@pytest.mark.parametrize("value, expected", [
    ("auto", "auto"), ("AUTO", "auto"), (" Loess ", "loess"), ("binned", "binned"),
    ("custom", "custom"), (" Custom ", "custom"),
])
def test_error_model_accepts_auto_the_function_names_and_custom(value, expected):
    assert checks._error_model_setting(value) == expected


@pytest.mark.parametrize("value", ["novaseq", "", None, True, 3, ["loess"]])
def test_error_model_refuses_everything_else(value):
    message = stop_text(checks._error_model_setting, value)
    assert "error_model must be one of auto, loess, binned, custom" in message
    assert f"got {value!r}" in message


# ── The learn_errors block itself ─────────────────────────────────────────────────

def test_a_block_whose_lines_are_all_commented_out_reads_as_null():
    # The premise of the stop below: YAML ignores the comment, so the key is left with
    # no value, which is null (None).
    text = "dada2:\n  learn_errors:\n    # accept_unconverged: [\"R2:0a55be94\"]\n"
    assert yaml.safe_load(text) == {"dada2": {"learn_errors": None}}


def test_missing_blocks_are_created_and_the_block_of_the_config_is_returned():
    # No dada2 block at all: both are created, empty, inside the config itself, so that
    # what the amplicon block writes into the returned dict (the normalised
    # accept_unconverged list) reaches every job that reads the config.
    amp_cfg = {}
    block = checks._learn_errors_block(amp_cfg)
    assert block == {} and amp_cfg == {"dada2": {"learn_errors": {}}}
    block["accept_unconverged"] = ["R2:0a55be94"]
    assert amp_cfg["dada2"]["learn_errors"]["accept_unconverged"] == ["R2:0a55be94"]
    # A block that is there is returned as it is, the same object.
    amp_cfg = {"dada2": {"max_ee": 2, "learn_errors": {"nbases": "1e8"}}}
    assert checks._learn_errors_block(amp_cfg) is amp_cfg["dada2"]["learn_errors"]


@pytest.mark.parametrize("amp_cfg, named, lost", [
    ({"dada2": None}, "amplicon.dada2 is empty in your run config", "max_ee, trunc_q and nbases"),
    ({"dada2": {"max_ee": 2, "learn_errors": None}},
     "amplicon.dada2.learn_errors is empty in your run config", "such as nbases, are lost"),
])
def test_a_block_left_empty_stops_with_its_name(amp_cfg, named, lost):
    message = stop_text(checks._learn_errors_block, amp_cfg)
    assert message.startswith(f"[MetaFlux] {named}")
    assert lost in message
    assert "or keep at least one setting under it" in message


# ── custom_error_function (experimental) ──────────────────────────────────────

BLOCK_KEY = "amplicon.dada2.learn_errors.custom_error_function"


def custom_stop(block, error_model):
    """Call the check of custom_error_function where it must stop; return the message."""
    with pytest.raises(SystemExit) as stop:
        checks._custom_error_function_setting(block, error_model)
    return str(stop.value)


@pytest.fixture
def r_file(tmp_path):
    """A small R file defining an error function, as a user would write one."""
    path = tmp_path / "my_errfun.R"
    path.write_text("my_errfun <- function(trans) dada2::loessErrfun(trans)\n")
    return path


@pytest.mark.parametrize("error_model", ["auto", "loess", "binned"])
@pytest.mark.parametrize("block", [None, "", {}, {"file": "", "name": ""},
                                   {"file": None, "name": None}, {"file": "  ", "name": ""}])
def test_custom_block_left_empty_means_no_function(block, error_model):
    # The template ships the block with both values empty: no function of the user's.
    # The block itself written as "" (custom_error_function: "") means the same.
    assert checks._custom_error_function_setting(block, error_model) is None


@pytest.mark.parametrize("block", [None, "", {}, {"file": "", "name": "my_errfun"}])
def test_custom_needs_a_file(block):
    message = custom_stop(block, "custom")
    assert message.startswith("[MetaFlux] amplicon.dada2.learn_errors.error_model is custom, so "
                              f"{BLOCK_KEY}.file must give the R file")
    assert "it is empty" in message
    assert "set error_model back to auto, loess or binned" in message


def test_custom_needs_a_name(r_file):
    message = custom_stop({"file": str(r_file), "name": ""}, "custom")
    assert f"{BLOCK_KEY}.name must give the name of your error function in {r_file}" in message
    assert "it is empty" in message


@pytest.mark.parametrize("given, missing", [("file", "name"), ("name", "file")])
def test_half_a_block_stops_under_the_other_settings(r_file, given, missing):
    block = {"file": "", "name": ""}
    block[given] = str(r_file) if given == "file" else "my_errfun"
    message = custom_stop(block, "loess")
    assert f"{BLOCK_KEY} sets a {given} but no {missing}" in message
    assert "Set both" in message and "or leave both empty" in message


@pytest.mark.parametrize("block", ["my_errfun.R", ["my_errfun.R", "my_errfun"], 3])
def test_the_block_must_be_a_block_of_keys(block):
    message = custom_stop(block, "custom")
    assert f"{BLOCK_KEY} must be a block with two keys" in message
    assert f"got {block!r}" in message


def test_file_and_name_must_be_text(r_file):
    message = custom_stop({"file": 3, "name": "my_errfun"}, "custom")
    assert message == f"[MetaFlux] {BLOCK_KEY}.file must be text, got 3"
    message = custom_stop({"file": str(r_file), "name": True}, "custom")
    assert message == f"[MetaFlux] {BLOCK_KEY}.name must be text, got True"


@pytest.mark.parametrize("error_model", ["custom", "loess"])
def test_a_missing_file_stops_with_the_resolved_path(tmp_path, monkeypatch, error_model):
    # A relative path is read from the folder snakemake is started in; the message names
    # the path as resolved, as written, and that folder.
    monkeypatch.chdir(tmp_path)
    message = custom_stop({"file": "funcs/missing.R", "name": "my_errfun"}, error_model)
    assert f"{BLOCK_KEY}.file: {tmp_path.resolve() / 'funcs' / 'missing.R'} is not an existing file" in message
    assert "written in the config as 'funcs/missing.R'" in message
    assert f"here {Path.cwd()}" in message


def test_a_folder_is_not_a_file(tmp_path):
    message = custom_stop({"file": str(tmp_path), "name": "my_errfun"}, "custom")
    assert "is not an existing file" in message


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0,
                    reason="root can read a file without read permission")
def test_an_unreadable_file_stops(r_file):
    r_file.chmod(0)
    try:
        message = custom_stop({"file": str(r_file), "name": "my_errfun"}, "custom")
    finally:
        r_file.chmod(0o644)
    assert f"{BLOCK_KEY}.file: {r_file.resolve()} cannot be read" in message


@pytest.mark.parametrize("name", [
    "1errfun",         # starts with a digit
    ".2x",             # a dot followed by a digit: R reads the number 0.2
    "_errfun",         # starts with an underscore
    "my-errfun",       # a dash, as in a file name
    "my errfun",       # a space inside
    "errfun()",        # brackets
    "naïve_errfun",    # a letter outside A-Z, a-z
    "function", "TRUE", "NA_real_", "in", "...", "..1",   # words R reserves
])
def test_a_name_r_cannot_give_a_function_stops(r_file, name):
    message = custom_stop({"file": str(r_file), "name": name}, "custom")
    assert f"{BLOCK_KEY}.name must be the name of an R function" in message
    # The message names the letters _is_r_name accepts, so an accented name is not told
    # that letters are allowed.
    assert "the letters a-z and A-Z (no accented letters)" in message
    assert f"got {name!r}" in message


@pytest.mark.parametrize("name", ["my_errfun", "gam.errfun2", "GAM_errfun", ".hidden", "._x",
                                  ".", "..x", "x.y_z"])
def test_names_r_can_give_a_function_pass(r_file, name):
    setting = checks._custom_error_function_setting({"file": str(r_file), "name": name}, "custom")
    assert setting["name"] == name


def test_a_relative_path_is_read_from_the_start_folder(tmp_path, monkeypatch):
    (tmp_path / "funcs").mkdir()
    (tmp_path / "funcs" / "gam.R").write_text("gam_errfun <- function(trans) trans\n")
    monkeypatch.chdir(tmp_path)
    setting = checks._custom_error_function_setting(
        {"file": " funcs/gam.R ", "name": " gam_errfun "}, "custom")
    assert setting == {"file": str((tmp_path / "funcs" / "gam.R").resolve()),
                       "file_as_configured": "funcs/gam.R", "name": "gam_errfun"}


@pytest.mark.parametrize("error_model", ["auto", "loess", "binned", "custom"])
def test_a_complete_block_is_returned_under_every_setting(r_file, error_model):
    # Under custom the models are learned with it; under the others the report draws it
    # as a comparison.
    setting = checks._custom_error_function_setting(
        {"file": str(r_file), "name": "my_errfun"}, error_model)
    assert setting == {"file": str(r_file.resolve()), "file_as_configured": str(r_file),
                       "name": "my_errfun"}


def test_the_shipped_example_passes_the_check():
    # The loessErrfun wrapper of the dry-run config 16S_custom_error_function, which the
    # documentation shows as the example.
    example = TEST_DIR / "custom_error_functions" / "loess_wrapper.R"
    setting = checks._custom_error_function_setting(
        {"file": str(example), "name": "loess_wrapper"}, "custom")
    assert setting["file"] == str(example.resolve())
    assert re.search(r"^loess_wrapper <- function\(trans\) \{$", example.read_text(), flags=re.MULTILINE)


# ── The selectable error functions, and the vendored DADA2 function ───────────

def registry_names():
    """The names of the registry entries in 60_error_functions.R, read between its
    marker lines: one entry per line written as `  name = list(` (two spaces)."""
    text = REGISTRY_FILE.read_text()
    start = text.index("# ---- registry: start ----")
    end = text.index("# ---- registry: end ----")
    block = text[start:end]
    return tuple(re.findall(r"^\s{2}(\w+)\s*=\s*list\(", block, flags=re.MULTILINE))


def test_the_three_lists_of_error_functions_agree():
    names = registry_names()
    assert names == ("loess", "binned")
    assert names == checks.ERROR_FUNCTION_NAMES
    tables = load_functions(CHOOSE_SCRIPT, ["SUITED_CLASS", "ON_OTHER_CLASS"])
    assert tuple(tables.SUITED_CLASS) == names
    assert tuple(tables.ON_OTHER_CLASS) == names
    assert set(tables.ON_OTHER_CLASS.values()) <= {"warn", "stop"}


def test_the_vendored_block_is_the_original():
    lines = VENDORED_FILE.read_text().split("\n")
    start = lines.index("# ---- vendored block: start ----")
    end = lines.index("# ---- vendored block: end ----")
    block = "".join(line + "\n" for line in lines[start + 1:end])
    header_text = "\n".join(lines[:start])
    recorded = re.search(r"md5 of the vendored block: ([0-9a-f]{32})", header_text).group(1)
    assert hashlib.md5(block.encode()).hexdigest() == recorded
    assert block.startswith("makeBinnedQualErrfun <- function(binnedQ) {\n")
    assert "Changes: none" in header_text


# ── Resources of learn_errors and error_model_report ──────────────────────────

def test_mem_mb_for_or_uses_the_named_value_or_the_fallback():
    settings = checks.mem_mb_for_or.__globals__
    settings["config"] = {"resources": {"mem_mb_default": 2000,
                                        "mem_mb": {"dada_seqtab": 16000, "learn_errors": 32000,
                                                   "dada_quality_plots": 4000}}}
    assert checks.mem_mb_for_or("learn_errors", 16000) == 32000
    # Not named: the fallback chosen by the rule, not mem_mb_default.
    assert checks.mem_mb_for_or("error_model_report", 4000) == 4000


def test_learning_falls_back_to_dada_seqtab_and_the_report_to_one_thread():
    # The lines of the two rules in 60_dada2.smk.
    learn = rule_text("learn_errors")
    assert 'threads: lambda wc: threads_for_or("learn_errors", threads_for("dada_seqtab"))' in learn
    assert 'mem_mb = lambda wc: mem_mb_for_or("learn_errors", mem_mb_for("dada_seqtab")),' in learn
    report = rule_text("error_model_report")
    assert 'threads: lambda wc: threads_for_or("error_model_report", 1)' in report
    assert ('mem_mb = lambda wc: mem_mb_for_or("error_model_report", mem_mb_for("dada_quality_plots")),'
            in report)
    # What they give with the template config/config.yaml, which names neither rule.
    lookups = load_functions(COMMON_SMK, ["threads_for", "threads_for_or", "mem_mb_for", "mem_mb_for_or"])
    with open(REPO_ROOT / "config" / "config.yaml") as handle:
        lookups.threads_for.__globals__["config"] = yaml.safe_load(handle)
    assert lookups.threads_for_or("learn_errors", lookups.threads_for("dada_seqtab")) == 8
    assert lookups.mem_mb_for_or("learn_errors", lookups.mem_mb_for("dada_seqtab")) == 16000
    assert lookups.threads_for_or("error_model_report", 1) == 1
    assert lookups.mem_mb_for_or("error_model_report", lookups.mem_mb_for("dada_quality_plots")) == 4000


# ── The files a run asks for ──────────────────────────────────────────────────

# The six files of the MultiQC section, by their key in _error_model_report_files.
MULTIQC_KEYS = ("mqc_summary", "mqc_fit_table", "mqc_fit_r1", "mqc_fit_r2",
                "mqc_convergence", "mqc_quality")
MULTIQC_IMAGE_KEYS = ("mqc_fit_r1", "mqc_fit_r2", "mqc_convergence", "mqc_quality")


def test_report_files_are_built_from_the_output_folder():
    files = checks._error_model_report_files(Path("/somewhere"))
    assert len(files) == 22
    assert len(set(files.values())) == 22                   # no two keys share a file
    assert all(str(path).startswith("/somewhere/stats/dada2/") for path in files.values())
    # The plots keep the names they had when dada_seqtab drew them.
    assert files["dada2_plot_r2"].name == "filtered_read_R2_error_plot.png"
    assert files["report_json"].name == "error_model_report.json"
    # The MultiQC section has a folder of its own, and MultiQC recognises its files by
    # the ending _mqc.tsv or _mqc.png.
    for key in MULTIQC_KEYS:
        assert files[key].parent == Path("/somewhere/stats/dada2/multiqc")
        assert files[key].name.endswith(("_mqc.tsv", "_mqc.png"))
    # Every other file sits directly in stats/dada2/, where MultiQC does not look.
    for key in set(files) - set(MULTIQC_KEYS):
        assert files[key].parent == Path("/somewhere/stats/dada2")


def test_targets_hold_the_models_the_record_the_report_and_the_verdict():
    targets = amplicon_targets("fixed")
    report = checks._error_model_report_files(Path("/out"))
    expected = {str(path) for path in report.values()} | {
        "/out/stats/dada2/error_model_choice.json",
        "/out/05.dada2/error_model_R1.rds",
        "/out/05.dada2/error_model_R2.rds",
        "/out/stats/dada2/error_model.json",
        "/out/stats/dada2/error_model_check.txt",
    }
    assert expected <= targets
    # The quality-values table is written by choose_error_model for people and read by
    # no step. As a target, deleting it would make choose_error_model run again and,
    # through the choice it rewrites, the learning of the error models and every step
    # after it; so the run does not ask for it.
    assert "/out/stats/dada2/quality_values.tsv" not in targets


def test_the_quality_values_have_a_file_of_their_own():
    # dada_filter writes the quality scan to stats/dada2/{sample}.quality_values.json, and
    # choose_error_model reads those files, not filter_stats.json, which keeps the layout
    # of v2.6.1 and is read by dada_seqtab only. In an output folder made before the file
    # existed, the missing input of choose_error_model (a step every run needs) makes
    # dada_filter run again under any --rerun-triggers; a field inside filter_stats.json
    # left that to the script's date, which only the code trigger checks.
    outputs = rule_section("dada_filter", "output")
    assert 'stats          = OUT / "stats" / "dada2" / "{sample}.filter_stats.json",' in outputs
    assert 'quality_values = OUT / "stats" / "dada2" / "{sample}.quality_values.json",' in outputs
    inputs = rule_section("choose_error_model", "input")
    assert ('quality_values = expand(OUT / "stats" / "dada2" / "{sample}.quality_values.json", '
            'sample=SAMPLES),') in inputs
    assert "filter_stats" not in inputs
    # Not a target: a deleted copy must not filter its sample again, which would make the
    # models be learned again through the new filtered reads.
    targets = amplicon_targets("fixed")
    assert not any(path.endswith(".quality_values.json") for path in targets)


# ── accept_unconverged reaches the check through the config ───────────────────

DADA2_RULES = WORKFLOW / "rules" / "amplicon" / "60_dada2.smk"


def rule_text(rule_name):
    """The whole text of one rule in 60_dada2.smk, from `rule <name>:` to the next rule."""
    text = DADA2_RULES.read_text()
    rule_start = text.index(f"\nrule {rule_name}:\n")
    next_rule = text.find("\nrule ", rule_start + 1)
    return text[rule_start:] if next_rule == -1 else text[rule_start:next_rule]


def rule_section(rule_name, section):
    """The text of one section (input, output, params, ...) of a rule in 60_dada2.smk:
    the lines after `    <section>:` up to the next section of the same rule."""
    match = re.search(rf"^    {section}:\n(.*?)(?=^    \w+:|\Z)", rule_text(rule_name),
                      flags=re.MULTILINE | re.DOTALL)
    return match.group(1) if match else ""


def test_accept_unconverged_reaches_the_check_through_the_config_only():
    # The check (60f_check_error_model.py) reads the accepted codes from the config
    # inside its job. Two lines make that work: the amplicon block of 00_common.smk
    # writes the normalised list back into the config, and the rule does not take the
    # codes as a param, so that adding or removing a code makes no step run again
    # (with a param, removing a code after an accepted run would rerun the check, which
    # would then stop the run again).
    common = COMMON_SMK.read_text()
    write_backs = re.findall(r'^\s+_learn_cfg\["accept_unconverged"\] = ACCEPT_UNCONVERGED\s*$',
                             common, flags=re.MULTILINE)
    assert len(write_backs) == 1
    params = rule_section("check_error_model", "params")
    assert "orientation" in params          # the section was found
    assert "accept_unconverged" not in params.lower()


def test_the_custom_function_file_is_an_input_where_it_should_be():
    # learn_errors takes the user's file (and its name) only under error_model: custom:
    # editing a function named for comparison only must never make the models be learned
    # again. The report takes the file whenever one is named, so that an edit draws the
    # function again. (The dry run of .test/config/16S_custom_error_function.yaml plans
    # both inputs; the CI does not compare the plans, hence this test.)
    learn_inputs = rule_section("learn_errors", "input")
    assert ('custom_function_file = [CUSTOM_ERROR_FUNCTION["file"]] if ERROR_MODEL_SETTING == "custom" else []'
            in learn_inputs)
    learn_params = rule_section("learn_errors", "params")
    assert 'custom_function_name               = CUSTOM_ERROR_FUNCTION["name"] if ERROR_MODEL_SETTING == "custom" else ""' in learn_params
    # The path as written in the config goes into the record; under the same condition,
    # so that a function named for comparison only never makes the learning run again.
    assert ('custom_function_file_as_configured = CUSTOM_ERROR_FUNCTION["file_as_configured"] '
            'if ERROR_MODEL_SETTING == "custom" else ""') in learn_params
    report_inputs = rule_section("error_model_report", "input")
    assert 'custom_function_file = [CUSTOM_ERROR_FUNCTION["file"]] if CUSTOM_ERROR_FUNCTION else []' in report_inputs
    report_params = rule_section("error_model_report", "params")
    assert 'custom_function_name               = CUSTOM_ERROR_FUNCTION["name"] if CUSTOM_ERROR_FUNCTION else ""' in report_params
    assert ('custom_function_file_as_configured = CUSTOM_ERROR_FUNCTION["file_as_configured"] '
            'if CUSTOM_ERROR_FUNCTION else ""') in report_params


# ── The conda environment of the user's own error function ────────────────────

ENVS = WORKFLOW / "envs"


def test_the_two_steps_run_in_dada2_custom_only_when_they_call_the_users_function():
    # learn_errors calls the user's function only under error_model: custom, the report
    # whenever one is named (the conditions of the custom_function_file inputs above).
    # Exactly then each runs in dada2_custom.yaml (dada2.yaml plus the R packages mgcv
    # and scam), otherwise in dada2.yaml, so a run without such a function never builds
    # the second environment, which conda cannot build as shipped on Apple Silicon Macs
    # or Linux on ARM (no r-scam there). The dry runs plan both environments, but the
    # runner does not check which rule gets which, hence this test.
    learn_conda = rule_section("learn_errors", "conda")
    assert ('"../../envs/dada2_custom.yaml" if ERROR_MODEL_SETTING == "custom" '
            'else "../../envs/dada2.yaml"') in learn_conda
    report_conda = rule_section("error_model_report", "conda")
    assert '"../../envs/dada2_custom.yaml" if CUSTOM_ERROR_FUNCTION else "../../envs/dada2.yaml"' in report_conda
    # Every other DADA2 step of the file stays in dada2.yaml, whatever the settings.
    for rule_name in ("dada_quality_plots", "dada_filter", "dada_seqtab"):
        assert rule_section(rule_name, "conda").strip() == '"../../envs/dada2.yaml"'


def pin_package_lines(pin_name):
    """The package lines of a pin file in workflow/envs/ (one download address per package
    build). The header lines are left out: they name the conda version that wrote the
    file, which may differ between the two files."""
    lines = (ENVS / pin_name).read_text().splitlines()
    return [line for line in lines if line.startswith("https://")]


def package_name(pin_line):
    """The package name of a pin line: .../linux-64/r-mgcv-1.9_4-r44h0e4624f_0.conda#<md5>
    gives r-mgcv (the file name without its version and build)."""
    file_name = pin_line.split("#")[0].split("/")[-1]
    return file_name.rsplit("-", 2)[0]


def test_dada2_custom_is_the_dada2_environment_plus_mgcv_and_scam():
    # dada2_custom.yaml must stay dada2.yaml plus r-mgcv and r-scam: the report refits
    # MetaFlux's own functions in it, and a model learned in it must rest on the same
    # software as one learned in dada2.yaml. On Linux its pin must hold every build of
    # dada2's pin, which is why that pin is rebuilt from dada2's and never solved
    # (CONTRIBUTING.md). The lines it adds are the two packages and nlme, which mgcv
    # needs; a refresh that brings in anything else fails here on purpose, so that the
    # new package is looked at before it ships.
    dada2_env = yaml.safe_load((ENVS / "dada2.yaml").read_text())
    custom_env = yaml.safe_load((ENVS / "dada2_custom.yaml").read_text())
    assert custom_env["channels"] == dada2_env["channels"]
    assert custom_env["dependencies"] == dada2_env["dependencies"] + ["r-mgcv", "r-scam"]

    dada2_lines = pin_package_lines("dada2.linux-64.pin.txt")
    custom_lines = pin_package_lines("dada2_custom.linux-64.pin.txt")
    missing_builds = [line for line in dada2_lines if line not in custom_lines]
    assert missing_builds == []
    added_lines = [line for line in custom_lines if line not in dada2_lines]
    assert sorted(package_name(line) for line in added_lines) == ["r-mgcv", "r-nlme", "r-scam"]



# ── The MultiQC section of the report ─────────────────────────────────────────

def multiqc_config():
    """workflow/multiqc_config.yaml as MultiQC reads it (YAML anchors resolved)."""
    with open(MULTIQC_CONFIG) as handle:
        return yaml.safe_load(handle)


def normalised(text):
    """The text with every run of spaces and line ends made one space."""
    return " ".join(text.split())


def test_every_figure_of_the_section_has_its_search_pattern_and_module():
    config = multiqc_config()
    files = checks._error_model_report_files(Path("/out"))
    for key in MULTIQC_IMAGE_KEYS:
        file_name = files[key].name
        section_id = file_name[:-len("_mqc.png")]          # MultiQC's id of the image
        assert config["sp"][section_id] == {"fn": file_name}
        entry = config["custom_data"][section_id]
        assert entry["parent_id"] == "metaflux_error_model"
        assert entry["parent_name"] == "DADA2 error model"
        assert entry["section_name"] and entry["description"]


def test_the_section_order_lists_the_module_and_its_six_sections():
    order = multiqc_config()["report_section_order"]
    files = checks._error_model_report_files(Path("/out"))
    image_sections = [files[key].name[:-len("_mqc.png")] for key in MULTIQC_IMAGE_KEYS]
    table_sections = ["metaflux_error_model_summary", "metaflux_error_model_goodness_of_fit"]
    for section_id in ["metaflux_error_model"] + image_sections + table_sections:
        assert section_id in order
    # Summary first, goodness of fit last; the module itself between the tool modules
    # and Software Versions.
    inside = sorted(image_sections + table_sections, key=lambda section: order[section]["order"])
    assert inside == ["metaflux_error_model_summary", "error_model_quality_values",
                      "error_model_convergence", "error_model_fit_R1", "error_model_fit_R2",
                      "metaflux_error_model_goodness_of_fit"]
    # The two tables written by the report script carry the ids the order names.
    script = REPORT_SCRIPT.read_text()
    for section_id in table_sections:
        assert f'"{section_id}"' in script


def report_statement():
    """Statement S as 60e_error_model_report.R holds it: the text between the quotes of
    REPORT_STATEMENT_TEXT, with every run of spaces and line ends made one space (the
    script does the same before using it)."""
    script = REPORT_SCRIPT.read_text()
    match = re.search(r'^REPORT_STATEMENT_TEXT <- "([^"]*)"', script, flags=re.MULTILINE)
    return normalised(match.group(1))


def test_the_statement_is_the_same_in_the_report_and_the_multiqc_config():
    statement = report_statement()
    assert statement.startswith("The comparison refits only the smoothing step:")
    assert statement.endswith("nor tell faulty reads from good ones.")
    for entry in multiqc_config()["custom_data"].values():
        if entry.get("parent_id") == "metaflux_error_model":
            assert normalised(entry["parent_description"]) == statement


def test_the_documentation_quotes_the_statement_word_for_word():
    # The page on the error model quotes the statement in its section on the report
    # (inside a note box, so indented and wrapped: line breaks and indentation do not
    # count). A change to the statement in the report script must reach the page too.
    page = REPO_ROOT / "docs" / "amplicon" / "error-model.md"
    assert report_statement() in normalised(page.read_text())


# ── Falco on the filtered reads ───────────────────────────────────────────────

def falco_stages_in_common_smk(remove_phix):
    """FALCO_STAGES as 00_common.smk sets it, with PhiX removal on or off.

    The line sits inside the amplicon block of 00_common.smk, which load_functions()
    cannot reach (it runs top-level definitions only), so the test finds the one line
    "FALCO_STAGES = ..." and evaluates its right-hand side: a plain list expression
    that reads nothing but REMOVE_PHIX."""
    text = COMMON_SMK.read_text()
    found = re.findall(r"^\s+FALCO_STAGES\s*=\s*(.+)$", text, flags=re.MULTILINE)
    assert len(found) == 1, "00_common.smk should set FALCO_STAGES on exactly one line"
    return eval(found[0], {"REMOVE_PHIX": remove_phix})


def test_falco_profiles_the_filtered_reads_with_and_without_phix_removal():
    # In pipeline order, the filtered reads last. These are also the two lists the
    # target tests use (amplicon_targets() in test_orientation_table.py).
    assert falco_stages_in_common_smk(True) == ["raw", "nophix", "stripped", "filtered"]
    assert falco_stages_in_common_smk(False) == ["raw", "stripped", "filtered"]


def test_targets_hold_a_falco_report_of_the_filtered_reads_of_each_sample_and_direction():
    expected = {"/out/stats/falco/S1_R1_filtered", "/out/stats/falco/S1_R2_filtered",
                "/out/stats/falco/S2_R1_filtered", "/out/stats/falco/S2_R2_filtered"}
    for remove_phix in (True, False):
        targets = amplicon_targets("fixed", remove_phix=remove_phix)
        filtered = {path for path in targets
                    if path.startswith("/out/stats/falco/") and path.endswith("_filtered")}
        assert filtered == expected


def multiqc_sample_name(name):
    """The name MultiQC shows for a sample of the given name, after the renames of
    multiqc_config.yaml: with sample_names_replace_regex on, MultiQC 1.35 applies every
    pattern of sample_names_replace in turn, in the order of the file, with Python's
    re.sub (clean_s_name in MultiQC's base_module.py)."""
    config = multiqc_config()
    assert config["sample_names_replace_regex"] is True
    for pattern, replacement in config["sample_names_replace"].items():
        name = re.sub(pattern, replacement, name)
    return name


def test_multiqc_gives_each_falco_stage_its_own_label():
    # Rule falco (20_qc.smk) runs Falco on a link named <sample>_R<r>_<stage>.fastq.gz;
    # Falco records that file name, and MultiQC takes it without .fastq.gz as the
    # sample's name. No later rename may change a label once given.
    assert multiqc_sample_name("S1_R1_raw") == "QC | raw reads | S1_R1"
    assert multiqc_sample_name("S1_R1_nophix") == "QC | PhiX-filtered reads | S1_R1"
    assert multiqc_sample_name("S1_R1_stripped") == "QC | primer-trimmed reads | S1_R1"
    assert multiqc_sample_name("S1_R2_filtered") == "QC | quality-filtered reads | S1_R2"
