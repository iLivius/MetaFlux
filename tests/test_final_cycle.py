"""The removal of the extra last sequencing cycle in primer trimming
(amplicon.cutadapt.trim_final_cycle): rule trim_primers in
workflow/rules/amplicon/30_preprocess.smk, its program
workflow/scripts/amplicon/30a_final_cycle.py, and the start-up checks
_trim_final_cycle_setting and _read_direction_truncation in
workflow/rules/shared/00_common.smk.

Why it matters: an Illumina run sequences one cycle more than its read length (2 x 151,
2 x 251, 2 x 301), and the calls of that extra cycle are poor. Where nothing truncates
the reads (ITS), a wrong last base stops DADA2 from merging the pairs; removing a good
base instead costs one base of overlap. Where DADA2 cuts the reads to a fixed length, the
cut decides, and removing the base there only moves the cut (18S test set), so auto
leaves those reads alone.

What is tested:
  - The decision of the program: with auto, the last base of an untruncated read
    direction is removed when the longest read is 76, 101, 151, 251 or 301 bp, and kept
    at any other length (300, 250, an uneven split such as 321, a file without reads);
    a truncated direction keeps it at every length; true always removes it, false never
    does, truncated or not.
  - Which read directions count as untruncated: both on ITS, whatever trunc_len says
    (60b_dada_filter.R sets truncLen to 0, 0); on the other markers only a direction
    with trunc_len.mode: manual and a cut of 0, read as pick_trunclen reads it.
  - How the program measures the longest read: in gzipped and in plain FASTQ files (the
    raw reads, which the rule reads when PhiX removal is off, may be either), in the
    first reads only (a longer read after them is not seen), in a file with fewer reads,
    and in an empty file.
  - The program as the rule runs it, with the Python interpreter: one word on standard
    output ("remove" or "keep"), which is all the rule's shell command reads, and one
    line on standard error, which the rule writes into its log; a mistyped setting or
    truncation word stops it with a message and no word.
  - The start-up check: true and false, and the words auto, true and false in any
    capitals, are accepted; anything else stops the run with the three values named. A
    value typed inline with --config arrives as text and is first read the way YAML
    reads a config file (_as_in_a_config_file), so "false" typed inline means false.
  - The rule: the cut goes into the 5' call of pass A and of pass B (which reads the
    input files again) and into neither 3' call, where it would remove a second base;
    the rule's own bash lines turn the word for R1 into --cut -1 and the word for R2
    into -U -1, never the other way round; the program is an input of the rule (that is
    what makes trimming run again in folders trimmed before it existed) and the checked
    setting and the truncation of each direction are params, passed to the program
    with its own read file; the log is opened by the program's first line and only added to
    after it, so Cutadapt's reports follow the decision instead of replacing it.
  - The template config/config.yaml sets auto, and 00_common.smk passes the setting
    through the check into TRIM_FINAL_CYCLE and the truncation of R1 and R2 into
    FINAL_CYCLE_TRUNCATION.

The program's functions and the start-up check are loaded with load_functions()
(helpers.py), which runs their definitions only; a stop is a SystemExit whose text the
tests read. The rule is read as text.
"""
import gzip
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from helpers import COMMON_SMK, REPO_ROOT, SCRIPTS, WORKFLOW, load_functions

PROGRAM = SCRIPTS / "amplicon" / "30a_final_cycle.py"
PREPROCESS_SMK = WORKFLOW / "rules" / "amplicon" / "30_preprocess.smk"

program = load_functions(PROGRAM, [
    "SCAN_READS", "FINAL_CYCLE_LENGTHS", "SETTINGS", "TRUNCATION_STATES",
    "open_fastq", "longest_read_length", "decide",
])
checks = load_functions(COMMON_SMK, ["_trim_final_cycle_setting", "_as_in_a_config_file",
                                     "_read_direction_truncation"])

STOP_FRAGMENT = "amplicon.cutadapt.trim_final_cycle must be auto, true or false"


def write_fastq(path, lengths, gzipped=True):
    """Write a FASTQ file with one read of each length in `lengths` (bases all A,
    qualities all I), gzipped unless gzipped is False. Returns the path as text, which
    is what the rule hands the program."""
    lines = []
    for number, length in enumerate(lengths):
        lines.append(f"@read{number}")
        lines.append("A" * length)
        lines.append("+")
        lines.append("I" * length)
    text = ""
    if lines:
        text = "\n".join(lines) + "\n"
    path = Path(path)
    if gzipped:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return str(path)


def run_program(setting, fastq_path, direction="R2", truncation="untruncated"):
    """Run the program as the rule does: python 30a_final_cycle.py <setting> <fastq>
    <direction> <truncation>. Returns the finished process (return code, stdout,
    stderr)."""
    return subprocess.run([sys.executable, "-I", str(PROGRAM), setting, fastq_path, direction,
                           truncation],
                          capture_output=True, text=True)


# ── The decision ──────────────────────────────────────────────────────────────

def test_the_lengths_auto_acts_on():
    # The read lengths of a run with one extra cycle, as the docs and the config list
    # them; a change here must be a change of those texts too.
    assert program.FINAL_CYCLE_LENGTHS == (76, 101, 151, 251, 301)
    assert program.SCAN_READS == 10_000


@pytest.mark.parametrize("longest", [76, 101, 151, 251, 301])
def test_auto_removes_the_last_base_at_the_length_of_a_run_with_an_extra_cycle(longest):
    assert program.decide("auto", longest, "untruncated") is True


@pytest.mark.parametrize("longest", [0, 75, 100, 150, 250, 281, 295, 300, 302, 321])
def test_auto_keeps_the_last_base_at_any_other_length(longest):
    # 300 and 250: runs set up without the extra cycle; 321 and 281: an uneven split of
    # a kit's cycles, which auto does not cover; 0: a file without reads.
    assert program.decide("auto", longest, "untruncated") is False


@pytest.mark.parametrize("longest", [0, 151, 250, 251, 300, 301])
def test_auto_keeps_the_last_base_of_truncated_reads_at_every_length(longest):
    # Where DADA2 cuts the reads to a fixed length, the cut decides (18S test set).
    assert program.decide("auto", longest, "truncated") is False


@pytest.mark.parametrize("truncation", ["truncated", "untruncated"])
@pytest.mark.parametrize("longest", [0, 250, 300, 301])
def test_true_and_false_depend_on_neither_length_nor_truncation(longest, truncation):
    assert program.decide("true", longest, truncation) is True
    assert program.decide("false", longest, truncation) is False


def test_the_truncation_words():
    assert program.TRUNCATION_STATES == ("truncated", "untruncated")


# ── Which read directions are untruncated (start-up) ──────────────────────────

@pytest.mark.parametrize("trunc_cfg", [
    {"mode": "auto", "manual_r1": None, "manual_r2": None},
    {"mode": "manual", "manual_r1": 260, "manual_r2": 240},
    {},
])
def test_its_is_untruncated_whatever_trunc_len_says(trunc_cfg):
    # 60b_dada_filter.R sets truncLen to 0, 0 for ITS, manual cuts included.
    for direction in ("r1", "r2"):
        assert checks._read_direction_truncation("ITS", trunc_cfg, direction) == "untruncated"


@pytest.mark.parametrize("marker", ["16S", "18S", "gyrB", "rpoB", "COI"])
def test_other_markers_are_truncated_in_auto_mode(marker):
    trunc_cfg = {"mode": "auto", "manual_r1": None, "manual_r2": None}
    for direction in ("r1", "r2"):
        assert checks._read_direction_truncation(marker, trunc_cfg, direction) == "truncated"


@pytest.mark.parametrize("manual_r1, manual_r2, expected_r1, expected_r2", [
    (0, 0, "untruncated", "untruncated"),
    (0, 240, "untruncated", "truncated"),
    ("0", 240, "untruncated", "truncated"),        # int(), as pick_trunclen reads it
    (260, 240, "truncated", "truncated"),
    (None, None, "truncated", "truncated"),        # pick_trunclen stops on missing cuts
    ("all", 0, "truncated", "untruncated"),        # not a number: pick_trunclen stops
])
def test_a_manual_cut_of_0_leaves_its_direction_untruncated(manual_r1, manual_r2,
                                                            expected_r1, expected_r2):
    trunc_cfg = {"mode": "manual", "manual_r1": manual_r1, "manual_r2": manual_r2}
    assert checks._read_direction_truncation("16S", trunc_cfg, "r1") == expected_r1
    assert checks._read_direction_truncation("16S", trunc_cfg, "r2") == expected_r2


def test_a_manual_cut_of_0_counts_only_in_manual_mode():
    # pick_trunclen reads manual_r1/manual_r2 only when mode is exactly "manual".
    trunc_cfg = {"mode": "auto", "manual_r1": 0, "manual_r2": 0}
    assert checks._read_direction_truncation("16S", trunc_cfg, "r1") == "truncated"


# ── The longest read ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("gzipped", [True, False])
def test_longest_read_of_a_gzipped_and_of_a_plain_file(tmp_path, gzipped):
    fastq = write_fastq(tmp_path / "S1_R1.fastq", [250, 301, 120], gzipped=gzipped)
    assert program.longest_read_length(fastq) == (301, 3)


def test_only_the_first_reads_are_read(tmp_path):
    # Three reads are read here; the longer fourth read is not seen.
    fastq = write_fastq(tmp_path / "S1_R1.fastq.gz", [300, 300, 300, 301])
    assert program.longest_read_length(fastq, max_reads=3) == (300, 3)


def test_a_file_without_reads(tmp_path):
    fastq = write_fastq(tmp_path / "S1_R1.fastq.gz", [])
    assert program.longest_read_length(fastq) == (0, 0)


# ── The program as the rule runs it ───────────────────────────────────────────

def test_the_program_prints_one_word_and_logs_why(tmp_path):
    fastq = write_fastq(tmp_path / "S1_R2.fastq.gz", [301, 290, 301])
    result = run_program("auto", fastq, "R2", "untruncated")
    assert result.returncode == 0, result.stderr
    assert result.stdout == "remove\n"
    assert result.stderr.startswith("[trim_primers] R2 final cycle: the last base of every "
                                    "read is removed")
    assert "the longest of the 3 reads of S1_R2.fastq.gz is 301 bp" in result.stderr
    assert "trim_final_cycle is auto, R2 is not truncated" in result.stderr


def test_the_program_keeps_the_base_of_truncated_reads_and_says_why(tmp_path):
    fastq = write_fastq(tmp_path / "S1_R1.fastq.gz", [301, 301])
    result = run_program("auto", fastq, "R1", "truncated")
    assert result.returncode == 0, result.stderr
    assert result.stdout == "keep\n"
    assert "R1 final cycle: the last base is kept" in result.stderr
    assert "DADA2 cuts R1 to a fixed length (truncLen)" in result.stderr
    result = run_program("true", fastq, "R1", "truncated")
    assert result.stdout == "remove\n"


def test_the_program_keeps_the_base_where_it_should(tmp_path):
    fastq_300 = write_fastq(tmp_path / "S1_R1.fastq.gz", [300, 300])
    fastq_301 = write_fastq(tmp_path / "S2_R1.fastq.gz", [301, 301])
    for setting, fastq in (("auto", fastq_300), ("false", fastq_301)):
        result = run_program(setting, fastq, "R1")
        assert result.returncode == 0, result.stderr
        assert result.stdout == "keep\n"
        assert "R1 final cycle: the last base is kept" in result.stderr
    result = run_program("true", fastq_300, "R1")
    assert result.stdout == "remove\n"
    assert "trim_final_cycle is true" in result.stderr


def test_a_mistyped_setting_stops_the_program(tmp_path):
    fastq = write_fastq(tmp_path / "S1_R1.fastq.gz", [301])
    result = run_program("always", fastq, "R1")
    assert result.returncode != 0
    assert result.stdout == ""
    assert "trim_final_cycle must be auto, true or false, got 'always'" in result.stderr


def test_a_mistyped_truncation_word_or_a_missing_word_stops_the_program(tmp_path):
    fastq = write_fastq(tmp_path / "S1_R1.fastq.gz", [301])
    result = run_program("auto", fastq, "R1", "untrimmed")
    assert result.returncode != 0
    assert result.stdout == ""
    assert "the truncation of R1 must be truncated or untruncated" in result.stderr
    # Three words, as the program took before the truncation was added.
    result = subprocess.run([sys.executable, "-I", str(PROGRAM), "auto", fastq, "R1"],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert result.stdout == ""
    assert "usage: 30a_final_cycle.py" in result.stderr


# ── The start-up check ────────────────────────────────────────────────────────

@pytest.mark.parametrize("value, expected", [
    (True, "true"), (False, "false"),
    ("auto", "auto"), (" Auto ", "auto"), ("TRUE", "true"), ("false", "false"),
])
def test_the_start_up_check_accepts(value, expected):
    assert checks._trim_final_cycle_setting(value) == expected


@pytest.mark.parametrize("value", [None, 1, 0, 2.5, "", "yes", "always", ["auto"]])
def test_the_start_up_check_refuses(value):
    with pytest.raises(SystemExit) as stop:
        checks._trim_final_cycle_setting(value)
    assert STOP_FRAGMENT in str(stop.value)


def test_a_value_typed_inline_means_what_it_means_in_a_config_file():
    # --config 'amplicon={cutadapt: {trim_final_cycle: false}}' hands over the text
    # "false"; 00_common.smk reads it back as YAML would, before any check.
    for typed, expected in (("false", "false"), ("true", "true"), ("auto", "auto")):
        value = checks._as_in_a_config_file(typed)
        assert checks._trim_final_cycle_setting(value) == expected


def test_the_template_and_the_wiring():
    template = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())
    assert template["amplicon"]["cutadapt"]["trim_final_cycle"] == "auto"
    common = COMMON_SMK.read_text()
    assert ('TRIM_FINAL_CYCLE = _trim_final_cycle_setting('
            '_cutadapt_cfg.get("trim_final_cycle", "auto"))') in common
    assert '"R1": _read_direction_truncation(AMPLICON_TYPE, _trunc_cfg, "r1"),' in common
    assert '"R2": _read_direction_truncation(AMPLICON_TYPE, _trunc_cfg, "r2"),' in common


# ── The rule ──────────────────────────────────────────────────────────────────

def trim_primers_text():
    """The whole text of rule trim_primers, from `rule trim_primers:` to the next rule
    (or the end of the file)."""
    text = PREPROCESS_SMK.read_text()
    start = text.index("\nrule trim_primers:\n")
    next_rule = text.find("\nrule ", start + 1)
    if next_rule == -1:
        return text[start:]
    return text[start:next_rule]


def cutadapt_calls(rule_text):
    """The text of every Cutadapt call of the rule, in order: from `cutadapt -j` up to
    the `2>&1` that ends its redirection to the log."""
    calls = []
    start = rule_text.find("cutadapt -j")
    while start != -1:
        end = rule_text.index("2>&1", start)
        calls.append(rule_text[start:end])
        start = rule_text.find("cutadapt -j", end)
    return calls


def test_the_cut_is_in_both_5_prime_calls_and_in_no_3_prime_call():
    calls = cutadapt_calls(trim_primers_text())
    five_prime_calls = [call for call in calls if "-g file:" in call]
    three_prime_calls = [call for call in calls if "-a file:" in call]
    # Pass A and pass B, each with a 5' and a 3' call.
    assert len(calls) == 4
    assert len(five_prime_calls) == 2
    assert len(three_prime_calls) == 2
    for call in five_prime_calls:
        assert "$final_cycle_cut" in call
    for call in three_prime_calls:
        assert "final_cycle" not in call


def final_cycle_options_block(rule_text):
    """The lines of the rule's shell command that turn the program's two words into
    Cutadapt options: from `final_cycle_cut=""` to the `fi` that closes the R2 test.
    They hold no Snakemake placeholder, so bash can run them as they are written."""
    start = rule_text.index('final_cycle_cut=""')
    r2_test = rule_text.index('if [ "$r2_final_cycle" = "remove" ]', start)
    # The first line that holds only `fi` after the R2 test closes it (a plain search
    # for "fi" would stop inside the word "final").
    closing_fi = re.search(r"^\s*fi\s*$", rule_text[r2_test:], re.MULTILINE)
    end = r2_test + closing_fi.end()
    return rule_text[start:end]


@pytest.mark.parametrize("r1_word, r2_word, expected_options", [
    ("remove", "remove", ["--cut", "-1", "-U", "-1"]),
    ("remove", "keep", ["--cut", "-1"]),
    ("keep", "remove", ["-U", "-1"]),
    ("keep", "keep", []),
])
def test_each_word_becomes_the_cut_of_its_own_read(r1_word, r2_word, expected_options):
    # The rule's own lines, run by bash in the strict mode of Snakemake's shell: the
    # word for R1 must give --cut -1 (Cutadapt's option for the first file) and the word
    # for R2 -U -1 (the second file), and "keep" nothing. printf prints the variable one
    # word per line, as Cutadapt receives it unquoted.
    script = ("set -euo pipefail\n"
              f"r1_final_cycle={r1_word}\n"
              f"r2_final_cycle={r2_word}\n"
              + final_cycle_options_block(trim_primers_text()) + "\n"
              + 'printf "%s\\n" $final_cycle_cut\n')
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == expected_options


def test_the_program_is_an_input_and_the_setting_a_param():
    text = trim_primers_text()
    assert ('final_cycle_program = workflow.basedir + "/scripts/amplicon/30a_final_cycle.py"'
            in text)
    assert PROGRAM.exists()
    assert "trim_final_cycle = TRIM_FINAL_CYCLE," in text
    assert 'r1_truncation = FINAL_CYCLE_TRUNCATION["R1"],' in text
    assert 'r2_truncation = FINAL_CYCLE_TRUNCATION["R2"],' in text
    # One call per read direction, each with its own input file and its own truncation.
    assert ("python {input.final_cycle_program} {params.trim_final_cycle} {input.r1} R1 "
            "{params.r1_truncation} 2> {log}" in text)
    assert ("python {input.final_cycle_program} {params.trim_final_cycle} {input.r2} R2 "
            "{params.r2_truncation} 2>> {log}" in text)


def test_the_log_is_opened_once_and_then_only_added_to():
    # Every redirection to the log, in order: the program's R1 line opens it (2>), and
    # everything after it must add to it (>> or 2>>), or the decision would be lost.
    redirections = re.findall(r"(2>>|2>|>>|>)\s*\{log\}", trim_primers_text())
    assert redirections[0] == "2>"
    for redirection in redirections[1:]:
        assert redirection in (">>", "2>>")
