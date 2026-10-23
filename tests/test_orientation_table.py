"""The orientation table of primer trimming under amplicon.primers.orientation: mixed
(rule trim_primers, workflow/rules/amplicon/30_preprocess.smk).

What is tested, and why it matters:
  - _amplicon_targets() in 00_common.smk, the list of files a plain run asks for,
    holds stats/orientation/{sample}.tsv for every sample under mixed, and nothing
    else changes between the two settings. The file has a second job besides showing
    how each sample's pairs split between the two orientations: a folder trimmed with
    mixed by MetaFlux v2.5.0 or earlier lacks it, and only because the run asks for
    it does primer trimming run again there (and with it the orientation fix), with
    or without Snakemake's records of that folder. Under fixed the list must stay as
    it was, or every fixed folder would run something again.
  - The small Python program inside the rule's shell command that writes the table
    from Cutadapt's JSON reports of the two 5' calls. It lives in the shell string,
    so the tests cut it out of the rule file (inline_program() in helpers.py) and run
    it with the interpreter, as the rule does. Its code must contain no curly braces
    and no single quotes, which Snakemake and bash would misread; a test checks that.

_amplicon_targets() reads the run's settings from names that 00_common.smk sets while
it loads the config (OUT, SAMPLES, ORIENTATION, ...). load_functions() runs only the
function's definition, so the tests set those names by hand, in the namespace the
function was loaded into, and give it a stand-in for Snakemake's expand() (the unit
tests do not need Snakemake).
"""
import itertools
import json
import subprocess
import sys
from pathlib import Path

from helpers import COMMON_SMK, WORKFLOW, inline_program, load_functions, write_text

PREPROCESS_SMK = WORKFLOW / "rules" / "amplicon" / "30_preprocess.smk"

TABLE_START = "python -c '"
TABLE_END = "' {output.json_pa5} $tmpdir/"
table_program = inline_program(PREPROCESS_SMK, TABLE_START, TABLE_END)

SAMPLES = ["S1", "S2"]
TABLES = {"/out/stats/orientation/S1.tsv", "/out/stats/orientation/S2.tsv"}


def expand_stand_in(pattern, **values):
    """Snakemake's expand() for the plain case _amplicon_targets() uses: the pattern
    filled in with every combination of the given values, as a list of paths."""
    names = list(values)
    filled = []
    for combination in itertools.product(*(values[name] for name in names)):
        filled.append(str(pattern).format(**dict(zip(names, combination))))
    return filled


def amplicon_targets(orientation, remove_phix=True, phylogeny=False):
    """The files a plain amplicon run of the two samples S1 and S2 asks for, writing to
    /out, with the given orientation setting, as a set of path strings."""
    loaded = load_functions(COMMON_SMK, ["_amplicon_targets"])
    settings = loaded._amplicon_targets.__globals__
    settings["OUT"] = Path("/out")
    settings["SAMPLES"] = SAMPLES
    settings["ORIENTATION"] = orientation
    settings["REMOVE_PHIX"] = remove_phix
    settings["FALCO_STAGES"] = ["raw", "nophix", "stripped"] if remove_phix else ["raw", "stripped"]
    settings["PHYLO_ENABLED"] = phylogeny
    # Off: the COI screens' tables are asked for in COI runs only, and none is tested here.
    settings["COI_SCREENS"] = False
    settings["expand"] = expand_stand_in
    return {str(path) for path in loaded._amplicon_targets()}


def test_mixed_asks_for_one_orientation_table_per_sample():
    targets = amplicon_targets("mixed")
    tables = {path for path in targets if "/stats/orientation/" in path}
    assert tables == TABLES


def test_fixed_asks_for_no_orientation_table():
    targets = amplicon_targets("fixed")
    assert not [path for path in targets if "/stats/orientation/" in path]


def test_the_tables_are_the_only_difference_between_the_two_settings():
    # With and without PhiX removal and the phylogeny, so that the comparison does not
    # rest on one combination of the other optional targets.
    for remove_phix, phylogeny in [(True, False), (False, True)]:
        mixed = amplicon_targets("mixed", remove_phix, phylogeny)
        fixed = amplicon_targets("fixed", remove_phix, phylogeny)
        assert mixed - fixed == TABLES
        assert fixed - mixed == set()


def test_program_avoids_what_snakemake_and_bash_would_misread():
    assert "{" not in table_program and "}" not in table_program   # Snakemake placeholders
    assert "'" not in table_program                                 # ends the bash string
    assert '"""' not in table_program                               # ends the rule's shell string


def write_report(path, pairs_in, pairs_kept):
    """A Cutadapt JSON report with the read counts the program reads (a real report has
    many more fields, which the program ignores)."""
    report = {"tag": "Cutadapt report",
              "read_counts": {"input": pairs_in, "output": pairs_kept}}
    return write_text(path, json.dumps(report))


def run_table_program(tmp_path, pass_a, pass_b, sample="S1"):
    """Run the program as the rule does: python <program> <pass A report> <pass B report>
    <sample> <table>. pass_a and pass_b are (pairs in, pairs kept). Returns the table's
    lines, split at the tabs."""
    program = write_text(tmp_path / "orientation_table.py", table_program)
    report_a = write_report(tmp_path / "passA_5prime.json", *pass_a)
    report_b = write_report(tmp_path / "passB_5prime.json", *pass_b)
    table = tmp_path / "stats" / "orientation" / f"{sample}.tsv"
    table.parent.mkdir(parents=True)
    result = subprocess.run([sys.executable, "-I", program, report_a, report_b, sample, str(table)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return [line.split("\t") for line in table.read_text().splitlines()]


def test_table_gives_the_pairs_each_pass_kept_and_their_shares(tmp_path):
    # 300 pairs went in; pass A kept 146 as sequenced, pass B 147 with mates swapped.
    lines = run_table_program(tmp_path, pass_a=(300, 146), pass_b=(300, 147))
    assert lines == [
        ["sample", "read_pairs_in", "kept_as_sequenced", "pct_as_sequenced",
         "kept_swapped", "pct_swapped"],
        ["S1", "300", "146", "48.67", "147", "49.00"],
    ]


def test_table_gives_na_shares_for_a_sample_without_pairs(tmp_path):
    lines = run_table_program(tmp_path, pass_a=(0, 0), pass_b=(0, 0))
    assert lines[1] == ["S1", "0", "0", "NA", "0", "NA"]
