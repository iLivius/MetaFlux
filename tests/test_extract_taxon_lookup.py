"""The taxon-name lookup of the read extraction (rule extract_taxon_reads,
workflow/rules/shotgun/55_extract_taxa.smk).

What is tested: shotgun.extract_taxa names taxa by their NCBI name, and KrakenTools
needs the taxid, so the rule's shell command first runs a small Python program that
finds the name in the sample's Kraken2 report. The program lives inside the rule's
shell string (it is passed to `python -c`), so the tests cut it out of the rule
file and run it with the interpreter, exactly as the rule does. The rules it must
follow:
  - a row counts only when its name, without its indentation, is exactly the
    requested name: Escherichia must not find Escherichia coli;
  - one matching row: its taxid is printed; no row: nothing is printed and the
    program ends normally (the rule then writes empty files);
  - several rows: a genus with a subgenus of the same name (Plasmodium, G and G1
    inside its branch) gives the genus, with a note in the log; any other case
    stops with an error and a message, because MetaFlux cannot tell which taxon
    is meant.
The second inline program of the rule, which cuts KrakenTools' progress counters
out of the log, is tested too. Both programs must contain no curly braces and no
single quotes, which Snakemake and bash would misread; a test checks that as well.
"""
import subprocess
import sys

from helpers import EXTRACT_SMK, inline_program, write_text
from make_test_kraken2_db import NODES, report_lines

LOOKUP_START = "taxid=$(python -c '"
LOOKUP_END = "' \"$target\" {input.report} {log})"
FILTER_START = "2>&1 | python -c '"
FILTER_END = "' >> {log}"

lookup_program = inline_program(EXTRACT_SMK, LOOKUP_START, LOOKUP_END)
filter_program = inline_program(EXTRACT_SMK, FILTER_START, FILTER_END)

# Reads on a few species of the test tree, written as the eight-column report
# Kraken2 produces for MetaFlux (--report-minimizer-data).
REPORT = "\n".join(report_lines(NODES, {562: 400, 4932: 100, 5833: 50, 900002: 300},
                                minimizer_columns=True, unclassified=20)) + "\n"


def run_lookup(tmp_path, target, report_text=REPORT):
    """Run the lookup as the rule does: python <program> <name> <report> <log>.
    Returns (stdout, exit code, stderr, log text). The rule creates the log before
    the lookup runs, so an empty log file is written first."""
    program = write_text(tmp_path / "lookup.py", lookup_program)
    report = write_text(tmp_path / "S1_report.txt", report_text)
    log = write_text(tmp_path / "S1_bacteria.log", "")
    result = subprocess.run([sys.executable, "-I", program, target, report, log],
                            capture_output=True, text=True)
    return result.stdout, result.returncode, result.stderr, (tmp_path / "S1_bacteria.log").read_text()


def test_programs_avoid_what_snakemake_and_bash_would_misread():
    for program in (lookup_program, filter_program):
        assert "{" not in program and "}" not in program   # Snakemake placeholders
        assert "'" not in program                           # ends the bash string
        assert '"""' not in program                         # ends the rule's shell string


def test_exact_name_gives_its_taxid_at_any_rank(tmp_path):
    stdout, code, stderr, log = run_lookup(tmp_path, "Escherichia coli")
    assert (stdout.strip(), code) == ("562", 0)
    stdout, code, stderr, log = run_lookup(tmp_path, "Bacteria")
    assert (stdout.strip(), code) == ("2", 0)
    assert log == ""                                        # nothing to note


def test_a_name_inside_a_longer_name_does_not_match(tmp_path):
    # "Escherichia" is the genus row, not the species row that contains the word.
    stdout, code, stderr, log = run_lookup(tmp_path, "Escherichia")
    assert (stdout.strip(), code) == ("561", 0)
    # A word that only occurs inside longer names matches no row.
    stdout, code, stderr, log = run_lookup(tmp_path, "coli")
    assert (stdout, code) == ("", 0)


def test_a_name_absent_from_the_report_prints_nothing(tmp_path):
    stdout, code, stderr, log = run_lookup(tmp_path, "Bacillus")
    assert (stdout, code, stderr) == ("", 0, "")


def test_genus_with_a_subgenus_of_the_same_name_gives_the_genus(tmp_path):
    # Plasmodium is a genus (G, 5820) with a subgenus Plasmodium (G1, 418103) inside
    # its branch; extracting the genus with its children covers both.
    stdout, code, stderr, log = run_lookup(tmp_path, "Plasmodium")
    assert (stdout.strip(), code) == ("5820", 0)
    assert "'Plasmodium' is the name of 2 rows" in log
    assert "extracts the rank G row (taxid 5820)" in log


def test_two_unrelated_rows_with_one_name_stop_the_job(tmp_path):
    # The same name on two branches (here a made-up genus name used twice): the
    # program exits with an error, and says so in the log and on stderr.
    two_branches = (
        "100.00\t20\t0\tR\t1\troot\n"
        "100.00\t20\t0\tR1\t131567\t  cellular organisms\n"
        " 50.00\t10\t0\tD\t2\t    Bacteria\n"
        " 50.00\t10\t10\tG\t111\t      Incertae sedis\n"
        " 50.00\t10\t0\tD\t2759\t    Eukaryota\n"
        " 50.00\t10\t10\tG\t222\t      Incertae sedis\n"
    )
    stdout, code, stderr, log = run_lookup(tmp_path, "Incertae sedis", two_branches)
    assert (stdout, code) == ("", 1)
    assert "is the name of 2 rows" in stderr
    assert "cannot tell which taxon to extract" in stderr
    assert "rank G taxid 111, rank G taxid 222" in log


def test_a_numbered_row_outside_the_genus_branch_still_stops(tmp_path):
    # A G1 row of the same name that is NOT inside the genus's branch is another
    # taxon, so the genus rule does not apply.
    report = (
        "100.00\t20\t0\tR\t1\troot\n"
        " 50.00\t10\t0\tG\t5820\t  Plasmodium\n"
        " 50.00\t10\t10\tS\t5833\t    Plasmodium falciparum\n"
        " 50.00\t10\t0\tG\t999\t  Othergenus\n"
        " 50.00\t10\t10\tG1\t998\t    Plasmodium\n"
    )
    stdout, code, stderr, log = run_lookup(tmp_path, "Plasmodium", report)
    assert (stdout, code) == ("", 1)
    assert "cannot tell which taxon to extract" in log


def test_log_filter_keeps_only_the_final_text_of_each_line(tmp_path):
    # KrakenTools redraws its counters with carriage returns; in a file every
    # redraw would stay. The filter keeps what a terminal shows when a line is done.
    program = write_text(tmp_path / "filter.py", filter_program)
    noisy = b"1000 reads processed\r2000 reads processed\r2500 reads processed\n5 reads written\nerror without newline"
    result = subprocess.run([sys.executable, "-I", program], input=noisy, capture_output=True)
    assert result.returncode == 0
    assert result.stdout == b"2500 reads processed\n5 reads written\nerror without newline\n"
