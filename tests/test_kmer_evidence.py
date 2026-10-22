"""The k-mer evidence gate: workflow/scripts/shotgun/40a_kmer_evidence.py, and the
check of its floor setting in 00_common.smk (_kmer_floor_setting).

What is tested: rule kmer_evidence reads one sample's Kraken2 report (8 columns, from
--report-minimizer-data) and judges every species by its count of distinct
minimizers. A species below the floor is removed from the report Bracken reads, and
its reads are handed to its parent. The floor is
shotgun.kmer_evidence.min_distinct_minimizers: a whole number, one floor for every
sample, or auto, a floor for each sample from its own depth:
    floor = min(6,500, (classified read pairs - Homo sapiens read pairs) / 1,000)
The tests run the whole script on small reports written on the spot and check:
  - with a number, the outputs are those of MetaFlux v2.5.0: the verdicts and their
    reasons, the gated report, the log line and the threshold line of the plot;
  - with auto: the floor from the classified read pairs, the human read pairs left
    out (the clade count of the Homo sapiens row, and a GTDB taxon that happens to
    have taxid 9606 not taken for human), the cap (also in a sample with human read
    pairs, which come off before the cap applies), a shallow sample whose floor is
    below 300, the boundary (a count equal to the floor is kept, one below it
    removed, also when the floor is not a whole number), the read floor and the 1 %
    safety rail still applied as with a number, and what the log, the evidence table
    and the plot say;
  - the setting check accepts auto, AUTO, " auto ", 0 and 300, and stops the run on
    -1, on other text and on values that YAML reads as true or false.

The plot. The script draws it with matplotlib, which the test environments (the CI
job and the metagenome environment) do not have. So every test here runs the script
with a small stand-in in its place (the `plot` fixture below): it draws nothing,
writes a placeholder file where the PNG would go, and keeps the horizontal lines the
script asks for, so that a test can read where the threshold line is and what its
legend says. The real figure is drawn in every pipeline run, in the python_utils
environment of the rule.
"""
import sys
import types
from pathlib import Path

import pytest

from helpers import COMMON_SMK, SCRIPTS, fake_snakemake, load_functions, run_script, write_text

SCRIPT = SCRIPTS / "shotgun" / "40a_kmer_evidence.py"

setting = load_functions(COMMON_SMK, ["_kmer_floor_setting"])

# The header of the evidence table, as MetaFlux v2.5.0 writes it. auto must not add,
# drop or rename a column.
EVIDENCE_HEADER = ("taxid\tname\trank\tclade_reads\tdirect_reads\ttotal_minimizers\t"
                   "distinct_minimizers\tduplicity\tminimizers_per_read\t"
                   "pct_of_classified\tdecision\treason")

# The stop message of the setting check, without the value it quotes at the end.
SETTING_STOP = ("shotgun.kmer_evidence.min_distinct_minimizers must be a whole number "
                ">= 0 (the default is 300) or the word auto")


# ── The matplotlib stand-in ───────────────────────────────────────────────────────

class StandInAxes:
    """Takes the place of a matplotlib Axes. Draws nothing; keeps each horizontal
    line as (height, legend label), which is where the script draws the floor. Only
    the calls the script makes exist here, so a new call in the script fails the
    tests with an AttributeError that names it."""

    def __init__(self):
        self.horizontal_lines = []

    def axhline(self, y, **style):
        self.horizontal_lines.append((y, style.get("label")))

    def set_title(self, *args, **style):
        pass

    def scatter(self, *args, **style):
        pass

    def axvline(self, *args, **style):
        pass

    def set_xscale(self, *args, **style):
        pass

    def set_yscale(self, *args, **style):
        pass

    def text(self, *args, **style):
        pass

    def set_xlabel(self, *args, **style):
        pass

    def set_ylabel(self, *args, **style):
        pass

    def legend(self, *args, **style):
        pass

    def grid(self, *args, **style):
        pass


class StandInFigure:
    """Takes the place of a matplotlib Figure: savefig writes a placeholder file, so
    the rule's third output exists as it would after a real run."""

    def tight_layout(self):
        pass

    def savefig(self, path, **options):
        Path(path).write_text("stand-in for the evidence plot\n")


@pytest.fixture(autouse=True)
def plot(monkeypatch):
    """Put the matplotlib stand-in in place for one test (autouse: every test of this
    file gets it, so none can reach a real matplotlib by accident) and return its
    axes, which keep what the script drew. monkeypatch puts back whatever was there
    before once the test ends."""
    axes = StandInAxes()

    def subplots(**options):
        return StandInFigure(), axes

    def close(figure):
        pass

    def use(backend):
        pass

    pyplot = types.ModuleType("matplotlib.pyplot")
    pyplot.subplots = subplots
    pyplot.close = close
    matplotlib = types.ModuleType("matplotlib")
    matplotlib.use = use
    matplotlib.pyplot = pyplot
    monkeypatch.setitem(sys.modules, "matplotlib", matplotlib)
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", pyplot)
    return axes


# ── Reports and runs ──────────────────────────────────────────────────────────────

def report_line(clade, direct, distinct, rank, taxid, depth, name, sample_total):
    """One line of a Kraken2 report with --report-minimizer-data: percentage of the
    sample, clade reads, reads on the node itself, total and distinct minimizers,
    rank code, taxid, and the name indented two spaces per level of the tree. Total
    minimizers are set to three times the distinct ones (a duplicity of 3)."""
    percentage = clade / sample_total * 100
    return (f"{percentage:6.2f}\t{clade}\t{direct}\t{3 * distinct}\t{distinct}\t{rank}\t"
            f"{taxid}\t{'  ' * depth}{name}")


def small_report(classified_pairs, species, human_pairs=0, human_distinct=50000):
    """The text of a small Kraken2 report, shaped like a real one:
        unclassified
        root
          Bacteria
            Bacillus (genus)
              the species, as (taxid, name, read_pairs, distinct_minimizers)
          Eukaryota > Homo > Homo sapiens > Homo sapiens neanderthalensis
                                            (only when human_pairs > 0)
    The root's clade reads are `classified_pairs`; the read pairs not under one of the
    species or under human stay on the root itself, as Kraken2 leaves the reads it
    cannot place lower. Every clade count is its own reads plus its children's, as in
    a real report. Homo sapiens is given many distinct minimizers by default, so that
    it passes every floor and does not get in the way of what a test checks.
    `human_pairs` is the clade count of the Homo sapiens row, the number auto leaves
    out. A tenth of those read pairs sit on the sub-rank row below it, so that the
    row's clade reads and its own reads differ, and a script that took the wrong one
    of the two fails the tests. (In the PlusPF reports of the test dataset, Homo
    sapiens has no row below it, so there the two are equal.)"""
    unclassified = 100000
    sample_total = unclassified + classified_pairs
    genus_reads = sum(read_pairs for taxid, name, read_pairs, distinct in species)
    root_own_reads = classified_pairs - genus_reads - human_pairs
    assert root_own_reads >= 0, "the species hold more read pairs than the sample"
    lines = [
        report_line(unclassified, unclassified, 0, "U", "0", 0, "unclassified", sample_total),
        report_line(classified_pairs, root_own_reads, 900000, "R", "1", 0, "root", sample_total),
        report_line(genus_reads, 0, 800000, "D", "2", 1, "Bacteria", sample_total),
        report_line(genus_reads, 0, 700000, "G", "1386", 2, "Bacillus", sample_total),
    ]
    for taxid, name, read_pairs, distinct in species:
        lines.append(report_line(read_pairs, read_pairs, distinct, "S", taxid, 3, name,
                                 sample_total))
    if human_pairs > 0:
        human_pairs_below = human_pairs // 10          # on the sub-rank row; see above
        lines.append(report_line(human_pairs, 0, human_distinct, "D", "2759", 1,
                                 "Eukaryota", sample_total))
        lines.append(report_line(human_pairs, 0, human_distinct, "G", "9605", 2, "Homo",
                                 sample_total))
        lines.append(report_line(human_pairs, human_pairs - human_pairs_below,
                                 human_distinct, "S", "9606", 3, "Homo sapiens",
                                 sample_total))
        lines.append(report_line(human_pairs_below, human_pairs_below, human_distinct // 10,
                                 "S1", "63221", 4, "Homo sapiens neanderthalensis",
                                 sample_total))
    return "\n".join(lines) + "\n"


def run_gate(tmp_path, report_text, floor, min_reads=0):
    """Run the script on one report, as rule kmer_evidence runs it for one sample:
    rank S, the safety rail at 1 %, the gate on (which only labels the log and the
    plot title; the gated report is written either way), the floor as 00_common.smk
    hands it over (an int, or the text "auto"). Returns the run's outputs:
      verdicts   {taxid: (decision, reason)}, read back from the evidence table
      header     the first line of the evidence table
      gated      the text of the gated report
      log        the text of the log"""
    report = write_text(tmp_path / "S1_report.txt", report_text)
    sm = fake_snakemake(
        input={"report": report},
        output={"evidence": str(tmp_path / "S1_species_evidence.tsv"),
                "gated": str(tmp_path / "S1_report_gated.txt"),
                "plot": str(tmp_path / "S1_evidence.png")},
        params={"sample": "S1", "rank": "S", "min_distinct": floor, "min_reads": min_reads,
                "protect_frac": 0.01, "enabled": True},
        log=[str(tmp_path / "S1.log")],
    )
    run_script(SCRIPT, sm)

    evidence_lines = (tmp_path / "S1_species_evidence.tsv").read_text().splitlines()
    verdicts = {}
    for line in evidence_lines[1:]:
        fields = line.split("\t")
        verdicts[fields[0]] = (fields[10], fields[11])
    return types.SimpleNamespace(
        verdicts=verdicts,
        header=evidence_lines[0],
        gated=(tmp_path / "S1_report_gated.txt").read_text(),
        log=(tmp_path / "S1.log").read_text(),
    )


def log_line(log_text, start):
    """The one line of a kmer_evidence log that starts with `start` (after the
    [kmer_evidence] tag)."""
    found = [line for line in log_text.splitlines()
             if line.startswith(f"[kmer_evidence] {start}")]
    assert len(found) == 1, f"expected one log line starting with {start!r}, got {found}"
    return found[0]


# ── A number: the floor behaves as in MetaFlux v2.5.0 ────────────────────────────

# A hand-written report with a strain under the species that fails, so that the whole
# subtree has to go, and the clade reads go to the genus's own count.
REPORT_WITH_STRAIN = (
    " 33.33\t50000\t50000\t0\t0\tU\t0\tunclassified\n"
    " 66.67\t100000\t99150\t3000000\t300000\tR\t1\troot\n"
    "  0.57\t850\t50\t8500\t2900\tD\t2\t  Bacteria\n"
    "  0.53\t800\t30\t8000\t2800\tG\t1386\t    Bacillus\n"
    "  0.33\t500\t500\t5000\t1200\tS\t1423\t      Bacillus subtilis\n"
    "  0.17\t250\t240\t2500\t120\tS\t1428\t      Bacillus thuringiensis\n"
    "  0.01\t10\t10\t100\t50\tS1\t1000001\t        Bacillus thuringiensis strain X\n"
    "  0.01\t20\t20\t200\t400\tS\t1396\t      Bacillus cereus\n"
)

# What the gate makes of it at the floor of 300: B. thuringiensis (120 distinct
# minimizers) and its strain are removed, and its 250 clade reads are added to the
# genus's own count (30 + 250 = 280). Every other line is unchanged.
GATED_WITH_STRAIN = (
    " 33.33\t50000\t50000\t0\t0\tU\t0\tunclassified\n"
    " 66.67\t100000\t99150\t3000000\t300000\tR\t1\troot\n"
    "  0.57\t850\t50\t8500\t2900\tD\t2\t  Bacteria\n"
    "  0.53\t800\t280\t8000\t2800\tG\t1386\t    Bacillus\n"
    "  0.33\t500\t500\t5000\t1200\tS\t1423\t      Bacillus subtilis\n"
    "  0.01\t20\t20\t200\t400\tS\t1396\t      Bacillus cereus\n"
)


def test_a_number_gives_the_outputs_of_v2_5_0(tmp_path, plot):
    run = run_gate(tmp_path, REPORT_WITH_STRAIN, floor=300)
    assert run.verdicts == {
        "1423": ("keep", "passed both floors"),
        "1428": ("prune", "distinct minimizers 120 < 300"),
        "1396": ("keep", "passed both floors"),
    }
    assert run.gated == GATED_WITH_STRAIN
    assert run.header == EVIDENCE_HEADER
    # The settings line reads as it always did, and a number writes no auto line.
    assert log_line(run.log, "sample=") == (
        "[kmer_evidence] sample=S1 rank=S min_distinct_minimizers=300 min_reads=0 "
        "never_prune_above_fraction=0.01")
    assert "auto floor" not in run.log
    # The plot: one threshold line, at exactly the configured number.
    assert plot.horizontal_lines == [(300, "threshold = 300")]
    assert type(plot.horizontal_lines[0][0]) is int


# ── auto: the floor from the sample's own depth ──────────────────────────────────

def test_auto_floor_is_one_per_1000_classified_pairs(tmp_path, plot):
    # 2,500,000 classified read pairs and no human row: the floor is 2,500.
    report = small_report(2_500_000, [
        ("1423", "Bacillus subtilis", 5000, 40000),
        ("1428", "Bacillus thuringiensis", 900, 2500),     # exactly on the floor
        ("1396", "Bacillus cereus", 900, 2499),            # one below it
    ])
    run = run_gate(tmp_path, report, floor="auto")
    assert run.verdicts == {
        "1423": ("keep", "passed both floors"),
        "1428": ("keep", "passed both floors"),
        "1396": ("prune", "distinct minimizers 2499 < 2500 (auto)"),
    }
    # The evidence table keeps its columns; the removed species leaves the gated report.
    assert run.header == EVIDENCE_HEADER
    assert "Bacillus cereus" not in run.gated
    assert "Bacillus thuringiensis" in run.gated
    # The log names the setting, then says how this sample's floor was worked out.
    assert "min_distinct_minimizers=auto" in log_line(run.log, "sample=")
    assert log_line(run.log, "auto floor:") == (
        "[kmer_evidence] auto floor: 2500000 classified read pairs, no Homo sapiens row "
        "in the report, so no read pairs left out: 2500000 / 1000 = 2500.000. A taxon "
        "at rank S needs at least 2500 distinct minimizers.")
    # The plot draws this sample's floor and says it is auto.
    assert plot.horizontal_lines == [(2500.0, "threshold = 2500 (auto)")]


def test_auto_leaves_the_human_read_pairs_out(tmp_path, plot):
    # 2,500,000 classified pairs, of which 300,000 are under Homo sapiens (270,000 on
    # the row itself and 30,000 on the row below it, see small_report): the floor is
    # (2,500,000 - 300,000) / 1,000 = 2,200. Not 2,500, which leaving nothing out would
    # give, and not 2,230, which taking the row's own reads instead of its clade
    # reads would give.
    report = small_report(2_500_000, [
        ("1423", "Bacillus subtilis", 900, 2300),          # kept only because of the subtraction
        ("1396", "Bacillus cereus", 900, 2199),
    ], human_pairs=300_000)
    run = run_gate(tmp_path, report, floor="auto")
    assert run.verdicts["1423"] == ("keep", "passed both floors")
    assert run.verdicts["1396"] == ("prune", "distinct minimizers 2199 < 2200 (auto)")
    assert log_line(run.log, "auto floor:") == (
        "[kmer_evidence] auto floor: 2500000 classified read pairs, less 300000 Homo "
        "sapiens read pairs (taxid 9606): 2200000 / 1000 = 2200.000. A taxon at rank S "
        "needs at least 2200 distinct minimizers.")
    assert plot.horizontal_lines == [(2200.0, "threshold = 2200 (auto)")]


def test_auto_takes_only_the_row_named_homo_sapiens_for_human(tmp_path):
    # In a GTDB-based Kraken2 index the taxids are not NCBI's: in
    # k2_gtdb_genome_reps_20250609 (GTDB r226), taxid 9606 is the archaeal genus
    # Thermococcus_A. Its 300,000 read pairs are not human and must stay in the
    # count, so the floor here is 2,500 (all classified pairs), not 2,200.
    sample_total = 2_600_000          # unclassified + classified
    report = "\n".join([
        report_line(100_000, 100_000, 0, "U", "0", 0, "unclassified", sample_total),
        report_line(2_500_000, 2_199_100, 900000, "R", "1", 0, "root", sample_total),
        report_line(900, 0, 800000, "D", "2", 1, "Bacteria", sample_total),
        report_line(900, 0, 700000, "G", "1386", 2, "Bacillus", sample_total),
        report_line(900, 900, 2300, "S", "1396", 3, "Bacillus cereus", sample_total),
        report_line(300_000, 0, 9000, "D", "2157", 1, "Archaea", sample_total),
        report_line(300_000, 300_000, 9000, "G", "9606", 2, "Thermococcus_A", sample_total),
    ]) + "\n"
    run = run_gate(tmp_path, report, floor="auto")
    assert run.verdicts["1396"] == ("prune", "distinct minimizers 2300 < 2500 (auto)")
    assert "no Homo sapiens row in the report" in log_line(run.log, "auto floor:")


def test_auto_floor_is_capped_at_6500(tmp_path, plot):
    # 10 million classified pairs would give 10,000; the cap holds the floor at 6,500.
    report = small_report(10_000_000, [
        ("1423", "Bacillus subtilis", 3000, 9000),         # removed by an uncapped floor
        ("1428", "Bacillus thuringiensis", 3000, 6500),
        ("1396", "Bacillus cereus", 3000, 6499),
    ])
    run = run_gate(tmp_path, report, floor="auto")
    assert run.verdicts == {
        "1423": ("keep", "passed both floors"),
        "1428": ("keep", "passed both floors"),
        "1396": ("prune", "distinct minimizers 6499 < 6500 (auto)"),
    }
    assert log_line(run.log, "auto floor:") == (
        "[kmer_evidence] auto floor: 10000000 classified read pairs, no Homo sapiens row "
        "in the report, so no read pairs left out: 10000000 / 1000 = 10000.000, above "
        "the cap of 6500, so the floor is 6500. A taxon at rank S needs at least 6500 "
        "distinct minimizers.")
    assert plot.horizontal_lines == [(6500, "threshold = 6500 (auto)")]


def test_auto_takes_the_human_read_pairs_off_before_the_cap(tmp_path, plot):
    # The order matters for a deep sample with much human DNA and host removal off.
    # 10 million classified pairs, 2 million of them human: (10,000,000 - 2,000,000)
    # / 1,000 = 8,000, above the cap, so the floor is 6,500. Capping first and taking
    # the human pairs off afterwards would give 6,500 - 2,000 = 4,500, and keep the
    # species with 6,499 distinct minimizers.
    report = small_report(10_000_000, [
        ("1428", "Bacillus thuringiensis", 3000, 6500),
        ("1396", "Bacillus cereus", 3000, 6499),
    ], human_pairs=2_000_000)
    run = run_gate(tmp_path, report, floor="auto")
    assert run.verdicts["1428"] == ("keep", "passed both floors")
    assert run.verdicts["1396"] == ("prune", "distinct minimizers 6499 < 6500 (auto)")
    assert log_line(run.log, "auto floor:") == (
        "[kmer_evidence] auto floor: 10000000 classified read pairs, less 2000000 Homo "
        "sapiens read pairs (taxid 9606): 8000000 / 1000 = 8000.000, above the cap of "
        "6500, so the floor is 6500. A taxon at rank S needs at least 6500 distinct "
        "minimizers.")
    assert plot.horizontal_lines == [(6500, "threshold = 6500 (auto)")]


def test_auto_floor_of_a_shallow_sample_is_below_300(tmp_path):
    # 150,000 classified pairs: the floor is 150, lower than the default 300. A
    # species with 200 distinct minimizers is kept with auto and removed at 300.
    report = small_report(150_000, [
        ("1423", "Bacillus subtilis", 400, 200),
        ("1396", "Bacillus cereus", 400, 149),
    ])
    run = run_gate(tmp_path / "auto", report, floor="auto")
    assert run.verdicts == {
        "1423": ("keep", "passed both floors"),
        "1396": ("prune", "distinct minimizers 149 < 150 (auto)"),
    }
    run_300 = run_gate(tmp_path / "fixed", report, floor=300)
    assert run_300.verdicts["1423"] == ("prune", "distinct minimizers 200 < 300")


def test_auto_boundary_a_count_on_the_floor_passes_one_below_fails(tmp_path):
    # A whole floor: 2,268,000 pairs -> 2,268. 2,268 passes, 2,267 does not.
    report = small_report(2_268_000, [
        ("1423", "Bacillus subtilis", 900, 2268),
        ("1396", "Bacillus cereus", 900, 2267),
    ])
    run = run_gate(tmp_path / "whole", report, floor="auto")
    assert run.verdicts["1423"] == ("keep", "passed both floors")
    assert run.verdicts["1396"] == ("prune", "distinct minimizers 2267 < 2268 (auto)")

    # A floor with a fraction: 2,268,300 pairs -> 2,268.3. The decision uses the
    # fraction, as the benchmark did, so 2,268 is just below and fails; the reason
    # shows the floor rounded UP (2,269), the smallest count that passes. A fraction
    # below one half is chosen on purpose: rounding to the nearest number would give
    # 2,268, so a floor rounded the wrong way, in the decision or in the reason,
    # fails this test.
    report = small_report(2_268_300, [
        ("1423", "Bacillus subtilis", 900, 2269),
        ("1396", "Bacillus cereus", 900, 2268),
    ])
    run = run_gate(tmp_path / "fraction", report, floor="auto")
    assert run.verdicts["1423"] == ("keep", "passed both floors")
    assert run.verdicts["1396"] == ("prune", "distinct minimizers 2268 < 2269 (auto)")
    assert "= 2268.300. A taxon at rank S needs at least 2269" in log_line(run.log, "auto floor:")


def test_auto_keeps_the_read_floor(tmp_path):
    # min_reads applies together with the auto floor, as with a number: 5,000
    # distinct minimizers clear the floor of 1,000, but 150 read pairs are below 200.
    report = small_report(1_000_000, [
        ("1423", "Bacillus subtilis", 250, 5000),
        ("1396", "Bacillus cereus", 150, 5000),
    ])
    run = run_gate(tmp_path, report, floor="auto", min_reads=200)
    assert run.verdicts["1423"] == ("keep", "passed both floors")
    assert run.verdicts["1396"] == ("prune", "clade reads 150 < 200")


def test_auto_keeps_the_safety_rail_on_all_classified_reads(tmp_path):
    # 1,000,000 classified pairs, 600,000 of them human: the floor is 400. Both
    # species are far below it. The rail still keeps a taxon holding more than 1 % of
    # the sample, and that share is of ALL classified reads, human ones included, as
    # with a number: 20,000 pairs (2 %) are kept, 8,000 pairs (0.8 %) are not, although
    # they would be 2 % of the non-human pairs.
    report = small_report(1_000_000, [
        ("1423", "Bacillus subtilis", 20_000, 100),
        ("1396", "Bacillus cereus", 8_000, 100),
    ], human_pairs=600_000)
    run = run_gate(tmp_path, report, floor="auto")
    assert run.verdicts["1423"] == ("keep", "protected: 2.00% of classified reads")
    assert run.verdicts["1396"] == ("prune", "distinct minimizers 100 < 400 (auto)")


# ── The setting check in 00_common.smk ───────────────────────────────────────────

def test_setting_accepts_auto_in_any_capitals_and_with_spaces():
    for value in ["auto", "AUTO", " auto ", "Auto"]:
        assert setting._kmer_floor_setting(value) == "auto"


def test_setting_hands_a_number_over_as_the_same_int_as_before():
    # The rule's params must not change for a number, or Snakemake would run the step
    # again: 0 and 300 come out as the ints 0 and 300, as in v2.5.0, and so do the
    # spellings v2.5.0 also read as 300.
    for value, expected in [(0, 0), (300, 300), (300.0, 300), ("300", 300)]:
        result = setting._kmer_floor_setting(value)
        assert result == expected
        assert type(result) is int


@pytest.mark.parametrize("value", [-1, "-1", "automatic", "high", "", True, False, 300.5, None])
def test_setting_stops_on_anything_else(value):
    with pytest.raises(SystemExit) as stop:
        setting._kmer_floor_setting(value)
    message = str(stop.value)
    assert SETTING_STOP in message
    assert message.endswith(f"got {value!r}")
