"""The ASV length window of rule dada_length_filter
(workflow/scripts/amplicon/60h_dada_length_filter.py).

What is tested, and why it matters:
  - ITS with extraction on (ITSx ran before the filter): the floor of the auto window is
    the probe's min, the shortest ITS region in the reference, without window_margin; the
    ceiling stays p95 + window_margin. Before, the floor was q1 - window_margin, 125 bp on
    the UNITE ITS2 reference, and it removed the yeasts with the shortest ITS2 of the ITS2
    test set (63-78 bp after ITSx: Candida haemulonii, Metschnikowia pulcherrima and a
    Dipodascus), with 64 % of the reads left after chimera removal. Checked for ITS2 and
    ITS1: the window, the ASVs kept and dropped (both edges are kept), the matching count
    table, the source written to asv_length_stats.json, the log line that names the rule,
    and the two window lines of the histogram. window_margin then moves only the ceiling.
  - Every other case keeps its window and its source text as before: 16S (a PCR-mode
    probe, with Metaxa2 extraction on), 18S, gyrB (a direct-mode probe, like ITS), ITS
    with extraction off (its ASVs then still carry their rRNA flanks), a manual range on
    an ITS run, and the fallback to expected_length when no probe ran.
  - Settings that leave nothing to size the window from stop the step with one ERROR line
    and exit 1, not a Python traceback: mode auto with neither a probe nor an expected
    length in trunclen.json (it is null when pick_trunclen could not resolve one), and a
    mode that is neither auto nor given a range.

The script is run whole with a stand-in snakemake object (run_script in helpers.py), on
small files written to a temporary folder: the ASV FASTA before and after extraction, the
count table, trunclen.json and a probe JSON as 40a_amplicon_probe.py writes it. The ITS
probes carry the statistics of the UNITE files MetaFlux downloads; the others are probes
cached by earlier test runs. Only the length of an ASV matters to the filter, so each test
ASV is a run of ACGT of the length the test needs. The script ends with sys.exit(main()),
so a run that passes ends in SystemExit with code 0, which run_filter checks.
The histogram is drawn with matplotlib, which the test environments (the CI job and the
metagenome environment) do not have. A stand-in takes its place, as in
test_kmer_evidence.py: it draws nothing, writes a placeholder file where the PNG would go,
and keeps the vertical lines the script draws at the two edges of the window.
"""
import json
import sys
import types
from pathlib import Path

import pytest

from helpers import SCRIPTS, fake_snakemake, run_script, write_text

SCRIPT = SCRIPTS / "amplicon" / "60h_dada_length_filter.py"

# ── Probe results ────────────────────────────────────────────────────────────────
# As 40a_amplicon_probe.py writes them (the fields the filter reads, and the reference
# tag the ITS log line names). The ITS values are those of the UNITE UCHIME files MetaFlux
# downloads; the ITS page of the manual gives the same table.
PROBE_ITS2 = {"amplicon_type": "ITS", "probe_mode": "direct",
              "reference_tag": "unite_uchime_ITS2", "n_amplicons": 154918,
              "min": 38, "q1": 175, "median": 204, "q3": 233, "p95": 293, "p99": 401,
              "max": 4461}
PROBE_ITS1 = {"amplicon_type": "ITS", "probe_mode": "direct",
              "reference_tag": "unite_uchime_ITS1", "n_amplicons": 160204,
              "min": 17, "q1": 174, "median": 200, "q3": 230, "p95": 293, "p99": 402,
              "max": 2427}
# The other markers: in each, min is far from q1 - 50, so a floor taken from min would
# show. 16S is 799F/1175R on SILVA 144 (the window [324, 435] of the 2.7.0 release
# notes); gyrB is read straight off DD7RZ8's amplicons, a direct-mode probe like ITS.
PROBE_16S = {"amplicon_type": "16S", "probe_mode": "pcr", "n_amplicons": 573114,
             "min": 38, "q1": 374, "median": 376, "q3": 379, "p95": 385, "p99": 694,
             "max": 1855}
PROBE_18S = {"amplicon_type": "18S", "probe_mode": "pcr", "n_amplicons": 75486,
             "min": 53, "q1": 379, "median": 381, "q3": 388, "p95": 494, "p99": 618,
             "max": 2187}
PROBE_GYRB = {"amplicon_type": "gyrB", "probe_mode": "direct", "n_amplicons": 118719,
              "min": 234, "q1": 244, "median": 247, "q3": 250, "p95": 250, "p99": 262,
              "max": 270}

# The source text of the ordinary auto window, which no case but ITS after ITSx may change.
Q1_SOURCE = "auto:probe_q1_p95+/-50"


# ── The matplotlib stand-in ───────────────────────────────────────────────────────

class StandInAxes:
    """Takes the place of a matplotlib Axes. Draws nothing; keeps each vertical line as
    (position, legend label), which is where the script marks the two window edges. Only
    the calls the script makes exist here, so a new call in the script fails the tests
    with an AttributeError that names it."""

    def __init__(self):
        self.vertical_lines = []

    def hist(self, *args, **style):
        pass

    def axvline(self, x, **style):
        self.vertical_lines.append((x, style.get("label")))

    def set_xlabel(self, *args, **style):
        pass

    def set_ylabel(self, *args, **style):
        pass

    def set_title(self, *args, **style):
        pass

    def legend(self, *args, **style):
        pass


class StandInFigure:
    """Takes the place of a matplotlib Figure: savefig writes a placeholder file, so the
    rule's histogram output exists as it would after a real run."""

    def tight_layout(self):
        pass

    def savefig(self, path, **options):
        Path(path).write_text("stand-in for the ASV length histogram\n")


@pytest.fixture(autouse=True)
def plot(monkeypatch):
    """Put the matplotlib stand-in in place for one test (autouse: every test of this file
    gets it, so none can reach a real matplotlib by accident) and return its axes, which
    keep the window lines the script drew. monkeypatch puts back whatever was there before
    once the test ends."""
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


# ── Input files and one run of the filter ─────────────────────────────────────────

def asv_fasta(lengths, extra_bases=0):
    """FASTA text with one ASV per length, named ASV_1, ASV_2, ... in the order given.
    `extra_bases` lengthens every ASV, as the rRNA flanks lengthen an ASV before ITSx or
    Metaxa2 trims them off."""
    lines = []
    for number, length in enumerate(lengths, start=1):
        total = length + extra_bases
        bases = ("ACGT" * (total // 4 + 1))[:total]
        lines.append(f">ASV_{number}")
        lines.append(bases)
    return "\n".join(lines) + "\n"


def count_table(n_asvs):
    """The name-keyed count table as dada_seqtab.R (or target_extract) writes it: a header
    of quoted ASV names after an empty first cell, then one row per sample. Two samples;
    ASV n has n reads in the first and 10 * n in the second."""
    header = ['""'] + [f'"ASV_{number}"' for number in range(1, n_asvs + 1)]
    first = ['"S1"'] + [str(number) for number in range(1, n_asvs + 1)]
    second = ['"S2"'] + [str(10 * number) for number in range(1, n_asvs + 1)]
    rows = [header, first, second]
    return "\n".join("\t".join(row) for row in rows) + "\n"


def run_filter(tmp_path, amp_type, extraction_enabled, lengths, probe=None,
               window_margin=50, manual_range=None, expected_length=None,
               mode="auto", expected_exit=0):
    """Run 60h_dada_length_filter.py once, as rule dada_length_filter would, on ASVs of
    the given `lengths` (the lengths the filter judges: after extraction when it ran).
    `probe` is the probe JSON as a dict, or None for a run without one (expected_length
    set by hand, so the rule hands an empty list). `expected_length` goes into
    trunclen.json, as pick_trunclen writes it (null for ITS).
    Returns (the filter_window of asv_length_stats.json, the names of the kept ASVs,
    the ASV columns of the kept count table, the log text). With an `expected_exit`
    other than 0 (a run that must stop), only the log text is filled in."""
    dada2_dir = tmp_path / "05.dada2"
    # Before extraction the ASVs still carry their rRNA flanks (140 bases here); when
    # extraction is off, the filter judges these very sequences.
    if extraction_enabled:
        seqs_pre = write_text(dada2_dir / "seqs.fasta", asv_fasta(lengths, extra_bases=140))
        seqs = write_text(dada2_dir / "seqs_extracted.fasta", asv_fasta(lengths))
    else:
        seqs_pre = write_text(dada2_dir / "seqs.fasta", asv_fasta(lengths))
        seqs = seqs_pre
    seqtab_names = write_text(dada2_dir / "seqtab_in_head_names.txt", count_table(len(lengths)))
    trunclen = {"amplicon_type": amp_type, "expected_length": expected_length}
    trunclen_json = write_text(tmp_path / "stats" / "trunclen.json", json.dumps(trunclen))
    probe_json = []
    if probe is not None:
        probe_json = [write_text(tmp_path / "refdb" / "cache" / "probe.json", json.dumps(probe))]

    stats_json = tmp_path / "stats" / "dada2" / "asv_length_stats.json"
    kept_fasta = dada2_dir / "seqs_lenfilt.fasta"
    kept_table = dada2_dir / "seqtab_lenfilt_head_names.txt"
    log_file = tmp_path / "logs" / "dada_length_filter.log"
    sm = fake_snakemake(
        input={"seqs_pre": seqs_pre, "seqs": seqs, "seqtab_names": seqtab_names,
               "trunclen_json": trunclen_json, "probe_json": probe_json},
        output={"seqs": str(kept_fasta), "seqtab_names": str(kept_table),
                "seqtab_seqs": str(dada2_dir / "seqtab_lenfilt_head_seqs.txt"),
                "stats_json": str(stats_json),
                "hist_png": str(tmp_path / "stats" / "dada2" / "asv_length_hist.png")},
        params={"mode": mode, "window_margin": window_margin,
                "manual_range": manual_range, "amp_type": amp_type,
                "extraction_enabled": extraction_enabled},
        log=[str(log_file)],
    )
    with pytest.raises(SystemExit) as finished:
        run_script(SCRIPT, sm)
    assert finished.value.code == expected_exit, log_file.read_text()
    if expected_exit != 0:
        return None, None, None, log_file.read_text()

    window = json.loads(stats_json.read_text())["filter_window"]
    kept_names = [line[1:] for line in kept_fasta.read_text().splitlines()
                  if line.startswith(">")]
    table_header = kept_table.read_text().splitlines()[0].split("\t")
    table_columns = [cell.strip('"') for cell in table_header[1:]]
    return window, kept_names, table_columns, log_file.read_text()


# ── ITS after ITSx: the floor is the probe's min ──────────────────────────────────

def test_its2_floor_is_the_shortest_reference_region(tmp_path, plot):
    """UNITE ITS2: the window is [38, 343] (it was [125, 343]). The ASVs of 63, 70 and
    78 bp stand for the short-ITS2 yeasts of the test set, and the one of 124 bp for an
    ITS2 just under the old floor: all are kept now. Below 38 and above 343 bp, ASVs
    are still dropped; 38 and 343 themselves are kept."""
    lengths = [30, 37, 38, 63, 70, 78, 124, 211, 343, 344]
    window, kept, columns, log = run_filter(tmp_path, "ITS", True, lengths, probe=PROBE_ITS2)

    assert window == {"min": 38, "max": 343, "source": "auto:ITS_probe_min_p95+50"}
    # ASV_3 to ASV_9 are the lengths 38 to 343; ASV_1 (30), ASV_2 (37), ASV_10 (344) go.
    assert kept == ["ASV_3", "ASV_4", "ASV_5", "ASV_6", "ASV_7", "ASV_8", "ASV_9"]
    assert columns == kept
    # The log names the rule, the reference and the floor q1 - window_margin would give.
    assert ("ITS floor rule: the floor is the probe's min, 38 bp (the shortest ITS region "
            "in unite_uchime_ITS2), without window_margin") in log
    assert "q1 - window_margin would give 125 bp" in log
    assert "Filter window: [38, 343]  (source: auto:ITS_probe_min_p95+50)" in log
    assert plot.vertical_lines == [(38, "filter min (38 bp)"), (343, "filter max (343 bp)")]


def test_its1_floor_is_the_shortest_reference_region(tmp_path, plot):
    """UNITE ITS1: the window is [17, 343] (it was [124, 343])."""
    lengths = [16, 17, 60, 123, 200, 344]
    window, kept, columns, log = run_filter(tmp_path, "ITS", True, lengths, probe=PROBE_ITS1)

    assert window == {"min": 17, "max": 343, "source": "auto:ITS_probe_min_p95+50"}
    assert kept == ["ASV_2", "ASV_3", "ASV_4", "ASV_5"]
    assert columns == kept
    assert "17 bp (the shortest ITS region in unite_uchime_ITS1)" in log
    assert "q1 - window_margin would give 124 bp" in log
    assert plot.vertical_lines == [(17, "filter min (17 bp)"), (343, "filter max (343 bp)")]


@pytest.mark.parametrize("margin, ceiling", [(0, 293), (100, 393)])
def test_its_window_margin_moves_only_the_ceiling(tmp_path, margin, ceiling):
    """For ITS after ITSx, window_margin pads the ceiling only; the floor stays the min."""
    lengths = [38, 200, ceiling, ceiling + 1]
    window, kept, columns, log = run_filter(tmp_path, "ITS", True, lengths,
                                            probe=PROBE_ITS2, window_margin=margin)

    assert window == {"min": 38, "max": ceiling,
                      "source": f"auto:ITS_probe_min_p95+{margin}"}
    assert kept == ["ASV_1", "ASV_2", "ASV_3"]


# ── Every other case: the window is unchanged ─────────────────────────────────────

@pytest.mark.parametrize("amp_type, probe, extraction_enabled, lengths, expected_window", [
    # 16S with Metaxa2 extraction on: extraction alone does not bring in the ITS rule.
    ("16S", PROBE_16S, True, [323, 324, 380, 435, 436], [324, 435]),
    ("18S", PROBE_18S, False, [328, 329, 400, 544, 545], [329, 544]),
    # gyrB has a direct-mode probe, like ITS, and no extractor.
    ("gyrB", PROBE_GYRB, False, [193, 194, 240, 300, 301], [194, 300]),
])
def test_other_markers_keep_the_q1_floor(tmp_path, plot, amp_type, probe,
                                         extraction_enabled, lengths, expected_window):
    """[q1 - 50, p95 + 50], the source text and the log lines are those of v2.6.1."""
    window, kept, columns, log = run_filter(tmp_path, amp_type, extraction_enabled,
                                            lengths, probe=probe)

    assert window == {"min": expected_window[0], "max": expected_window[1],
                      "source": Q1_SOURCE}
    assert kept == ["ASV_2", "ASV_3", "ASV_4"]
    assert columns == kept
    assert "ITS floor rule" not in log
    assert (f"Probe stats: q1={probe['q1']}, p95={probe['p95']}, "
            f"median={probe['median']}, n={probe['n_amplicons']}") in log
    assert plot.vertical_lines == [
        (expected_window[0], f"filter min ({expected_window[0]} bp)"),
        (expected_window[1], f"filter max ({expected_window[1]} bp)"),
    ]


def test_its_without_extraction_keeps_the_q1_floor(tmp_path):
    """extraction.enabled: false on ITS: the filter judges whole ASVs, which still carry
    their 5.8S and LSU flanks, so the ITSx reasoning does not apply and the window stays
    [q1 - 50, p95 + 50] = [125, 343]."""
    lengths = [124, 125, 209, 343, 344]
    window, kept, columns, log = run_filter(tmp_path, "ITS", False, lengths, probe=PROBE_ITS2)

    assert window == {"min": 125, "max": 343, "source": Q1_SOURCE}
    assert kept == ["ASV_2", "ASV_3", "ASV_4"]
    assert "ITS floor rule" not in log


def test_manual_range_wins_on_an_its_run(tmp_path):
    """length_filter.range is used as given, also for ITS after ITSx."""
    lengths = [40, 100, 300, 301]
    window, kept, columns, log = run_filter(tmp_path, "ITS", True, lengths,
                                            probe=PROBE_ITS2, manual_range=[100, 300])

    assert window == {"min": 100, "max": 300, "source": "manual:config_range"}
    assert kept == ["ASV_2", "ASV_3"]
    assert "ITS floor rule" not in log


def test_fallback_without_a_probe_is_unchanged(tmp_path):
    """No probe (expected_length set by hand): expected_length -15 % / +15 % from
    trunclen.json, as before."""
    lengths = [339, 340, 400, 460, 461]
    window, kept, columns, log = run_filter(tmp_path, "16S", True, lengths,
                                            expected_length=400)

    assert window == {"min": 340, "max": 460,
                      "source": "auto:fallback_expected_length±15pct"}
    assert kept == ["ASV_2", "ASV_3", "ASV_4"]


# ── Settings that leave nothing to size the window from ────────────────────────────


def test_no_probe_and_no_expected_length_stops_with_the_settings_to_change(tmp_path):
    """mode auto, no probe, and trunclen.json's expected_length null (pick_trunclen could
    not resolve one): one ERROR line that names what to set, exit 1."""
    log = run_filter(tmp_path, "16S", True, [340, 400], expected_length=None,
                     expected_exit=1)[3]
    assert ("[dada_length_filter] ERROR: length_filter.mode is auto, but nothing can size "
            "the window") in log
    assert "length_filter.range" in log
    assert "Traceback" not in log


def test_a_mode_without_a_range_stops_with_an_error_line(tmp_path):
    log = run_filter(tmp_path, "16S", True, [340, 400], mode="manual", expected_exit=1)[3]
    assert ("[dada_length_filter] ERROR: length_filter.mode='manual' is not 'auto', and no "
            "length_filter.range is set.") in log
