"""Manual trunc_len mode and the length filter's expected_length fallback.

What is tested:
  - workflow/scripts/amplicon/50a_pick_trunclen.py in manual mode records the
    resolved expected_length (and expected_length_source) in trunclen.json for
    integer, [min, max] range, and 'auto' probe inputs, leaves it null for ITS
    or when expected_length is null, and still validates non-null values early.
  - workflow/scripts/amplicon/60d_dada_length_filter.py in auto mode without a
    probe JSON (manual expected_length) reads expected_length from that
    manual-mode trunclen.json and builds the expected +/- 15% length window,
    instead of crashing with KeyError after DADA2.
  - When trunclen.json has no expected_length (null or absent), auto length
    filtering without a probe JSON raises a ValueError explaining what to set.
"""
import json
import sys
import types
from pathlib import Path

import pytest

from helpers import SCRIPTS, fake_snakemake, run_script, write_text

PICK_SCRIPT = SCRIPTS / "amplicon" / "50a_pick_trunclen.py"
FILTER_SCRIPT = SCRIPTS / "amplicon" / "60d_dada_length_filter.py"


class StandInAxes:
    """Minimal matplotlib Axes stand-in so 60d_dada_length_filter.py runs
    without matplotlib installed, matching tests/test_kmer_evidence.py."""

    def hist(self, *args, **kwargs):
        pass

    def axvline(self, *args, **kwargs):
        pass

    def set_xlabel(self, *args, **kwargs):
        pass

    def set_ylabel(self, *args, **kwargs):
        pass

    def set_title(self, *args, **kwargs):
        pass

    def legend(self, *args, **kwargs):
        pass


class StandInFigure:
    """Minimal matplotlib Figure stand-in that writes a placeholder PNG."""

    def tight_layout(self):
        pass

    def savefig(self, path, **kwargs):
        Path(path).write_text("stand-in histogram\n")


@pytest.fixture
def standin_matplotlib(monkeypatch):
    """Install a lightweight matplotlib stand-in for scripts that draw plots."""
    pyplot = types.ModuleType("matplotlib.pyplot")
    pyplot.subplots = lambda **kwargs: (StandInFigure(), StandInAxes())
    pyplot.close = lambda fig: None
    matplotlib = types.ModuleType("matplotlib")
    matplotlib.use = lambda backend: None
    matplotlib.pyplot = pyplot
    monkeypatch.setitem(sys.modules, "matplotlib", matplotlib)
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", pyplot)


def run_manual_pick(tmp_path, expected_length, amp_type="rpoB", probe_json=None):
    """Run 50a_pick_trunclen.py in manual mode and return the parsed trunclen.json."""
    out_json = tmp_path / "stats" / "trunclen.json"
    sm = fake_snakemake(
        input={
            "falco_r1": [],
            "falco_r2": [],
            "probe_json": [str(probe_json)] if probe_json is not None else [],
        },
        output={"json": str(out_json)},
        params={
            "amplicon_type": amp_type,
            "expected_length": expected_length,
            "min_overlap": 12,
            "mode": "manual",
            "q_threshold": 20,
            "q_floor": 15,
            "resolve_policy": "raise_trunc",
            "manual_r1": 250,
            "manual_r2": 200,
            "probe_stat": "p95",
            "min_read_coverage_pct": 95,
        },
        log=[str(tmp_path / "logs" / "pick_trunclen.log")],
    )
    with pytest.raises(SystemExit) as exited:
        run_script(PICK_SCRIPT, sm)
    assert exited.value.code == 0
    return out_json, json.loads(out_json.read_text())


def run_length_filter_auto(tmp_path, trunclen_json, amp_type="rpoB"):
    """Run 60d_dada_length_filter.py in auto mode without a probe JSON."""
    seqs_path = write_text(
        tmp_path / "05.dada2" / "seqs.fasta",
        f">ASV_1\n{'A' * 387}\n>ASV_2\n{'C' * 100}\n",
    )
    seqtab_path = write_text(
        tmp_path / "05.dada2" / "seqtab_head_names.txt",
        '""\t"ASV_1"\t"ASV_2"\n"S1"\t100\t10\n',
    )
    stats_out = tmp_path / "stats" / "dada2" / "asv_length_stats.json"
    sm = fake_snakemake(
        input={
            "seqs_pre": seqs_path,
            "seqs": seqs_path,
            "seqtab_names": seqtab_path,
            "trunclen_json": str(trunclen_json),
            "probe_json": [],
        },
        output={
            "seqs": str(tmp_path / "05.dada2" / "seqs_lenfilt.fasta"),
            "seqtab_names": str(tmp_path / "05.dada2" / "seqtab_lenfilt_head_names.txt"),
            "seqtab_seqs": str(tmp_path / "05.dada2" / "seqtab_lenfilt_head_seqs.txt"),
            "stats_json": str(stats_out),
            "hist_png": str(tmp_path / "stats" / "dada2" / "asv_length_hist.png"),
        },
        params={
            "mode": "auto",
            "window_margin": 50,
            "manual_range": None,
            "amp_type": amp_type,
            "extraction_enabled": False,
        },
        log=[str(tmp_path / "logs" / "dada_length_filter.log")],
    )
    with pytest.raises(SystemExit) as exited:
        run_script(FILTER_SCRIPT, sm)
    assert exited.value.code == 0
    return json.loads(stats_out.read_text())


def test_manual_trunclen_feeds_auto_length_filter_fallback(tmp_path, standin_matplotlib):
    # Manual trunc_len + manual expected_length + auto length_filter (no probe JSON):
    # pick_trunclen must record expected_length so dada_length_filter can compute
    # its +/- 15% fallback window instead of failing with KeyError.
    out_json, trunclen = run_manual_pick(tmp_path, expected_length=387)
    assert trunclen["r1"] == 250
    assert trunclen["r2"] == 200
    assert trunclen["expected_length"] == 387
    assert trunclen["expected_length_source"] == "manual:int"

    stats = run_length_filter_auto(tmp_path, out_json)
    assert stats["filter_window"] == {
        "min": 328,
        "max": 446,
        "source": "auto:fallback_expected_length±15pct",
    }
    assert stats["kept_asvs"]["n"] == 1
    assert stats["n_dropped_by_length_filter"] == 1


def test_manual_trunclen_resolves_range_and_auto_probe(tmp_path):
    _, range_tl = run_manual_pick(tmp_path / "range", expected_length=[115, 150], amp_type="18S")
    assert range_tl["expected_length"] == 150
    assert range_tl["expected_length_source"] == "manual:range_max"

    probe_path = tmp_path / "probe.json"
    probe_path.write_text(json.dumps({
        "n_amplicons": 100, "min": 350, "q1": 380, "median": 395,
        "q3": 405, "p95": 420, "p99": 430, "max": 440,
    }))
    _, auto_tl = run_manual_pick(tmp_path / "auto", expected_length="auto", amp_type="16S",
                                 probe_json=probe_path)
    assert auto_tl["expected_length"] == 420
    assert auto_tl["expected_length_source"] == "auto:p95"


def test_manual_trunclen_allows_null_expected_length_and_its(tmp_path):
    _, null_tl = run_manual_pick(tmp_path / "null", expected_length=None, amp_type="rpoB")
    assert null_tl["expected_length"] is None
    assert null_tl["expected_length_source"] is None

    _, its_tl = run_manual_pick(tmp_path / "its", expected_length=300, amp_type="ITS")
    assert its_tl["expected_length"] is None
    assert its_tl["expected_length_source"] == "not_applicable:ITS_truncLen_overridden"


def test_length_filter_auto_fallback_rejects_missing_expected_length(tmp_path, standin_matplotlib):
    trunclen_json = write_text(
        tmp_path / "stats" / "trunclen.json",
        json.dumps({"r1": 250, "r2": 200, "mode": "manual"}),
    )
    with pytest.raises(ValueError, match="expected_length"):
        run_length_filter_auto(tmp_path, trunclen_json)
