#!/usr/bin/env python3
"""ASV length distribution statistics and length-based filter.

Runs AFTER target_extract (Metaxa2 / ITSx) when extraction is enabled, so the
filter sees post-extraction lengths that are directly comparable to the probe.
When extraction is disabled, runs directly on the post-DADA2 seqs.fasta.

The script always also reads the pre-extraction FASTA (seqs.fasta) so the
stats JSON can show both distributions side-by-side — useful for evaluating
how much Metaxa2/ITSx trimmed each ASV.

Window source priority:
  1. manual_range in config (explicit [min, max])
  2. auto + probe JSON available → [probe_q1 - margin, probe_p95 + margin]
  3. auto + no probe JSON → [expected_length * 0.85, expected_length * 1.15]
"""
import json
import math
import sys
from collections import Counter
from pathlib import Path
from statistics import mean, median, stdev

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── FASTA / seqtab helpers ────────────────────────────────────────────────────

def parse_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current_id = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            elif current_id is not None:
                chunks.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(chunks)
    return seqs


def write_fasta(seqs: dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for asv_id, seq in seqs.items():
            fh.write(f">{asv_id}\n{seq}\n")


def filter_seqtab_by_ids(in_path: Path, kept_ids: set[str], out_path: Path,
                         log_fn) -> tuple[list[str], list[str], list[list[str]]]:
    """Subset an ID-keyed seqtab to columns whose header is in kept_ids.

    Returns (sample_names, kept_asv_ids_in_order, per_sample_count_rows) for
    use by the seq-keyed regeneration step.
    """
    with in_path.open() as fh:
        raw_header = fh.readline().rstrip("\n")
        data_lines = fh.readlines()
    header = raw_header.split("\t")
    # R write.table wraps every column name in double quotes; strip them for matching.
    header_clean = [h.strip('"') for h in header]
    keep_idx = [0] + [i for i, h in enumerate(header_clean[1:], 1) if h in kept_ids]

    kept_asv_ids = [header_clean[i] for i in keep_idx[1:]]
    sample_names: list[str] = []
    count_rows: list[list[str]] = []

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("\t".join(header[i] for i in keep_idx) + "\n")
        for line in data_lines:
            cols = line.rstrip("\n").split("\t")
            fh.write("\t".join(cols[i] for i in keep_idx) + "\n")
            sample_names.append(cols[0])
            count_rows.append([cols[i] for i in keep_idx[1:]])
    log_fn(f"[dada_length_filter] {out_path.name}: kept {len(kept_asv_ids)} / {len(header)-1} columns")
    return sample_names, kept_asv_ids, count_rows


def write_seq_keyed_seqtab(out_path: Path, sample_names: list[str],
                           kept_asv_ids: list[str], count_rows: list[list[str]],
                           asv_to_seq: dict[str, str]) -> None:
    """Rebuild a sequence-keyed seqtab from the ID-keyed counts and the FASTA mapping."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        header_cells = [""] + [f'"{asv_to_seq[a]}"' for a in kept_asv_ids]
        fh.write("\t".join(header_cells) + "\n")
        for sample, row in zip(sample_names, count_rows):
            fh.write("\t".join([sample] + row) + "\n")


# ── Statistics ────────────────────────────────────────────────────────────────

def pct(sorted_vals: list[int], p: float) -> int:
    if not sorted_vals:
        return 0
    idx = min(int(p / 100 * len(sorted_vals)), len(sorted_vals) - 1)
    return sorted_vals[idx]


def compute_stats(lengths: list[int]) -> dict:
    if not lengths:
        return {"n": 0}
    s = sorted(lengths)
    return {
        "n": len(s),
        "min": s[0], "max": s[-1],
        "q1": pct(s, 25), "median": pct(s, 50), "q3": pct(s, 75),
        "p95": pct(s, 95), "p99": pct(s, 99),
        "mean": round(mean(lengths), 2),
        "stdev": round(stdev(lengths), 2) if len(lengths) > 1 else 0.0,
    }


# ── Window determination ──────────────────────────────────────────────────────

def determine_window(mode, manual_range, probe_json_path, trunclen_json_path,
                     window_margin, log_fn) -> tuple[int, int, str]:
    margin = int(window_margin)

    if manual_range is not None and len(manual_range) == 2:
        min_len, max_len = int(manual_range[0]), int(manual_range[1])
        return min_len, max_len, "manual:config_range"

    if mode == "auto":
        pj = None
        if probe_json_path:
            pj_path = Path(probe_json_path) if not isinstance(probe_json_path, list) else (
                Path(probe_json_path[0]) if probe_json_path else None
            )
            if pj_path and pj_path.exists():
                pj = json.loads(pj_path.read_text())

        if pj is not None:
            min_len = int(pj["q1"]) - margin
            max_len = int(pj["p95"]) + margin
            source = f"auto:probe_q1_p95+/-{margin}"
            log_fn(f"[dada_length_filter] Probe stats: q1={pj['q1']}, p95={pj['p95']}, "
                   f"median={pj['median']}, n={pj['n_amplicons']}")
            return min_len, max_len, source

        # Fallback: expected_length from trunclen.json ± 15%
        tl = json.loads(Path(trunclen_json_path).read_text())
        exp = int(tl["expected_length"])
        min_len = math.floor(exp * 0.85)
        max_len = math.ceil(exp * 1.15)
        return min_len, max_len, "auto:fallback_expected_length±15pct"

    raise ValueError(
        f"length_filter.mode={mode!r} is not 'auto', and no length_filter.range is set. "
        "Set mode to 'auto' or 'manual', and if 'manual', a [min, max] range."
    )


# ── Histogram ─────────────────────────────────────────────────────────────────

def write_histogram(lengths_pre: list[int], lengths_post: list[int],
                    min_len: int, max_len: int, extraction_enabled: bool,
                    out_png: Path, log_fn) -> None:
    """Plot the filter-source distribution; overlay pre-extraction when extraction is on."""
    fig, ax = plt.subplots(figsize=(12, 5))
    all_vals = lengths_pre + lengths_post
    lo, hi = min(all_vals), max(all_vals)
    bin_width = max(1, (hi - lo) // 60)
    bins = list(range(lo - bin_width, hi + 2 * bin_width, bin_width))

    kept = [l for l in lengths_post if min_len <= l <= max_len]

    if extraction_enabled:
        ax.hist(lengths_pre, bins=bins, color="lightgray", edgecolor="white", linewidth=0.4,
                label=f"pre-extraction (n={len(lengths_pre)})")
    ax.hist(lengths_post, bins=bins, color="steelblue", edgecolor="white", linewidth=0.4,
            label=f"filter source (n={len(lengths_post)})")
    ax.hist(kept, bins=bins, color="seagreen", edgecolor="white", linewidth=0.4, alpha=0.7,
            label=f"kept (n={len(kept)})")
    ax.axvline(min_len, color="red", linestyle="--", linewidth=1.5,
               label=f"filter min ({min_len} bp)")
    ax.axvline(max_len, color="red", linestyle=":",  linewidth=1.5,
               label=f"filter max ({max_len} bp)")
    ax.set_xlabel("ASV length (bp)")
    ax.set_ylabel("Number of ASVs")
    title = "ASV length distribution"
    if extraction_enabled:
        title += " — pre-extraction vs filter source"
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    log_fn(f"[dada_length_filter] Histogram saved: {out_png}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    sm = snakemake  # noqa: F821

    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = log_path.open("w")
    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)
        print(msg, file=sys.stderr, flush=True)

    seqs_pre_in   = Path(sm.input.seqs_pre)
    seqs_in       = Path(sm.input.seqs)
    seqtab_names  = Path(sm.input.seqtab_names)
    trunclen_json = sm.input.trunclen_json

    seqs_out          = Path(sm.output.seqs)
    seqtab_names_out  = Path(sm.output.seqtab_names)
    seqtab_seqs_out   = Path(sm.output.seqtab_seqs)
    stats_json_out    = Path(sm.output.stats_json)
    hist_png_out      = Path(sm.output.hist_png)

    p = sm.params
    mode               = p.mode
    window_margin      = p.window_margin
    manual_range       = p.manual_range  # None or [min, max]
    amp_type           = p.amp_type
    extraction_enabled = bool(p.extraction_enabled)

    # probe_json: empty list when WANTS_PROBE=False, single-element list otherwise
    probe_json_path = None
    if hasattr(sm.input, "probe_json") and sm.input.probe_json:
        pj = sm.input.probe_json
        probe_json_path = pj[0] if isinstance(pj, list) else pj

    log(f"[dada_length_filter] amp_type={amp_type}  mode={mode}  "
        f"window_margin={window_margin}  extraction_enabled={extraction_enabled}")

    # ── Length distributions: pre-extraction and filter source ────────────────
    fasta_pre = parse_fasta(seqs_pre_in)
    pre_lens = [len(s) for s in fasta_pre.values()]
    stats_pre = compute_stats(pre_lens)
    log(f"[dada_length_filter] Pre-extraction ({len(fasta_pre)} ASVs from {seqs_pre_in.name}): "
        f"min={stats_pre['min']}, q1={stats_pre['q1']}, median={stats_pre['median']}, "
        f"q3={stats_pre['q3']}, p95={stats_pre['p95']}, max={stats_pre['max']}, "
        f"mean={stats_pre['mean']}, stdev={stats_pre['stdev']}")

    fasta_post = parse_fasta(seqs_in)
    post_lens = [len(s) for s in fasta_post.values()]
    stats_post = compute_stats(post_lens)
    log(f"[dada_length_filter] Filter source ({len(fasta_post)} ASVs from {seqs_in.name}): "
        f"min={stats_post['min']}, q1={stats_post['q1']}, median={stats_post['median']}, "
        f"q3={stats_post['q3']}, p95={stats_post['p95']}, max={stats_post['max']}, "
        f"mean={stats_post['mean']}, stdev={stats_post['stdev']}")

    if extraction_enabled:
        n_pre = len(fasta_pre)
        n_post = len(fasta_post)
        log(f"[dada_length_filter] Extraction kept {n_post} / {n_pre} ASVs "
            f"(dropped: {n_pre - n_post} ASVs ITSx/Metaxa2 could not extract)")

    # ── Window determination ──────────────────────────────────────────────────
    min_len, max_len, source = determine_window(
        mode, manual_range, probe_json_path, trunclen_json, window_margin, log
    )
    log(f"[dada_length_filter] Filter window: [{min_len}, {max_len}]  (source: {source})")

    # ── Apply the filter to the filter-source sequences ───────────────────────
    kept = {asv_id: seq for asv_id, seq in fasta_post.items()
            if min_len <= len(seq) <= max_len}
    dropped = len(fasta_post) - len(kept)
    log(f"[dada_length_filter] Kept {len(kept)} / {len(fasta_post)} ASVs  ({dropped} dropped)")

    if len(kept) == 0:
        log("[dada_length_filter] ERROR: all ASVs filtered out — widen the window or check the probe")
        log_fh.close()
        return 1

    # ── Stats JSON: kept distribution + side-by-side pre/post ─────────────────
    kept_lens = [len(seq) for seq in kept.values()]
    stats_kept = compute_stats(kept_lens)
    result: dict = {
        "amp_type":  amp_type,
        "extraction_enabled": extraction_enabled,
        "filter_window": {"min": min_len, "max": max_len, "source": source},
        "n_dropped_by_length_filter": dropped,
        "pre_extraction_asvs": {
            **stats_pre,
            "length_counts": dict(sorted(Counter(pre_lens).items())),
        },
        "filter_source_asvs": {
            **stats_post,
            "length_counts": dict(sorted(Counter(post_lens).items())),
        },
        "kept_asvs": {
            **stats_kept,
            "length_counts": dict(sorted(Counter(kept_lens).items())),
        },
    }
    stats_json_out.parent.mkdir(parents=True, exist_ok=True)
    stats_json_out.write_text(json.dumps(result, indent=2) + "\n")

    # ── Write outputs ─────────────────────────────────────────────────────────
    write_fasta(kept, seqs_out)
    write_histogram(pre_lens, post_lens, min_len, max_len,
                    extraction_enabled, hist_png_out, log)

    kept_ids = set(kept.keys())
    sample_names, kept_asv_ids, count_rows = filter_seqtab_by_ids(
        seqtab_names, kept_ids, seqtab_names_out, log
    )
    write_seq_keyed_seqtab(seqtab_seqs_out, sample_names, kept_asv_ids,
                           count_rows, kept)
    log(f"[dada_length_filter] Wrote seq-keyed seqtab → {seqtab_seqs_out.name}")

    log(f"[dada_length_filter] DONE. {len(kept)} ASVs retained in [{min_len}, {max_len}] bp window.")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
