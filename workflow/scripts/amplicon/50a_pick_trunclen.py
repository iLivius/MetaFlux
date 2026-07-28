#!/usr/bin/env python3
"""Auto-pick truncLen from falco's per-base sequence quality table.

What truncLen does, and the trade-off
-------------------------------------
``truncLen`` is one length per read direction; any bases past it are discarded.
That single number is pulled two opposite ways:

    - Trim shorter -> cuts off the low-quality 3' tail of the read, which is
      what DADA2's error model wants (cleaner bases).
    - Trim too short -> R1 and R2 stop overlapping in the middle, so DADA2 can
      no longer merge them into one amplicon, and a pair that won't merge is
      dropped entirely.

So the cut should remove the bad tail but stop while the forward and reverse
reads still reach each other.

How the cut is chosen
---------------------
Falco generates a FastQC-compatible ``fastqc_data.txt``. For each direction we
read the ``>>Per base sequence quality`` section, take the "Lower Quartile"
(Q1, 25th percentile) at each position, aggregate across samples (median per
bin), and cut at the first position where the aggregated Q1 drops below
``q_threshold`` — i.e. where quality falls off.

Bin notation: positions 1..9 are reported as single integers, then ``10-14``,
``15-19``, ... in fixed-width bins. ``truncLen`` is set to ``(bin_start - 1)``
so we keep cycles up to (but not including) the first sub-threshold bin.

Keeping the pair mergeable
--------------------------
The two cuts must still leave enough overlap to merge:

    truncR1 + truncR2 >= amplicon_length + min_overlap

If the quality-based cuts break this, ``resolve_policy`` decides what to do:
    - raise_trunc: extend the cuts past the quality drop to recover the missing
      bases, splitting the deficit between R1/R2 by how much room each still has
      (its headroom); error if even the combined headroom is not enough.
    - relax_q    : lower the Q threshold by 1 (down to ``q_floor``) and re-cut,
      until the overlap holds.
    - error      : abort and report the deficit.

ITS amplicons
-------------
The picker still runs on ITS — it computes the Q-based cuts and writes them to
``trunclen.json`` — but ``dada_filter`` overrides ``truncLen`` to ``c(0, 0)`` for
ITS: a fixed-length cut would discard genuinely short fungal amplicons, since ITS
length varies. The 3' low-quality tail is instead trimmed adaptively there by
``truncQ``, with ``maxEE``/``minLen`` as the quality gate — so if ITS R2 quality
is poor, the lever is a stricter ``truncQ``/``maxEE``, not a fixed ``truncLen``.
The cuts written here are kept for reference/QC only and do not affect ITS
filtering (this also means ``manual_r1``/``manual_r2`` have no effect for ITS).

Manual mode bypasses Q analysis and uses ``manual_r1``/``manual_r2`` from config.
"""
import json
import math
import sys
from pathlib import Path
from statistics import median


# ── Read the quality table from a falco report ───────────────────
def parse_qbins(fastqc_data_path: Path) -> list[tuple[int, int, float]]:
    """Return [(bin_start, bin_end, q1), ...] from a falco fastqc_data.txt."""
    bins: list[tuple[int, int, float]] = []
    in_block = False
    with fastqc_data_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">>Per base sequence quality"):
                in_block = True
                continue
            if not in_block:
                continue
            if line.startswith(">>END_MODULE"):
                break
            if line.startswith("#"):
                continue  # header row
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            base_field = parts[0]
            q1 = float(parts[3])  # "Lower Quartile"
            if "-" in base_field:
                lo, hi = base_field.split("-", 1)
                bin_start, bin_end = int(lo), int(hi)
            else:
                bin_start = bin_end = int(base_field)
            bins.append((bin_start, bin_end, q1))
    if not bins:
        raise ValueError(f"No 'Per base sequence quality' bins found in {fastqc_data_path}")
    return bins


# ── Read how long the reads in one sample actually are ───────────
def parse_length_distribution(fastqc_data_path: Path) -> list[tuple[int, float]]:
    """Return [(length, n_reads), ...] from a falco fastqc_data.txt.

    Falco's ">>Sequence Length Distribution" block says how many reads ended at
    each length. We need it because the quality table alone is misleading: falco
    reports a quality bin for every position up to the LONGEST read in the file,
    even when only a handful of reads get that far. After primer trimming most
    reads stop a little short of the full read length, so trusting the quality
    table's last position means trusting a position backed by a few reads.

    The length field is usually a single integer, but falco bins it ("270-274")
    once the spread is wide, exactly like the quality table. A binned entry is
    counted at its LOWEST length: we only want to claim a read reached a position
    if it certainly did.
    """
    counts: list[tuple[int, float]] = []
    in_block = False
    with fastqc_data_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">>Sequence Length Distribution"):
                in_block = True
                continue
            if not in_block:
                continue
            if line.startswith(">>END_MODULE"):
                break
            if line.startswith("#"):
                continue  # header row
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            length_field = parts[0]
            n_reads = float(parts[1])
            if "-" in length_field:
                lo, _hi = length_field.split("-", 1)
                length = int(lo)          # conservative: shortest length in the bin
            else:
                length = int(length_field)
            counts.append((length, n_reads))
    return counts


# ── How far do most of this sample's reads actually reach? ────────
def read_coverage_cap(length_counts: list[tuple[int, float]], min_pct: float) -> int | None:
    """Longest position at least `min_pct` percent of the sample's reads still reach.

    Input is parse_length_distribution's output for one sample. Walking from the
    longest length downwards and accumulating reads, we stop as soon as the running
    total covers the requested percentage — that length is the cap.

    Why this matters: DADA2's filterAndTrim DISCARDS any read shorter than truncLen
    rather than padding it. So a truncLen above this cap does not merely trim fewer
    bases, it throws away every read that fell short. On a library whose quality
    never degrades, the quality-based cut has nothing to stop it landing there.

    Returns None if the file had no length data, in which case the caller leaves
    the sample uncapped and behaves as before.
    """
    if not length_counts:
        return None
    total = sum(n for _length, n in length_counts)
    if total <= 0:
        return None
    needed = total * min_pct / 100.0
    running = 0.0
    for length, n in sorted(length_counts, reverse=True):
        running += n
        if running >= needed:
            return length
    # Reached here only through float rounding; the shortest length covers everything.
    return min(length for length, _n in length_counts)


# ── Drop quality bins past a sample's usable read length ─────────
def cap_bins_at(bins: list[tuple[int, int, float]], cap: int | None) -> list[tuple[int, int, float]]:
    """Trim a sample's quality bins so none extends past `cap`.

    A bin straddling the cap is shortened to end on it; bins entirely beyond it are
    dropped. `cap=None` means "no length data for this sample" and returns the bins
    untouched.
    """
    if cap is None:
        return bins
    out: list[tuple[int, int, float]] = []
    for bin_start, bin_end, q1 in bins:
        if bin_start > cap:
            continue
        out.append((bin_start, min(bin_end, cap), q1))
    return out


# ── Convert binned quality scores to one value per cycle position ─
def expand_bins_to_positions(bins: list[tuple[int, int, float]]) -> dict[int, float]:
    """Expand binned Q1 to a per-position dict: {pos: q1_of_containing_bin}."""
    out: dict[int, float] = {}
    for bin_start, bin_end, q1 in bins:
        for pos in range(bin_start, bin_end + 1):
            out[pos] = q1
    return out


# ── Summarise quality across all samples at each cycle position ──
def aggregate_q1_across_samples(per_sample_bins: list[list[tuple[int, int, float]]],
                                 min_coverage_pct: float = 100.0
                                 ) -> list[tuple[int, int, float]]:
    """Aggregate Q1 per position across samples by median.

    Bin layouts can differ at the tail (variable read lengths across samples).
    Expand each sample's bins to per-position Q1, then take the median per
    position across samples that have data there. Re-emit as (pos, pos, q1)
    tuples (one-position "bins") — keeps the downstream first_drop_trunclen
    logic the same.

    Coverage gate: only keep positions covered by at least `min_coverage_pct`
    of SAMPLES. Default 100 means every sample must reach the position, so the
    ceiling is the shortest sample's usable read length.

    Note what this gate does NOT do: it counts samples, not reads. Lowering it
    admits positions that some samples never reach, which RAISES the ceiling —
    the opposite of sacrificing short-read outliers. Reads that are too short
    within a sample are handled before this function is called, by capping each
    sample's bins with read_coverage_cap(); that is the gate that protects
    against truncLen landing above the bulk of the data.
    """
    if not per_sample_bins:
        raise ValueError("No samples to aggregate.")
    per_sample_pos = [expand_bins_to_positions(b) for b in per_sample_bins]
    n_samples = len(per_sample_pos)
    # 100% gate → need data from ALL samples; lower gates round up to nearest sample.
    # math.ceil, not int(): int() truncates DOWN (e.g. 95% of 21 samples = 19.95 →
    # int() gives 19, but "at least 95%" needs ceil() = 20) — only matters once a
    # caller passes a non-default min_coverage_pct, which none does today.
    min_n = max(1, math.ceil(min_coverage_pct * n_samples / 100))
    if min_coverage_pct >= 100:
        min_n = n_samples
    all_positions = sorted({pos for s in per_sample_pos for pos in s})
    out: list[tuple[int, int, float]] = []
    for pos in all_positions:
        q1s = [s[pos] for s in per_sample_pos if pos in s]
        if len(q1s) >= min_n:
            out.append((pos, pos, median(q1s)))
    return out


# ── Find the position where quality first falls below the threshold ─
def first_drop_trunclen(bins: list[tuple[int, int, float]], threshold: float) -> int:
    """truncLen = last cycle still above threshold.

    If Q1 < threshold first occurs in bin (start, end), truncate at start - 1.
    If Q1 never drops below threshold, truncate at the last bin's end.
    """
    for bin_start, _bin_end, q1 in bins:
        if q1 < threshold:
            return bin_start - 1
    return bins[-1][1]


# ── Convert the config amplicon length setting to a single number ─
def resolve_expected_length(expected_length, amp_type, probe_json_path, log_fn,
                            probe_stat: str = "p95") -> tuple[int, str]:
    """Validate and resolve amplicon.expected_length from config to a single int.

    Accepts:
      - int / float                       e.g. 338
      - [min, max] list                   e.g. [320, 360]  (uses max for overlap constraint)
      - "auto"                            reads the amplicon_probe JSON and picks
                                          the stat resolved per-marker upstream (probe_stat).
    Returns (resolved_int, source_label) for inclusion in trunclen.json.
    """
    if expected_length is None:
        raise ValueError(
            "amplicon.expected_length is empty/null in config. Set an integer (e.g. 338), "
            "a [min, max] range, or 'auto' to invoke the amplicon_probe rule."
        )
    if isinstance(expected_length, bool):
        raise ValueError(
            f"amplicon.expected_length must be int / [min,max] / 'auto', got bool: {expected_length}"
        )
    if isinstance(expected_length, str):
        if expected_length.strip().lower() == "auto":
            if probe_json_path is None or not Path(probe_json_path).exists():
                raise FileNotFoundError(
                    f"amplicon.expected_length='auto' but probe JSON is missing: {probe_json_path}. "
                    f"Snakemake should have produced this via the amplicon_probe rule; "
                    f"check the pipeline DAG and refdb/cache/."
                )
            probe = json.loads(Path(probe_json_path).read_text())
            # Which statistic to read is already resolved per-marker upstream:
            # 00_common.smk maps the marker pack's probe_stat_key onto
            # amplicon.probe_length_stat.<key> and passes it in as probe_stat.
            # ITS is the only marker that returns before this point (see the
            # amp_type == "ITS" branch in main(), above) — every other marker
            # reaches here, pcr-mode (16S/18S/rpoB) or direct-mode (gyrB) alike.
            # 00_common.smk's AMPLICON_TYPE != "ITS" guard on probe_length_stat
            # is what guarantees probe_stat is a valid key by the time we get here.
            stat = probe_stat
            if stat not in probe:
                raise KeyError(
                    f"Probe JSON does not contain stat '{stat}'. "
                    f"Available keys: {list(probe.keys())}"
                )
            resolved = int(probe[stat])
            log_fn(
                f"[pick_trunclen] auto-resolved expected_length = {resolved} bp "
                f"(via {stat} of probe distribution: n={probe['n_amplicons']}, "
                f"min={probe['min']}, q1={probe['q1']}, median={probe['median']}, "
                f"q3={probe['q3']}, p95={probe['p95']}, p99={probe['p99']}, max={probe['max']})"
            )
            return resolved, f"auto:{stat}"
        raise ValueError(
            f"amplicon.expected_length string must be 'auto' (got: {expected_length!r}). "
            f"Use an integer or [min, max] list otherwise."
        )
    if isinstance(expected_length, list):
        if len(expected_length) != 2:
            raise ValueError(
                f"amplicon.expected_length list must have exactly 2 elements [min, max]; "
                f"got {len(expected_length)} element(s): {expected_length}"
            )
        for v in expected_length:
            if not isinstance(v, (int, float)) or isinstance(v, bool):
                raise ValueError(
                    f"amplicon.expected_length list elements must be numeric; got {v!r}"
                )
        if any(v <= 0 for v in expected_length):
            raise ValueError(
                f"amplicon.expected_length values must be positive; got {expected_length}"
            )
        if expected_length[0] > expected_length[1]:
            raise ValueError(
                f"amplicon.expected_length [min, max] must satisfy min <= max; got {expected_length}"
            )
        return int(max(expected_length)), "manual:range_max"
    if isinstance(expected_length, (int, float)):
        v = int(expected_length)
        if v <= 0:
            raise ValueError(
                f"amplicon.expected_length must be positive; got {expected_length}"
            )
        if v < 50 or v > 10_000:
            raise ValueError(
                f"amplicon.expected_length={v} is outside the plausible amplicon range "
                f"(50-10000 bp). Check the value in config.amplicon.expected_length."
            )
        return v, "manual:int"
    raise ValueError(
        f"amplicon.expected_length must be int / [min,max] / 'auto', got "
        f"{type(expected_length).__name__}: {expected_length!r}"
    )


# ── Save a result dictionary to a JSON file ──────────────────────
def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")


def main() -> int:
    # ── Inputs, parameters, and log file ────────────────────────
    sm = snakemake  # noqa: F821 (injected by Snakemake)

    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = log_path.open("w")
    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)
        print(msg, file=sys.stderr, flush=True)

    falco_r1 = [Path(p) for p in sm.input.falco_r1]
    falco_r2 = [Path(p) for p in sm.input.falco_r2]
    out_json = Path(sm.output.json)

    p = sm.params
    amp_type        = p.amplicon_type
    expected_length = p.expected_length
    min_overlap     = int(p.min_overlap)
    mode            = p.mode
    q_threshold     = float(p.q_threshold)
    q_floor         = float(p.q_floor)
    resolve_policy  = p.resolve_policy
    manual_r1       = p.manual_r1
    manual_r2       = p.manual_r2
    probe_stat  = str(p.probe_stat)
    min_read_coverage_pct = float(p.min_read_coverage_pct)

    log(f"[pick_trunclen] amp_type={amp_type} mode={mode} q_threshold={q_threshold}")
    log(f"[pick_trunclen] R1 falco files: {len(falco_r1)}; R2 falco files: {len(falco_r2)}")

    # ── Manual mode: use the lengths provided directly in config ─
    if mode == "manual":
        if manual_r1 is None or manual_r2 is None:
            log("[pick_trunclen] ERROR: trunc_len.mode=manual but manual_r1/manual_r2 are not set")
            return 2
        log(f"[pick_trunclen] Manual mode: r1={manual_r1} r2={manual_r2}")
        write_json(out_json, {
            "r1": int(manual_r1), "r2": int(manual_r2),
            "mode": "manual", "amplicon_type": amp_type,
            "resolved_via": "manual",
        })
        return 0

    # ── Auto mode: parse quality files and find the primary cuts ─
    bins_r1 = [parse_qbins(p) for p in falco_r1]
    bins_r2 = [parse_qbins(p) for p in falco_r2]

    # Cap each sample at the length most of its reads still reach, BEFORE looking
    # at quality. Falco reports a quality bin for every position up to the longest
    # read in a file, so without this step the ceiling is set by however few reads
    # happened to run to full length — and since filterAndTrim discards reads
    # shorter than truncLen, a cut placed there loses nearly the whole library.
    if not 0 < min_read_coverage_pct <= 100:
        log(f"[pick_trunclen] ERROR: trunc_len.min_read_coverage_pct must be >0 and <=100, "
            f"got {min_read_coverage_pct}")
        return 2

    caps_r1 = [read_coverage_cap(parse_length_distribution(p), min_read_coverage_pct) for p in falco_r1]
    caps_r2 = [read_coverage_cap(parse_length_distribution(p), min_read_coverage_pct) for p in falco_r2]
    bins_r1 = [cap_bins_at(b, c) for b, c in zip(bins_r1, caps_r1)]
    bins_r2 = [cap_bins_at(b, c) for b, c in zip(bins_r2, caps_r2)]

    # A sample whose bins were emptied has no usable positions left, which would break
    # the aggregation below. In practice this needs reads shorter than the first quality
    # bin, so it means the input is not what the run assumes — say the wrong files, or
    # reads already trimmed to almost nothing. Fail here with a readable reason rather
    # than further down with an index error.
    for label, caps, capped, paths in (("R1", caps_r1, bins_r1, falco_r1),
                                       ("R2", caps_r2, bins_r2, falco_r2)):
        for cap, s, path in zip(caps, capped, paths):
            if not s:
                log(f"[pick_trunclen] ERROR: {label} sample {path.parent.name} has no quality "
                    f"positions left after capping at {cap} bp "
                    f"({min_read_coverage_pct}% read-coverage). Check the primer-trimmed reads.")
                return 2

    # Report which sample sets each ceiling, and by how much the read-coverage cap
    # pulled it in — a large gap here is the signature of the failure described above.
    def cap_report(caps, bins, falco_paths, label):
        rows = []
        for cap, s, path in zip(caps, bins, falco_paths):
            if not s:
                continue
            raw_max = max(b[1] for b in parse_qbins(path))
            rows.append((cap if cap is not None else raw_max, raw_max, path.parent.name))
        if not rows:
            return
        rows.sort(key=lambda x: x[0])
        capped, raw_max, name = rows[0]
        log(f"[pick_trunclen] {label} ceiling {capped} bp, set by {name} "
            f"(its longest read is {raw_max} bp; {min_read_coverage_pct}% of its reads reach {capped} bp)")
        if raw_max - capped >= 10:
            log(f"[pick_trunclen] NOTE: {label} read-coverage cap pulled the ceiling in by "
                f"{raw_max - capped} bp. Without it the cut could have landed above most reads.")

    cap_report(caps_r1, bins_r1, falco_r1, "R1")
    cap_report(caps_r2, bins_r2, falco_r2, "R2")

    agg_r1 = aggregate_q1_across_samples(bins_r1)
    agg_r2 = aggregate_q1_across_samples(bins_r2)

    primary_r1 = first_drop_trunclen(agg_r1, q_threshold)
    primary_r2 = first_drop_trunclen(agg_r2, q_threshold)
    max_r1 = agg_r1[-1][1]
    max_r2 = agg_r2[-1][1]
    log(f"[pick_trunclen] Primary cuts at Q1>={q_threshold}: R1={primary_r1}, R2={primary_r2}")
    log(f"[pick_trunclen] Ceilings after the read-coverage and all-samples gates: "
        f"R1={max_r1}, R2={max_r2}")

    if primary_r1 < 50 or primary_r2 < 50:
        log(f"[pick_trunclen] WARNING: very short primary cuts (R1={primary_r1}, R2={primary_r2}); inspect falco reports.")

    # ── ITS: overlap constraint not applicable — truncLen is overridden to c(0,0) ─
    # Skip resolve_expected_length entirely (not used for ITS) and write the
    # primary cuts so dada_filter can read them before overriding to (0, 0).
    if amp_type == "ITS":
        log("[pick_trunclen] ITS mode: skipping overlap constraint (truncLen will be overridden to c(0,0) in dada_filter)")
        write_json(out_json, {
            "r1": int(primary_r1),
            "r2": int(primary_r2),
            "mode": "auto",
            "amplicon_type": amp_type,
            "q_threshold_requested": q_threshold,
            "q_threshold_used": q_threshold,
            "expected_length": None,
            "expected_length_source": "not_applicable:ITS_truncLen_overridden",
            "min_overlap": min_overlap,
            "overlap_slack": None,
            "resolved_via": "primary",
            "aggregation": "median_across_samples_per_bin",
        })
        log(f"[pick_trunclen] DONE. truncLen=(R1={primary_r1}, R2={primary_r2}) — will be overridden to c(0,0) for ITS")
        return 0

    # ── Non-ITS (16S/18S/gyrB/rpoB): resolve the expected amplicon length and check overlap ──
    probe_json_path = None
    if hasattr(sm.input, "probe_json") and sm.input.probe_json:
        # Snakemake gives a list (or string) here depending on how it was declared.
        pj = sm.input.probe_json
        probe_json_path = pj[0] if isinstance(pj, list) else pj

    exp_len, exp_len_source = resolve_expected_length(
        expected_length, amp_type, probe_json_path, log,
        probe_stat=probe_stat,
    )
    log(f"[pick_trunclen] expected_length (for overlap): {exp_len} (source: {exp_len_source}), min_overlap: {min_overlap}")

    min_total = exp_len + min_overlap
    slack = (primary_r1 + primary_r2) - min_total
    log(f"[pick_trunclen] Overlap: r1+r2={primary_r1 + primary_r2}, need>={min_total}, slack={slack}")

    # ── Auto mode: enforce the overlap constraint if violated ────
    cut_r1, cut_r2 = primary_r1, primary_r2
    q_used = q_threshold
    resolved_via = "primary"

    if slack < 0:
        log(f"[pick_trunclen] Constraint violated by {-slack} bp; resolve_policy={resolve_policy}")
        if resolve_policy == "raise_trunc":
            headroom_r1 = max_r1 - cut_r1
            headroom_r2 = max_r2 - cut_r2
            deficit = -slack
            if headroom_r1 + headroom_r2 < deficit:
                log(f"[pick_trunclen] ERROR: raise_trunc cannot satisfy overlap "
                    f"(deficit={deficit}, headroom R1={headroom_r1}+R2={headroom_r2})")
                return 3

            # Which read gives up the deficit is decided by the QUALITY of the bases
            # about to be added, not by how many are left.
            #
            # Splitting by headroom alone (ceiling - cut) is actively perverse: headroom
            # is largest exactly where the cut landed early, i.e. where quality collapsed
            # soonest. On typical Illumina that is R2, so the read with the worst tail
            # would be extended furthest — pushing the recovery deepest into the bases
            # least worth keeping, which then cost reads at maxEE anyway.
            #
            # Instead, walk the deficit one base at a time and each time extend whichever
            # read has the better aggregated Q1 at the position it would gain. Ties go to
            # R1, which is the stronger read on Illumina chemistry. A read stops being a
            # candidate once its own headroom is used up.
            q_by_pos_r1 = {pos: q for pos, _end, q in agg_r1}
            q_by_pos_r2 = {pos: q for pos, _end, q in agg_r2}

            bump_r1 = bump_r2 = 0
            gained_q_r1: list[float] = []
            gained_q_r2: list[float] = []
            for _ in range(deficit):
                can_r1 = bump_r1 < headroom_r1
                can_r2 = bump_r2 < headroom_r2
                if not can_r1 and not can_r2:
                    break
                # Quality of the next base each read would gain; missing positions sort last.
                next_q_r1 = q_by_pos_r1.get(cut_r1 + bump_r1 + 1, float("-inf")) if can_r1 else float("-inf")
                next_q_r2 = q_by_pos_r2.get(cut_r2 + bump_r2 + 1, float("-inf")) if can_r2 else float("-inf")
                if next_q_r1 >= next_q_r2:
                    gained_q_r1.append(next_q_r1)
                    bump_r1 += 1
                else:
                    gained_q_r2.append(next_q_r2)
                    bump_r2 += 1

            cut_r1 += bump_r1
            cut_r2 += bump_r2
            resolved_via = "raise_trunc"

            def _mean_q(vals: list[float]) -> str:
                finite = [v for v in vals if v != float("-inf")]
                return f"{sum(finite) / len(finite):.1f}" if finite else "n/a"

            log(f"[pick_trunclen] Raised: R1 +{bump_r1} (mean Q1 of added bases "
                f"{_mean_q(gained_q_r1)}), R2 +{bump_r2} (mean Q1 {_mean_q(gained_q_r2)}); "
                "allocation is quality-weighted, not headroom-weighted")

        elif resolve_policy == "relax_q":
            q_try = q_threshold
            ok = False
            while q_try > q_floor:
                q_try -= 1
                c1 = first_drop_trunclen(agg_r1, q_try)
                c2 = first_drop_trunclen(agg_r2, q_try)
                if c1 + c2 >= min_total:
                    cut_r1, cut_r2 = c1, c2
                    q_used = q_try
                    resolved_via = "relax_q"
                    log(f"[pick_trunclen] Relaxed to Q1>={q_try}: R1={c1}, R2={c2}")
                    ok = True
                    break
            if not ok:
                log(f"[pick_trunclen] ERROR: relax_q reached floor Q={q_floor} without satisfying overlap")
                return 4

        elif resolve_policy == "error":
            log(f"[pick_trunclen] ERROR: overlap constraint violated (deficit={-slack} bp at Q>={q_threshold})")
            return 5

        else:
            log(f"[pick_trunclen] ERROR: unknown resolve_policy={resolve_policy!r}")
            return 6

    # ── Write the final truncLen decision to JSON ────────────────
    final_slack = (cut_r1 + cut_r2) - min_total
    result = {
        "r1": int(cut_r1),
        "r2": int(cut_r2),
        "mode": "auto",
        "amplicon_type": amp_type,
        "q_threshold_requested": q_threshold,
        "q_threshold_used": q_used,
        "expected_length": exp_len,
        "expected_length_source": exp_len_source,
        "min_overlap": min_overlap,
        "overlap_slack": final_slack,
        "resolved_via": resolved_via,
        "aggregation": "median_across_samples_per_bin",
    }
    write_json(out_json, result)
    log(f"[pick_trunclen] DONE. truncLen=(R1={cut_r1}, R2={cut_r2}) via {resolved_via}; final slack={final_slack}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
