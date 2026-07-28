#!/usr/bin/env python3
"""Aggregate per-stage read counts for the shotgun pipeline into a wide TSV.

Columns (left to right):
  sample → raw → [nophix] → [no_host] → trimmed → classified → unclassified

All counts are read PAIRS throughout:
  raw          : pairs entering the pipeline
                 source: BBDuk PhiX stats #Total/2 (if phix removal enabled)
                         else fastp summary.before_filtering.total_reads/2
  nophix       : pairs surviving PhiX removal (only if REMOVE_PHIX)
                 source: BBDuk PhiX stats (#Total - #Matched) / 2
  no_host      : pairs surviving host decontamination (only if HOST_GENOMES set)
                 source: BBMap dehost stats 'surviving_pairs' (kept R1 records)
  trimmed      : pairs passing fastp adapter/quality filter
                 source: fastp filtering_result.passed_filter_reads / 2
  classified   : pairs classified by Kraken2
                 source: Kraken2 report rank-R (root) read count
  unclassified : pairs not classified by Kraken2
                 source: Kraken2 report rank-U read count

Both BBDuk and fastp count R1 and R2 individually; dividing by 2 gives pairs.
Kraken2 already counts paired reads as one unit; no division needed.
"""
import json
import sys
from pathlib import Path


def log(msg: str, fh) -> None:
    print(msg, file=fh, flush=True)
    print(msg, file=sys.stderr, flush=True)


def parse_bbduk_stats(path: Path) -> tuple[int, int]:
    """Return (total_pairs, matched_pairs) from a BBDuk stats file."""
    total = matched = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith("#Total"):
                total = int(line.split()[1])
            elif line.startswith("#Matched"):
                matched = int(line.split()[1])
    return total // 2, matched // 2


def parse_host_surviving(path: Path) -> int:
    """Return surviving (non-host) pairs from the BBMap dehost stats file.
    The decontam_host rule writes a single 'surviving_pairs\\t<N>' line."""
    with path.open() as fh:
        for line in fh:
            if line.startswith("surviving_pairs"):
                return int(line.split()[1])
    return 0


def parse_fastp_counts(path: Path) -> tuple[int, int]:
    """Return (before_pairs, passed_pairs) from a fastp JSON."""
    d = json.loads(path.read_text())
    before = d["summary"]["before_filtering"]["total_reads"] // 2
    passed = d["filtering_result"]["passed_filter_reads"] // 2
    return before, passed


def parse_kraken2_report(path: Path) -> tuple[int, int]:
    """Return (classified_pairs, unclassified_pairs) from a Kraken2 report.
    Rank U = unclassified; rank R = root of classified tree.

    Column indices assume the 8-column format produced by --report-minimizer-data
    (35_kraken2.smk always passes this flag): it inserts two minimizer-count
    columns before rank/taxid/name, shifting rank from the standard report's
    index 3 to index 5 here. Without --report-minimizer-data, rank would be at
    cols[3] instead — do not "fix" this back to 3 without checking that flag."""
    classified = unclassified = 0
    with path.open() as fh:
        for line in fh:
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            rank = cols[5].strip()
            reads = int(cols[1])
            if rank == "U":
                unclassified = reads
            elif rank == "R":
                classified = reads
                break  # root row is always line 2; no need to read further
    return classified, unclassified


def _index(paths: list[str], strip_suffix: str) -> dict[str, Path]:
    """Turn a flat list of per-sample file paths into a {sample_id: path} lookup,
    recovering the sample id by stripping the stage's fixed filename suffix
    (e.g. "_fastp.json") from each file's basename. Snakemake's expand() gives
    us these paths in SAMPLES order already, but indexing by sample id makes
    the per-sample loop below (and any missing-file bug) trivial to reason
    about, instead of relying on every input list staying in lock-step order."""
    return {Path(p).name.replace(strip_suffix, ""): Path(p) for p in paths}


# ── Main ──────────────────────────────────────────────────────────────────────
# Paths, sample list, and the two decontamination toggles all come from rule
# aggregate_read_counts in 80_stats.smk. output.counts (stats/read_tracking.txt)
# is a terminal deliverable listed in _shotgun_targets() (00_common.smk) —
# requested directly by `rule all`; nothing further in the pipeline reads it.
sm = snakemake  # noqa: F821

log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log_fh = log_path.open("w")

# remove_phix / remove_host mirror REMOVE_PHIX / HOST_GENOMES from 00_common.smk:
# whether those decontamination steps ran at all for this run, which decides
# whether the "nophix" / "no_host" columns exist below.
samples     = list(sm.params.samples)
remove_phix = bool(sm.params.remove_phix)
remove_host = bool(sm.params.remove_host)

# Build one {sample_id: path} lookup per pipeline stage (see _index above), so
# the per-sample loop below can fetch each sample's file by name instead of
# assuming every input list is in the same order. The phix/host indexes are
# only built when that stage actually ran; otherwise they stay empty and are
# never consulted (guarded by remove_phix / remove_host below).
fastp_idx   = _index(list(sm.input.fastp_jsons),    "_fastp.json")
kraken_idx  = _index(list(sm.input.kraken_reports), "_report.txt")
phix_idx    = _index(list(sm.input.phix_stats),     "_dephix_stats.txt") if remove_phix else {}
host_idx    = _index(list(sm.input.host_stats),     "_dehost_stats.txt") if remove_host else {}

# Column header for the output TSV — must match the column order the per-sample
# loop below appends to `row`, one column per pipeline stage that actually ran.
stages = ["sample", "raw"]
if remove_phix:
    stages.append("nophix")
if remove_host:
    stages.append("no_host")
stages += ["trimmed", "classified", "unclassified"]

rows: list[list[str]] = []

# One row per sample, columns built in the same order as `stages` above.
for s in samples:
    # fastp's own "before" count is always read, even when PhiX removal ran:
    # it's only actually used for the "raw" column when PhiX removal did NOT
    # run (the else branch below); when PhiX removal did run, "raw" instead
    # comes from BBDuk's own total further down, since fastp then sees reads
    # only after PhiX has already been stripped out.
    fastp_before, fastp_passed = parse_fastp_counts(fastp_idx[s])

    row: list[str] = [s]

    if remove_phix:
        # "raw" = every read pair BBDuk saw; "nophix" = same minus PhiX matches.
        phix_total, phix_matched = parse_bbduk_stats(phix_idx[s])
        row.append(str(phix_total))
        row.append(str(phix_total - phix_matched))
    else:
        # No PhiX step, so fastp is the first thing to see the raw reads —
        # its "before filtering" count IS the raw count here.
        row.append(str(fastp_before))

    if remove_host:
        # Pairs left after BBMap's host-genome decontamination step.
        row.append(str(parse_host_surviving(host_idx[s])))

    # "trimmed" = pairs that passed fastp's adapter/quality filter, regardless
    # of whether PhiX/host removal ran before it.
    row.append(str(fastp_passed))

    # "classified" / "unclassified" = Kraken2's own split of the trimmed reads.
    classified, unclassified = parse_kraken2_report(kraken_idx[s])
    row.append(str(classified))
    row.append(str(unclassified))

    rows.append(row)

out_path = Path(sm.output.counts)
out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as out:
    out.write("\t".join(stages) + "\n")
    for row in rows:
        out.write("\t".join(row) + "\n")

log(f"[aggregate_read_counts] Written {out_path} for {len(samples)} sample(s)", log_fh)
log_fh.close()
