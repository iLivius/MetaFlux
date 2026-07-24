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
    return {Path(p).name.replace(strip_suffix, ""): Path(p) for p in paths}


# ── Main ──────────────────────────────────────────────────────────────────────
sm = snakemake  # noqa: F821

log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log_fh = log_path.open("w")

samples     = list(sm.params.samples)
remove_phix = bool(sm.params.remove_phix)
remove_host = bool(sm.params.remove_host)

fastp_idx   = _index(list(sm.input.fastp_jsons),    "_fastp.json")
kraken_idx  = _index(list(sm.input.kraken_reports), "_report.txt")
phix_idx    = _index(list(sm.input.phix_stats),     "_dephix_stats.txt") if remove_phix else {}
host_idx    = _index(list(sm.input.host_stats),     "_dehost_stats.txt") if remove_host else {}

stages = ["sample", "raw"]
if remove_phix:
    stages.append("nophix")
if remove_host:
    stages.append("no_host")
stages += ["trimmed", "classified", "unclassified"]

rows: list[list[str]] = []

for s in samples:
    fastp_before, fastp_passed = parse_fastp_counts(fastp_idx[s])

    row: list[str] = [s]

    if remove_phix:
        phix_total, phix_matched = parse_bbduk_stats(phix_idx[s])
        row.append(str(phix_total))
        row.append(str(phix_total - phix_matched))
    else:
        row.append(str(fastp_before))

    if remove_host:
        row.append(str(parse_host_surviving(host_idx[s])))

    row.append(str(fastp_passed))

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
