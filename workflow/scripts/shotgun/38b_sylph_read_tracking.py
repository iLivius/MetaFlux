#!/usr/bin/env python3
"""Per-stage read counts for a profiler: sylph run (rule sylph_read_tracking,
38_sylph.smk), written to stats/read_tracking.txt.

Columns (left to right):
  sample → raw → [nophix] → [no_host] → trimmed

These are the cleaning columns of the read-tracking table the Kraken2 rule
aggregate_read_counts writes (90a_aggregate_read_counts.py), from the same files and
with the same arithmetic. That rule also adds Kraken2's classified/unclassified
columns; under profiler: sylph Kraken2 does not run, its rule file is not read, and
this script writes the table instead. All counts are read PAIRS:
  raw      : pairs entering the pipeline
             source: BBDuk PhiX stats #Total/2 (if PhiX removal ran)
                     else fastp summary.before_filtering.total_reads/2
  nophix   : pairs surviving PhiX removal (only if REMOVE_PHIX)
             source: BBDuk PhiX stats (#Total - #Matched) / 2
  no_host  : pairs surviving host removal (only if HOST_GENOMES is set)
             source: BBMap dehost stats 'surviving_pairs' (kept R1 records)
  trimmed  : pairs passing fastp's adapter/quality filter
             source: fastp filtering_result.passed_filter_reads / 2
             (the column total of 03.abundance/sylph_table.tsv)
"""
import json
import sys
from pathlib import Path


def log(msg: str, fh) -> None:
    print(msg, file=fh, flush=True)
    print(msg, file=sys.stderr, flush=True)


# ── Parsers ───────────────────────────────────────────────────────────────────
# The four functions below are copied word for word from
# 90a_aggregate_read_counts.py, which must stay untouched so that finished Kraken2
# runs are not re-run (Snakemake re-runs a script rule whose script has changed).
# Keep the two copies in step: a change there must be made here too. The smoke tests
# check that the two scripts give the same cleaning columns on the same files.

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


def _index(paths: list[str], strip_suffix: str) -> dict[str, Path]:
    """Turn a flat list of per-sample file paths into a {sample_id: path} lookup,
    recovering the sample id by stripping the stage's fixed filename suffix
    (e.g. "_fastp.json") from each file's basename. Snakemake's expand() gives
    us these paths in SAMPLES order already, but indexing by sample id makes
    the per-sample loop below (and any missing-file bug) trivial to reason
    about, instead of relying on every input list staying in lock-step order."""
    return {Path(p).name.replace(strip_suffix, ""): Path(p) for p in paths}


# ── Main ──────────────────────────────────────────────────────────────────────
# Paths, the sample list and the two decontamination switches come from rule
# sylph_read_tracking (38_sylph.smk). The output, stats/read_tracking.txt, is a final
# file listed in _shotgun_targets() (00_common.smk); nothing in the pipeline reads it.
sm = snakemake  # noqa: F821

log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log_fh = log_path.open("w")

# remove_phix / remove_host mirror REMOVE_PHIX / HOST_GENOMES from 00_common.smk:
# whether those cleaning steps ran at all, which decides whether the "nophix" /
# "no_host" columns exist.
samples     = list(sm.params.samples)
remove_phix = bool(sm.params.remove_phix)
remove_host = bool(sm.params.remove_host)

# One {sample: file} lookup per cleaning stage. The PhiX and host lookups are built
# only when that stage ran; otherwise they stay empty and are never used.
fastp_idx = _index(list(sm.input.fastp_jsons), "_fastp.json")
if remove_phix:
    phix_idx = _index(list(sm.input.phix_stats), "_dephix_stats.txt")
else:
    phix_idx = {}
if remove_host:
    host_idx = _index(list(sm.input.host_stats), "_dehost_stats.txt")
else:
    host_idx = {}

# Column header: one column per cleaning stage that actually ran, in the same order
# the loop below fills each row.
stages = ["sample", "raw"]
if remove_phix:
    stages.append("nophix")
if remove_host:
    stages.append("no_host")
stages.append("trimmed")

rows: list[list[str]] = []

for s in samples:
    # fastp's "before" count is the raw count only when PhiX removal did NOT run;
    # otherwise fastp saw the reads after PhiX was removed, and "raw" comes from
    # BBDuk's own total instead (the same rule as 90a_aggregate_read_counts.py).
    fastp_before, fastp_passed = parse_fastp_counts(fastp_idx[s])

    row: list[str] = [s]

    if remove_phix:
        # "raw" = every read pair BBDuk saw; "nophix" = the same minus PhiX matches.
        phix_total, phix_matched = parse_bbduk_stats(phix_idx[s])
        row.append(str(phix_total))
        row.append(str(phix_total - phix_matched))
    else:
        # No PhiX step: fastp was the first to see the raw reads.
        row.append(str(fastp_before))

    if remove_host:
        # Pairs left after BBMap's host removal.
        row.append(str(parse_host_surviving(host_idx[s])))

    # "trimmed" = pairs that passed fastp, whichever steps ran before it. These are
    # the read pairs sylph sketched, and the column total of sylph_table.tsv.
    row.append(str(fastp_passed))

    rows.append(row)

out_path = Path(sm.output.counts)
out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as out:
    out.write("\t".join(stages) + "\n")
    for row in rows:
        out.write("\t".join(row) + "\n")

log(f"[sylph_read_tracking] Written {out_path} for {len(samples)} sample(s)", log_fh)
log_fh.close()
