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
             source: BBDuk PhiX log 'Input:' reads/2 (if PhiX removal ran)
                     else fastp summary.before_filtering.total_reads/2
  nophix   : pairs surviving PhiX removal (only if REMOVE_PHIX)
             source: BBDuk PhiX log 'Result:' reads/2 (the reads BBDuk kept)
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
# 90a_aggregate_read_counts.py, the script that writes this same table when Kraken2
# runs (rule aggregate_read_counts, 80_stats.smk). They are copied, not imported,
# because 90a builds its whole table as soon as it is loaded: importing it here would
# run that too.
# Keep the two copies word for word in step: a change to a parser there must be made
# here too, so that profiler: sylph and profiler: kraken2 write the same cleaning
# columns from the same files. The smoke tests check that the two scripts give the
# same cleaning columns on the same files.

def parse_bbduk_log(path: Path) -> tuple[int, int]:
    """Return (pairs_in, pairs_kept) from the log of rule decontam_phix.

    BBDuk ends its log with a short summary that counts READS (both mates):
        Input:          11836374 reads     1730990881 bases.
        ...
        Result:         11836370 reads (100.00%)     1730990288 bases (100.00%)
    "Input" is what BBDuk read, "Result" is what it wrote out. BBDuk removes
    the whole pair when either mate matches PhiX, so both counts are even and
    half of each is a number of read pairs.

    The stats file (stats=) is not used for this: its #Matched counts the reads
    that matched PhiX, so a pair with PhiX in one mate only would count as half
    a pair, and halving #Matched under-counts the pairs BBDuk removed."""
    reads_in = None
    reads_kept = None
    with path.open() as fh:
        for line in fh:
            if line.startswith("Input:"):
                reads_in = int(line.split()[1])
            elif line.startswith("Result:"):
                reads_kept = int(line.split()[1])
    # Both lines are printed only when BBDuk finished; a log without them
    # should stop the run here rather than put a made-up number in the table.
    if reads_in is None or reads_kept is None:
        raise ValueError(f"no 'Input:' and 'Result:' lines in BBDuk log {path}")
    return reads_in // 2, reads_kept // 2


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
    # removesuffix() takes the suffix off the end of the file name only.
    # replace() would also cut it out of the middle of a sample name: the PhiX
    # log of a sample called "Plot.log2" is "Plot.log2.log", which
    # replace(".log", "") turns into "Plot2", so the loop below could not find it.
    return {Path(p).name.removesuffix(strip_suffix): Path(p) for p in paths}


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
    phix_idx = _index(list(sm.input.phix_logs), ".log")
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
        # "raw" = every read pair BBDuk read; "nophix" = the pairs it kept. So
        # raw - nophix counts every pair BBDuk removed: pairs with PhiX in either
        # mate, plus any pair with a read shorter than 10 bp, which BBDuk drops
        # by default (rare in raw reads).
        phix_in, phix_kept = parse_bbduk_log(phix_idx[s])
        row.append(str(phix_in))
        row.append(str(phix_kept))
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
