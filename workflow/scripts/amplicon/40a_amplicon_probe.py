#!/usr/bin/env python3
"""Amplicon length probe — two modes depending on amplicon type.

PCR mode (16S): two-pass cutadapt against SILVA.
  Pass 1: ``-g file:fwd`` with ``--discard-untrimmed``
          → keep reference sequences with the forward primer at 5', trim it off.
  Pass 2: ``-a file:rev_rc`` with ``--discard-untrimmed``
          → of those, keep sequences with revcomp(rev) at 3', trim it off.
          The survivors are the in-silico amplicon bodies.

Direct mode (ITS): no primer trimming.
  The reference FASTA is already the pre-extracted ITS1 or ITS2 subregion
  from UNITE UCHIME; lengths are read directly without running cutadapt.
  The reference FASTA is gzip-copied to the amplicons output for consistency.

Both modes compute the same length distribution statistics and write a JSON:
  {n_amplicons, min, q1, median, q3, p95, p99, max, mean, stdev, probe_mode, ...}

Downstream ``pick_trunclen`` reads the probe_length_stat from config to select
which statistic to use as expected_length.
"""
import gzip
import json
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

# ── Inputs, parameters, and log file ────────────────────────────
sm = snakemake  # noqa: F821 (injected by Snakemake)

ref_fasta = Path(sm.input.ref_fasta)
fwd = str(sm.input.fwd)
rev = str(sm.input.rev)
out_json = Path(sm.output.json)
out_amplicons = Path(sm.output.amplicons)

p = sm.params
probe_mode  = str(p.probe_mode)
max_err     = float(p.max_err)
pcr_min_len = int(p.pcr_min_len)
amp_type    = str(p.amp_type)
primer_hash = str(p.primer_hash)
ref_tag     = str(p.ref_tag)

threads = int(sm.threads)
log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log = log_path.open("w")


# ── Run an external tool and log its output ──────────────────────
def run(cmd: list[str], step: str) -> None:
    log.write(f"\n[amplicon_probe] {step}: {' '.join(cmd)}\n")
    log.flush()
    res = subprocess.run(cmd, stdout=log, stderr=log)
    if res.returncode != 0:
        log.write(f"[amplicon_probe] {step} failed (exit {res.returncode})\n")
        log.close()
        sys.exit(res.returncode)


# ── Read the lengths of all sequences in a FASTA file ───────────
def iter_fasta_lengths(path: Path):
    """Yield sequence lengths from a (gzipped) FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    seq_len = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_len:
                    yield seq_len
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if seq_len:
            yield seq_len


# ── Compute a percentile from a sorted list of lengths ──────────
def pct(sorted_lengths: list[int], p: float) -> int:
    if not sorted_lengths:
        return 0
    idx = min(int(p / 100 * len(sorted_lengths)), len(sorted_lengths) - 1)
    return sorted_lengths[idx]


out_amplicons.parent.mkdir(parents=True, exist_ok=True)

# ── Two-pass in-silico PCR (16S) ─────────────────────────────────
if probe_mode == "pcr":
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        rev_rc = tmp / "rev_rc.fa"
        pass1 = tmp / "pass1.fa.gz"

        # Revcomp the rev primer file → rev_rc.fa
        log.write(f"\n[amplicon_probe] revcomp rev primer: seqtk seq -r {rev} > {rev_rc}\n")
        log.flush()
        with rev_rc.open("w") as out:
            res = subprocess.run(["seqtk", "seq", "-r", rev], stdout=out, stderr=log)
        if res.returncode != 0:
            log.write(f"seqtk failed (exit {res.returncode})\n")
            log.close()
            sys.exit(res.returncode)

        # Pass 1: keep only sequences that start with the forward primer
        run([
            "cutadapt", "-j", str(threads),
            "-g", f"file:{fwd}",
            "-e", str(max_err),
            "--discard-untrimmed",
            "-o", str(pass1),
            str(ref_fasta),
        ], "pass1 (5' fwd)")

        # Pass 2: of those, keep only sequences that end with the reverse primer
        run([
            "cutadapt", "-j", str(threads),
            "-a", f"file:{rev_rc}",
            "-e", str(max_err),
            "--discard-untrimmed",
            "--minimum-length", str(pcr_min_len),
            "-o", str(out_amplicons),
            str(pass1),
        ], "pass2 (3' rev_rc)")

    lengths = list(iter_fasta_lengths(out_amplicons))

    if not lengths:
        log.write(
            "\n[amplicon_probe] ERROR: no reference sequences had both primers. "
            "Check primer sequences against the reference DB orientation and "
            "primer error tolerance.\n"
        )
        log.close()
        sys.exit(2)

# ── Direct length measurement from pre-extracted subregion (ITS) ─
else:
    log.write(
        f"\n[amplicon_probe] direct mode: reading lengths from {ref_fasta} "
        "(UNITE UCHIME pre-extracted ITS subregion, no cutadapt)\n"
    )
    log.flush()

    lengths = list(iter_fasta_lengths(ref_fasta))

    if not lengths:
        log.write(
            "\n[amplicon_probe] ERROR: no sequences found in reference FASTA. "
            f"Check that {ref_fasta} exists and is not empty.\n"
        )
        log.close()
        sys.exit(2)

    # Copy the reference FASTA to the amplicons output (gzipped) so the
    # output file declared in the rule always exists regardless of mode.
    log.write(f"[amplicon_probe] copying {ref_fasta} → {out_amplicons}\n")
    log.flush()
    if str(ref_fasta).endswith(".gz"):
        shutil.copy2(ref_fasta, out_amplicons)
    else:
        with open(ref_fasta, "rb") as src, gzip.open(out_amplicons, "wb") as dst:
            shutil.copyfileobj(src, dst)

# ── Length distribution statistics ──────────────────────────────
sorted_lengths = sorted(lengths)
result = {
    "amplicon_type": amp_type,
    "probe_mode": probe_mode,
    "primer_hash": primer_hash,
    "reference_tag": ref_tag,
    "reference_fasta": str(ref_fasta),
    "n_amplicons": len(lengths),
    "min": min(lengths),
    "q1": pct(sorted_lengths, 25),
    "median": pct(sorted_lengths, 50),
    "q3": pct(sorted_lengths, 75),
    "p95": pct(sorted_lengths, 95),
    "p99": pct(sorted_lengths, 99),
    "max": max(lengths),
    "mean": round(statistics.mean(lengths), 2),
    "stdev": round(statistics.stdev(lengths), 2) if len(lengths) > 1 else 0.0,
}
if probe_mode == "pcr":
    result["cutadapt_max_error_rate"] = max_err
    result["cutadapt_min_length"] = pcr_min_len

# ── Write results and close log ──────────────────────────────────
out_json.parent.mkdir(parents=True, exist_ok=True)
out_json.write_text(json.dumps(result, indent=2) + "\n")

log.write(
    f"\n[amplicon_probe] DONE. n={len(lengths)} sequences; "
    f"median={result['median']} p95={result['p95']} q1={result['q1']} q3={result['q3']}\n"
)
log.close()
