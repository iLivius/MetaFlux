#!/usr/bin/env python3
"""Estimate the expected amplicon length before any real reads are trimmed.

Why this exists: two downstream steps need a target length for this marker —
``pick_trunclen`` (does R1+R2 still overlap enough after quality trimming?)
and ``dada_length_filter`` (is this ASV a plausible amplicon or junk?). Rather
than asking the user to know that number in advance, we measure it from the
reference database itself, in silico (computationally, without a real PCR
reaction) — either by cutting the reference the same way the real primers
would, or by reading it directly when the reference is already amplicon-sized.
Which approach applies is set per marker pack (``probe_mode: pcr | direct``,
see workflow/markers/*.yaml) and passed in via ``snakemake.params.probe_mode``.

Input: ``ref_fasta`` — a full-length or pre-trimmed reference FASTA fetched/
prepared by the rules in 10_refdb.smk (SILVA for 16S, SILVA-Euk for 18S,
UNITE UCHIME's ITS1/ITS2 subregion for ITS, DD7RZ8 for gyrB, FROGS for rpoB).
``fwd``/``rev`` are this run's primer FASTA files (only used in PCR mode).

PCR mode (16S, 18S, rpoB): two-pass cutadapt against the full-length reference.
  Pass 1: ``-g file:fwd`` with ``--discard-untrimmed``
          → keep reference sequences containing the forward primer, trim it off.
  Pass 2: ``-a file:rev_rc`` with ``--discard-untrimmed``
          → of those, keep sequences containing revcomp(rev) — the reverse
          complement, i.e. the reverse primer read on the opposite strand,
          which is how it actually appears at the 3' end of a PCR product —
          and trim it off.
          The survivors are the in-silico amplicon bodies: what would be left
          of each reference sequence if it were actually PCR-amplified with
          these two primers.

Direct mode (ITS, gyrB): no primer trimming. The reference is ALREADY an
  amplicon-length sequence, whatever marker/DB it comes from — UNITE UCHIME's
  pre-extracted ITS1/ITS2 subregion for ITS, or DD7RZ8's pre-trimmed in-silico
  amplicons for gyrB. Lengths are read directly without running cutadapt, and
  the reference FASTA is gzip-copied to the amplicons output for consistency
  (so the rule's declared output always exists, regardless of which branch ran).

Both modes compute the same length distribution statistics and write a JSON:
  {n_amplicons, min, q1, median, q3, p95, p99, max, mean, stdev, probe_mode, ...}
plus amp_type/primer_hash/reference_tag, which together identify exactly which
inputs produced this result — the probe JSON is cached under refdb/cache/ under
a filename built from those same three values, so a rerun with unchanged
marker/reference/primers skips this rule entirely.

Output: the JSON above (``snakemake.output.json``) is read next by
``pick_trunclen`` (50a_pick_trunclen.py), which picks one statistic
(config amplicon.probe_length_stat) as the expected amplicon length for the
R1+R2 overlap check, and by ``dada_length_filter`` (60d_dada_length_filter.py),
which uses q1/p95 to size the ASV length-filter window. The amplicons FASTA
(``snakemake.output.amplicons``) is kept for inspection only — nothing reads
it back in.
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
max_err     = float(p.max_err)        # cutadapt -e: fraction of primer bases allowed to mismatch
pcr_min_len = int(p.pcr_min_len)      # PCR-mode floor (bp); rejects spurious near-zero-length "hits"
amp_type    = str(p.amp_type)
# primer_hash/ref_tag identify exactly which primers + reference produced this
# result. They go into the output JSON below and into the probe cache filename
# (set in 40_probe.smk), so a rerun with the same marker/reference/primers can
# reuse the cached probe instead of redoing the cutadapt passes.
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

# ── Two-pass in-silico PCR (16S, 18S, rpoB) ──────────────────────
if probe_mode == "pcr":
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        rev_rc = tmp / "rev_rc.fa"
        pass1 = tmp / "pass1.fa.gz"

        # Revcomp (reverse-complement: reverse the sequence and swap each base for its
        # pair, A<->T and C<->G) the rev primer file → rev_rc.fa. A PCR reverse primer
        # binds the template's other strand, so on the same strand as the forward
        # primer/reference it always appears as its reverse complement — that is the
        # form cutadapt needs to search for in pass 2 below.
        log.write(f"\n[amplicon_probe] revcomp rev primer: seqtk seq -r {rev} > {rev_rc}\n")
        log.flush()
        with rev_rc.open("w") as out:
            res = subprocess.run(["seqtk", "seq", "-r", rev], stdout=out, stderr=log)
        if res.returncode != 0:
            log.write(f"seqtk failed (exit {res.returncode})\n")
            log.close()
            sys.exit(res.returncode)

        # Pass 1: keep only sequences where the forward primer matches anywhere (unanchored,
        # i.e. searched at any position, not just the very first base). Reference records
        # are full-length genes (e.g. the whole 16S rRNA gene), so the primer's true binding
        # site sits somewhere in the middle of the record, never at position 1 — an anchored
        # search would never match. --discard-untrimmed drops any reference sequence the
        # primer wasn't found in at all, and trims off everything before the match.
        run([
            "cutadapt", "-j", str(threads),
            "-g", f"file:{fwd}",
            "-e", str(max_err),
            "--discard-untrimmed",
            "-o", str(pass1),
            str(ref_fasta),
        ], "pass1 (5' fwd)")

        # Pass 2: of those, keep only sequences where the reverse primer matches anywhere
        # (unanchored, same reasoning as pass 1) and trim off everything after it — what's
        # left between the two cuts is the in-silico amplicon. --minimum-length drops any
        # survivor shorter than pcr_min_len, which only happens if the two primers matched
        # implausibly close together (a spurious hit, not a real amplicon).
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

# ── Direct length measurement (ITS, gyrB): reference is already amplicon-length ─
# (ITS: UNITE UCHIME pre-extracted subregion; gyrB: DD7RZ8 pre-trimmed amplicons)
else:
    log.write(
        f"\n[amplicon_probe] direct mode: reading lengths from {ref_fasta} "
        "(reference is already amplicon-length, no cutadapt)\n"
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
