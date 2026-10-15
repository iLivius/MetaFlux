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

PCR mode (16S, 18S, rpoB): two-pass cutadapt against the full-length reference
(three passes for rpoB).
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
  Pass 3 (only when the marker pack sets ``probe_exclude_second_site: true``,
          today rpoB): ``-g`` with each forward primer, required at its own full
          length, with ``--discard-trimmed``
          → keep only the amplicon bodies that do NOT hold a second full-length
          forward-primer site. The others are left out of the length statistics;
          the amplicons FASTA still contains them.
          Why: a primer can match a reference twice, at its real site and at a
          weaker look-alike stretch. With a degenerate primer cutadapt keeps the
          LEFTMOST of the two matches, even when the later one fits better (here
          the later one is often a perfect match). For the reverse primer (pass
          2) that is the site nearest the forward primer, the right one, so no
          body holds a second reverse site. For the forward primer it is the
          site farther from the reverse primer, so the product starts at the
          wrong site and comes out too long.
          Univ_rpoB_deg_F also matches, with the 4 mismatches that max_error_rate
          0.2 allows, a stretch about 920 bp upstream of its real site in 7,289 of
          the 42,995 FROGS genes the two primers amplify (17%, mostly
          Gammaproteobacteria and Actinomycetes). Those products, most of them
          1.0-2.0 kb long, pulled p95 to 1,332 bp; without them p95 is 395 bp.
          Which of the two sites a PCR would use cannot be told from the sequence
          alone (the inner site is usually the better match, not always), so such
          a product is left out rather than cut back to the inner site.
          16S and 18S leave this pass off: the same thing happens there, but
          rarely (0.5 to 0.7% of the products), and their probe results are kept
          as they are.

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
marker/reference/primers skips this rule entirely. When pass 3 runs, the JSON
also holds n_excluded_second_site, and the reference tag in the filename ends
in "_excl-second-site" (set in 00_common.smk).

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
# Reference lineages left out of the length statistics (marker pack field
# probe_exclude_taxa; empty for most markers). For 16S this is Chloroplast,
# Mitochondria and Eukaryota: organelle and eukaryote SSU amplified in silico by
# bacterial primers are often much longer (insertions), and SILVA 144 carries
# about eight times more of them than 138.2 — enough to move p95 from 386 to 397 bp
# for 799F/1175R, and with it truncLen, although these are exactly the sequences
# the taxonomy filter removes later. The amplicons FASTA still contains them.
exclude_taxa = {str(t) for t in (p.exclude_taxa or [])}
# Whether pass 3 runs: leave out of the statistics every product that holds a second
# full-length forward-primer site (marker pack field probe_exclude_second_site; true
# for rpoB only, see the module docstring). Read in 00_common.smk, which also stops
# the run when the field is set on a direct-mode marker.
exclude_second_site = bool(p.exclude_second_site)

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


# ── Read header + length of every sequence in a FASTA file ─────
def iter_fasta_records(path: Path):
    """Yield (header, sequence length) for each record of a (gzipped) FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header, seq_len = None, 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None and seq_len:
                    yield header, seq_len
                header, seq_len = line[1:].rstrip("\n"), 0
            else:
                seq_len += len(line.strip())
        if header is not None and seq_len:
            yield header, seq_len


def is_excluded(header: str) -> bool:
    """True if any ';'-separated rank of the header is one of exclude_taxa.

    The SILVA DADA2 headers ARE the lineage (">Bacteria;...;Mitochondria;"), and
    cutadapt keeps them on the amplicons. A SILVA 144 "--other" suffix is ignored,
    as in the taxonomy filter.
    """
    for rank in header.split(";"):
        rank = rank.strip()
        if rank.endswith("--other"):
            rank = rank[: -len("--other")]
        if rank in exclude_taxa:
            return True
    return False


def lengths_for_stats(path: Path) -> tuple[list[int], int]:
    """Lengths of the records that count for the statistics, and how many were left out."""
    kept, n_excluded = [], 0
    for header, length in iter_fasta_records(path):
        if exclude_taxa and is_excluded(header):
            n_excluded += 1
        else:
            kept.append(length)
    return kept, n_excluded


# ── Read the primer sequences of a primer FASTA file ────────────
def read_primer_sequences(path: Path) -> list[str]:
    """Every sequence of a primer FASTA file, upper case, in file order.

    Used by pass 3 only, which hands each forward primer to cutadapt on its own so
    that each one can be required at its own full length. A primer file may hold
    several primers (a primer mix); each of them counts as a primer site.
    """
    sequences: list[str] = []
    current: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                current = []
            else:
                current.append(line)
    if current:
        sequences.append("".join(current).upper())
    return sequences


# ── Count the records of a FASTA file ───────────────────────────
def count_fasta_records(path: Path) -> int:
    """Number of (non-empty) records in a (gzipped) FASTA file."""
    n = 0
    for _header, _length in iter_fasta_records(path):
        n += 1
    return n


# ── Compute a percentile from a sorted list of lengths ──────────
def pct(sorted_lengths: list[int], p: float) -> int:
    if not sorted_lengths:
        return 0
    idx = min(int(p / 100 * len(sorted_lengths)), len(sorted_lengths) - 1)
    return sorted_lengths[idx]


out_amplicons.parent.mkdir(parents=True, exist_ok=True)

# Filled by pass 3 only (pcr mode, packs with probe_exclude_second_site: true): how
# many in-silico products there were, and how many of them pass 3 left out.
n_products = 0
n_second_site = 0

# ── Two-pass in-silico PCR (16S, 18S, rpoB), plus pass 3 for rpoB ─
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

        # Pass 3 (packs with probe_exclude_second_site: true, i.e. rpoB): find the
        # amplicon bodies that still hold a whole forward-primer site, and keep only
        # the others for the statistics. Why, with the rpoB numbers: module docstring.
        #   -g          : every forward primer as a 5' search, as in pass 1, each
        #                 passed on its own so that each gets its own min_overlap.
        #   min_overlap : the primer's own length, so only a whole primer site counts.
        #                 With cutadapt's default of 3 bases, a body that merely starts
        #                 with the last bases of the forward primer would count too.
        #   -e          : the same mismatch rate as pass 1, so a second site is one
        #                 that pass 1 would have taken as a primer site.
        #   --discard-trimmed : write only the bodies in which no primer was found.
        # The reverse primer is not searched again: pass 2 already keeps the reverse
        # site nearest the forward primer, so a body never holds a second one.
        # The kept bodies go to a temporary file; out_amplicons keeps every body.
        stats_fasta = out_amplicons
        if exclude_second_site:
            single_site = tmp / "single_site.fa.gz"
            pass3_cmd = ["cutadapt", "-j", str(threads)]
            for primer in read_primer_sequences(Path(fwd)):
                pass3_cmd += ["-g", f"{primer};min_overlap={len(primer)}"]
            pass3_cmd += [
                "-e", str(max_err),
                "--discard-trimmed",
                "-o", str(single_site),
                str(out_amplicons),
            ]
            run(pass3_cmd, "pass3 (second forward-primer site)")
            n_products = count_fasta_records(out_amplicons)
            n_second_site = n_products - count_fasta_records(single_site)
            stats_fasta = single_site

        # Read the lengths here, while the temporary folder (and pass 3's file in it)
        # still exists.
        lengths, n_excluded = lengths_for_stats(stats_fasta)

    if not lengths:
        if n_products > 0 and n_second_site == n_products:
            # Pass 3 left nothing: every product holds a second forward-primer site.
            log.write(
                f"\n[amplicon_probe] ERROR: each of the {n_products} in-silico "
                "products holds a second full-length forward-primer site, so none is "
                "left to measure the amplicon length from. Set amplicon.expected_length "
                "to a number, and length_filter.mode: manual with a range, for these "
                "primers.\n"
            )
        else:
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

    lengths, n_excluded = lengths_for_stats(ref_fasta)

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
# Which reference lineages were left out of the statistics above (none for most
# markers), and how many amplicons that removed.
result["excluded_taxa"] = sorted(exclude_taxa)
result["n_excluded"] = n_excluded
# How many products pass 3 left out because they hold a second forward-primer site.
# Written only when pass 3 ran (rpoB), so the JSON of every other marker stays as it was.
if exclude_second_site:
    result["n_excluded_second_site"] = n_second_site

# ── Write results and close log ──────────────────────────────────
out_json.parent.mkdir(parents=True, exist_ok=True)
out_json.write_text(json.dumps(result, indent=2) + "\n")

if exclude_taxa:
    log.write(f"\n[amplicon_probe] left out of the statistics: {n_excluded} sequences "
              f"from {', '.join(sorted(exclude_taxa))}\n")
if exclude_second_site:
    log.write(f"\n[amplicon_probe] left out of the statistics: {n_second_site} of "
              f"{n_products} products, which hold a second full-length "
              f"forward-primer site\n")
log.write(
    f"\n[amplicon_probe] DONE. n={len(lengths)} sequences; "
    f"median={result['median']} p95={result['p95']} q1={result['q1']} q3={result['q3']}\n"
)
log.close()
