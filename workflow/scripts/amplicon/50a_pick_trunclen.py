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

How the cut is chosen (method: quality, the default)
----------------------------------------------------
Falco generates a FastQC-compatible ``fastqc_data.txt``. For each direction we
read the ``>>Per base sequence quality`` section, take the "Lower Quartile"
(Q1, 25th percentile) at each position, take the median across samples at each
position, and cut at the first position where that median drops below
``q_threshold`` — i.e. where quality falls off.

Position groups: Falco reports positions 1..9 one by one, then groups of positions
(``10-14``, ``15-19``, ... on 2 x 250 and 2 x 300 reads), each with the mean of its
per-position values. Each group is expanded back to one value per position, so
``truncLen`` is the last position before the first group below the threshold.

Before quality is read, each sample is capped at the length that
``min_read_coverage_pct`` of its reads reach (the read-length limit). Reads shorter
than half the expected amplicon length are left out of that count: they come from a
short off-target product, not from the target (see offtarget_floor below).

What the cut is based on (logged): the kind of quality values the reads carry
(continuous on a MiSeq, a few "bins" on most newer instruments; see
scan_quality_values below), what ``q_threshold`` means on them, and whether each cut
came from quality or from the read-length limit.

Keeping the pair mergeable
--------------------------
The two cuts must still leave enough overlap to merge:

    truncR1 + truncR2 >= amplicon_length + min_overlap

If the quality-based cuts break this, ``resolve_policy`` decides what to do:
    - raise_trunc: extend the cuts past the quality drop to recover the missing
      bases, one base at a time, each time on the read whose next base has the
      better median Q1 (ties go to R1); error if even both read-length limits
      together are not enough.
    - relax_q    : lower the Q threshold by 1 (down to ``q_floor``) and re-cut,
      until the overlap holds.
    - error      : abort and report the deficit.

The other method (method: simulate)
-----------------------------------
Instead of reading a quality summary, simulate filterAndTrim on a sample of real
read pairs (the idea behind FIGARO): for each read, find the longest truncLen at
which it would still pass the run's truncQ, maxN (0) and maxEE, then choose the
R1/R2 pair of cuts that keeps the most pairs while still leaving
expected length + min_overlap + simulate_extra_overlap bases (of the cut pairs within
SIMULATE_TOLERANCE_PCT of the most, the longest). It reads the same
quality values the filter will, so binning does not matter to it, and it sees what
no per-position summary can (reads with an N at a position that varies from tile to
tile). The quality rule's cut is always computed and logged next to it.
Each of the three stops prints the same numbers (overlap_stop_text): how far the
reads reach after primer removal, what merging needs (expected amplicon length +
min_overlap), where the expected length came from (amplicon.expected_length, or
the probe statistic named in amplicon.probe_length_stat.<key>), and what to check.

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
import bisect
import gzip
import itertools
import json
import math
import sys
import time
from pathlib import Path
from statistics import median


# ── Settings of the checks added for binned data and short products ──
# How many reads to read from the start of each primer-trimmed file to see which
# quality values it carries. Every value a file uses appears within its first 1,000
# reads (measured on the MiSeq and NextSeq test files), so 10,000 is plenty for the
# values themselves. It is NOT enough to judge quality: NextSeq files are written
# tile by tile, so the start of a file is one tile, not the whole run.
QUALITY_SCAN_READS = 10_000

# The classes of quality values: the same thresholds as rule choose_error_model
# (60c_choose_error_model.py), so that this step and the error model act on the same
# answer. Q0-Q2 are the no-call values of every Illumina table (an N gets Q2), so they
# are not counted. Fewer than 3 values above Q2: no usable quality information;
# 3 to 10 with no two neighbours (such as 29 and 30): binned; 20 or more: continuous;
# anything else: undetermined.
NO_CALL_MAX = 2
MIN_QUALITY_VALUES = 3
MAX_BINNED_VALUES = 10
MIN_CONTINUOUS_VALUES = 20

# Short off-target products (rule of change b). A read shorter than this share of the
# expected amplicon length is left out when the read-length limit is computed.
# Why one half:
#   - Cutadapt removes the other primer from a read that runs through the end of its
#     product, so a primer-trimmed read shorter than the sequenced length IS its whole
#     product. A read under half the expected length therefore comes from a product
#     less than half as long as the target: a primer dimer or an unspecific product,
#     not natural length variation of the marker (ITS, the one marker that varies that
#     much, never uses this limit).
#   - Such a product could only give an ASV far below the ASV length window
#     (probe q1 - 50 bp), which removes it later anyway.
#   - Target reads that do NOT run through (the amplicon is longer than the read) are
#     always longer than this: R1 + R2 must exceed the expected length for the pair to
#     merge at all, so each read reaches more than half of it in any design that works.
OFFTARGET_FRACTION_OF_EXPECTED = 0.5
# ... unless the short reads are most of the sample: then they are not a minority
# product but the sample itself (wrong primers, wrong expected length), and the
# old limit, over all reads, is kept for that sample, with a warning.
MAX_OFFTARGET_SHARE_PCT = 50.0

# method: simulate. How many read pairs to take from each sample, spread evenly over
# the whole file. 20,000 pairs measure the share of pairs kept to within about
# +/- 0.7 points per sample (95 % range), and much better pooled over samples.
# The same number from every sample means every sample has the same say in the
# choice, however deep it was sequenced (a sample with fewer pairs gives all it has).
SIMULATE_PAIRS_PER_SAMPLE = 20_000
# method: simulate. Cut pairs that keep almost as many pairs as the best one count as
# equally good, and among them the one with the longest R1 + R2 is taken: in points
# of the share of sampled pairs kept. Why: the share kept is often nearly flat over a
# wide range of cuts (on a MiSeq COI library, R1 34 / R2 183 kept 99.59 % and 57 / 200
# kept 99.28 %), so the strict maximum is decided by differences smaller than the
# sampling noise (about 0.1-0.3 points with 20,000 pairs per sample) and lands on
# lopsided cuts, such as a 42 bp R1, that leave most of the amplicon read by one read
# only. The longer total adds overlap, the stretch where mergePairs checks R1
# against R2 and rejects pairs that disagree, at a cost of at most this many points.
SIMULATE_TOLERANCE_PCT = 0.5


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


# ── How many reads a file holds, from its falco report ───────────
def parse_total_sequences(fastqc_data_path: Path) -> int | None:
    """The number of reads in the file, from the "Total Sequences" line of falco's
    ">>Basic Statistics" block (falco counts every read of the file).

    Used by method: simulate to spread its sample of read pairs evenly over the whole
    file without reading the file twice. R1 and R2 of a sample hold the same number of
    reads (Cutadapt keeps or drops whole pairs), so the R1 report is enough.
    Returns None if the line is missing; the caller then takes every pair.
    """
    with fastqc_data_path.open() as fh:
        for line in fh:
            if line.startswith("Total Sequences"):
                return int(line.rstrip("\n").split("\t")[1])
            if line.startswith(">>END_MODULE"):
                break
    return None


# ── Open a FASTQ file and walk through its reads ─────────────────
def open_fastq(path):
    """Open a FASTQ file for reading in bytes, gzipped (03.stripped/*.fastq.gz) or not.

    gzip is recognised by its first two bytes, not by the file name, so a test file
    written without compression works too.
    """
    path = Path(path)
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return path.open("rb")


def fastq_records(fh):
    """Yield (sequence, quality) of each read of an open FASTQ file, as bytes.

    A FASTQ record is four lines: header, bases, "+", quality letters. Only the bases
    and the quality letters are needed here.
    """
    while True:
        header = fh.readline()
        if not header:
            return
        seq = fh.readline().rstrip(b"\r\n")
        fh.readline()
        qual = fh.readline().rstrip(b"\r\n")
        yield seq, qual


# ── Change (a): which quality values do the reads carry? ─────────
def scan_quality_values(fastq_path, max_reads: int = QUALITY_SCAN_READS) -> dict[int, int]:
    """Count how many bases carry each quality value in the first reads of a file.

    Input : one primer-trimmed FASTQ file (03.stripped/{sample}_R1/R2.fastq.gz, the
            reads DADA2 will be handed).
    Output: {Phred value: number of bases}, used by main() to tell continuous MiSeq
            values from binned ones and to explain what q_threshold means on them.
    Quality letters are read as Phred+33, the encoding of every current Illumina FASTQ
    (DADA2's own reader assumes the same for these files).
    """
    letters: dict[int, int] = {}
    with open_fastq(fastq_path) as fh:
        for n_read, (_seq, qual) in enumerate(fastq_records(fh)):
            if n_read >= max_reads:
                break
            # bytes.count runs in C, so counting per distinct letter of the read is
            # much faster than looping over every base in Python.
            for letter in set(qual):
                letters[letter] = letters.get(letter, 0) + qual.count(letter)
    return {letter - 33: n for letter, n in letters.items()}


def quality_value_class(values_above_q2) -> str:
    """Say what kind of quality values a set of values above Q2 is.

    The same rule, word for word, as classify() in 60c_choose_error_model.py, which
    picks DADA2's error function from the same values:
      "none"         fewer than 3 values (rewritten qualities, no information);
      "binned"       3 to 10 values, no two of them neighbours (no instrument table has
                     neighbouring values);
      "continuous"   20 or more (a MiSeq writes 2 and every value from 7 to 38);
      "undetermined" anything else (for example continuous values cut by a high
                     trunc_q, which keep neighbours).
    """
    values = sorted(set(values_above_q2))
    count = len(values)
    if count < MIN_QUALITY_VALUES:
        return "none"
    if count <= MAX_BINNED_VALUES:
        for lower, upper in zip(values, values[1:]):
            if upper - lower == 1:
                return "undetermined"
        return "binned"
    if count < MIN_CONTINUOUS_VALUES:
        return "undetermined"
    return "continuous"


def values_text(values) -> str:
    """A set of quality values as a short text: "4 values (2, 9, 20, 34)", or
    "33 values (2-38)" when there are more than ten."""
    values = sorted(values)
    if not values:
        return "no values"
    if len(values) > 10:
        return f"{len(values)} values ({values[0]}-{values[-1]})"
    return f"{len(values)} values ({', '.join(str(v) for v in values)})"


def threshold_meaning(q_threshold: float, values_above_q2, value_class: str) -> tuple[str, list[str]]:
    """What "median Q1 below q_threshold" asks of the reads, and warnings about it.

    Input : the threshold from the config, the quality values above Q2 found by
            scan_quality_values (pooled over the samples of one read direction), and
            their class from quality_value_class.
    Output: (one sentence for the log, a list of warning sentences).

    Why it matters: Q1 at a position is the lowest value that at least a quarter of
    the bases there are at or below. On binned data a single position can only hold
    one of the few values, so the threshold can only ask "are a quarter of the bases
    at one of the values below it?". Every threshold between two neighbouring values
    asks the same question, and a threshold EQUAL to a value sits on the edge between
    two questions: one unit higher already counts that value as low. This is what
    happens with the default 20 on NextSeq 1000/2000 600-cycle data (9, 20, 34).
    Falco's 5-position group means and the median across samples can still fall
    between the values, so the cut can also move a little within one question.
    """
    t = q_threshold
    values = sorted(set(values_above_q2))
    if value_class == "continuous":
        return (f"the cut falls where at least a quarter of the bases are below Q{t:g}", [])
    if value_class != "binned":
        return ("not explained: the values are neither continuous nor clearly binned", [])

    below = [v for v in values if v < t]
    warnings: list[str] = []
    if not below:
        meaning = (f"no value above the no-call values lies below {t:g}, so only a quarter "
                   f"of no-calls (N) could set a quality cut")
        warnings.append(f"q_threshold {t:g} is at or below the lowest quality value "
                        f"({values[0]}); quality will practically never set the cut.")
        return meaning, warnings
    if t > values[-1]:
        meaning = (f"every value lies below {t:g}, so the median Q1 is below it at the first "
                   f"position: the quality cut is 0 and the final cuts come from "
                   f"resolve_policy alone")
        warnings.append(f"q_threshold {t:g} is above the highest quality value "
                        f"({values[-1]}); the quality profile cannot set the cut.")
        return meaning, warnings
    meaning = (f"the cut falls where at least a quarter of the bases carry {below[-1]} or "
               f"lower (the values below {t:g}: {', '.join(str(v) for v in below)})")
    if t in values:
        above = [v for v in values if v > t]
        warnings.append(
            f"q_threshold {t:g} equals one of the quality values. A threshold one unit "
            f"higher ({t + 1:g}) would count {t:g} as low too and ask whether a quarter of "
            f"the bases carry {t:g} or lower, a different question that can move the cut by "
            f"tens of bases (the values above it: {', '.join(str(v) for v in above)}).")
    return meaning, warnings


# ── Change (b): reads too short to be the target do not set the limit ─
def offtarget_floor(exp_len: int) -> int:
    """The shortest read that can still be part of the target amplicon: half the
    expected amplicon length, rounded up (see OFFTARGET_FRACTION_OF_EXPECTED for why).

    Input : the expected amplicon length from resolve_expected_length (probe p95 or the
            user's number). Output: a length in bp, used by read_coverage_cap_with_floor.
    """
    return math.ceil(exp_len * OFFTARGET_FRACTION_OF_EXPECTED)


def read_coverage_cap_with_floor(length_counts: list[tuple[int, float]], min_pct: float,
                                 floor: int | None) -> tuple[int | None, float, bool]:
    """A sample's read-length limit, leaving out reads shorter than `floor`.

    Input : the sample's read-length distribution (parse_length_distribution), the
            share of reads the limit must reach (min_read_coverage_pct) and the floor
            from offtarget_floor (None for ITS: the old limit over all reads).
    Output: (cap, share of the sample's reads below the floor in %, applied?).
            The cap goes to cap_bins_at like the old one; the share and the flag are
            logged per sample by main().

    Why: in a sample with, say, 20 % of its reads in a 158 bp off-target product, the
    length 95 % of ALL reads reach falls inside that product, and because the run takes
    the strictest sample, one such sample caps both reads of the whole run too short to
    overlap. Those short reads are lost at any cut that can merge, so the limit should
    describe the other reads. A read-length bin of falco's ("214-218") counts at its
    shortest length, as in parse_length_distribution.
    The floor is NOT applied when more than MAX_OFFTARGET_SHARE_PCT of the reads lie
    below it (then the "short product" is the sample itself), and changes nothing when
    no read lies below it: such samples get exactly the old cap.
    """
    total = sum(n for _length, n in length_counts)
    if floor is None or total <= 0:
        return read_coverage_cap(length_counts, min_pct), 0.0, False
    kept = [(length, n) for length, n in length_counts if length >= floor]
    share_out = 100.0 * (total - sum(n for _length, n in kept)) / total
    if share_out == 0.0 or share_out > MAX_OFFTARGET_SHARE_PCT:
        return read_coverage_cap(length_counts, min_pct), share_out, False
    return read_coverage_cap(kept, min_pct), share_out, True


# ── Change (c), method: simulate — what filterAndTrim would do to one read ─
# Error probability of each quality letter, looked up by its byte value (Phred+33):
# a base with quality Q is wrong with probability 10^(-Q/10). Letters below "!" never
# occur; they get 0 so that the table can be indexed by any byte.
ERROR_BY_LETTER = [10.0 ** (-(letter - 33) / 10.0) if letter >= 33 else 0.0
                   for letter in range(256)]


def longest_passing_cut(seq: bytes, qual: bytes, trunc_q: int, max_ee: float,
                        max_len: float | None = None) -> int:
    """The longest truncLen at which this read would still pass filterAndTrim.

    Input : one read (bases and quality letters) and the run's filter settings, the same
            ones dada_filter (60b_dada_filter.R) hands to filterAndTrim: trunc_q,
            this direction's max_ee, max_len (None = no limit), maxN fixed at 0.
    Output: a length M. The read passes at every truncLen from min_len up to M and
            fails above it; 0 means it fails at every cut.

    What filterAndTrim does to a read, in its own order, and why one number is enough:
      1. maxLen, on the whole read: a longer read is dropped at any cut.
      2. truncQ: the read is cut just before its first base with quality <= trunc_q.
      3. truncLen L: a read now shorter than L is dropped, otherwise cut to L.
      4. maxN = 0: dropped if an N remains in the first L bases.
      5. maxEE: dropped if the expected errors of the first L bases (the sum of their
         error probabilities) exceed max_ee.
    Every condition gets harder as L grows (expected errors only add up), so each read
    passes exactly up to the smallest of: the bases before its first low-quality base,
    the bases before its first N, and the longest start of the read whose expected
    errors stay within max_ee. The same reasoning reproduced DADA2's own reads.out on
    26 of 26 MiSeq samples when it was first checked.
    """
    if max_len is not None and len(seq) > max_len:
        return 0
    # Step 2: position of the first base with quality <= trunc_q (bytes.find runs in C).
    usable = len(qual)
    for letter in range(33, 34 + trunc_q):
        pos = qual.find(bytes([letter]), 0, usable)
        if pos != -1:
            usable = pos
    # Step 4: bases before the first N.
    first_n = seq.find(b"N")
    if first_n != -1 and first_n < usable:
        usable = first_n
    # Step 5: running sum of error probabilities over the usable bases, left to right as
    # DADA2 adds them; bisect counts how many starts of the read stay within max_ee.
    running_ee = list(itertools.accumulate(ERROR_BY_LETTER[letter] for letter in qual[:usable]))
    return bisect.bisect_right(running_ee, max_ee)


def simulate_sample(r1_path, r2_path, n_reads: int | None, n_wanted: int, trunc_q: int,
                    max_ee_r1: float, max_ee_r2: float, max_len: float | None) -> list[tuple[int, int]]:
    """The longest passing cuts (R1, R2) of read pairs spread evenly over one sample.

    Input : the sample's primer-trimmed R1 and R2 files (03.stripped/), the number of
            reads in them (parse_total_sequences; None = take every pair), how many
            pairs to take, and the filter settings.
    Output: one (M_R1, M_R2) per sampled pair (longest_passing_cut of each read); a pair
            passes filterAndTrim at (L1, L2) exactly when L1 <= M_R1 and L2 <= M_R2.
            main() pools these over the samples into survival_table.

    The pairs are spread over the WHOLE file, never taken from its start: a NextSeq
    writes its reads tile by tile, and in one soil R2 file the first 100,000 reads carry
    no N while 19.9 % of all its reads do. Pair j is read number floor(j * n / n_wanted),
    so exactly n_wanted pairs are taken, evenly spaced from the first read to the last.
    R1 and R2 are read side by side: DADA2 pairs them by their place in the two files.
    """
    pairs: list[tuple[int, int]] = []
    take_all = n_reads is None or n_reads <= n_wanted
    j = 0
    next_pick = 0
    with open_fastq(r1_path) as fh1, open_fastq(r2_path) as fh2:
        for index, (rec1, rec2) in enumerate(zip(fastq_records(fh1), fastq_records(fh2))):
            if not take_all:
                if j >= n_wanted:
                    break
                if index < next_pick:
                    continue
                j += 1
                next_pick = (j * n_reads) // n_wanted
            m1 = longest_passing_cut(rec1[0], rec1[1], trunc_q, max_ee_r1, max_len)
            m2 = longest_passing_cut(rec2[0], rec2[1], trunc_q, max_ee_r2, max_len)
            pairs.append((m1, m2))
    return pairs


def survival_table(pairs: list[tuple[int, int]], longest: int) -> list[list[int]]:
    """kept[L1][L2] = how many of the sampled pairs pass at truncLen (L1, L2).

    Input : the (M_R1, M_R2) of every sampled pair (simulate_sample, all samples
            pooled) and the longest cut to tabulate (the longest M seen).
    Output: a (longest + 2) x (longest + 2) table of counts, read by best_simulated_cuts
            and by main() for the quality rule's cut.

    How: first count the pairs at each exact (M_R1, M_R2), then add up from the long
    end in both directions, so each cell holds the pairs with M_R1 >= L1 and
    M_R2 >= L2 (a "2-D running total"; the last row and column stay 0).
    """
    size = longest + 2
    kept = [[0] * size for _ in range(size)]
    for m1, m2 in pairs:
        kept[m1][m2] += 1
    for l1 in range(longest, -1, -1):
        row, row_below = kept[l1], kept[l1 + 1]
        for l2 in range(longest, -1, -1):
            row[l2] += row[l2 + 1] + row_below[l2] - row_below[l2 + 1]
    return kept


def kept_at(kept: list[list[int]], l1: int, l2: int) -> int:
    """Pairs kept at (l1, l2) from survival_table; 0 beyond the longest cut tabulated."""
    size = len(kept)
    if l1 >= size or l2 >= size:
        return 0
    return kept[l1][l2]


def best_simulated_cuts(kept: list[list[int]], need: int, min_len: int,
                        tolerance: float = 0.0) -> tuple[int, int, int] | None:
    """The R1/R2 cuts to use: the longest pair among those that keep (nearly) the most
    sampled pairs while L1 + L2 >= need.

    Input : survival_table's counts, the length the two cuts must reach together
            (expected length + min_overlap + simulate_extra_overlap), min_len (a cut
            below it would fail filterAndTrim's minLen), and the tolerance in pairs
            (SIMULATE_TOLERANCE_PCT of the sampled pairs; 0 = the strict maximum).
    Output: (L1, L2, pairs kept), or None when no cut pair reaches `need`.

    Two steps:
      1. The most pairs any cut pair reaching `need` keeps. Fewer bases never lose a
         pair, so for each L1 it is enough to look at the shortest L2 that reaches.
      2. Among all cut pairs reaching `need` that keep at least that many minus the
         tolerance, the one with the longest R1 + R2 (the most overlap for mergePairs
         to check), then the longer R1 (on Illumina the better read).
    """
    longest = len(kept) - 2
    shortest_cut = max(1, min_len)
    most = None
    for l1 in range(shortest_cut, longest + 1):
        l2 = max(shortest_cut, need - l1)
        if l2 <= longest and (most is None or kept[l1][l2] > most):
            most = kept[l1][l2]
    if most is None:
        return None
    best = None
    for l1 in range(shortest_cut, longest + 1):
        for l2 in range(max(shortest_cut, need - l1), longest + 1):
            if kept[l1][l2] < most - tolerance:
                break                    # longer L2 only keeps fewer pairs
            candidate = (l1 + l2, l1)
            if best is None or candidate > best[0]:
                best = (candidate, l1, l2)
    _rank, l1, l2 = best
    return l1, l2, kept[l1][l2]


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


# ── Read a setting that must be a whole number ───────────────────
def whole_number(value) -> int | None:
    """Return `value` as an int when it is a whole number, else None.

    Accepts 20, 20.0 and the text "20" (an inline --config value can arrive as text);
    refuses true/false (Python would count them as 1 and 0), fractions, other text
    and an empty value. Used for amplicon.trunc_len.simulate_extra_overlap; main()
    stops with a readable message on None.
    """
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        try:
            return int(value.strip())
        except ValueError:
            return None
    return None


# ── Save a result dictionary to a JSON file ──────────────────────
def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")


# ── Say where the expected amplicon length came from (for a stop text) ─
def describe_expected_length(exp_len: int, exp_len_source: str, probe_json_path,
                             probe_stat_key: str) -> tuple[str, str]:
    """Two sentences for overlap_stop_text below: what the expected length is, and
    which setting changes it.

    The overlap check measures the reads against one number, the expected amplicon
    length, and a user who meets the stop needs to know where that number came from
    before deciding whether it or the reads are at fault.

    Input : the length and its source label, both from resolve_expected_length above
            ("auto:<stat>", "manual:int" or "manual:range_max"); the probe JSON (only
            read for "auto"; written by rule amplicon_probe, 40a_amplicon_probe.py); and
            the key of amplicon.probe_length_stat this marker reads (from the marker
            pack via 50_trunclen.smk; 16S for 16S and for gyrB).
    Output: (what_it_is, what_to_change), two sentences of plain text.
    """
    if exp_len_source.startswith("auto:"):
        stat = exp_len_source.split(":", 1)[1]
        setting = f"amplicon.probe_length_stat.{probe_stat_key}"
        probe = json.loads(Path(probe_json_path).read_text())
        if probe.get("probe_mode") == "pcr":
            measured = f"the {probe['n_amplicons']:,} in-silico PCR products of these primers"
        else:
            measured = (f"the {probe['n_amplicons']:,} reference sequences, which are already "
                        f"cut to the amplicon")
        what_it_is = (
            f"{exp_len} bp, the {stat} of the lengths of {measured} "
            f"(amplicon.expected_length: auto; {setting}: {stat}). The probe measured "
            f"q1 {probe['q1']}, median {probe['median']}, q3 {probe['q3']}, "
            f"p95 {probe['p95']} and p99 {probe['p99']} bp."
        )
        # A lower statistic is offered only down to the median: below it the check
        # would be sized for fewer than half of the amplicons.
        order = ["min", "q1", "median", "q3", "p95", "p99", "max"]
        lower = [s for s in ("q3", "median")
                 if stat in order and order.index(s) < order.index(stat)]
        what_to_change = (f"If the amplicon between your primers is shorter than {exp_len} bp, "
                          f"set amplicon.expected_length to its length in bp")
        if lower:
            offered = ", ".join(f"{s}: {probe[s]} bp" for s in lower)
            what_to_change += (f", or choose a lower statistic of the probe under {setting} "
                               f"({offered})")
        what_to_change += "."
    else:
        if exp_len_source == "manual:range_max":
            where = "the upper end of the range set in amplicon.expected_length"
        else:
            where = "the number set in amplicon.expected_length"
        what_it_is = (f"{exp_len} bp, {where} (amplicon.probe_length_stat is not used "
                      f"with a length of your own).")
        what_to_change = (
            f"If the amplicon between your primers is shorter than {exp_len} bp, correct "
            f"amplicon.expected_length, or set it to auto to have it measured from the "
            f"reference."
        )
    return what_it_is, what_to_change


# ── The text of a stop at the overlap check ─────────────────────────
def overlap_stop_text(first_line: str, exp_len: int, exp_len_source: str, probe_json_path,
                      probe_stat_key: str, min_overlap: int, max_r1: int, max_r2: int,
                      min_read_coverage_pct: float, policy_hint: str = "") -> str:
    """The message pick_trunclen prints when no pair of cuts passes the overlap check.

    Used by the three resolve_policy branches of main() (raise_trunc, relax_q, error),
    which each say in `first_line` how far their cuts got. The rest is the same for
    all three, so that every stop gives the numbers behind the check and the settings
    that feed it:
      - the read length available: how far the reads reach after primer removal,
        R1 + R2 (the ceilings left by the read-coverage cap; the log lines above name
        the sample that sets each one);
      - what merging needs: the expected amplicon length plus amplicon.min_overlap;
      - where the expected length came from (describe_expected_length above);
      - what to check, including the case where the expected length is right and the
        reads are simply too short, which no setting fixes.
    `policy_hint` is an extra line for relax_q and error when the reads themselves
    would be long enough (raise_trunc would then pass).
    Output: one text of several lines; main() writes it to the log and to stderr.
    """
    min_total = exp_len + min_overlap
    available = max_r1 + max_r2
    what_it_is, what_to_change = describe_expected_length(
        exp_len, exp_len_source, probe_json_path, probe_stat_key)
    lines = [
        f"[pick_trunclen] ERROR: {first_line}",
        f"  Read length available: after primer removal the reads reach at most "
        f"R1 {max_r1} + R2 {max_r2} = {available} bp (the length that "
        f"{min_read_coverage_pct:g}% of each sample's reads reach, "
        f"amplicon.trunc_len.min_read_coverage_pct; the lines above name the sample that "
        f"sets each).",
        f"  Merging needs: {exp_len} bp (expected amplicon length) + {min_overlap} bp of "
        f"overlap (amplicon.min_overlap) = {min_total} bp.",
        f"  Expected amplicon length: {what_it_is}",
        "  What to check:",
        f"  - {what_to_change}",
    ]
    if policy_hint:
        lines.append(f"  - {policy_hint}")
    if available < min_total:
        lines.append(
            f"  - If {exp_len} bp is right, these reads are too short for R1 and R2 to "
            f"overlap by {min_overlap} bp, and no truncation setting makes the pairs merge."
        )
    return "\n".join(lines)


# ── The quality rule's overlap check: keep the pair mergeable ───────
def resolve_overlap(resolve_policy: str, primary_r1: int, primary_r2: int, max_r1: int,
                    max_r2: int, agg_r1, agg_r2, min_total: int, q_threshold: float,
                    q_floor: float, stop_args: dict, log) -> dict:
    """Apply resolve_policy when the quality cuts are too short to overlap.

    Input : the quality cuts (primary_r1/r2), the read-length limits (max_r1/r2), the
            median Q1 per position (agg_r1/r2, from aggregate_q1_across_samples), the
            length the two cuts must reach together (expected length + min_overlap),
            the thresholds, the shared arguments of overlap_stop_text, and the log.
    Output: a dict with the final cuts (cut_r1, cut_r2), the threshold used (q_used),
            how they were reached (resolved_via: primary, raise_trunc or relax_q), and,
            when no pair of cuts can merge, exit_code (3, 4, 5 or 6) and stop_text.
            exit_code 0 means the cuts are fine.

    main() acts on the result: under method: quality a stop ends the run as before;
    under method: simulate the quality rule only serves as the comparison logged next
    to the simulated cut, so its stop is reported, not acted on. That is why this
    function returns the stop text instead of writing it as an error itself.
    """
    cut_r1, cut_r2 = primary_r1, primary_r2
    result = {"cut_r1": cut_r1, "cut_r2": cut_r2, "q_used": q_threshold,
              "resolved_via": "primary", "exit_code": 0, "stop_text": ""}
    slack = (primary_r1 + primary_r2) - min_total
    if slack >= 0:
        return result

    # One line for relax_q and error: true when the reads themselves are long enough,
    # so that raise_trunc (which extends the cuts up to the ceilings) would pass.
    raise_trunc_would_pass = (max_r1 + max_r2) >= min_total
    raise_trunc_hint = (
        f"The reads themselves reach far enough ({max_r1 + max_r2} bp): "
        f"amplicon.trunc_len.resolve_policy: raise_trunc (the default) extends the cuts "
        f"past the quality drop instead."
    )

    log(f"[pick_trunclen] Constraint violated by {-slack} bp; resolve_policy={resolve_policy}")
    if resolve_policy == "raise_trunc":
        headroom_r1 = max_r1 - cut_r1
        headroom_r2 = max_r2 - cut_r2
        deficit = -slack
        if headroom_r1 + headroom_r2 < deficit:
            # Even cut at the ceilings the reads do not reach: say so with the
            # numbers and the settings behind them, then stop (exit code 3).
            short_by = min_total - (max_r1 + max_r2)
            result["exit_code"] = 3
            result["stop_text"] = overlap_stop_text(
                f"R1 and R2 cannot be cut so that the pairs still merge. Even at their "
                f"full usable length (resolve_policy: raise_trunc) the reads reach "
                f"{max_r1 + max_r2} bp, {short_by} bp short of the {min_total} bp "
                f"that merging needs.",
                **stop_args)
            return result

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

        result["cut_r1"] = cut_r1 + bump_r1
        result["cut_r2"] = cut_r2 + bump_r2
        result["resolved_via"] = "raise_trunc"

        def _mean_q(vals: list[float]) -> str:
            finite = [v for v in vals if v != float("-inf")]
            return f"{sum(finite) / len(finite):.1f}" if finite else "n/a"

        log(f"[pick_trunclen] Raised: R1 +{bump_r1} (mean Q1 of added bases "
            f"{_mean_q(gained_q_r1)}), R2 +{bump_r2} (mean Q1 {_mean_q(gained_q_r2)}); "
            "allocation is quality-weighted, not headroom-weighted")
        return result

    if resolve_policy == "relax_q":
        q_try = q_threshold
        # The cuts at the last threshold tried, for the stop text: the primary
        # cuts when the loop below never runs (q_threshold already at the floor).
        c1, c2 = primary_r1, primary_r2
        while q_try > q_floor:
            q_try -= 1
            c1 = first_drop_trunclen(agg_r1, q_try)
            c2 = first_drop_trunclen(agg_r2, q_try)
            if c1 + c2 >= min_total:
                result.update(cut_r1=c1, cut_r2=c2, q_used=q_try, resolved_via="relax_q")
                log(f"[pick_trunclen] Relaxed to Q1>={q_try}: R1={c1}, R2={c2}")
                return result
        # Even at the lowest threshold allowed the quality cuts do not reach:
        # say how far they got, with the numbers behind the check (exit code 4).
        result["exit_code"] = 4
        result["stop_text"] = overlap_stop_text(
            f"R1 and R2 cannot be cut so that the pairs still merge. With the "
            f"quality threshold lowered to Q{q_try:g} (amplicon.trunc_len.q_floor: "
            f"{q_floor:g}; resolve_policy: relax_q) the cuts reach R1 {c1} + R2 {c2} "
            f"= {c1 + c2} bp, {min_total - (c1 + c2)} bp short of the {min_total} bp "
            f"that merging needs.",
            policy_hint=raise_trunc_hint if raise_trunc_would_pass else "",
            **stop_args)
        return result

    if resolve_policy == "error":
        # This policy never moves the cuts: report the quality cuts with the numbers
        # behind the check and stop (exit code 5).
        result["exit_code"] = 5
        result["stop_text"] = overlap_stop_text(
            f"R1 and R2 cannot be cut so that the pairs still merge. The quality cuts "
            f"at Q{q_threshold:g} (amplicon.trunc_len.q_threshold) reach R1 {primary_r1} "
            f"+ R2 {primary_r2} = {primary_r1 + primary_r2} bp, {-slack} bp short of the "
            f"{min_total} bp that merging needs, and resolve_policy: error stops here.",
            policy_hint=raise_trunc_hint if raise_trunc_would_pass else "",
            **stop_args)
        return result

    result["exit_code"] = 6
    result["stop_text"] = f"[pick_trunclen] ERROR: unknown resolve_policy={resolve_policy!r}"
    return result


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
    # Which key of amplicon.probe_length_stat this marker reads (16S, 18S, rpoB, COI;
    # gyrB reads 16S). Only used to name that setting in the stop texts below.
    probe_stat_key = str(p.probe_stat_key)
    min_read_coverage_pct = float(p.min_read_coverage_pct)

    log(f"[pick_trunclen] amp_type={amp_type} mode={mode} q_threshold={q_threshold}")
    log(f"[pick_trunclen] R1 falco files: {len(falco_r1)}; R2 falco files: {len(falco_r2)}")

    # ── Manual mode: use the lengths provided directly in config ─
    if mode == "manual":
        if manual_r1 is None or manual_r2 is None:
            log("[pick_trunclen] ERROR: trunc_len.mode=manual but manual_r1/manual_r2 are not set")
            return 2
        log(f"[pick_trunclen] Manual mode: r1={manual_r1} r2={manual_r2}")
        # The expected amplicon length is recorded in manual mode too. The manual cuts
        # do not use it, but rule dada_length_filter does: with length_filter.mode: auto
        # and no probe (amplicon.expected_length set to a number or a [min, max] range),
        # its window is that length -15 % / +15 %, read from this file. Without it, that
        # combination crashed in dada_length_filter with a bare KeyError, after the
        # whole DADA2 run. The length is resolved exactly as in auto mode below
        # (resolve_expected_length: a number as it is, a range as its max, 'auto' as the
        # probe's statistic). A value that cannot be resolved does not stop this step,
        # because the cuts do not need it: it is written as null with a note, and
        # dada_length_filter stops with its own message only if it needs the value.
        # ITS never has one (see the ITS branch of auto mode).
        manual_exp_len = None
        manual_exp_source = None
        if amp_type == "ITS":
            manual_exp_source = "not_applicable:ITS_truncLen_overridden"
        else:
            manual_probe_json = None
            if hasattr(sm.input, "probe_json") and sm.input.probe_json:
                pj = sm.input.probe_json
                manual_probe_json = pj[0] if isinstance(pj, list) else pj
            try:
                manual_exp_len, manual_exp_source = resolve_expected_length(
                    expected_length, amp_type, manual_probe_json, log,
                    probe_stat=probe_stat,
                )
                log(f"[pick_trunclen] expected_length: {manual_exp_len} (source: {manual_exp_source}); "
                    f"the manual cuts do not use it, the length filter can")
            except (ValueError, FileNotFoundError, KeyError) as problem:
                manual_exp_source = "unresolved"
                log(f"[pick_trunclen] NOTE: amplicon.expected_length could not be resolved, so it is "
                    f"recorded as null ({problem}). The manual cuts do not need it; "
                    f"length_filter.mode: auto without a probe would.")
        write_json(out_json, {
            "r1": int(manual_r1), "r2": int(manual_r2),
            "mode": "manual", "amplicon_type": amp_type,
            "expected_length": manual_exp_len,
            "expected_length_source": manual_exp_source,
            "resolved_via": "manual",
        })
        return 0

    # ── Auto mode: check the settings of the method and of the filter ─
    # method: quality (the default) cuts from falco's quality profile; simulate
    # simulates filterAndTrim on sampled read pairs (see the module docstring). Both
    # are checked here, like resolve_policy, so a typo stops with a readable message.
    method = str(p.method).strip().lower()
    if method not in ("quality", "simulate"):
        log(f"[pick_trunclen] ERROR: amplicon.trunc_len.method must be quality or simulate, "
            f"got {p.method!r}")
        return 2
    # The safety margin the simulate method adds to the overlap rule, in bp.
    extra_overlap = whole_number(p.simulate_extra_overlap)
    if extra_overlap is None or extra_overlap < 0:
        log(f"[pick_trunclen] ERROR: amplicon.trunc_len.simulate_extra_overlap must be a "
            f"whole number of bases >= 0, got {p.simulate_extra_overlap!r}")
        return 2
    # The filter settings dada_filter will hand to filterAndTrim (amplicon.dada2), which
    # the simulate method replays read by read. max_ee is [R1, R2] in the template; a
    # single number counts for both reads, as in filterAndTrim.
    trunc_q = int(p.trunc_q)
    max_ee = p.max_ee if isinstance(p.max_ee, (list, tuple)) else [p.max_ee, p.max_ee]
    max_ee_r1, max_ee_r2 = float(max_ee[0]), float(max_ee[-1])
    min_len = int(p.min_len)
    max_len = None if p.max_len is None else float(p.max_len)
    log(f"[pick_trunclen] method={method}")

    # ── Change (a): which quality values do the reads carry? ─────
    # Read the first QUALITY_SCAN_READS reads of every primer-trimmed file, pool the
    # values per read direction, and say in the log what kind they are and what
    # q_threshold means on them. Nothing here changes a cut: it explains it.
    stripped_r1 = [Path(x) for x in sm.input.stripped_r1]
    stripped_r2 = [Path(x) for x in sm.input.stripped_r2]
    for label, paths in (("R1", stripped_r1), ("R2", stripped_r2)):
        pooled: dict[int, int] = {}
        for path in paths:
            for value, n_bases in scan_quality_values(path).items():
                pooled[value] = pooled.get(value, 0) + n_bases
        values = sorted(pooled)
        above_q2 = [v for v in values if v > NO_CALL_MAX]
        value_class = quality_value_class(above_q2)
        log(f"[pick_trunclen] {label} quality values (first {QUALITY_SCAN_READS:,} reads of "
            f"each of the {len(paths)} primer-trimmed files): {values_text(values)}, "
            f"{value_class}")
        meaning, warnings = threshold_meaning(q_threshold, above_q2, value_class)
        log(f"[pick_trunclen] {label} q_threshold {q_threshold:g} on these values: {meaning}")
        for warning in warnings:
            log(f"[pick_trunclen] WARNING: {label} {warning}")

    # ── Resolve the expected amplicon length (all markers but ITS) ─
    # Needed before the read-length limits below (change b) as well as by the overlap
    # check. ITS has no single expected length (its truncLen is overridden to 0 by
    # dada_filter), so it keeps the old limits and skips the overlap check.
    probe_json_path = None
    exp_len = exp_len_source = None
    if amp_type != "ITS":
        if hasattr(sm.input, "probe_json") and sm.input.probe_json:
            # Snakemake gives a list (or string) here depending on how it was declared.
            pj = sm.input.probe_json
            probe_json_path = pj[0] if isinstance(pj, list) else pj
        exp_len, exp_len_source = resolve_expected_length(
            expected_length, amp_type, probe_json_path, log,
            probe_stat=probe_stat,
        )
        log(f"[pick_trunclen] expected_length (for overlap): {exp_len} (source: {exp_len_source}), min_overlap: {min_overlap}")

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

    # Change (b): reads shorter than half the expected amplicon length belong to a short
    # off-target product and are left out of each sample's limit (offtarget_floor and
    # read_coverage_cap_with_floor explain why). Each sample's limit and the share of
    # its reads left out are logged, so a reader can see which samples carry the product.
    floor = offtarget_floor(exp_len) if exp_len is not None else None
    caps_r1: list[int | None] = []
    caps_r2: list[int | None] = []
    for label, falco_paths, caps in (("R1", falco_r1, caps_r1), ("R2", falco_r2, caps_r2)):
        if floor is not None:
            log(f"[pick_trunclen] {label} read-length limit per sample: the length "
                f"{min_read_coverage_pct:g}% of its reads reach, not counting reads shorter "
                f"than {floor} bp (half the expected amplicon length: a short off-target "
                f"product, not the target)")
        for path in falco_paths:
            cap, share_out, applied = read_coverage_cap_with_floor(
                parse_length_distribution(path), min_read_coverage_pct, floor)
            caps.append(cap)
            if floor is None:
                continue
            name = path.parent.name
            if applied:
                log(f"[pick_trunclen]   {name}: {cap} bp ({share_out:.2f}% of its reads are "
                    f"shorter than {floor} bp and were left out)")
            elif share_out > MAX_OFFTARGET_SHARE_PCT:
                log(f"[pick_trunclen] WARNING: {name}: {share_out:.1f}% of its reads are "
                    f"shorter than {floor} bp, more than half of the sample, so they are not "
                    f"treated as an off-target product: its limit, {cap} bp, counts all its "
                    f"reads. Check its primers and amplicon.expected_length.")
            else:
                log(f"[pick_trunclen]   {name}: {cap} bp (no read shorter than {floor} bp)")
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

    # Change (a), continued: say what each quality cut is based on. A cut equal to the
    # read-length limit means the median Q1 never fell below q_threshold within it, so
    # quality decided nothing and the filtering is left to max_ee and trunc_q; up to
    # v2.6.1 the log showed this only as two lines with the same numbers.
    for label, primary, ceiling in (("R1", primary_r1, max_r1), ("R2", primary_r2, max_r2)):
        if primary >= ceiling:
            log(f"[pick_trunclen] {label} quality cut {primary} bp is the read-length limit: "
                f"the median Q1 never fell below {q_threshold:g} within it, so quality did "
                f"not set this cut")
        else:
            log(f"[pick_trunclen] {label} quality cut {primary} bp is set by quality: the "
                f"median Q1 falls below {q_threshold:g} at position {primary + 1} "
                f"(read-length limit {ceiling} bp)")

    if primary_r1 < 50 or primary_r2 < 50:
        log(f"[pick_trunclen] WARNING: very short primary cuts (R1={primary_r1}, R2={primary_r2}); inspect falco reports.")

    # ── ITS: overlap constraint not applicable — truncLen is overridden to c(0,0) ─
    # Skip resolve_expected_length entirely (not used for ITS) and write the
    # primary cuts so dada_filter can read them before overriding to (0, 0).
    if amp_type == "ITS":
        if method == "simulate":
            # No single cut is applied to ITS, so there is nothing to simulate; the
            # quality cuts are written for the record, as under method: quality.
            log("[pick_trunclen] method=simulate does not apply to ITS (truncLen is "
                "overridden to 0); writing the quality cuts for the record")
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

    # ── Non-ITS (16S/18S/gyrB/rpoB/COI): check the overlap of the quality cuts ──
    min_total = exp_len + min_overlap
    slack = (primary_r1 + primary_r2) - min_total
    log(f"[pick_trunclen] Overlap: r1+r2={primary_r1 + primary_r2}, need>={min_total}, slack={slack}")

    # The arguments every stop text below shares (overlap_stop_text above): the
    # expected length and where it came from, the overlap asked for, and how far the
    # reads reach (the ceilings max_r1/max_r2 left by the read-coverage cap).
    stop_args = dict(exp_len=exp_len, exp_len_source=exp_len_source,
                     probe_json_path=probe_json_path, probe_stat_key=probe_stat_key,
                     min_overlap=min_overlap, max_r1=max_r1, max_r2=max_r2,
                     min_read_coverage_pct=min_read_coverage_pct)

    # The quality rule's final cuts, after resolve_policy (resolve_overlap above).
    quality = resolve_overlap(resolve_policy, primary_r1, primary_r2, max_r1, max_r2,
                              agg_r1, agg_r2, min_total, q_threshold, q_floor,
                              stop_args, log)

    # ── method: quality — write the quality rule's cuts, or stop as before ──
    if method == "quality":
        if quality["exit_code"] != 0:
            log(quality["stop_text"])
            return quality["exit_code"]
        cut_r1, cut_r2 = quality["cut_r1"], quality["cut_r2"]
        final_slack = (cut_r1 + cut_r2) - min_total
        result = {
            "r1": int(cut_r1),
            "r2": int(cut_r2),
            "mode": "auto",
            "amplicon_type": amp_type,
            "q_threshold_requested": q_threshold,
            "q_threshold_used": quality["q_used"],
            "expected_length": exp_len,
            "expected_length_source": exp_len_source,
            "min_overlap": min_overlap,
            "overlap_slack": final_slack,
            "resolved_via": quality["resolved_via"],
            "aggregation": "median_across_samples_per_bin",
        }
        write_json(out_json, result)
        log(f"[pick_trunclen] DONE. truncLen=(R1={cut_r1}, R2={cut_r2}) via {quality['resolved_via']}; final slack={final_slack}")
        return 0

    # ── method: simulate — replay filterAndTrim on sampled read pairs ──
    # For every sample: take SIMULATE_PAIRS_PER_SAMPLE pairs spread evenly over its
    # files, and find for each read the longest cut at which it would still pass the
    # filter (simulate_sample). Pool all samples into one table of pairs kept at every
    # (R1, R2) cut (survival_table) and take the best cut pair that leaves
    # expected length + min_overlap + simulate_extra_overlap bases (best_simulated_cuts).
    start_time = time.time()
    per_sample: list[tuple[str, list[tuple[int, int]]]] = []
    for r1_path, r2_path, falco_path in zip(stripped_r1, stripped_r2, falco_r1):
        n_reads = parse_total_sequences(falco_path)
        pairs = simulate_sample(r1_path, r2_path, n_reads, SIMULATE_PAIRS_PER_SAMPLE,
                                trunc_q, max_ee_r1, max_ee_r2, max_len)
        per_sample.append((falco_path.parent.name.replace("_R1_stripped", ""), pairs))
    pooled_pairs = [pair for _name, pairs in per_sample for pair in pairs]
    if not pooled_pairs:
        log("[pick_trunclen] ERROR: method=simulate found no read pairs in the primer-trimmed files.")
        return 2
    longest = max(max(m1, m2) for m1, m2 in pooled_pairs)
    kept = survival_table(pooled_pairs, longest)
    n_sampled = len(pooled_pairs)
    log(f"[pick_trunclen] Simulated filterAndTrim (truncQ {trunc_q}, maxN 0, maxEE "
        f"{max_ee_r1:g}/{max_ee_r2:g}) on {n_sampled:,} read pairs from {len(per_sample)} "
        f"samples (up to {SIMULATE_PAIRS_PER_SAMPLE:,} per sample, spread over each whole "
        f"file) in {time.time() - start_time:.0f} s")

    def pct(n: int) -> str:
        return f"{100.0 * n / n_sampled:.2f}%"

    # How the margin changes the answer, for the reader who wants to set it: the best
    # cut pair and the share of pairs it keeps at a few margins around the one used.
    tolerance = n_sampled * SIMULATE_TOLERANCE_PCT / 100.0
    log(f"[pick_trunclen] Cut pairs within {SIMULATE_TOLERANCE_PCT:g} points of the most pairs "
        f"kept count as equally good; of those, the longest R1 + R2 is taken")
    for margin in sorted({0, 10, 20, 30, 40, extra_overlap}):
        option = best_simulated_cuts(kept, min_total + margin, min_len, tolerance)
        mark = " (used)" if margin == extra_overlap else ""
        if option is None:
            log(f"[pick_trunclen]   extra overlap {margin} bp: no cut pair reaches {min_total + margin} bp{mark}")
        else:
            log(f"[pick_trunclen]   extra overlap {margin} bp: best R1={option[0]}, R2={option[1]} "
                f"keeps {pct(option[2])} of the sampled pairs{mark}")

    best = best_simulated_cuts(kept, min_total + extra_overlap, min_len, tolerance)
    margin_used = extra_overlap
    if best is None or best[2] == 0:
        # The reads cannot carry the margin: fall back to the overlap rule alone, so the
        # pairs still merge, and say so.
        best = best_simulated_cuts(kept, min_total, min_len, tolerance)
        margin_used = 0
        if best is not None and best[2] > 0:
            log(f"[pick_trunclen] WARNING: no cut pair keeps any read pair with "
                f"{extra_overlap} bp of extra overlap; using the overlap rule alone "
                f"({min_total} bp).")
    if best is None or best[2] == 0:
        longest_r1 = max(m1 for m1, _m2 in pooled_pairs)
        longest_r2 = max(m2 for _m1, m2 in pooled_pairs)
        log(overlap_stop_text(
            f"R1 and R2 cannot be cut so that the pairs still merge. In the simulation "
            f"(method: simulate) the longest cuts at which any sampled read still passes "
            f"the filter are R1 {longest_r1} and R2 {longest_r2} bp, and no cut pair that "
            f"reaches the {min_total} bp merging needs keeps a single read pair.",
            **stop_args))
        return 3
    cut_r1, cut_r2, n_kept = best

    # Always log the quality rule's answer next to the simulated one, with the share of
    # the same sampled pairs it would keep, so the two can be compared on every run.
    if quality["exit_code"] == 0:
        q_r1, q_r2 = quality["cut_r1"], quality["cut_r2"]
        log(f"[pick_trunclen] Quality rule for comparison: R1={q_r1}, R2={q_r2} via "
            f"{quality['resolved_via']}, keeps {pct(kept_at(kept, q_r1, q_r2))} of the "
            f"sampled pairs")
    else:
        first_line = quality["stop_text"].splitlines()[0].replace("[pick_trunclen] ERROR: ", "")
        log(f"[pick_trunclen] Quality rule for comparison: would have stopped "
            f"(exit code {quality['exit_code']}): {first_line}")
    log(f"[pick_trunclen] Simulated cut: R1={cut_r1}, R2={cut_r2}, keeps {pct(n_kept)} of "
        f"the sampled pairs (R1 + R2 = {cut_r1 + cut_r2} bp, need >= {min_total} + "
        f"{margin_used} bp extra overlap)")

    # Per sample, the share of its sampled pairs kept at each cut: shows whether one
    # sample pays for the choice (for example one with its N bases at an early cycle).
    for name, pairs in per_sample:
        if not pairs:
            continue
        n_sim = sum(1 for m1, m2 in pairs if m1 >= cut_r1 and m2 >= cut_r2)
        line = (f"[pick_trunclen]   {name}: {len(pairs):,} pairs sampled, "
                f"{100.0 * n_sim / len(pairs):.1f}% kept at the simulated cut")
        if quality["exit_code"] == 0:
            n_q = sum(1 for m1, m2 in pairs if m1 >= q_r1 and m2 >= q_r2)
            line += f", {100.0 * n_q / len(pairs):.1f}% at the quality rule's cut"
        log(line)

    # Same keys as under method: quality, so dada_filter, dada_length_filter and the
    # report read the file unchanged. No quality threshold set this cut, hence null.
    final_slack = (cut_r1 + cut_r2) - min_total
    result = {
        "r1": int(cut_r1),
        "r2": int(cut_r2),
        "mode": "auto",
        "amplicon_type": amp_type,
        "q_threshold_requested": q_threshold,
        "q_threshold_used": None,
        "expected_length": exp_len,
        "expected_length_source": exp_len_source,
        "min_overlap": min_overlap,
        "overlap_slack": final_slack,
        "resolved_via": "simulate",
        "aggregation": "simulated_filterAndTrim_on_sampled_pairs",
    }
    write_json(out_json, result)
    log(f"[pick_trunclen] DONE. truncLen=(R1={cut_r1}, R2={cut_r2}) via simulate; final slack={final_slack}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
