# Amplicon length and truncation

Two length decisions shape an amplicon run, and MetaFlux measures both rather than
asking for them to be typed in:

- **`truncLen`** — where each read is cut before DADA2 sees it. Chosen from the
  per-base quality profile of the run itself.
- **the ASV length window** — which assembled ASVs are kept at the end. Chosen from
  the length distribution of the amplicon that the supplied primers actually produce
  against the marker's reference database.

The two are linked: the same expected amplicon length that keeps the truncation
mergeable also sizes the ASV window. This page walks through both, in the order the
pipeline computes them.

ITS behaves differently at almost every step below. Its exceptions are flagged as they
come up and collected in [ITS: the exception](#its-the-exception).

## What `truncLen` does, and why it is a compromise

`truncLen` is one length per read direction. In DADA2's `filterAndTrim` it does two
things at once:

- every read is cut to that length, and
- every read **shorter** than that length is discarded.

That single number is pulled in two opposite directions:

- **Trim shorter** and the low-quality 3′ tail goes away. This is what DADA2's error
  model wants, because it learns a substitution rate for each of the twelve nucleotide
  transitions as a function of quality score, and a noisy tail both raises those rates
  and generates spurious uniques.
- **Trim too short** and R1 and R2 stop overlapping in the middle of the amplicon.
  `mergePairs` then cannot join them, and a pair that does not merge is discarded as a
  whole.

So the cut has to remove the bad tail but stop while the forward and reverse reads
still reach each other. MetaFlux picks the first part from the data and then enforces
the second part as a hard constraint.

## Choosing the cut from the quality profile

The rule is `pick_trunclen` (`workflow/scripts/amplicon/50a_pick_trunclen.py`). It runs
after primer trimming and after per-stage QC, and it writes a single small file,
`stats/trunclen.json`, which every per-sample `dada_filter` job then reads.

**What it takes in.** Falco produces a FastQC-compatible `fastqc_data.txt` for each
sample and direction at each QC stage. `pick_trunclen` reads only the *stripped* stage —
the reads as they exist after Cutadapt has removed the primers, which is exactly what
DADA2 will be handed:

```text
stats/falco/{sample}_R1_stripped/fastqc_data.txt
stats/falco/{sample}_R2_stripped/fastqc_data.txt
```

**What it does with them.** From each file it parses the `>>Per base sequence quality`
block and keeps one column: the **Lower Quartile**, the 25th percentile of quality at
that position, usually written Q1. Using Q1 rather than the mean is deliberate — the
mean stays flattering long after the worst quarter of reads has gone bad, and it is the
bad quarter that generates spurious variants.

Falco reports quality for positions 1 to 9 individually and for groups of positions
beyond that. The group width depends on the longest read in the file: 5 positions
(`10-14`, `15-19`, …) when the longest read is 140 to 334 bp, as in 2 × 250 and 2 × 300
runs after primer trimming. Each group reports the mean of its per-position values. The
picker expands each group back to one value per position, so samples with different
group layouts (different read lengths) can be compared position by position.

**Aggregating across samples.** At each cycle position the per-sample Q1 values are
combined by taking the **median across samples**. One truncation is chosen for the whole
run, so one sample with an unusually poor tail cannot drag the cut inwards for
everything else.

**Two gates then decide how far the cut is even allowed to go.** They answer different
questions and both are needed.

*Within a sample — how far do most of its reads actually reach?* Falco reports quality
for every position up to the **longest** read in a file, however few reads get that
far. This matters because `filterAndTrim` **discards** any read shorter than `truncLen`
instead of padding it: a cut placed where only a few reads reach discards almost the
whole library. So before quality is
consulted at all, each sample's quality table is cut off at the length that
`amplicon.trunc_len.min_read_coverage_pct` percent of its reads still reach — 95% by
default, read from Falco's "Sequence Length Distribution" module.

*Across samples — does every sample reach this position?* A position is then only
considered if **all** samples have data there, so the run-wide ceiling is the smallest of
the per-sample ceilings and no sample is cut beyond its own range. This gate is fixed at
100% in the code and is not exposed as a config key.

The next section explains why the read-coverage gate is applied first.

### Why the read-coverage gate exists

If this part of truncation goes wrong, the run loses almost all of its reads. Three facts
combine to cause it.

**One: reads are not all the same length once the primer is off.** A library that came
off the sequencer at a uniform 301 bp does not stay uniform. The primer does not sit at
exactly the same offset in every read, and 3′ quality trimming removes different amounts
from different reads. What comes out of primer trimming is ragged: most reads clustered
tightly around one length, with a thin tail running longer.

**Two: Falco reports quality for every position that any read reached.** If a single
read out of two hundred thousand happens to be 289 bp long, Falco's per-base quality
table has a row for position 289. That table says nothing about *how many* reads got
that far — it is a quality profile, not a census.

**Three: DADA2 discards short reads rather than padding them.** `filterAndTrim` with
`truncLen = N` cuts every read at position N **and throws away outright any read that
never reached N**.

Together, these produce the failure. On a clean run, per-base quality
never drops below `q_threshold`, so the picker has no quality signal to cut on and simply
runs to the end of the table. The end of the table was written by the longest read in the
file. Fact three then deletes everything shorter.

!!! example "What that looks like on real data"

    `Surface_dolomite` (run DRR316764) of the 18S test set
    ([PRJDB12231](../about/test-datasets.md)), R1 after primer trimming: **211,974 reads**.
    Reading Falco's length distribution from the longest length downwards:

    | Read length | reads ending exactly here | reads reaching **at least** here |
    |---|--:|--:|
    | 289 bp | 1 | 1 |
    | 288 bp | 2 | 3 |
    | 285 bp | 16 | 23 |
    | 282 bp | 850 | 908 |
    | **281 bp** | **91,042** | 91,950 |
    | 276 bp | 25,516 | 122,513 |
    | 275 bp | 37,199 | 159,712 |
    | 272 bp | 135 | 162,376 |
    | **271 bp** | **49,512** | **211,888** |

    Read from the quality tables alone, the run-wide ceiling came out at **288 bp** —
    set by a different sample, since this one's own table runs to 289. Either way, only
    **three reads here are 288 bp or longer**, so a cut there discards 211,971 of
    211,974. Without this gate (MetaFlux v2.1.0 and earlier) the cuts landed at 288/239:
    only 9 of the 908,768 primer-trimmed read pairs in the four samples are long enough
    for them, and the run finished with 3 ASVs.

**What the gate does instead.** Before quality is consulted at all, each sample is asked
one question: *what is the longest position that at least `min_read_coverage_pct` of my
reads still reach?* Ninety-five percent of 211,974 is about 201,375. Walking up the
right-hand column above: at 272 bp only 162,376 reads (76.6%) qualify, which is not
enough; at 271 bp 211,888 do. So 271 becomes this sample's ceiling and the quality table
is truncated there. The picker now cannot choose anything beyond it, whatever the quality profile
says.

Two properties are worth noticing:

- **The share of reads kept is usually much higher than the setting.** Because read
  lengths cluster rather than spread evenly, the 95% request here retains 99.96% of the
  sample.
- **It is per sample, and then the strictest one wins.** The four 18S samples capped at
  275, 271, 276 and 271, so the run-wide R1 ceiling is 271. No sample is ever cut past
  its own limit.

**Reading it in the log.** `logs/pick_trunclen.log` names the sample that set each
ceiling and flags when the gate did real work:

```text
[pick_trunclen] R1 ceiling 271 bp, set by Surface_dolomite_R1_stripped (its longest read is 289 bp; 95.0% of its reads reach 271 bp)
[pick_trunclen] NOTE: R1 read-coverage cap pulled the ceiling in by 18 bp. Without it the cut could have landed above most reads.
```

The `NOTE` fires when the gap is 10 bp or more; read it rather than skip it. A large gap means the read lengths are ragged — normal after primer trimming,
and extreme for markers whose amplicons genuinely vary in length: on the ITS test set
([PRJNA862334](../about/test-datasets.md)) the gap is 89 bp (run SRR20648167: longest
read 298 bp, 95% of its reads reach 209 bp).

**When to change it.** Rarely. A higher percentage gives a shorter cut, because it asks
for a position that more reads reach; a lower percentage gives a longer cut.

| | Effect on the cut | Effect on reads |
|---|---|---|
| **Raise** it (e.g. 99) | **shorter** | fewer discarded |
| **Lower** it (e.g. 90) | **longer** | more discarded |

Measured on `Surface_dolomite` R1, so the shape is visible:

| `min_read_coverage_pct` | 50 | 80 | 90 | **95** | 98 | 99 | 100 |
|---|--:|--:|--:|--:|--:|--:|--:|
| resulting cap | 276 | 271 | 271 | **271** | 271 | 271 | **79** |

Between roughly 80 and 99 the value makes little difference, because read lengths
cluster rather than spread evenly. `100` does not switch the gate off: it requires every
read to reach the cut, so the ceiling drops to the shortest read in the file. On this
sample that is 79 bp, which then fails the overlap constraint and stops the run.

**Short off-target products are left out of the count.** If more than 5% of a sample's
reads belonged to a product much shorter than the amplicon, the 95% length would fall
within that product, and because the run-wide ceiling is set by the strictest sample, a
single such sample would cap both reads too short to overlap. In a 78-sample 16S soil
study sequenced in eight NextSeq 1000/2000 runs, 13 samples had 5 to 40% of their reads in
a product of 150 to 165 bp, with the mode at 158 bp. Counted over all reads, each of them
capped both reads at 158 bp, and the picker stopped at the overlap check.

The gate therefore counts only reads at least **half the expected amplicon length** long.
The reasoning:

- Cutadapt removes the other primer from a read that runs through the end of its
  product, so a primer-trimmed read shorter than its sequenced length is its whole
  product. A read under half the expected length comes from a product less than half as
  long as the target: a primer dimer or an unspecific product, not the natural length
  variation of the marker.
- Such a product could only give an ASV far below the
  [ASV length window](#the-asv-length-filter), which removes it later anyway.
- Target reads that do not run through (amplicon longer than the read) are always
  longer than the floor: R1 and R2 together must exceed the expected length for the pair
  to merge at all.

On the soil study (expected length 431 bp, floor 216 bp) the gate left out 0 to 40.7% of
each sample's reads, the 13 samples capped at 272 to 281 bp like the others, and the cuts
are 272/280. On the eight MiSeq test libraries no `truncLen` changed (the largest share
left out was 2.7%, in one COI library).

Two exceptions. When more than half of a sample's reads lie below the floor, they are the
sample rather than a minority product: that sample keeps the limit over all its reads,
and the log warns (check its primers and `amplicon.expected_length`). ITS has no single
expected length, so its gate still counts every read.

The log lists every sample's limit and the share of its reads left out:

```text
[pick_trunclen] R1 read-length limit per sample: the length 95% of its reads reach, not counting reads shorter than 216 bp (half the expected amplicon length: a short off-target product, not the target)
[pick_trunclen]   CABI-SG-cDNA-2.1_R1_stripped: 281 bp (40.13% of its reads are shorter than 216 bp and were left out)
```

!!! warning "There is no value that disables this gate"

    `100` is the strictest setting, not the most permissive. No value restores the
    behaviour of MetaFlux v2.1.0 and earlier, which had no gate; to reproduce it, use
    `mode: manual` and set the lengths yourself.

The gate runs under `mode: auto` only — in `manual` mode the lengths come straight from
`manual_r1` / `manual_r2` and no quality table is consulted at all, so the value has no
effect. (The key itself still has to be present and numeric even for a manual run, which
the shipped template guarantees.) On the 18S test set the default picks 271/239 unaided,
against 270/239 worked out by hand from a retention curve: 828,018 read pairs pass the
length filter, against 828,170 with the hand-picked lengths (measured in v2.2.0), a
difference of 0.02%.

**The cut itself.** Walking the aggregated profile from position 1, the cut is placed at
the last cycle still at or above `q_threshold` — that is, one position before the first
one that drops below it. If quality never drops below the threshold, the cut lands at
the last position that survived the coverage gate. A warning is logged if either cut
comes out below 50 bp, which usually means the falco reports are worth a look.

!!! warning "Binned quality scores"

    The default `q_threshold` was set on MiSeq data, which carry 33 distinct Phred values
    (2 and every integer from 7 to 38). The other instruments in the table bin quality
    scores. The NovaSeq 6000, NextSeq 1000/2000, NovaSeq X and MiSeq i100 use three bins plus
    a no-call value. The values depend on the
    instrument, its control software and, on the NextSeq 1000/2000, the reagent kit:

    | Instrument | Q values in the FASTQ |
    |:---|:---|
    | MiSeq (original) | 2 and every integer from 7 to 38: continuous |
    | HiSeq 2000/2500 | binned by default since control software HCS 2.2.38, for example 6, 15, 22, 27, 33, 37, 40 |
    | NextSeq 500/550 | 2, 14, 21, 27, 32, 36 |
    | NovaSeq 6000 | 2, 11, 25, 37 (12, 23, 37 on early software) |
    | NextSeq 1000/2000 | 12, 26, 34 (standard 50-300-cycle kits); 9, 20, 34 (standard 600-cycle kit); 12, 24, 40 (XLEAP-SBS) |
    | NovaSeq X / X Plus | 12, 20, 37 on control software before v1.2; 2, 12, 24, 40 on v1.2; 2, 9, 24, 40 (or 41) on later versions |
    | MiSeq i100 | 2, 12, 24, 38, or 2, 9, 23, 38 on later control software |

    Values were taken from Illumina's documentation, except for the MiSeq, the NextSeq
    500/550 (four runs of one NextSeq 500) and the current NovaSeq 6000, which were read from
    FASTQ files.

    The same values decide how DADA2 learns its error model: with `error_model: auto`,
    MetaFlux counts the distinct values in the filtered reads and takes DADA2's binned fit
    for 3 to 10 values above Q2 with gaps between them (no table above has two neighbouring
    values), anchored on the values found in the reads, never on this
    table (see [DADA2 error model](error-model.md#continuous-and-binned-quality-scores),
    and the setting `learn_errors.error_model` in
    [the configuration reference](../reference/configuration.md#dada2-tier-1-operational-parameters)).

    On binned data, the lower quartile of a single position in one sample can only take one
    of these values. It drops to a lower value only when at least a quarter of the bases at
    that position are at or below that value. Marginal calls are the bases in the lowest bin above the no-call
    value. Before binning, they are Q17 or lower on the NextSeq 1000/2000, NovaSeq X and MiSeq
    i100, and below about Q15 on the NovaSeq 6000. With the default `q_threshold: 20`, the
    cut is placed where at least a quarter of the bases are marginal calls. On MiSeq data, the
    same setting places the cut where at least a quarter of the bases are Q19 or lower. On
    the NextSeq 1000/2000, NovaSeq X and MiSeq i100, all calls from Q18 to Q29 receive the
    same value. Calls at Q18 and Q19 are therefore no longer counted as marginal, and no
    threshold distinguishes a tail at Q20 from one at Q29.

    **The cut can end up at the read-length limit.** Two 16S data sets with the
    standard 600-cycle chemistry (values 9, 20, 34) were tested. The first was the soil study
    above without its 13 samples with the short product (65 samples). The second was 31,514
    read pairs from four tiles of a NextSeq 2000 mock community run. In both, the median lower
    quartile across samples did not fall below 20 within the length reached by 95% of the
    reads. The cut was placed at that limit, and quality filtering was left to
    `max_ee` and `trunc_q`; see
    [automatic truncLen on very clean runs](../troubleshooting.md#amplicon-automatic-trunclen-can-drop-all-reads-on-very-clean-runs).
    On the MiSeq test data, quality set the R2 cut for 16S, 18S, ITS, gyrB and rpoB, and
    the R1 cut for 16S, rpoB and one of the three COI libraries. Every other cut, including
    R2 in all three COI libraries, was the read-length limit.

    **`q_threshold` selects a bin boundary.** A threshold above the low value and up to the
    middle value places the cut where a quarter of the bases are marginal calls. A
    threshold above the middle value places it where a quarter of the bases are below Q30.
    On a current NovaSeq 6000, all values from 12 to 25 give the first result; on NextSeq
    1000/2000 600-cycle data, all values from 10 to 20; on NovaSeq X data from control
    software before v1.2, all values from 13 to 20. On the last two, the default of 20 lies
    on the boundary. On MiSeq 16S, 18S and rpoB reads re-encoded with the 600-cycle values, a
    threshold of 21 shortened the quality cuts by 15 to 50 bp; on the soil study, it shortened
    the R2 cut by 1 bp. Falco's group means and the median across samples also produce
    intermediate values. On the soil study, thresholds from 21 to 34 moved the R2 cut from
    279 to 234 bp. A threshold above
    the top value fails at the first position: the quality cut is 0, and the final cuts
    come from `resolve_policy` alone.

    **What the log says.** `pick_trunclen` reads the first 10,000 reads of every
    primer-trimmed file and logs, per read direction, the quality values found, their
    kind (binned, continuous, by the same rule as the choice of the error function), what
    `q_threshold` means on them, and whether each cut came from quality or from the
    read-length limit. It warns when `q_threshold` equals one of the binned values, as the
    default does on 9, 20, 34 data, or lies above the highest. None of this changes a cut.
    On the NextSeq 2000 mock community run:

    ```text
    [pick_trunclen] R1 quality values (first 10,000 reads of each of the 1 primer-trimmed files): 3 values (9, 20, 34), binned
    [pick_trunclen] R1 q_threshold 20 on these values: the cut falls where at least a quarter of the bases carry 9 or lower (the values below 20: 9)
    [pick_trunclen] WARNING: R1 q_threshold 20 equals one of the quality values. A threshold one unit higher (21) would count 20 as low too and ask whether a quarter of the bases carry 20 or lower, a different question that can move the cut by tens of bases (the values above it: 34).
    [pick_trunclen] R1 quality cut 303 bp is the read-length limit: the median Q1 never fell below 20 within it, so quality did not set this cut
    ```

    The first 10,000 reads show every value a file uses, but not its quality: a NextSeq
    writes its files tile by tile. To judge the cut on the reads themselves, see
    [the simulate method](#the-simulate-method).

    **`relax_q` has little effect.** It lowers the threshold one unit at a time down to
    `q_floor` (15 by default). No current three-bin table has a value from 15 to 19.
    Lowering the threshold to the floor therefore moves the cut only past intermediate values
    from Falco's group means or the median across samples. If the overlap still fails at the
    floor, the run stops.
    On binned data, keep `resolve_policy: raise_trunc`.

### Config keys read by `pick_trunclen`

| Key | Default | What it does |
|---|---|---|
| `amplicon.trunc_len.mode` | `auto` | `auto` = derive the cuts from the quality profile; `manual` = use the two values below and skip the analysis entirely. |
| `amplicon.trunc_len.q_threshold` | `20` | The Q1 value the cut is placed at. Quality falling below this is what "the tail starts here" means. |
| `amplicon.trunc_len.q_floor` | `15` | The lowest threshold `relax_q` is allowed to fall back to. Ignored by the other policies. |
| `amplicon.trunc_len.resolve_policy` | `raise_trunc` | What to do when the quality-based cuts break the overlap constraint. See below. |
| `amplicon.trunc_len.manual_r1` | unset | R1 cut in bp, used only when `mode: manual`. |
| `amplicon.trunc_len.manual_r2` | unset | R2 cut in bp, used only when `mode: manual`. |
| `amplicon.trunc_len.method` | `quality` | `auto` only. `quality` = the rule described on this page; `simulate` = replay `filterAndTrim` on sampled read pairs, see [the simulate method](#the-simulate-method). |
| `amplicon.trunc_len.simulate_extra_overlap` | `20` | `simulate` only: bp of overlap required on top of `expected_length + min_overlap`. |
| `amplicon.dada2.trunc_q`, `max_ee`, `filter.min_len`, `filter.max_len` | see [DADA2](../reference/configuration.md) | Read by `simulate` only, which replays the filter with the run's own settings. |
| `amplicon.expected_length` | `auto` | Amplicon length used by the overlap constraint. See [where the expected length comes from](#where-the-expected-amplicon-length-comes-from). |
| `amplicon.min_overlap` | `12` | Overlap the two cuts must leave, in bp. Design constraint only. |

!!! note "`amplicon.min_overlap` is not the same setting as `dada2.merge.min_overlap`"

    `amplicon.min_overlap` (default `12`) is used **only** in the truncation arithmetic
    below — it decides where the reads are cut. `amplicon.dada2.merge.min_overlap`
    (also `12`, DADA2's own default) is the threshold `mergePairs` applies at runtime to
    each individual pair. They now hold the same value, so the truncation is sized for
    exactly the overlap the merge step will demand; raising the design constraint above
    the runtime threshold is what buys a margin.

In `manual` mode the quality analysis is skipped completely: the two configured values
are written straight to `stats/trunclen.json` and used as they are. No overlap check is
performed in manual mode — the responsibility moves to whoever set the numbers.

## Keeping the pair mergeable

For the two reads to still meet in the middle:

```text
truncLen_R1 + truncLen_R2  >=  expected_amplicon_length + min_overlap
```

The difference between the left and the right side is recorded in
`stats/trunclen.json` as `overlap_slack`. Positive slack means the cuts are comfortable;
negative slack means the quality-based cuts, taken at face value, would produce
unmergeable pairs.

When the slack is negative, `amplicon.trunc_len.resolve_policy` decides what happens:

| Policy | What it does | When it gives up |
|---|---|---|
| `raise_trunc` *(default)* | Extends the cuts past the quality drop until the deficit is covered, taking each additional base from whichever read has the **better aggregated Q1 at the position it would gain**. See below. | If the reads, cut at their ceilings, still do not reach the required length, the run stops. |
| `relax_q` | Lowers `q_threshold` by 1, re-cuts both directions at the new threshold, and repeats until the constraint holds. | If the threshold reaches `q_floor` without the constraint holding, the run stops. |
| `error` | Changes nothing and aborts immediately, reporting the deficit in bp. | Always — that is the point of it. |

The trade-off between the first two is worth stating plainly. `raise_trunc` keeps the
quality bar where it was and buys the missing bases by accepting a longer, slightly
worse tail. `relax_q` keeps the cut at a genuine quality transition but accepts a lower
bar for what counts as good. Neither is universally right; `raise_trunc` is the default
because losing a whole pair is not recoverable, whereas extra low-quality bases are
usually corrected by DADA2's denoising step before merging. `mergePairs` itself repairs
nothing: it works on already-denoised sequences, resolves the residual disagreements in
favour of the more abundant parent, and rejects outright any pair exceeding
`max_mismatch`.

!!! note "How `raise_trunc` decides which read gives up the bases"

    The deficit is walked one base at a time, and each base is taken from whichever read
    has the better aggregated Q1 at the position it would gain. Ties go to R1, and a
    read stops being a candidate once it reaches its own ceiling.

    Splitting by headroom instead — the distance from each cut to its ceiling — reads as
    the obvious approach and is the wrong one. Headroom is largest exactly where the cut
    landed earliest, which is where quality collapsed soonest. On typical Illumina that
    is R2, so a headroom split sends most of the extension into the worst tail on the
    run: the bases least worth keeping, which then cost reads at `maxEE` anyway. Taking
    the good bases first spends the same number of bases and keeps more of them useful.

    The log records the split and the mean quality of what each read gained. On the 16S
    test set ([PRJNA305879](../about/test-datasets.md)) the whole deficit goes to R1:

    ```text
    [pick_trunclen] Raised: R1 +30 (mean Q1 of added bases 21.8), R2 +0 (mean Q1 n/a); allocation is quality-weighted, not headroom-weighted
    ```

Each of the three stops prints the numbers behind the check and the settings that feed
it: how far the reads reach after primer removal, what merging needs, and where the
expected length came from. A 2×300 bp run whose expected length is the 582 bp of the
example above stops with:

```text
[pick_trunclen] ERROR: R1 and R2 cannot be cut so that the pairs still merge. Even at their full usable length (resolve_policy: raise_trunc) the reads reach 562 bp, 32 bp short of the 594 bp that merging needs.
  Read length available: after primer removal the reads reach at most R1 303 + R2 259 = 562 bp (the length that 95% of each sample's reads reach, amplicon.trunc_len.min_read_coverage_pct; the lines above name the sample that sets each).
  Merging needs: 582 bp (expected amplicon length) + 12 bp of overlap (amplicon.min_overlap) = 594 bp.
  Expected amplicon length: 582 bp, the p95 of the lengths of the 607,963 in-silico PCR products of these primers (amplicon.expected_length: auto; amplicon.probe_length_stat.16S: p95). The probe measured q1 403, median 422, q3 427, p95 582 and p99 640 bp.
  What to check:
  - If the amplicon between your primers is shorter than 582 bp, set amplicon.expected_length to its length in bp, or choose a lower statistic of the probe under amplicon.probe_length_stat.16S (q3: 427 bp, median: 422 bp).
  - If 582 bp is right, these reads are too short for R1 and R2 to overlap by 12 bp, and no truncation setting makes the pairs merge.
```

Under `relax_q` and `error` the first line says how far the quality cuts got instead, and
a further line names `raise_trunc` when the reads themselves would reach far enough.

Whichever route was taken is recorded as `resolved_via` in `stats/trunclen.json`
(`primary`, `raise_trunc`, `relax_q` or `manual`), along with `q_threshold_requested`
and `q_threshold_used` — those two differ only when `relax_q` fired.

## The simulate method

`amplicon.trunc_len.method: simulate` chooses the cuts from the reads themselves instead
of from a quality summary; this is the idea behind FIGARO. The default stays `quality`.

**When to use it.** On binned quality values the quality rule often cannot fire and
leaves the cut at the read-length limit. There `simulate` is the recommended method: on
the binned mock data it found the same true sequences and kept 6-7% more reads in ASVs.
On MiSeq data the default `quality` keeps the cuts of earlier versions.

**What it does.**

1. From each sample it takes 20,000 read pairs spread evenly over the **whole** file,
   never its start: a NextSeq writes its reads tile by tile, and in one soil R2 file the
   first 100,000 reads carried no N while 19.9% of all its reads did.
2. For each read it finds the longest `truncLen` at which the read still passes
   `filterAndTrim` with the run's own settings: the read is cut before its first base at
   `trunc_q` or below, must hold no N (maxN 0), and the expected errors of its first
   `truncLen` bases must stay within `max_ee` (`max_len`, when set, is applied first).
   Every one of these conditions only gets harder as the cut grows, so one number per
   read is enough: a pair passes at (R1, R2) exactly when both cuts are within its two
   numbers.
3. Over all sampled pairs it counts the pairs kept at every R1/R2 cut pair and takes the
   pair that keeps the most while `R1 + R2 >= expected_length + min_overlap +
   simulate_extra_overlap`. Of the cut pairs within 0.5 points of the most, it takes the
   one with the longest `R1 + R2`, then the longer R1.

Each sample gives the same number of pairs, so a deeply sequenced sample does not decide
for the others. The method reads the same quality values as the filter, so binned values
need no translation, and it sees what no per-position summary shows: whole tiles with an
N at a cycle that varies from tile to tile. ITS stays with the quality rule, because
`dada_filter` applies no fixed cut to ITS.

**Why the two extra rules.** The share of pairs kept is often almost flat over a wide
range of cuts, and the strict maximum then follows differences below the sampling noise
(about 0.1 to 0.3 points). On a MiSeq COI library with a 205 bp amplicon, R1 34 / R2 183
kept 99.59% of the pairs and 57 / 200 kept 99.28%; the strict maximum with the 20 bp
margin was a 42 bp R1. The rule is also pulled towards R2 because R2 has the larger
`max_ee` (5 against 2). Short, lopsided cuts leave most of the amplicon read by one read
only, while the overlap is where `mergePairs` checks R1 against R2 and rejects pairs that
disagree. `simulate_extra_overlap` (20 bp) and the 0.5-point tolerance both buy back
overlap at a small cost in reads. `0` keeps the most pairs.

**What it logs.** The share of sampled pairs kept by the best cut pair at extra overlaps
of 0, 10, 20, 30 and 40 bp, the quality rule's cut with the share it keeps, the
simulated cut, and both shares per sample:

```text
[pick_trunclen] Simulated filterAndTrim (truncQ 2, maxN 0, maxEE 2/5) on 1,545,764 read pairs from 78 samples (up to 20,000 per sample, spread over each whole file) in 247 s
[pick_trunclen] Quality rule for comparison: R1=272, R2=280 via primary, keeps 72.13% of the sampled pairs
[pick_trunclen] Simulated cut: R1=206, R2=264, keeps 77.19% of the sampled pairs (R1 + R2 = 470 bp, need >= 443 + 20 bp extra overlap)
```

The simulation predicts the filter closely: on the real runs below, the share of pairs
that `filterAndTrim` kept differed from the predicted share by at most
0.42 points.

**Measured on the test data.** The cuts of the quality rule and of `simulate`, with the
share of the sampled read pairs each keeps:

| Library | `quality` | kept | `simulate` | kept |
|---|--:|--:|--:|--:|
| 16S V5-V7, MiSeq (6 samples) | 244/154 | 61.8% | 216/204 | 64.9% |
| 18S V4, MiSeq (4) | 271/239 | 95.2% | 271/259 | 93.8% |
| ITS2, MiSeq (6) | 209/204 | – | not applied (ITS) | – |
| gyrB, MiSeq (5) | 150/194 | 95.1% | 137/185 | 95.8% |
| rpoB, MiSeq (5) | 254/189 | 89.8% | 236/197 | 91.7% |
| COI BF3/BR2, MiSeq (1) | 249/271 | 71.6% | 184/270 | 84.6% |
| COI Leray-XT, MiSeq (1) | 258/269 | 84.7% | 178/192 | 96.9% |
| COI fwhF2/fwhR2n, MiSeq (1) | 203/204 | 88.4% | 63/198 | 99.0% |
| rpoB, MiSeq re-encoded with 9/20/34 (5) | 264/189 | 80.6% | 235/197 | 88.0% |
| 16S V3-V4 mock, NextSeq 2000, 9/20/34 (1) | 303/259 | 84.2% | 236/253 | 89.8% |
| 16S V3-V4 soil, NextSeq 1000/2000, 9/20/34 (78) | 272/280 | 72.1% | 206/264 | 77.2% |

On every library but 18S, `simulate` keeps more pairs; on 18S the quality cut already
leaves only 4 bp of overlap to spare, and the 20 bp margin costs 1.5 points. The run time
is a few seconds per library; the 78 soil samples took about 4 minutes, nearly all of it
spent decompressing the reads.

**ASVs on mock communities.** Full runs up to the ASV table, with both methods and
otherwise identical settings (error function `auto`). An ASV identical to a true sequence
of the mock counts as a truth match; every other ASV counts as spurious.

| Mock (true sequences) | Method | Cuts | Reads in ASVs | ASVs | True sequences found | Reads in them | Spurious ASVs | Reads in them |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| Zymo D6331, 16S V3-V4, NextSeq 2000 (47) | `quality` | 303/259 | 180,514 | 33 | 18 | 93.38% | 15 | 6.62% |
| | `simulate` | 236/253 | 191,516 (+6.1%) | 35 | 18 | 93.40% | 17 | 6.60% |
| rpoB, MiSeq re-encoded 9/20/34 (19) | `quality` | 264/189 | 96,163 | 54 | 19 | 96.86% | 35 | 3.14% |
| | `simulate` | 235/197 | 103,255 (+7.4%) | 40 | 19 | 97.50% | 21 | 2.50% |
| rpoB, MiSeq (19) | `quality` | 254/189 | 105,593 | 47 | 19 | 97.42% | 28 | 2.58% |
| | `simulate` | 236/197 | 107,548 (+1.9%) | 44 | 19 | 97.48% | 25 | 2.52% |
| COI fwhF2/fwhR2n, MiSeq (374 barcodes) | `quality` | 203/204 | 1,199,935 | 890 | 362 | 98.06% | 528 | 1.94% |
| | `simulate` | 63/198 | 1,329,856 (+10.8%) | 850 | 365 | 98.09% | 485 | 1.91% |

How the truth was set: for Zymo D6331, the stretch between the primers of every 16S copy
with exact primer sites; for rpoB, one sequence per isolate, taken from the runs (the
publication gives no sequences), so only the spurious ASVs can be compared, not which
true sequences are found; for COI, the 374 reference barcodes of the mock plus the
extra bulk-sample barcodes, an ASV counting as true when it occurs exactly within one.

On all four mocks `simulate` kept more reads, found the same or more true sequences, put
the same or a larger share of the reads into them, and gave fewer spurious ASVs on three
(Zymo: two more, together 0.004% more of the reads). That includes the COI library with
the lopsided 63/198 cut. Not yet tested: a 16S MiSeq mock, an ITS run (not affected), and
data from a NovaSeq 6000 or an XLEAP-SBS instrument.

## Where the expected amplicon length comes from

`amplicon.expected_length` is what the overlap constraint measures against. It accepts
three forms:

| Form | Example | Value used | `expected_length_source` |
|---|---|---|---|
| Single integer, in bp | `338` | the integer itself | `manual:int` |
| `[min, max]` range | `[320, 360]` | the **maximum** — the constraint has to hold for the longest amplicons, or those are the ones that fail to merge | `manual:range_max` |
| `auto` | `auto` | one statistic of the measured probe distribution | `auto:<stat>` |

The integer form is validated: it must be positive and fall between 50 and 10000 bp, a
sanity range that catches a length typed in the wrong unit. The range form must have
exactly two positive numbers with `min <= max`.

`auto` is the interesting case. It makes MetaFlux measure the amplicon instead of being
told about it, and that measurement comes from the length probe.

## The amplicon length probe

The `amplicon_probe` rule (`workflow/scripts/amplicon/40a_amplicon_probe.py`) runs
**only** when `amplicon.expected_length: auto`. With a manual length there is nothing to
measure, so the rule never runs. Whether that also saves a download depends on the
marker: for ITS and 18S the probe substrate (UNITE UCHIME, SILVA-Eukaryotic) is used by
nothing else and is not fetched; for 16S, gyrB, rpoB and COI it is downloaded either way.
For gyrB, rpoB, COI and 16S on SILVA 138.2 the probe substrate is the taxonomy reference
itself; for 16S on SILVA 144 it is SILVA's own download, from which MetaFlux builds the
training set.

Its job is to answer one question: given *these* primers and *this* reference database,
how long is the amplicon, and how much does that length vary? It answers it in one of
two ways, chosen by the marker pack, not by the user.

### `probe_mode: pcr` — in-silico PCR

Used when the reference database holds full-length sequences (16S, 18S, rpoB, COI).
Cutadapt
is run twice against the reference FASTA:

1. The reverse primer file is reverse-complemented with `seqtk seq -r`, because in a
   reference sequence written 5′→3′ the reverse primer appears as its reverse
   complement at the 3′ end.
2. **Pass 1** — `-g file:<fwd>` with `--discard-untrimmed`: keep only reference
   sequences in which the forward primer is found *anywhere*, and cut off the primer
   together with everything upstream of it.
3. **Pass 2** — `-a file:<rev_rc>` with `--discard-untrimmed`: of those survivors, keep
   only the ones in which the reverse-complemented reverse primer is found, and cut it
   off together with everything downstream. A `--minimum-length 30` floor is applied
   here.

The primers are deliberately left unanchored — no `^` or `$` — and that is what makes the
in-silico PCR work: a full-length reference is cut down to the amplicon window wherever
the primer sites happen to sit in it. A V3–V4 or V4 forward primer binds several hundred
bases into a SILVA entry that is around 1450 bp long, so requiring the primers at the
sequence ends would match almost nothing, and what did match would come out near
full-length rather than as the ~400 bp distribution the probe actually measures.

What comes out is the set of in-silico amplicon bodies — the reference sequences reduced
to exactly the stretch the primer pair would amplify. Their lengths are the probe
distribution. The primer mismatch tolerance is `amplicon.cutadapt.max_error_rate`
(default `0.2`), the same value used for trimming the real reads.

If no reference sequence carries both primers, the rule fails with a message pointing at
primer orientation and the error tolerance — a genuine failure, since sizing a window
from zero amplicons is meaningless.

#### A primer that matches twice

A primer can match a reference twice: at its real site and at a weaker look-alike
stretch elsewhere in the gene. With a degenerate primer, cutadapt keeps the leftmost of
the two matches, even when the later one fits better. For the reverse primer that is the
site nearest the forward primer, the right one. For the forward primer it is the site
farther from the reverse primer, so the product starts too early and comes out too long.

In rpoB this is common. With the `Univ_rpoB_deg` primers on the FROGS RefSeq rpoB release
MetaFlux uses (`rpob_refseq_cc_20240707`), the forward primer also matches, with the 4
mismatches that `max_error_rate: 0.2` allows, a stretch
about 920 bp upstream of its real site in 7,289 of the 42,995 products (17%), mostly
Gammaproteobacteria and Actinomycetes. Most of those products are 1.0–2.0 kb long; they
put p95 at 1,332 bp, so a fully automatic run would check the read overlap, and size the
ASV length window, for an amplicon more than three times too long; with 2×300 bp reads
`pick_trunclen` stops there, since no pair of truncation lengths spans 1,332 bp. The
rpoB marker pack therefore sets `probe_exclude_second_site: true`, which adds a
**pass 3**: every amplicon body is searched again for the forward primer, at its full
length and with the same mismatch rate, and the bodies that still hold a forward-primer
site are left out of the statistics. The probe then gives p95 395 bp (median 386 bp).
Such a body is left out rather than cut back to its inner site, because the sequence
alone does not tell which of the two sites a PCR would use.

In 16S it depends on the forward primer, and the 16S marker pack sets the field too. The
17-mer 341F (`CCTACGGGNGGCWGCAG`, the V3–V4 primer of Illumina's 16S protocol) also
matches, with 2 mismatches, `CCTACGGGAGAAAGCAG` in V2 (*E. coli* positions 185–201) of
*Pseudomonas*, *Acinetobacter*, *Stenotrophomonas*, *Aeromonas*, *Shewanella* and their
relatives. On SILVA 144, 6.3% of the 341F/805R products start there and still hold the
real site: 583 bp long instead of 427 bp. p95 leaves only 5% of the products above it,
so it lands among them, at 582 bp, and a 2×300 bp run, whose reads reach about 562 bp
after primer removal, stops at `pick_trunclen`. With pass 3 the probe gives p95 429 bp. For 799F/1175R
(V5–V7) the share is 0.69% on SILVA 144 (p95 386 → 385 bp) and 0.61% on 138.2 (p95 384 bp
either way), for 515F/806R (V4) 0.15% (p95 255 bp either way).

18S does not set the field. The same thing happens there, but rarely: 0.5% of the 18S
products with the TAReuk V4 primers. Left out, they would move p95 by 7 bp.

### `probe_mode: direct` — lengths read straight off the reference

Used when the reference database already ships amplicon-length sequences (ITS, gyrB).
No cutadapt, no primers: the lengths are read directly from the reference FASTA, and the
reference is gzip-copied to the rule's amplicon output so the declared output file exists
in both modes.

For these two markers it is the only option that works. UNITE's UCHIME
release contains pre-extracted ITS1/ITS2 subregions, and the common ITS primers (ITS1F,
ITS4) bind in the flanking rRNA regions that those sequences no longer contain — an
in-silico PCR would find nothing to cut against. DD7RZ8 ships primer-trimmed in-silico
gyrB amplicons rather than full genes, with the same consequence.

### Probe wiring per marker

| Marker | `probe_mode` | Probe substrate | Reference tag in the cache filename | `probe_length_stat` key read |
|---|---|---|---|---|
| 16S | `pcr` | SILVA trainset (chloroplast, mitochondrial and eukaryote products, and products that hold a second forward-primer site, left out of the statistics) | `silva_v144_toGenus_excl-chloroplast-eukaryota-mitochondria_excl-second-site` (`silva_v138.2_…` with release 138.2) | `16S` |
| ITS | `direct` | UNITE UCHIME, pre-extracted ITS1 or ITS2 | `unite_uchime_ITS1` / `unite_uchime_ITS2` | none — ITS never resolves a stat |
| 18S | `pcr` | SILVA-Eukaryotic 18S v132 | `silva_euk_18s_v132` | `18S` |
| gyrB | `direct` | DD7RZ8 v6 | `gyrb_dd7rz8_v6` | `16S` (gyrB reuses that key) |
| rpoB | `pcr` | FROGS RefSeq rpoB (products that hold a second forward-primer site left out of the statistics) | `rpob_refseq_cc_20240707_excl-second-site` | `rpoB` |
| COI | `pcr` | MIDORI2 GB273 UNIQ NUC_SP, the converted SINTAX file | `midori2_gb273_uniq_nuc_sp_co1` | `COI` |

For 16S, the in-silico products of chloroplast, mitochondrial and eukaryote references
are left out of the statistics (marker pack field `probe_exclude_taxa`, written into the
cache filename). Organelle SSU is often much longer than bacterial 16S between the same
primers, and SILVA 144 holds eight times more of it than 138.2: for 799F/1175R the
p95 would be 390 bp with it and is 385 bp without (384 bp with 138.2 either way).
Since `expected_length` feeds the overlap constraint, those extra 5 bp can push
truncLen further into the quality-poor 3' end, where it costs reads, for amplicons the
taxon filter removes later anyway.

For 18S the probe substrate and the taxonomy reference are deliberately different
databases: SILVA-Eukaryotic is probed, PR2 classifies. See
[18S rRNA](markers/18S.md).

### What the probe writes

Both modes compute the same statistics and write the same JSON:

```text
amplicon_type, probe_mode, primer_hash, reference_tag, reference_fasta,
n_amplicons, min, q1, median, q3, p95, p99, max, mean, stdev,
excluded_taxa, n_excluded
```

`pcr` mode adds `cutadapt_max_error_rate` and `cutadapt_min_length`. `excluded_taxa`
names the reference lineages left out of the statistics (the pack's
`probe_exclude_taxa`: Chloroplast, Eukaryota and Mitochondria for 16S, none for the
other markers) and `n_excluded` how many sequences that left out; `n_amplicons` and the
statistics count the rest. With `probe_exclude_second_site: true` (16S and rpoB) the JSON also
holds `n_excluded_second_site`, the number of products pass 3 left out because they hold
a second forward-primer site.

`amplicon.probe_length_stat.<key>` selects which one of these becomes
`expected_length`. Accepted values are `min`, `q1`, `median`, `q3`, `p95`, `p99` and
`max`; the shipped config sets `p95` for `16S`, `18S`, `rpoB` and `COI`. A percentile
near the top of the distribution is the safe choice for the overlap constraint, because
the constraint has to hold for the longer amplicons in the pool — sizing it on the median
would leave the long tail unmergeable.

!!! warning "The key must exist for the marker being run"

    Every marker except ITS *reads* a `probe_length_stat` key when
    `expected_length: auto`. The key is *required* for every marker except ITS, whether
    or not `expected_length` is `auto` — a 16S, 18S, gyrB, rpoB or COI run with a
    manual integer length and no matching key still stops at start-up, with the exact
    key to add in the message. Note that gyrB reads `probe_length_stat.16S`, not a key of its
    own.

### Probe caching

Probe results are cached in the repository, not in the output directory, so they survive
reruns and are shared between output directories:

```text
refdb/cache/probe_{type}_{reference_tag}_{primer_hash}.json
refdb/cache/probe_{type}_{reference_tag}_{primer_hash}.amplicons.fa.gz
```

`primer_hash` is the first 12 hex characters of the SHA-256 of the forward primer file,
a `|` separator, and the reverse primer file. Because the hash is taken over the file
contents, editing a primer sequence changes the filename and the probe re-runs; running
the same primers again hits the same cache entry and the probe is skipped. The reference
tag is in the name for the same reason — a different database is a different measurement.

The `.amplicons.fa.gz` file holds the in-silico amplicons themselves (in `direct` mode, a
gzipped copy of the reference). It is kept for inspection: it is the direct answer to
"which references did my primers actually match, and how long were the products".

## The ASV length filter

The second length decision happens at the far end of the pipeline, in the
`dada_length_filter` rule (`workflow/scripts/amplicon/60h_dada_length_filter.py`). By
then DADA2 has produced ASVs and, for 16S and ITS, Metaxa2 or ITSx has optionally
trimmed each ASV down to the marker region. The filter runs **after** extraction, so the
ASV lengths it sees are directly comparable to the probe distribution — both describe
the same stretch of sequence.

The window is chosen by the first of these that applies:

| Priority | Condition | Window | `source` recorded |
|---|---|---|---|
| 1 | `length_filter.range` is set to `[min, max]` | exactly that range | `manual:config_range` |
| 2 | `length_filter.mode: auto` and a probe JSON exists | `[probe_q1 − window_margin, probe_p95 + window_margin]`; for ITS with extraction on, `[probe_min, probe_p95 + window_margin]` | `auto:probe_q1_p95+/-<margin>`; for ITS with extraction on, `auto:ITS_probe_min_p95+<margin>` |
| 3 | `length_filter.mode: auto`, no probe JSON, and `trunc_len.mode: auto` — the only mode that records `expected_length` in `stats/trunclen.json` | `expected_length × 0.85` to `expected_length × 1.15`, read from `stats/trunclen.json` | `auto:fallback_expected_length±15pct` |

`amplicon.length_filter.window_margin` defaults to `50` bp and is added on both sides
(for ITS with extraction on, at the top only; see below).
The window is deliberately asymmetric in its sources: **q1** at the bottom and **p95** at
the top, so it spans the bulk of the real amplicon distribution.

P95 rather than P99 at the top sets a tighter ceiling. The ITS probe is the case to watch:
it reads lengths straight off pre-extracted ITS1/ITS2 subregions, so its upper tail is real
length variation between fungal lineages, and a p95 ceiling will clip genuinely long ITS
ASVs, biased towards Basidiomycota. For ITS, check `stats/dada2/asv_length_hist.png`
against the window and set an explicit `length_filter.range` if long-ITS taxa matter.

For ITS with extraction on, the floor is the probe's `min` (38 bp for ITS2, 17 bp for
ITS1), without margin: each ASV reaches the filter as the ITS region ITSx cut out of it,
so a short region is a short ITS, not a fragment. `q1 − window_margin` would put the ITS2
floor at 125 bp, which removed three yeasts of the [ITS test set](../about/test-datasets.md),
64 % of its reads after chimera removal (see
[the ITS page](markers/its.md#region-extraction-with-itsx)). With
`extraction.enabled: false` the window is built as for the other markers.

!!! warning "A manual `expected_length` with a manual truncation leaves route 3 nothing to read"

    Route 3 takes `expected_length` out of `stats/trunclen.json`, and only
    `trunc_len.mode: auto` writes that field — in `manual` mode `pick_trunclen` records
    just `r1`, `r2`, `mode`, `amplicon_type` and `resolved_via`. A run that combines
    `trunc_len.mode: manual` with a manual `expected_length` and the shipped
    `length_filter.mode: auto` therefore has neither a probe JSON (route 2) nor an
    `expected_length` to fall back on, and it fails in `dada_length_filter` — after the
    whole DADA2 chain and any extraction have already finished. Unlike the ITS case
    below, nothing catches this combination at start-up. Set `length_filter.range`
    explicitly whenever `expected_length` is not `auto`, or keep `trunc_len.mode: auto`.

**What the filter produces.** ASVs outside the window are dropped from the FASTA and
from the count table, and both keyed forms of the table are rewritten:

| Output | Contents |
|---|---|
| `05.dada2/seqs_lenfilt.fasta` | the surviving ASV sequences — the input to taxonomy |
| `05.dada2/seqtab_lenfilt_head_names.txt` | counts, columns keyed by `ASV_#` — also read by taxonomy |
| `05.dada2/seqtab_lenfilt_head_seqs.txt` | the same counts, columns keyed by the full sequence; a terminal deliverable, not consumed downstream |
| `stats/dada2/asv_length_stats.json` | length distributions at three stages — pre-extraction, filter source, kept — plus the window and its source |
| `stats/dada2/asv_length_hist.png` | the filter-source and kept distributions as a histogram, with the window boundaries drawn in; the pre-extraction distribution is overlaid as well when extraction ran |

The stats JSON reports the pre-extraction distribution alongside the post-extraction one
even when extraction is enabled, which makes it easy to see how much Metaxa2 or ITSx
actually trimmed off each ASV. For 16S the two usually match closely, because a 16S
amplicon sits entirely inside the gene and Metaxa2 has almost nothing to trim; for ITS
they typically diverge, because ITSx strips off the flanking 5.8S/SSU/LSU rRNA fragments
that surround the true ITS region — a large gap here reflects real trimming, not a
problem with the filter.

If the window removes every ASV, the rule fails rather than writing an empty table —
the message says to widen the window or check the probe.

## ITS: the exception

ITS length varies genuinely and substantially between fungal taxa. A fixed truncation
would not be trimming noise, it would be discarding real short amplicons. So for ITS,
`dada_filter` overrides `truncLen` to `c(0, 0)` — no fixed cut at all.

What this means in practice:

- `pick_trunclen` **still runs** for ITS. It computes the quality-based cuts and writes
  them to `stats/trunclen.json` for reference and QC, then skips the overlap constraint
  entirely and records `expected_length: null` with
  `expected_length_source: not_applicable:ITS_truncLen_overridden`.
- The override is unconditional, so it also wins over `mode: manual`. `manual_r1` and
  `manual_r2` have no effect whatsoever for ITS.
- The 3′ tail is **not** left untrimmed. `truncQ` (default `2`) still trims every read
  adaptively at its first base at or below that quality, and `maxEE` (default `[2, 5]`
  for R1/R2) discards reads with too many expected errors. The lever for a poor ITS R2 is
  therefore a stricter `truncQ` or `maxEE`, never a fixed `truncLen`.
- `minLen` for ITS is auto-derived from the probe distribution, using the statistic named
  in `amplicon.dada2.filter.min_len_stat`, when that is set (it is `null` by default).
  The ITS probe measures extracted ITS1/ITS2 lengths, which are below per-read length
  and therefore a sensible per-read floor. When no probe JSON exists, ITS falls back to
  the fixed `amplicon.dada2.filter.min_len` (default `20`). Every other marker always
  uses that config value, because their probe measures the full amplicon — longer than a
  single read, so using it as a per-read floor would drop everything.

!!! warning "ITS with `length_filter.mode: auto` requires `expected_length: auto`"

    Because the ITS branch always writes `expected_length: null` to
    `stats/trunclen.json`, the length filter's fallback route (priority 3 above) has
    nothing to read. Without `expected_length: auto` there is no probe JSON either, so
    neither route works. MetaFlux catches this combination at start-up rather than after
    the full DADA2 run: either set `expected_length: auto`, or set
    `length_filter.mode: manual` with an explicit `length_filter.range`.

See [ITS](markers/its.md) for the rest of the ITS path.

## Reading the decision afterwards

`stats/trunclen.json` is small and worth opening after every run. In `auto` mode it
holds:

| Field | Meaning |
|---|---|
| `r1`, `r2` | the truncation lengths actually handed to `filterAndTrim` |
| `mode` | `auto` or `manual` |
| `amplicon_type` | the marker the run used |
| `q_threshold_requested` / `q_threshold_used` | differ only when `relax_q` fired; `q_threshold_used` is `null` under `method: simulate` |
| `expected_length` | the resolved amplicon length (`null` for ITS) |
| `expected_length_source` | `manual:int`, `manual:range_max`, `auto:<stat>`, or the ITS placeholder |
| `min_overlap` | the design overlap the constraint used |
| `overlap_slack` | bases of margin left over; `null` for ITS |
| `resolved_via` | `primary`, `raise_trunc`, `relax_q`, `simulate` or `manual` |
| `aggregation` | `median_across_samples_per_bin`, or `simulated_filterAndTrim_on_sampled_pairs` under `method: simulate` |

In `manual` mode only `r1`, `r2`, `mode`, `amplicon_type` and `resolved_via` are written,
since nothing else was computed.

The full reasoning — the quality values of the reads and what `q_threshold` means on
them, each sample's read-length limit and the share of its reads left out, per-direction
primary cuts and what set each of them, coverage ceilings, the sample driving each
ceiling, the overlap arithmetic, any policy that fired and, under `method: simulate`, the
simulated cut next to the quality rule's — is in `logs/pick_trunclen.log`.

## When `auto` goes wrong

The read-coverage gate guards against one failure mode by default; lowering
`min_read_coverage_pct` weakens it.

On a very clean run, per-base quality can stay above `q_threshold` right to the last
cycle. The picker then has no quality signal to cut on and takes the ceiling instead.
If that ceiling were the maximum read length, it would sit above where most reads end
after primer trimming — and since `filterAndTrim` discards every read shorter than
`truncLen`, the DADA2 step fails with `Error: No reads passed the filter`.

This is what `min_read_coverage_pct` prevents. The 18S test set
([PRJDB12231](../about/test-datasets.md)) is a real instance: its R1 quality never falls
below Q20, its longest reads reach 289 bp while 70% of R1 reads stop at exactly 281, and
only 9 of its 908,768 primer-trimmed read pairs are long enough for the 288/239 cut that
MetaFlux v2.1.0 and earlier chose. With the default 95% cap the ceiling comes in at 271,
`auto` returns 271/239 unaided, and 90.9% of the raw read pairs reach the final table.

It can still bite in two situations:

- **`min_read_coverage_pct` lowered.** A lower percentage permits a longer cut, which
  is the direction that brings this failure back. (Raising it shortens the cut instead,
  and `100` collapses the ceiling onto the shortest read in the file rather than
  disabling anything.)
- **Read lengths so ragged that even the 95% length is unrepresentative** — possible on
  heavily quality-trimmed input, or where a marker's amplicons vary enormously in length.

In either case the symptom is the same and so is the fix: read the post-trimming
"Sequence Length Distribution" from the `QC | primer-trimmed reads | {sample}_R{1,2}`
entries in MultiQC (or the Falco reports under `stats/falco/{sample}_R{1,2}_stripped/`),
switch to `mode: manual`, and set `manual_r1` / `manual_r2` a little below the length
most reads still reach — while keeping
`truncLen_R1 + truncLen_R2 >= amplicon_length + min_overlap` so the pairs still merge.
Set the two directions independently: usually only one of them is stuck on its ceiling,
and overriding the other throws away good sequence for nothing.
Re-running is cheap: trimming and QC are already done, so only DADA2 and the steps after
it repeat.

Full symptom and fix in [Troubleshooting](../troubleshooting.md).

## Related pages

| Page | Read it for |
|---|---|
| [Amplicon overview](overview.md) | the whole DADA2 path, step by step |
| [DADA2 error model](error-model.md) | how the quality values of the filtered reads choose the error function, and the report on the error models |
| [Markers](markers/index.md) | what each marker's reference databases are and how they differ |
| [Configuration](../reference/configuration.md) | every config key and its default |
| [Output files](../reference/output.md) | what else lands in `05.dada2/` and `stats/` |
| [Troubleshooting](../troubleshooting.md) | the `auto` truncation caveat and other practical failures |
