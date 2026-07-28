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
  `mergePairs` then cannot join them, and a pair that will not merge is dropped
  entirely — not shortened, dropped.

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

Falco reports quality in bins: positions 1 to 9 individually, then `10-14`, `15-19`, and
so on in fixed-width bins. Each bin is expanded back to one value per cycle position, so
that samples whose bin layout differs at the tail (different read lengths) can still be
compared position by position.

**Aggregating across samples.** At each cycle position the per-sample Q1 values are
combined by taking the **median across samples**. One truncation is chosen for the whole
run, so one sample with an unusually poor tail cannot drag the cut inwards for
everything else.

**Two gates then decide how far the cut is even allowed to go.** They answer different
questions and both are needed.

*Within a sample — how far do most of its reads actually reach?* Falco reports a quality
bin for every position up to the **longest** read in a file, however few reads get that
far. Trusting that would be a mistake, because `filterAndTrim` **discards** any read
shorter than `truncLen` instead of padding it: a cut placed where only a handful of reads
reach does not trim a few bases, it throws away the library. So before quality is
consulted at all, each sample's bins are cut off at the length that
`amplicon.trunc_len.min_read_coverage_pct` percent of its reads still reach — 95% by
default, read from Falco's "Sequence Length Distribution" module.

*Across samples — does every sample reach this position?* A position is then only
considered if **all** samples have data there, so the run-wide ceiling is the smallest of
the per-sample ceilings and no sample is cut beyond its own range. This gate is fixed at
100% in the code and is not exposed as a config key.

The order matters, and the next section is why.

### Why the read-coverage gate exists

This is the one part of truncation worth understanding properly, because when it goes
wrong it does not go slightly wrong — it costs the entire run. Three separate facts
combine into the problem, and none of them is obviously dangerous on its own.

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
never reached N**. There is no partial credit.

Put those together and the failure assembles itself. On a clean run, per-base quality
never drops below `q_threshold`, so the picker has no quality signal to cut on and simply
runs to the end of the table. The end of the table was written by the longest read in the
file. Fact three then deletes everything shorter.

!!! example "What that looks like on real data"

    `Surface_dolomite`, the 18S test set, R1 after primer trimming: **211,974 reads**.
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
    211,974. That is not hypothetical: before this gate existed the four-sample 18S run
    kept 9 read pairs out of 908,768 and finished with 3 ASVs.

**What the gate does instead.** Before quality is consulted at all, each sample is asked
one question: *what is the longest position that at least `min_read_coverage_pct` of my
reads still reach?* Ninety-five percent of 211,974 is about 201,375. Walking up the
right-hand column above: at 272 bp only 162,376 reads (76.6%) qualify, which is not
enough; at 271 bp 211,888 do. So 271 becomes this sample's ceiling and the quality table
is truncated there. The picker now cannot choose anything beyond it, whatever the quality profile
says.

Two properties are worth noticing:

- **You usually keep far more than you asked for.** The setting is a floor, not a target.
  Because read lengths cluster rather than spreading evenly, the 95% request here
  actually retains 99.96% of the sample.
- **It is per sample, and then the strictest one wins.** The four 18S samples capped at
  275, 271, 276 and 271, so the run-wide R1 ceiling is 271. No sample is ever cut past
  its own limit.

**Reading it in the log.** `logs/pick_trunclen.log` names the sample that set each
ceiling and flags when the gate did real work:

```text
[pick_trunclen] R1 ceiling 271 bp, set by Surface_dolomite_R1_stripped (its longest read is 289 bp; 95.0% of its reads reach 271 bp)
[pick_trunclen] NOTE: R1 read-coverage cap pulled the ceiling in by 18 bp. Without it the cut could have landed above most reads.
```

The `NOTE` fires when the gap is 10 bp or more, and is worth reading rather than
skipping. A large gap means the read lengths are ragged — normal after primer trimming,
and extreme for markers whose amplicons genuinely vary in length, where the ITS test set
shows an 89 bp gap.

**When to change it.** Rarely — and note which way it runs, because it is the opposite
of what "raise the percentage to be more permissive" would suggest. The value is a
*guarantee about reads*, not an allowance about length: asking for a higher percentage
asks for a position more reads reach, which is necessarily **further in**.

| | Effect on the cut | Effect on reads |
|---|---|---|
| **Raise** it (e.g. 99) | **shorter** | fewer discarded |
| **Lower** it (e.g. 90) | **longer** | more discarded |

Measured on `Surface_dolomite` R1, so the shape is visible:

| `min_read_coverage_pct` | 50 | 80 | 90 | **95** | 98 | 99 | 100 |
|---|--:|--:|--:|--:|--:|--:|--:|
| resulting cap | 276 | 271 | 271 | **271** | 271 | 271 | **79** |

Two things fall out of that. Between roughly 80 and 99 the value barely matters, because
read lengths cluster rather than spreading evenly — so there is little to gain by tuning
it. And **`100` is a cliff, not an off switch**: it demands that *every* read reach the
cut, so the ceiling collapses onto the single shortest read in the file. On this sample
that is 79 bp, which then fails the overlap constraint and stops the run.

!!! warning "There is no value that disables this gate"

    `100` does not restore the old behaviour — it is the most aggressive setting
    available, not the least. If you want the old behaviour for comparison, there is no
    config value for it; use `mode: manual` and set the lengths yourself.

The gate runs under `mode: auto` only — in `manual` mode the lengths come straight from
`manual_r1` / `manual_r2` and no quality table is consulted at all, so the value has no
effect. (The key itself still has to be present and numeric even for a manual run, which
the shipped template guarantees.) On the 18S set the default picks 271/239 unaided —
against 270/239 worked out by hand from a retention curve — the run keeps 828,018 reads
where the hand-picked pair kept 828,170, a difference of 0.02%.

**The cut itself.** Walking the aggregated profile from position 1, the cut is placed at
the last cycle still at or above `q_threshold` — that is, one position before the first
one that drops below it. If quality never drops below the threshold, the cut lands at
the last position that survived the coverage gate. A warning is logged if either cut
comes out below 50 bp, which usually means the falco reports are worth a look.

!!! warning "Quality-score binning changes how much these settings can do"

    The method above assumes quality scores form a continuous scale, which was true of
    four-channel SBS instruments — the original MiSeq and HiSeq — where roughly forty
    distinct Phred values appear. Current two-channel instruments do not report a
    continuum. Their real-time analysis software collapses every base call into three or
    four bins:

    | Instrument | Q values actually emitted |
    |---|---|
    | MiSeq, HiSeq 2500 (four-channel SBS) | ~40 distinct values — effectively continuous |
    | NovaSeq 6000 (RTA3) | 12, 23, 37 |
    | NextSeq 1000/2000 (RTA3) | 12, 26, 34 |
    | NovaSeq X / X Plus (RTA4) | 2, 9 (or 12), 24, 40 |
    | MiSeq i100 (RTA4) | 2, 9 (or 12), 23 (or 24), 38 |

    The cut is still placed correctly: the step from the high bin down to the low bin is
    a genuine quality collapse, and that is what the picker finds. What is lost is
    resolution, and it has two practical consequences on binned data.

    **`q_threshold` becomes a choice between bins, not a value.** On a NovaSeq 6000 every
    threshold from 13 to 23 selects exactly the same cut, because no base carries a score
    in between. Moving it from 20 to 18 or 22 changes nothing.

    **`relax_q` has almost nothing to work with.** It steps the threshold down one point
    at a time as far as `q_floor` (15 by default), and on a binned instrument that whole
    range sits inside a single bin — so the policy either changes nothing and then stops
    the run at the floor, or crosses a bin edge and changes the cut abruptly. On
    two-channel data, prefer `resolve_policy: raise_trunc`.

    A related effect: XLEAP-SBS chemistry puts the top bin at Q40, so the aggregated
    lower quartile may never fall below the threshold at all and the cut lands at the
    read length. That is the same situation described in
    [automatic truncLen on very clean runs](../troubleshooting.md).

### Config keys read by `pick_trunclen`

| Key | Default | What it does |
|---|---|---|
| `amplicon.trunc_len.mode` | `auto` | `auto` = derive the cuts from the quality profile; `manual` = use the two values below and skip the analysis entirely. |
| `amplicon.trunc_len.q_threshold` | `20` | The Q1 value the cut is placed at. Quality falling below this is what "the tail starts here" means. |
| `amplicon.trunc_len.q_floor` | `15` | The lowest threshold `relax_q` is allowed to fall back to. Ignored by the other policies. |
| `amplicon.trunc_len.resolve_policy` | `raise_trunc` | What to do when the quality-based cuts break the overlap constraint. See below. |
| `amplicon.trunc_len.manual_r1` | unset | R1 cut in bp, used only when `mode: manual`. |
| `amplicon.trunc_len.manual_r2` | unset | R2 cut in bp, used only when `mode: manual`. |
| `amplicon.expected_length` | `auto` | Amplicon length used by the overlap constraint. See [where the expected length comes from](#where-the-expected-amplicon-length-comes-from). |
| `amplicon.min_overlap` | `12` | Overlap the two cuts must leave, in bp. Design constraint only. |

!!! note "`amplicon.min_overlap` is not the same knob as `dada2.merge.min_overlap`"

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
| `raise_trunc` *(default)* | Extends the cuts past the quality drop until the deficit is covered, taking each additional base from whichever read has the **better aggregated Q1 at the position it would gain**. See below. | If the combined headroom of R1 and R2 is smaller than the deficit, the run stops with the deficit and both headrooms reported. |
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

    The log records the split and the mean quality of what each read gained:

    ```text
    [pick_trunclen] Raised: R1 +20 (mean Q1 of added bases 20.0), R2 +10 (mean Q1 12.0);
    allocation is quality-weighted, not headroom-weighted
    ```

Whichever route was taken is recorded as `resolved_via` in `stats/trunclen.json`
(`primary`, `raise_trunc`, `relax_q` or `manual`), along with `q_threshold_requested`
and `q_threshold_used` — those two differ only when `relax_q` fired.

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
nothing else and is not fetched; for 16S, gyrB and rpoB the probe substrate is the same
file as the taxonomy reference, so it is downloaded either way.

Its job is to answer one question: given *these* primers and *this* reference database,
how long is the amplicon, and how much does that length vary? It answers it in one of
two ways, chosen by the marker pack, not by the user.

### `probe_mode: pcr` — in-silico PCR

Used when the reference database holds full-length sequences (16S, 18S, rpoB). Cutadapt
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

### `probe_mode: direct` — lengths read straight off the reference

Used when the reference database already ships amplicon-length sequences (ITS, gyrB).
No cutadapt, no primers: the lengths are read directly from the reference FASTA, and the
reference is gzip-copied to the rule's amplicon output so the declared output file exists
in both modes.

This is not a shortcut, it is the only correct option for these two. UNITE's UCHIME
release contains pre-extracted ITS1/ITS2 subregions, and the common ITS primers (ITS1F,
ITS4) bind in the flanking rRNA regions that those sequences no longer contain — an
in-silico PCR would find nothing to cut against. DD7RZ8 ships primer-trimmed in-silico
gyrB amplicons rather than full genes, with the same consequence.

### Probe wiring per marker

| Marker | `probe_mode` | Probe substrate | Reference tag in the cache filename | `probe_length_stat` key read |
|---|---|---|---|---|
| 16S | `pcr` | SILVA trainset | `silva_v138.2_toGenus` | `16S` |
| ITS | `direct` | UNITE UCHIME, pre-extracted ITS1 or ITS2 | `unite_uchime_ITS1` / `unite_uchime_ITS2` | none — ITS never resolves a stat |
| 18S | `pcr` | SILVA-Eukaryotic 18S v132 | `silva_euk_18s_v132` | `18S` |
| gyrB | `direct` | DD7RZ8 v6 | `gyrb_dd7rz8_v6` | `16S` (gyrB reuses that key) |
| rpoB | `pcr` | FROGS RefSeq rpoB | `rpob_refseq_cc_20240707` | `rpoB` |

For 18S the probe substrate and the taxonomy reference are deliberately different
databases: SILVA-Eukaryotic is probed, PR2 classifies. See
[18S rRNA](markers/18S.md).

### What the probe writes

Both modes compute the same statistics and write the same JSON:

```text
amplicon_type, probe_mode, primer_hash, reference_tag, reference_fasta,
n_amplicons, min, q1, median, q3, p95, p99, max, mean, stdev
```

`pcr` mode adds `cutadapt_max_error_rate` and `cutadapt_min_length`.

`amplicon.probe_length_stat.<key>` selects which one of these becomes
`expected_length`. Accepted values are `min`, `q1`, `median`, `q3`, `p95`, `p99` and
`max`; the shipped config sets `p95` for `16S`, `18S` and `rpoB`. A percentile near the
top of the distribution is the safe choice for the overlap constraint, because the
constraint has to hold for the longer amplicons in the pool — sizing it on the median
would leave the long tail unmergeable.

!!! warning "The key must exist for the marker being run"

    Every marker except ITS *reads* a `probe_length_stat` key when
    `expected_length: auto`. The key is *required* for every marker except ITS at
    startup, whether or not `expected_length` is `auto` — a 16S, 18S, gyrB or rpoB run
    with a manual integer length and no matching key still stops at parse time, with
    the exact key to add in the message. Note that gyrB reads `probe_length_stat.16S`,
    not a key of its own.

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
`dada_length_filter` rule (`workflow/scripts/amplicon/60d_dada_length_filter.py`). By
then DADA2 has produced ASVs and, for 16S and ITS, Metaxa2 or ITSx has optionally
trimmed each ASV down to the marker region. The filter runs **after** extraction, so the
ASV lengths it sees are directly comparable to the probe distribution — both describe
the same stretch of sequence.

The window is chosen by the first of these that applies:

| Priority | Condition | Window | `source` recorded |
|---|---|---|---|
| 1 | `length_filter.range` is set to `[min, max]` | exactly that range | `manual:config_range` |
| 2 | `length_filter.mode: auto` and a probe JSON exists | `[probe_q1 − window_margin, probe_p95 + window_margin]` | `auto:probe_q1_p95+/-<margin>` |
| 3 | `length_filter.mode: auto`, no probe JSON, and `trunc_len.mode: auto` — the only mode that records `expected_length` in `stats/trunclen.json` | `expected_length × 0.85` to `expected_length × 1.15`, read from `stats/trunclen.json` | `auto:fallback_expected_length±15pct` |

`amplicon.length_filter.window_margin` defaults to `50` bp and is added on both sides.
The window is deliberately asymmetric in its sources: **q1** at the bottom and **p95** at
the top, so it spans the bulk of the real amplicon distribution.

P95 rather than P99 at the top sets a tighter ceiling. The ITS probe is the case to watch:
it reads lengths straight off pre-extracted ITS1/ITS2 subregions, so its upper tail is real
length variation between fungal lineages, and a p95 ceiling will clip genuinely long ITS
ASVs, biased towards Basidiomycota. For ITS, check `stats/dada2/asv_length_hist.png`
against the window and set an explicit `length_filter.range` if long-ITS taxa matter.

!!! warning "A manual `expected_length` with a manual truncation leaves route 3 nothing to read"

    Route 3 takes `expected_length` out of `stats/trunclen.json`, and only
    `trunc_len.mode: auto` writes that field — in `manual` mode `pick_trunclen` records
    just `r1`, `r2`, `mode`, `amplicon_type` and `resolved_via`. A run that combines
    `trunc_len.mode: manual` with a manual `expected_length` and the shipped
    `length_filter.mode: auto` therefore has neither a probe JSON (route 2) nor an
    `expected_length` to fall back on, and it fails in `dada_length_filter` — after the
    whole DADA2 chain and any extraction have already finished. Unlike the ITS case
    below, nothing catches this combination at parse time. Set `length_filter.range`
    explicitly whenever `expected_length` is not `auto`, or keep `trunc_len.mode: auto`.

**What the filter produces.** ASVs outside the window are dropped from the FASTA and
from the count table, and both keyed forms of the table are rewritten:

| Output | Contents |
|---|---|
| `5.dada2/seqs_lenfilt.fasta` | the surviving ASV sequences — the input to taxonomy |
| `5.dada2/seqtab_lenfilt_head_names.txt` | counts, columns keyed by `ASV_#` — also read by taxonomy |
| `5.dada2/seqtab_lenfilt_head_seqs.txt` | the same counts, columns keyed by the full sequence; a terminal deliverable, not consumed downstream |
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
    neither route works. MetaFlux catches this combination at startup rather than after
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
| `q_threshold_requested` / `q_threshold_used` | differ only when `relax_q` fired |
| `expected_length` | the resolved amplicon length (`null` for ITS) |
| `expected_length_source` | `manual:int`, `manual:range_max`, `auto:<stat>`, or the ITS placeholder |
| `min_overlap` | the design overlap the constraint used |
| `overlap_slack` | bases of margin left over; `null` for ITS |
| `resolved_via` | `primary`, `raise_trunc`, `relax_q` or `manual` |
| `aggregation` | `median_across_samples_per_bin` |

In `manual` mode only `r1`, `r2`, `mode`, `amplicon_type` and `resolved_via` are written,
since nothing else was computed.

The full reasoning — per-direction primary cuts, coverage ceilings, the sample driving
each ceiling, the overlap arithmetic, and any policy that fired — is in
`logs/pick_trunclen.log`.

## When `auto` goes wrong

One failure mode is worth knowing about, both because it is guarded against by default
and because the guard can be turned off.

On a very clean run, per-base quality can stay above `q_threshold` right to the last
cycle. The picker then has no quality signal to cut on and takes the ceiling instead.
If that ceiling were the maximum read length, it would sit above where most reads end
after primer trimming — and since `filterAndTrim` discards every read shorter than
`truncLen`, the DADA2 step fails with `Error: No reads passed the filter`.

This is what `min_read_coverage_pct` prevents. The 18S test set is a real instance: its
R1 quality never falls below Q20, its longest reads reach 289 bp while half stop at 281,
and a cut at 288 kept 9 read pairs out of 908,768. With the default 95% cap the ceiling
comes in at 271, `auto` returns 271/239 unaided, and the run keeps 90.9% of its reads.

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
| [Markers](markers/index.md) | what each marker's reference databases are and how they differ |
| [Configuration](../reference/configuration.md) | every config key and its default |
| [Output files](../reference/output.md) | what else lands in `5.dada2/` and `stats/` |
| [Troubleshooting](../troubleshooting.md) | the `auto` truncation caveat and other practical failures |
