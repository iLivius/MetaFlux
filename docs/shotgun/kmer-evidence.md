# Species-level k-mer evidence

!!! note "Read this after [Choosing confidence and threshold](confidence-and-threshold.md)"

    Bracken's depth-scaled read threshold, described on that page, has the largest effect
    on precision in shotgun mode; set it first. This gate catches one narrow kind of error
    that the threshold cannot detect. It is off by default; see
    [Why the default is `enabled: false`](#why-the-default-is-enabled-false).

Kraken2's read count tells you how many reads it chose to label with a taxon. It does
not tell you how much sequence *unique to that organism* actually turned up in the
sample — and those are different questions.

- A species that is genuinely present contributes many **different** species-specific
  sequences.
- A spurious call is usually a small set of sequences — conserved, repetitive, or
  shared with a relative that really is abundant — matched over and over.

In a read count the two look identical. Running Kraken2 with `--report-minimizer-data`
adds the number that separates them. MetaFlux always does.

## The two extra columns

The report becomes eight columns instead of six, with two inserted before the rank code:

| Column | Contents |
|--:|---|
| 4 | Number of minimizers in the read data associated with this taxon |
| 5 | An estimate of the number of **distinct** minimizers in the read data associated with this taxon |

Column 5 is an estimate, not an exact count — Kraken2 uses a HyperLogLog sketch rather
than storing every minimizer. Both columns are **clade-cumulative**, like the read
counts.

!!! warning "What "associated with this taxon" actually means"

    The columns do not count the minimizers carried by the reads assigned to the taxon.

    Kraken2 attributes each minimizer **occurrence** to whichever taxon that minimizer
    maps to **in the database** (its lowest common ancestor), pooled over every read in
    the sample, regardless of where each read finally got assigned.

    This is visible in any report: taxa with two or three assigned reads can carry
    hundreds of thousands of minimizer counts. In SRR39106521, one of the test samples
    described [below](#what-it-does-on-real-data), the genus row *Streptomyces* has
    **3 clade reads and 343,359 total minimizers**, about 500 times the 690 k-mer
    positions that three pairs of 149-base reads can hold. Counts that high cannot come
    from the assigned reads alone; they are expected when minimizers are attributed by
    LCA.

    So column 5 for a species means: **how many distinct database minimizers specific to
    that species were seen anywhere in this sample.** That is a statement about how much
    species-specific sequence is present — the minimizer counterpart of what KrakenUniq
    calls *unique k-mers*, so the published thresholds below are a fair starting point;
    how well they convert to minimizers is discussed
    [below](#choosing-a-threshold). It is **not** a genome-coverage measure, and
    describing it as "how much of the genome the reads covered" would be wrong.

!!! note "The flag itself changes no classification"

    `--report-minimizer-data` is purely additive reporting: with or without it, Kraken2
    writes a byte-for-byte identical per-read output and the same values in every shared
    report column.

## Reading the numbers

Two species from SRR39106521, with similar read counts:

| Species | Reads | Distinct minimizers | Total minimizers | Duplicity |
|---|--:|--:|--:|--:|
| *Pseudomonas aeruginosa* | 4,775 | 68,276 | 221,570 | 3.2 |
| *Pigmentiphaga* sp. CHJ604 | 4,129 | 929 | 226,726 | 244.1 |

The two also have almost the same total minimizer count. The first is supported by
68,276 distinct species-specific minimizers. The second is supported by 929, each one
seen about 240 times.

**Duplicity** (total ÷ distinct) is the most directly readable of the three numbers,
because it does not require knowing the sequencing depth. It also barely tracks
abundance: across the 2,418 species rows in the Kraken2 reports of the four test
samples, Pearson *r* against log reads is **+0.17** for log duplicity and **+0.71** for
the log distinct count. That makes it the better number to *eyeball*. It is nonetheless a
**worse filter** (see [What does not work](#what-does-not-work)), so MetaFlux reports it
and thresholds the distinct count.

## What MetaFlux writes

Every shotgun run produces these, per sample, whether or not the gate is switched on:

```text
02.classification/kraken2/{sample}_species_evidence.tsv   one row per taxon at the Bracken rank
02.classification/kraken2/{sample}_report_gated.txt       the report with failing taxa removed
stats/kmer_evidence/{sample}_evidence.png                 reads vs distinct minimizers, log-log
```

The gated report is written even when the gate is off, precisely so you can diff it
against the Kraken2 report and see what enabling the gate would do to *your* data before
trusting it.

The evidence table carries, per taxon: taxid, name, rank, clade and direct reads, total
and distinct minimizers, duplicity, minimizers per read, percent of classified reads,
and the verdict with its reason.

## The gate

```yaml
shotgun:
  kmer_evidence:
    enabled: false                     # report only; pipeline unchanged
    min_distinct_minimizers: 300       # the floor; see the unit conversion below
    min_reads: 0                       # extra read floor, applied together with the above
    never_prune_above_fraction: 0.01   # safety rail
```

### Where it acts

```text
kraken2  →  kmer_evidence  →  bracken  →  finalize_otu_table
```

The gate sits **between** Kraken2 and Bracken. It does not replace or skip Bracken — it
hands Bracken a cleaned report instead of the raw one. Because the correction happens
upstream, it flows all the way through: **`otu_table.biom` and `otu_table.tsv` become the
corrected tables.** There is no separate filtered file to reconcile.

Switching the gate off restores the previous behaviour exactly: Bracken reads Kraken2's
own report and nothing downstream changes.

### How a taxon is removed

Two details must be handled correctly; otherwise Bracken redistributes the reads of the
pruned report wrongly, without an error.

**Whole subtree, never a bare row.** Bracken rebuilds the taxonomy tree from the
*indentation* of the report, attaching each row to the last row one level above it.
Delete a species row and leave its strain (`S1`) children behind, and those children get
re-parented onto whichever species happens to precede them, whose read count is then
inflated without any warning. So when a taxon fails, everything indented beneath it goes
too.

**The reads are credited to the parent, not discarded.** A pruned taxon's whole clade total
is added to its parent's own-node count, so Bracken redistributes those reads among the
siblings that passed. The parent's clade total is untouched, so the report's arithmetic
still balances: for every node, `clade == direct + sum(children clade)`, and no column
other than a parent's direct count ever changes.

!!! warning "Two things this does not guarantee"

    **Not all reads land on close relatives.** The rationale for crediting the parent is
    that a spurious call usually borrows evidence from a neighbour. That holds for most
    of the reads that move, not all: part of them cross genus or family boundaries.

    **Some reads are still lost.** If *every* species under a genus fails, Bracken has
    nothing left at the target rank to redistribute to, and drops those reads, usually a
    small fraction of the sample. Crediting the parent is much better than deleting rows
    outright, but it is not conservation.

## Choosing a threshold

Three published sources give thresholds for this kind of filter, all expressed in
KrakenUniq's unique k-mers:

| Source | What it reports |
|---|---|
| [KrakenUniq](https://link.springer.com/article/10.1186/s13059-018-1568-0) (Breitwieser, Baker & Salzberg, 2018) | On a synthetic mix of reads from bacterial isolates (1–34 million read pairs), the best threshold (highest recall at a false discovery rate of 5 % or less) rose by about 2,000 unique k-mers per million reads sequenced: about 3,700 at 1 million, 21,000 at 10 million. For ten patient biopsies, mostly host DNA, a fixed 1,000 unique k-mers together with at least 10 reads removed the background calls and kept every confirmed pathogen. |
| [Oskolkov (2026)](https://doi.org/10.3389/fmicb.2026.1603339), *Front. Microbiol.* 17:1603339 | Simulated ancient metagenomes (fragments of about 42 bases, 0.1–0.7 million reads): about 1,000 unique k-mers gave the best F1, and a read-count filter added nothing. The best value rose with depth, by about 200 unique k-mers per 100,000 reads, a trend the paper calls relatively weak. Run on the same data, Kraken2 needed two to three times more distinct minimizers than KrakenUniq needed unique k-mers. |
| [aMeta](https://link.springer.com/article/10.1186/s13059-023-03083-9) (Pochon et al., 2023) | Defaults of 1,000 unique k-mers and 200 assigned reads, set from more than 1,200 ancient libraries; users are invited to try 500–1,500 k-mers depending on depth. The 2026 study above comes from the same group and reuses aMeta's simulated data. |

!!! warning "Published thresholds count k-mers; Kraken2 reports minimizers"

    **In plain terms.** Kraken2 reads a sequence as overlapping 35-base k-mers, but it
    does not look each one up. For every k-mer it picks one 31-base piece inside it, the
    *minimizer*, and stores and counts only that. Neighbouring k-mers usually pick the
    same piece, so a region of a genome gives about three times fewer distinct minimizers
    than distinct k-mers. The studies above counted unique k-mers and settled on about
    1,000. In a Kraken2 report the same amount of sequence shows up as about 333
    distinct minimizers. The default floor, **300**, is that figure rounded down, because
    the conversion is only approximate. It was not tuned on MetaFlux's own test data, and
    it is a low, cautious value rather than an exact equivalent. When Oskolkov (2026) ran
    both tools on the same simulated data, the best Kraken2 floor was two to three times
    *higher* than the best KrakenUniq one, where the conversion predicts three times
    lower, and MetaFlux's own benchmark finds useful floors in the thousands (see
    [Scaling the floor with depth](#scaling-the-floor-with-depth)).

    **The arithmetic**, for a database built with other settings. How many k-mers share
    a minimizer depends on the database's `k` and `l`. Each k-mer of length `k`
    spans `w = k − l + 1` candidate minimizer positions of length `l`, and the standard
    expected density of selected minimizers is `2 / (w + 1)`. For the usual nucleotide
    settings (`k = 35`, `l = 31`) that is `w = 5` and a density of `1/3` — roughly **one
    distinct minimizer per three k-mers**.

    So, counting sequence alone, a published *1,000 unique k-mers* corresponds to about
    **333** distinct minimizers here.

    Check `k` and `l` on the first line of your database's `inspect.txt` and recompute if
    they differ. (The `2/(w+1)` density is the expectation for a random minimizer scheme;
    Kraken2's spaced-seed masking shifts it slightly, so treat 333 as an approximate
    conversion, not an exact one.)

### What it does on real data

The numbers below come from the four metagenomes of BioProject
[PRJNA1477390](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1477390) (runs
SRR39106521–SRR39106524), microbial consortia enriched on naproxen from wastewater
treatment plants (see [Shotgun test dataset](../about/test-datasets.md#shotgun-test-dataset)).
They were processed with the MetaFlux v2.5.0 defaults (Kraken2 confidence 0.15, Bracken
`threshold: auto`) and PlusPF `k2_pluspf_20260626`, once with the gate off and once with
it on at the default floor of 300. The floor is `min_distinct_minimizers`: the number of
distinct minimizers a species needs in the Kraken2 report to stay in the table.

| | SRR39106521 | SRR39106522 | SRR39106523 | SRR39106524 |
|---|--:|--:|--:|--:|
| Species in `otu_table.tsv`, gate off | 99 | 80 | 116 | 120 |
| Species, gate on (floor 300) | 84 | 66 | 94 | 104 |
| Genera, gate off → on | 39 → 34 | 31 → 28 | 46 → 39 | 52 → 44 |
| Read pairs moved to another species | 4.9 % | 1.3 % | 5.8 % | 1.8 % |
| Change in total read pairs | −0.01 % | +0.01 % | −0.01 % | −0.09 % |
| Genomes detected by sylph in the same reads | 10 | 9 | 14 | 16 |
| 16S rRNA clusters at 97 % in the same reads, no reference | 14 | 13 | 11 | 15 |

- **The gate removes 14 to 22 species per sample and almost no reads.** The reads of a
  removed species go back to its parent taxon, and Bracken shares them among the
  relatives that passed.
- **The removed species are not only rare ones.** Bracken spreads the reads that Kraken2
  could only place at genus level or above over the species below, and it had given the
  removed species 10 to 27 times more reads than Kraken2 assigned to them. In every
  sample one of the 20 most abundant species is removed. The largest case is
  *Pseudomonas* sp. JS425 in SRR39106521: 372 read pairs from Kraken2, carrying 182
  distinct minimizers, raised by Bracken to 69,992 read pairs, the eighth most abundant
  species in the sample. Its read count passes the Bracken threshold; only the
  distinct-minimizer count shows how little species-specific sequence is behind it (more
  [below](#unchanged-f1-does-not-mean-an-unchanged-table)).
- **Fewer species is the right direction for these samples.** Two independent counts
  from the same reads put these cultures at about 10 to 16 members each. sylph, run by
  MetaFlux in `both` mode against `gtdb-r232-c200-dbv1` and
  `fungi-refseq-2025-10-11-c200-dbv1`, detects 9 to 16 genomes per sample. The 16S rRNA
  gene fragments in the reads were extracted with Metaxa2 2.2.3, and the fragments that
  cover the same region of the gene (the region with the most fragments, 284 to 504 per
  sample) were clustered at 97 % identity: 11 to 15 clusters per sample (13 to 18 at
  99 %; single-fragment clusters excluded). This count does not use a reference
  database. Both counts are lower bounds: sylph misses organisms without a closely
  related reference genome, and the 16S count detects only members represented by at
  least two fragments in that region, roughly those above 0.5 % of the 16S gene copies.
- **The gate does not bring the list down to that size.** Most of what remains is a
  database problem. The most abundant organism in SRR39106521 and SRR39106523,
  *Pusillimonas noertemannii*, has no genome in PlusPF. Its reads are reported under
  several *Achromobacter* species, which together hold more than half of each of these
  two tables; the largest of them have tens of thousands of distinct minimizers (see the
  warning on missing organisms in [Databases](databases.md#choosing-an-index)). A k-mer
  floor cannot correct that.

The same four samples at other floors: species left in `otu_table.tsv`, and the share of
the gate-off table held by the species that the floor removes.

| Floor | SRR39106521 | SRR39106522 | SRR39106523 | SRR39106524 | Read pairs held by removed species |
|---|--:|--:|--:|--:|--:|
| none (gate off) | 99 | 80 | 116 | 120 | 0 |
| 200 | 92 | 73 | 105 | 112 | 0.9–4.7 % |
| **300 (default)** | **84** | **66** | **94** | **104** | **1.3–5.8 %** |
| 1,000 | 55 | 44 | 61 | 80 | 8–21 % |
| 3,000 | 29 | 26 | 31 | 54 | 13–42 % |

A lower floor keeps the species with 200 to 300 distinct minimizers, 7 to 11 more per
sample. These samples have no known composition, so they cannot show whether those
species are real. A higher floor removes many more species but also moves a large share
of the table: at 3,000, 13 to 42 % of the read pairs end up under a different species. On
benchmark communities of known composition, a fixed floor of 10,000 already lowered some
samples' F1 by removing true species, and lower fixed floors did the same once shallower
sequencing was modelled (3,000 at 30 % of the depth, 1,000 at 10 %). A floor scaled to
the depth of the sample avoids this; see the next section.

### Scaling the floor with depth

A species that is really present leaves more distinct minimizers the deeper a sample is
sequenced, and so, up to a limit, does a spurious call that borrows sequence from an
abundant relative. A fixed floor is therefore lenient on deep samples and strict on
shallow ones. KrakenUniq and Oskolkov (2026) reported the same for unique k-mers (table
above), but neither gives a depth rule for Kraken2.

MetaFlux's benchmark was used to fit one. It covers 42 samples of known composition: one
real mock community (ZymoBIOMICS D6331, SRA run SRR38976715, BioProject PRJNA1474435)
and 41 simulated CAMI II samples from the marine, rhizosphere and strain-madness sets
(see [Benchmark data](../about/test-datasets.md#benchmark-data)), with 2.6–6.5 million
classified read pairs each, classified with three Kraken2 indices: `k2_pluspf_20260626`
(at confidence 0.15 and 0.20), `k2_gtdb_genome_reps_20250609` (GTDB r226 representative
genomes) and `k2_NCBI_reference_20251007` (one genome per NCBI species). A floor of **one
distinct minimizer per 1,000 classified read pairs** gave the best result that was safe
everywhere:

- species-level F1 rose by 0.066 on average (0.036 with PlusPF), against 0.0020 for the
  fixed default of 300 (0.0003 with PlusPF);
- it removed 2,288 of 13,115 false species calls and 6 of 13,419 true ones (none with
  PlusPF), and lowered no sample's F1;
- it stayed safe when shallower sequencing was modelled, down to 1 % of the reads, while
  every fixed floor tested, from 300 to 10,000, lowered some samples' F1 at some depth.

This is the same order of magnitude as KrakenUniq's slope (about 670 distinct minimizers
per million read pairs once converted), but the number itself comes from MetaFlux's
benchmark. As a starting point:

| Classified read pairs in the shallowest sample, host excluded | `min_distinct_minimizers` |
|--:|--:|
| 500,000 | 500 |
| 1 million | 1,000 |
| 2 million | 2,000 |
| 5 million | 5,000 |
| 6.5 million or more | 6,500 |

Each sample's classified count is in `stats/read_tracking.txt` (column `classified`) and
in its Bracken log (`[bracken] auto threshold: 5e-05 x N classified read pairs`).
MetaFlux applies one floor to every sample in a run, so size it on the shallowest sample;
the deeper samples then get a more lenient floor than they could take. Below about
500,000 classified pairs the scaled floor gains little. Use it even when it comes out
lower than the default of 300: in benchmark samples thinned to 26,000–195,000 classified
pairs, the fixed floor of 300 lowered some samples' F1 and the scaled floor lowered none.

Limits of this rule:

- **It rests on limited data.** Three quarters of the F1 gain came from the one real mock
  community; on the simulated CAMI samples it was 0.01–0.03. The shallower depths were
  modelled from the observed counts, not sequenced.
- **Do not go above 1.5 per 1,000 pairs.** At that rate rare true species start to be
  removed (92 of 13,419 true calls), and from 2 per 1,000 some samples' F1 falls even at
  full depth.
- **Do not go above 6,500 for deeper samples.** The benchmark stops at 6.5 million
  classified pairs. A floor that keeps growing with depth eventually exceeds what a real
  species with a small species-specific part of its genome can reach; in the model this
  began at about 10 million pairs.
- **Exclude host reads.** PlusPF contains the human genome and host removal is off by
  default, so host reads would raise the floor without adding any microbial evidence.
  Subtract the read pairs Kraken2 assigned to the host (for human, the *Homo sapiens*
  row of the report), or remove the host first ([Decontamination](decontamination.md#host-removal)).
- **Library layout.** The rule assumes paired 2×150 bp Illumina reads with little mate
  overlap and few duplicates. Shorter reads, single-end data, overlapping mates and
  duplicate-rich libraries give fewer distinct minimizers per read pair, so the same
  floor is stricter for them; use a lower value.
- **It changes the table more than 300 does.** On the test dataset above (1.7–3.5 million
  classified pairs), a floor of 1,673, sized on the shallowest sample with its host reads
  excluded, keeps 33–68 species per sample, against 66–104 at 300, and the species it
  removes hold 10–33 % of the read pairs, which move to relatives. That is still well
  above the 9–16 genomes sylph finds and the 11–15 16S clusters. Compare `otu_table.tsv` with the gate off and on
  before relying on it.

## Why the default is `enabled: false`

The gate is off by default because a correctly scaled
[Bracken threshold](confidence-and-threshold.md) already removes most of what it would
remove. It has been run on 46 samples: the D6331 mock, the four PRJNA1477390 metagenomes
and 41 CAMI II samples. In gate-on/gate-off runs of CAMI II samples at several depths
with PlusPF `k2_pluspf_20251015`, the gate changed species-level F1 by a mean of
**+0.0009** at the depth-scaled threshold, and by zero on the strain-madness and
rhizosphere samples at every depth tested. At the former fixed `-t 10` the mean was
**+0.0079**, about nine times more. So the gate mainly helps when the Bracken threshold is
not scaled to depth; with a scaled threshold it adds little.

!!! note "Which database build those two numbers come from"

    Both figures come from gate-on/gate-off runs on PlusPF `k2_pluspf_20251015` (15
    sample-depth combinations at the scaled threshold, 11 at the fixed one). The
    `k2_pluspf_20260626` benchmark runs had the gate off; removing from their Bracken
    output the species the gate would prune at the default floor of 300 changes mean
    species-level F1 by +0.0003.

### Unchanged F1 does not mean an unchanged table

The two filters ask different questions. The threshold asks *how many reads*; the gate
asks *how much distinct sequence*. A taxon can pass one and fail the other, and an
aggregate score of 250+ species barely notices when a single one flips.

*Pseudomonas* sp. JS425 in the four test samples, with the Bracken threshold scaled to
depth (`threshold_alpha` 5×10⁻⁵):

| Sample | Classified pairs | Bracken `-t` | Kraken2 pairs | Distinct minimizers | Passes `-t` | After Bracken, gate off | Gate |
|---|--:|--:|--:|--:|:--|--:|:--|
| SRR39106521 | 2,268,774 | 113 | 372 | 182 | yes | 69,992 | removed |
| SRR39106522 | 3,458,847 | 173 | 103 | 181 | no | – | removed |
| SRR39106523 | 1,673,810 | 84 | 301 | 181 | yes | 51,418 | removed |
| SRR39106524 | 2,177,429 | 109 | 123 | 181 | yes | 8,935 | removed |

In SRR39106521 the 372 read pairs pass the threshold, but they rest on 182 distinct
minimizers, each hit about 60 times, and Bracken raises the call to 3 % of the sample.
The threshold keeps JS425 in three samples and drops it from the fourth only because the
cut-off falls differently there; in a presence/absence plot that would look like a real
difference between samples. The gate removes it from all four.

These samples have no known composition, so they do not enter the F1 scores above. That
is why the benchmark shows almost no gain from the gate while the table changes.

**When to consider turning it on:** communities the database covers poorly, where reads
from organisms absent from the database are assigned to their relatives. The evidence
for this is limited. It rests on cases like JS425 in the wastewater test samples, which
have no known composition; on the CAMI II rhizosphere samples, the poorly covered
benchmark set with a known composition, the gate hardly changed F1. For gut, clinical
or other samples close to RefSeq, leave it off. Either way the gated report is written
even when the gate is off, so compare it with the Kraken2 report before deciding.

## Validation status

What the gate was tested against, and what those tests could and could not show.

### What has been checked

The gate was checked on one real mock community (Zymo D6331, below) and on
[CAMI II](https://www.nature.com/articles/s41592-022-01431-4)'s marine, rhizosphere and
strain-madness gold-standard sets (41 samples): in gate-on/gate-off runs with PlusPF
`k2_pluspf_20251015`, and by applying it to the Bracken calls of the
`k2_pluspf_20260626` runs. These results were scored with the benchmark's own scorer. On
the two `k2_pluspf_20260626` runs (Kraken2 confidence 0.15 and 0.20) that scorer was
checked against OPAL 1.0.13, the CAMI assessment tool: given the same input preparation,
OPAL's species-level counts of true positives, false positives and false negatives equal
the scorer's in every sample. The full methodology and the
numbers that matter for choosing a threshold are on the
[confidence and threshold](confidence-and-threshold.md) page. The result for this gate:
once the Bracken threshold is scaled to depth it adds at most 0.002 F1 (PlusPF
`k2_pluspf_20251015`), but it still removes calls like JS425 above, which a read-count
threshold cannot detect.

One real mock community: the [ZymoBIOMICS Gut Microbiome
Standard (D6331)](https://files.zymoresearch.com/protocols/_d6331_zymobiomics_gut_microbiome_standard.pdf),
SRA run [SRR38976715](https://www.ncbi.nlm.nih.gov/sra/SRR38976715) (BioProject
PRJNA1474435; see [Benchmark data](../about/test-datasets.md#benchmark-data)), with about
5.3 million classified read pairs.

The standard is sold as 21 strains. A species-level profiler can only ever see **17
species**, because five of the 21 are *Escherichia coli* strains at 2.8% each that collapse
into one *E. coli* row at 14%. That ceiling is arithmetic, not a reference-database gap.
Zymo lists the standard's *Fusobacterium* as *F. nucleatum* (strain 2/1/50A), but sylph
(GTDB r232), MetaPhlAn 4 and Kraken2 with GTDB r226 or PlusPF `k2_pluspf_20260626` place
most of these reads in *F. animalis*, a former subspecies of *F. nucleatum* now ranked as
a species. It is therefore scored as *F. animalis*, and calls of *F. nucleatum*,
*F. polymorphum* and *F. vincentii* count as false positives. *Homo sapiens*, from human
DNA in the sample, is not scored.

The run was processed twice through the pipeline with PlusPF `k2_pluspf_20251015`,
Kraken2 confidence 0.15 and Bracken at the former fixed `-t 10`: once with the gate off
and once with it on at a floor of 333. A floor of 300 gives the same species list on this
run. At the current default threshold (`-t 265` for this run), only 2 of the 34 false
positives in the table below have enough reads to pass (*F. nucleatum* and
*F. polymorphum*), as does the unscored *Homo sapiens*; the gate removes none of them,
and the threshold itself drops *Salmonella enterica*.

| | Gate off | Gate on |
|---|--:|--:|
| True members recovered | 15 / 17 | **15 / 17** |
| False-positive species | 34 | **22** |
| Precision | 0.306 | **0.405** |
| F1 | 0.455 | **0.556** |

- **No true member was removed.** *Salmonella enterica*, 0.01% of the standard and 134
  classified reads, has 3,547 distinct minimizers at a duplicity of 1.5, more than ten
  times the floor, and is kept. The gate depends on this difference: a rare species can
  still have many distinct minimizers.
- The two members never recovered (*Enterococcus faecalis* at 0.001%, *Clostridium
  perfringens* at 0.0001%) have one Kraken2 read each, below Bracken's `-t 10`, so neither
  reaches the table in either run. Both have species entries in the database, so they were
  missed for lack of reads.
- **Quantitative accuracy did not change.** With the human row left out, L1 distance to
  the theoretical composition was 47.15 without the gate and 47.14 with it, Shannon *H*
  2.145 and 2.144, and the Bray–Curtis dissimilarity between the two profiles 0.00009.
  The 82 taxa the gate removed at 333 held 426 of the 5.3 million reads (81 taxa and 423
  reads at 300). Any ordination or differential-abundance test would give the same answer
  from both tables.
- **Richness, on the other hand, improved substantially** — which is where the gate did its
  work on this dataset:

| Observed richness | Truth | Gate off | Gate on |
|---|--:|--:|--:|
| Species | 17 | 49 (+188%) | **37 (+118%)** |
| Genera | 17 | 25 (+47%) | **18 (+6%)** |

On this mock the gate changed which species are listed, not their abundances.
If you report how many taxa you detected, or work with presence/absence, that matters. If
you only ever compare relative abundances between samples, on data like this it will not
change your conclusions.

### Why richness is still over-estimated with the gate on

With the gate on, 37 species other than human are still reported against 17 in the
standard. Most of the 22 false positives are not the kind of call the gate is designed to
remove.

Of the 22 false-positive species that survive, **19 are congeners of a genuine member** —
and the other three are defensible. The table also lists *Homo sapiens*, which is not
scored:

| Survivor | Reads after Bracken | Distinct minimizers | What it is |
|---|--:|--:|---|
| *Fusobacterium nucleatum*, *F. vincentii*, *F. polymorphum* | 180,835 | 947–29,678 | Close relatives of the true member, *F. animalis*; all four were subspecies of *F. nucleatum* until *F. animalis*, *F. vincentii* and *F. polymorphum* were raised to species rank. With `k2_pluspf_20251015`, Kraken2 gives *F. nucleatum* more reads than *F. animalis* (16,709 against 6,062) |
| *Homo sapiens* | 965 | 10,062 | Trace human DNA from handling; present in the sample, a contaminant that is not scored |
| *Shigella dysenteriae* | 839 | 1,117 | *Shigella* sits phylogenetically inside *E. coli*; the genomes are near-identical |
| *Sphingomonas paucimobilis*, *Cutibacterium acnes* | 51 | 464–779 | Common extraction-kit and skin contaminants; plausibly present in the sample |
| 16 further congeners of *Veillonella*, *Fusobacterium*, *Saccharomyces*, *Escherichia*, *Bacteroides*, *Candida*, *Clostridioides* | 1,874 | 400–5,380 | Close relatives of an abundant true member |

All these survivors carry **large** distinct-minimizer counts — *Saccharomyces paradoxus*
has 5,380 on just 73 reads, at a duplicity of 1.7. That is broad, un-repetitive evidence,
and it is genuinely in the read data: *S. cerevisiae* is a real member of the standard,
and its sequence matches a great deal of what the database files under its sister
species.

So there are two different ways a species-level call goes wrong, and only one is an
evidence problem:

- **Thin evidence** — a few conserved minimizers hit over and over. Low distinct count,
  high duplicity. **The gate removes these.**
- **Mis-attributed evidence** — a genuinely present organism whose sequence matches a close
  relative's database entry. High distinct count, low duplicity. **The gate does not remove
  these**, and should not: the sequence is in the sample, and only the species label is
  wrong.

The residual inflation is almost entirely the second kind, which is a property of how
finely the reference database splits closely related organisms — not something any
k-mer counting method can detect. The practical consequence is visible in the table above:
at **genus** rank, where those distinctions disappear, the gate brings richness from +47%
to **+6%** of the true value.

### What this does not show

!!! warning "On the mock the gate removed no true species; on CAMI II it changed F1 very little"

    Every organism in D6331 has a species-level entry in the database used. When that is
    true, classification is already close to correct, most false positives are a thin
    tail of taxa with few reads, and there is little for an evidence filter to repair,
    which is what the CAMI II benchmark found: with the Bracken threshold scaled to depth,
    this gate adds at most 0.002 F1 across marine, rhizosphere and strain-madness samples
    (PlusPF `k2_pluspf_20251015`).

    The gate is still useful in cases like JS425 above: a taxon with an ordinary read
    count but few distinct minimizers. A read-count threshold cannot detect this, however
    it is scaled.

The CAMI II benchmark described on the [confidence and threshold](confidence-and-threshold.md)
page covers what one mock cannot: larger and more even communities (marine, 256–381
species per sample), several samples at several depths, communities with genomes absent
from the reference (rhizosphere: only 165 of its 307 gold-standard species, 54 %, are in
PlusPF `k2_pluspf_20260626`), and sample types beyond gut and wastewater. It agrees with
the D6331 result that the gate removes no true species at full sequencing depth, and it
shows that once the Bracken threshold is scaled to depth the gate changes F1 very little.

Still untested: everything above used short-read Illumina data. Long-read or hybrid data,
and sample types with sparser reference coverage than rhizosphere, have not been tried. If
you switch the gate on, check the evidence tables for what it removed.

### Limitation: very close relatives

Very close relatives defeat the method, and no k-mer evidence measure fixes this.
*Shigella dysenteriae* is not in D6331, but 13 reads were assigned to it — carrying 1,117
distinct minimizers, so it passes the gate, and Bracken raised it to 839 reads (834 with
the gate off). *Shigella* and *E. coli* are similar enough that minimizers the database
marks as *Shigella*-specific are genuinely present in *E. coli* reads. That is a property
of the reference taxonomy, not something this filter can detect.

### Threshold stability

A fixed threshold should give a species the same verdict in every sample, but it does not
always. Of 620 species found in more than one of the four PRJNA1477390 test samples,
**123 (20 %)** got different verdicts at the floor of 300:

- **62 % differ mainly in depth**: the species had at least five times fewer reads in one
  sample, and with few reads it cannot reach enough distinct minimizers. This change of
  verdict is probably correct.
- **38 % have read counts less than five-fold apart.** 20 of these 47 species lie close to
  the floor, between 200 and 450 distinct minimizers, in every sample. *Sphingomonas
  paucimobilis* is one: 155 to 199 read pairs in each sample, 350 distinct minimizers in
  SRR39106524 (kept) and 277 to 295 in the other three (removed). The other 27 lie further
  from the floor in at least one sample, sometimes much further: *Sphingopyxis
  lindanitolerans* has 93 to 105 distinct minimizers in three samples and 10,626 in the
  fourth, on 95 to 282 read pairs.

Any fixed cut-off behaves like this, whatever its value. Check the evidence tables for
species whose verdict differs between samples.

## What does not work

Three approaches were tested on the four PRJNA1477390 test samples and rejected. They are
described here because each seems reasonable before it is measured.

**A threshold on distinct minimizers alone is partly an abundance filter.** The distinct
count correlates with read count (Pearson *r* = +0.71 on log-log across the 2,418 species
rows of the four Kraken2 reports), so some of what it removes, it removes for being rare.
This is why the published thresholds pair the k-mer count with an explicit read floor —
`min_reads` exists for exactly that, and is 0 by default only because
[Oskolkov (2026)](https://doi.org/10.3389/fmicb.2026.1603339) found that filtering on
unique k-mers alone was enough to reach the best F1, with no read-count filter needed.

**Normalising to minimizers-per-read inverts the filter.** That ratio *falls* as
abundance rises (*r* = −0.72): across the four samples its median is 84 for single-read
species and 1.4 for species with 1,000 reads or more. A "ratio below 5" rule on
SRR39106521 removed 80 % of the reads that had reached species rank, including the two
species with the most reads, while keeping all 160 single-read species.

**Thresholding duplicity removes the most abundant species.** Across all species rows
duplicity hardly changes with abundance (*r* = +0.17), but the most deeply sequenced
species have high values, because each of their minimizers is seen many times.
*Achromobacter denitrificans*, the largest species in all four samples, has a duplicity of
23 to 74 on 42,000 to 1.1 million distinct minimizers. A duplicity cut-off of 10 prunes
about as many species rows as the floor of 300 (65 % against 66 %), and its verdicts
differ between samples less often (14 % of species against 20 %). But it would remove 73
to 93 % of each sample's table, *A. denitrificans* included, where the floor of 300
removes 1.3 to 5.8 %. With the gate's 1 % safety rail, which would keep
*A. denitrificans*, it would still remove 16 to 65 %. Duplicity is therefore reported, not
thresholded.

**Flagging outliers against a fitted trend line finds no distinct population.** Fitting
`log(distinct minimizers)` against `log(reads)` and taking residuals looks principled,
but the residuals are essentially normal: in SRR39106521 a "2 SD below the line" rule
flagged 14 of 578 species, where a normal distribution predicts 13.1. There is no
separable cluster to find, so the rule just removes the bottom few percent of *any*
sample. Worse, because the line is fitted per sample, the same species can fall on
opposite sides of it in different samples of one study.

The conclusion for the design: use one fixed floor for every sample of a run, so that
every sample is judged against the same number, and treat the trend line as a *ranking
for the eye* rather than a test.
