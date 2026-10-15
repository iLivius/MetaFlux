# Species-level k-mer evidence

!!! note "Read this after [Choosing confidence and threshold](confidence-and-threshold.md)"

    The single biggest precision lever in shotgun mode is Bracken's depth-scaled read
    threshold, covered on that page. Set that first. This gate is a backstop for a narrow
    class of error the threshold cannot see — valuable when it fires, but off by default;
    see [Why the default is `enabled: false`](#why-the-default-is-enabled-false).

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

    The obvious reading — *the minimizers carried by the reads that were assigned to
    this taxon* — is wrong, and it changes how you interpret the numbers.

    Kraken2 attributes each minimizer **occurrence** to whichever taxon that minimizer
    maps to **in the database** (its lowest common ancestor), pooled over every read in
    the sample, regardless of where each read finally got assigned.

    The give-away in any real report is taxa with a couple of assigned reads carrying
    hundreds of thousands of minimizer counts. In one test sample a *Streptomyces* row
    has **2 clade reads and 345,108 total minimizers** — 737× more than the ~468 k-mer
    positions two read pairs physically contain. Under the per-read reading that is
    impossible; under LCA attribution it is ordinary.

    So column 5 for a species means: **how many distinct database minimizers specific to
    that species were seen anywhere in this sample.** That is a statement about how much
    species-specific sequence is present — which is exactly the quantity KrakenUniq
    calls *unique k-mers*, so the published thresholds below do transfer once the units
    are converted. It is **not** a genome-coverage measure, and describing it as "how
    much of the genome the reads covered" would be wrong.

!!! note "The flag itself changes no classification"

    `--report-minimizer-data` is purely additive reporting: with or without it, Kraken2
    writes a byte-for-byte identical per-read output and the same values in every shared
    report column.

## Reading the numbers

Two species from one real sample, with similar read counts:

| Species | Reads | Distinct minimizers | Total minimizers | Duplicity |
|---|--:|--:|--:|--:|
| *Achromobacter xylosoxidans* | 7,756 | 224,500 | 980,110 | 4.4 |
| *Novosphingobium* sp. RL4 | 4,597 | 426 | 196,778 | 461.9 |

Similar abundance, utterly different evidence. The first is supported by 224,500
distinct species-specific minimizers. The second is supported by 426, each one seen
about 460 times.

**Duplicity** (total ÷ distinct) is the most directly readable of the three numbers,
because it does not require knowing the sequencing depth. It also barely tracks
abundance — measured Pearson *r* = **+0.20** against log reads, versus **+0.71** for the
raw distinct count. That makes it the better number to *eyeball*. It is nonetheless a
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
    min_distinct_minimizers: 333       # the floor; see the unit conversion below
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

Two details have to be right or the pruned report quietly corrupts Bracken's arithmetic.

**Whole subtree, never a bare row.** Bracken rebuilds the taxonomy tree from the
*indentation* of the report, attaching each row to the last row one level above it.
Delete a species row and leave its strain (`S1`) children behind, and those children get
re-parented onto whichever species happens to precede them — silently inflating an
innocent neighbour. So when a taxon fails, everything indented beneath it goes too.

**The reads go back to the parent, not in the bin.** A pruned taxon's whole clade total
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

!!! warning "Those numbers are k-mers. Kraken2 reports minimizers. Not the same unit."

    **In plain terms.** Kraken2 reads a sequence as overlapping 35-base k-mers, but it
    does not look each one up. For every k-mer it picks one 31-base piece inside it, the
    *minimizer*, and stores and counts only that. Neighbouring k-mers usually pick the
    same piece, so a stretch of genome gives about three times fewer distinct minimizers
    than distinct k-mers. The studies above counted unique k-mers and settled on about
    1,000. In a Kraken2 report the same amount of sequence shows up as about **333**
    distinct minimizers, and that conversion is where the default comes from. It was not
    tuned on MetaFlux's own test data, and it is a cautious value rather than an
    equivalent: when Oskolkov (2026) ran both tools on the same simulated data, the best
    Kraken2 floor was two to three times *higher* than the best KrakenUniq one, and
    MetaFlux's own benchmark finds useful floors in the thousands (see
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
    Kraken2's spaced-seed masking shifts it slightly, so treat 333 as an order-of-
    magnitude conversion, not a precise one.)

### What it does on real data

The numbers below come from the four metagenomes of BioProject
[PRJNA1477390](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1477390) (runs
SRR39106521–SRR39106524), microbial consortia enriched on naproxen from wastewater
treatment plants. They were processed with the default settings and the PlusPF index
`k2_pluspf_20260626`, once with the gate off and once with it on. The floor
is `min_distinct_minimizers`: the number of distinct minimizers a species needs in the
Kraken2 report to stay in the table.

| | SRR39106521 | SRR39106522 | SRR39106523 | SRR39106524 |
|---|--:|--:|--:|--:|
| Species in `otu_table.tsv`, gate off | 99 | 80 | 116 | 120 |
| Species, gate on (floor 333) | 84 | 65 | 93 | 101 |
| Genera, gate off → on | 39 → 34 | 31 → 28 | 46 → 39 | 52 → 43 |
| Read pairs moved to another species | 4.9 % | 1.6 % | 5.8 % | 2.1 % |
| Change in total read pairs | −0.01 % | +0.01 % | −0.01 % | −0.09 % |
| Genomes detected by sylph in the same reads | 10 | 9 | 14 | 16 |

- **The gate removes 15 to 23 species per sample and almost no reads.** The reads of a
  removed species go back to its parent taxon, and Bracken shares them among the
  relatives that passed.
- **The removed species are not only rare ones.** Bracken spreads the reads that Kraken2
  could only place at genus level or above over the species below, and it had given the
  removed species 11 to 27 times more reads than Kraken2 assigned to them. In every
  sample one of the 20 most abundant species is removed. The largest case is
  *Pseudomonas* sp. JS425 in SRR39106521: 372 read pairs from Kraken2, carrying 182
  distinct minimizers, raised by Bracken to 69,992 read pairs, the eighth most abundant
  species in the sample. Its read count passes the Bracken threshold; only the
  distinct-minimizer count shows how little of its genome is behind it (more
  [below](#but-no-change-in-f1-is-not-no-change-in-the-table)).
- **Fewer species is the right direction for these samples.** sylph, run by MetaFlux on
  the same reads in `both` mode, detects 9 to 16 genomes per sample.
- **The gate does not bring the list down to that size.** Most of what remains is a
  database problem. The most abundant organism in SRR39106521 and SRR39106523,
  *Pusillimonas noertemannii*, has no genome in PlusPF, and its reads are reported under
  several *Achromobacter* species, each with tens of thousands of distinct minimizers
  (see the warning on missing organisms in [Databases](databases.md#choosing-an-index)).
  A k-mer floor cannot correct that.

The same four samples at other floors: species left in `otu_table.tsv`, and the share of
the gate-off table held by the species that the floor removes.

| Floor | SRR39106521 | SRR39106522 | SRR39106523 | SRR39106524 | Read pairs held by removed species |
|---|--:|--:|--:|--:|--:|
| none (gate off) | 99 | 80 | 116 | 120 | 0 |
| 200 | 92 | 73 | 105 | 112 | 0.9–4.7 % |
| **333 (default)** | **84** | **65** | **93** | **101** | **1.6–5.8 %** |
| 1,000 | 55 | 44 | 61 | 80 | 8–21 % |
| 3,000 | 29 | 26 | 31 | 54 | 13–42 % |

A lower floor keeps the species with 200 to 333 distinct minimizers, 8 to 12 more per
sample. These samples have no known composition, so they cannot show whether those
species are real. A higher floor removes many more species but also moves a large share
of the table: at 3,000, 13 to 42 % of the read pairs end up under a different species. On
benchmark communities of known composition, a fixed floor of 10,000 already removed true
species from some samples, and lower fixed floors did the same once shallower sequencing
was modelled (3,000 at 30 % of the depth, 1,000 at 10 %). A floor scaled to the depth of
the sample avoids this; see the next section.

### Scaling the floor with depth

A species that is really present leaves more distinct minimizers the deeper a sample is
sequenced, and so, up to a limit, does a spurious call that borrows sequence from an
abundant relative. A fixed floor is therefore lenient on deep samples and strict on
shallow ones. KrakenUniq reported the same for unique k-mers (table above), but neither
publication gives a rule for Kraken2.

MetaFlux's benchmark was used to fit one. It covers 42 samples of known composition: one
real mock community (ZymoBIOMICS D6331, SRA run SRR38976715) and 41 simulated CAMI II
samples from the marine, rhizosphere and strain-madness sets, with 2.6–6.5 million
classified read pairs each, classified with three Kraken2 indices: `k2_pluspf_20260626`
(at confidence 0.15 and 0.20), `k2_gtdb_genome_reps_20250609` (GTDB r226 representative
genomes) and `k2_NCBI_reference_20251007` (one genome per NCBI species). A floor of **one
distinct minimizer per 1,000 classified read pairs** gave the best result that was safe
everywhere:

- species-level F1 rose by 0.064 on average (0.029 with PlusPF), against 0.0025 for the
  fixed 333 (0.0004 with PlusPF);
- it removed 2,286 of 13,114 false species calls and 6 of 13,419 true ones (none with
  PlusPF), and lowered no sample's F1;
- it stayed safe when shallower sequencing was modelled, down to 1 % of the reads, while
  every fixed floor from 333 to 10,000 lowered some samples' F1 at some depth.

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
500,000 classified pairs the gate changes little, whatever the floor.

Limits of this rule:

- **It rests on limited data.** Three quarters of the F1 gain came from the one real mock
  community; on the simulated CAMI samples it was 0.01–0.03. The shallower depths were
  modelled from the observed counts, not sequenced.
- **Do not go above 1.5 per 1,000 pairs.** At that rate rare true species start to be
  removed (92 of 13,419 true calls), and from 2 per 1,000 the F1 falls.
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
- **It changes the table more than 333 does.** On the test dataset above (1.7–3.5 million
  classified pairs), a floor of 1,674 sized on the shallowest sample keeps 33–68 species
  per sample, against 65–101 at 333, and the species it removes hold 10–33 % of the read
  pairs, which move to relatives. That is still well above the 9–16 genomes sylph finds.
  Compare `otu_table.tsv` with the gate off and on before relying on it.

## Why the default is `enabled: false`

Not because the gate is untested — it has been run across 45 real and simulated samples
spanning four environments — but because a correctly scaled
[Bracken threshold](confidence-and-threshold.md) already removes most of what it would
remove. Measured at the depth-scaled threshold, the gate changes species-level F1 by a
mean of **+0.0009**, and by *exactly* zero on the strain-madness and rhizosphere sets at
every depth tested. Left at the old fixed `-t 10` it was worth about nine times more
(mean +0.0079) — which is the clearest statement of its role: **it is a safety net for an
untuned pipeline, not a second independent filter.**

!!! note "Which database build those two numbers come from"

    Both figures come from gate-on/gate-off comparisons on one PlusPF build. The
    mechanism does not depend on the database, so they are not expected to shift much
    with a newer one.

### But "no change in F1" is not "no change in the table"

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

**When to turn it on:** environmental or otherwise poorly-covered communities, where reads
from organisms absent from the database have to land somewhere. For gut, clinical, or
anything close to RefSeq, leave it off. Either way the preview gated report is written
even when disabled, so `diff` it against the Kraken2 report before deciding.

## Validation status

What the gate was tested against, and what those tests could and could not show.

### What has been checked

The gate was checked on one real mock community (Zymo D6331, below, kept for the worked
detail) and on [CAMI II](https://www.nature.com/articles/s41592-022-01431-4)'s marine,
rhizosphere and strain-madness gold-standard sets — 41 samples, two PlusPF builds,
scored with OPAL. The full methodology and the numbers that matter for choosing a
threshold live on the [confidence and threshold](confidence-and-threshold.md) page; the
finding specific to *this* gate is that it contributes essentially nothing (≤0.002 F1)
once the Bracken threshold is scaled correctly, but still catches the JS425-style case
above that a read-count floor cannot.

One real mock community: the [ZymoBIOMICS Gut Microbiome
Standard (D6331)](https://files.zymoresearch.com/protocols/_d6331_zymobiomics_gut_microbiome_standard.pdf),
a shotgun run of ~5.3 M classified read pairs, processed twice through the whole pipeline
with the gate off and on.

The standard is sold as 21 strains. A species-level profiler can only ever see **17
species**, because five of the 21 are *Escherichia coli* strains at 2.8% each that collapse
into one *E. coli* row at 14%. That ceiling is arithmetic, not a reference-database gap.

| | Gate off | Gate on |
|---|--:|--:|
| True members recovered | 15 / 17 | **15 / 17** |
| False-positive species | 35 | **23** |
| Precision | 0.300 | **0.395** |
| F1 | 0.448 | **0.545** |

- **No true member was removed.** *Salmonella enterica* — 0.01% of the standard, just 134
  classified reads — carries 3,547 distinct minimizers at a duplicity of 1.5, ten times
  clear of the threshold, and survives. That is the distinction the whole method rests on:
  *rare* is not the same as *poorly evidenced*.
- The two members never recovered (*Enterococcus faecalis* at 0.001%, *Clostridium
  perfringens* at 0.0001%) are absent from the Kraken2 report in **both** runs. Both have
  species entries in the database, so this is a depth limit, not a reference gap and not
  the gate.
- **Quantitative accuracy did not change.** L1 distance to the theoretical composition was
  43.44 either way, Shannon *H* 2.146 either way, and Bray–Curtis between the two profiles
  0.00009 — numerically indistinguishable. The 82 taxa the gate removed held 426 reads
  between them out of 5.3 million. Any ordination or differential-abundance test would
  give the same answer from both tables.
- **Richness, on the other hand, improved substantially** — which is where the gate did its
  work on this dataset:

| Observed richness | Truth | Gate off | Gate on |
|---|--:|--:|--:|
| Species | 17 | 50 (+194%) | **38 (+124%)** |
| Genera | 17 | 26 (+53%) | **19 (+12%)** |

The split is clean: **the gate cleans the name list, it does not touch the abundances.**
If you report how many taxa you detected, or work with presence/absence, that matters. If
you only ever compare relative abundances between samples, on data like this it will not
change your conclusions.

### Why richness is still over-estimated with the gate on

38 species against a true 17 is still nearly double, and it is worth being precise about
what the extra 21 are, because they are not what the gate is built to catch.

Of the 23 non-truth species that survive, **19 are congeners of a genuine member** — and
all four of the remainder are defensible:

| Survivor | Reads | Distinct minimizers | What it is |
|---|--:|--:|---|
| *Fusobacterium animalis*, *vincentii*, *polymorphum* | 82,228 | 947–13,157 | Former **subspecies of *F. nucleatum***, elevated to species rank. Not errors — the mock's genuine member, split by a taxonomy change |
| *Homo sapiens* | 965 | 10,062 | Trace human DNA from handling. Really there |
| *Shigella dysenteriae* | 839 | 1,117 | *Shigella* sits phylogenetically inside *E. coli*; the genomes are near-identical |
| *Sphingomonas paucimobilis*, *Cutibacterium acnes* | 51 | 464–779 | Textbook extraction-kit and skin contaminants. Plausibly really there |
| 14 further congeners of *Veillonella*, *Bacteroides*, *Escherichia*, *Saccharomyces*, *Candida*, *Clostridioides* | ~1,900 | 400–5,380 | Close relatives of an abundant true member |

The pattern is the point. These survivors all carry **large** distinct-minimizer counts —
*Saccharomyces paradoxus* has 5,380 on just 73 reads, at a duplicity of 1.7. That is broad,
un-repetitive evidence, and it is genuinely in the read data: *S. cerevisiae* is a real
member of the standard, and its sequence matches a great deal of what the database files
under its sister species.

So there are two different ways a species-level call goes wrong, and only one is an
evidence problem:

- **Thin evidence** — a few conserved minimizers hit over and over. Low distinct count,
  high duplicity. **The gate removes these.**
- **Mis-attributed evidence** — a genuinely present organism whose sequence matches a close
  relative's database entry. High distinct count, low duplicity. **The gate cannot remove
  these, and arguably should not: the sequence really is present. The label is what is
  wrong, not the evidence.**

The residual inflation is almost entirely the second kind, which is a property of how
finely the reference database splits closely related organisms — not something any
k-mer counting method can detect. The practical consequence is visible in the table above:
at **genus** rank, where those distinctions disappear, the gate brings richness from +53%
to **+12%** of the true value.

### What this does not show

!!! warning "A mock community proved the gate safe. Bigger, gold-standard communities showed it is not the main event."

    Every organism in D6331 has a species-level entry in the database used. When that is
    true, classification is already close to correct, false positives are a thin tail of
    near-zero-read taxa, and there is little for an evidence filter to repair — which is
    exactly what the CAMI II benchmark found: with the Bracken threshold scaled to depth,
    this gate adds at most 0.002 F1 across marine, rhizosphere and strain-madness samples.

    The JS425 case above is the one place it still earns its keep: a taxon whose *read
    count* looks unremarkable but whose *distinct sequence support* does not. A read-count
    floor, however well scaled, cannot see that distinction by construction.

The CAMI II benchmark described on the [confidence and threshold](confidence-and-threshold.md)
page covers what one mock cannot: larger and more even communities (marine, 256–478
species per sample), several samples at several depths, communities with genomes absent
from the reference (rhizosphere, where 52% of gold species are missing from the
database), and sample types beyond gut and wastewater. It confirms the D6331 result and
narrows the gate's role: safe, rarely decisive.

What genuinely remains open: everything above used short-read Illumina data. Long-read
or hybrid assemblies, and sample types with even sparser reference coverage than
rhizosphere, are untested. Treat the gate as a tool to run with the evidence tables open
beside it, not as a setting to switch on and forget.

### A limitation worth knowing about

Very close relatives defeat the method, and no k-mer evidence measure fixes this.
*Shigella dysenteriae* is not in D6331, but 13 reads were assigned to it — carrying 1,117
distinct minimizers, so it passes the gate, and Bracken raised it to 834 reads. *Shigella*
and *E. coli* are similar enough that minimizers the database marks as *Shigella*-specific
are genuinely present in *E. coli* reads. The evidence is real; the label is wrong. That is
a property of the reference taxonomy, not something this filter can detect.

### Threshold stability

An absolute threshold gives a species the same verdict in every sample — up to a point.
Of 594 species seen in more than one test sample, **119 (20%)** got inconsistent verdicts.
Splitting them:

- **55% are depth-driven** — the organism was genuinely 5× rarer in one sample, and too
  few reads simply cannot accumulate enough distinct minimizers. Arguably correct.
- **45% are knife-edge** — similar read counts, verdict still flips, because the species
  sits within a whisker of the line. *Novosphingobium* sp. RL4 above is one: 426 distinct
  minimizers in one sample (kept) and 293–301 in the other three (pruned), despite
  duplicity of 460–800 marking it as spurious in all four.

This is inherent to any hard cut-off, not a defect of this particular number. It is the
main reason to read the evidence tables rather than trust a threshold blindly.

## What does not work

Three approaches were tested against real data and rejected. Recorded because each looks
reasonable until it is measured.

**A threshold on distinct minimizers alone is partly an abundance filter.** The distinct
count correlates with read count (Pearson *r* = +0.71 on log-log across 2,312 species
calls), so some of what it removes, it removes for being rare. This is why the published
thresholds pair the k-mer count with an explicit read floor — `min_reads` exists for
exactly that, and is 0 by default only because
[Oskolkov (2026)](https://doi.org/10.3389/fmicb.2026.1603339) found that filtering on
unique k-mers alone was enough to reach the best F1, with no read-count filter needed.

**Normalising to minimizers-per-read inverts the filter.** That ratio *falls* as
abundance rises (*r* = −0.74) — median 92 for single-read species, 0.68 for species with
thousands of reads. A "ratio below 5" rule on one test sample removed 80% of the reads
that had reached species rank, including the two most abundant species in the sample,
while retaining every one of the 158 single-read species. It removes the data and keeps
the noise.

**Thresholding duplicity is more stable-looking but less stable in practice.** Duplicity
is nearly independent of abundance, which makes it attractive. But as a *filter* it is
worse: at a comparable pruning rate it gave 28.5% inconsistent verdicts across samples
against 20% for the distinct count. It is only consistent at cut-offs so permissive they
remove almost nothing. Reported, not thresholded.

**Flagging outliers against a fitted trend line finds no distinct population.** Fitting
`log(distinct minimizers)` against `log(reads)` and taking residuals looks principled,
but the residuals are essentially normal: a "2 SD below the line" rule flagged 13 species
where a pure Gaussian predicts 12.4. There is no separable cluster to find, so the rule
just removes the bottom few percent of *any* sample. Worse, because the line is fitted per
sample, the same species can fall on opposite sides of it in different samples of one
study.

The conclusion for the design: use one fixed floor for every sample of a run, so a
species gets the same verdict in every sample, and treat the trend line as a *ranking for
the eye* rather than a test.
