# Species-level k-mer evidence

!!! note "Read this after [Choosing confidence and threshold](confidence-and-threshold.md)"

    The single biggest precision lever in shotgun mode is Bracken's depth-scaled read
    threshold, covered on that page. Benchmarked across a real mock community and three
    CAMI II environments, this gate adds **at most 0.002 F1** on top of a correctly
    scaled threshold — it is not redundant, but it is a backstop for a narrow class of
    error, not the main event. Read the threshold page first if you have not.

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

    `--report-minimizer-data` is purely additive reporting. Verified directly: the same
    5,000 read pairs against the same database twice, identical settings, differing only
    in this flag, gave a byte-for-byte identical per-read output and identical values in
    every shared report column.

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
02b.evidence/{sample}_species_evidence.tsv   one row per taxon at the Bracken rank
02b.evidence/{sample}_report_gated.txt       the report with failing taxa removed
stats/kmer_evidence/{sample}_evidence.png    reads vs distinct minimizers, log-log
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
kraken2  →  kmer_evidence  →  bracken  →  kraken-biom  →  finalize_otu_table
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
still balances: for every node, `clade == direct + sum(children clade)`. Verified on the
test data — zero violations, zero indentation gaps, and no column other than a parent's
direct count ever changes.

!!! warning "Two things this does not guarantee"

    **Not all reads land on close relatives.** The rationale for crediting the parent is
    that a spurious call usually borrows evidence from a neighbour. That holds for most
    of the traffic, not all: on four test samples **70%** of the reads gained by
    surviving taxa can be matched to a loss inside the same genus; the other **30%**
    crosses genus or family boundaries.

    **Some reads are still lost.** If *every* species under a genus fails, Bracken has
    nothing left at the target rank to redistribute to, and drops those reads. Measured:
    **−0.16%** of reads, with 27–38 genera per sample emptied completely. Crediting the
    parent is much better than deleting rows outright, but it is not conservation.

## Choosing a threshold

There is published guidance, and it is worth reading before inventing a number.

| Source | Recommendation |
|---|---|
| [KrakenUniq](https://link.springer.com/article/10.1186/s13059-018-1568-0) (Breitwieser & Salzberg, 2018) | ~2,000 unique k-mers per million reads; 1,000 for pathogen discovery in patient samples |
| [Refining filtering criteria of the Kraken family](https://pmc.ncbi.nlm.nih.gov/articles/PMC13223122/) (2026) | ~1,000 unique k-mers optimal by F1 across three ancient-metagenome dataset types; read-count filtering "largely unnecessary" when the unique k-mer count is high; suggests scaling with depth, ~200 per 100,000 reads |
| [aMeta](https://link.springer.com/article/10.1186/s13059-023-03083-9) (Pochon et al., 2023) | Ships 1,000 unique k-mers **and** 200 assigned reads as user-configurable defaults |

!!! warning "Those numbers are k-mers. Kraken2 reports minimizers. Not the same unit."

    Minimizers are a *sample* of the k-mers, so the same amount of matched reference
    sequence yields a much smaller number.

    How much smaller depends on the database's `k` and `l`. Each k-mer of length `k`
    spans `w = k − l + 1` candidate minimizer positions of length `l`, and the standard
    expected density of selected minimizers is `2 / (w + 1)`. For the usual nucleotide
    settings (`k = 35`, `l = 31`) that is `w = 5` and a density of `1/3` — roughly **one
    distinct minimizer per three k-mers**.

    So a published *1,000 unique k-mers* is about **333** here. Applying the published
    number directly to the minimizer column makes the filter roughly three times harsher
    than its authors intended.

    Check `k` and `l` on the first line of your database's `inspect.txt` and recompute if
    they differ. (The `2/(w+1)` density is the expectation for a random minimizer scheme;
    Kraken2's spaced-seed masking shifts it slightly, so treat 333 as an order-of-
    magnitude conversion, not a precise one.)

### What that does to real data

Four wastewater-enrichment metagenomes (naproxen-degrading consortia from a treatment plant), 2,312 species calls pooled. **Read the denominators** — most reads
in a metagenome stop at genus or above and are never candidates for a species-level
filter, so the two columns differ by a factor of ~30:

| Threshold | Species kept | % of species-rank reads kept | % of all classified reads kept |
|---|--:|--:|--:|
| none | 2312 (100%) | 100.0% | 100.00% |
| ≥ 200 | 1050 (45.4%) | 99.6% | 99.85% |
| **≥ 333** (converted, **default**) | **754 (32.6%)** | **98.7%** | **99.54%** |
| ≥ 500 | 541 (23.4%) | 97.8% | 99.20% |
| ≥ 1000 (published number, unconverted) | 345 (14.9%) | 96.4% | 98.70% |

Only 36.5% of classified reads reach species rank at all. The converted threshold has the
shape a useful filter should have: it drops two thirds of the species *names* while
keeping 98.7% of the reads that carried a species label.

### What it does to the final table

The same four samples, gate off versus on:

| | Gate off | Gate on |
|---|--:|--:|
| Rows in `otu_table.tsv` | 357 | 233 |
| Total counts | 9,287,679 | 9,271,415 (−0.18%) |
| Genera represented | 121 | 71 |

Bray–Curtis dissimilarity between the two profiles is **0.026–0.098** at species level and
**0.005–0.031** at genus level. The genus-level picture barely moves, which is what a
filter that corrects *names* rather than *biology* should look like.

!!! danger "It is not only a tail-trim — check before you believe it"

    The per-sample log reports the pruned reads using Kraken2's own counts. Bracken
    redistributes higher-rank reads down onto species, so a pruned taxon may have been
    holding far more of the **final** table than it holds in the Kraken2 report — on the
    test data, **6–19× more**. Two or three of every sample's top 20 species disappear.

    The clearest case, and the reason to keep this gate on even though the [depth-scaled
    Bracken threshold](confidence-and-threshold.md) does most of the precision work
    elsewhere: *Pseudomonas* sp. JS425 had 373 Kraken2 reads on 189 distinct minimizers
    (each position hit ~62 times) — comfortably above any sane *read-count* floor — and
    Bracken inflated it to **64,961 reads, 3% of the sample and rank 9 in the community**.
    A 174× amplification of a thinly-evidenced call that a threshold on read count alone
    cannot see, because the read count itself looked fine. This is a narrow failure mode,
    not the typical one, but when it happens it is large.

    Run once with the gate off, once with it on, and compare `03.abundance` before
    committing to a threshold.

!!! danger "The depth-scaling rule does not extrapolate"

    The published thresholds come from ancient-DNA and clinical-pathogen work, mostly at
    shallower depth than a modern shotgun run. The depth-scaling suggestion in particular
    should be treated with caution: extrapolated to a multi-million-pair sample it implies
    a threshold roughly ten times the fixed one, far outside the range it was fitted on.
    MetaFlux uses the fixed converted value and leaves the scaling to you.

## Validation status

The default is `enabled: false`. Not because the gate is untested — it has since been run
across 45 real and simulated samples spanning four environments — but because, as the
benchmark below shows, most of what it would change on a well-covered sample is now
already handled by the [depth-scaled Bracken threshold](confidence-and-threshold.md), and
turning on a second filter you have not looked at is still not something to do by default.

### What has been checked

The original single-mock validation below (Zymo D6331) still stands and is kept for the
worked detail. It has since been extended to [CAMI
II](https://www.nature.com/articles/s41592-022-01431-4)'s marine, rhizosphere, and
strain-madness gold-standard sets — 41 samples, two Kraken2 database builds eight months
apart, scored against OPAL — which is where the size and depth questions raised at the
end of this section were actually answered. The full methodology and the numbers that
matter for choosing a threshold live on the [confidence and
threshold](confidence-and-threshold.md) page; the finding specific to *this* gate is that
it contributes essentially nothing (≤0.002 F1) once the Bracken threshold is scaled
correctly, but still catches the JS425-style case above that a read-count floor cannot.

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
    exactly what the CAMI II follow-up found: with the Bracken threshold scaled to depth,
    this gate adds at most 0.002 F1 across marine, rhizosphere and strain-madness samples.

    The JS425 case above is the one place it still earns its keep: a taxon whose *read
    count* looks unremarkable but whose *distinct sequence support* does not. A read-count
    floor, however well scaled, cannot see that distinction by construction.

The four gaps originally listed here — larger and more even communities, several samples
at several depths, communities containing genomes absent from the reference, and sample
types beyond gut and wastewater — have since been closed by the CAMI II benchmark
described on the [confidence and threshold](confidence-and-threshold.md) page: 41 samples,
marine (256–478 species/sample) and rhizosphere (52% of gold species absent from the
database — the exact "genomes absent from the reference" case), across a depth range and
two database builds. The result did not overturn the D6331 finding, it confirmed the
gate's role is narrower than that finding alone suggested: safe, rarely decisive.

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
exactly that, and is 0 by default only because the 2026 paper found read filtering
largely unnecessary once the k-mer count is high.

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

The lesson driving the design: use an **absolute** threshold, so a species gets the same
verdict in every run, and treat the trend line as a *ranking for the eye* rather than a
test.
