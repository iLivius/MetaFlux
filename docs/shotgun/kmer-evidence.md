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
    min_distinct_minimizers: 300       # the floor: a whole number, or auto (see below)
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
    nothing left at the target rank to redistribute to, and drops those reads; the same
    happens one level up when every species of a family or an order fails. At the
    default floor of 300 that is a small fraction of the sample: on the
    [test dataset](#what-it-does-on-real-data) the table lost at most 0.09 % of its read
    pairs. With `auto` it was 2.7 to 15.1 % of the same tables, because a whole order
    failed in two of them (see [Limits of the rule](#limits-of-the-rule)). Crediting the
    parent is much better than deleting rows outright, but it is not conservation.

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
[PRJNA1477390](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1477390) (runs
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
| Read pairs held by the removed species | 4.9 % | 1.3 % | 5.8 % | 1.8 % |
| Change in total read pairs | −0.01 % | +0.01 % | −0.01 % | −0.09 % |
| Genomes detected by sylph in the same reads | 10 | 9 | 14 | 16 |
| Species above 1 % of `otu_table.tsv`, gate off | 21 | 10 | 23 | 12 |
| Abundant members counted from the 16S rRNA gene, no reference | 7 | 6 | 4 | 10 |

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
- **Fewer species is the right direction for these samples.** sylph, run by MetaFlux in
  `both` mode against `gtdb-r232-c200-dbv1` and `fungi-refseq-2025-10-11-c200-dbv1`,
  detects 9 to 16 genomes per sample; it misses organisms without a closely related
  reference genome. A count that needs no reference database comes from the 16S rRNA
  gene fragments in the reads, extracted with Metaxa2 2.2.3: fragments from the same part
  of the gene are compared with one another and grouped at 97 % identity. It sees only
  the more abundant members, roughly those above 1–2 % of the 16S gene copies, and finds
  4 to 10 per sample, mostly from the same families that sylph and Kraken2 report. For that
  abundant part of the community Kraken2 reports 10 to 23 species, those above 1 % of
  its table, because it splits the reads of one organism over several related
  species.
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
species are real. A higher floor removes many more species, and those species hold a
large share of the table: at 3,000, 13 to 42 % of the read pairs. Part of them moves to
other species and part leaves the table (measured with `auto`, whose floors on these
samples run from 1,673 to 3,458: see [Limits of the rule](#limits-of-the-rule)). On
benchmark communities of known composition, a fixed floor of 10,000 already lowered some
samples' F1 by removing true species, and lower fixed floors did the same once shallower
sequencing was modelled (3,000 at 30 % of the depth, 1,000 at 10 %). A floor scaled to
the depth of the sample avoids removing true species that way; see the next section.

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
benchmark.

#### `min_distinct_minimizers: auto`

`auto` applies this rule for you, sample by sample:

```yaml
shotgun:
  kmer_evidence:
    enabled: true
    min_distinct_minimizers: auto
```

After Kraken2 has run, rule `kmer_evidence` works out each sample's floor from that
sample's own Kraken2 report:

```text
classified pairs   clade reads of the root row
human pairs        clade reads of the Homo sapiens row (0 if the report has none)
floor              (classified pairs − human pairs) / 1,000, at most 6,500
```

| Classified read pairs, human excluded | Floor set by `auto` |
|--:|--:|
| 150,000 | 150 |
| 500,000 | 500 |
| 1 million | 1,000 |
| 2 million | 2,000 |
| 5 million | 5,000 |
| 6.5 million or more | 6,500 |

- **One floor per sample.** With a number, MetaFlux applies that one floor to every
  sample of the run. With `auto`, each sample is judged against its own floor, so a deep
  sample gets a higher floor than a shallow one in the same run. The four samples of the
  [test dataset](#what-it-does-on-real-data) get floors of 1,673 to 3,458.
- **A species stays when its count reaches the floor**: distinct minimizers ≥ the
  unrounded floor (1,672.619 in the example below), as in the benchmark, and it must also
  clear `min_reads`. The 1 % safety rail works as with a number; its share is of all
  classified reads, human ones included.
- **Human reads are left out.** PlusPF contains the human genome and host removal is off
  by default, so human reads would raise the floor without adding any microbial
  evidence. This is MetaFlux's addition to the rule: in the benchmark, the host reads
  left in a sample counted as classified. The row is recognised by its taxid (9606) and
  its name, so a database whose taxids are not NCBI's (a GTDB-based one) is not affected.
  Reads of other hosts are not subtracted. With PlusPF they mostly stay unclassified; a
  database that holds other host genomes classifies them, and they then raise the floor:
  remove the host first ([Decontamination](decontamination.md#host-removal)) or set a
  number. The *Homo sapiens* row itself is still judged against the floor: in all four
  samples of the test dataset it fails (605 to 1,743 distinct minimizers) and leaves the
  table with its read pairs, where at 300 it stays.
- **The cap of 6,500.** The benchmark stops at 6.5 million classified pairs. A floor that
  keeps growing with depth eventually exceeds what a real species with a small
  species-specific part of its genome can reach; in the model this began at about 10
  million pairs.
- **Shallow samples.** Below about 500,000 classified pairs the scaled floor gains little,
  and below 300,000 `auto` sets it lower than the default of 300. That is intended: in
  benchmark samples thinned to 26,000–195,000 classified pairs, the fixed floor of 300
  lowered some samples' F1 and the scaled floor lowered none. On the one benchmark
  sample really run at 0.76 million classified pairs, `auto` (floor 759) and 300 removed
  the same three false species
  ([below](#measured-at-known-depths-and-on-deeper-samples)).
- **Where to see the floor.** Each sample's log, `logs/kmer_evidence/{sample}.log`, shows
  how its floor was worked out:

    ```text
    [kmer_evidence] auto floor: 1673810 classified read pairs, less 1191 Homo sapiens read pairs (taxid 9606): 1672619 / 1000 = 1672.619. A taxon at rank S needs at least 1673 distinct minimizers.
    ```

    The evidence table gives the floor in its `reason` column (for *Sphingobium
    wenxiniae* in that sample, `distinct minimizers 1664 < 1673 (auto)`), and the plot
    draws it, labelled `(auto)`.
    Both show it rounded up: the counts are whole numbers, so "at least 1,672.619" and
    "at least 1,673" are the same test.

To size a number by hand instead, for instance to judge every sample of a run against the
same floor, divide the classified read pairs of the shallowest sample, less its human
read pairs, by 1,000. Each sample's classified count is in `stats/read_tracking.txt`
(column `classified`), and its human read pairs are in the *Homo sapiens* row of its
Kraken2 report. The deeper samples are then judged against a lower floor than `auto`
would set them: more lenient where their extra depth is spread over the community, but
not where it sits on one dominant organism, whose reads raise `auto`'s floor for the
other species without adding to their evidence (see the end of this page). On the test
dataset the single floor of 1,673 keeps 41, 33, 45 and 68 species where `auto` keeps 35,
25, 45 and 61.

#### Why `auto` is not the default

`auto` is opt-in, and the default stays 300:

- **Its gain rests mostly on one real mock community.** Three quarters of the F1 gain came
  from the ZymoBIOMICS D6331 sample; on the three simulated CAMI II sets (41 samples) the
  mean gain per set was 0.01–0.03.
- **Part of its safety was modelled, not measured.** For the benchmark samples the
  shallower depths were modelled from the observed counts, not sequenced. One of them
  was then subsampled to three shallower depths and run for real with PlusPF, and
  behaved largely as the model predicted
  ([below](#measured-at-known-depths-and-on-deeper-samples)); the other 41 rest on the
  model. No benchmark sample is deeper than 6.5 million classified pairs, which is why
  there is a cap.
- **It changes the table much more than 300 does, and read pairs leave the table**: 2.7
  to 15.1 % of the read pairs on the test dataset, against at most 0.09 % at 300, most of
  them from an organism with no close reference genome (the last point under
  [Limits of the rule](#limits-of-the-rule)).

#### Measured at known depths and on deeper samples

What `auto` does on samples sequenced at known depths, and on samples deeper than any
in the benchmark:

**Measured at four depths.** CAMI II marine sample 0, one of the benchmark samples, was
subsampled to 2, 5 and 10 million read pairs and run together with its full 16.6
million, with the settings of its benchmark run (PlusPF `k2_pluspf_20260626`, Kraken2
confidence 0.15): once with the gate off, once at 300 and once with `auto`. After
quality filtering, Kraken2 classified 0.76, 1.90, 3.79 and 6.31 million pairs. Species
were scored as in the benchmark, against the sample's 252 gold-standard species:

| Classified read pairs | 0.76 M | 1.90 M | 3.79 M | 6.31 M |
|---|--:|--:|--:|--:|
| Floor set by `auto` | 759 | 1,896 | 3,788 | 6,309 |
| True / false species, gate off | 217 / 71 | 216 / 69 | 216 / 68 | 216 / 68 |
| True / false species, 300 | 217 / 68 | 216 / 67 | 216 / 66 | 216 / 67 |
| True / false species, `auto` | 217 / 68 | 216 / 64 | 216 / 63 | 216 / 61 |
| Recall, every arm | 0.861 | 0.857 | 0.857 | 0.857 |
| Precision, gate off → 300 → `auto` | 0.754 → 0.761 → 0.761 | 0.758 → 0.763 → 0.771 | 0.761 → 0.766 → 0.774 | 0.761 → 0.763 → 0.780 |
| Species-level F1, gate off | 0.804 | 0.804 | 0.806 | 0.806 |
| F1 change, 300 | +0.0045 | +0.0030 | +0.0030 | +0.0015 |
| F1 change, `auto` | +0.0045 | +0.0076 | +0.0076 | +0.0107 |
| True / false species, `auto`, as the depth model predicted | 216 / 66 | 216 / 64 | 216 / 63 | the model's starting point |

- **No true species was removed**, at any depth, by either floor: recall is the same in
  the three arms of every rung. Relative to its floor, the true species closest to
  removal was *Aliivibrio wodanis* at full depth, 7,263 distinct minimizers against
  6,309 (1.15 times the floor); at 0.76 million the weakest true species,
  *Synechococcus* sp. WH 8020, had 1,186 against 759.
- **The depth model held where it was really tested.** The model starts from the
  full-depth run and thins its calls to a shallower depth, so the 6.31 M column is its
  starting point, not a test, and the two middle rungs are the real tests: at 1.90 and
  3.79 million it gave the true and false species of `auto` exactly, and the true
  species kept the distinct minimizers it predicted for them (median within 1 %). At
  0.76 million it was off in two ways. Species near Bracken's read threshold come and go
  by chance when reads are subsampled (six crossed it at 0.76 million, one at 1.90),
  which the model cannot see and no floor can correct; and three false species had 1.6
  to 3.6 times fewer distinct minimizers than the model gave them, so at that depth 300
  removed as many false species as `auto` did.
- **The read pairs in the table hardly changed**: by less than 0.05 %, with either floor,
  at every depth.

This is one simulated sample on one index: it measures the depths the benchmark only
modelled, for this sample, and says nothing about the other 41 benchmark samples.

**A deep mock community.** MBARC-26 is a mock of 23 bacteria and 3 archaea with known
genomes ([Singer et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5037974/)), here
in the library sequenced by
[Hardwick et al. (2018)](https://doi.org/10.1038/s41467-018-05555-0) (SRA run
SRR6392275). With PlusPF `k2_pluspf_20260626` it has 51.9 million classified read pairs,
eight times the deepest benchmark sample, so `auto` sets the cap, 6,500; without the cap
the floor would be 51,851. The gate was run on the existing Kraken2 report and Bracken
was not run again, so what follows is which species each floor keeps, not how the read
counts move.

| Floor | Members removed, of 26 | Other species rows of the report kept, of 729 |
|---|---|--:|
| 300 | *Nocardiopsis dassonvillei* | 105 |
| `auto` (the cap, 6,500) | *N. dassonvillei* | 3 |
| Without the cap (51,851) | *N. dassonvillei*, *Escherichia coli* | 1 |

- *N. dassonvillei* has 15 read pairs from Kraken2 in this library and 264 distinct
  minimizers, so every floor removes it from the report: the one true member a floor
  removed in these tests, and a very thin one. *E. coli* K-12 has 713 read pairs and
  9,450 distinct minimizers, a member that is there and that a floor growing without
  limit would remove from the report; the cap keeps it.
- **The final table keeps the same species at every floor**, the same as with the gate
  off: Bracken's read threshold at this depth (`-t 2593`) already leaves out both thin
  members, and the other 24 members, and *Lysobacter enzymogenes*, which is not in the
  mock (154,092 distinct minimizers), pass every floor.
- The library has 2×125 bp reads, shorter than the 2×150 bp the rule assumes (see
  [Limits of the rule](#limits-of-the-rule)), so the floor is relatively stricter here
  than on the benchmark; it still removed no member that Bracken keeps.

**Deep, diverse soils.** Seven forest soils from the Schönbuch, Germany (BioProject
PRJNA717813; the O, A and B horizons of a Cambisol, a Podsol and a Stagnosol), with 50
to 108 million read pairs after quality filtering, of which PlusPF `k2_pluspf_20260626`
classified only 4.0 to 8.6 million. `auto` set floors of 3,953 to 6,500, the cap in four
of them. Their composition is not known, so this shows how much `auto` removes, not
whether it is right. As for the mock, the gate was run on the existing reports and
Bracken was not run again.

- **Most species rows of the report go, few species of the table.** `auto` removed 59 to
  76 % of the species rows of each Kraken2 report (300: 9 to 13 %); those rows held 0.13
  to 0.22 % of the classified reads. Of the 39 to 69 species per soil in the final table
  with the gate off, it removed none in three soils and one to three in the other four,
  all of them kept at 300: *Kribbella* sp. NBC_01484 in four soils (425 to 856 read
  pairs from Kraken2, 2,300 to 3,703 distinct minimizers) and two *Bradyrhizobium*
  species in one.
- **Those species held up to 5.6 % of a table, and that is the most a table can lose.**
  Bracken had raised *Kribbella* sp. NBC_01484 from its 425 to 856 read pairs to
  67,000–172,000; the three species removed from one soil held 5.6 % of its table, the
  one species removed from each of the other three 1.1 to 2.6 %. A species that passes
  the gate never loses read pairs when its relatives are removed, so the loss is at most
  this share; how much of it moves to relatives and how much leaves the table was not
  measured.
- **The cap changed no species of the table.** Without it the four capped soils would get
  floors of 7,007 to 8,607, and 42 to 232 more species rows of their reports would go,
  none of them in the final table.

#### Limits of the rule

These hold with `auto` and with a number sized by hand:

- **It rests on limited data** (see [Why `auto` is not the default](#why-auto-is-not-the-default)).
- **Do not go above 1.5 per 1,000 pairs** when you size a number by hand; `auto` uses 1
  per 1,000. At 1.5 rare true species start to be removed (92 of 13,419 true calls), and
  from 2 per 1,000 some samples' F1 falls even at full depth. On the CAMI II sample
  measured at full depth (above), 1.5 per 1,000 would set the floor at 9,462 and remove a
  true species, *Aliivibrio wodanis* (7,263 distinct minimizers).
- **Do not go above 6,500**, by hand either, for the reason `auto` stops there (the cap,
  above).
- **Host reads.** `auto` leaves out the *Homo sapiens* read pairs only. When you size a
  number by hand, subtract the read pairs Kraken2 assigned to the host (for human, the
  *Homo sapiens* row of the report), or remove the host first
  ([Decontamination](decontamination.md#host-removal)).
- **Library layout.** The rule assumes paired 2×150 bp Illumina reads with little mate
  overlap and few duplicates. Shorter reads, single-end data, overlapping mates and
  duplicate-rich libraries give fewer distinct minimizers per read pair, so the same
  floor is stricter for them; use a lower number rather than `auto`.
- **The floor follows the depth of the whole sample, not the reads of each species.** In
  a sample dominated by one organism, the other species face a higher floor without
  having more reads of their own; an example is at the end of this page.
- **It changes the table more than 300 does, and read pairs leave the table.** On the
  test dataset above (1.7–3.5 million classified pairs), `auto` keeps 25–61 species per
  sample, against 66–104 at 300, and 2.7–15.1 % of the read pairs leave the table,
  against at most 0.09 % at 300:

    | | SRR39106521 | SRR39106522 | SRR39106523 | SRR39106524 |
    |---|--:|--:|--:|--:|
    | Classified read pairs | 2,268,774 | 3,458,847 | 1,673,810 | 2,177,429 |
    | *Homo sapiens* read pairs, left out | 1,843 | 1,721 | 1,191 | 1,451 |
    | Floor set by `auto` | 2,267 | 3,458 | 1,673 | 2,176 |
    | Species in `otu_table.tsv`, gate off | 99 | 80 | 116 | 120 |
    | Species, 300 | 84 | 66 | 94 | 104 |
    | Species, one floor of 1,673 for every sample | 41 | 33 | 45 | 68 |
    | Species, `auto` | 35 | 25 | 45 | 61 |
    | Genera, gate off → `auto` | 39 → 12 | 31 → 10 | 46 → 18 | 52 → 23 |
    | Read pairs held by the species `auto` removes | 29 % | 13 % | 33 % | 14 % |
    | Change in total read pairs, 300 | −0.01 % | +0.01 % | −0.01 % | −0.09 % |
    | Change in total read pairs, `auto` | −15.1 % | −11.1 % | −4.2 % | −2.7 % |

    - Of the read pairs the removed species held, 15 to 87 % went to species that
      stayed; the rest left the table, because no species of their genus, family or
      order was left to take them.
    - Most of what left belongs to one order, Sphingomonadales. In SRR39106521 and
      SRR39106522 every one of its 24 and 23 species fails (the strongest,
      *Novosphingobium pentaromativorans*, has 2,174 and 2,156 distinct minimizers
      against floors of 2,267 and 3,458), and the whole order, 16 and 11 % of the two
      tables, is dropped; the rest of each table changes by less than 1 %. Where one
      species of the order passes, it collects most of the order's reads: in
      SRR39106523 *N. pentaromativorans* goes from 40,710 to 283,996 read pairs, on
      about the same distinct-minimizer count as in the two samples that drop it.
    - The order seems to be a real member of these communities without a close
      reference genome, not a spurious call. That is a reading of two indirect signs,
      not a measured truth: the 16S rRNA gene fragments in the same reads, classified
      by Metaxa2 against SILVA, put 13–21 % of them in the order, and sylph finds no
      genome of the order in any sample. Kraken2 then spreads its reads over many named
      species, none with much species-specific sequence, and the floor removes them
      all.
    - `auto` also removes *Homo sapiens* from all four samples (605 to 1,743 distinct
      minimizers, below its sample's floor in all four); at 300 the row stays. A user
      who reads human contamination off `otu_table.tsv` loses it.
    - The 25–61 species left are still well above the 9–16 genomes sylph finds. Compare
      `otu_table.tsv` with the gate off and on before relying on `auto`, and do not take
      the estimate of dropped read pairs in `logs/kmer_evidence/{sample}.log` for the
      measured change: with `auto` it was 3.5 and 3.9 times too low in the two samples
      that lose a whole order, because it only sees the reads Kraken2 assigned, not
      those Bracken moves.

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

The conclusion for the design: judge every species against a floor that is not fitted to
the evidence of the other species in its sample, and treat the trend line as a *ranking
for the eye* rather than a test. The floor is one fixed number for every sample of a run
or, with `auto`, a number set by the read pairs Kraken2 classified in the sample, human
ones left out. With `auto` the same species can therefore get different verdicts in
samples of different depth. That is intended when the whole community was sequenced
deeper: a species that is really present then leaves more distinct minimizers (see
[Scaling the floor with depth](#scaling-the-floor-with-depth)). But the floor follows
the depth of the whole sample, not the reads of each species, so a species with the same
reads in a deeper sample, where it is a smaller share, faces a higher floor. In the test
dataset *Novosphingobium pentaromativorans* has 2,156 to 2,195 distinct minimizers in
all four samples, on 2,955 to 3,550 read pairs from Kraken2; `auto` keeps it in the two
samples with floors of 1,673 and 2,176 and removes it from the two with floors of 2,267
and 3,458. In the deepest of the four, 57 % of the classified pairs belong to one
species, *Achromobacter denitrificans*, which raises the floor for everything else
without adding to their evidence.
