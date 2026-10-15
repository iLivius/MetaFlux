# Choosing confidence and threshold

Two Kraken2/Bracken settings affect how many false species calls reach the final table:
the Kraken2 confidence score and the Bracken read threshold. In the benchmark on this
page, changing confidence had little effect on accuracy. The threshold had the largest
effect of any shotgun setting tested, and a fixed value suits only one sequencing depth:
at Bracken's default of 10, accuracy fell the more deeply a sample was sequenced.

- `shotgun.kraken.confidence` — **leave it at 0.15.**
- `shotgun.bracken.threshold` — **set it to `auto`** (the default).

This page is the benchmark behind that recommendation: what was tested, how, and what it
did and did not show. If you only want the setting, skip to
[The recommendation](#the-recommendation).

## Effect on a mock community

Kraken2 classifies a read to the taxon its k-mers best match, however thin that match
is. Bracken then re-estimates abundance from those calls. Neither step, by default, asks
*how much* evidence a taxon actually has before trusting it — Bracken's own floor
(`-t`, ten reads) has been ten reads since its first release, at any sequencing depth.

Species-level results on a mock community of known composition, the ZymoBIOMICS Gut
Microbiome Standard (D6331), SRA run
[SRR38976715](https://www.ncbi.nlm.nih.gov/sra/SRR38976715) (BioProject PRJNA1474435),
classified with PlusPF `k2_pluspf_20260626` at Kraken2 confidence 0.15, the v2.5.0
defaults:

| | Species reported | Precision | Recall | **F1** |
|---|--:|--:|--:|--:|
| Bracken's default, `-t 10` | 54 | 0.278 | 0.882 | 0.423 |
| MetaFlux default, `threshold: auto` (`-t 265`) | 20 | 0.700 | 0.824 | **0.757** |

The mock's *Fusobacterium*, listed by Zymo as *F. nucleatum*, is scored as
*F. animalis*, so calls of *F. nucleatum*, *F. polymorphum* and *F. vincentii* count as
false positives ([Benchmark data](../about/test-datasets.md#benchmark-data) gives the
reason).

Recall fell slightly, from 0.88 to 0.82: *Salmonella enterica*, with 135 Kraken2 reads,
fell below the threshold. Precision rose from 0.28 to 0.70. At `-t 10`, 39 of the 54
species reported were not in the mock. Bracken reports a species only when Kraken2 gave
it at least `-t` reads, so each of them had at least ten Kraken2 reads. Most of them
(29 of the 39) are in the same genus as a member of the mock and share conserved
regions with it; database cross-contamination is another source. Five of the six false
species left at `-t 265` are *Fusobacterium* species closely related to the mock's
*F. animalis*; the sixth is human.

## Kraken2 confidence — keep it at 0.15

Raising the confidence score removes far fewer weakly supported species than the Bracken
threshold does, because Kraken2 does not reject a read that falls short of the
confidence threshold: it moves the call **up** the taxonomy tree until the score clears
the bar, and only reports "unclassified" if even the root fails. Raising confidence
therefore mostly moves reads to coarser ranks, which leaves more reads for Bracken to
redistribute.

The measured consequence: raising confidence lowers the read count of species-level
calls, and a rare species can lose a large share of its few reads. On the same D6331
reads, *Salmonella enterica*, 0.01 % of the mock's DNA, had 135 Kraken2 reads at
confidence 0.15 and 85 at 0.20 (PlusPF `k2_pluspf_20260626`; 134 and 87 with
`k2_pluspf_20251015`).

Raising confidence and applying a read-count floor both reduce the same small read
counts, so together they remove rare species that are present. [Wright, Comeau & Langille
(2023)](https://doi.org/10.1099/mgen.0.000949) swept confidence 0.00–1.00 on simulated
data: 0.60 gave the highest mean F1, while 0.15 gave the lowest L1 distance, that is,
the most accurate relative abundances. 0.15 is the source of MetaFlux's default, and
nothing here changes it.

## Bracken threshold

### Why a fixed number is wrong

The right floor depends on how many reads a sample has. A species with 200 reads is
better supported in a sample of 2 million classified read pairs than in one of
15 million: `auto` sets `-t 100` in the first, which keeps it, and `-t 750` in the
second, which removes it. A fixed threshold is therefore right for one depth only. In
the first benchmark behind the default, with PlusPF `k2_pluspf_20251015`, one CAMI II
marine sample (sample 0) was subsampled to 2, 5 and 10 million read pairs and also used
at its full 16.6 million. At Bracken's default of 10, F1 was lower the deeper the sample
was sequenced:

| classified pairs | F1 at fixed `-t 10` |
|--:|--:|
| 1.8 M | 0.740 |
| 4.5 M | 0.667 |
| 9.0 M | 0.612 |
| 15.0 M | 0.576 |

With the community and all other settings unchanged, F1 fell from 0.740 to 0.576 (22 %)
as depth increased. A depth-scaled threshold over the same range stayed at 0.82–0.83.

### The rule

```
threshold = max(threshold_min, threshold_alpha × classified_read_pairs)
```

**Units matter: classified read *pairs*, not raw reads, not trimmed reads.** Kraken2 in
paired mode counts a pair as one unit, and the classified fraction depends on how well
the database covers the community: with PlusPF `k2_pluspf_20260626`, 99.8 % of the
D6331 mock's read pairs were classified, against 47–75 % for the four wastewater
metagenomes of the [shotgun test dataset](../about/test-datasets.md#shotgun-test-dataset),
at similar depth (5.3 million read pairs entering Kraken2 for the mock, 3.5–4.7 million
for the wastewater metagenomes). Raw input reads would therefore be the wrong
denominator.

```yaml
shotgun:
  bracken:
    threshold: auto            # or a fixed integer
    threshold_alpha: 5.0e-05
    threshold_min: 10
```

With `auto`, MetaFlux reads the classified-pair count from Kraken2's own report after it
finishes — the count does not exist until then — and computes the threshold per sample.
The decision is written to that sample's Bracken log; for the D6331 run it reads:

```
[bracken] auto threshold: 5e-05 x 5292994 classified read pairs -> -t 265
```

### Where alpha comes from, and what is and is not new here

`-t` is Bracken's own parameter, and Bracken does not scale it: it takes a flat integer
and does nothing else with it. A cut-off that grows with depth is not a new idea.
Benchmarks usually filter on relative abundance, a fixed share of the reads (for example
[Ye et al. 2019](https://www.sciencedirect.com/science/article/pii/S0092867419307755)),
and KrakenUniq found that the best unique-k-mer threshold rose roughly linearly with
sequencing depth ([Breitwieser et al. 2018](https://doi.org/10.1186/s13059-018-1568-0)).

MetaFlux adds a measured constant and computes the threshold for each sample
automatically:

- a depth exponent of **0.967**, fitted in the first benchmark on one CAMI II marine
  sample classified at four depths (1.8 to 15 million classified pairs, the series
  above), consistent with a threshold directly proportional to depth
- a coefficient, **alpha = 5×10⁻⁵**, chosen in the first benchmark (PlusPF
  `k2_pluspf_20251015`) and re-fitted on all 41 CAMI II samples with
  `k2_pluspf_20260626`, where the median optimum was 6.0×10⁻⁵; on the D6331 mock the
  optimum is several times higher (next section)
- the `auto` computation itself, run per sample inside the `bracken` rule

None of this is standard terminology. "Alpha" is the name given here to the
proportionality constant; if you cite this elsewhere, describe it as a depth-scaled
Bracken threshold with a benchmark-fitted coefficient, not as a published Bracken feature.

### How firmly alpha is established

| Community (samples) | Species per sample | Optimal alpha |
|---|--:|--:|
| ZymoBIOMICS D6331, run SRR38976715 | 17 | 2.3×10⁻⁴ – 3.5×10⁻⁴ |
| CAMI II strain-madness, samples 0–9 | 20 | 4.1×10⁻⁵ – 1.2×10⁻⁴ |
| CAMI II marine, samples 0–9 | 256–381 | 5.2×10⁻⁵ – 2.1×10⁻⁴ |
| CAMI II rhizosphere, samples 0–20 | 73–126 | 1.5×10⁻⁵ – 2.6×10⁻⁴ |

The optimal alpha is the `-t` with the highest species F1, divided by the sample's
classified read pairs. D6331: the run in [Effect on a mock
community](#effect-on-a-mock-community) (PlusPF `k2_pluspf_20260626`, confidence 0.15),
where F1 reaches 0.903 for any `-t` from 1,191 to 1,832. CAMI II rows: the range of the
per-sample optima in the 41-sample re-fit with PlusPF `k2_pluspf_20260626`.

Pooling the 41 CAMI II samples gives a median optimal alpha of **6.0×10⁻⁵**, 1.2× the
default. One constant works for these samples despite the spread because F1 changes
little over a wide range of alpha: using 5×10⁻⁵ everywhere costs a mean of **0.023 F1**
against tuning each sample to its own best (median cost 0.014). In the first benchmark,
Bracken's default of 10 cost up to **0.34 F1**.

The D6331 mock is an exception: there the default costs 0.146 F1 against the optimum
(0.757 against 0.903).

!!! note "Adjusting alpha for a known community"

    A smaller `threshold_alpha` is a looser filter, a larger one a stricter filter. The
    best value differs between communities and between samples of one community (table
    above), and it does not follow the number of species: the optimum of the 17-species
    D6331 mock is 4.5 to 7 times the default, while the optima of the 20-species
    strain-madness samples range from 0.8 to 2.4 times it. The default is recommended for
    routine use.

## The recommendation

```yaml
shotgun:
  kraken:
    confidence: 0.15            # unchanged; do not raise it alongside the threshold
  bracken:
    threshold: auto
    threshold_alpha: 5.0e-05
    threshold_min: 10
```

Worked examples: 2M classified pairs → `-t` ≈ 100 · 5M → ≈ 250 · 15M → ≈ 750.

## What this does — and does not — fix

The threshold raises precision substantially, lowers recall slightly and changes
abundance accuracy little.

On the D6331 run, the L1 distance between the Bracken profile and the genomic-DNA
composition in Zymo's manual (the sum of absolute differences, in percentage points) was
40.88 at `-t 10` and 40.87 at `-t 265`, and the Bray–Curtis dissimilarity between the
two profiles was 0.0009. The threshold changes **which taxa are named**, not **how much of
each is reported**. If your downstream work is ordination or differential abundance on
already-good tables, expect little change. If it is species counts, richness, or
presence/absence, expect a large one.

**It cannot rescue a database that does not cover your community.** Only 165 of the 307
true species of the CAMI II rhizosphere (plant-associated) set, samples 0–20, are in
PlusPF at all: 54 % over the whole set and 42 % per sample on average, the same in
`k2_pluspf_20251015` and `k2_pluspf_20260626`. No parameter setting recovers the rest,
and species F1 stays near 0.3 at every setting tested. The gap lies in the reference
database, and a database update alone should not be expected to close it.

## The benchmark behind these numbers

**Real data with known composition:** the [ZymoBIOMICS Gut Microbiome
Standard (D6331)](https://files.zymoresearch.com/protocols/_d6331_zymobiomics_gut_microbiome_standard.pdf),
SRA run SRR38976715 (BioProject PRJNA1474435; see
[Benchmark data](../about/test-datasets.md#benchmark-data)), 17 species after
strain-collapse, with its *Fusobacterium* scored as *F. animalis*; classified with PlusPF
`k2_pluspf_20251015` and `k2_pluspf_20260626`.

**Simulated data with gold-standard truth:** [CAMI
II](https://www.nature.com/articles/s41592-022-01431-4) marine samples 0–9, rhizosphere
(plant-associated) samples 0–20 and strain-madness samples 0–9 (every strain-madness
sample holds the same 20 species by design), 41 samples in all. All 41 were classified
with PlusPF `k2_pluspf_20260626`; five of them were also classified with
`k2_pluspf_20251015` for the first choice of alpha.

**Real data with no known composition:** the four wastewater metagenomes of the
[shotgun test dataset](../about/test-datasets.md#shotgun-test-dataset) (BioProject
PRJNA1477390, runs SRR39106521–SRR39106524). Without a truth table they give no accuracy
figures; this page uses them for the classified fractions under [The rule](#the-rule),
and [What it does on real data](kmer-evidence.md#what-it-does-on-real-data) shows the
k-mer evidence gate on them.

**Scorer check:** the benchmark scorer was checked against
[OPAL](https://github.com/CAMI-challenge/OPAL) 1.0.13, the CAMI assessment tool, on the
Bracken species tables of the two runs with PlusPF `k2_pluspf_20260626` (Kraken2
confidence 0.15 and 0.20), which cover the 41 CAMI II samples and the D6331 mock. With
its inputs prepared as the scorer prepares them, OPAL's species-level true positives,
false positives and false negatives equal the scorer's in every sample. Run on the gold
standards as distributed, OPAL differs only in the taxa those preparation rules touch:
gold taxids that NCBI has merged since the CAMI release, placeholder names such as
"uncultured bacterium", which the scorer does not score, and *F. animalis*, which the
scorer counted as the mock's *F. nucleatum* in this check; the D6331 figures on this
page score *F. animalis* as the true member instead. The first benchmark, the 41-sample
re-fit behind the CAMI II rows of the alpha table and the `k2_pluspf_20251015` runs
were not part of this check.

**Comparator:** [MetaPhlAn 4](https://github.com/biobakery/MetaPhlAn), database
`mpa_vJan26_CHOCOPhlAnSGB_202605`, on the same trimmed reads of the D6331 run. MetaPhlAn
reported 14 of the 17 species and no false ones (F1 0.90). MetaFlux at its v2.5.0
defaults (PlusPF `k2_pluspf_20260626`) reported the same 14 species and six false ones,
five of them *Fusobacterium* species related to the mock's *F. animalis* (F1 0.76).

## References

- Lu, J. et al. (2017). Bracken: estimating species abundance in metagenomics data.
  *PeerJ Computer Science* 3:e104. <https://peerj.com/articles/cs-104/>
- Ye, S. H. et al. (2019). Benchmarking Metagenomics Tools for Taxonomic Classification.
  *Cell* 178:779–794. <https://www.sciencedirect.com/science/article/pii/S0092867419307755>
- Breitwieser, F. P., Baker, D. N. & Salzberg, S. L. (2018). KrakenUniq: confident and
  fast metagenomics classification using unique k-mer counts. *Genome Biology* 19:198.
  <https://doi.org/10.1186/s13059-018-1568-0>
- Wright, R. J., Comeau, A. M. & Langille, M. G. I. (2023). From defaults to databases:
  parameter and database choice dramatically impact the performance of metagenomic
  taxonomic classification tools. *Microbial Genomics* 9(3).
  <https://doi.org/10.1099/mgen.0.000949>
- Meyer, F. et al. (2022). Critical Assessment of Metagenome Interpretation — the second
  round of challenges (CAMI II). *Nature Methods* 19:429–440.
  <https://www.nature.com/articles/s41592-022-01431-4>
- Meyer, F. et al. OPAL: taxonomic profiling assessment.
  <https://github.com/CAMI-challenge/OPAL>
- Zymo Research. ZymoBIOMICS Gut Microbiome Standard (D6331), Instruction Manual v1.2.0.
  <https://files.zymoresearch.com/protocols/_d6331_zymobiomics_gut_microbiome_standard.pdf>

See also [Species-level k-mer evidence](kmer-evidence.md) for the optional gate that sits
between Kraken2 and Bracken — the same benchmark found it contributes little once the
threshold above is scaled correctly, but it catches a narrow class of error the threshold
cannot.
