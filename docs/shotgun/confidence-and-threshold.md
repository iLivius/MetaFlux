# Choosing confidence and threshold

Shotgun mode has two Kraken2/Bracken settings that decide how much of the final table is
real. One of them barely matters. The other is the single biggest accuracy lever in the
whole shotgun path, and MetaFlux's own default for it was wrong until this was measured.

- `shotgun.kraken.confidence` — **leave it at 0.15.**
- `shotgun.bracken.threshold` — **the one that matters.** Set it to `auto`.

This page is the benchmark behind that recommendation: what was tested, how, and what it
did and did not show. If you only want the setting, skip to
[The recommendation](#the-recommendation).

## The problem in one picture

Kraken2 classifies a read to the taxon its k-mers best match, however thin that match
is. Bracken then re-estimates abundance from those calls. Neither step, by default, asks
*how much* evidence a taxon actually has before trusting it — Bracken's own floor
(`-t`, ten reads) has been ten reads since its first release, at any sequencing depth.

Run a real Kraken2/Bracken pipeline against a mock community of known composition and
the effect is not subtle:

| | Species reported | Precision | Recall | **F1** |
|---|--:|--:|--:|--:|
| MetaFlux defaults | 47 | 0.319 | 0.882 | 0.469 |
| MetaFlux, `threshold: auto` | 15 | 0.933 | 0.824 | **0.875** |

Recall — finding the organisms that are really there — was never the problem: 0.88 either
way. **Precision was catastrophic**, because roughly two of every three species named
were phantoms: conserved regions shared with an abundant relative, database
cross-contamination, or reads that Kraken2 placed at genus level and Bracken then pushed
down onto whichever species happened to be listed first.

## Kraken2 confidence — keep it at 0.15

Confidence is a stringency knob, and raising it looks like the obvious fix. It is not,
for a mechanistic reason worth understanding before touching it: Kraken2 does not reject
a read that falls short of the threshold — it moves the call **up** the taxonomy tree
until the score clears the bar, and only reports "unclassified" if even the root fails.
Raising confidence therefore mostly pushes reads to coarser ranks, handing Bracken more
work, not less.

The measured consequence: raising confidence **shrinks the read count of every taxon**,
rare ones hardest. On a real mock community, *Salmonella enterica* — present at 0.01% —
had:

| confidence | 0.05 | 0.15 | 0.20 | 0.30 | 0.50 |
|---|--:|--:|--:|--:|--:|
| reads | 362 | 134 | 87 | 34 | 10 |

Push confidence up **and** apply a read-count floor, and the two filters compound against
exactly the rare true organism you most want to keep. That is how a plausible-looking
tuning choice quietly deletes real biology. [Wright, Comeau & Langille
(2023)](https://doi.org/10.1099/mgen.0.000949) swept confidence 0.00–1.00 on simulated
data and found two F1 optima, 0.60 and 0.15, with 0.15 favouring accurate relative
abundance — the source of MetaFlux's existing default. Nothing here changes it.

## Bracken threshold — the game changer

### Why a fixed number is wrong

The right floor depends on how many reads you have. Fifty reads is solid evidence in a
2-million-pair sample and noise in a 15-million-pair one. A **fixed** threshold therefore
gets worse the deeper you sequence:

| classified pairs | F1 at fixed `-t 10` |
|--:|--:|
| 1.8 M | 0.740 |
| 4.5 M | 0.667 |
| 9.0 M | 0.612 |
| 15.0 M | 0.576 |

Same community, same everything else — F1 falls by a quarter purely because the sample
was sequenced deeper. A depth-scaled threshold, over the identical range, held flat at
~0.82–0.83.

### The rule

```
threshold = max(threshold_min, threshold_alpha × classified_read_pairs)
```

**Units matter: classified read *pairs*, not raw reads, not trimmed reads.** Kraken2 in
paired mode counts a pair as one unit, and the classified fraction can differ 2× between
samples of similar depth depending on how well the database covers that community — using
raw input reads as the denominator would be a different, wrong number.

```yaml
shotgun:
  bracken:
    threshold: auto            # was: a fixed integer
    threshold_alpha: 5.0e-05
    threshold_min: 10
```

With `auto`, MetaFlux reads the classified-pair count from Kraken2's own report after it
finishes — the count does not exist until then — and computes the threshold per sample.
The decision is written to that sample's Bracken log:

```
[bracken] auto threshold: 5e-05 x 5294079 classified read pairs -> -t 265
```

### Where alpha comes from, and what is and is not new here

`-t` is Bracken's own parameter. The *principle* that it should scale with depth is
already published: [Ye et al. (2019)](https://www.sciencedirect.com/science/article/pii/S0092867419307755)
state that a 10× deeper sample should use a roughly 10× higher threshold for comparable
results. Bracken's own source contains no such logic — it takes a flat integer and does
nothing else with it.

What MetaFlux adds is the measured constant and an implementation that applies it without
you having to think about it:

- a fitted depth exponent of **0.967** (marine community, four depths spanning 8×) —
  confirmation that "scale with depth" really does mean direct proportionality, not
  something looser
- a coefficient, **alpha = 5×10⁻⁵**, fitted across every benchmarked community and
  confirmed stable across **two Kraken2 database builds eight months apart**
- the `auto` computation itself, run per sample inside the `bracken` rule

None of this is standard terminology. "Alpha" is simply the name given here to the
proportionality constant; if you cite this elsewhere, describe it as a depth-scaled
Bracken threshold with a benchmark-fitted coefficient, not as a published Bracken feature.

### How firmly alpha is established

| Community | Species | Optimal alpha |
|---|--:|--:|
| ZymoBIOMICS D6331 mock | 17 | 2.4×10⁻⁵ |
| CAMI II strain-madness | 20 | 4.1–1.2×10⁻⁴ (10 samples) |
| CAMI II marine | 256–478 | 5.2–2.1×10⁻⁴ (10 samples) |
| CAMI II rhizosphere | 73–307 | 1.5×10⁻⁵–2.6×10⁻⁴ (21 samples) |

Pooling every sample from every community and refitting from scratch gives a median
optimal alpha of **6.0×10⁻⁵** — a 1.2× shift from the value shipped, on an independent,
8× larger set of samples. The reason one fixed constant still works despite that spread
is that the optimum is a **broad plateau, not a sharp peak**: across 41 CAMI II samples,
using 5×10⁻⁵ everywhere costs a mean of **0.023 F1** against tuning each sample to its own
individual best (median cost 0.014), while leaving the threshold at the historic default
of 10 costs up to **0.34**. Being off by 2× on alpha is nearly free; not scaling at all is
not.

!!! note "If you know your community, a rough refinement is possible — rarely worth it"

    Simple, defined mixtures tolerate a smaller alpha (looser filter); complex,
    hundred-plus-species communities want a larger one (stricter filter). The gain from
    hand-picking is small enough, and the plateau wide enough, that the universal default
    is recommended for routine use.

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

**Precision, substantially. Recall, not touched. Abundance accuracy, not at all.**

L1 distance to the true composition and Bray–Curtis dissimilarity were measured flat
across every threshold tested, on every dataset — both real and simulated. This fixes
**which taxa are named**, not **how much of each is reported**. If your downstream work
is ordination or differential abundance on already-good tables, expect no change. If it
is species counts, richness, or presence/absence, expect a large one.

**It cannot rescue a database that does not cover your community.** CAMI II's
rhizosphere (soil-associated) set has only ~52% of its true species present in the
PlusPF database used here **at all** — no parameter setting recovers them, and F1 caps
near 0.3 regardless of tuning. This was checked against two PlusPF builds eight months
apart (2025-10 and 2026-06): the newer, larger database did not move the ceiling
(0.521 → 0.518). Environmental coverage gaps are a reference-database problem, not
a parameter problem, and a database refresh alone should not be expected to close them.

## The benchmark behind these numbers

**Real data with known composition:** the [ZymoBIOMICS Gut Microbiome
Standard (D6331)](https://files.zymoresearch.com/protocols/_d6331_zymobiomics_gut_microbiome_standard.pdf),
17 species after strain-collapse.

**Simulated data with gold-standard truth:** [CAMI
II](https://www.nature.com/articles/s41592-022-01431-4) marine (10/10 samples), rhizosphere
/ plant-associated (21/21), and strain-madness (10/100 — every sample of that set shares
the same 20 species by design, so ten are as informative as all of them). 41 samples,
630+ parameter combinations, two independent Kraken2 PlusPF builds.

**Real data with no truth, used for concordance only:** four wastewater
naproxen-degrading enrichment cultures.

**Independent scorer cross-check:** every CAMI-derived precision/recall/F1 figure was
re-scored with [OPAL](https://github.com/CAMI-challenge/OPAL), the CAMI-challenge scoring
tool, against the official gold standards. Every value matched the pipeline's own scorer
to at least four significant figures.

**Comparator:** [MetaPhlAn 4](https://github.com/biobakery/MetaPhlAn) scored on the same
samples. It out-performs tuned MetaFlux on the gut-derived Zymo mock (F1 0.90 vs 0.88) —
its marker database is built for exactly that niche — but **tuned MetaFlux outperforms it
on the CAMI marine community** (F1 0.83 vs 0.81), with better recall throughout. Neither
tool is uniformly better; which wins depends on how well the reference database matches
the sample.

## References

- Lu, J. et al. (2017). Bracken: estimating species abundance in metagenomics data.
  *PeerJ Computer Science* 3:e104. <https://peerj.com/articles/cs-104/>
- Ye, S. H. et al. (2019). Benchmarking Metagenomics Tools for Taxonomic Classification.
  *Cell* 178:779–794. <https://www.sciencedirect.com/science/article/pii/S0092867419307755>
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
