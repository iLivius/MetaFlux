# Validation of the phylogeny defaults

The [Phylogeny](../amplicon/phylogeny.md) page quotes measured figures for how much the
module's choices matter and for how much a short amplicon fragment can and cannot
recover. This page is where those figures come from. Two experiments produced them —
a **sensitivity analysis on the real 16S test run** and a **full-gene check** against
reference sequences — and they are different kinds of evidence, so each is described
with what it can and cannot support. The trees, alignments and analysis scripts behind
every number are kept with the project's benchmark material and are available on
request.

!!! note "Terms"
    Patristic distance, Robinson–Foulds distance, UFBoot, likelihood mapping and the
    rest are defined in the [glossary](../amplicon/phylogeny.md#glossary) on the
    Phylogeny page.


!!! note "These experiments predate the v2.4.0 input-order change"
    Both were run before MetaFlux started sorting the phylogeny input by sequence and
    building the tree on one thread. The sequences are the same and the comparisons
    between settings are unaffected, but the alignment rows are in a different order, so
    a fresh run of today's default does **not** reproduce the "T01" tree below: it finds
    its own optimum on the same data (total tree length 16.2 rather than 18.6, seven
    long-branch tips either way, six of them the same). Read the tables as *"changing
    this setting moved the output by this much"*, which is what they were built to
    answer, not as a description of the tree your own run will produce.

## The sensitivity analysis on the real test run

**Design.** One real dataset — the 16S test set, 211 ASVs after contaminant filtering,
6 samples, V5–V7 — and the pipeline's own default run as the reference point. Then one
thing changed at a time, everything else held fixed (same ASV set, seed 42; 4 threads for IQ-TREE, RAxML-NG
clamped to its own recommendation of 1, FastTree single-threaded as in the pipeline):

| Axis | Variants compared to the default |
|---|---|
| Alignment strategy | L-INS-i instead of FFT-NS-2 |
| Alignment masking | trimAl `-automated1`, trimAl `-gappyout`, AliFilter defaults — applied offline to the pipeline's alignment, tree rebuilt with the identical IQ-TREE command |
| Backend and model | IQ-TREE with ModelFinder; FastTree GTR+Γ; FastTree JC (the model QIIME 2 and ampliseq actually end up using); RAxML-NG GTR+Γ; RAxML-NG with MOOSE |

Ten trees in total, plus one IQ-TREE run with support values and one likelihood-mapping
run on the default alignment. In statistics this is a one-factor-at-a-time **sensitivity
analysis**: it answers "how much do the outputs move when this choice changes?", not
"which choice is correct".

**What was measured, and why that measure.**

- *Topology* — normalized Robinson–Foulds distance, RF for short (Robinson & Foulds 1981): the share
  of internal branches two trees do not have in common, 0 for identical trees, 1 when
  they share none.
- *Distances* — Pearson correlation between the two trees' patristic distance matrices
  (the branch-length distance between every pair of tips; the cophenetic correlation of
  Sokal & Rohlf 1962), computed over the full matrix as `ape::cophenetic.phylo()`
  returns it. This is the quantity UniFrac and Faith's PD actually consume, so
  agreement here matters more than agreement on topology.
- *Scale* — total tree length, because Faith's PD is by definition a sum of branch
  lengths (Faith 1992): a method that shrinks the tree shrinks every PD value with it.
- *Downstream stability* — a per-sample PD proxy (total branch length of the subtree
  spanning the ASVs present in that sample), correlated between trees across the six
  samples.
- *Taxonomic coherence* — the fraction of SINTAX genera with at least two ASVs that form
  a single clade, under the same rooting for every tree. The only external check
  available without a reference phylogeny, and a rough one: it inherits classifier
  errors, and some genera are genuinely not monophyletic in 16S.
- *Support* — ultrafast bootstrap (Hoang et al. 2018) and SH-aLRT (Guindon et al. 2010),
  read at IQ-TREE's documented single-gene thresholds.
- *Signal* — likelihood mapping (Strimmer & von Haeseler 1997): draw 5,000 random
  four-taxon subsets and ask, for each, whether the data clearly favour one of the three
  possible arrangements.
- *Model selection* — ModelFinder (Kalyaanamoorthy et al. 2017) by BIC (Schwarz 1978).
  The models named here: GTR (Tavaré 1986), TPM3u (the three-substitution-type
  family of Kimura 1981), +G gamma rates (Yang 1994), +R FreeRate (Yang 1995; Soubrier
  et al. 2012).

**Results.** Every tree is compared with the pipeline default, T01 (FFT-NS-2 alignment,
IQ-TREE, GTR+F+G4, seed 42). Robinson–Foulds is normalized (0 = same topology). The
patristic correlation and the per-sample PD correlation are Pearson's *r*.

| Tree | One thing changed | RF vs default | Patristic *r* | Tree length | PD *r* (6 samples) | Genera monophyletic |
|---|---|---|---|---|---|---|
| T01 | — (default) | 0 | 1 | 18.61 | 1 | 11 / 19 |
| T02 | alignment: L-INS-i | 0.457 | 0.954 | 17.58 | 0.997 | 11 / 19 |
| T03 | mask: trimAl `-automated1` | 0.606 | 0.914 | 13.44 | 0.997 | 10 / 19 |
| T04 | mask: trimAl `-gappyout` | 0.620 | 0.912 | 13.39 | 0.997 | 10 / 19 |
| T05 | mask: AliFilter defaults | 0.389 | 0.942 | 17.46 | 0.9995 | 11 / 19 |
| T06 | model: ModelFinder (picked TPM3u+R5) | 0.125 | 0.998 | 15.02 | 0.9999 | 11 / 19 |
| T07 | backend: FastTree GTR+Γ | 0.413 | 0.861 | 15.30 | 0.993 | 12 / 19 |
| T08 | backend: FastTree JC+Γ | 0.389 | 0.851 | 14.66 | 0.989 | 12 / 19 |
| T09 | backend: RAxML-NG GTR+Γ | 0.577 | 0.885 | 17.06 | 0.999 | 11 / 19 |
| T10 | backend: RAxML-NG MOOSE (picked GTR+FU+R5) | 0.466 | 0.916 | 15.06 | 0.999 | 10 / 19 |

What the table says, in order of how much each choice moved the output:

- *Masking* is the only choice that changes the scale: both trimAl masks removed 24 % of
  the columns and with them 28 % of total tree length (18.61 → 13.4), and every
  per-sample PD value shrank with it. The two trimAl trees are nearly the same tree
  (RF 0.154, *r* 0.997 between them). AliFilter removed almost nothing and changed
  almost nothing. No mask improved taxonomic coherence (10 / 19 genera vs 11 / 19).
- *Backend* moves topology and distances more than any setting within a backend:
  FastTree's two models give essentially the same tree (RF 0.111, *r* 0.986 between
  T07 and T08), yet both sit farthest from the IQ-TREE default in patristic distance
  (*r* 0.85–0.86). RAxML-NG under the same GTR+Γ model differs from IQ-TREE in topology
  (RF 0.577) but less in distances (*r* 0.885).
- *Model selection* changes the least: ModelFinder's TPM3u+R5 tree is the closest to
  the default in topology (RF 0.125) and in distances (*r* 0.998), while its total
  length is 19 % shorter — the FreeRate model rescales the branches, so absolute PD
  values are not comparable across models even when the tree is.
- *Alignment strategy* (T02) moves topology a lot (46 % of internal branches differ)
  but distances little (*r* 0.954). One tree per setting cannot say how much of that is
  the tier and how much is ordinary run-to-run variation of the tree search; the
  full-gene check below, where seeds were replicated, puts numbers on both.
- *Per-sample PD* barely notices any of it: *r* ≥ 0.989 against the default for every
  tree, and the ranking of the six samples is preserved by the two IQ-TREE trees
  (Spearman 1.0) and nearly so by the rest (Spearman 0.83–0.94). With six samples this
  is a coarse check.

**Support and signal on the default alignment.** The support run (T11: default tree,
`-B 1000 --alrt 1000`) puts 208 internal branches on the tree; 77 (37 %) reach
UFBoot ≥ 95, 141 (68 %) reach SH-aLRT ≥ 80, and 67 (32 %) meet both — the documented
read for a single-gene tree. Likelihood mapping on 5,000 random quartets finds
77.9 % fully resolved, 9.5 % partly resolved and 12.7 % unresolved. Both numbers say the
same thing: a ~440-column V5–V7 alignment supports the broad structure of the tree and
leaves about a third of the fine splits undecided, which is the expected amount for a
fragment this short.

**Model selection, for the record.** ModelFinder (BIC) chose TPM3u+R5 on the default
alignment; RAxML-NG's MOOSE chose GTR+FU+R5. Both prefer a five-category FreeRate
description of among-site rate variation and they disagree only on the substitution
matrix, which is why the module keeps the pinned GTR+F+G4 default and makes selection
available rather than automatic (see the Phylogeny page).

**Long-branch screen on the default tree.** Pendant edges: median 4.6 × 10⁻³, 75th
percentile 0.052, longest 0.657. Flagging at 5× the median would mark 81 of 211 tips
and at 20× still 33, which is why the QC report uses the boxplot fence Q3 + 3×IQR
(about 0.21 here: 0.208 with R's default quantile definition, 0.216 with the one the QC
report uses): 7 tips flagged either way, each the only ASV of its lineage in the run.

**What it is not.** It is not an accuracy benchmark: there is no known true tree for
environmental ASVs, so nothing in it says which tree is *right*. It is one dataset, one
amplicon region, one tree per setting, and six samples (so the PD correlations rest on six
points). Run-to-run variation of the tree search was not part of the design, but it was
measured afterwards on the same 211 ASVs: two IQ-TREE searches that differed only in the
order of two input sequences (see *Reproducibility* on the Phylogeny page) gave trees with
RF 0.37 between them, patristic *r* 0.93, tree lengths 18.0 and 18.6, and per-sample PD
*r* 0.999. That is the noise floor for the table above. The RF column (0.13–0.62) sits
largely inside it, so topological differences between settings should not be read as
effects of the settings; the masking effect on tree length (−28 %) and the backend effect
on distances (*r* 0.85–0.89) lie outside it, and the ModelFinder result (*r* 0.998) lies
well within it. The masks were produced with the same
binaries and flags the pipeline would use, but outside the pipeline. Read its numbers as
*"this choice moved the output by this much on this data"*. The question it cannot answer
— how far any fragment tree is from the best tree the same organisms' whole genes can
give — is what the next experiment is for.

## The full-gene check

**The idea.** There is no true tree for environmental ASVs, but there is the standard
substitute used by Janssen et al. (2018): take organisms whose *whole* 16S gene is known,
cut the amplicon out of each gene *in silico*, build a tree from the fragments exactly as
the pipeline would, and ask how well it recovers the tree built from the whole genes. The
full-gene tree is not "the truth" either — it is a maximum-likelihood estimate from a
1.4 kb gene — but it is the best-informed tree those organisms can give, and it is the
right yardstick for the question "what does a 376 bp fragment lose?".

**Design.**

- **Reference genes:** 250 full-length sequences from SILVA 138.2 (median 1,444 bp;
  19 phyla, 18 bacterial and 1 archaeal), chosen to mirror the test run: up to 3 sequences for each of the 55 genera
  SINTAX found there, 1 per family for ASVs that only reached family (171 taxon-matched),
  the rest random SILVA prokaryotes (79: 78 bacteria, one archaeon). Every tree in the experiment has exactly these 250 tips.
- **Fragments:** cut out with MetaFlux's own in-silico PCR (two-pass cutadapt, 20% primer
  mismatch allowed) — **V5–V7** with the test run's primers 799F/1175R (median 376 bp,
  the same length the pipeline's probe reports for this pair) and, as a second region,
  **V4** with the EMP primers 515F-Y/806R-B (median 253 bp).
- **Reference tree:** full genes, MAFFT L-INS-i (2,623 columns), IQ-TREE GTR+F+G4,
  seed 42. Rebuilt with seeds 43 and 44 to measure its own search noise, plus one
  ModelFinder run and one support run on the full genes.
- **Fragment grid, V5–V7:** the same axes as the sensitivity analysis, now with **three
  seeds** for each backend whose search starts from random choices — IQ-TREE GTR+F+G4, IQ-TREE ModelFinder, RAxML-NG GTR+Γ
  (3 seeds each), FastTree GTR+Γ and JC (deterministic, 1 each), the L-INS-i alignment
  and the trimAl-masked alignment (3 seeds each) — and V4 at pipeline defaults (3 seeds).
- **Scored against the reference:** normalized Robinson–Foulds (RF) distance (topology), patristic correlation
  (distances), tree-length ratio (branch-length compression), and genus/family monophyly
  — this time against SILVA's own taxonomy, a genuine external check rather than SINTAX
  calls on ASVs.

**Results** (mean over seeds; reference = full gene, GTR+F+G4, seed 42):

| Tree | Robinson–Foulds (RF) distance to reference | Patristic r | Tree length ÷ reference | Genus monophyletic | Family monophyletic |
|---|--:|--:|--:|--:|--:|
| Reference rebuilt, seeds 43/44 *(seed-to-seed variation — the noise floor)* | 0.02–0.04 | 0.999–1.000 | 0.98 | 84% | 89% |
| Reference, ModelFinder (→ GTR+F+R7) | 0.06 | 0.999 | **0.75** | 84% | 89% |
| **V5–V7, pipeline default** (FFT-NS-2, IQ-TREE GTR+F+G4) | 0.58 | 0.856 | 1.00 | 56% | 67% |
| V5–V7, pipeline default with support (`-B 1000 --alrt 1000`) | 0.56 | 0.863 | 1.02 | 59% | 68% |
| V5–V7, IQ-TREE ModelFinder (→ GTR+F+I+R5) | 0.59 | 0.862 | **0.67** | 55% | 67% |
| V5–V7, RAxML-NG GTR+Γ | 0.58 | 0.831 | 1.02 | 57% | 64% |
| V5–V7, FastTree GTR+Γ | 0.57 | 0.864 | **0.65** | 56% | 64% |
| V5–V7, FastTree JC | 0.59 | 0.868 | **0.62** | 56% | 64% |
| V5–V7, L-INS-i alignment | 0.60 | 0.878 | 1.11 | 55% | 64% |
| V5–V7, trimAl `-gappyout` | 0.59 | 0.872 | 0.96 | 56% | 61% |
| V4, pipeline default | 0.60 | 0.867 | 0.75 | 59% | 72% |

Support (UFBoot ≥ 95 *and* SH-aLRT ≥ 80): **64%** of internal branches on the full gene,
**38%** on the V5–V7 fragment.

**Reading it.**

- **The fragment loses most of the topology and keeps a fair amount of the distances.**
  Every V5–V7 tree sits in one band — about 58% of internal branches differ from the
  full-gene tree (Robinson–Foulds distance 0.56–0.60), while patristic distances correlate at 0.83–0.88. For
  scale: the full-gene reference rebuilt with another seed differs from itself by 0.02–0.06, and the same *fragment* alignment rebuilt with another seed by 0.22–0.28 —
  so roughly a third of the fragment's disagreement is the tree search wandering, and
  the rest is information the fragment does not have. That band *is* the framing paragraph at the top of
  the [Phylogeny](../amplicon/phylogeny.md) page, measured.
- **None of the settings tested escapes the band.** Aligner, mask, backend, model: a
  spread of 0.04 in RF and 0.05 in patristic r — smaller than the seed-to-seed variation
  of a single setting — against a gap of ~0.55 to the full gene. The pipeline's choices
  are not where the distance to the full-gene tree is decided; the fragment is. The
  defaults stand.
- **V4 is not a better region than V5–V7** on this measure (RF 0.60, r 0.87).
- **Roughly a third of the taxonomic coherence is lost with the fragment** (a fifth to
  two fifths, depending on the measure): genera monophyletic 84% → 55–59%, families
  89% → 61–72%; well-supported branches 64% → 38%.
- **Branch-length scale depends on the rate model and on the region — and neither
  tells you the true scale.** V5–V7 fragment trees under GTR+F+G4 come out the same
  total length as the full-gene GTR+F+G4 tree (ratio ≈ 1.0); the V4 fragment under the
  same model does not (0.75). FreeRate models and FastTree make the V5–V7 tree 33–38%
  shorter — and the 1.4 kb *reference itself* loses 25% when ModelFinder picks GTR+F+R7
  for it, so on this data the shrinkage is a rate-model effect, not a fragment effect.
  (Medlar et al. 2014 describe a separate short-fragment underestimation; the V4 row is
  the only place this grid shows something like it.) Either way: absolute PD must never
  be compared between runs built under different models or from different regions.
- **Model selection at 250 references** lands in the GTR family (GTR+F+R7 on the full
  gene, GTR+F+I+R5 on the fragment), with FreeRate preferred over gamma both times — the
  same rate-model preference the sensitivity analysis found on the real ASVs.

**Limits.** One gene set of 250, one taxonomic mix (built to resemble one plant-associated
community), one region pair; the reference is itself an estimate (only 64% of its
branches are well supported, so some of the "disagreement" is the reference's own
uncertainty); and the fragments are cut from clean full-length references, so they carry
no sequencing error — a real ASV set adds sequencing error on top of the
fragment's own information loss, so there is no reason to expect it to do better. What this experiment
does establish, and the sensitivity analysis could not, is the *size* of the gap between every fragment tree in this grid and the full-gene
tree, and that it dwarfs every pipeline choice tested.

## References

The papers behind every measure and model named here — Robinson & Foulds 1981, Sokal &
Rohlf 1962, Faith 1992, Strimmer & von Haeseler 1997, Hoang et al. 2018, Guindon et al.
2010, Schwarz 1978, Tavaré 1986, Kimura 1981, Yang 1994 and 1995, Soubrier et al. 2012,
Kalyaanamoorthy et al. 2017, Janssen et al. 2018, Medlar et al. 2014, Apprill et al. 2015
— are listed on the [Citation](citation.md) page.
