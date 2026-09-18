# Phylogeny (16S)

Optional. Builds an approximate, same-region 16S marker-fragment phylogeny from the
eligible ASVs, for exploratory phylogenetic-diversity analysis (e.g. Faith's PD,
UniFrac) performed in your own downstream tools. It is not a species tree or a resolved
deep phylogeny: short single-region fragments carry limited deep signal, so treat
topology beyond close relatives, and any single long branch, with caution.

Off by default. Turn it on with:

```yaml
amplicon:
  phylogeny:
    enabled: true
```

## Glossary

The terms this page leans on, in plain words. Skip it if you build trees for a living.

| Term | What it means here |
|---|---|
| **Alignment** | The ASV sequences written one under the other, with gaps inserted so that equivalent positions line up in the same column. Every tree program reads this, not the raw sequences. |
| **ASV** | Amplicon sequence variant — DADA2's unit of output: one exact denoised sequence. Each ASV is one **tip** of the tree. |
| **Backend** | The program that builds the tree from the alignment: IQ-TREE, FastTree or RAxML-NG. |
| **Branch length** | The length of one line segment in the tree, in expected substitutions per site — how much sequence change separates its two ends. Not time. |
| **BIC** | Bayesian information criterion. A score for choosing between models: lower is better, and it charges a penalty for every extra parameter, so a model has to earn its complexity. |
| **Faith's PD** | Phylogenetic diversity: the total branch length of the part of the tree spanned by the ASVs present in a sample. As usually computed (`picante::pd`, default `include.root = TRUE`) the path up to the root is included too, which is why the root position matters. A sum of branch lengths, so anything that scales the tree scales PD. |
| **FFT-NS-2 / L-INS-i** | MAFFT's two alignment strategies used here. FFT-NS-2 is fast: it builds the alignment in one pass, adding sequences group by group along a rough guide tree. L-INS-i then repeatedly re-aligns to improve it — more accurate, but documented only up to ~200 sequences. |
| **Γ, +G, +G4, +R, +I** | Ways of letting different alignment positions evolve at different speeds (in rRNA, stems barely change, loops change fast). **+G4** / **Γ**: speeds drawn from a gamma distribution in 4 categories (RAxML-NG and FastTree write it as +G or Γ; FastTree's version uses 20 categories). **+R5**: five speed categories estimated freely (FreeRate). **+I**: a fraction of positions that never change. |
| **In-silico PCR** | Cutting the amplicon out of a full-length reference gene on the computer, by finding the two primer sites. MetaFlux does this to size amplicons; the full-gene check below does it to make fragments of known origin. |
| **Likelihood mapping** | A test of how much phylogenetic signal an alignment carries: draw many random **quartets** (sets of four sequences) and ask, for each, whether the data clearly favour one of the three possible ways to connect them. |
| **Long branch** | A tip joined to the rest of the tree by an unusually long branch — what a contaminant, an off-target sequence, or a genuinely isolated lineage looks like. Inflates PD, distorts unweighted UniFrac. |
| **Maximum likelihood (ML)** | The way all three backends choose a tree: the tree and branch lengths under which the observed sequences are most probable, given a substitution model. FastTree is an *approximate* ML method. |
| **Midpoint rooting** | Placing the root halfway along the longest tip-to-tip path. Convenient, but it can move when tips are pruned — whenever the longest path changes — and it is thrown off by a single long branch. |
| **Model (substitution model)** | The assumed rules of sequence change. **JC** — every change equally likely; **GTR** — each of the six kinds of change (A↔C, A↔G, …) has its own rate; **TPM3u** — a GTR relative in which those six kinds are grouped into three rate parameters; **+F** — base frequencies counted from the data. The site-speed parts (+G4, +R5, +I) are the row above. |
| **ModelFinder / MOOSE / MFP** | Automatic model choice in IQ-TREE (ModelFinder, requested with `MFP`) and in RAxML-NG (MOOSE, requested with `DNA`). Both try many models and keep the best by BIC. |
| **Monophyletic** | A group of tips that sit together on the tree as one complete clade, with nothing else mixed in. "Are the ASVs of the same genus monophyletic?" is a rough check of whether the tree reflects known taxonomy. |
| **Newick** | The text format trees are written in: nested parentheses with branch lengths, e.g. `((ASV_1:0.01,ASV_2:0.02):0.05,ASV_3:0.03);`. |
| **Outgroup** | A lineage known to sit outside the group of interest. Used to root a tree; a paralogue such as *parE* in a gyrB dataset behaves like one whether you want it or not. |
| **Outlier fence (Q3 + 3×IQR)** | The standard boxplot rule for calling a value an outlier: the upper quartile plus three times the interquartile range. Used by the long-branch report because it does not care how small the median is. |
| **Patristic distance** | The distance between two tips measured *along the tree*: add up the branch lengths on the path from one to the other. It is built from the same branch lengths UniFrac and PD use, so two trees whose patristic distances agree give near-identical diversity numbers — which is why it is the comparison that matters most on this page. |
| **Pendant edge** | The single branch that connects a tip to the rest of the tree. Its length is the most direct measure of how isolated that ASV is. |
| **Pruning** | Removing tips from a tree. Where a removed tip's parent is left with just one branch in and one out (a "degree-2 node"), that node is dissolved by adding its two branch lengths together, so distances among the remaining tips do not change. |
| **Robinson–Foulds (RF) distance** | How different two tree shapes are: the number of internal branches present in one tree but not the other. Normalized to 0 (identical) – 1 (nothing in common). Compares only the branching pattern, not branch lengths. |
| **Root / unrooted** | A rooted tree has a designated oldest point; an unrooted one has only the branching pattern and lengths. MetaFlux exports unrooted trees, because the right root depends on which tips you keep. |
| **Root-to-tip distance** | The sum of branch lengths from a tip up to the root. Needs a root, so the QC report measures it on a temporary midpoint-rooted copy; the exported tree stays unrooted. Good at spotting a whole clade hanging off one long **stem** (the branch leading to a clade), poor at spotting a single long branch. |
| **Saturation** | So many repeated changes at the same positions that their history is erased — distances stop growing with divergence. The reason rpoB is excluded. |
| **Seed** | The starting number for the random choices a tree search makes. Same seed, same input file, same settings and one thread → the same tree. Change the thread count, or the order of the input sequences, and the search settles on a different, about equally good tree (see Reproducibility below). `amplicon.seed` is reused here. |
| **Site pattern** | A distinct column in the alignment. IQ-TREE's work per step scales with the number of these, which is why a ~440-column 16S alignment is "short". |
| **Split** | The two groups of tips an internal branch separates. Support values and RF distance are both about splits. |
| **Support (UFBoot, SH-aLRT)** | Numbers on internal branches saying how sure the data are about that split. UFBoot (ultrafast bootstrap): how often the branch reappears when the alignment columns are resampled, read as trustworthy at ≥ 95. SH-aLRT: a per-branch likelihood test, read at ≥ 80. FastTree's "SH-like" supports are on a 0–1 scale and are not bootstraps. |
| **Tip / internal branch / topology** | A tip (leaf) is an ASV. Internal branches connect groups of tips. The topology is the branching pattern alone, ignoring lengths. |
| **Tree length** | The sum of all branch lengths. The scale on which PD is measured. |
| **UniFrac** | A distance between two samples based on how much of the tree's branch length they share (unweighted: presence/absence; weighted: abundance-weighted). Root-dependent. |

## What MetaFlux does and does not do

MetaFlux aligns the ASVs, infers a tree, checks it, and writes it out. That is the whole
module. It computes **no diversity statistics** — no UniFrac, no Faith's PD, no
ordination, no PERMANOVA. Those belong in your R session, on your own filtered feature
set, and the [Using the tree in R](#using-the-tree-in-r) section below is the handover
point.

The value here is not that MetaFlux can build a tree you could not build yourself. It is
that the build becomes part of the versioned, seeded, resource-managed run, with every
resolved setting recorded, instead of a script you re-write for each project.

## 16S only — and why the other markers are not "coming soon"

A de novo tree is only meaningful for a marker that can be **globally aligned across
divergent taxa**. 16S rRNA meets that test. The others each fail it for a specific,
published reason, and those reasons are not the same as each other:

- **ITS** is non-coding and hypervariable, with roughly 20-fold length variation across
  fungi (full-length ITS runs from ~250 bp in some Saccharomycetales to ~1500 bp). It
  "is not amenable to robust multiple alignments and phylogenetic reconstruction much
  beyond the genus level" (Tedersoo et al. 2022), and QIIME 2's own fungal ITS tutorial
  recommends non-phylogenetic metrics such as Bray–Curtis or Jaccard instead. Use those.
  The old workaround, ghost-tree, is unmaintained.
- **gyrB** amplicons co-amplify the paralogue *parE* — up to 35% of reads overall, and
  60–95% within some families (Poirier et al. 2018). A paralogue is an **out-group**, so
  it enters a tree as a deeply divergent clade on a long stem: it inflates Faith's PD
  directly and distorts UniFrac while looking exactly like real phylogenetic structure.
  MetaFlux tags the paralogue so `discard: [other]` can remove it, but that is not
  applied by default.
- **rpoB** shows mutational saturation at **all three** codon positions at phylum and
  domain scope; the marker's reference study for microbial ecology switched to
  amino-acid alignment for exactly this reason (Case et al. 2007).
- **18S** is the closest candidate and the one we would revisit first, being an rRNA
  gene like 16S. But branch lengths from short 18S amplicons are documented to be
  underestimated (Medlar et al. 2014), V9 in particular carries little phylogenetic
  signal, and we found no study validating UniFrac or Faith's PD from short 18S
  amplicons. It is **deferred pending validation**, not excluded on principle.

Protein-coding markers (gyrB, rpoB, and COI when it arrives) would additionally need
**codon-aware alignment** and a codon or amino-acid substitution model — a plain
nucleotide alignment of a coding sequence inserts frameshift-breaking gaps that become
spurious branch length. That is a separate module, not a configuration option.

Enabling `phylogeny` for any marker other than 16S is forced off with a warning when the
config is read, before anything runs, and the run continues normally.

!!! note "This is a deliberate difference from other pipelines"
    QIIME 2 and nf-core/ampliseq are marker-agnostic in code. An ampliseq ITS run with a
    metadata file will build a MAFFT/FastTree phylogeny from bare ITS1 or ITS2 fragments
    and report weighted UniFrac, unweighted UniFrac and Faith's PD, by default, with no
    warning. We would rather not produce a number you cannot defend.

## Where the ASVs come from

The tree is built from the ASV set in `6.taxonomy/` — **after** the contaminant filter —
and never from `5.dada2/seqs_lenfilt.fasta`.

This is the single most important design decision in the module. MetaFlux's contaminant
filter (`amplicon.taxonomy.filter`: chloroplast, mitochondria, wrong domain, off-target)
runs *inside* the taxonomy step. The stage-60 FASTA therefore still contains exactly the
off-target ASVs that a de novo tree renders as long branches — and long branches are what
distort these metrics. Faith's PD is an absolute sum of branch length, so a long branch
adds to it directly; unweighted UniFrac is more sensitive still. In one published case, a
de novo tree gave three low-abundance archaeal ASVs an artificially long branch and
produced a **fake cluster separation** in unweighted UniFrac that vanished when the tree
was built properly.

Consuming the filtered set removes that problem outright, and it guarantees
the tree's tips and the abundance table you pair it with in R are the same feature set.

!!! note "The ASV set can shift slightly between runs, for a reason upstream of this module"
    With `taxonomy.method: sintax` on more than one thread, VSEARCH races several threads
    on one random-number stream, so per-rank confidence values shift between runs — on the
    16S test set by up to about 13 percentage points (an order call at 0.72 in one run,
    0.85 in the next), enough to add or drop up to three trailing ranks of a lineage near
    the cutoff. Between any two complete runs of that set from identical configs, 10–17 of
    211 ASVs changed their taxonomy string, with **every read count identical**, and the
    tree input stayed the same 211 ASVs each time. That is luck, not a guarantee: the
    `discard` tokens (`o__Chloroplast`, `f__Mitochondria`) sit exactly at the ranks whose
    confidence drifts, so a mitochondrial ASV called to family in one run and only to order
    in the next is discarded in the first run and reaches the tree in the second. Set
    `resources.threads.assign_taxonomy: 1` if you need the ASV set itself to be
    reproducible.

!!! warning "The filter has to be doing something"
    If `amplicon.taxonomy.filter` is disabled or its `keep`/`discard` lists are empty,
    the tree is built from an unfiltered ASV set and the protection above does not apply.
    MetaFlux warns, before anything runs, when you enable phylogeny in that state. Check the
    long-branch report either way.

## Choosing a backend

```yaml
    backend: iqtree     # iqtree | fasttree | raxml-ng
```

| Backend | What it is | When to pick it |
|---|---|---|
| **`iqtree`** *(default)* | Full maximum likelihood, IQ-TREE 3. Real model selection available. | The default for good reason — the best accuracy/effort trade-off at amplicon scale. |
| `fasttree` | Approximate ML (heuristic search, ML branch lengths). ~93 s for 10,000 ASVs on one core. | Very large ASV sets, or when you want the same engine QIIME 2 and ampliseq use. |
| `raxml-ng` | Full ML alternative. | Cross-checking a result against a second full-ML implementation. |

All three write the same canonical outputs, so switching backends does not change how you
use the result downstream. Their native files are kept in `7.phylogeny/<backend>/`, so
switching leaves the previous backend's output in place for comparison (the run log
points this out so you do not mistake it for current output).

### The model

```yaml
    model: auto         # auto | an explicit model string
```

`auto` resolves to **GTR+F+G4** for IQ-TREE, **GTR+Γ** for RAxML-NG, and **GTR+Γ** for
FastTree.

IQ-TREE's own default is `-m MFP`, which runs ModelFinder on every execution — 22 base
models crossed with frequency and rate-heterogeneity variants. MetaFlux pins a model
instead, because MFP makes runtime depend on the data and re-does on every run a decision
that does not change; IQ-TREE's documentation recommends exactly this select-once-then-pin
pattern for repeated analyses.

Model selection is still available, and whatever it picks is recorded in
`phylogeny.params.json`:

- `model: MFP` — IQ-TREE's ModelFinder
- `model: DNA` — RAxML-NG's MOOSE
- FastTree has **no model selection at all**. Its only nucleotide models are Jukes-Cantor
  and GTR, so only `gtr` (the default) and `jc` are accepted; anything else is a hard
  error rather than a silently ignored setting.

!!! info "Worth knowing about FastTree elsewhere"
    QIIME 2 invokes FastTree with no model flags at all, so the default 16S tree from
    QIIME 2 and ampliseq is a **Jukes-Cantor** tree, and there is no supported way to
    change it there. MetaFlux passes `-gtr` explicitly.

!!! info "What ModelFinder actually picks on the 16S test set — and whether it matters"
    Run once on the real 211-ASV test dataset ([PRJNA305879](../about/test-datasets.md),
    V5–V7), ModelFinder selects **TPM3u+R5** by BIC: a *simpler* substitution matrix
    than GTR (three exchange-rate classes instead of six) but a *richer* site-rate
    model (five free rate categories instead of a 4-category gamma). RAxML-NG's MOOSE
    agrees on the rate side, fitting GTR+R5. So on this dataset the two selectors agree on replacing the 4-category gamma with
    FreeRate and disagree on the matrix (TPM3u vs GTR) — the rate model, not the
    exchange-rate matrix, is where both depart from GTR+F+G4.

    Measured consequence of pinning GTR+F+G4 anyway, on that dataset: patristic
    distances between the two trees correlate at r = 0.998, per-sample PD at
    r = 0.9999 with identical sample ranking — the numbers your diversity analysis
    consumes are effectively unchanged. What does shift is the absolute scale (total
    tree length 18.6 under GTR+F+G4 vs 15.0 under TPM3u+R5, ~20%) and ~12% of
    internal branches — largely the poorly supported ones. Two practical rules
    follow: never compare absolute PD between runs built under different models, and
    if you want the selected model, run `model: MFP` once and then pin what it picked.

### Branch support

```yaml
    support: false
```

Off by default: neither UniFrac nor Faith's PD reads support values, and computing them
dominates the runtime. `true` adds `-B 1000 --alrt 1000` for IQ-TREE (read as UFBoot ≥ 95
and SH-aLRT ≥ 80 for a single-gene tree like this one), keeps FastTree's SH-like local
supports (0–1 scale — these are **not** bootstrap values, and the resampling that
produces them is seeded from `amplicon.seed`), and is not available for RAxML-NG, where
it is forced off with a warning.

!!! warning "Support values do not survive rerooting"
    In Newick a support value is stored as a label on the node below the edge it
    describes. When you reroot the tree in R, the edge above a node changes for every
    node between the old and new root, so labels can end up describing splits they were
    never computed for. Read support values off the unrooted tree as exported; treat
    them with suspicion on any rerooted copy you make yourself.

## The alignment

The aligner is **MAFFT**, fixed and not selectable. What you can choose is the strategy:

```yaml
    aligner_strategy: auto      # auto | linsi | fftns2
```

`auto` applies MetaFlux's own rule at MAFFT's documented ~200-sequence ceiling for its
accuracy strategies:

| ASVs | Strategy | Flags |
|---|---|---|
| < 200 | L-INS-i | `--localpair --maxiterate 1000 --threadit 0` |
| ≥ 200 | FFT-NS-2 | `--retree 2 --maxiterate 0` |

MAFFT's own `--auto` is deliberately **not** used. Its decision thresholds are not
published anywhere in MAFFT's documentation — they exist only in the source — and they
switch on dataset size. Adding a single ASV can cross a boundary and silently change the
alignment algorithm between two runs of the same study. MetaFlux's rule is documented,
uses MAFFT's own published ceiling, and the resolved choice is written into
`phylogeny.params.json`.

Set `fftns2` explicitly if you want two runs of different sizes to be strictly comparable.

`--adjustdirection` is never used: it prefixes reverse-complemented sequences with `_R_`,
which would break the ASV-ID match against your abundance table with no error raised.
Orientation is settled upstream, at primer trimming.

### No alignment masking

MetaFlux does not mask, trim, or filter alignment columns, and there is no option to.
Three reasons:

- Automated alignment filtering **frequently worsens** single-gene phylogenetic inference
  (Tan et al. 2015) — and a short 16S fragment is exactly that case.
- IQ-TREE explicitly instructs that constant sites must **not** be removed, because they
  inform branch-length estimation. Branch lengths are the quantity this module exports.
- The QIIME/ampliseq mask you might expect this to match is effectively a **pass-through**
  at its shipped defaults: `max_gap_frequency=1.0` removes nothing ever, and
  `min_conservation` is computed over non-gap characters only, so a column that is 95%
  gaps and 5% `A` scores 1.0 and is retained.

So "no masking" is closer to what those pipelines effectively do than it appears, and it
is honest about it.

This was also checked empirically on the 211-ASV test set, by masking the pipeline's own
alignment offline and rebuilding the tree (IQ-TREE, GTR+F+G4, same seed): trimAl
`-automated1` and `-gappyout` each removed 24% of columns, which cut **28% of total tree
length** (Faith's PD scales directly with that), moved the topology further from every
other tree in the comparison than any other single choice did, and left genus-level
coherence no better (10/19 genera monophyletic vs 11/19 unmasked — a difference of one genus). AliFilter was
gentler (14% of columns) but showed no gain either — and it is not installable from
bioconda, which rules it out as a MetaFlux option regardless. Nothing in the measurement
argues for adding a masking step.

## Rooting — the tree is exported unrooted, and only unrooted

`7.phylogeny/asv_16s.unrooted.nwk` is **the** deliverable, and there is deliberately no
option to produce a rooted copy.

Faith's PD and UniFrac are both root-dependent. In practice you will prune the tree
before computing them — decontam against your negative-control samples, abundance
filtering of rare ASVs — and any root placed on the **full** ASV set stops being valid
the moment the first tip is removed (a midpoint in particular moves). Midpoint rooting
is also unstable under rate variation and long branches (Mai et al. 2017), which is
precisely the situation a short-amplicon tree can present. A rooted file that is only
correct until you do the pruning you will certainly do would be a trap, not a
convenience.

So rooting happens in your R session, after pruning, immediately before the metric —
see [Using the tree in R](#using-the-tree-in-r).

## Tip labels

Tips are bare ASV IDs — `ASV_1`, `ASV_2`, … — with no taxonomy appended, ever. Beyond
tidiness there are two hard reasons: FastTree truncates a name at the first space unless
names are quoted, and both FastTree and RAxML-NG reject the `:,()` characters a lineage
string is full of.

Join taxonomy back on in R, keyed on the ASV ID, from `6.taxonomy/asv_table.txt`.

!!! note "ASV IDs are stable within a run, not across runs"
    IDs are assigned by decreasing abundance when the sequence table is built. `ASV_7` in
    two different runs is not the same organism. This holds even for two runs of the same
    reads with the same configuration: ASVs with exactly equal total abundance are numbered
    in the order they are first met in the reads, and that order is not fixed, because the
    phiX-removal step (bowtie2, multithreaded) writes its reads in a run-dependent order. In
    three complete test runs, two ASVs with 2 reads each swapped IDs once. Within a run the
    table, the taxonomy and the tree always agree with each other; only cross-run
    comparisons must be keyed on the sequence, never on the ID.

## Checks, and when no tree is built

`phylo_export` refuses to write the canonical tree unless it passes every check:

- every ASV in `asv_table.txt` appears **exactly once** as a tip, and there are no extra
  tips;
- the Newick parses;
- every branch length is present, finite, and non-negative.

**Fewer than 4 eligible ASVs**: no tree is built, and this is recorded rather than
treated as a failure. Below four tips there is only one possible unrooted topology, so
there is nothing to infer. The run completes normally and
`stats/phylogeny/phylogeny_qc.json` says why:

```json
{
  "skipped": true,
  "reason": "Only 3 eligible ASV(s) ... at least 4 are needed.",
  "n_eligible": 3
}
```

### The long-branch report

This is the number worth looking at. `phylogeny_qc.json` flags tips sitting on unusually
long branches, by two measures:

- **Pendant edge** — the length of the single branch leading to that tip. This is the
  **primary flag**. It needs no root, which makes it the right measure for an unrooted
  tree, and it is the direct expression of what a contaminant looks like: a sequence far
  from everything else, joined to the rest of the tree by one long branch.
- **Root-to-tip** on a midpoint-rooted copy, reported alongside. Better at spotting a
  whole off-target *clade* hanging off one long stem, where each member's own pendant
  edge is short.

The two use different threshold rules because they behave differently, and both rules
were tuned on data rather than reasoned out:

- Root-to-tip distances share a root and, on a tree with a clear centre, cluster, so **1.5× the median** works
  (roughly the ratio in the published artifact case). But this measure is nearly blind
  to a *single* long branch: on a synthetic test with one planted contaminant, midpoint
  rooting split the offending branch across the root and the contaminant sat at 1.05×
  the median — invisible — while its pendant edge was 22.7× the median. How many tips
  it flags also depends on the backend and on which optimum the search found: 9 to 58 of
  the same 211 ASVs across the trees checked for this release, against 6 to 7 by the
  pendant fence. Read the pendant flags first.
- Pendant edges cannot use a multiple of the median at all: a real ASV set carries many
  near-identical sequences whose pendant edges are effectively zero, so the median is
  close to zero and even generous multiples of it flag far too much of the tree (measured on
  the 16S test set: 81 of 211 tips at 5× the median, still 33 at 20×). The rule is instead the standard boxplot outlier fence,
  **Q3 + 3×IQR**, which flagged 7 tips there. Each is the only ASV of its lineage in the
  run — two are the run's sole Acidobacteriota and Bdellovibrionota — but two others are
  1–3-read ASVs classified no deeper than phylum and a third (5 reads) reaches only
  order, and nothing outside the run was checked, so "divergent singleton" and "artefact"
  cannot be told apart from the tree alone.

The five longest tips by each measure are always listed, whether or not they cross a
threshold.

Flagged tips are **reported, never removed**: dropping them would desynchronise the tree
from your abundance table. And a flag is not an accusation — in an environmental sample,
a long pendant edge is often a genuine divergent lineage with a single representative,
which is consistent with what the test set shows. Check what the tip was classified
as, consider whether the contaminant filter should have caught it, and re-run your
analysis without it to see whether your conclusion depends on it.

## Using the tree in R

MetaFlux stops at the tree; this is where you pick it up. Below is a complete worked
example — from the three files a run leaves behind through to a UniFrac ordination —
that was written and run against MetaFlux's own 16S test output, so the code is known
to work rather than merely to look plausible.

### Which packages, and why not phyloseq

The example uses **`ape`** (read and prune the tree), **`phangorn`** (midpoint rooting),
**`rbiom`** (UniFrac) and **`vegan`** (ordination, PERMANOVA), plus **`decontam`** if you
sequenced negative controls and **`rtk`** if you want rarefied non-phylogenetic indices.

It deliberately does **not** use `phyloseq`. That is not a judgement on anyone's past
work — it was the right tool for years — but three things have changed:

- **It is in maintenance mode.** Bioconductor's own successor is the
  `TreeSummarizedExperiment` / `mia` family, and the *Orchestrating Microbiome Analysis*
  book is explicitly positioned as the successor to the phyloseq tutorials.
- **Its object model buys you nothing here.** MetaFlux writes a plain table and a plain
  Newick file. Assembling them into an S4 object first, only to pull them apart again,
  adds ceremony between you and your data — and hides exactly the step where the tree and
  the table can silently stop matching.
- **Its UniFrac is slow.** It is implemented in R; `rbiom`'s is compiled, and the gap is
  orders of magnitude on a real ASV set.

If you want a maintained *framework* rather than plain matrices, use
`TreeSummarizedExperiment` + `mia`: it holds the tree in `rowTree()` and has the same
metrics. The example below stays framework-free on purpose, so that every object in it is
a matrix or a data frame you can print and inspect.

### The example

```r
# ─────────────────────────────────────────────────────────────────────────────
# From a MetaFlux 16S run to a UniFrac ordination.
#
# Reads the three files a MetaFlux amplicon run leaves behind, filters the ASVs
# the way a real study does, prunes and roots the tree, and ends at a UniFrac
# ordination with a PERMANOVA. Every step that can silently corrupt the result
# is checked rather than assumed.
#
# Packages (all CRAN except decontam, which is Bioconductor):
#   ape       read, prune and inspect the tree
#   phangorn  midpoint rooting
#   rbiom     UniFrac and Faith's PD (fast C++; reads a plain matrix + tree)
#   vegan     ordination and PERMANOVA
#   decontam  optional, only if you sequenced negative controls
# ─────────────────────────────────────────────────────────────────────────────

library(ape)
library(phangorn)
library(rbiom)
library(vegan)

# vegan and rbiom both export rarefy(), and phangorn and vegan both export
# diversity(); whichever package is attached last silently wins. Every rbiom call
# below is written as rbiom::fn() so this example cannot break on library() order.

out_dir <- "path/to/your/metaflux_output"          # <- edit this

# ── 1. Read what MetaFlux wrote ──────────────────────────────────────────────
# asv_table.txt: rows are ASVs, columns are your samples plus one `taxonomy`
# column. R's read.table handles MetaFlux's quoting as-is; check.names = FALSE
# keeps sample names such as "WT-1" from being mangled into "WT.1".
asv_table <- read.table(
  file.path(out_dir, "6.taxonomy", "asv_table.txt"),
  sep = "\t", header = TRUE, row.names = 1,
  quote = "\"", comment.char = "", check.names = FALSE
)

# Counts and taxonomy travel together in that file; split them apart, because
# every numeric step below needs a pure integer matrix.
taxonomy <- asv_table[, "taxonomy", drop = TRUE]
names(taxonomy) <- rownames(asv_table)
counts <- as.matrix(asv_table[, setdiff(colnames(asv_table), "taxonomy"), drop = FALSE])

# The tree. MetaFlux exports it unrooted on purpose: you are about to remove
# tips, and a root chosen before that stops being the right root afterwards.
tree <- read.tree(file.path(out_dir, "7.phylogeny", "asv_16s.unrooted.nwk"))

cat(sprintf("read %d ASVs x %d samples; tree has %d tips\n",
            nrow(counts), ncol(counts), length(tree$tip.label)))

# ── 2. Check the invariant before trusting anything ──────────────────────────
# Every phylogenetic metric silently returns a number even when the tree and the
# table describe different features — it is just the wrong number. MetaFlux
# guarantees these two sets match; verify it anyway, because you are about to
# subset both and that is where they drift apart.
stopifnot(setequal(rownames(counts), tree$tip.label))
stopifnot(!any(duplicated(tree$tip.label)))

# ── 3. Remove contaminants (only if you sequenced negative controls) ─────────
# Skip this block if you did not. decontam's prevalence method compares how
# often each ASV appears in true samples versus in blanks; threshold 0.5 is the
# "more prevalent in controls than in samples" rule, which is the stricter and
# more usual choice for a well-designed blank set.
#
# is_control <- colnames(counts) %in% c("BLANK1", "BLANK2")
# library(decontam)
# verdict <- isContaminant(t(counts), neg = is_control, method = "prevalence",
#                          threshold = 0.5)
# cat(sprintf("decontam flagged %d of %d ASVs\n", sum(verdict$contaminant), nrow(verdict)))
# counts <- counts[!verdict$contaminant, !is_control, drop = FALSE]

# ── 4. Filter rare ASVs ──────────────────────────────────────────────────────
# A plain, defensible rule: keep an ASV seen at least `min_reads` times in at
# least `min_samples` samples. Tune to your design — this is a study decision,
# not a technical default. Filtering before the metrics matters because a
# singleton on a long branch moves Faith's PD and unweighted UniFrac more than
# any abundant ASV does.
min_reads   <- 2
min_samples <- 2
keep <- rowSums(counts >= min_reads) >= min_samples
cat(sprintf("abundance filter keeps %d of %d ASVs (%.1f%% of reads)\n",
            sum(keep), length(keep), 100 * sum(counts[keep, ]) / sum(counts)))
counts <- counts[keep, , drop = FALSE]

# ── 5. Prune the tree to the ASVs that survived ──────────────────────────────
# Prune, do not rebuild. Dropping a tip removes it and merges the two branches
# that met at its parent, so distances among the remaining tips are unchanged.
# The pruned tree is exactly the full tree restricted to the ASVs you kept, and
# its branch lengths were estimated from more data than a tree re-inferred from
# the survivors alone would have had.
tree <- keep.tip(tree, rownames(counts))
stopifnot(setequal(rownames(counts), tree$tip.label))
cat(sprintf("tree pruned to %d tips\n", length(tree$tip.label)))

# ── 6. Root the tree — now, not earlier ─────────────────────────────────────
# Unweighted UniFrac and Faith's PD both depend on where the root sits. On the
# MetaFlux 16S test data, moving the root changed unweighted UniFrac distances
# by up to 7e-3; weighted UniFrac was unaffected. Midpoint rooting is the usual
# choice when no outgroup is available. Do it here, after pruning, because the
# midpoint of the full tree is not the midpoint of this one.
#
# With a real outgroup, prefer it:  tree <- root(tree, outgroup = "ASV_57", resolve.root = TRUE)
tree <- midpoint(tree)
stopifnot(is.rooted(tree))

# ── 7. Even out sequencing depth: rarefaction, done properly ────────────────
# Sequencing depth is not biology. A sample sequenced ten times deeper shows more
# lineages, so UniFrac (especially unweighted) and Faith's PD both track depth
# unless you account for it. There are three defensible positions, and it is
# worth knowing which one this code takes:
#
#   1. Do not rarefy; normalise instead (DESeq2/edgeR-style, or CSS). Argued by
#      McMurdie & Holmes (2014). Their target was RARE-FYING — one subsample,
#      data thrown away.
#   2. Rarefy once. Simple, and what most tutorials show. It is also the version
#      (1) rightly criticises: one draw of the dice decides your distances.
#   3. RAREFACTION proper — subsample many times, compute the statistic on each
#      draw, average the results. Schloss (2024) re-examined the 2014 critique
#      and found it does not carry over to this version, which controls uneven
#      effort better than the normalisation methods do.
#
# This example takes position 3. Nothing below throws data away permanently: the
# subsampling happens inside the loop, once per iteration.
n_iters <- 100
depth   <- min(colSums(counts))
cat(sprintf("sample depths: %s\nrarefying to %d, %d iterations\n",
            paste(sprintf("%s=%d", colnames(counts), colSums(counts)), collapse = ", "),
            depth, n_iters))

# One subsample of every sample to `depth`, without replacement.
rarefy_once <- function(m, depth) {
  out <- apply(m, 2, function(x) {
    drawn <- sample(rep(seq_along(x), x), depth)
    as.integer(table(factor(drawn, levels = seq_along(x))))
  })
  rownames(out) <- rownames(m)
  out
}

# ── 8. Alpha diversity: Faith's PD, exactly, without iterating ──────────────
# For phylogenetic diversity you do not need to simulate at all. A branch is
# missing from a subsample only if NONE of the reads below it were drawn, and
# for m reads drawn without replacement from N that probability is
# hypergeometric. So the expected PD after rarefying has a closed form
#
#     E[PD] = sum over branches of  L_b * (1 - C(N - k_b, m) / C(N, m))
#
# with k_b the reads below branch b — Hurlbert's rarefaction argument applied to
# branches instead of species (Nipperess & Matsen 2013). This returns the exact
# value that averaging infinitely many rarefaction iterations would converge to:
# no seed, no iteration count, no Monte Carlo noise. Verified against a
# 400-iteration simulation on MetaFlux's test data, and it uses the same
# root-inclusive convention as rbiom's `faith`.
expected_faith_pd <- function(tree, x, depth) {
  total_reads <- sum(x)
  if (total_reads < depth) return(NA_real_)          # too shallow to rarefy to `depth`
  descendants  <- phangorn::Descendants(tree, seq_len(max(tree$edge)), type = "tips")
  child_node   <- tree$edge[, 2]
  reads_below  <- vapply(child_node,
                         function(nd) sum(x[tree$tip.label[descendants[[nd]]]]),
                         numeric(1))
  # lchoose keeps this stable for large N; when N - k < m the term is 0 and the
  # branch is certain to be present, which exp(-Inf) gives for free.
  p_present <- 1 - exp(lchoose(total_reads - reads_below, depth) -
                       lchoose(total_reads, depth))
  sum(tree$edge.length * p_present)
}

faith_pd <- vapply(colnames(counts),
                   function(s) expected_faith_pd(tree, counts[, s], depth),
                   numeric(1))
print(round(faith_pd, 4))

# Non-phylogenetic indices have no such shortcut, so those do need iterating.
# The rtk package does this efficiently — but note it takes no tree and offers
# no phylogenetic measure, which is why Faith's PD is handled above instead.
#
# library(rtk)
# r <- rtk(counts, repeats = n_iters, depth = depth, margin = 2)
# shannon_mean <- sapply(get.diversity(r, div = "shannon"), mean)

# ── 9. Beta diversity: distances averaged over rarefaction iterations ───────
# Compute the distance matrix on each subsample and average the matrices — not
# the other way round. Averaging the rarefied tables first would just hand back
# something close to the original table and defeat the point. Averaging
# distances is what mothur's `dist.shared(..., subsample=T, iters=)` does for
# Bray-Curtis, and nothing about UniFrac forbids the same treatment: the mean of
# distance matrices is still symmetric, still zero on the diagonal, and still
# satisfies the triangle inequality.
#
# Measured on MetaFlux's 16S test data: two single rarefactions produced
# unweighted UniFrac distances differing by up to 0.109, while two independent
# 100-iteration averages differed by at most 0.016. That factor of ~7 is the
# whole argument for doing it this way.
#
# All the metrics are computed on the SAME subsample within each iteration, so
# that any later comparison between them reflects the metrics themselves and not
# two different draws of the dice.
average_dists <- function(counts, tree, metrics, depth, iters) {
  totals <- setNames(vector("list", length(metrics)), metrics)
  for (i in seq_len(iters)) {
    rarefied <- rarefy_once(counts, depth)
    biom_i   <- rbiom::as_rbiom(rarefied, tree = tree)
    for (m in metrics) {
      d <- as.matrix(rbiom::bdiv_distmat(biom_i, bdiv = m))
      totals[[m]] <- if (is.null(totals[[m]])) d else totals[[m]] + d
    }
  }
  lapply(totals, function(x) as.dist(x / iters))
}

set.seed(42)   # the subsampling is random; fix it so the run is repeatable
dists <- average_dists(counts, tree,
                       c("unweighted_unifrac", "weighted_unifrac", "bray"),
                       depth, n_iters)
uw <- dists$unweighted_unifrac
wt <- dists$weighted_unifrac
bc <- dists$bray          # Bray-Curtis: no tree involved, the non-phylogenetic view

# ── 10. Ordination ──────────────────────────────────────────────────────────
# PCoA (classical multidimensional scaling) is the standard partner for a
# distance matrix. The eigenvalues tell you how much of the structure the first
# two axes actually show — quote it on the axis labels, never omit it.
pcoa <- cmdscale(uw, k = 2, eig = TRUE)
var_explained <- 100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0])
plot(pcoa$points, pch = 19, cex = 1.4,
     xlab = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
     ylab = sprintf("PCoA 2 (%.1f%%)", var_explained[2]),
     main = "Unweighted UniFrac, PCoA")
text(pcoa$points, labels = rownames(pcoa$points), pos = 3, cex = 0.7)

# ── 11. Test a hypothesis ───────────────────────────────────────────────────
# PERMANOVA asks whether group centroids differ. It is sensitive to differences
# in within-group spread as well, so always run betadisper alongside: a
# significant adonis2 with a significant betadisper may only mean one group is
# more variable than the other.
#
# group <- factor(c("treated", "treated", "control", "control", "treated", "control"))
# print(adonis2(uw ~ group, permutations = 999))
# print(permutest(betadisper(uw, group), permutations = 999))

# ── 12. Is the tree telling you anything Bray-Curtis doesn't? ───────────────
# The question worth asking before you report a phylogenetic result: does it
# actually differ from the answer you would have got without a tree? Two ways to
# ask, and they answer slightly different questions.
#
# Mantel compares the DISTANCE MATRICES directly — do the two metrics rank the
# sample pairs the same way? Procrustes compares the ORDINATIONS — after
# rotating, scaling and translating one onto the other, how well do the points
# superimpose? Procrustes is the more relevant test when what you publish is a
# PCoA plot, because that is exactly the object it compares.
cat("\n--- Mantel: do the distance matrices agree? ---\n")
print(mantel(uw, bc, permutations = 999))   # unweighted UniFrac vs Bray-Curtis
print(mantel(wt, bc, permutations = 999))   # weighted   UniFrac vs Bray-Curtis

cat("\n--- Procrustes: do the ordinations agree? ---\n")
# Rotate the Bray-Curtis ordination onto each UniFrac ordination. `protest` adds
# a permutation test; its "Correlation in a symmetric Procrustes rotation" is the
# number to quote (the `t0` field), and 1 means the two are the same picture.
pr_uw <- protest(cmdscale(uw, k = 2), cmdscale(bc, k = 2), permutations = 999)
pr_wt <- protest(cmdscale(wt, k = 2), cmdscale(bc, k = 2), permutations = 999)
print(pr_uw)
print(pr_wt)

# Per-sample residuals show WHERE the two views disagree: a large residual means
# that sample sits in a different place depending on whether the tree was used.
cat("\nper-sample Procrustes residuals (unweighted UniFrac vs Bray-Curtis):\n")
print(round(residuals(pr_uw), 3))

# The two ordinations superimposed. Arrows run from the Bray-Curtis position to
# the UniFrac one, so a long arrow is a sample the tree moved.
plot(pr_uw, main = "Bray-Curtis -> unweighted UniFrac")
```

### What the checks in it are for

**The invariant check (step 2) is the one not to skip.** Every phylogenetic metric will
happily return a number when the tree and the table describe different feature sets — it
is simply the wrong number, and nothing warns you. MetaFlux guarantees the two match on
output; the check exists because *you* are about to subset both, and that is where they
drift apart.

**Prune, don't rebuild.** Dropping a tip removes it and merges the two branches that met
at its parent, so patristic distances among the remaining tips are unchanged — verified
on this module's own output. UniFrac or Faith's PD on the pruned tree therefore equals
what you would get from the full tree restricted to those features. A tree re-inferred
from only the survivors would differ, and having been estimated from fewer sequences,
would generally have *worse*-informed branch lengths. Rebuild only if you removed a very
large fraction of the ASVs, or if the removed set was specifically what distorted the
topology; there is no benchmarked threshold, so it is a judgement call.

**Root after pruning, not before.** This is the step people get wrong, and it is why
MetaFlux exports the tree unrooted. Measured on the 16S test data: moving the root
changed unweighted UniFrac distances by up to 7×10⁻³, while weighted UniFrac was
unaffected to numerical precision. So the root position matters for exactly the metric
that is already the most artefact-prone.

**Watch the namespace collisions.** `vegan` and `rbiom` both export `rarefy()`, and
`phangorn` and `vegan` both export `diversity()`. Whichever package is attached last
wins, silently. The example calls every `rbiom` function as `rbiom::fn()` for that reason.

**On rarefaction — the choice this example makes, and why.** Sequencing depth is not
biology, and both unweighted UniFrac and Faith's PD track it. Three positions are
defensible, and it matters which one you are taking:

| Position | What it does | Where it stands |
|---|---|---|
| Do not rarefy; normalise | DESeq2/edgeR-style or CSS scaling | Argued by McMurdie & Holmes (2014) — whose target was *rarefying*, not rarefaction |
| Rarefy once | One subsample, the rest discarded | What most tutorials show, and exactly what the 2014 critique is about |
| **Rarefaction proper** | Subsample many times, compute the statistic each time, average | What this example does. Schloss (2024) re-examined the 2014 critique and found it does not carry to this version, which controls uneven effort better than the normalisation methods |

The distinction the 2014 paper's title obscured is between *rarefying* (one draw, data
thrown away) and *rarefaction* (many draws, averaged, nothing permanently discarded).
The example subsamples inside the loop, so no data is lost.

**For beta diversity, average the distance matrices, not the tables.** Compute UniFrac
on each subsample and average the resulting matrices — the order matters, because
averaging the rarefied tables first would hand back something close to the original
table. This is what mothur's `dist.shared(..., subsample=T, iters=)` does for
Bray–Curtis, and nothing about UniFrac forbids the same treatment: an average of
distance matrices is still symmetric, still zero on the diagonal, and still satisfies
the triangle inequality. Measured on the 16S test data: two single rarefactions gave
unweighted UniFrac distances differing by up to **0.109**, while two independent
100-iteration averages differed by at most **0.016** — a sevenfold reduction, which is
the whole argument for doing it this way. Weighted UniFrac is far less affected
(0.0098 → 0.0008), as expected for an abundance-weighted metric.

**For Faith's PD, do not iterate at all — there is a closed form.** A branch is missing
from a subsample only if none of the reads below it were drawn, and for *m* reads drawn
without replacement from *N* that probability is hypergeometric. So

```
E[PD] = sum over branches b of  L_b * ( 1 - C(N - k_b, m) / C(N, m) )

  L_b = length of branch b        k_b = reads below branch b
  N   = reads in the sample       m   = rarefaction depth
```

— Hurlbert's rarefaction argument
applied to branches instead of species (Nipperess & Matsen 2013). This is the exact
value that averaging infinitely many iterations converges to: no seed, no iteration
count, no Monte Carlo noise. The example implements it in eight lines. It was checked
against a 400-iteration simulation on the test data (agreement to within simulation
error) and uses the same root-inclusive convention as `rbiom`'s `faith`.

!!! note "rtk does not do Faith's PD"
    `rtk` is the efficient choice for *multiple rarefaction of non-phylogenetic indices*
    — it returns richness, Shannon, Simpson, inverse Simpson, Chao1 and evenness, in
    compiled code. It takes no tree argument and offers no phylogenetic measure, so it
    cannot produce Faith's PD. Use the closed form above for PD, and `rtk` for the rest.

**Report both UniFracs.** Unweighted asks which lineages are present, weighted asks how
abundant they are. Unweighted is the sensitive one — it is where a single long branch or
a contaminant does its damage, and it is the metric the long-branch report exists to
protect. If the two disagree, the difference lives in the rare taxa, which is where a de
novo tree from a short marker is least trustworthy.

**Check the tree against Bray–Curtis before you believe it.** The last step of the
example asks the question that decides whether enabling this module changed anything:
does the phylogenetic result differ from the one you would have got without a tree?
Two tests, answering slightly different questions — Mantel compares the distance
matrices (do the metrics rank sample pairs the same way?), Procrustes compares the
ordinations after rotating one onto the other (do the *plots* say the same thing?).
Procrustes is the more relevant of the two when a PCoA is what you publish.

On the 16S test data, with all three metrics computed on the same rarefaction draws:

| Comparison | Mantel *r* | Procrustes correlation | *p* |
|---|--:|--:|--:|
| Weighted UniFrac vs Bray–Curtis | 0.911 | 0.862 | 0.003 |
| Unweighted UniFrac vs Bray–Curtis | 0.031 | 0.238 | 0.98 |

That contrast is the useful part, and it generalises further than this small dataset
does. **Weighted UniFrac reproduces Bray–Curtis almost exactly**, because both are
dominated by the abundant taxa — and among abundant, well-sampled lineages the tree
adds little that abundance did not already say. **Unweighted UniFrac gives a completely
different picture**, because it is driven by which rare lineages are present, and it is
their placement that the tree supplies.

So the practical rule: if your conclusion rests on **weighted** UniFrac, the tree is
largely confirming what Bray–Curtis already told you, and little hangs on the
phylogeny being right. If it rests on **unweighted** UniFrac, it rests entirely on the
tree — and specifically on the placement of rare, long-branch tips, which is what a de
novo tree from a short marker is least able to get right. That is the case where you
should read `stats/phylogeny/phylogeny_qc.json`'s long-branch report before believing
the result, and re-run the analysis without the flagged tips to see whether the
conclusion survives.

The per-sample Procrustes residuals tell you *where* the two views disagree, and
`plot(pr_uw)` draws the two ordinations superimposed with an arrow per sample: a long
arrow is a sample the tree moved.

!!! warning "Six samples is an illustration, not a result"
    The test dataset has six samples, so every *p*-value above is bounded below by
    1/720 ≈ 0.0014 and none of it is evidence about anything biological. The numbers
    are there to show the machinery working and the contrast between the two UniFracs;
    run the same comparison on your own study before drawing conclusions from it.

**PERMANOVA needs replication.** `adonis2` on a handful of samples cannot produce a
meaningful *p*-value: with two samples per group there are only a few hundred distinct
permutations. Always pair it with `betadisper` — a significant `adonis2` alongside a
significant `betadisper` may mean only that one group is more variable than the other,
not that the centroids differ.

## What is recorded

`7.phylogeny/phylogeny.params.json` holds everything resolved at run time: the resolved
aligner strategy, the resolved (or selected) model, effective thread counts, tool
versions, the seed, SHA-256 checksums of both input tables and the extracted FASTA, and
any `extra_args`. That record is what lets you pin an automatic choice explicitly later.

## Reproducibility, and what is not yet verified

`amplicon.seed` seeds every backend — there is no separate phylogeny seed.

**Exactly reproducible only for the same input file, the same seed and one thread.**
Everything else about the tree search is a heuristic that settles into one of several
nearly equally good trees, and which one it finds depends on things that look
irrelevant. All of the following was measured on the 211-ASV test set with IQ-TREE 3.1.3
(GTR+F+G4, seed 42):

- **Same alignment file, 1 thread:** two runs gave byte-identical Newick files.
- **Same alignment file, 4 threads (the shipped default): not reproducible.** On one
  input file, three runs reached the same log-likelihood and tree length and one of
  them differed from the other two on 12 of 416 splits, every one a branch of
  2 × 10⁻⁶ or less that `ape::di2multi(tree, tol = 1e-5)` removes. On another input
  file (the same sequences, two of them in a different order) two runs ended in
  different optima: log-likelihoods 0.1 apart, tree lengths 18.0 and 18.7, and RF 0.08
  and 0.12 against the single-threaded tree of the same file (patristic r 0.996 and
  0.984). A multithreaded search visits candidate trees in an order that depends on
  thread timing, and can settle in a different local optimum; `di2multi()` does not
  remove that kind of difference.
- **Same sequences, two of them in a different order:** a different tree. Three complete
  runs of the pipeline from raw reads, identical configuration, gave identical ASV
  sequences and counts, but in one run two ASVs with the same total abundance (2 reads
  each) received each other's IDs, because the read order after phiX removal is not
  fixed between runs (see the note on ASV IDs above). That swapped two records in the
  alignment input, and IQ-TREE converged to a different local optimum: 154 of 416 splits
  differ (normalized Robinson–Foulds 0.37), patristic correlation r = 0.93, total tree
  length 18.0 against 18.6 — and a *better* log-likelihood (−15 485 against −15 509). At
  1 thread the same reordering gave the same picture (RF 0.35, r = 0.94). Per-sample
  Faith's PD hardly noticed: r = 0.999 between the two trees, sample ranking identical.
- **Thread count alone:** the reordered input at 1 and at 4 threads differed on 34 splits
  (RF 0.08, r = 0.996), so the thread count also moves the search, but less than the
  input order does.

What follows from this: read the topology as *one* of several trees the data support
about equally well, not as *the* tree; compare trees by patristic correlation or by
per-sample PD, never with `diff`, and use `ape::dist.topo()` only after `di2multi()`;
and to reproduce a tree exactly between runs, the whole path to it must be identical —
the same eligible ASV set and the same ASV numbering (neither is guaranteed between runs
with the shipped defaults; see the two notes above on `sintax` threads and on ASV IDs),
`resources.threads.phylo_tree: 1` (140 s against 90 s at 4 threads on 211 ASVs, hours
on thousands) and the same seed. The distances and diversity values you actually use
are the stable part.
- Single-threaded FastTree with `support: false` uses no randomness at all and is fully
  deterministic. Its parallel build, `FastTreeMP`, is documented non-deterministic and is
  **never** used by MetaFlux.
- RAxML-NG's default seed is the wall clock, so MetaFlux always passes `--seed`
  explicitly.

Two questions that used to sit here as "not yet verified" now have measured answers,
from the benchmark on the 211-ASV test set described on the
[Phylogeny validation](../about/validation.md) page:

- **The L-INS-i / FFT-NS-2 tier boundary does change the tree.** Same data, same
  backend and model, different tier: 46% of internal branches differ. Part of that is the tree search
  itself — on the reference fragments, where seeds were replicated, the *same*
  alignment rebuilt with another seed already differs by 22–28%, and the tier switch
  by 34–39% at a matched seed — so the tier is a real effect, but the 46% figure (one
  tree per setting) does not separate the two. Patristic
  distances stay far more stable (r = 0.95) and per-sample PD barely moves
  (r = 0.997), but topology is not comparable across the boundary. If you will
  compare runs whose ASV counts straddle 200, pin `aligner_strategy: fftns2`.
- **How much signal a ~440-column alignment of ~376 bp fragments carries**: likelihood mapping (`iqtree3 --lmap`, 5,000
  quartets, GTR+F+G4) fully resolves **77.9%** of quartets, leaving ~22% partly or
  fully unresolved; and in a support run (`-B 1000 --alrt 1000`) only **32%** of
  internal branches meet the documented UFBoot ≥ 95 + SH-aLRT ≥ 80 read. Both numbers
  say the same thing as the framing paragraph: fine structure and distances are
  usable, deep topology is not to be trusted.

## How these defaults were checked

Two experiments, written up in full on the
[Phylogeny validation](../about/validation.md) page:

- **A sensitivity analysis on the real 211-ASV test run** — one setting changed at a
  time (alignment strategy, masking, backend, model) and the output compared to the
  pipeline default. Result: no choice moves the distances your diversity metrics use by
  more than a few percent, and masking removes 28% of tree length for no measurable gain.
- **A full-gene check on 250 SILVA reference genes** — the V5–V7 amplicon cut out of each
  full-length gene *in silico*, the fragment tree compared to the full-gene tree. Result:
  every fragment tree, under every setting, differs from the full-gene tree on ~58% of
  internal branches (Robinson–Foulds distance) while patristic distances correlate at 0.83–0.88 — and the spread
  between settings (0.04) is smaller than the search noise of a single setting.

Read those numbers as *"this choice moved the output by this much on this data"*, not as
accuracy claims; the validation page says exactly what each experiment can and cannot
support.

## Advanced: `extra_args`

```yaml
    extra_args:
      mafft: ""
      iqtree: ""
      fasttree: ""
      raxml_ng: ""
```

Appended verbatim to the tool's command line and recorded in `phylogeny.params.json`.
**Not validated.** A flag that contradicts the settings above — `-T AUTO`, say — is
passed straight through and is your responsibility. This is the one route by which the
machine-dependent behaviour MetaFlux otherwise blocks can re-enter a run. Normally leave
these empty.

## Resources

Threads and memory come from the shared `resources` block, like every other rule:
`phylo_input`, `phylo_align`, `phylo_tree`, `phylo_export`, `phylo_qc`. See
[Configuration](../reference/configuration.md#threads). The defaults are deliberately
modest, and `phylo_tree` behaves differently per backend — IQ-TREE parallelises across
alignment columns and a 16S alignment is short; RAxML-NG clamps to its own recommendation
because it hard-errors when given too many threads; FastTree is pinned single-threaded.

## References

Janssen et al. 2018, *mSystems* 3:e00021-18 (full-gene design) · Medlar et al. 2014, *BMC Evol Biol* 14:235 (fragment branch-length underestimation) · Apprill et al. 2015, *Aquat Microb Ecol* 75:129 and Parada et al. 2016 (V4 primers 806R-B / 515F-Y) · Tan et al. 2015, *Syst. Biol.* 64:778 · Mai et al. 2017, *PLoS ONE* 12:e0182238 ·
Tedersoo et al. 2022, *Mol Ecol* 31:2769 · Poirier et al. 2018, *PLoS ONE* 13:e0204629 ·
Case et al. 2007, *Appl Environ Microbiol* 73:278. Tool citations and the benchmark's methods references (Robinson & Foulds, Sokal &
Rohlf, Faith, Strimmer & von Haeseler, Hoang, Guindon, Schwarz, Tavaré, Kimura, Yang,
Soubrier, Kalyaanamoorthy) are on the [Citation](../about/citation.md) page.
