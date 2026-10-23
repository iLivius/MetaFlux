# Amplicon mode

Amplicon mode takes raw paired-end Illumina metabarcoding reads and returns a table of
exact amplicon sequence variants (ASVs) with a taxonomic assignment for each one. The
denoising engine is DADA2; six markers are supported — 16S rRNA, ITS, 18S rRNA, gyrB,
rpoB and COI — and the marker is what decides which reference databases are used, whether a
region extractor runs, and which classifiers are available.

Two config keys define the run:

```yaml
mode: amplicon
amplicon:
  type: 16S          # 16S | ITS | 18S | gyrB | rpoB | COI
```

Everything else is either a path (input reads, output directory, primer FASTAs) or
tuning. Marker-specific facts — database URLs, rank models, which extractor applies —
live in the marker packs under `workflow/markers/<type>.yaml` and are never set by hand.

## The steps at a glance

The pipeline runs as nine stages, plus an optional tenth. Output directories are numbered
in pipeline order under `out_dir`, so the directory listing reads like the workflow.

| # | Step | Tool | Snakemake rule | Writes to |
|---|------|------|----------------|-----------|
| 1 | PhiX removal *(optional)* | bowtie2 | `rm_phix` | `02.no_phix/` |
| 2 | Primer trimming | Cutadapt | `trim_primers` | `03.stripped/` |
| 3 | Per-stage QC | Falco | `falco` | `stats/falco/` |
| 4 | Amplicon-length probe | Cutadapt (or direct read) | `amplicon_probe` | `refdb/cache/` |
| 5 | Truncation-length selection | in-house | `pick_trunclen` | `stats/trunclen.json` |
| 6 | ASV inference | DADA2 | `dada_filter`, `dada_seqtab` | `04.filtered/`, `05.dada2/` |
| 7 | Marker-region extraction *(optional)* | Metaxa2 or ITSx | `target_extract` | `05.dada2/` |
| 8 | ASV length filter | in-house | `dada_length_filter` | `05.dada2/`, `stats/dada2/` |
| 9 | Taxonomy assignment | DADA2 RDP or VSEARCH SINTAX | `assign_taxonomy` | `06.taxonomy/` |
| 10 | [Phylogeny](phylogeny.md) *(optional, 16S only)* | MAFFT + IQ-TREE, FastTree or RAxML-NG | `phylo_input`, `phylo_align`, `phylo_tree`, `phylo_export`, `phylo_qc` | `07.phylogeny/`, `stats/phylogeny/` |

Two more rules close the run: `aggregate_read_counts` builds `stats/read_tracking.txt`,
and `multiqc` builds `multiqc/multiqc_report.html`. One small rule runs before any of the
taxonomy: `check_taxonomy_reference` checks the reference the classifier will read
against the marker pack and writes `stats/taxonomy_reference_check.json` (see
[step 9](#9-taxonomy-assignment)).

Stage 10 is off by default and runs only for 16S. It takes the contaminant-filtered ASVs
from `06.taxonomy/`, aligns them, and writes an unrooted tree for phylogenetic-diversity
analysis you perform yourself — MetaFlux computes no diversity statistics. See
[Phylogeny](phylogeny.md).

## Input

Amplicon mode requires local paired FASTQs named `{sample}_R1.fastq.gz` and
`{sample}_R2.fastq.gz` in `input.fastq_dir`. Sample names are discovered from those
filenames at start-up, when Snakemake reads the config and works out which jobs it will
need. SRA-style `{sample}_1/_2` naming is
accepted by shotgun mode but rejected here with an explicit error, because the two
conventions cannot be told apart safely once mixed.

The raw files are never modified. `link_reads` symlinks them into `01.reads/` purely so the
inputs of a run are visible next to its outputs; every rule that needs raw reads reads
them from `fastq_dir` directly.

## 1. PhiX removal *(optional)*

Illumina runs are routinely spiked with the PhiX174 genome as a sequencing control, and
those reads have no place in a community profile. `rm_phix` maps each pair against a
bowtie2 index of the PhiX reference and keeps the pairs that fail to align *concordantly*
(`--un-conc-gz`), which are written to `02.no_phix/`. Failing to align concordantly is not
the same as not matching at all: a pair in which only one mate aligned to PhiX is also
written out and kept, so a small PhiX residue can survive. The shotgun path is stricter —
BBDuk matches by k-mer and, with its default `removeifeitherbad=t`, drops a pair as soon
as either mate matches.

The reference is fetched once from NCBI (`references.phix.fetch_url`, RefSeq assembly
GCF_000819615.1) by `fetch_phix` and indexed by `build_phix_index`; both are cached under
`refdb/phix/`, so this cost is paid on the first run only.

Toggled with `amplicon.decontamination.remove_phix` (shipped default `true`). When it is
`false` the step is skipped entirely: `02.no_phix/` is not created, Cutadapt reads the raw
FASTQs instead, and the `nophix` QC stage, its read count and its read-tracking column all
disappear rather than showing up empty.

The bowtie2 alignment summary is written to `logs/rm_phix/{sample}_phixFilter.log` and is
picked up by MultiQC, so the PhiX fraction per sample is visible in the report.

## 2. Primer trimming

DADA2 requires primer-free reads: universal primers carry degenerate positions, and DADA2
reads that synthetic variation as real biological variation. This inflates the ASV count
and, because the spurious variants share long identical stretches, it interferes with
chimera detection — unremoved primers are the usual cause of an implausibly high chimeric
fraction. `trim_primers` runs Cutadapt on the output of step 1 (or on the raw reads when
PhiX removal is off) and writes to `03.stripped/`.

Trimming happens in two calls, because a read can carry primer sequence at both ends:

| Call | Cutadapt options | Purpose |
|------|------------------|---------|
| 5′ | `-g file:fwd -G file:rev`, `--discard-untrimmed` (and `--overlap 10` with `orientation: mixed`, see [below](#the-mixed-orientation-second-pass)) | Removes the forward primer from R1 and the reverse primer from R2. Pairs without a recognisable primer are dropped. |
| 3′ | `-a file:rev_rc -A file:fwd_rc` | Removes read-through: when the amplicon is shorter than the read, the sequencer runs past the far end into the reverse complement of the other primer. Reads without read-through pass through untouched. |

The two reverse-complemented primer files are produced once per run by `revcomp_primers`
(seqtk) into `_aux/primers/`.

Both calls use `amplicon.cutadapt.max_error_rate` (default `0.2`, i.e. up to 20% of the
matched primer length may differ), `amplicon.cutadapt.min_length` (default `50`, dropping
anything shorter after trimming) and `--pair-filter=any`, which discards both mates when
either one fails a filter.

Degeneracy costs nothing from that error budget: Cutadapt reads IUPAC codes in the primer
natively, so an `R` in the primer matches A or G for free. The budget is there for
sequencing errors in the first cycles and for genuine primer/template mismatch, which is
common with universal primers. Cutadapt's own default is `0.1`; `0.2` is deliberately more
permissive here because the 5′ pass runs with `--discard-untrimmed`, so a primer that
fails to match costs the whole pair. The per-sample match rates in `stats/cutadapt/` — or
the Cutadapt section of the MultiQC report — show what the tolerance is actually buying.

There is deliberately **no quality trimming here**: DADA2 learns its error model from the
run's own quality scores, and pre-trimming on quality would distort it by removing exactly
the low-quality bases the model needs in order to fit.

### The mixed-orientation second pass

Some library preparations (and most amplicon data pooled from several runs) deliver reads
in both orientations: in some pairs R1 carries the forward primer, in others it carries
the reverse one. Because the 5′ call uses `--discard-untrimmed`, a single fixed-orientation
pass would silently throw away every swapped pair — often half the data.

Setting `amplicon.primers.orientation: mixed` adds a second Cutadapt pass over the same
input with the primers swapped (`-g rev -G fwd`, then `-a fwd_rc -A rev_rc`). A pair this
pass keeps still has its mates the other way round, its R1 being the reverse read, so the
two mates are swapped back: the pair's R2 goes to the R1 file and its R1 to the R2 file,
after the pairs of the first pass. A swapped pair then reaches DADA2 as the pairs of the
first pass do, with its forward read in R1, starting right after the forward primer, and
its reverse read in R2.

The mates are swapped, not reverse-complemented. Reverse-complementing each mate, which
is what MetaFlux did up to v2.5.0, puts both reads on the right strand but leaves R1
covering the 3′ end of the amplicon and R2 the 5′ end. `mergePairs`, with
`amplicon.dada2.merge.trim_overhang` at its default `true`, then trims away whatever one
read extends past the other and keeps only the stretch where the two overlap: on a
418-bp COI amplicon in which every R1 carries the reverse primer, none of the merged
reads came within 9 bp of the full length that way, against 99.4 % with the mates
swapped (99.1 % at exactly 418 bp).

**With `mixed`, both 5′ calls need 10 primer bases.** The two passes read the same
input, and each keeps a pair when its 5′ call finds a primer on both mates. By default
Cutadapt counts a 5′ primer as found as soon as a read starts with the primer's last three
bases (a minimum overlap of 3, because a read may begin partway into its primer), or when
the whole primer turns up anywhere in the read within the error rate. With some primer
pairs the start of each primer resembles the end of the other. The COI primer BF3 ends in
HCG, which the start of a BR2 read matches whenever the read's degenerate third base is a
G (TCG), and BR2 ends in YCA, which the start of a BF3 read matches whenever its third
base is an A (CCA). In the two BF3/BR2 libraries below, 42 % and 61 % of the BR2 reads
start with TCG, and 31 % and 28 % of the BF3 reads with CCA. A pair whose two mates both
start that way passes the call meant for the other orientation too, with only those three
bases cut, and up to v2.5.0 such a pair was written twice: once as it should be, and once
with its mates the wrong way round and most of a primer still on. If DADA2 merged that
wrong copy, the result was the reverse complement of the amplicon, or of part of it, with
primer remnants on.

So with `mixed`, both 5′ calls need at least 10 bases of their primer (Cutadapt's
`--overlap 10`). A read that starts with its primer has all of it, while a chance match
of 10 bases or more is rare. The 3′ calls keep Cutadapt's default, because a read that
runs through a short amplicon into the far primer can legitimately end a few bases into
it. With `fixed` nothing changes: its one pass keeps Cutadapt's default and gives the same
reads as in v2.5.0.

Measured with the rule's own commands at the default settings, each library against
`fixed` run on the same reads with the two files in their right orientation. The columns
give the pairs `fixed` keeps, the pairs `mixed` writes twice with Cutadapt's 3-base
default and with 10 bases, and the pairs that `fixed` keeps and `mixed` with 10 bases
loses:

| Primers | Reads | `fixed` keeps | Twice, 3 bases | Twice, 10 bases | Lost |
|---|---|---|---|---|---|
| COI BF3/BR2 | insect mock (SRR8082166), every R1 carries BR2 | 910,014 | 120,930 (13.3 %) | 1 | 509 (0.06 %) |
| COI BF3/BR2 | another library of the same mock (SRR8082198), R1 carries BF3 | 676,377 | 115,376 (17.1 %) | 0 | 201 (0.03 %) |
| COI fwhF2/fwhR2n | the same mock (SRR8082148), every R1 carries fwhR2n | 1,383,136 | 1,875 (0.14 %) | 0 | 1 |
| COI mlCOIintF-XT/jgHCO2198 | the same mock (SRR8082178), R1 carries mlCOIintF-XT | 1,089,685 | 0 | 0 | 10 (0.001 %) |
| rpoB, Univ_rpoB_deg | [rpoB test set](../about/test-datasets.md), 5 mocks | 157,641 | 487 (0.31 %) | 0 | 141 (0.09 %) |
| gyrB, F64/R353 | [gyrB test set](../about/test-datasets.md), 5 mocks | 318,145 | 113 (0.04 %) | 32 | 136 (0.04 %) |
| 18S, TAReuk454FWD1/TAReukREV3 | [18S test set](../about/test-datasets.md), 4 samples | 908,768 | 2 | 0 | 52 (0.006 %) |
| 16S, 799F/1175R | [16S test set](../about/test-datasets.md), 6 samples | 347,792 | 0 | 0 | 521 (0.15 %) |
| ITS, ITS3/ITS4 | [ITS test set](../about/test-datasets.md), 6 samples | 455,764 | 0 | 0 | 517 (0.11 %) |

The pairs still written twice with 10 bases are pairs in which the other orientation's
primer matches 10 or more of its bases, in 31 of the 33 as a whole primer found inside
one of the reads; the gyrB primers are highly degenerate.

The lost pairs are what the 10 bases cost. In all but 2 of the 2,088 lost pairs, `fixed`
found the primer of a mate from only 3 to 9 of its bases at the start of the read, and
removed just those bases. To check whether such a mate is a cleanly trimmed amplicon read,
a simple test was used: a trimmed read counts as amplicon-like when at least 10 of the
reads that `fixed` keeps from the same sample and in the same direction (R1 or R2) start
with the same 15 bases. Depending on the sample and direction, 93 to 99.8 % of the reads
that `fixed` keeps pass this test, but the mate cut short passes it in only 69 (3 %) of
the lost pairs. In the others it still carries most of a primer that differs from the
given one by more than the error rate allows, or is a read of something other than the
amplicon, such as a low-complexity read. The other mate, by contrast, passes the test in
1,792 of the 2,086 lost pairs (86 %).

What the swap means for DADA2: in a truly mixed library, the R1 file holds sequencer
read 1 of the pairs that came in the right way round and sequencer read 2 of the swapped
ones, and the R2 file the other way round. DADA2 learns one error model for all the R1
reads of a run and one for all the R2 reads ([step 6](#6-asv-inference)), because read 2
is usually the noisier of the two; in a mixed library each model is fitted to a blend of
both reads instead, and the truncation lengths ([step 5](#5-truncation-length-selection))
are chosen from the same blended quality profiles.

Measured on a genuinely mixed library: two BF3/BR2 libraries of the COI insect mock,
sequenced in opposite orientations and joined as one sample of 1.6 million pairs (57 % of
them reversed), run with `mixed`, against the two libraries run separately with `fixed`
(the reversed one with its files swapped) and their tables merged by sequence. The blended
run recovered 363 of the 374 specimens and the separate runs 364 (the one the blend missed
has 3 reads); 600 kept ASVs are shared, of 654 and 664; 99.2 % of the reads left after
chimera removal are the full 418 bp in both; and the two tables are almost identical, with
a Bray-Curtis dissimilarity of 0.011 on specimen proportions (a 0-to-1 measure of how
different two tables are, 0 for identical ones). The blended run keeps 7 % fewer reads,
all lost at DADA2's quality filter: the truncation lengths are chosen on the blended
quality profile, 268 / 270 bases here, where each library on its own got 249 / 271 and
260 / 233, and the longer cuts let more pairs fail the expected-error limit. The error
models themselves converged as quickly as on a single library. This is an upper bound of
what the blend costs, because the two libraries also come from two sequencing runs, so two
error profiles share one file, which a mixed library from one run would not have. With the
trimming of v2.5.0 and earlier the same library gave no result at all: DADA2's error model
had not converged after 172 rounds and 10.5 hours when the run was stopped, because half
of the reads in each file were reverse-complemented, with their quality running backwards,
and 7 % still began with most of a primer (see
[Troubleshooting](../troubleshooting.md#amplicon-dada2s-error-model-never-converges)).

A library in which every pair is swapped needs no second pass: swap the two input files
of each sample (R1 ↔ R2) and keep `fixed`. With `mixed`, the second pass would keep the
pairs that `fixed` keeps on the swapped files, less those the 10 bases cost (509 of
910,014 in SRR8082166, above), and the first pass would add almost nothing (1 pair
there). But Cutadapt would read every pair twice, and the Cutadapt section of the MultiQC
report, which shows the first pass alone, would show nearly every pair as discarded. With
the files swapped, Cutadapt runs one pass and its report describes all the reads.

Measured on the 16S test set with the two files of every sample swapped, `mixed` gave
the taxonomy, the ASV names and the tree of the `fixed` run on the original files. The
trimmed reads differed only by what the 10 bases cost and what the swap brought in: 517
pairs dropped and 4 trimmed 3 bases shorter (together 0.15 % of the pairs), and 22 pairs
that came the other way round in the original files, so `fixed` had dropped them and the
swap put them the right way round. Nearly all of these fail DADA2's quality filter
anyway, and the final ASV table differed by 7 of 193,097 reads, in two ASVs. With the
mates of every second pair swapped instead, the trimmed files held the same pairs in
another order, the first pass's pairs first, and DADA2 gave the same counts; only ASVs
with exactly the same counts can then be numbered in another order. With
`method: sintax`, an ASV that changes its number also changes its place in the
classifier's input, and VSEARCH bases each ASV's random bootstrap draws on the seed and
on that place, so a lineage near the confidence cutoff can gain or lose a rank: here 5 of
226 ASVs, each of 2 to 6 reads, and their counts did not change. Compare such runs by
sequence, not by ASV number.

!!! tip

    Leave `orientation: fixed` unless the data really is mixed. The symptom of a wrong
    setting is a large drop in reads at the primer-trimming stage — visible in
    `stats/read_tracking.txt` and in the Cutadapt section of the MultiQC report.

Per-sample Cutadapt statistics are written as JSON to `stats/cutadapt/` — one file for the
5′ call and one for the 3′ call of the first pass. The second pass writes its Cutadapt
reports to `logs/trim_primers/{sample}.log`, so with `mixed` the Cutadapt section of
the MultiQC report describes the first pass alone, and there the pairs that only the second
pass keeps count as discarded. `stats/03.stripped_reads.counts` and the `stripped` column of
`stats/read_tracking.txt` count the pairs of both passes, so the rare pair kept by both
counts twice.

With `mixed`, `stats/orientation/{sample}.tsv` shows how the pairs of each sample split
between the two orientations: the read pairs that went into primer trimming, the pairs the
first pass kept as they were sequenced, the pairs the second pass kept and swapped, and
each as a percentage of the pairs that went in. The counts are those of the two 5′ calls,
where the orientation of a pair is decided, taken from Cutadapt's own reports. The 3′ calls
after them can still drop a pair that ends up shorter than `min_length`, so the two kept
columns can add up to a little more than the pairs in `03.stripped/`; the difference is the
pairs those calls dropped. The rare pair kept by both passes counts in both columns, as it
counts twice in `03.stripped/`. Nothing in the pipeline reads the file; the
[output reference](../reference/output.md#stats) lists its columns.

The file also makes sure that the swap and the 10-base rule reach output folders trimmed
before they existed. A folder trimmed with `mixed` by MetaFlux v2.5.0 or earlier has no
such file, and a run with `mixed` asks for it, so primer trimming runs again for every
sample, and so does every step after it: the trimmed reads with their QC and read counts,
the truncation lengths, DADA2, the region extraction and the phylogeny where they run,
the taxonomy, the read tracking and the MultiQC report are made again. The PhiX removal
and the QC and counts of the raw reads are kept, and so are the reference files and the
length probe when the run starts from the MetaFlux folder that made the output folder. A
fresh clone, an unpacked release or a second worktree looks for the probe, and by default
the references, in its own `refdb/` folder, so from there the probe is computed again and
the references are fetched again, unless the old `refdb/` was carried over.

Primer trimming runs again whichever MetaFlux folder the run starts from, and whatever
`--rerun-triggers` says. Without the file it would not always happen. Snakemake notices a
changed command only by comparing it with its record of the step (what it wrote down when
the step finished), kept in `.snakemake/` inside the MetaFlux folder it runs from, and a
fresh clone, an unpacked release or a second worktree has no records of the output
folder. From there, DADA2 and the steps after it would still run again, because their
scripts are newer than the old outputs, but on the reads the old version trimmed, and
without a word. A folder trimmed with `fixed` has no such file and is not asked for one,
so nothing runs again because of it.

## 3. Per-stage QC

`falco` profiles reads at each stage of preprocessing: `raw`, `nophix` (only when PhiX
removal ran) and `stripped`. Falco is a drop-in FastQC replacement that reads FASTQ much
faster and writes the same `fastqc_data.txt`, `fastqc_report.html` and `summary.txt`,
one directory per sample, direction and stage under `stats/falco/`.

Running the same QC three times is what makes it possible to see *where* reads were lost
and *how* the quality profile changed, instead of only seeing the final, post-trimming
snapshot. Each Falco call
is fed a stage-renamed symlink of its input, because MultiQC takes the sample name from the
filename recorded inside `fastqc_data.txt` — without that, the three stages of one sample
would collapse into a single entry.

The `stripped`-stage tables are not only a report: `pick_trunclen` reads their per-base
quality columns in the next step.

## 4. Amplicon-length probe

Two later decisions need to know how long the amplicon should be: the truncation lengths
(step 5) must leave R1 and R2 overlapping, and the ASV length filter (step 8) needs a
window of plausible lengths. Rather than hard-coding a number per primer pair, MetaFlux
measures it from the reference database with the primers actually in use.

`amplicon_probe` runs only when `amplicon.expected_length: auto`. It works in one of two
modes, fixed by the marker pack:

| Mode | Markers | What it does |
|------|---------|--------------|
| `pcr` | 16S, 18S, rpoB | In-silico PCR. Two Cutadapt passes over a full-length reference: find the forward primer wherever it sits and cut it away along with everything before it, then, of those sequences, find the reverse complement of the reverse primer and cut it away along with everything after it. What survives between the two primer sites is the predicted amplicon. |
| `direct` | ITS, gyrB | The reference is already cut to amplicon length — UNITE's pre-extracted UCHIME ITS1/ITS2 subregions, or DD7RZ8's pre-trimmed gyrB amplicons — so lengths are read straight off it. In-silico PCR would fail here: there are no primer sites left to cut against. |

Either way, the step writes a JSON holding the length distribution: `n_amplicons`, `min`,
`q1`, `median`, `q3`, `p95`, `p99`, `max`, `mean` and `stdev`. For 16S, chloroplast,
mitochondrial and eukaryote references are left out of these statistics, because their
in-silico products are often longer and the taxon filter removes them later anyway; the
JSON records what was left out (`excluded_taxa`, `n_excluded`). For rpoB, products that
still hold a second forward-primer site are left out too (`n_excluded_second_site`):
their forward primer matched a look-alike stretch upstream of its real site. Which
single statistic feeds the truncation constraint is set by
`amplicon.probe_length_stat.<key>` (the shipped config uses `p95` for 16S, 18S and rpoB).

The result is cached in the repository under `refdb/cache/`, not under `out_dir`, and the
filename encodes the marker, the reference (for 16S including the lineages left out, for
rpoB the second-site rule), and a short hash of both primer FASTAs. Same marker, same
reference, same primers means the cached probe is reused — across reruns and across different output directories. Change a primer and the hash changes, so the probe is
recomputed automatically.

Details, including which statistic to pick and when `auto` is the wrong choice, are on the
[amplicon length and truncation](length-and-truncation.md) page.

## 5. Truncation-length selection

DADA2's `truncLen` cuts every read to a fixed length and **discards any read shorter than
it**. That one number is pulled in opposite directions: cutting shorter removes the
low-quality 3′ tail the error model dislikes, but cutting too short leaves R1 and R2 no
longer overlapping, and pairs that cannot be merged are lost.

`pick_trunclen` writes the decision to `stats/trunclen.json`, which `dada_filter` and
`dada_length_filter` read. In `auto` mode it takes the lower-quartile (Q1) quality at each
cycle from the Falco tables of the `stripped` stage, takes the median across samples per
position, and cuts at the last position still at or above `trunc_len.q_threshold`
(default `20`) — that is, one position before the first drop. Only positions covered by
every sample are considered, so the cut can never exceed the read length of the shortest
sample.

It then checks the merge constraint:

```text
truncLen_R1 + truncLen_R2  >=  expected_length + min_overlap
```

`min_overlap` defaults to `12` bp. If the quality-based cuts break the constraint,
`trunc_len.resolve_policy` decides: `raise_trunc` (the default) extends the cuts past the
quality drop, splitting the shortfall between R1 and R2 according to how much room each
still has; `relax_q` lowers the quality threshold one point at a time down to
`trunc_len.q_floor` (default `15`) and re-cuts; `error` stops and reports the shortfall.

In `manual` mode the quality analysis is skipped and `trunc_len.manual_r1` /
`manual_r2` are used as given.

!!! warning "Clean runs can defeat `auto`"

    If quality never drops below the threshold, the cut lands at the full read length —
    but after primer trimming most reads are a little shorter than that, and
    `filterAndTrim` will then discard almost all of them. The failure appears at step 6 as
    *"No reads passed the filter"*. The fix is to switch to `mode: manual` with values a
    little below the length most reads still reach. See
    [Troubleshooting](../troubleshooting.md).

ITS is the exception: the cuts are still computed and recorded, but `truncLen` is
overridden to `c(0, 0)` at the filtering step, because fungal ITS amplicons vary genuinely
in length and a fixed cut would discard the short ones. Quality control for ITS comes from
`truncQ` and `maxEE` instead. This override also wins over `manual` mode, so `manual_r1` /
`manual_r2` have no effect for ITS.

## 6. ASV inference

This is the DADA2 core: three rules, run in order.

**`dada_quality_plots`** aggregates the per-base quality profile across all primer-trimmed
samples and saves it as PNG and PDF (`stats/dada2/stripped_read_R{1,2}_qual_plot.*`).
Diagnostic only — it does not feed the truncation decision, which comes from Falco via
`pick_trunclen`.

**`dada_filter`** runs DADA2's `filterAndTrim` per sample on the reads from `03.stripped/`,
using the truncation lengths from `stats/trunclen.json`. Filtered reads go to `04.filtered/`,
and a small per-sample JSON of in/out counts goes to `stats/dada2/`, where `dada_seqtab`
collects it for read tracking. The gates:

| Parameter | Config key | Shipped default |
|-----------|------------|-----------------|
| Max expected errors `[R1, R2]` | `amplicon.dada2.max_ee` | `[2, 5]` |
| Adaptive 3′ quality trim | `amplicon.dada2.trunc_q` | `2` |
| Minimum read length | `amplicon.dada2.filter.min_len` | `20` |
| Maximum read length | `amplicon.dada2.filter.max_len` | `null` (no limit) |

Ambiguous bases are not tolerated (`maxN = 0`), and DADA2's own PhiX screen is always left
off here (`rm.phix = FALSE`, hardcoded) — removing PhiX is step 1's job, when step 1 is
enabled. For ITS, and only when a probe has run (`expected_length: auto`), `min_len` is
instead derived from the probe distribution using `amplicon.dada2.filter.min_len_stat`
(default `q1`) — the extracted ITS subregions are shorter than a read, so a probe-derived
floor is meaningful there, whereas for the other markers the probe describes the full
amplicon and would drop every read. Without a probe JSON, ITS falls back to the config
`min_len` like every other marker.

**`dada_seqtab`** runs the denoising chain once over all samples together:

1. `learnErrors` on R1 and R2 separately — builds a run-specific error model from a
   random subsample of reads (`amplicon.dada2.learn_errors.nbases`, default `1e8`).
   Diagnostic plots go to `stats/dada2/filtered_read_R{1,2}_error_plot.*`.
2. `derepFastq` then `dada` — the denoising itself, which decides whether a rare sequence
   is a real variant or a sequencing error of an abundant one.
   `amplicon.dada2.dada.omega_a` (default `1e-40`) sets that sensitivity: it is the
   abundance p-value a unique sequence must beat to be split off as its own ASV, so lower
   values are stricter and yield fewer ASVs. The default is deliberately conservative;
   *raising* it (say to `1e-20`) splits off more ASVs, at the cost of more false
   positives. `amplicon.dada2.pool` (default `false`) controls whether samples are
   processed independently, pseudo-pooled, or fully pooled.
3. `mergePairs` — joins R1 and R2 into full amplicons, requiring
   `amplicon.dada2.merge.min_overlap` bp of overlap (default `12`) with at most
   `max_mismatch` mismatches (default `0`) — both DADA2's own defaults.
4. `makeSequenceTable` then `removeBimeraDenovo` — builds the ASV × sample count table and
   removes chimeras (`amplicon.dada2.chimera.method`, default `consensus`).

!!! note "Two `min_overlap` keys, two different jobs"

    `amplicon.min_overlap` is a *design* constraint used only when choosing truncation
    lengths — it makes sure the reads are cut so they can still meet. `amplicon.dada2.merge.min_overlap`
    is the *runtime* threshold `mergePairs` applies to each individual pair.

The outputs land in `05.dada2/`: `seqs.fasta` (ASV sequences), `seqtab_head_names.txt`
(counts keyed by `ASV_#` ID), `seqtab_head_seqs.txt` (the same counts keyed by the full
sequence) and `read.counts`.

`amplicon.seed` (default `42`) fixes the random-number draw in `learnErrors`, which is what
makes a rerun on the same input reproduce the same ASVs.

## 7. Marker-region extraction *(optional)*

What this step does depends on the marker. An ITS ASV carries fragments of the flanking
rRNA genes — 5.8S, SSU, LSU — outside the marker region proper; ITSx trims them, which
makes ASV lengths directly comparable to the probe distribution and stops near-identical
sequences from being split by their flanks alone. A 16S amplicon, by contrast, lies
entirely inside the 16S gene, since every standard primer binds within it, so Metaxa2 has
little or nothing to trim. There it acts as a filter instead, dropping ASVs it cannot
recognise as SSU rRNA — off-target amplification — which is why the pre- and
post-extraction length distributions normally coincide for 16S.

`target_extract` runs when `amplicon.extraction.enabled` is `true` **and** the marker has an
extractor:

| Marker | Extractor | Notes |
|--------|-----------|-------|
| 16S | Metaxa2 | Run with `-t all` (all domains). Every extracted sequence is kept; pruning is left to the taxonomy filter in step 9, which also uses Metaxa2's own call to remove mitochondria and chloroplasts (`metaxa2_origin.tsv`). |
| ITS | ITSx | Run with `-t all --only_full F`. Partial detections must be kept: a targeted ITS1 or ITS2 amplicon never contains the full ITS region. The region kept is the one named by `amplicon.its_region`. |
| 18S, gyrB, rpoB, COI | *none* | No extractor exists in MetaFlux for these. Extraction is forced off with a warning if the config still asks for it. |

Both extractors write `seqs_extracted.fasta` and `seqtab_extracted_head_names.txt` into
`05.dada2/`, alongside the extractor's own report (`metaxa2_extraction.results.txt` or
`itsx_extraction.summary.txt`) and its full raw output directory, kept for inspection.
Metaxa2 also writes `metaxa2_origin.tsv`, its origin call for each ASV, which the
taxonomy filter reads in step 9. ASVs the extractor could not resolve are dropped at
this point.

ITSx also collapses ASVs: distinct ASVs whose extracted subregions come out identical —
they differed only in the flanks that were just removed — become a single representative,
and their per-sample counts are summed. The mapping is recorded in
`itsx_collapse_map.tsv`.

## 8. ASV length filter

Off-target amplification and mis-merged pairs tend to show up as ASVs of implausible
length. `dada_length_filter` keeps only ASVs inside a length window and writes
`seqs_lenfilt.fasta`, `seqtab_lenfilt_head_names.txt` and
`seqtab_lenfilt_head_seqs.txt` to `05.dada2/`.

It is not a backstop for chimeras. A bimera carries both parents' primer sites and spans
the same region, so it comes out the same length as a real amplicon and passes the window
untouched. Chimeras are removed only by `removeBimeraDenovo` at step 6, and the chimeric
fraction in `stats/read_tracking.txt` is what to watch.

The filter runs *after* extraction, so the lengths it sees are the trimmed marker region —
directly comparable to what the probe measured. The window comes from, in order of
priority:

1. `amplicon.length_filter.range`, when set — an explicit `[min, max]`.
2. `mode: auto` with a probe available — `[probe q1 − window_margin, probe p95 + window_margin]`,
   with `window_margin` defaulting to `50` bp. P95 rather than P99 sets a tighter ceiling.
   The ITS probe is the case to watch: it reads lengths straight off pre-extracted
   ITS1/ITS2 subregions, so its upper tail is real length variation between fungal
   lineages, and a p95 ceiling will clip genuinely long ITS ASVs, Basidiomycota most
   of all. For ITS, check
   `stats/dada2/asv_length_hist.png` against the window and set an explicit
   `length_filter.range` if long-ITS taxa matter.
3. `mode: auto` with no probe — a fallback of ±15% around `expected_length`.

Two diagnostics go to `stats/dada2/`. `asv_length_stats.json` always records three
distributions — pre-extraction, filter source, and kept. `asv_length_hist.png` overlays the
pre-extraction series only when extraction ran; when no extractor applies, the filter
source is simply the DADA2 output. The gap between the pre-extraction and filter-source
distributions is exactly how much Metaxa2 or ITSx trimmed; for ITS that gap is usually
large, for 16S the two normally coincide. If the window ends up excluding everything, the
rule stops with an error rather than writing an empty table.

The window logic is covered in full on the
[amplicon length and truncation](length-and-truncation.md) page.

## 9. Taxonomy assignment

`assign_taxonomy` classifies the length-filtered ASVs and writes the final tables to
`06.taxonomy/`. Two implementations exist, and one is chosen at start-up based on
`amplicon.taxonomy.method`:

| Method | Tool | Threshold key | Available for |
|--------|------|---------------|---------------|
| `rdp` | DADA2 `assignTaxonomy` (naive Bayesian RDP classifier), plus `addSpecies` for 16S | `min_boot` (default `80`) | 16S, ITS, 18S, gyrB, rpoB |
| `sintax` | `vsearch --sintax` | `sintax_cutoff` (default `0.8`, roughly equivalent to `min_boot: 80`) | 16S, ITS, 18S, COI |

gyrB and rpoB are RDP-only: their reference releases ship a DADA2 training set and no
SINTAX build, so `method: sintax` is rejected at start-up with a clear message rather than
failing later. COI is SINTAX-only, for the opposite reason: a DADA2 classifier trained on
MIDORI2 would need about 225 GB of memory at species level, so `method: rdp` is rejected
in the same way.

The reference database follows from the marker — SILVA for 16S, UNITE for ITS, PR2 for 18S,
DD7RZ8 for gyrB, FROGS RefSeq for rpoB, MIDORI2 for COI — and is fetched on first use. Some of them arrive
by way of a local conversion step rather than a plain download: for 16S,
`convert_silva_train` writes SILVA 144's training set with its kingdom names without
SILVA's `--other` suffix and `convert_silva_sintax` builds the SINTAX database from it;
for 18S under `sintax`, `convert_pr2_sintax` moves the domain of PR2's UTAX file to the
SINTAX letter `d:`; and the rpoB reference is rebuilt from the FROGS release by
`convert_rpob_to_dada2`. See the [marker pages](markers/index.md) for what each one
contains.

Before any classification, `check_taxonomy_reference` reads the exact file the
classifier will use and checks it against the marker pack: the number of ranks (or the
SINTAX letters and their order), the rank prefixes, and the names in the top ranks. A
wrong file, such as another SILVA release's training set pinned in the run config,
stops the run there, before a name can land in the wrong column. A file that passes
gets a record, `stats/taxonomy_reference_check.json`, and only then does
`assign_taxonomy` start. See [the reference check](markers/index.md#the-reference-check-reference_values).

Both paths write the same four tables; the `sintax` path additionally leaves
`_sintax_raw.tsv`, the raw VSEARCH output, for troubleshooting.

| File | Rows | Columns |
|------|------|---------|
| `asv_table.txt` | `ASV_#` IDs | per-sample counts + taxonomy |
| `asv_table_seqs.txt` | ASV sequences | per-sample counts + taxonomy |
| `taxon_seq_table.txt` | `ASV_#` IDs | one column per rank + the sequence |
| `asv_table_unfiltered.txt` | `ASV_#` IDs, every ASV that reached the filter | per-sample counts + taxonomy + `metaxa2_origin` + `filter_status` |

The taxonomy string is built from the marker's rank model, which follows the ranks of
the reference in use: `d__` is always the domain and `k__` always the kingdom, and a
string has a kingdom only where the reference defines one. So 16S on SILVA 144 has eight
ranks (Domain, Kingdom, Phylum … Species) and on SILVA 138.2 seven, with no kingdom; ITS
has seven starting at the kingdom, because UNITE has no domain; rpoB has seven starting
at the domain; 18S has PR2's nine ranks under `method: rdp` and eight under
`method: sintax`, where PR2's UTAX file merges Division and Subdivision into one field.
gyrB is seven ranks deep too, but its first rank is the DD7RZ8 paralog tag — `gyrB` for
the true gene, `other` for a co-amplified paralog — so a gyrB lineage has no domain or
kingdom at all. See the [gyrB page](markers/gyrb.md), and
[taxonomy strings in the final tables](../reference/output.md#taxonomy-strings-in-the-final-tables)
for every marker's slots side by side.

The last thing the step does is apply the keep/discard contaminant filter
(`amplicon.taxonomy.filter`), which matches whole rank segments of each ASV's lineage.
`discard` also checks the classifier's best hit below the reported ranks, so an ASV whose
order or family call fell under the cutoff can still be removed
(`filter.discard_min_support`). The filter has no built-in default for the lists: left
unset, both resolve to `[]` and match nothing.
The lists have to match the marker being run. A token that this run's strings can never
contain (`k__Bacteria`, or any `d__` token on ITS) stops the run at start-up, before any
step; a token that the marker does have but the study does not want, such as a 16S
`keep: [d__Bacteria, d__Archaea]` on an 18S run, passes and removes nearly every ASV.
For 16S the filter also removes, by default, the ASVs that Metaxa2
called mitochondria or chloroplast in step 7 (`filter.organelles`), whatever the lists say;
every ASV, removed ones included, stays listed in `asv_table_unfiltered.txt`. It is
described on the [keeping and discarding taxa](taxon-filter.md) page.

## Read tracking and the QC report

Two rules close the run.

`aggregate_read_counts` joins the per-stage counts into `stats/read_tracking.txt`: one row
per sample, one column per stage, in pipeline order — `raw`, `nophix` (only when PhiX
removal ran), `stripped`, `filtered`, `denoised`, `merged`, `non_chimeric`,
`post_extraction` (only when extraction ran), `post_length_filter` and
`post_taxonomy_filter`. This is the first file to open when a run produces fewer ASVs than
expected: it shows which step consumed the reads.

`multiqc` collects the Falco reports, the Cutadapt JSONs and the bowtie2 PhiX logs into
`multiqc/multiqc_report.html`.

## What changes between markers

Everything above is the same for all six markers except the four columns below, and the
two COI screens (a reading-frame check and an off-target guard against bacterial COX1)
that run before the taxon filter in COI runs.

| Marker | Probe substrate and mode | Extractor | Taxonomy reference | Methods |
|--------|--------------------------|-----------|--------------------|---------|
| [16S](markers/16S.md) | SILVA, in-silico PCR | Metaxa2 | SILVA (+ species assignment) | rdp, sintax |
| [ITS](markers/its.md) | UNITE UCHIME ITS1/ITS2, direct | ITSx | UNITE | rdp, sintax |
| [18S](markers/18S.md) | SILVA-Euk, in-silico PCR | none | PR2 | rdp, sintax |
| [gyrB](markers/gyrb.md) | DD7RZ8, direct | none | DD7RZ8 | rdp only |
| [rpoB](markers/rpob.md) | FROGS RefSeq, in-silico PCR | none | FROGS RefSeq | rdp only |
| [COI](markers/coi.md) | MIDORI2, in-silico PCR | none | MIDORI2 (converted) | sintax only |

Start with the [marker pages](markers/index.md) before setting up a run: each one covers
the primers the reference was built for, the filter tokens that make sense for it, and the
quirks that follow from its database.

## Where to go next

- [Amplicon length and truncation](length-and-truncation.md) — steps 4, 5 and 8 in depth.
- [Keeping and discarding taxa](taxon-filter.md) — the step-9 contaminant filter.
- [Markers](markers/index.md) — per-marker references, primers and defaults.
- [Output files](../reference/output.md) — what every file under `out_dir` contains.
- [Configuration](../reference/configuration.md) — the full config key reference.
- [Troubleshooting](../troubleshooting.md) — the failure modes worth recognising early.
