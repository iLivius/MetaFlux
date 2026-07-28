# Amplicon mode

Amplicon mode takes raw paired-end Illumina metabarcoding reads and returns a table of
exact amplicon sequence variants (ASVs) with a taxonomic assignment for each one. The
denoising engine is DADA2; five markers are supported — 16S rRNA, ITS, 18S rRNA, gyrB and
rpoB — and the marker is what decides which reference databases are used, whether a
region extractor runs, and which classifiers are available.

Two config keys define the run:

```yaml
mode: amplicon
amplicon:
  type: 16S          # 16S | ITS | 18S | gyrB | rpoB
```

Everything else is either a path (input reads, output directory, primer FASTAs) or
tuning. Marker-specific facts — database URLs, rank models, which extractor applies —
live in the marker packs under `workflow/markers/<type>.yaml` and are never set by hand.

## The steps at a glance

The pipeline runs as nine stages. Output directories are numbered in pipeline order under
`out_dir`, so the directory listing reads like the workflow.

| # | Step | Tool | Snakemake rule | Writes to |
|---|------|------|----------------|-----------|
| 1 | PhiX removal *(optional)* | bowtie2 | `rm_phix` | `2.no_phix/` |
| 2 | Primer trimming | Cutadapt | `trim_primers` | `3.stripped/` |
| 3 | Per-stage QC | Falco | `falco` | `stats/falco/` |
| 4 | Amplicon-length probe | Cutadapt (or direct read) | `amplicon_probe` | `refdb/cache/` |
| 5 | Truncation-length selection | in-house | `pick_trunclen` | `stats/trunclen.json` |
| 6 | ASV inference | DADA2 | `dada_filter`, `dada_seqtab` | `4.filtered/`, `5.dada2/` |
| 7 | Marker-region extraction *(optional)* | Metaxa2 or ITSx | `target_extract` | `5.dada2/` |
| 8 | ASV length filter | in-house | `dada_length_filter` | `5.dada2/`, `stats/dada2/` |
| 9 | Taxonomy assignment | DADA2 RDP or VSEARCH SINTAX | `assign_taxonomy` | `6.taxonomy/` |

Two more rules close the run: `aggregate_read_counts` builds `stats/read_tracking.txt`,
and `multiqc` builds `multiqc/multiqc_report.html`.

## Input

Amplicon mode requires local paired FASTQs named `{sample}_R1.fastq.gz` and
`{sample}_R2.fastq.gz` in `input.fastq_dir`. Sample names are discovered from those
filenames at parse time — the moment before the run starts, when Snakemake reads the
config and works out which jobs it will need. SRA-style `{sample}_1/_2` naming is
accepted by shotgun mode but rejected here with an explicit error, because the two
conventions cannot be told apart safely once mixed.

The raw files are never modified. `link_reads` symlinks them into `1.reads/` purely so the
inputs of a run are visible next to its outputs; every rule that needs raw reads reads
them from `fastq_dir` directly.

## 1. PhiX removal *(optional)*

Illumina runs are routinely spiked with the PhiX174 genome as a sequencing control, and
those reads have no place in a community profile. `rm_phix` maps each pair against a
bowtie2 index of the PhiX reference and keeps the pairs that fail to align *concordantly*
(`--un-conc-gz`), which are written to `2.no_phix/`. Failing to align concordantly is not
the same as not matching at all: a pair in which only one mate aligned to PhiX is also
written out and kept, so a small PhiX residue can survive. The shotgun path is stricter —
BBDuk matches by k-mer and, with its default `removeifeitherbad=t`, drops a pair as soon
as either mate matches.

The reference is fetched once from NCBI (`references.phix.fetch_url`, RefSeq assembly
GCF_000819615.1) by `fetch_phix` and indexed by `build_phix_index`; both are cached under
`refdb/phix/`, so this cost is paid on the first run only.

Toggled with `amplicon.decontamination.remove_phix` (shipped default `true`). When it is
`false` the step is skipped entirely: `2.no_phix/` is not created, Cutadapt reads the raw
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
PhiX removal is off) and writes to `3.stripped/`.

Trimming happens in two calls, because a read can carry primer sequence at both ends:

| Call | Cutadapt options | Purpose |
|------|------------------|---------|
| 5′ | `-g file:fwd -G file:rev`, `--discard-untrimmed` | Removes the forward primer from R1 and the reverse primer from R2. Pairs without a recognisable primer are dropped. |
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
input with the primers swapped (`-g rev -G fwd`, then `-a fwd_rc -A rev_rc`). Its output is
reverse-complemented with seqtk so it matches the orientation of the first pass, and then
concatenated onto it. Both 5′ calls require a recognisable primer on each mate and drop
the pair otherwise, so a correctly oriented pair is claimed only by the first pass and a
swapped pair only by the second — nothing is counted twice.

!!! tip

    Leave `orientation: fixed` unless the data really is mixed. The symptom of a wrong
    setting is a large drop in reads at the primer-trimming stage — visible in
    `stats/read_tracking.txt` and in the Cutadapt section of the MultiQC report.

Per-sample Cutadapt statistics are written as JSON to `stats/cutadapt/` — one file for the
5′ call and one for the 3′ call of the first pass. The swap-pass calls are logged but not
reported separately.

## 3. Per-stage QC

`falco` profiles reads at each stage of preprocessing: `raw`, `nophix` (only when PhiX
removal ran) and `stripped`. Falco is a drop-in FastQC replacement that reads FASTQ much
faster and writes the same `fastqc_data.txt`, `fastqc_report.html` and `summary.txt`,
one directory per sample, direction and stage under `stats/falco/`.

Running the same QC three times is what makes it possible to see *where* reads were lost
and *how* the quality profile changed, instead of only seeing the endpoint. Each Falco call
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
`q1`, `median`, `q3`, `p95`, `p99`, `max`, `mean` and `stdev`. Which single statistic
feeds the truncation constraint is set by `amplicon.probe_length_stat.<key>` (the shipped
config uses `p95` for 16S, 18S and rpoB).

The result is cached in the repository under `refdb/cache/`, not under `out_dir`, and the
filename encodes the marker, the reference, and a short hash of both primer FASTAs. Same
marker, same reference, same primers means the cached probe is reused — across reruns and
across different output directories. Change a primer and the hash changes, so the probe is
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

**`dada_filter`** runs DADA2's `filterAndTrim` per sample on the reads from `3.stripped/`,
using the truncation lengths from `stats/trunclen.json`. Filtered reads go to `4.filtered/`,
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

The outputs land in `5.dada2/`: `seqs.fasta` (ASV sequences), `seqtab_head_names.txt`
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
| 16S | Metaxa2 | Run with `-t all` (all domains). Every extracted sequence is kept; domain-level pruning is left to the taxonomy filter in step 9. |
| ITS | ITSx | Run with `-t all --only_full F`. Partial detections must be kept: a targeted ITS1 or ITS2 amplicon never contains the full ITS region. The region kept is the one named by `amplicon.its_region`. |
| 18S, gyrB, rpoB | *none* | No extractor exists in MetaFlux for these. Extraction is forced off with a warning if the config still asks for it. |

Both extractors write `seqs_extracted.fasta` and `seqtab_extracted_head_names.txt` into
`5.dada2/`, alongside the extractor's own report (`metaxa2_extraction.results.txt` or
`itsx_extraction.summary.txt`) and its full raw output directory, kept for inspection.
ASVs the extractor could not resolve are dropped at this point.

ITSx has one extra behaviour worth knowing: distinct ASVs whose extracted subregions come
out identical — they differed only in the flanks that were just removed — are collapsed
into a single representative and their per-sample counts summed. The mapping is recorded in
`itsx_collapse_map.tsv`.

## 8. ASV length filter

Off-target amplification and mis-merged pairs tend to show up as ASVs of implausible
length. `dada_length_filter` keeps only ASVs inside a length window and writes
`seqs_lenfilt.fasta`, `seqtab_lenfilt_head_names.txt` and
`seqtab_lenfilt_head_seqs.txt` to `5.dada2/`.

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
`6.taxonomy/`. Two implementations exist and one is chosen when the workflow is parsed,
based on `amplicon.taxonomy.method`:

| Method | Tool | Threshold key | Available for |
|--------|------|---------------|---------------|
| `rdp` | DADA2 `assignTaxonomy` (naive Bayesian RDP classifier), plus `addSpecies` for 16S | `min_boot` (default `80`) | all five markers |
| `sintax` | `vsearch --sintax` | `sintax_cutoff` (default `0.8`, roughly equivalent to `min_boot: 80`) | 16S, ITS, 18S |

gyrB and rpoB are RDP-only: their reference releases ship a DADA2 training set and no
SINTAX build, so `method: sintax` is rejected at startup with a clear message rather than
failing later.

The reference database follows from the marker — SILVA for 16S, UNITE for ITS, PR2 for 18S,
DD7RZ8 for gyrB, FROGS RefSeq for rpoB — and is fetched on first use. Two of them arrive
by way of a local conversion step rather than a plain download: the 16S SINTAX database is
built from the SILVA DADA2 trainset by `convert_silva_sintax`, and the rpoB reference is
rebuilt from the FROGS release by `convert_rpob_to_dada2`. See the
[marker pages](markers/index.md) for what each one contains.

Both paths write the same three tables; the `sintax` path additionally leaves
`_sintax_raw.tsv`, the raw VSEARCH output, for troubleshooting.

| File | Rows | Columns |
|------|------|---------|
| `asv_table.txt` | `ASV_#` IDs | per-sample counts + taxonomy |
| `asv_table_seqs.txt` | ASV sequences | per-sample counts + taxonomy |
| `taxon_seq_table.txt` | `ASV_#` IDs | one column per rank + the sequence |

The taxonomy string is built from the marker's rank model, so it has the right depth for
the reference in use — the Linnaean seven for 16S, ITS and rpoB, and nine ranks for
PR2/18S under `method: rdp` (PR2's UTAX file merges Division and Subdivision into one
field, so `method: sintax` yields eight). gyrB is seven ranks deep too, but its first rank
is the DD7RZ8 paralog tag — `gyrB` for the true gene, `other` for a co-amplified paralog —
rather than a kingdom, so a gyrB lineage carries no `k__` segment at all. See the
[gyrB page](markers/gyrb.md).

The last thing the step does is apply the keep/discard contaminant filter
(`amplicon.taxonomy.filter`), which matches whole rank segments of each ASV's lineage. It
has no built-in default: left unset, both lists resolve to `[]` and nothing is filtered.
The lists have to match the marker being run — a 16S keep list applied to an 18S run would
discard every ASV. It is described on the
[keeping and discarding taxa](taxon-filter.md) page.

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

Everything above is the same for all five markers except the four columns below.

| Marker | Probe substrate and mode | Extractor | Taxonomy reference | Methods |
|--------|--------------------------|-----------|--------------------|---------|
| [16S](markers/16S.md) | SILVA, in-silico PCR | Metaxa2 | SILVA (+ species assignment) | rdp, sintax |
| [ITS](markers/its.md) | UNITE UCHIME ITS1/ITS2, direct | ITSx | UNITE | rdp, sintax |
| [18S](markers/18S.md) | SILVA-Euk, in-silico PCR | none | PR2 | rdp, sintax |
| [gyrB](markers/gyrb.md) | DD7RZ8, direct | none | DD7RZ8 | rdp only |
| [rpoB](markers/rpob.md) | FROGS RefSeq, in-silico PCR | none | FROGS RefSeq | rdp only |

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
