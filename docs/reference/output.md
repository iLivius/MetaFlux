# Output files

Everything a run produces goes under the directory named in `output.out_dir`.
Two things are deliberately written elsewhere: reference databases, which are
cached under `references.refdb_root` (`refdb/` by default) so that every run
reuses the same downloads, and — in shotgun BioProject mode — the FASTQs fetched
from SRA, which land in `input.fastq_dir`.

Inside `out_dir` the directories are numbered, with two digits, in the order the
stages run, so a plain `ls` reads like the pipeline itself. Both modes also write
`logs/` (one log per rule, one per sample where the rule runs per sample),
`stats/read_tracking.txt` and `multiqc/multiqc_report.html`.

!!! warning "One output directory per mode"

    Amplicon and shotgun runs use different numbered directories, but they share
    the names `stats/read_tracking.txt` and `multiqc/multiqc_report.html`. Pointing
    both modes at the same `out_dir` overwrites those two files. Give each mode its
    own `out_dir`.

---

## Amplicon mode

### Layout

```text
out_dir/
├── 01.reads/         # symlinks to the input FASTQs (visibility only)
├── 02.no_phix/       # PhiX-filtered reads      (absent when remove_phix: false)
├── 03.stripped/      # primer-trimmed reads
├── 04.filtered/      # DADA2 filterAndTrim output
├── 05.dada2/         # ASV sequences and count tables at each stage
├── 06.taxonomy/      # asv_table.txt, asv_table_seqs.txt, taxon_seq_table.txt,
│                     # asv_table_unfiltered.txt (every ASV, removed ones included)
│                     # (+ _sintax_raw.tsv on the sintax path)
├── 07.phylogeny/     # alignment + unrooted tree  (only when phylogeny.enabled: true)
├── stats/            # read counts, falco/, cutadapt/, trunclen.json, dada2/,
│                     # taxonomy_reference_check.json
│                     # (+ phylogeny/ when phylogeny.enabled: true)
├── multiqc/          # multiqc_report.html + multiqc_data/
├── logs/             # per-rule logs
├── benchmarks/       # runtime/memory for rm_phix and trim_primers
└── _aux/primers/     # reverse-complemented primers used by Cutadapt
```

Three files carry the result of an amplicon run: `06.taxonomy/asv_table.txt`
(ASV × sample counts plus a taxonomy string), `05.dada2/seqs_lenfilt.fasta` (the
sequences those ASV IDs refer to), and `stats/read_tracking.txt` (how many reads
survived each stage, per sample). Next to them, `06.taxonomy/asv_table_unfiltered.txt`
is the same table before the taxon filter, with the step that removed each ASV.
Everything else is intermediate data or QC.

### Reads: `01.reads/` → `04.filtered/`

| Directory | Files | Rule | What happened |
|---|---|---|---|
| `01.reads/` | `{sample}_R1.fastq.gz`, `_R2` | `link_reads` | Symlinks pointing back at `input.fastq_dir`. Nothing reads from here — the downstream rules open the originals. The links exist so the run directory shows which files went in. |
| `02.no_phix/` | `{sample}_R1.fastq.gz`, `_R2` | `rm_phix` | Bowtie2 against the fetched PhiX reference; the pairs that do *not* align are kept. Present only when `amplicon.decontamination.remove_phix` is `true` (the default). |
| `03.stripped/` | `{sample}_R1.fastq.gz`, `_R2` | `trim_primers` | Cutadapt removes the primers at the 5' end and any read-through of the opposite primer at the 3' end. With `primers.orientation: mixed` a second pass recovers reads that arrived the other way round and reverse-complements them back. |
| `04.filtered/` | `{sample}_R1_filt.fastq.gz`, `_R2_filt` | `dada_filter` | DADA2 `filterAndTrim` with the `truncLen` decided in `stats/trunclen.json`. This is the last per-sample FASTQ; from here DADA2 works on all samples at once. |

### Inside `05.dada2/`

The directory looks busy, but most of it is the same two artifacts — the **ASV
sequences** and the **ASV count table** — written out again after each processing
step. The table reads left to right in pipeline order. Within the count tables,
`head_names` columns are `ASV_#` identifiers and `head_seqs` columns are the full
ASV sequence; the counts are identical, only the labels differ.

| Artifact | DADA2 | + Extraction *(optional)* | + Length filter |
|---|---|---|---|
| **ASV sequences** (FASTA) | `seqs.fasta` | `seqs_extracted.fasta` | `seqs_lenfilt.fasta` |
| **Counts** — `ASV_#` keyed | `seqtab_head_names.txt` | `seqtab_extracted_head_names.txt` | `seqtab_lenfilt_head_names.txt` |
| **Counts** — sequence keyed | `seqtab_head_seqs.txt` | — | `seqtab_lenfilt_head_seqs.txt` |
| **Extras** | `read.counts` | Metaxa2 / ITSx results and working dir | — |

The rules behind the columns are `dada_seqtab`, `target_extract` and
`dada_length_filter`, in that order.

The last populated column is the final result. `assign_taxonomy` reads two of these
files: the length-filtered FASTA `seqs_lenfilt.fasta` and the `ASV_#`-keyed
`seqtab_lenfilt_head_names.txt`. For 16S with extraction on it also reads Metaxa2's
call for each ASV, `metaxa2_origin.tsv` (see the extraction trail below), for the
organelle filter. (The read-tracking step also opens
`read.counts` and the `head_names` seqtabs, but only to count reads per stage.)
The sequence-keyed tables — `seqtab_head_seqs.txt` and
`seqtab_lenfilt_head_seqs.txt` — are terminal deliverables: the same counts
relabelled by sequence, convenient for merging runs or feeding another tool, but
no rule in the workflow reads them.

With extraction switched off, the middle column does not exist and the length
filter runs straight on the DADA2 output. Extraction is off whenever
`amplicon.extraction.enabled` is `false`, and is forced off (with a warning) for
markers that have no region extractor — 18S, gyrB and rpoB.

**The extras.** `read.counts` is a TSV written by `dada_seqtab`: one row per
sample, columns `stripped`, `filtered`, `denoised`, `merged`, `non_chimeric`. It
is the raw material for the DADA2 part of `stats/read_tracking.txt`.

The extraction step leaves a per-extractor trail:

| Marker | Extractor | Files in `05.dada2/` |
|---|---|---|
| 16S | Metaxa2 | `metaxa2_extraction.results.txt` (Metaxa2's own summary), `metaxa2_origin.tsv` (Metaxa2's origin call per extracted ASV, read by the organelle filter) plus the raw `seqs.*` working files under `metaxa2/` |
| ITS | ITSx | `itsx_extraction.summary.txt`, `itsx_collapse_map.tsv`, plus the raw `seqs.*` working files under `_itsx_tmp/` |

`itsx_collapse_map.tsv` is worth knowing about. Two ASVs can differ only in the
flanking 5.8S/LSU sequence; once ITSx trims those flanks away the remaining ITS1
or ITS2 region is identical, so the ASVs are merged into one and their counts are
added. The map records which original `ASV_#` went into which representative, so
that merge is auditable rather than silent.

### `06.taxonomy/`

Three tables, all written by `assign_taxonomy`, all built from the same filtered
set of ASVs. They differ only in what labels the rows. A fourth,
`asv_table_unfiltered.txt`, is written by the same rule and holds every ASV that
reached the filter, the removed ones included (see below).

With `taxonomy.method: sintax` — the shipped default — the directory also holds
`_sintax_raw.tsv`, the raw VSEARCH `--tabbedout` result the parser was fed. It is
kept for troubleshooting and no rule reads it. Above one thread its lines come in the
order VSEARCH finishes them, so sort it before comparing two runs; the tables built from
it do not depend on that order.

| File | Rows | Columns |
|---|---|---|
| `asv_table.txt` | `ASV_#` | one per sample, then `taxonomy` |
| `asv_table_seqs.txt` | the ASV sequence | one per sample, then `taxonomy` |
| `taxon_seq_table.txt` | `ASV_#` | one per taxonomic rank, then `sequence` |
| `asv_table_unfiltered.txt` | `ASV_#`, every ASV that reached the filter | one per sample, then `taxonomy`, `metaxa2_origin`, `filter_status` |

`asv_table_unfiltered.txt` is `asv_table.txt` before the filter, with two more
columns. `metaxa2_origin` is Metaxa2's call for the ASV (`mitochondria`,
`chloroplast`, `bacteria`, `archaea`, `eukaryota`, `uncertain`), or `NA` when
Metaxa2 did not run (markers other than 16S, or extraction off) or made no call.
`filter_status` is `kept` or the first filter step that removed the ASV, in the order
the steps run: `removed_keep`, `removed_discard`, `removed_best_hit`,
`removed_organelle`. The `kept` rows are exactly `asv_table.txt`, and the read
counts of all rows add up to the `post_length_filter` column of
`stats/read_tracking.txt`. Use it when the removed ASVs matter, for instance to
compare how many organelle reads two primer pairs pick up — see
[keeping and discarding taxa](../amplicon/taxon-filter.md#the-table-with-everything-asv_table_unfilteredtxt).
Nothing in the pipeline reads it; the tree and every other table use the filtered set.

The rank columns in `taxon_seq_table.txt` come from the marker pack, so they
match the reference's actual depth: `Domain, Kingdom, Phylum … Species` for 16S on
SILVA 144, `Domain, Phylum … Species` for 16S on SILVA 138.2 and for rpoB,
`Kingdom … Species` for ITS, `Gene, Phylum … Species` for gyrB, and for 18S PR2's
nine ranks
(`Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species`)
under `method: rdp`, or eight under `method: sintax`, where PR2's UTAX file joins
Division and Subdivision into one `Division_Subdivision` column. The same column can
therefore sit at a different position from one marker or release to the next: pick the
columns by name (see [reading the tables in R](#reading-the-tables-in-r)).

!!! note "File format"

    The four tables written by `assign_taxonomy` follow R's
    `write.table(col.names = NA)` conventions, and the rdp and the sintax path write
    them in exactly the same format: tab-separated, the header line starts with an
    empty quoted cell (`""`) for the row-name column, and text values (row names,
    taxonomy strings, rank names) are wrapped in double quotes. A missing value is an
    unquoted `NA`, as R writes it: a rank the classifier did not reach in
    `taxon_seq_table.txt`, and a missing `metaxa2_origin`.
    `read.table(..., sep = "\t", header = TRUE, row.names = 1)` in R and
    `pd.read_csv(..., sep = "\t", index_col = 0)` in pandas both read them without
    extra arguments, `NA` as a missing value.

The contaminant filter (`amplicon.taxonomy.filter.keep` / `.discard`, the best-hit
check and, for 16S, the organelle filter) has already been applied when these tables
are written, so the ASVs it removed are absent here, except from
`asv_table_unfiltered.txt` — see [keeping and discarding taxa](../amplicon/taxon-filter.md).

### `07.phylogeny/`

Only when `amplicon.phylogeny.enabled: true` (16S only). The whole directory is absent
otherwise. See [Phylogeny](../amplicon/phylogeny.md).

| File | What it is | What reads it |
|---|---|---|
| `asv_16s.unrooted.nwk` | **The deliverable.** Newick tree of the contaminant-filtered ASVs, tips labelled `ASV_N`. Unrooted on purpose — root it in R after pruning. | Nothing in MetaFlux. Your R session. |
| `asv_16s.aln.fasta` | The MAFFT alignment the tree was built from. Retained so it can be inspected or reused. | `phylo_tree`, `phylo_export`, `phylo_qc`. |
| `asv_16s.fasta` | The extracted ASV sequences, bare `ASV_N` headers — the module's input. | `phylo_align`. |
| `phylogeny.params.json` | Everything resolved at run time: aligner strategy, model, thread counts, tool versions, seed, input checksums, `extra_args`. | `phylo_qc`; and you, when you need to pin a choice. |
| `input.json` | Eligible ASV count, per-ASV lengths, input checksums, and two ordered ID lists: `asv_ids` as the taxonomy table lists them, and `fasta_order`, the sequence-sorted order the alignment was actually built in. Also what decides whether a tree is attempted. | `phylo_qc`, `phylo_export`. |
| `aln_run.json` | Per-step provenance for the alignment, merged into `phylogeny.params.json`. | `phylo_export`. |
| `<backend>/` | The backend's own files — `iqtree/`, `fasttree/` or `raxml-ng/` — plus `tree_run.json`. Inspection only. | `phylo_export` reads the tree and `tree_run.json`; the rest is for you. |

Two things to know about what is and is not here:

- **With fewer than 4 eligible ASVs, only `asv_16s.fasta` and `input.json` are written.**
  No tree is built and the run still succeeds — `stats/phylogeny/phylogeny_qc.json`
  records why. Below four tips there is only one possible unrooted topology.
- **Switching `backend` leaves the previous backend's directory in place**, deliberately,
  so results can be compared. The `phylo_export` log notes when it finds one. The
  canonical files at the top level always come from the backend configured *now*.

!!! note "Rebuilding a deleted tree"
    `stats/phylogeny/phylogeny_qc.json` is the file `rule all` requests; the tree and
    alignment are built as its dependencies. So deleting `asv_16s.unrooted.nwk` on its own
    does not trigger a rebuild. Delete the QC report as well, or run
    `snakemake -R phylo_export`.

### `stats/`

| Path | Contents |
|---|---|
| `read_tracking.txt` | The per-stage summary, described below. |
| `01.raw_reads.counts` | R1 read count per sample, `{sample} : N`, with a trailing total. |
| `02.nophix_reads.counts` | Same, after PhiX removal. Only when `remove_phix: true`. |
| `03.stripped_reads.counts` | Same, after primer trimming. |
| `trunclen.json` | The truncation decision: what `truncLen` was chosen for R1 and R2 and how. See [amplicon length and truncation](../amplicon/length-and-truncation.md). |
| `taxonomy_reference_check.json` | What the reference check saw in the file the classifier reads: the file, the number of records, the ranks per record, the SINTAX letters (under `sintax`) and the names found in each rank the marker pack has a rule for. Written by `check_taxonomy_reference` before any classification, and only when the file passes; `assign_taxonomy` waits for it, and nothing reads its content. See [the reference check](../amplicon/markers/index.md#the-reference-check-reference_values). |
| `falco/{sample}_R{1,2}_{stage}/` | Falco QC per sample, per read direction, per stage (`raw`, `nophix`, `stripped`), each holding `fastqc_data.txt`, `fastqc_report.html` and `summary.txt`. |
| `cutadapt/{sample}.passA_5prime.cutadapt.json` and `.passA_3prime.` | Cutadapt's own report for the 5' and 3' trimming passes. These feed MultiQC. The extra swap pass run under `orientation: mixed` writes no JSON and does not appear in MultiQC. |
| `phylogeny/phylogeny_qc.json` | Only when `phylogeny.enabled: true`. Tip accounting, alignment shape, and the long-branch report — or a record of why no tree was built. Written on every enabled run, including a skipped one. Nothing reads it back; it is for you. See [Phylogeny](../amplicon/phylogeny.md#the-long-branch-report). |
| `dada2/{sample}.filter_stats.json` | Per-sample `filterAndTrim` result. |
| `dada2/stripped_read_R{1,2}_qual_plot.png` and `.pdf` | Aggregate quality profile of the primer-trimmed reads — the picture behind the `truncLen` choice. |
| `dada2/filtered_read_R{1,2}_error_plot.png` and `.pdf` | DADA2 error-model fit. A bad fit here explains a bad ASV set downstream. |
| `dada2/asv_length_stats.json` and `asv_length_hist.png` | ASV length distributions, described below. |

**`read_tracking.txt`** is one row per sample, one column per stage, in pipeline
order:

```text
raw → [nophix] → stripped → filtered → denoised → merged → non_chimeric
    → [post_extraction] → post_length_filter → post_taxonomy_filter
```

`nophix` appears only when PhiX removal ran, `post_extraction` only when the
extraction step ran. The last column counts reads still present after the
contaminant filter, so the drop between `post_length_filter` and
`post_taxonomy_filter` is exactly what the taxon filter removed: the keep/discard
lists, the best-hit check and, for 16S, the organelle filter.
`06.taxonomy/asv_table_unfiltered.txt` shows which of them removed each ASV. A sample
that collapses at one particular column points straight at the step to inspect.

**`asv_length_stats.json`** and the matching histogram record ASV lengths at
three points: `pre_extraction_asvs` (straight out of DADA2), `filter_source_asvs`
(what the length filter was actually given — the extracted sequences when
extraction ran, otherwise the same as pre-extraction), and `kept_asvs` (what
survived the length window). The JSON also records the window itself
(`filter_window`, with `min`, `max` and where those came from) and
`n_dropped_by_length_filter`. This is most informative for ITS, where ITSx trims
a lot; for 16S the first two distributions usually look the same.

### `multiqc/`

`multiqc_report.html` collects the Falco reports for all three read stages, the
Cutadapt JSONs and the Bowtie2 PhiX logs into one page. `multiqc_data/` holds the
parsed numbers behind it. Falco writes FastQC-format output, so MultiQC files it
under a "FastQC" heading. MultiQC renames the three Falco stages to
`QC | raw reads`, `QC | PhiX-filtered reads` and `QC | primer-trimmed reads`; the
renaming rules are in the shipped MultiQC config.

### `logs/`, `benchmarks/` and `_aux/`

`logs/` holds the captured stdout/stderr of every rule, laid out as
`logs/<rule>/<sample>.log` for per-sample rules and `logs/<rule>.log` for the
rest. It is the first place to look when a rule fails.

`benchmarks/` holds Snakemake's runtime and memory measurements, but only for the
two rules that declare them: `rm_phix` and `trim_primers`.

`_aux/primers/` holds the reverse-complemented primer FASTAs that Cutadapt needs
for 3' read-through trimming. They are generated from the primer files given in
`amplicon.primers`, and are of no interest after the run.

### Outside `out_dir`: the probe cache

When `amplicon.expected_length` is `auto`, the amplicon-length probe writes its
result to `refdb/cache/` inside the MetaFlux checkout, not to `out_dir`. This one
path is fixed: unlike the reference databases it does not follow
`references.refdb_root`.

```text
refdb/cache/probe_{type}_{reference-tag}_{primer-hash}.json
refdb/cache/probe_{type}_{reference-tag}_{primer-hash}.amplicons.fa.gz
```

The filename carries the marker, the reference it was probed against (for 16S also
the lineages left out of the statistics; for rpoB `_excl-second-site`, because products
that hold a second forward-primer site are left out), and a hash of both primer
sequences, so the probe is skipped on every later run with the same marker and primers —
including runs writing to a different `out_dir`. Changing a primer changes the hash and
the probe runs again.

---

## Shotgun mode

### Layout

```text
out_dir/
├── 01.preprocessing/    # QC reports and decontamination stats (reads are temporary)
├── 02.classification/
│   ├── kraken2/         # Kraken2 reports, per-read output, classified reads, k-mer
│   │                    # evidence (profiler kraken2 | both)
│   └── sylph/           # sylph sketches and sylph-tax files (profiler sylph | both)
├── 03.abundance/        # Bracken reports + otu_table.tsv / otu_table.biom; sylph tables
├── 04.extracted_reads/  # per-taxon reads (only when extract_taxa is set)
├── stats/               # read_tracking.txt; sylph/ summary (profiler sylph | both)
├── multiqc/             # multiqc_report.html + multiqc_data/
├── samples.tsv          # sample → source → FASTQ manifest
└── logs/                # per-rule logs
```

The result of a shotgun run is `03.abundance/otu_table.tsv` (taxon × sample
abundances) with its BIOM twin `otu_table.biom`, the per-sample Bracken reports
next to them, and `stats/read_tracking.txt`. With `shotgun.profiler: sylph` the
Kraken2 and Bracken files are not made and the result is
`03.abundance/sylph_table.tsv` instead; with `both` there are both (see
[below](#with-shotgunprofiler-sylph-or-both)).

### `01.preprocessing/`

The cleaned FASTQs here are **temporary** — Snakemake removes a file marked that
way as soon as every rule that needed it has finished: `_dephix` goes once host
removal (or fastp) has read it, `_dehost` once fastp has, and `_trim` once Kraken2
has, or `sylph_sketch`, or both, depending on `shotgun.profiler`. What persists is the
QC and the counts:

| File | Written by | Contents |
|---|---|---|
| `{sample}_fastp.html` / `_fastp.json` | `trim_adapters` | fastp's adapter and quality trimming report. The JSON is also read by MultiQC, by the read-tracking step, and by Bracken (to pick the k-mer distribution matching the post-trim read length). |
| `{sample}_dephix_stats.txt` | `decontam_phix` | BBDuk PhiX-removal counts. Only when `shotgun.decontamination.remove_phix: true`. |
| `{sample}_dehost_stats.txt` | `decontam_host` | Surviving non-host pairs after BBMap. Only when `shotgun.decontamination.host_genomes` is non-empty. |
| `refs/bbmap_index/` | `build_host_index` | The BBMap index of the host reference, built once and reused for every sample. Host removal only. |

The temporary FASTQs are named for the stage they came off:
`{sample}_dephix_R{1,2}.fastq.gz` (BBDuk, PhiX) → `{sample}_dehost_R{1,2}.fastq.gz`
(BBMap, host) → `{sample}_trim_R{1,2}.fastq.gz` (fastp). Each stage is optional;
the chain simply skips what is switched off. The concatenated host reference
`refs/host_combined.fna.gz` is temporary too — it is deleted once the BBMap index
has been built from it. See [decontamination](../shotgun/decontamination.md).

### `02.classification/`

One subfolder per profiler: `kraken2/` when Kraken2 runs (`shotgun.profiler: kraken2`
or `both`) and `sylph/` when sylph runs (`sylph` or `both`; see
[below](#with-shotgunprofiler-sylph-or-both)). `kraken2/` holds:

| File | Contents |
|---|---|
| `{sample}_report.txt` | Kraken2's tree report. |
| `{sample}_output.txt` | Kraken2's per-read assignment — one line per read pair. Large, and only useful for read extraction and troubleshooting. |
| `{sample}_R1.fastq.gz`, `_R2` | The reads Kraken2 managed to classify. Kept, because per-taxon extraction reads them. |
| `{sample}_species_evidence.tsv`, `{sample}_report_gated.txt` | The k-mer evidence of every taxon at the Bracken rank, and the Kraken2 report without the taxa that fail it. Written on every run; Bracken reads the gated report only when `shotgun.kmer_evidence.enabled` is `true`. See [species-level k-mer evidence](../shotgun/kmer-evidence.md). |

!!! warning "The Kraken2 report has extra columns"

    MetaFlux runs Kraken2 with `--report-minimizer-data`, which inserts two
    minimizer-count columns into the report. The result has 8 columns rather than
    the standard 6, with rank, taxid and name pushed to the right. Anything that
    parses these reports by column position needs to know that — including the
    workflow's own read-tracking script, which reads rank from column 6.

### `03.abundance/`

Bracken re-estimates abundance from the Kraken2 report at the rank set in
`shotgun.bracken.tax_lev`, correcting the read counts that Kraken2 leaves stranded
at higher ranks.

| File | Contents |
|---|---|
| `{sample}_report.txt` | The Bracken-corrected report, same layout as a Kraken2 report. Feeds the OTU table and MultiQC. Since v2.5.0 MetaFlux uses Bracken 3.1p1, bioconda's build of 3.1 with a fix for the report: plain 3.1 writes `0` in the third column (reads assigned directly to the taxon) on almost every row when `tax_lev` is not `S`. The second column, the clade count, is correct in both, and it is the only count MetaFlux uses. |
| `{sample}_output.txt` | Bracken's own per-taxon table for that sample. |
| `otu_table.biom` | Cross-sample abundance table in BIOM (HDF5) format. |
| `otu_table.tsv` | The same table as tab-separated text. |

Rule `finalize_otu_table` builds both OTU tables straight from the per-sample Bracken
reports, from the same rows and numbers, so they always agree. Each row is one taxon at
`tax_lev`, each column one sample, and each value is the clade count of that taxon's row
in the sample's Bracken report (0 where the sample does not have it). Each row carries
a lineage of eight slots taken from the reports' own rank codes (see
[below](#taxonomy-strings-in-the-final-tables)). Rows come in the order the taxa are
first met, reading the samples in the order of the table's columns, which are sorted by
sample name (as in `samples.tsv`). The rule then applies the optional
`shotgun.taxonomy_filter` keep/discard lists. In the BIOM file the sample metadata
field `Id` holds the plain sample name.

### With `shotgun.profiler: sylph` or `both`

sylph writes its own files, next to the Kraken2 ones under `both` and instead of them
under `sylph`. The two sets are never merged.

| File | Contents |
|---|---|
| `03.abundance/sylph_table.tsv` | Read pairs per genome and sample, estimated from genome coverage, plus a last row `unassigned`. Each column adds up to `trimmed` in `stats/read_tracking.txt`, except by a few pairs per million in a sample sylph explains to 100 % (see [Reading the numbers](../shotgun/sylph.md#reading-the-numbers)). |
| `03.abundance/sylph_table.biom` | The same table in BIOM (HDF5) format. |
| `03.abundance/sylph_table.genome_copies.tsv` | Each genome's share of genome copies (%) among the detected genomes; same rows, no `unassigned`. |
| `03.abundance/sylph_genomes.tsv` | One row per genome: species name, whether it is a placeholder, the FASTA header of the reference genome, the placeholder-genus flag, best coverage and ANI. |
| `03.abundance/sylph_profile.tsv` | sylph's own output, unchanged: one row per sample and genome. |
| `stats/sylph/sylph_summary.tsv` | One row per sample: read pairs after QC, the share the detected genomes explain, genomes found, flagged genomes. `sylph_summary_mqc.tsv` next to it is the copy MultiQC reads. |
| `02.classification/sylph/` | The sketches (`sketches/{sample}.paired.sylsp`, kept), sylph-tax's `taxprof/{sample}.sylphmpa` files, and sylph's own log lines for each sample. |

How to read these numbers is explained on the [sylph](../shotgun/sylph.md) page.

### `04.extracted_reads/`

Only produced when `shotgun.extract_taxa` lists at least one taxon name. Files are
`{sample}_{taxon}_R1.fastq.gz` and `_R2`, where `{taxon}` is the configured name
lowercased with spaces turned into underscores — `Escherichia coli` becomes
`escherichia_coli`. Reads assigned to that taxon *and everything below it* are
pulled out of the classified reads. A taxon absent from a sample yields an empty
pair of files rather than a failure. See
[per-taxon read extraction](../shotgun/read-extraction.md).

### `stats/read_tracking.txt`, `samples.tsv` and `multiqc/`

`read_tracking.txt` is one row per sample, one column per stage:

```text
sample → raw → [nophix] → [no_host] → trimmed → classified → unclassified
```

`nophix` appears only when PhiX removal ran, `no_host` only when host removal ran.
**Every number is read pairs**, not reads: BBDuk and fastp count the two mates
separately and the script halves them, while Kraken2 already counts a pair once.
`classified` and `unclassified` come from the Kraken2 report's root and
unclassified rows. Under `shotgun.profiler: sylph`, Kraken2 does not run and the table
ends at `trimmed`.

`samples.tsv` is a small manifest with the columns `sample`, `source`,
`bioproject`, `r1`, `r2` — handy for a downstream R or Python script that should
not have to re-parse the config or re-glob the FASTQ directory. `source` is
`local` or `bioproject`.

`multiqc/multiqc_report.html` collects fastp, Kraken2, Bracken and — when PhiX
removal ran — the BBDuk statistics. BBMap host removal contributes nothing to the
report; its counts are in `stats/read_tracking.txt` and its mapping rate in
`logs/decontam_host/{sample}.log`. When sylph runs, the report also has a sylph table
section, built from `stats/sylph/sylph_summary_mqc.tsv`; under `profiler: sylph` it
has no Kraken2 or Bracken sections.

---

## Taxonomy strings in the final tables

Both modes write taxonomy the same way: ranks joined by `;` with no space, each rank
carrying its prefix. Which ranks a table has follows one rule, the same for every marker
and every profiler:

- **`d__` is always the domain slot and `k__` always the kingdom slot.** A table has a
  kingdom only where its reference defines one; MetaFlux never makes one up.
- **Kraken2 tables follow Kraken2's report codes**: always eight slots,
  `d__ k__ p__ c__ o__ f__ g__ s__`, on any database.
- **Amplicon tables and sylph tables follow their reference files.**

| Output | Slots | Count |
|---|---|--:|
| 16S, SILVA 144 | `d__ k__ p__ c__ o__ f__ g__ s__` (`s__` from `addSpecies`, so under `method: rdp` only) | 8 |
| 16S, SILVA 138.2 | `d__ p__ c__ o__ f__ g__ s__` | 7 |
| ITS, UNITE | `k__ p__ c__ o__ f__ g__ s__` (UNITE has no domain) | 7 |
| 18S, PR2, `method: rdp` | `d__ sg__ dv__ sbd__ c__ o__ f__ g__ s__` | 9 |
| 18S, PR2, `method: sintax` | `d__ sg__ dv__ c__ o__ f__ g__ s__` (`dv__` holds Division and Subdivision joined by `-`) | 8 |
| gyrB, DD7RZ8 | the gene tag (`gyrB` or `other`, no prefix), then `p__ c__ o__ f__ g__ s__` | 7 |
| rpoB, FROGS | `d__ p__ c__ o__ f__ g__ s__` | 7 |
| Kraken2 (`otu_table`), any database | `d__ k__ p__ c__ o__ f__ g__ s__` (`k__` empty on a GTDB database, which has no kingdom rank, and for taxa without a kingdom, most protists among them) | 8 |
| sylph (`sylph_table`) | `d__ p__ c__ o__ f__ g__ s__` (GTDB has no kingdom) | 7 |

In the amplicon tables, ranks the classifier could not resolve are simply left out, so a
string can end early (`d__Bacteria;k__Pseudomonadati;p__Pseudomonadota`). The shotgun
tables always write every slot and leave empty any slot they have no name for, and every
slot below `shotgun.bracken.tax_lev` (`…;g__Bacteroides;s__` in a genus-level table).

**Values in `d__` that are not domains.** Three kinds:

- **A virus realm, in the Kraken2 table.** A virus has no domain. ICTV, which names virus
  taxa, made the realm the top rank of virus taxonomy, in place of the domain, and NCBI
  keeps domains and realms in one field of its taxonomy, so the realm sits where a
  cellular organism's domain sits, with the viral kingdom in `k__`:
  `d__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;c__Caudoviricetes;o__;f__Peduoviridae;g__Peduovirus;s__Peduovirus P37`.
  A virus that NCBI places in no realm keeps an empty `d__`, and so does anything else
  without a domain, plasmids for example. Realm names follow the NCBI taxonomy the
  database was built with. `d__Viruses` never appears, also on a database built when
  NCBI still ranked Viruses a superkingdom.
- **`d__Unassigned`**, on the `unassigned` row of the sylph table.
- **PR2's organelle tags**, in 18S tables: `d__Eukaryota:plas` (plastid), `:mito`,
  `:nucl` (nucleomorph), `:apic` (apicoplast) and `:chro` (chromatophore).

**`k__` in ITS tables** holds whatever UNITE files at its top rank: kingdoms
(`k__Fungi`, `k__Viridiplantae`, `k__Metazoa`), protist groups (`k__Alveolata`,
`k__Stramenopila`, `k__Rhizaria`, …) and placeholders (`k__Protista`, `k__unidentified`,
`k__Eukaryota_kgd_Incertae_sedis`).

**With `shotgun.profiler: both`** the same organism reads differently in the two tables:
`otu_table.tsv` has eight slots and `sylph_table.tsv` seven, and the two databases can
name a taxon differently. From one run:

```text
otu_table.tsv    d__Archaea;k__Methanobacteriati;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter smithii
sylph_table.tsv  d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanocatella;s__Methanocatella smithii
```

No MetaFlux step joins the two tables. Compare them by rank name, not by column
position.

**Other tools** may write the same lineage another way: kraken-biom and MetaPhlAn, for
example, put the domain under `k__` (`k__Bacteria`). Filter tokens and scripts written
for their tables need the `d__` prefix in MetaFlux's.

With SILVA 144 and `silva_other: verbatim` (the default) a 16S rank name can end in
`--other` (`g__Bradyrhizobium--other`); the species slot is still written with the plain
genus (`g__Lactobacillus--other;s__Lactobacillus iners`). The kingdom never carries the
suffix. See
[16S rRNA](../amplicon/markers/16S.md#the-other-suffix-amplicontaxonomysilva_other).

The exact form of the species slot follows the reference database, not the
pipeline:

| Source | Species slot | Example |
|---|---|---|
| 16S, SILVA via `rdp` (`addSpecies`) | `s__Genus species` (a binomial, space-separated) | `s__Bacillus subtilis` |
| 16S via `sintax` | none: SILVA's SINTAX file stops at genus | — |
| ITS, UNITE, both methods | UNITE's own `s__Genus_species` (underscore) | `s__Agaricus_bisporus` |
| 18S, PR2 | PR2's own `s__Genus_species` | `s__Unruhdinium_kevei` |
| gyrB (DD7RZ8), rpoB (FROGS) | the reference's own `s__Genus_species` | `s__Staphylococcus_simiae` |
| Shotgun `otu_table` | the whole species name, as the Kraken2 database writes it | `s__Escherichia coli`, `s__Peduovirus P37`, `s__uncultured Martelella sp.` |
| Shotgun `sylph_table` | the GTDB species name exactly as GTDB writes it, placeholders included | `s__Escherichia coli_E`, `s__JABGUH01 sp018646675` |

In the amplicon tables the binomial is only assembled when the reference hands back a
bare, lowercase epithet (`subtilis`, from `addSpecies`; or, under `rdp`, the epithet
that DADA2 cuts out of UNITE's species, see [ITS](../amplicon/markers/its.md)). A
reference that already returns a complete name — PR2's `Unruhdinium_kevei`, for
instance — is passed through untouched, so the genus is never glued on twice.

This happens in the taxonomy string only. In `taxon_seq_table.txt`, the `Species` column
of 16S under `rdp` holds `addSpecies`' bare epithet (`subtilis`), and the binomial is
only in the string; [Reading the tables in R](#reading-the-tables-in-r) shows how to
build it from the table. For every other marker and method the `Species` cell holds the
same name as the `s__` slot.

In the sylph table, a lineage table that writes the domain as `k__` (FungiRefSeq's
`k__Eukaryota`) is rewritten to `d__`, so every row uses the same seven prefixes. See
[sylph](../shotgun/sylph.md#what-you-get).

## Reading the tables in R

The tables follow the ranks of their reference (above), so one rank can sit in a
different column from one marker, release or profiler to the next: `tax[, 2]` is the
kingdom on SILVA 144 and the phylum on SILVA 138.2. Pick ranks by name,
`tax[, "Phylum"]`, which is the same column in every table that has a phylum.

**Amplicon tables.** `taxon_seq_table.txt` has one column per rank, already named after
it, so it is the taxonomy table to hand to phyloseq, mia or any other package; there is
no need to parse the `taxonomy` strings. A rank column with no value in any row comes
into R as `logical` rather than as text, because `read.table` sees only `NA` in it;
`Species` on 16S under `sintax` is always such a column. `colClasses = "character"`
reads every column as text:

```r
tax <- read.table("06.taxonomy/taxon_seq_table.txt", sep = "\t", header = TRUE,
                  row.names = 1, colClasses = "character")
tax[, "Phylum"]
```

On 16S under `rdp`, `Species` holds only the epithet (`subtilis`). This line writes the
binomial of the taxonomy string into it, after dropping a trailing `--other` from the
genus (SILVA 144 with `silva_other: verbatim`, the default). Run it on 16S `rdp` tables
only: in the other tables the cell already holds the whole name.

```r
tax[, "Species"] <- ifelse(is.na(tax[, "Species"]), NA,
                           paste(sub("--other$", "", tax[, "Genus"]), tax[, "Species"]))
```

Do not parse the taxonomy strings with phyloseq's `parse_taxonomy_greengenes`. It knows
only the one-letter prefixes `k__` to `s__`, so it names the `d__` column `Rank1`, and it
misreads the longer 18S prefixes: `sg__TSAR` comes out as a column named `Species`
holding `sTSAR`.

**Shotgun BIOM files** (`otu_table.biom`, `sylph_table.biom`):

- phyloseq's `import_biom()` with `parse_taxonomy_greengenes` names the first column
  `Rank1`, for the same reason; rename it:

    ```r
    library(phyloseq)
    ps <- import_biom("03.abundance/otu_table.biom",
                      parseFunction = parse_taxonomy_greengenes)
    colnames(tax_table(ps))[1] <- "Domain"
    ```

- mia names the columns after the prefixes itself (`domain`, `kingdom`, `phylum`, …)
  when asked to:

    ```r
    tse <- mia::importBIOM("03.abundance/otu_table.biom",
                           rank.from.prefix = TRUE, prefix.rm = TRUE)
    ```

- ampvis2 accepts a `Domain` or a `Kingdom` column, but stops with "Cannot have both
  Domain and Kingdom columns at the same time in taxonomy." on a table that has both, as
  every Kraken2 table and every 16S table on SILVA 144 does. Drop one of the two columns
  before loading.

---

## See also

- [Configuration](configuration.md) — the keys referred to above.
- [Running MetaFlux](running.md) — what a re-run rebuilds, and what it leaves alone.
- [Troubleshooting](../troubleshooting.md) — what an empty or collapsing table usually means.
