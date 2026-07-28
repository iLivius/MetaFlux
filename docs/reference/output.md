# Output files

Everything a run produces goes under the directory named in `output.out_dir`.
Two things are deliberately written elsewhere: reference databases, which are
cached under `references.refdb_root` (`refdb/` by default) so that every run
reuses the same downloads, and — in shotgun BioProject mode — the FASTQs fetched
from SRA, which land in `input.fastq_dir`.

Inside `out_dir` the directories are numbered in the order the stages run, so a
plain `ls` reads like the pipeline itself. Both modes also write `logs/` (one log
per rule, one per sample where the rule runs per sample), `stats/read_tracking.txt`
and `multiqc/multiqc_report.html`.

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
├── 1.reads/          # symlinks to the input FASTQs (visibility only)
├── 2.no_phix/        # PhiX-filtered reads      (absent when remove_phix: false)
├── 3.stripped/       # primer-trimmed reads
├── 4.filtered/       # DADA2 filterAndTrim output
├── 5.dada2/          # ASV sequences and count tables at each stage
├── 6.taxonomy/       # asv_table.txt, asv_table_seqs.txt, taxon_seq_table.txt
│                     # (+ _sintax_raw.tsv on the sintax path)
├── stats/            # read counts, falco/, cutadapt/, trunclen.json, dada2/
├── multiqc/          # multiqc_report.html + multiqc_data/
├── logs/             # per-rule logs
├── benchmarks/       # runtime/memory for rm_phix and trim_primers
└── _aux/primers/     # reverse-complemented primers used by Cutadapt
```

Three files carry the result of an amplicon run: `6.taxonomy/asv_table.txt`
(ASV × sample counts plus a taxonomy string), `5.dada2/seqs_lenfilt.fasta` (the
sequences those ASV IDs refer to), and `stats/read_tracking.txt` (how many reads
survived each stage, per sample). Everything else is intermediate data or QC.

### Reads: `1.reads/` → `4.filtered/`

| Directory | Files | Rule | What happened |
|---|---|---|---|
| `1.reads/` | `{sample}_R1.fastq.gz`, `_R2` | `link_reads` | Symlinks pointing back at `input.fastq_dir`. Nothing reads from here — the downstream rules open the originals. The links exist so the run directory shows which files went in. |
| `2.no_phix/` | `{sample}_R1.fastq.gz`, `_R2` | `rm_phix` | Bowtie2 against the fetched PhiX reference; the pairs that do *not* align are kept. Present only when `amplicon.decontamination.remove_phix` is `true` (the default). |
| `3.stripped/` | `{sample}_R1.fastq.gz`, `_R2` | `trim_primers` | Cutadapt removes the primers at the 5' end and any read-through of the opposite primer at the 3' end. With `primers.orientation: mixed` a second pass recovers reads that arrived the other way round and reverse-complements them back. |
| `4.filtered/` | `{sample}_R1_filt.fastq.gz`, `_R2_filt` | `dada_filter` | DADA2 `filterAndTrim` with the `truncLen` decided in `stats/trunclen.json`. This is the last per-sample FASTQ; from here DADA2 works on all samples at once. |

### Inside `5.dada2/`

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

The last populated column is the final result. `assign_taxonomy` reads exactly two
of these files: the length-filtered FASTA `seqs_lenfilt.fasta` and the
`ASV_#`-keyed `seqtab_lenfilt_head_names.txt`. (The read-tracking step also opens
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

| Marker | Extractor | Files in `5.dada2/` |
|---|---|---|
| 16S | Metaxa2 | `metaxa2_extraction.results.txt` (Metaxa2's own summary) plus the raw `seqs.*` working files under `metaxa2/` |
| ITS | ITSx | `itsx_extraction.summary.txt`, `itsx_collapse_map.tsv`, plus the raw `seqs.*` working files under `_itsx_tmp/` |

`itsx_collapse_map.tsv` is worth knowing about. Two ASVs can differ only in the
flanking 5.8S/LSU sequence; once ITSx trims those flanks away the remaining ITS1
or ITS2 region is identical, so the ASVs are merged into one and their counts are
added. The map records which original `ASV_#` went into which representative, so
that merge is auditable rather than silent.

### `6.taxonomy/`

Three tables, all written by `assign_taxonomy`, all built from the same filtered
set of ASVs. They differ only in what labels the rows.

With `taxonomy.method: sintax` — the shipped default — the directory also holds
`_sintax_raw.tsv`, the raw VSEARCH `--tabbedout` result the parser was fed. It is
kept for troubleshooting and no rule reads it.

| File | Rows | Columns |
|---|---|---|
| `asv_table.txt` | `ASV_#` | one per sample, then `taxonomy` |
| `asv_table_seqs.txt` | the ASV sequence | one per sample, then `taxonomy` |
| `taxon_seq_table.txt` | `ASV_#` | one per taxonomic rank, then `sequence` |

The rank columns in `taxon_seq_table.txt` come from the marker pack, so they
match the reference's actual depth: `Kingdom … Species` for 16S, ITS and rpoB,
`Gene, Phylum … Species` for gyrB, and PR2's nine ranks
(`Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species`)
for 18S.

!!! note "File format"

    These follow R's `write.table(col.names = NA)` conventions, on both the rdp
    and the sintax path: tab-separated, the header line starts with an empty cell
    for the row-name column, and text values (row names, taxonomy strings, rank
    names) are wrapped in double quotes. `read.table(..., sep = "\t", header = TRUE,
    row.names = 1)` in R and `pd.read_csv(..., sep = "\t", index_col = 0)` in
    pandas both read them without extra arguments.

The contaminant filter (`amplicon.taxonomy.filter.keep` / `.discard`) has already
been applied when these tables are written, so the ASVs it removed are absent
here — see [keeping and discarding taxa](../amplicon/taxon-filter.md).

### `stats/`

| Path | Contents |
|---|---|
| `read_tracking.txt` | The per-stage summary, described below. |
| `1.raw_reads.counts` | R1 read count per sample, `{sample} : N`, with a trailing total. |
| `2.nophix_reads.counts` | Same, after PhiX removal. Only when `remove_phix: true`. |
| `3.stripped_reads.counts` | Same, after primer trimming. |
| `trunclen.json` | The truncation decision: what `truncLen` was chosen for R1 and R2 and how. See [amplicon length and truncation](../amplicon/length-and-truncation.md). |
| `falco/{sample}_R{1,2}_{stage}/` | Falco QC per sample, per read direction, per stage (`raw`, `nophix`, `stripped`), each holding `fastqc_data.txt`, `fastqc_report.html` and `summary.txt`. |
| `cutadapt/{sample}.passA_5prime.cutadapt.json` and `.passA_3prime.` | Cutadapt's own report for the 5' and 3' trimming passes. These feed MultiQC. The extra swap pass run under `orientation: mixed` writes no JSON and does not appear in MultiQC. |
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
`post_taxonomy_filter` is exactly what the keep/discard lists removed. A sample
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

The filename carries the marker, the reference it was probed against, and a hash
of both primer sequences, so the probe is skipped on every later run with the same
marker and primers — including runs writing to a different `out_dir`. Changing a
primer changes the hash and the probe runs again.

---

## Shotgun mode

### Layout

```text
out_dir/
├── 01.preprocessing/    # QC reports and decontamination stats (reads are temporary)
├── 02.classification/   # Kraken2 reports, per-read output, classified reads
├── 03.abundance/        # Bracken reports + otu_table.tsv / otu_table.biom
├── 04.extracted_reads/  # per-taxon reads (only when extract_taxa is set)
├── stats/               # read_tracking.txt
├── multiqc/             # multiqc_report.html + multiqc_data/
├── samples.tsv          # sample → source → FASTQ manifest
└── logs/                # per-rule logs
```

The result of a shotgun run is `03.abundance/otu_table.tsv` (taxon × sample
abundances) with its BIOM twin `otu_table.biom`, the per-sample Bracken reports
next to them, and `stats/read_tracking.txt`.

### `01.preprocessing/`

The cleaned FASTQs here are **temporary** — Snakemake removes a file marked that
way as soon as every rule that needed it has finished: `_dephix` goes once host
removal (or fastp) has read it, `_dehost` once fastp has, and `_trim` once Kraken2
has. What persists is the QC and the counts:

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

| File | Contents |
|---|---|
| `{sample}_report.txt` | Kraken2's tree report. |
| `{sample}_output.txt` | Kraken2's per-read assignment — one line per read pair. Large, and only useful for read extraction and troubleshooting. |
| `{sample}_R1.fastq.gz`, `_R2` | The reads Kraken2 managed to classify. Kept, because per-taxon extraction reads them. |

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
| `{sample}_report.txt` | The Bracken-corrected report, same layout as a Kraken2 report. Feeds the OTU table, per-taxon extraction, and MultiQC. |
| `{sample}_output.txt` | Bracken's own per-taxon table for that sample. |
| `otu_table.biom` | Cross-sample abundance table in BIOM (HDF5) format. |
| `otu_table.tsv` | The same table as tab-separated text. |

The two OTU tables are rendered from the same filtered table by
`finalize_otu_table`, so they always agree. That rule also tidies up what
`kraken-biom` produces: it strips the `_report` suffix from sample names,
normalises the taxonomy strings (see below), and applies the optional
`shotgun.taxonomy_filter` keep/discard lists. The untidied intermediate
`otu_table.raw.biom` is temporary and does not survive the run.

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
unclassified rows.

`samples.tsv` is a small manifest with the columns `sample`, `source`,
`bioproject`, `r1`, `r2` — handy for a downstream R or Python script that should
not have to re-parse the config or re-glob the FASTQ directory. `source` is
`local` or `bioproject`.

`multiqc/multiqc_report.html` collects fastp, Kraken2, Bracken and — when PhiX
removal ran — the BBDuk statistics. BBMap host removal contributes nothing to the
report; its counts are in `stats/read_tracking.txt` and its mapping rate in
`logs/decontam_host/{sample}.log`.

---

## Taxonomy strings in the final tables

Both modes write taxonomy the same way structurally: ranks joined by `;` with no
space, each rank carrying its prefix (`k__`, `p__`, `c__`, `o__`, `f__`, `g__`,
`s__`; 18S uses the PR2 rank set, `d__`, `sg__`, `dv__` and so on). Ranks the
classifier could not resolve are simply left out, so a string can end early.

The exact form of the species slot follows the reference database, not the
pipeline:

| Source | Species slot | Example |
|---|---|---|
| SILVA / 16S via `rdp` | `s__Genus species` (a binomial, space-separated) | `s__Bacillus subtilis` |
| UNITE / ITS | UNITE's own `s__Genus_species` (underscore), emitted verbatim | `s__Genus_species` |
| Shotgun `otu_table` | the same binomial reconstruction as SILVA, where the reference supplies a bare epithet | `s__Bacillus subtilis` |
| `sintax` path | whatever label the SINTAX database carries | — |

The binomial is only assembled when the reference hands back a bare, lowercase
epithet (`subtilis`). A reference that already returns a complete name — PR2's
`Unruhdinium_kevei`, for instance — is passed through untouched, so the genus is
never glued on twice.

---

## See also

- [Configuration](configuration.md) — the keys referred to above.
- [Running MetaFlux](running.md) — what a re-run rebuilds, and what it leaves alone.
- [Troubleshooting](../troubleshooting.md) — what an empty or collapsing table usually means.
