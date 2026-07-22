# MetaFlux

```
______  ___    _____       _______________
___   |/  /______  /______ ___  ____/__  /___  _____  __
__  /|_/ /_  _ \  __/  __ `/_  /_   __  /_  / / /_  |/_/
_  /  / / /  __/ /_ / /_/ /_  __/   _  / / /_/ /__>  <
/_/  /_/  \___/\__/ \__,_/ /_/      /_/  \__,_/ /_/|_|

v1.0.0
```

**Unified Illumina 16S/ITS amplicon and shotgun taxonomic profiling workflow.**

![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A59.0-brightgreen.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![DOI](https://img.shields.io/badge/DOI-pending%20release-lightgrey.svg)

## Synopsis

MetaFlux is a [Snakemake](https://snakemake.github.io/) workflow that takes raw
paired-end Illumina reads through to a taxonomically annotated abundance table. It
runs in one of two modes, selected by a single config key:

- **`amplicon`** — a [DADA2](https://benjjneb.github.io/dada2/)-based pipeline for
  16S rRNA and ITS metabarcoding: PhiX removal, primer trimming, exact ASV
  inference, marker-region extraction, length filtering, and taxonomic assignment
  against SILVA (16S) or UNITE (ITS).
- **`shotgun`** — a [Kraken2](https://github.com/DerrickWood/kraken2) +
  [Bracken](https://github.com/jenniferlu717/Bracken) pipeline for shotgun
  metagenomic taxonomic profiling: decontamination, quality trimming, read
  classification, abundance re-estimation, and OTU-table construction.

Both modes run from the same command and the same config file, and produce a single
[MultiQC](https://multiqc.info/) report.

`MetaFlux` belongs to the [BioFlux](https://github.com/stars/iLivius/lists/bioflux) family of pipelines.

## Table of Contents

- [Quick Start](#quick-start)
- [Rationale](#rationale)
- [Description](#description)
  - [Amplicon mode](#amplicon-mode)
  - [Shotgun mode](#shotgun-mode)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running MetaFlux](#running-metaflux)
- [Output](#output)
- [Troubleshooting](#troubleshooting)
- [Acknowledgements](#acknowledgements)
- [Citation](#citation)
- [References](#references)
- [License](#license)

## Quick Start

```bash
# 1. Clone
git clone https://github.com/iLivius/MetaFlux.git
cd MetaFlux

# 2. Create the launcher environment (Snakemake + conda support)
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake

# 3. Edit config/config.yaml — set `mode`, input/output paths, and references
#    (amplicon: primer FASTAs; shotgun: a Kraken2 database)

# 4. Dry-run to check the plan, then execute
snakemake --sdm conda --cores 16 --configfile config/config.yaml --dry-run
snakemake --sdm conda --cores 16 --configfile config/config.yaml
```

`--sdm conda` (software-deployment-method) tells Snakemake to create each rule's
own pinned conda environment automatically on first run. The mode and input paths
are set in the config before launching — see [Configuration](#configuration).

[↑ Back to top](#table-of-contents)

## Rationale

Over many years of amplicon-sequencing and metagenomics work, a set of methods and
parameter choices accumulates into an in-house pipeline.

MetaFlux re-implements that methodology as a more modern, modular Snakemake workflow,
adding several capabilities the previous implementation never had:

1. **Data-driven read truncation.** The DADA2 `truncLen` can be chosen
   automatically from the per-base quality profile, with an overlap-aware safeguard
   so forward and reverse reads still merge — instead of being picked by eye.
2. **Agnostic amplicon-length inference.** The expected amplicon length is derived
   from the reference database using the supplied primers (in-silico PCR against SILVA, or
   the UNITE pre-extracted subregions for ITS), and that distribution then drives an
   automatic ASV length filter — no hard-coded numbers.
3. **Orientation-aware primer trimming.** Cutadapt strips primers and, when reads
   arrive in mixed orientation, recovers and re-orients the swapped reads.
4. **Fine-grained resource control.** Every rule's CPU and RAM are set in one place,
   so the same workflow scales from a laptop to an HPC cluster unchanged.
5. **One pipeline, two modes.** Amplicon and shotgun processing share a single
   reproducible entry point and infrastructure.

Each step runs in its own pinned conda environment, the workflow is resumable after
an interruption, and every parameter lives in a single version-controlled config
file — so a run can be reproduced exactly.

[↑ Back to top](#table-of-contents)

## Description

The `mode` key in the config decides which pipeline runs. A few steps are shared by
both modes — read QC and its MultiQC report, read tracking, reference handling, and
the per-rule resource settings — while everything else is specific to the selected
mode.

```
                 ┌──────────────┬───────────────┐
   config: mode  │   amplicon   │    shotgun    │
                 └──────┬───────┴───────┬───────┘
        raw paired-end Illumina reads (R1/R2)
                        │               │
        bowtie2 PhiX*   │               │  BBDuk PhiX*  +  BBMap host*
        Cutadapt primer │               │  fastp trim
        DADA2 ASVs      │               │  Kraken2 classify
        Metaxa2 / ITSx* │               │  Bracken abundance
        length filter   │               │  kraken-biom OTU table
        SILVA / UNITE   │               │  KrakenTools per-taxon reads*
        taxonomy        │               │
                        ▼               ▼
              ASV table + taxonomy   OTU table + reports
                        └───────┬───────┘
                     read tracking + MultiQC

   * optional steps, toggled in config
```

### Amplicon mode

DADA2-based 16S/ITS metabarcoding. Steps (each maps to a numbered output directory):

1. **PhiX removal** *(optional)* — [bowtie2](https://github.com/BenLangmead/bowtie2)
   maps against the PhiX reference; unmapped reads are kept. Toggle with
   `amplicon.decontamination.remove_phix`.
2. **Primer trimming** — [Cutadapt](https://cutadapt.readthedocs.io/) removes 5′
   primers and read-through 3′ adapters. With `orientation: mixed`, a second pass
   recovers reads sequenced in the opposite orientation and re-orients them.
3. **Per-stage QC** — [Falco](https://github.com/smithlabcode/falco) profiles reads
   at the raw, PhiX-filtered, and primer-trimmed stages.
4. **Amplicon-length probe** — in-silico PCR of the supplied primers against SILVA (16S)
   or direct measurement of the UNITE ITS1/ITS2 subregions (ITS), yielding a length
   distribution used downstream.
5. **Truncation-length selection** — `auto` picks `truncLen` from the per-base Q1
   profile under a merge-overlap constraint; `manual` uses fixed values.
6. **ASV inference** — DADA2 `filterAndTrim` → error learning → denoising →
   `mergePairs` → chimera removal.
7. **Marker extraction** *(optional)* — [Metaxa2](https://microbiology.se/software/metaxa2/) (16S)
   or [ITSx](https://microbiology.se/software/itsx/) (ITS) isolates the target region. Toggle
   with `amplicon.extraction.enabled`.
8. **ASV length filter** — keeps ASVs within a window derived from the probe
   distribution (`auto`) or a manual range; the before/after length distributions
   are logged to `stats/dada2/`.
9. **Taxonomy** — DADA2 RDP naive Bayesian classifier *or*
   [VSEARCH](https://github.com/torognes/vsearch) SINTAX (selectable), against SILVA
   (16S) or UNITE (ITS), followed by a configurable include/exclude contaminant filter.

### Shotgun mode

Kraken2 + Bracken taxonomic profiling. Steps:

1. **Input** — local paired FASTQs, or automatic NCBI SRA retrieval from a
   **BioProject** accession (resolves PAIRED runs via e-utils, fetches with
   `prefetch` + `fasterq-dump`).
2. **PhiX removal** *(optional)* — [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
   against its bundled PhiX reference (`shotgun.decontamination.remove_phix`).
3. **Host/contaminant removal** *(optional)* —
   [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) maps against
   one or more user-supplied reference genomes (e.g. a masked human genome);
   unmapped reads are kept (`shotgun.decontamination.host_genomes`).
4. **Adapter/quality trimming** — [fastp](https://github.com/OpenGene/fastp).
5. **Classification** — [Kraken2](https://github.com/DerrickWood/kraken2) against a
   user-supplied database (e.g. PlusPF).
6. **Abundance re-estimation** — [Bracken](https://github.com/jenniferlu717/Bracken).
7. **OTU table** — [kraken-biom](https://github.com/smdabdoub/kraken-biom) writes a
   combined BIOM + TSV table.
8. **Per-taxon read extraction** *(optional)* —
   [KrakenTools](https://github.com/jenniferlu717/KrakenTools) pulls reads assigned
   to taxa of interest.

[↑ Back to top](#table-of-contents)

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/iLivius/MetaFlux.git
cd MetaFlux
```

### 2. Install Snakemake

Only Snakemake (≥ 9) needs to be installed by hand; conda installs every other tool,
one environment per rule, when the workflow runs with `--sdm conda`.

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

### 3. Reference databases

| Database                           | Mode               | Provisioning                                     |
|------------------------------------|--------------------|--------------------------------------------------|
| PhiX                               | both               | **Auto-fetched** (NCBI) and indexed on first run |
| SILVA (toGenus + species + SINTAX) | amplicon (16S)     | **Auto-fetched** from Zenodo on first run        |
| UNITE (general + UCHIME + SINTAX)  | amplicon (ITS)     | **Auto-fetched** from PlutoF on first run        |
| Kraken2 database                   | shotgun            | **User-supplied** (see below)                    |
| Host genome(s)                     | shotgun (optional) | **User-supplied** (path or URL)                  |

The amplicon reference databases download automatically the first time they are
needed and are cached under `refdb/`. **Shotgun mode requires a user-supplied Kraken2
database** — it is large and version-sensitive. Pre-built indexes, with sizes, build
dates, and md5 checksums, are listed in Ben Langmead's [catalog](https://benlangmead.github.io/aws-indexes/k2).
A standard choice is the **PlusPF** index. A single `wget` stream works but is slow on
a ~100 GB file; a multi-connection downloader is faster and resumable.

Pick the current build from the catalog:

```bash
mkdir -p /path/to/db && cd /path/to/db

# substitute the latest dated build from the catalog
aria2c -x16 -s16 https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260226.tar.gz

# verify integrity against the md5 listed in the catalog
md5sum k2_pluspf_20260226.tar.gz

# contains hash.k2d, opts.k2d, taxo.k2d, *.kmer_distrib
tar -xzf k2_pluspf_20260226.tar.gz
```

> If the AWS CLI is installed, `aws s3 cp --no-sign-request <url>` is the fastest option (parallel multipart). Plain `wget -c <url>` also works and resumes a partial download.

> **NOTE:** the full PlusPF index is ~100 GB. With memory-mapping enabled (default),
> Kraken2 needs only a modest private workspace and lets the OS cache the database;
> without it, plan for RAM ≥ the `hash.k2d` size. Smaller capped indexes
> (e.g. `pluspf_16gb`) work too — MetaFlux sizes the Kraken2 memory request from the
> actual index automatically.

For **host removal** (optional), the BBTools authors' masked human reference is a
good choice — masking ribosomal/low-complexity and cross-kingdom-homologous regions
avoids removing real microbial reads:

```bash
# hg19 + viral masking, ~889 MB
wget -O /path/to/db/human_virus_masked.fasta.gz \
  "https://zenodo.org/records/4116107/files/human_virus_masked.fasta.gz"
```

[↑ Back to top](#table-of-contents)

## Configuration

All behaviour is controlled by `config/config.yaml`. Paths may be absolute or
relative; **relative paths resolve against the directory snakemake is invoked from**,
not the config file's location — absolute paths are safest.

Each section is tagged `[shared]`, `[amplicon]`, or `[shotgun]` so it is clear which
mode consumes it.

### Mode selection

```yaml
mode: amplicon          # amplicon | shotgun
```

> **NOTE:** `mode` must match the data. If it is set to one mode while the
> input/output paths point at the other mode's data, MetaFlux runs the wrong
> analysis — it does not always error, because both modes accept `{sample}_R1/_R2`
> filenames. The mode can also be overridden without editing the file:
> `--config mode=shotgun`.

### Input and output `[shared]`

```yaml
input:
  fastq_dir: /path/to/fastq_dir          # local paired FASTQs
  # Shotgun-only: fetch a BioProject from NCBI SRA instead of local files
  bioproject: null                       # e.g. PRJNA603575 (null = use local FASTQs)
  accession_list: null                   # optional subset, e.g. [SRR…, SRR…]

output:
  out_dir: /path/to/output
```

- **Amplicon** filenames must be `{sample}_R1.fastq.gz` / `{sample}_R2.fastq.gz`.
- **Shotgun** accepts `{sample}_R1/_R2` or SRA-style `{sample}_1/_2` (auto-detected).
- When `bioproject` is set, MetaFlux resolves all PAIRED runs for that accession;
  `accession_list` restricts this to specific runs (e.g. to keep only WGS runs from a
  BioProject that also contains amplicon or RNA-Seq runs).

### References `[shared / amplicon / shotgun]`

```yaml
references:
  phix:                                   # [shared] — only amplicon fetches/indexes it
    fasta: refdb/phix/phix.fna
    fetch_url: ftp://ftp.ncbi.nlm.nih.gov/.../GCF_000819615.1_..._genomic.fna.gz
  silva:                                  # [amplicon, 16S] — auto-fetched
    train:   refdb/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz
    species: refdb/silva/silva_v138.2_assignSpecies.fa.gz
    sintax:  refdb/silva/silva_nr99_v138.2_sintax.fa.gz
    fetch_url_train:   https://zenodo.org/...
    fetch_url_species: https://zenodo.org/...
  unite:                                  # [amplicon, ITS] — auto-fetched
    fasta:   refdb/unite/...fasta
    sintax:  refdb/unite/...fasta.gz
    fetch_url: https://...
  kraken_db: /path/to/db/k2_pluspf        # [shotgun] — user-supplied
```

### Amplicon parameters `[amplicon]`

```yaml
amplicon:
  type: 16S                               # 16S | ITS
  its_region: ITS2                        # ITS1 | ITS2 (ignored for 16S)
  seed: 100                               # fixes learnErrors/assignTaxonomy RNG fully; sintax needs threads:1 too — see Troubleshooting

  decontamination:
    remove_phix: true                     # bowtie2 PhiX removal; false = skip

  primers:
    fwd: /path/to/fwd_primer.fasta
    rev: /path/to/rev_primer.fasta
    orientation: fixed                    # fixed | mixed (mixed enables the swap pass)

  expected_length: auto                   # int | [min, max] | "auto"
  probe_length_stat: { 16S: p95 }         # statistic used from the probe distribution
  min_overlap: 12                         # required forward/reverse overlap, in bp

  trunc_len:
    mode: manual                          # auto | manual
    q_threshold: 20                       # Q1 cutoff for auto
    resolve_policy: raise_trunc           # raise_trunc | relax_q | error
    manual_r1: 260                        # used when mode == manual
    manual_r2: 240

  length_filter:
    mode: auto                            # auto | manual
    window_margin: 50                     # bp added around the probe q1/p95

  extraction:
    enabled: true                         # Metaxa2 (16S) / ITSx (ITS)

  taxonomy:
    method: sintax                        # rdp | sintax
    min_boot: 80                          # rdp bootstrap minimum
    sintax_cutoff: 0.8                    # sintax confidence (≈ min_boot 80)
    filter: { enabled: true }             # contaminant include/exclude filter
```

See [Troubleshooting](#troubleshooting) for the `auto` truncation caveat.

### Shotgun parameters `[shotgun]`

```yaml
shotgun:
  decontamination:
    remove_phix: true                     # BBDuk PhiX removal; false = skip
    host_genomes: []                      # list of paths/URLs to FASTA(.gz); [] = skip
    host_min_id: 0.95                     # BBMap minimum identity for host removal
  fastp:
    min_read_length: 100                  # --length_required after trimming
  kraken:
    confidence: 0.15                      # Kraken2 confidence threshold
    hit_groups: 3
    memory_mapping: true                  # mmap the DB (low RAM) vs load into RAM
  bracken:
    tax_lev: S                            # taxonomic level for re-estimation
    threshold: 10
  extract_taxa: []                        # optional taxa to pull reads for
```

`host_genomes` is a YAML **list** — each reference is its own `- ` item (a local
path *or* a URL); `[]` means skip host removal. Multiple entries are concatenated and
indexed once.

### Resources `[shared]`

Per-rule CPU and RAM live in one block, with global fallbacks. These can be tuned to
the target machine or cluster; the values feed straight into Snakemake's cluster
executor.

```yaml
resources:
  threads_default: 4
  mem_mb_default: 2000
  threads:
    kraken2: 16
    decontam_host: 16
    # … one entry per rule
  mem_mb:
    kraken2: 20000
    assign_taxonomy: 16000
    # … one entry per rule
```

[↑ Back to top](#table-of-contents)

## Running MetaFlux

Run from the repository root.

```bash
# Always preview the plan first
snakemake --sdm conda --cores 16 --configfile config/config.yaml --dry-run

# Execute
snakemake --sdm conda --cores 16 --configfile config/config.yaml
```

Useful patterns:

```bash
# Override the mode without editing the config
snakemake --sdm conda --cores 16 --config mode=shotgun

# Limit concurrent jobs (e.g. memory-bound steps)
snakemake --sdm conda --jobs 4 --cores 16

# Resume after an interruption
snakemake --sdm conda --cores 16 --rerun-incomplete

# Cluster execution: supply a profile (e.g. SLURM) — per-rule threads/mem_mb are honoured
snakemake --sdm conda --profile <slurm-profile>
```

> **NOTE (shotgun, BioProject mode):** set `input.bioproject` to an NCBI accession
> and MetaFlux fetches the PAIRED runs into `fastq_dir` automatically. The first
> Kraken2 call reads the whole database from disk (slow when cold, fast once the OS
> has cached it).

[↑ Back to top](#table-of-contents)

## Output

Outputs are written under `out_dir` with numbered, stage-ordered directories.

### Amplicon mode

```
out_dir/
├── 1.reads/              # input symlinks (visibility)
├── 2.no_phix/            # PhiX-filtered reads          (omitted if remove_phix: false)
├── 3.stripped/           # primer-trimmed reads
├── 4.filtered/           # DADA2 filterAndTrim output
├── 5.dada2/              # ASV sequences + tables, extraction, length-filtered seqtab
├── 6.taxonomy/           # asv_table.txt, asv_table_seqs.txt, taxon_seq_table.txt
├── stats/                # read_tracking.txt, falco/, cutadapt/, trunclen.json, dada2/ (error/quality plots + ASV length stats)
├── multiqc/              # multiqc_report.html
├── logs/                 # per-rule logs
└── benchmarks/           # per-rule runtime/memory benchmarks
```

Key files: `6.taxonomy/asv_table.txt` (ASV × sample counts + taxonomy, keyed by
`ASV_#` ID) and `6.taxonomy/asv_table_seqs.txt` (the same table keyed by the ASV
sequence), `5.dada2/seqs.fasta` (ASV sequences), `stats/read_tracking.txt`
(reads surviving each stage, per sample).

**Inside `5.dada2/`.** It looks busy, but most of it is just two artifacts — the **ASV
sequences** and the **ASV count table** — saved at each processing stage. The table
below reads left-to-right in pipeline order; within the count table, `head_names`
columns are `ASV_#` IDs and `head_seqs` columns are the full ASV sequence (same counts,
relabeled).

| Artifact                    | DADA2                   | + Extraction *(optional)*         | + Length filter                 |
|-----------------------------|-------------------------|-----------------------------------|---------------------------------|
| **ASV sequences** (FASTA)   | `seqs.fasta`            | `seqs_extracted.fasta`            | `seqs_lenfilt.fasta`            |
| **Counts** — `ASV_#` keyed  | `seqtab_head_names.txt` | `seqtab_extracted_head_names.txt` | `seqtab_lenfilt_head_names.txt` |
| **Counts** — sequence keyed | `seqtab_head_seqs.txt`  | —                                 | `seqtab_lenfilt_head_seqs.txt`  |
| **Extras**                  | `read.counts`           | Metaxa2 / ITSx results + dir      | —                               |

The last populated column is the final result: `seqtab_lenfilt_head_seqs.txt`
(length-filtered, sequence-keyed) is what the taxonomy step reads. With extraction
disabled, the middle column is absent and the length filter runs on the DADA2 output.

The length filter also writes `stats/dada2/asv_length_stats.json` +
`asv_length_hist.png`: ASV lengths at three stages — pre-extraction, post-extraction,
kept — showing how much Metaxa2/ITSx trimmed. Most useful for ITS; for 16S the first
two stages usually match.

### Shotgun mode

```
out_dir/
├── 01.preprocessing/     # (dephix / dehost stats), fastp JSON/HTML
├── 02.classification/    # Kraken2 {sample}_report.txt / _output.txt + classified reads
├── 03.abundance/         # Bracken reports + otu_table.tsv / otu_table.biom
├── 04.extracted_reads/   # KrakenTools per-taxon reads (if extract_taxa set)
├── stats/                # read_tracking.txt
├── multiqc/              # multiqc_report.html
├── samples.tsv           # sample → source → FASTQ manifest
└── logs/                 # per-rule logs
```

Key files: `03.abundance/otu_table.tsv` (taxon × sample abundances),
`03.abundance/{sample}_report.txt` (Bracken), `stats/read_tracking.txt`
(raw → trimmed → classified/unclassified, per sample).

**Inside `01.preprocessing/`.** The cleaned reads themselves are temporary — Snakemake
deletes the `_dephix` / `_dehost` / `_trim` FASTQs once Kraken2 has read them — so what
persists here is per-sample QC and stats, not bulky sequence files:

- `{sample}_fastp.html` / `_fastp.json` — fastp trimming report (the JSON also feeds MultiQC).
- `{sample}_dephix_stats.txt` — BBDuk PhiX-removal counts *(only if `remove_phix: true`)*.
- `{sample}_dehost_stats.txt` — BBMap host-removal counts *(only if `host_genomes` is set)*.
- `refs/` — the concatenated host reference and its BBMap index, built once and reused *(host removal only)*.

The suffix marks the stage the reads came off: `_dephix` (BBDuk PhiX) → `_dehost`
(BBMap host) → `_trim` (fastp). With decontamination off, only the fastp report remains.

[↑ Back to top](#table-of-contents)

## Troubleshooting

### Amplicon: automatic `truncLen` can drop all reads on very clean runs

In `amplicon.trunc_len.mode: auto`, MetaFlux truncates where the per-base lower-
quartile (Q1) quality first falls below `q_threshold`. If quality never drops (clean
MiSeq runs often stay high to the last cycle), `truncLen` lands at the maximum read
length — but after primer trimming most reads end a little earlier, and DADA2's
`filterAndTrim` **discards any read shorter than `truncLen`**, which can drop almost
everything. The symptom is a failure at the DADA2 step:

```
Error: No reads passed the filter ...
```

**Fix:** set the truncation length by hand. Look at the read lengths after primer
trimming — the `*_R{1,2}_stripped` "Sequence Length Distribution" in MultiQC/Falco —
then switch to `mode: manual` and set `manual_r1`/`manual_r2` a little below the length
most reads still reach, so `filterAndTrim` stops discarding them. Keep
`truncLen_R1 + truncLen_R2 ≥ amplicon_length + min_overlap` so the pairs still overlap.
Re-running picks up from there: the trimming and QC are already done, so only DADA2 and
the steps after it run again, and trying a different value is fast.

### Shotgun: BBDuk PhiX removal fails on variable-length (pre-trimmed) reads

The BBDuk PhiX step (`decontam_phix`) reads the R1 and R2 files in parallel — one read
from each at a time, assuming they stay paired — and expects **uniform-length** reads
(as they come off the sequencer). If the input was adapter/quality-trimmed *before*
deposition so read lengths vary — common for SRA/BioProject data, especially when R1
and R2 of a pair differ — the reads from the two files stop matching up and BBDuk
aborts. This is not corrupt data (the FASTQs are intact and correctly paired) — it is a
BBDuk limitation with non-uniform lengths.

The BBMap host step (`decontam_host`) is **not** affected: it aligns reads, so ragged
lengths are no problem. To work around the PhiX step:

- **Skip it** — `remove_phix: false`. In shotgun data PhiX is a negligible spike-in
  fraction, and Kraken2 will not classify it as anything real.
- **Remove PhiX with BBMap instead** — keep `remove_phix: false` and list a phiX174
  FASTA under `host_genomes`; BBMap drops it at the host step, ragged lengths and all.
  The NCBI URL (the same phiX174 the amplicon path fetches) can go in directly:

  ```yaml
  decontamination:
    remove_phix: false
    host_genomes:
      - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz
  ```

- **Use raw, untrimmed reads** if available — MetaFlux trims with fastp anyway.

### Amplicon: ASV counts or taxonomy calls differ between identical re-runs

Three steps in the amplicon path draw on a random-number generator: DADA2's
`learnErrors(randomize=TRUE)` (which subsamples reads to build the error model),
DADA2's `assignTaxonomy` (rdp method — bootstrap confidence), and VSEARCH's
`--sintax` (sintax method — bootstrap confidence, seeded from a random source by
default). Without a fixed seed, each of these draws a different sample every run,
so a handful of ASVs or taxonomy calls near a decision threshold can flip between
two runs on the exact same input — this is expected DADA2/VSEARCH behaviour, not
a bug or corrupted data.

`amplicon.seed` fixes the RNG for all three, but **only `learnErrors` and
`assignTaxonomy` (rdp) become fully reproducible at any thread count** — verified
by running the amplicon test set twice at 8 threads with the same seed:
`seqs.fasta`, ASV counts, and rdp-path taxonomy came back byte-identical.

**`--sintax` (the sintax method) does not**: VSEARCH runs its bootstrap
confidence step across threads that share one RNG stream, so which thread
consumes which random draw still depends on scheduling, seed or no seed. Verified
the same way — with `amplicon.seed` fixed, 8-threaded `assign_taxonomy` reruns
still shifted a handful of per-rank confidence values (e.g. a genus call at 0.99
in one run, 1.00 in the next), while single-threaded reruns came back identical.
If byte-reproducible sintax output matters more than speed, set
`resources.threads.assign_taxonomy: 1`; otherwise treat sub-percent confidence
drift near `sintax_cutoff` as expected run-to-run noise, same as before the seed
was added.

### Config (YAML) errors

Snakemake reports any config problem with one generic message
(*"Config file is not valid JSON or YAML…"*). This almost always means a YAML
indentation/list-syntax issue or a stray tab — not a problem with the parameter
values. Validate directly to get the exact line:

```bash
python -c "import yaml; yaml.safe_load(open('config/config.yaml'))"
```

A *"expected block end, but found block sequence start"* message points to a
misplaced/over-indented `- ` list item (commonly in `host_genomes`).

[↑ Back to top](#table-of-contents)

## Acknowledgements

Developed at the [AIT Austrian Institute of Technology](https://www.ait.ac.at/).
MetaFlux consolidates and modernises methods refined across many amplicon and
metagenomics collaborations, and is part of the **BioFlux** family of workflows.

## Citation

Publications that use MetaFlux should cite this repository. A formal release with a
Zenodo DOI is forthcoming:

> Antonielli, L. (2026). *MetaFlux: a unified Illumina 16S/ITS amplicon and shotgun
> taxonomic profiling workflow.* Zenodo. DOI: pending release.

The underlying tools should also be cited (see [References](#references)).

## References

1. Köster, J. & Rahmann, S. (2012). Snakemake — a scalable bioinformatics workflow engine. *Bioinformatics*.
2. Callahan, B. J., et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*.
3. Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*.
4. Langmead, B. & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*.
5. de Sena Brandine, G. & Smith, A. D. (2019). Falco: high-speed FastQC emulation for fastq files. *F1000Research*.
6. Bengtsson-Palme, J., et al. (2015). Metaxa2: improved identification and taxonomic classification of small and large subunit rRNA in metagenomic data. *Molecular Ecology Resources*.
7. Bengtsson-Palme, J., et al. (2013). ITSx: improved software detection and extraction of ITS1 and ITS2. *Methods in Ecology and Evolution*.
8. Rognes, T., et al. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ*.
9. Quast, C., et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*.
10. Abarenkov, K., et al. (2024). UNITE general FASTA release for eukaryotes. *UNITE Community*.
11. Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*.
12. Bushnell, B. BBTools (BBDuk, BBMap). *DOE Joint Genome Institute*.
13. Wood, D. E., Lu, J. & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*.
14. Lu, J., et al. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*.
15. Lu, J., et al. (2022). Metagenome analysis using the Kraken software suite (KrakenTools). *Nature Protocols*.
16. Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*.

## License

MetaFlux is released under the [MIT License](LICENSE). Third-party tools invoked by
the workflow are distributed under their own licenses (a mix of MIT, BSD, and
GPL/LGPL); see each tool's repository.

[↑ Back to top](#table-of-contents)
