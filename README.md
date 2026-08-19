# MetaFlux

```
______  ___    _____       _______________
___   |/  /______  /______ ___  ____/__  /___  _____  __
__  /|_/ /_  _ \  __/  __ `/_  /_   __  /_  / / /_  |/_/
_  /  / / /  __/ /_ / /_/ /_  __/   _  / / /_/ /__>  <
/_/  /_/  \___/\__/ \__,_/ /_/      /_/  \__,_/ /_/|_|

v2.3.1
```

**Unified short-read multi-marker amplicon and shotgun taxonomic profiling workflow.**

![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A59.0-brightgreen.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![DOI](https://img.shields.io/badge/DOI-pending%20release-lightgrey.svg)

## 📖 Documentation

**Full documentation: [iLivius.github.io/MetaFlux](https://iLivius.github.io/MetaFlux/)**

| | |
|---|---|
| [Installation](https://iLivius.github.io/MetaFlux/getting-started/installation/) | environments, reference databases, Kraken2 index |
| [Quick start](https://iLivius.github.io/MetaFlux/getting-started/quick-start/) | first run, start to finish |
| [Amplicon mode](https://iLivius.github.io/MetaFlux/amplicon/overview/) | the DADA2 pipeline, step by step |
| [Markers](https://iLivius.github.io/MetaFlux/amplicon/markers/) | 16S · ITS · 18S · gyrB · rpoB, one page each |
| [Shotgun mode](https://iLivius.github.io/MetaFlux/shotgun/overview/) | Kraken2 + Bracken profiling |
| [Configuration reference](https://iLivius.github.io/MetaFlux/reference/configuration/) | every key, default and trade-off |
| [Output files](https://iLivius.github.io/MetaFlux/reference/output/) | what each run writes, and what reads it |
| [Troubleshooting](https://iLivius.github.io/MetaFlux/troubleshooting/) | the failures worth knowing about in advance |
| [Test datasets](https://iLivius.github.io/MetaFlux/about/test-datasets/) | the mock communities each marker is validated against |

This README is deliberately short: it covers what MetaFlux is and how to get a first run
going. Everything else — parameter reference, per-marker guidance, output semantics and
the methodological caveats — lives on the documentation site.

## Synopsis

MetaFlux is a [Snakemake](https://snakemake.github.io/) workflow that takes raw
paired-end short reads through to a taxonomically annotated abundance table. It
runs in one of two modes, selected by a single config key:

- **`amplicon`** — a [DADA2](https://benjjneb.github.io/dada2/)-based pipeline for
  multi-marker metabarcoding (16S rRNA, ITS, 18S rRNA, gyrB, rpoB): PhiX removal,
  primer trimming, exact ASV inference, optional marker-region extraction, length
  filtering, and taxonomic assignment against SILVA (16S), UNITE (ITS), PR2 (18S),
  DD7RZ8 (gyrB), or FROGS (rpoB).
- **`shotgun`** — a [Kraken2](https://github.com/DerrickWood/kraken2) +
  [Bracken](https://github.com/jenniferlu717/Bracken) pipeline for shotgun
  metagenomic taxonomic profiling: decontamination, quality trimming, read
  classification, abundance re-estimation, and OTU-table construction.

Both modes run from the same command and the same config file, and produce a single
[MultiQC](https://multiqc.info/) report.

```
                 ┌──────────────┬───────────────┐
   config: mode  │   amplicon   │    shotgun    │
                 └──────┬───────┴───────┬───────┘
        raw paired-end short reads (R1/R2)
                        │               │
        bowtie2 PhiX*   │               │  BBDuk PhiX*  +  BBMap host*
        Cutadapt primer │               │  fastp trim
        DADA2 ASVs      │               │  Kraken2 classify
        Metaxa2 / ITSx* │               │  Bracken abundance
        length filter   │               │  kraken-biom OTU table
        multi-marker    │               │  KrakenTools per-taxon reads*
        taxonomy        │               │
                        ▼               ▼
              ASV table + taxonomy   OTU table + reports
                        └───────┬───────┘
                     read tracking + MultiQC

   * optional steps, toggled in config
```

What the workflow adds over running these tools by hand: `truncLen` chosen from the
run's own quality profile under a merge-overlap constraint; expected amplicon length
measured from the reference database with the supplied primers rather than hard-coded;
orientation-aware primer trimming; per-rule CPU and memory in one place; and one
reproducible entry point for both modes.
See [Rationale](https://iLivius.github.io/MetaFlux/about/rationale/).

`MetaFlux` belongs to the [BioFlux](https://github.com/stars/iLivius/lists/bioflux) family of pipelines.

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

Only Snakemake (≥ 9) is installed by hand. `--sdm conda`
(software-deployment-method) tells Snakemake to build each rule's own conda
environment on first run, so every other tool is provisioned automatically.

Always dry-run first: it validates the config, resolves sample names and reference
paths, and prints the plan without executing anything.

The shipped `config/config.yaml` is a template whose `/path/to/...` entries are
placeholders. A common pattern is to copy it to `config/config.local.yaml`
(git-ignored) and run with `--configfile config/config.local.yaml`, leaving the
template untouched.
Full details in the [configuration reference](https://iLivius.github.io/MetaFlux/reference/configuration/).

## Reference databases

Amplicon reference databases download automatically the first time they are needed and
are cached under `refdb/`. Shotgun mode needs a Kraken2 index supplied by hand.

| Database | Mode | Provisioning |
|---|---|---|
| PhiX | both | auto-fetched (NCBI) for amplicon; BBDuk's bundled copy for shotgun |
| SILVA (toGenus + species + SINTAX) | amplicon (16S) | auto-fetched on first run |
| UNITE (general + UCHIME + SINTAX) | amplicon (ITS) | auto-fetched on first run |
| PR2 (DADA2 + UTAX) + SILVA-Euk | amplicon (18S) | auto-fetched on first run |
| DD7RZ8 (DADA2 trainset) | amplicon (gyrB) | auto-fetched on first run |
| FROGS RefSeq (converted to DADA2) | amplicon (rpoB) | auto-fetched and converted on first run |
| **Kraken2 / Bracken index** | shotgun | **user-supplied** |
| Host genome(s) | shotgun (optional) | user-supplied, path or URL |

Choosing and downloading a Kraken2 index, and the memory it needs, are covered in
[reference databases](https://iLivius.github.io/MetaFlux/shotgun/databases/).

## Citation

Publications that use MetaFlux should cite this repository. A formal release with a
Zenodo DOI is forthcoming:

> Antonielli, L. (2026). *MetaFlux: a unified short-read multi-marker amplicon and
> shotgun taxonomic profiling workflow.* Zenodo. DOI: pending release.

MetaFlux is a wrapper around published tools and reference databases, and those do the
actual work — cite them too. The full reference list is on the
[citation page](https://iLivius.github.io/MetaFlux/about/citation/).

## Acknowledgements

Developed at the [AIT Austrian Institute of Technology](https://www.ait.ac.at/).
MetaFlux consolidates and modernises methods refined across many amplicon and
metagenomics collaborations, and is part of the **BioFlux** family of workflows.

Portions of this codebase were developed with the assistance of Claude Code.

## License

MetaFlux is released under the [MIT License](LICENSE). Third-party tools invoked by
the workflow are distributed under their own licenses (a mix of MIT, BSD, and
GPL/LGPL); see each tool's repository.
