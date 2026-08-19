# Overview

MetaFlux is a [Snakemake](https://snakemake.github.io/) workflow that takes raw
paired-end short reads through to a taxonomically annotated abundance table. It
covers two kinds of study — marker-gene metabarcoding and shotgun metagenomics — and a
single key in the config file, `mode`, decides which of the two runs. This site
documents version **v2.3.1**.

MetaFlux belongs to the [BioFlux](https://github.com/stars/iLivius/lists/bioflux)
family of pipelines.

## The two modes

|                     | `amplicon`                                                                                                          | `shotgun`                                                                                       |
|---------------------|---------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| **Input**           | PCR products from one marker gene                                                                                     | untargeted whole-genome sequencing reads                                                          |
| **Core tools**      | [DADA2](https://benjjneb.github.io/dada2/)                                                                            | [Kraken2](https://github.com/DerrickWood/kraken2) + [Bracken](https://github.com/jenniferlu717/Bracken) |
| **Markers**         | 16S rRNA, ITS, 18S rRNA, gyrB, rpoB                                                                                   | not applicable                                                                                    |
| **Main steps**      | PhiX removal, primer trimming, exact ASV inference, optional marker-region extraction, length filtering, taxonomy      | decontamination, quality trimming, read classification, abundance re-estimation, OTU-table construction |
| **Reference data**  | SILVA (16S), UNITE (ITS), PR2 (18S), DD7RZ8 (gyrB), FROGS (rpoB) — fetched automatically on first use                  | a Kraken2 database, supplied by the user                                                          |
| **Final table**     | `6.taxonomy/asv_table.txt`                                                                                            | `03.abundance/otu_table.tsv`                                                                      |

**Amplicon mode** resolves sequences exactly rather than clustering them at an
arbitrary similarity threshold: DADA2 infers amplicon sequence variants (ASVs), which
are then length-filtered and classified against the reference database that matches the
marker. See the [amplicon overview](amplicon/overview.md) and the
[per-marker pages](amplicon/markers/index.md).

**Shotgun mode** skips amplification altogether and assigns whole-genome reads directly
against a Kraken2 index, with Bracken re-estimating abundances at the chosen taxonomic
level. See the [shotgun overview](shotgun/overview.md).

## One command, one config, one report

Both modes are launched the same way, from the repository root:

```bash
snakemake --sdm conda --cores 16 --configfile config/config.yaml
```

`--sdm conda` (short for `--software-deployment-method`) tells Snakemake to build each
rule's own conda environment on first run, so only Snakemake itself (version 9 or
newer) has to be installed by hand.

Everything else — mode, input and output paths, marker, primers, databases, per-rule CPU
and memory — lives in `config/config.yaml`. Both modes finish by writing read-tracking
statistics and a single [MultiQC](https://multiqc.info/) report under
`multiqc/multiqc_report.html`, so a run can be checked end to end from one HTML file.

!!! note "What has to be prepared in advance"

    Amplicon runs need only the primer FASTAs: the marker reference databases download
    themselves the first time they are needed and are cached under `refdb/`. Shotgun
    runs need a Kraken2 database on disk before the first launch — these are large and
    version-sensitive, so MetaFlux never fetches one automatically. See
    [Installation](getting-started/installation.md) and
    [Reference databases](shotgun/databases.md).

## Where to go next

| Page                                                        | Read it for                                                                 |
|-------------------------------------------------------------|-----------------------------------------------------------------------------|
| [Installation](getting-started/installation.md)              | cloning the repository, the Snakemake launcher environment, reference databases |
| [Quick start](getting-started/quick-start.md)                | the shortest path from a clone to a finished run                            |
| [Choosing a mode](getting-started/choosing-a-mode.md)        | which mode a given dataset belongs in, and what happens if the wrong one is set |
| [Amplicon mode](amplicon/overview.md)                        | the DADA2 path step by step, marker by marker                               |
| [Shotgun mode](shotgun/overview.md)                          | the Kraken2 + Bracken path step by step                                     |
| [Configuration](reference/configuration.md)                  | every config key, what consumes it, and its default                         |
| [Output files](reference/output.md)                          | what each output directory holds and which files are the deliverables       |
| [Running MetaFlux](reference/running.md)                     | dry runs, resuming, cluster execution                                       |
| [Troubleshooting](troubleshooting.md)                        | the failure modes that come up in practice, and how to get past them        |
| [Rationale](about/rationale.md)                              | why the workflow makes the choices it does                                  |
| [Citation and references](about/citation.md)                 | how to cite MetaFlux and the tools it calls                                 |

MetaFlux is released under the MIT License. The tools it invokes carry their own
licenses.
