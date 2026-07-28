# Installation

MetaFlux is a Snakemake workflow, not a compiled program: the clone contains rule
files, scripts and marker definitions, and nothing needs to be installed
system-wide. Only Snakemake has to be present before the first run. Every tool the
workflow actually calls — DADA2, Cutadapt, bowtie2, Kraken2, Bracken and the rest —
is installed by conda while the workflow runs.

## 1. Clone the repository

```bash
git clone https://github.com/iLivius/MetaFlux.git
cd MetaFlux
```

The clone is the working copy. Four places in it matter:

| Path | What is in it |
|------|---------------|
| `config/config.yaml` | The configuration template — every parameter of a run |
| `workflow/rules/` | The rule files, split into `shared/`, `amplicon/`, `shotgun/` |
| `workflow/envs/` | One conda environment specification per tool stack |
| `workflow/markers/` | One file per amplicon marker (`16S`, `ITS`, `18S`, `gyrB`, `rpoB`), declaring that marker's reference databases and download URLs |

A fifth directory, `refdb/`, appears on the first run: it is where downloaded
reference databases are cached. It is git-ignored, so it never ends up in a commit.

## 2. Create the launcher environment

Snakemake 9 or newer:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

This environment is only a launcher. It reads the config, plans the run and starts
jobs; it does not contain any of the analysis software. There is no reason to add
DADA2 or Kraken2 to it, and doing so does not help — the rules will not see them.

## 3. What `--sdm conda` does

Every MetaFlux run is started with `--sdm conda`. `--sdm` is short for
`--software-deployment-method`, the modern spelling of the deprecated `--use-conda`
flag:

```bash
snakemake --sdm conda --cores 16 --configfile config/config.yaml
```

Each rule names an environment file in `workflow/envs/` — the DADA2 rules point at
`dada2.yaml`, the Kraken2 rule at `kraken2.yaml`, and so on. With `--sdm conda`,
Snakemake builds each of those environments the first time a rule that needs it is
about to run, then reuses them on every later run. The environments are created
under `.snakemake/conda` in the directory snakemake was started from, which is also
git-ignored.

Two practical consequences:

- **The first run is slow before it is even doing science.** Solving and downloading
  several conda environments takes minutes; afterwards they are already on disk and
  the workflow starts immediately.
- **Tools never have to coexist.** The R/DADA2 stack, Cutadapt, Kraken2 and MultiQC
  each get their own environment, so an upgrade to one cannot break another. The
  software set is declared in the repository, but only partly pinned: some tools
  are fixed to an exact version, others only to a minimum, and a few are left to
  the solver.

!!! note

    Only the environments the selected mode actually needs are built. The Snakefile
    includes `workflow/rules/amplicon/*.smk` **or** `workflow/rules/shotgun/*.smk`
    depending on `mode`, so an amplicon run never solves the Kraken2 environment and
    a shotgun run never solves the DADA2 one.

## 4. Reference databases

Taxonomic assignment needs a reference. MetaFlux splits these into two groups: the
marker databases, which it downloads itself, and the shotgun databases, which are
too large and too version-sensitive to fetch behind the scenes and must be supplied
by hand.

| Database                           | Mode               | Provisioning                                                           |
|------------------------------------|--------------------|------------------------------------------------------------------------|
| PhiX                               | both               | **Auto-fetched** (NCBI) for amplicon; BBDuk's bundled copy for shotgun |
| SILVA (toGenus + species + SINTAX) | amplicon (16S)     | **Auto-fetched + converted** — toGenus and species from Zenodo, SINTAX built locally from the trainset |
| UNITE (general + UCHIME + SINTAX)  | amplicon (ITS)     | **Auto-fetched** from PlutoF on first run                              |
| PR2 (DADA2 + UTAX) + SILVA-Euk     | amplicon (18S)     | **Auto-fetched** from GitHub/Zenodo on first run                       |
| DD7RZ8 (DADA2 trainset)            | amplicon (gyrB)    | **Auto-fetched** from Recherche Data Gouv                              |
| FROGS RefSeq (converted to DADA2)  | amplicon (rpoB)    | **Auto-fetched + converted** on first run                              |
| Kraken2 database                   | shotgun            | **User-supplied** (see below)                                          |
| Host genome(s)                     | shotgun (optional) | **User-supplied** (path or URL)                                        |

### Marker databases (fetched automatically)

Each marker file in `workflow/markers/` declares the database files that marker
needs, where they live under `refdb_root`, and the URL to get them from. The
workflow turns a plain download straight into a fetch rule; references that arrive
inside an archive, or in a format DADA2 or VSEARCH cannot read, keep a dedicated
rule that unpacks or reshapes them.

Only the databases of the marker set in `amplicon.type` are downloaded. Running 16S
does not pull UNITE or PR2.

| Marker | Files created under `refdb_root` | Where they come from |
|--------|----------------------------------|----------------------|
| 16S | `silva/silva_nr99_v138.2_toGenus_trainset.fa.gz`, `silva/silva_v138.2_assignSpecies.fa.gz`, `silva/silva_nr99_v138.2_sintax.fa.gz` | Zenodo; the SINTAX file is not downloaded but built locally from the trainset (`convert_silva_sintax`) |
| ITS | `unite/sh_general_release_dynamic_s_all_19.02.2025.fasta`, `unite/utax_reference_dataset_all_19.02.2025.fasta.gz`, `unite/uchime_its1.fasta`, `unite/uchime_its2.fasta` | PlutoF; the general release arrives as a `.tgz` and the UCHIME subregions as a `.zip`, both unpacked by `fetch_unite` / `fetch_uchime` |
| 18S | `silva_euk/silva_132.18s.dada2.fa.gz`, `pr2/pr2_SSU_dada2.fasta.gz`, `pr2/pr2_SSU_UTAX.fasta.gz` | Zenodo (SILVA-Euk v132) and the PR2 v5.1.1 GitHub release |
| gyrB | `gyrb/train_set_gyrB_v6.fa.gz` | Recherche Data Gouv (INRAE Dataverse), DOI 10.57745/DD7RZ8 |
| rpoB | `rpob/rpob_refseq_cc_20240707_dada2.fa.gz` | FROGS rpoB databank (INRAE Toulouse); downloaded as a `.tar.gz` and rewritten into a DADA2 trainset by `convert_rpob_to_dada2` |

`refdb_root` defaults to `refdb`, i.e. inside the clone:

```yaml
references:
  refdb_root: refdb            # where fetched databases are cached
```

Point it somewhere else — a shared filesystem, say — and every marker database
follows. The `references` block lists no marker databases at all, because the marker
files already know their own paths. An entry here is only needed to **relocate** a
single database, for example when a copy already exists on the machine:

```yaml
references:
  silva: {train: /shared/dbs/silva_nr99_v138.2_toGenus_trainset.fa.gz}
```

Downloads happen once. A database is a normal Snakemake output: present on disk, its
fetch rule does not run again, so a second project reusing the same `refdb_root`
starts classifying immediately.

!!! note "The probe cache is not under `refdb_root`"

    The amplicon-length probe result is cached separately, always in `refdb/cache/`
    inside the repository, whatever `refdb_root` is set to. For 16S, 18S and rpoB the
    probe is an in-silico PCR of the primers against a full-length reference; for ITS
    and gyrB it is a direct length measurement of a reference that is already
    amplicon-length. Its filenames carry the marker, the reference version and a
    12-character fingerprint of the two primer FASTA files, so a change of primers
    produces a different name and the probe is recomputed rather than silently
    reused. See [amplicon length and truncation](../amplicon/length-and-truncation.md).

### PhiX

PhiX is spiked into most Illumina runs and both modes can remove it, but they get
the reference differently:

- **Amplicon** — `fetch_phix` downloads the NCBI phiX174 genome
  (`GCF_000819615.1`) to `refdb/phix/phix.fna` and `build_phix_index` builds the
  bowtie2 index next to it. Both run automatically on the first run with
  `amplicon.decontamination.remove_phix: true`.
- **Shotgun** — nothing is downloaded. BBDuk ships `phix174_ill.ref.fa.gz` inside
  its own conda package and MetaFlux uses that.

### Kraken2 database (shotgun, required)

Shotgun mode has no default database and does not fetch one: Kraken2 indexes are
large, dated, and the choice of index determines what can be detected at all. The
path goes in the config:

```yaml
references:
  kraken_db: /path/to/db/k2_pluspf   # [shotgun] — user-supplied
```

This must be the **directory** holding the index files (`hash.k2d`, `opts.k2d`,
`taxo.k2d` and the `*.kmer_distrib` files Bracken needs), not one of the files. If
the path does not exist, MetaFlux stops before any job runs, with:

```text
[MetaFlux] mode=shotgun requires kraken_db but not found: /path/to/db/k2_pluspf
```

Pre-built indexes, with sizes, build dates and md5 checksums, are listed in Ben
Langmead's [catalog](https://benlangmead.github.io/aws-indexes/k2). A standard
choice is the **PlusPF** index. A single `wget` stream works but is slow on a
~100 GB file; a multi-connection downloader is faster and resumable.

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

`aria2c` is an external download helper, not part of MetaFlux. If the AWS CLI is
installed, `aws s3 cp --no-sign-request <url>` is the fastest option (parallel
multipart). Plain `wget -c <url>` also works and resumes a partial download.

!!! warning "Memory, not just disk"

    The full PlusPF index is ~100 GB. With memory-mapping enabled
    (`shotgun.kraken.memory_mapping: true`, the default), Kraken2 reads the database
    through the operating system's file cache: the memory request is a fixed ~20 GB
    of private workspace per process, and parallel jobs share the cached database
    instead of each loading a copy. With memory-mapping off, each process loads the
    whole hash privately and MetaFlux sizes the request from the actual `hash.k2d`
    file plus 10%. Index size therefore matters only in that second case: a capped
    index such as `pluspf_16gb` lowers the request when memory-mapping is off, but
    with it on the request stays at the fixed ~20 GB. Details in
    [shotgun reference databases](../shotgun/databases.md).

### Host genome (shotgun, optional)

Host removal is off unless `shotgun.decontamination.host_genomes` lists something.
Each entry is a local FASTA path or a URL; several entries are concatenated and
indexed once.

For human material, use a **masked** reference. Masking ribosomal and low-complexity
regions, and sequences homologous to microbes, plants and fungi, stops real
microbial reads from being thrown away with the host. The masked hg19 + viral
reference (Handley 2020) is a good choice:

```bash
# hg19 + viral masking, ~889 MB
wget -O /path/to/db/human_virus_masked.fasta.gz \
  "https://zenodo.org/records/4116107/files/human_virus_masked.fasta.gz"
```

A URL can be given directly in `host_genomes`, but downloading once to a local path
and pointing at the file avoids re-fetching ~900 MB on every fresh run. See
[decontamination](../shotgun/decontamination.md).

## Next

- [Quick start](quick-start.md) — the shortest path from a clone to a first result.
- [Choosing a mode](choosing-a-mode.md) — amplicon or shotgun, and what each needs.
- [Configuration](../reference/configuration.md) — every config key in detail.
