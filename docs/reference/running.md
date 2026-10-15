# Running MetaFlux

MetaFlux is launched with `snakemake`. Everything that varies between runs lives
in the config file, so the command line stays short and looks the same for both
modes.

## Run from the repository root

The Snakefile resolves `configfile: "config/config.yaml"` against the directory the
command was launched from, not against its own location. Two more lookups behave
the same way: `references.refdb_root` defaults to the relative path `refdb`, and
user-supplied marker packs are read from `config/markers/`. Launching from
somewhere else silently changes where all three point.

Paths inside the config (`input.fastq_dir`, `output.out_dir`, `references.kraken_db`,
the primer FASTAs) are best given as absolute paths, which makes them immune to this.

## The two commands

```bash
# Always preview the plan first
snakemake --sdm conda --cores 16 --configfile config/config.yaml --dry-run

# Execute
snakemake --sdm conda --cores 16 --configfile config/config.yaml
```

A dry run (`--dry-run`, or `-n`) builds the job list and prints it without running
anything or writing any file. Do one before every run: MetaFlux does most of its
checks while the workflow is read, before any job starts. A dry run therefore
catches, among others:

- `mode` set to something other than `amplicon` or `shotgun`;
- an `amplicon.type` that matches no marker pack;
- primer FASTAs that do not exist (amplicon), or a `kraken_db` that does not exist (shotgun);
- a config still using the retired `include_pattern` / `exclude_pattern` filter keys;
- a keep/discard token that can never match the run's taxonomy strings, such as
  `k__Bacteria` from a config written for MetaFlux v2.4.0 or earlier (now `d__Bacteria`);
- a marker pack in `config/markers/` written for an older MetaFlux;
- ITS combined with `length_filter.mode: auto` but a manual `expected_length`, which would otherwise crash *after* the full DADA2 run.

While reading the workflow, MetaFlux also prints a banner line to stderr:

```text
[MetaFlux v2.5.0] mode=amplicon, input=local, 24 sample(s) in /data/fastq
```

Before it, an amplicon run prints a line saying whether the organelle filter is on, and
a 16S run one more naming the SILVA release in use — the quickest check that
`references.silva.release` and `amplicon.taxonomy.silva_other` were picked up:

```text
[MetaFlux] organelle filter: on — ASVs that Metaxa2 calls mitochondria or chloroplast are removed
[MetaFlux] 16S reference: SILVA 144 (silva_other: verbatim; trainset silva_nr99_v144_toGenus_7rank_other-verbatim_trainset.fa.gz)
```

Checking that the sample count is the expected one costs nothing and catches a
mistyped `fastq_dir` or a mate-naming mismatch before a long run starts.

!!! note "`--configfile` is optional"

    The Snakefile already loads `config/config.yaml`. Passing `--configfile` makes
    the choice explicit, and is how a different config — one per marker, or one per
    project — is selected.

## `--sdm conda` and `--conda-prefix`

`--sdm conda` is short for `--software-deployment-method conda`. It tells Snakemake
to build each rule's own conda environment from the YAML files in
`workflow/envs/` and to run that rule inside it. Only Snakemake itself has to be
installed by hand; DADA2, Cutadapt, Kraken2, Bracken and the rest are installed by
the workflow on first use. Leaving the flag out makes every rule run against
whatever happens to be on `PATH`, which is exactly what the per-rule environments
exist to avoid.

The environments are built once and then reused. By default they go into
`.snakemake/conda/` inside the working directory, which is why running from the
repository root keeps reusing the same ones. `--conda-prefix` moves them somewhere
else:

```bash
snakemake --sdm conda --conda-prefix ~/.metaflux-envs --cores 16
```

That is useful when several checkouts, several users, or several working
directories should share one set of environments instead of each building its own —
they are a few gigabytes in total, and building them takes minutes.

!!! tip "Pre-build the environments"

    `--conda-create-envs-only` builds every environment and stops. Handy before a
    cluster run, so the compute nodes are not each waiting on conda, and on
    machines where the network is only reachable from the login node.

## Useful patterns

```bash
# Override the mode without editing the config
snakemake --sdm conda --cores 16 --config mode=shotgun

# Cap concurrent memory (e.g. memory-bound steps)
snakemake --sdm conda --cores 16 --resources mem_mb=32000

# Resume after an interruption
snakemake --sdm conda --cores 16 --rerun-incomplete

# Cluster execution: supply a profile (e.g. SLURM) — per-rule threads/mem_mb are honoured
snakemake --sdm conda --profile <slurm-profile>
```

**Overriding the mode.** `--config mode=shotgun` replaces the top-level `mode` key
for that invocation only; the file on disk is untouched and every other setting
still comes from it. The shotgun section of the config has to be filled in for this
to work — `references.kraken_db` in particular, which is checked while the workflow
is read. Nested settings (`amplicon.type`, for instance) can't be reached with a dotted
key — `--config` only accepts plain top-level identifiers, so overriding one
means replacing its whole parent block with inline YAML
(`--config 'amplicon={type: ITS}'`) rather than a single `key=value` pair.
Keeping one config file per marker and selecting it with `--configfile` is
clearer.

!!! warning "Give each mode its own `out_dir`"

    Both modes write `stats/read_tracking.txt` and `multiqc/multiqc_report.html`.
    Running the second mode into the first mode's output directory overwrites them.
    The same applies to running two markers into one directory: the ASV files in
    `05.dada2/` are not named after the marker.

**Limiting concurrency.** `--cores` is the total CPU budget on the machine; a rule
asking for more threads than that is scaled down to fit. Lowering `--cores` is
therefore the way to run fewer heavy jobs side by side. Memory is capped
separately: every heavy rule declares a `mem_mb` — from `resources.mem_mb.*` in the
config, except Kraken2 and `sylph_profile`, whose figures are computed while the
workflow is read — so `--resources mem_mb=32000` lets Snakemake start jobs only while
their declared memory still fits in 32 GB. This matters most for Kraken2 with
`shotgun.kraken.memory_mapping: false`, where every concurrent process loads its own
full copy of the database, and for the DADA2 steps, which hold all samples in
memory at once.

`--jobs` is not a second lever here. On a local run Snakemake treats it as an alias
for `--cores`, and passing both switches the core budget off rather than adding a
job cap on top of it — `--jobs 4 --cores 16` would run four jobs at once whatever
their thread counts add up to. It belongs on cluster runs only, see below.

**Resuming.** An interrupted run leaves two kinds of debris. Files that a job was
part-way through writing are marked incomplete, and `--rerun-incomplete` rebuilds
them instead of trusting them. A run killed outright (`Ctrl-C` twice, a node
failure, an OOM kill) can also leave the working directory locked, which Snakemake
reports on the next launch; `snakemake --unlock` clears it. Everything that
finished cleanly is kept and not recomputed.

**Cluster execution.** A profile holds the submission settings for the scheduler
(partition, account, how `threads` and `mem_mb` map to `sbatch` arguments), so the
MetaFlux command line does not change. The per-rule values come from
`resources.threads.<rule>` and `resources.mem_mb.<rule>` in the config, falling back
to `threads_default` (4) and `mem_mb_default` (2000) for any rule not listed.
Kraken2 is the exception: its memory request is computed while the workflow is read,
from `shotgun.kraken.memory_mapping` and the actual size of the database's
`hash.k2d` — 20000 MB with memory mapping on (the database itself is shared through
the operating system's file cache), otherwise the hash size plus 10 % headroom. With
a profile, `--jobs` becomes the cap on jobs queued with the scheduler at once —
the one setting where it does something `--cores` does not.

## What happens on a re-run

Snakemake only rebuilds what is missing or out of date. Re-running a completed
workflow with no changes does nothing at all and says so. Three practical
consequences:

- **A deleted or edited output is rebuilt, and only what is needed to rebuild it
  runs.** Deleting `06.taxonomy/asv_table.txt` re-runs `assign_taxonomy`; the DADA2
  steps upstream are already satisfied and are skipped.
- **Changing the config can be enough to mark a step out of date.** Current
  Snakemake versions consider more than file timestamps: a changed rule parameter, a
  changed script, or a changed conda environment also count. Because the contaminant
  keep/discard lists reach `assign_taxonomy` as rule parameters, editing them
  re-runs the taxonomy step alone and leaves the denoising untouched. Adding
  `--rerun-triggers mtime` restricts the comparison to timestamps only.
- **A dry run shows what would run.** `-n` prints the jobs a re-run would execute.

!!! note "Re-running an output folder made with v2.4.0"

    **Every step runs again.** The output folders have new names: the amplicon folders
    `1.reads` to `7.phylogeny` are now `01.reads` to `07.phylogeny`, and the read-count
    files in `stats/` take the same two digits (`01.raw_reads.counts` …); in shotgun mode
    the Kraken2 files move from `02.classification/` into `02.classification/kraken2/`,
    together with the k-mer evidence files of `02b.evidence/`. The few steps whose files
    kept their names run again too, because MetaFlux now ships a pin file next to every
    conda environment (the exact package builds it was tested with; see
    [Installation](../getting-started/installation.md#3-what-sdm-conda-does)) and
    Snakemake counts that file as part of the software a step ran with: the QC of the raw
    reads in `stats/`, the primer reverse complements in `_aux/`, the reference files
    already in `refdb/` (downloaded and built again, same content) and, in BioProject
    mode, the download of the SRA runs. Only `samples.tsv` is kept, the one step without
    a conda environment whose file kept its name. To keep the SRA runs and the reference
    files already on disk, start the run once with `--rerun-triggers mtime params input
    code`: it leaves the software-environment trigger out, and everything that moved runs
    again anyway. The old folders and count files, the old files directly in
    `02.classification/`, and `03.abundance/otu_table.raw.biom` (v2.4.0's intermediate
    table, which a run made with `--notemp` keeps) are left as they are and can be
    deleted.

    **Results that change on purpose:**

    - Taxonomy strings: `d__` is always the domain slot and `k__` the kingdom slot (see
      [taxonomy strings](output.md#taxonomy-strings-in-the-final-tables)). 16S moves to
      SILVA 144, the new default: `d__Bacteria;k__Pseudomonadati;p__…`. On SILVA 138.2
      (`references: {silva: {release: "138.2"}}`, the release v2.4.0 used) and on rpoB:
      `d__Bacteria;p__…`. ITS, 18S and gyrB keep their ranks. The shotgun `otu_table`
      has eight slots on every Kraken2 database, a virus carries its realm in `d__`,
      species names are whole, and in `otu_table.biom` the sample `Id` is the plain
      sample name (it was `<sample>_report`).
    - 16S on SILVA 138.2: the length probe now leaves chloroplast, mitochondrial and
      eukaryote references out. For 799F/1175R, 515F/806R and 341F/805R the probe's p95
      does not change, so neither do `truncLen` and the denoised ASVs and counts; the
      `auto` ASV length window can move by a base (for 799F/1175R, the primers of the
      16S [test dataset](../about/test-datasets.md), its lower edge moves from 324 to
      325 bp).
    - 18S with `sintax`: almost every ASV whose taxonomy was empty now carries its
      domain, nearly always `d__Eukaryota` (158 of 159 on the 18S
      [test dataset](../about/test-datasets.md)).
    - ITS with `rdp`: the species reads `Genus_epithet`, as with `sintax`. All taxonomy
      tables have one format under both methods (`""` as first header cell, `NA` for a
      rank not reached), and the rank cells of `taxon_seq_table.txt` hold bare names.
    - The taxon filter can remove more organelle ASVs (best-hit check, organelle
      filter). A lineage near the cutoff can gain or lose a rank (VSEARCH 2.32 under
      `sintax`; under `rdp`, as between two identical runs), so an ASV at the cutoff of a
      `keep` token's rank can be kept by one version and removed by the other.
    - Shotgun: Kraken2 2.17.1, Bracken 3.1p1 and bbmap 39.52 are pinned, so a folder made
      with other releases can get slightly different results; at a `tax_lev` other than
      `S` the OTU table is now built correctly (see
      [Troubleshooting](../troubleshooting.md#shotgun-empty-or-nearly-empty-otu-table-with-tax_lev-other-than-s));
      `nophix` in the read tracking now counts exactly the read pairs BBDuk kept, so
      `raw` minus `nophix` is the pairs it removed (a pair with PhiX in one mate only
      used to count as half a pair).

    **Stops at start-up:** a filter token such as `k__Bacteria` (the message names the
    one to write, `d__Bacteria`); `shotgun.bracken.tax_lev: D`, refused like `K`; a
    marker pack in `config/markers/` written for an older MetaFlux (copy the shipped pack
    again and re-apply your change). See [Troubleshooting](../troubleshooting.md).

    **Reference files:** the SILVA 144 files are downloaded and built once, and the 18S
    PR2 SINTAX file is built once (neither existed in v2.4.0); the rpoB reference is
    built again under a new name (the FROGS rpoB archive is downloaded again); the SILVA
    138.2 SINTAX file is rebuilt under its old name, with the same content, so a refdb
    folder shared with a v2.4.0 checkout keeps working for both.

Two caches deliberately sit outside this machinery and survive across runs:
reference databases under `refdb/`, and the amplicon-length probe result under
`refdb/cache/` inside the MetaFlux checkout — that second path is fixed and does
not follow `references.refdb_root`. The probe's filename encodes the marker, its
reference (for 16S also the lineages left out of the statistics; for rpoB
`_excl-second-site`, because products that hold a second forward-primer site are left
out), and a hash of both primer sequences, so it is reused whenever those are unchanged
and recomputed when a primer changes. See [output files](output.md) for what a run
writes where.

## Two mode-specific notes

!!! note "Shotgun, BioProject mode"

    Setting `input.bioproject` to an NCBI accession makes MetaFlux resolve the
    project's PAIRED runs through NCBI e-utils and download them into `fastq_dir`
    before anything else. Single-end runs are skipped — they do not fit the rest of
    the pipeline. The run list is cached as `{accession}_runinfo.csv` in `fastq_dir`,
    so later invocations do not need the network for that step, and
    `input.accession_list` narrows the set to named runs.

!!! note "Shotgun, first Kraken2 call"

    The first Kraken2 job reads the whole database off disk. That is slow while the
    cache is cold and fast once the operating system has the database in memory, so
    a first sample that seems to hang is usually just loading.

## See also

- [Configuration](configuration.md) — every key referred to here.
- [Output files](output.md) — what the run produces and which file feeds which step.
- [Troubleshooting](../troubleshooting.md) — YAML errors, empty tables, and the `truncLen` trap.
