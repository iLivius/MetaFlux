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
anything or writing any file. It is worth the ten seconds every time, because
MetaFlux does most of its validation while the workflow is being read — before any
job starts. A dry run therefore catches, among others:

- `mode` set to something other than `amplicon` or `shotgun`;
- an `amplicon.type` that matches no marker pack;
- primer FASTAs that do not exist (amplicon), or a `kraken_db` that does not exist (shotgun);
- a config still using the retired `include_pattern` / `exclude_pattern` filter keys;
- ITS combined with `length_filter.mode: auto` but a manual `expected_length`, which would otherwise crash *after* the full DADA2 run.

The same parse step prints a banner line to stderr:

```text
[MetaFlux v2.1.0] mode=amplicon, input=local, 24 sample(s) in /data/fastq
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
is read. Nested settings (`amplicon.type`, for instance) are awkward to override
this way; keeping one config file per marker and selecting it with `--configfile`
is clearer.

!!! warning "Give each mode its own `out_dir`"

    Both modes write `stats/read_tracking.txt` and `multiqc/multiqc_report.html`.
    Running the second mode into the first mode's output directory overwrites them.
    The same applies to running two markers into one directory: the ASV files in
    `5.dada2/` are not named after the marker.

**Limiting concurrency.** `--cores` is the total CPU budget on the machine; a rule
asking for more threads than that is scaled down to fit. Lowering `--cores` is
therefore the way to run fewer heavy jobs side by side. Memory is capped
separately: every heavy rule declares a `mem_mb` — from `resources.mem_mb.*` in the
config, except Kraken2, whose figure is computed while the workflow is read — so
`--resources mem_mb=32000` lets Snakemake start jobs only while their declared memory
still fits in 32 GB. This matters most for Kraken2 with
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
  runs.** Deleting `6.taxonomy/asv_table.txt` re-runs `assign_taxonomy`; the DADA2
  steps upstream are already satisfied and are skipped.
- **Changing the config can be enough to mark a step out of date.** Current
  Snakemake versions consider more than file timestamps: a changed rule parameter, a
  changed script, or a changed conda environment also count. Because the contaminant
  keep/discard lists reach `assign_taxonomy` as rule parameters, editing them
  re-runs the taxonomy step alone and leaves the denoising untouched. Adding
  `--rerun-triggers mtime` restricts the comparison to timestamps only.
- **A dry run answers the question directly.** `-n` prints exactly which jobs a
  re-run would execute, which is faster than reasoning about it.

Two caches deliberately sit outside this machinery and survive across runs:
reference databases under `refdb/`, and the amplicon-length probe result under
`refdb/cache/` inside the MetaFlux checkout — that second path is fixed and does
not follow `references.refdb_root`. The probe's filename encodes the marker, its
reference, and a hash of both primer sequences, so it is reused whenever those
are unchanged and recomputed when a primer changes. See
[output files](output.md) for what a run writes where.

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
