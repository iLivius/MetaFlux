# Quick start

The shortest path from an empty directory to a taxonomically annotated table. Five
steps: clone, create the launcher environment, edit the config, preview the plan,
run it.

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

Run these from the repository root: Snakemake looks for `workflow/Snakefile` there,
and relative paths in the config resolve against the directory snakemake is invoked
from, not against the config file. Absolute paths are safest.

Steps 1 and 2 are covered in more detail under [installation](installation.md),
including where the reference databases come from and what `--sdm conda` builds.

## Step 3 in detail — the four lines that define a run

`config/config.yaml` is long, but most of it is tuning that can stay as shipped.
Four settings decide what the run actually is:

| Setting | What it decides |
|---------|-----------------|
| `mode` | `amplicon` or `shotgun` — see [choosing a mode](choosing-a-mode.md) |
| `amplicon.type` | `16S`, `ITS`, `18S`, `gyrB` or `rpoB` (amplicon runs only) |
| `input.fastq_dir` | The directory holding the paired FASTQ files |
| `output.out_dir` | Where every result is written |

Then whatever the chosen mode needs before it can start: the two primer FASTA files
for amplicon (`amplicon.primers.fwd` / `.rev`), or the Kraken2 database directory
for shotgun (`references.kraken_db`).

Input files are found by pattern, not listed by hand. MetaFlux globs `fastq_dir` for
`{sample}_R1` / `{sample}_R2` pairs — shotgun also accepts the SRA-style
`{sample}_1` / `{sample}_2` — and the part before the mate suffix becomes the sample
name used in every output table. All files must share one extension, and in practice
that extension has to be `fastq` or `fastq.gz`. The startup check also lets `fq` and
`fq.gz` through, but every rule that reads a raw file resolves its path as
`.fastq.gz` or `.fastq`, so a directory of `.fq.gz` files passes validation and then
fails later on an input the rules cannot resolve.

!!! tip "Keep the template, edit a copy"

    The workflow always loads `config/config.yaml` as a base layer, and a file passed
    with `--configfile` is applied on top of it. So copying the template and editing
    the copy works cleanly:

    ```bash
    cp config/config.yaml config/config.local.yaml
    snakemake --sdm conda --cores 16 --configfile config/config.local.yaml
    ```

    `.gitignore` tracks only `config/config.yaml`; any other `config/*.yaml` stays
    out of git, so a run config with real paths in it is never committed by accident.

## Step 4 in detail — why the dry-run matters

`--dry-run` asks Snakemake to work out everything it would do and then stop without
running a single job. It costs seconds. Skipping it costs whatever the workflow
manages to do before it hits the problem — and the first thing a real run does is
build conda environments and download reference databases, which can take a long
time before anything fails.

Nearly all of MetaFlux's own validation happens while the plan is being built, so a
dry-run exercises it. What is checked before the first job starts:

- `mode` is one of `amplicon` / `shotgun`.
- `fastq_dir` yields at least one `{sample}_R1` file (or `{sample}_1` in shotgun
  mode), all with one consistent extension, and sample names without characters
  Snakemake cannot use as wildcards (a wildcard is the placeholder — here
  `{sample}` — that a rule is instantiated over, once per sample).
- Amplicon: `amplicon.type` matches a marker; both primer FASTA files exist;
  `taxonomy.method` is available for that marker (`gyrB` and `rpoB` have no SINTAX
  reference, so only `rdp` works for them); the statistic the length probe needs is
  present in `probe_length_stat`.
- Shotgun: the `kraken_db` directory exists.

Failures are reported plainly, for example:

```text
[MetaFlux] no *_1.* FASTQs found in /path/to/fastq_dir
[MetaFlux] mode=amplicon requires {sample}_R1/_R2 naming; found {sample}_1/_2 in /path/to/fastq_dir
[MetaFlux] mode=amplicon requires primer FASTA but not found: /path/to/fwd_primer.fasta
[MetaFlux] mode=shotgun requires kraken_db but not found: /path/to/kraken2_db
```

The first message names the fallback pattern, `_1`, because MetaFlux looks for
`_R1/_R2` first and only reports the failure after the SRA-style names have also
come up empty — it means neither convention matched anything in that directory.

Two things that check does not do. MetaFlux creates `fastq_dir` if it is missing,
because BioProject mode downloads into it, so a mistyped path produces the
`no *_1.* FASTQs found` message above rather than a missing-directory error. And
mates are matched from R1 only — a missing R2 is reported later by Snakemake as a
missing input file, not by MetaFlux.

The dry-run also prints a one-line summary of what MetaFlux understood, which is
worth reading before committing cores to a run:

```text
[MetaFlux v2.2.0] mode=amplicon, input=local, 24 sample(s) in /path/to/fastq_dir
```

A sample count that does not match the sequencing run usually means a naming
mismatch — an unexpected mate suffix, or files sitting one directory deeper than
`fastq_dir`. Below that summary, Snakemake lists the jobs it would run and how many
times each. A job count of zero means everything is already up to date; large
unexpected counts usually mean a config change invalidated an earlier stage.

!!! warning

    The dry-run cannot catch everything. It checks paths, names and parameter
    validity — not whether `mode` matches the biology of the data. Both modes accept
    `{sample}_R1/_R2` files, so amplicon reads pushed through shotgun mode (or the
    reverse) will run without complaint and produce a meaningless table. See
    [choosing a mode](choosing-a-mode.md).

## Step 5 — the run

```bash
snakemake --sdm conda --cores 16 --configfile config/config.yaml
```

`--cores` is the total number of CPUs Snakemake may use at once across all
concurrent jobs; the split between rules comes from the `resources` block in the
config. The first run also builds the conda environments and downloads the reference
databases the selected marker needs, so it is considerably slower than the second.

The run is resumable. If it is interrupted, repeating the same command picks up from
the last completed step; add `--rerun-incomplete` if Snakemake reports outputs left
in an incomplete state. Changing a parameter re-runs only the steps downstream of
it — trying a different truncation length, for example, does not repeat primer
trimming.

## Where the results are

Everything lands under `out_dir` in numbered, stage-ordered directories. The two
files most runs end at:

| Mode | Final table | Read tracking |
|------|-------------|---------------|
| amplicon | `6.taxonomy/asv_table.txt` | `stats/read_tracking.txt` |
| shotgun | `03.abundance/otu_table.tsv` (plus `otu_table.biom`) | `stats/read_tracking.txt` |

Both modes also write `multiqc/multiqc_report.html`, and `read_tracking.txt` shows
how many reads survived each stage per sample — the first thing to look at when a
table comes out thinner than expected. The full layout is in
[output files](../reference/output.md).

## Next

- [Choosing a mode](choosing-a-mode.md) — amplicon or shotgun.
- [Configuration](../reference/configuration.md) — every key, with its default.
- [Running MetaFlux](../reference/running.md) — cluster profiles, job limits, useful flags.
- [Troubleshooting](../troubleshooting.md) — the failures that actually happen.
