# `.test/`: data and configs for the CI dry runs

This folder holds what `.github/workflows/ci.yml` needs to dry-run MetaFlux in every
mode without downloading anything: a few hundred read pairs per sample, the primer
files, stand-in databases that satisfy the start-up checks, and one run config per
case. Nothing here is meant for a real run, and nothing here is a reference.

## Why dry runs

`snakemake -n` reads the workflow, runs every start-up check (the marker pack, the
reference rules, the filter tokens, the database files, the sample names) and lists
the jobs it would run, without running one. That takes a few seconds per config and
needs no conda environment, no reference database and no Kraken2 index, so it can
run on GitHub's free runners for every push. It proves that every marker, taxonomy
method, SILVA release and shotgun profiler still builds its plan, and that the
start-up checks still stop what they must stop. It does not prove that the tools
produce the right results: that is what the real test datasets are for
(`docs/about/test-datasets.md`).

## Running it yourself

From the repository root, with Snakemake on your PATH (any way of installing it will
do; the CI uses `pip install snakemake==9.23.1`):

```bash
bash .test/run_dry_runs.sh
```

One line per config says PASS or FAIL with the number of jobs planned or the stop
message; the full Snakemake output of each config is kept under `.test/out/dry_run_logs/`
(git-ignored). The exit code is the number of failures. To dry-run one config by hand:

```bash
snakemake -n --cores 1 --configfile .test/config/16S_rdp_silva144.yaml
```

## What is in here

| Path | What it is |
|---|---|
| `config/*.yaml` | One config per case that must pass. Each sets only the keys that differ from the template `config/config.yaml`, which Snakemake loads underneath it as in a real run; the header comment of each file says what the case covers. |
| `config/expected_stop/*.yaml` | Configs that must stop at start-up. Each carries a `# EXPECT:` line with the fragment of the stop message the runner looks for. |
| `data/<marker>/` | The first 300 read pairs of two samples of each public test dataset (16S, ITS, 18S, gyrB, rpoB) and of two wastewater samples for shotgun mode, as `{sample}_R1/_R2.fastq.gz`. Enough for Snakemake to find the samples; far too few for any result to mean anything. |
| `primers/<marker>/` | The forward and reverse primer FASTA of each marker, as used on the full test datasets. |
| `stubs/kraken2_db/` | A stand-in Kraken2 database folder (see below). |
| `stubs/sylph/`, `stubs/sylph_taxonomy/` | Stand-in sylph database files and sylph-tax lineage tables (see below). |
| `stubs/host/host_standin.fa.gz` | A made-up 300 bp "host genome", so that the host-removal rules appear in one shotgun plan. |
| `run_dry_runs.sh` | The runner the CI calls; its header says exactly what it does. |

## The matrix

Amplicon, one config per marker, method and SILVA release that MetaFlux offers:

| Config | Marker | Method | SILVA release | Also covers |
|---|---|---|---|---|
| `16S_rdp_silva144` | 16S | rdp | 144 (default) | Metaxa2 extraction, PhiX removal, the organelle filter |
| `16S_sintax_silva144` | 16S | sintax | 144 | the phylogeny stage with IQ-TREE |
| `16S_rdp_silva138` | 16S | rdp | 138.2 | the phylogeny stage with FastTree |
| `16S_sintax_silva138` | 16S | sintax | 138.2 | |
| `ITS_rdp`, `ITS_sintax` | ITS (ITS2) | rdp, sintax | – | ITSx extraction, the UNITE downloads |
| `18S_rdp`, `18S_sintax` | 18S | rdp, sintax | – | the PR2 downloads and the PR2 SINTAX conversion |
| `gyrB_rdp` | gyrB | rdp | – | fixed amplicon length, manual truncation and length window |
| `rpoB_rdp` | rpoB | rdp | – | the FROGS archive download and conversion |

gyrB and rpoB have no SINTAX reference, so `method: sintax` is refused for them at
start-up; that refusal is one of the expected-stop cases rather than a matrix cell.

Shotgun, one config per profiler:

| Config | Profiler | Also covers |
|---|---|---|
| `shotgun_kraken2` | Kraken2 + Bracken | the k-mer evidence gate, the taxon filter, read extraction |
| `shotgun_sylph` | sylph | two databases with built-in sylph-tax taxonomies |
| `shotgun_both` | both | host removal with the stand-in host genome, the taxon filter, read extraction |

Expected stops (`config/expected_stop/`):

| Config | What must stop the run | Expected message fragment |
|---|---|---|
| `16S_old_tokens` | filter tokens that write the domain under `k__`, as configs up to v2.4.0 did | `write d__Bacteria` |
| `gyrB_sintax` | `method: sintax` on a marker without a SINTAX reference | `method: sintax is not available` |
| `shotgun_tax_lev_D` | Bracken at the domain rank | `tax_lev: D is not allowed` |

## The stand-in databases

A dry run never opens a database, but MetaFlux's start-up checks do look at a few
files, and the stubs give them exactly what they ask for:

- **Kraken2** (`stubs/kraken2_db/`). At start-up MetaFlux checks that the folder
  exists and reads `taxo.k2d`, Kraken2's own copy of the taxonomy, to learn which taxa
  the database ranks as a domain, a kingdom or a virus realm (the filter tokens are
  checked against them, and the OTU-table step fills its `d__` and `k__` slots from
  them). `taxo.k2d` here is a real file in kraken2-build's layout, holding a cut-down
  taxonomy of about forty taxa, written by `tests/make_test_kraken2_db.py`:

  ```bash
  python tests/make_test_kraken2_db.py .test/stubs/kraken2_db
  ```

  `hash.k2d` and `opts.k2d` are empty files: MetaFlux reads the size of `hash.k2d`
  only to size the Kraken2 job's memory when memory mapping is off (the default is
  on), and a dry run never asks Kraken2 to load them. Kraken2 cannot run on this
  folder, which is the point: it is a few kilobytes instead of tens of gigabytes. A
  test (`tests/test_kraken2_taxonomy_file.py`) checks that the committed `taxo.k2d` is
  what the script writes, so the two cannot drift apart.
- **sylph** (`stubs/sylph/`). MetaFlux checks that each database file exists, ends in
  `.syldb` or `.syl2db`, reads the release and the sketch density from its name
  (`gtdb-r232-c200-dbv1.syldb` means GTDB r232 at c 200) and takes its size for the
  memory request. The stubs are small text files under the real file names; sylph
  never opens them in a dry run.
- **sylph-tax lineage tables** (`stubs/sylph_taxonomy/`). For a built-in taxonomy name
  (`GTDB_r232`, `FungiRefSeq-latest`) MetaFlux checks that the table sylph-tax would
  read is present in `references.sylph_taxonomy_dir`. The stubs are tiny gzip files
  under the real table names.
- **Host genome** (`stubs/host/`). A local path in `host_genomes` is only used when
  the index-building job runs, so any file would do for a dry run; the stub is a
  made-up 300 bp FASTA.

## Adding a case

Copy the closest config, change its header comment, `out_dir` and the keys the case
needs, and run `bash .test/run_dry_runs.sh`. For a case that must stop at start-up,
put it under `config/expected_stop/` with a `# EXPECT:` line quoting a fragment of
the stop message, as the three there do. Keep the data as it is: a few hundred read
pairs per sample are enough for Snakemake to find the samples, and the repository
should stay small.
