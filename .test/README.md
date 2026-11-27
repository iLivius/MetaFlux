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
| `data/<marker>/` | The first 300 read pairs of two samples of each public test dataset (16S, ITS, 18S, gyrB, rpoB; for COI, two slices of 300 pairs of the one mock run) and of two wastewater samples for shotgun mode, as `{sample}_R1/_R2.fastq.gz`. Enough for Snakemake to find the samples; far too few for any result to mean anything. |
| `data/16S_digit_names/` | Links to the two 16S samples of `data/16S/`, named `01` and `002`, for the expected stop `16S_digit_sample_names`. |
| `primers/<marker>/` | The forward and reverse primer FASTA of each marker, as used on the full test datasets. |
| `stubs/kraken2_db/` | A stand-in Kraken2 database folder (see below). |
| `stubs/sylph/`, `stubs/sylph_taxonomy/` | Stand-in sylph database files and sylph-tax lineage tables (see below). |
| `stubs/host/host_standin.fa.gz` | A made-up 300 bp "host genome", so that the host-removal rules appear in one shotgun plan. |
| `custom_error_functions/loess_wrapper.R` | A minimal error function of your own (experimental `error_model: custom`): it hands DADA2's transition counts to `loessErrfun` and returns its table unchanged. Used by `16S_custom_error_function` and shown as the example on the documentation page of the error model. |
| `custom_error_functions/binned_from_counts.R` | An error function of your own with a fit of its own: DADA2's binned fit, with every quality score of the counts at which all four true bases have bases as a bin. It carries an unchanged copy of `makeBinnedQualErrfun` (DADA2 1.38.0, under DADA2's license), which DADA2 1.34.0 lacks. Shown as an example on the documentation page of the error model; no test or dry run uses it. |
| `custom_error_functions/monotone_gam.R` | An error function of your own built on a monotone GAM (R package scam): per kind of misread, an error rate that may only fall as the quality score rises, with fall-backs for thin counts. Needs scam (`workflow/envs/dada2_custom.yaml`). Shown as an example on the documentation page of the error model; no test or dry run uses it. |
| `custom_error_functions/custom_function_checks.R` | The R checks of `error_model: custom` on made-up tables (base R only): each broken table must stop with its reason, attributes are dropped, a rate below 1e-7 warns, the random-number state and kind are put back, and the round keeper writes each warning once. Run by `tests/test_custom_error_function_r.py` when `METAFLUX_RSCRIPT` names an Rscript; not part of the dry runs. |
| `r_checks/error_model_checks.R` | R checks of the error-model scripts on made-up input: `learn_errors` sorts the filtered files in the same order in every locale, `dada_seqtab` pairs each R1 file with the R2 file of the same name and stops when one is missing (this case needs the R packages pacman and dada2), and the texts of the flat-rounds log line, of the binned fit's stop on quality scores outside its anchors and of the report's titles with many anchors. Run by `tests/test_error_model_r.py` when `METAFLUX_RSCRIPT` names an Rscript; not part of the dry runs. |
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
| `COI_sintax` | COI | sintax | – | the downloads of MIDORI2 and of the bacterial COI homologues, the MIDORI2 conversion, the COI screens, `orientation: mixed` |
| `16S_error_model_settings` | 16S | sintax | 144 | the settings of DADA2's error learning written as text, in capitals or as a single code (`max_consist: "200"`, `max_flat_rounds: 0`, `error_model: " LOESS "`, `accept_unconverged: "R2:0A55BE94"`), which the start-up checks accept and normalise |
| `16S_custom_error_function` | 16S | sintax | 144 | the experimental `error_model: custom` with `custom_error_function` naming `custom_error_functions/loess_wrapper.R` by a relative path: the start-up check accepts the block, and the plan is built (that the file is an input of `learn_errors` and `error_model_report`, and that both run in `workflow/envs/dada2_custom.yaml`, is checked on the rule text by `tests/test_error_model_settings.py`) |

gyrB and rpoB have no SINTAX reference, so `method: sintax` is refused for them at
start-up, and COI has no rdp trainset, so `method: rdp` is refused for it; those
refusals are expected-stop cases rather than matrix cells.

Shotgun, one config per profiler, and one for the `auto` floor of the k-mer evidence
gate:

| Config | Profiler | Also covers |
|---|---|---|
| `shotgun_kraken2` | Kraken2 + Bracken | the k-mer evidence gate, the taxon filter, read extraction |
| `shotgun_sylph` | sylph | two databases with built-in sylph-tax taxonomies |
| `shotgun_both` | both | host removal with the stand-in host genome, the taxon filter, read extraction |
| `shotgun_kmer_auto` | Kraken2 + Bracken | the k-mer evidence gate with `min_distinct_minimizers: auto`, at `tax_lev: G`, so the start-up warning for a gate below species is printed for `auto` |

Expected stops (`config/expected_stop/`):

| Config | What must stop the run | Expected message fragment |
|---|---|---|
| `16S_old_tokens` | filter tokens that write the domain under `k__`, as configs up to v2.4.0 did | `write d__Bacteria` |
| `gyrB_sintax` | `method: sintax` on a marker without a SINTAX reference | `method: sintax is not available` |
| `COI_rdp` | `method: rdp` on a marker without an rdp trainset | `method: rdp is not available` |
| `shotgun_tax_lev_D` | Bracken at the domain rank | `tax_lev: D is not allowed` |
| `shotgun_kmer_floor_word` | a k-mer evidence floor that is neither a whole number nor `auto` | `min_distinct_minimizers must be a whole number >= 0 (the default is 300) or the word auto` |
| `16S_max_consist_zero` | a cap of 0 rounds on DADA2's error learning, which would hand DADA2's rough starting estimate to the denoising | `max_consist must be a whole number of at least 1` |
| `16S_accept_unconverged_true` | `accept_unconverged: true`, an on/off switch where the setting takes the codes of accepted error models | `accept_unconverged takes the codes` |
| `16S_max_flat_rounds_fraction` | a fraction (`2.5`) as the number of flat rounds of the early stop | `max_flat_rounds must be a whole number >= 0` |
| `16S_error_model_word` | an instrument name (`novaseq`) as the error function of DADA2's error learning, which takes `auto`, `loess`, `binned`, `scam` or `custom` | `error_model must be one of auto, loess, binned, scam` |
| `16S_custom_error_function_missing_file` | `error_model: custom` with an error-function file that does not exist (a typo in its name) | `is not an existing file` |
| `16S_custom_error_function_bad_name` | `error_model: custom` with a function name R cannot give a function (`loess-wrapper`, with a dash) | `custom_error_function.name must be the name of an R function` |
| `16S_learn_errors_block_empty` | a `learn_errors:` block whose only line (`accept_unconverged`) is commented out, which YAML reads as empty and which would replace the template's whole block | `amplicon.dada2.learn_errors is empty in your run config` |
| `16S_two_amplicon_blocks` | a second `amplicon:` block pasted at the end of the run config, which YAML would keep instead of the first, losing the marker, the primers and the taxonomy filter without a word | `the key 'amplicon' is written twice in the same block (lines 18 and 29)` |
| `ITS_trim_final_cycle_word` | a word other than `auto`, `true` or `false` (`always`) for the removal of the extra last sequencing cycle in primer trimming | `trim_final_cycle must be auto, true or false` |
| `16S_digit_sample_names` | amplicon samples named `01` and `002`, which R and pandas read back as numbers | `mode: amplicon cannot use these sample names: ['002', '01']` |

Settings typed inline with `--config` (after the two groups above, in `run_dry_runs.sh`):

| Check | What is typed inline | What must happen |
|---|---|---|
| `phylogeny_off`, `phylogeny_on` | `amplicon={phylogeny: {enabled: false}}` over a config that turns the phylogeny on, and `true` over one that leaves it off | the plan holds no phylogeny job, and phylogeny jobs, as the same values written in a config file would give |
| `accept_code` | `amplicon={dada2: {learn_errors: {accept_unconverged: [R2:0a55be94]}}}` | the dry run plans the run: a model code typed inline is accepted at start-up (how a code is normalised is checked by `tests/test_error_model_settings.py`, and how inline text becomes a list by `tests/test_inline_config.py`) |
| `max_consist_zero` | `amplicon={dada2: {learn_errors: {max_consist: 0}}}` | the run stops at start-up with `max_consist must be a whole number of at least 1`, as in a config file |

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
the stop message, as the ones there do. Keep the data as it is: a few hundred read
pairs per sample are enough for Snakemake to find the samples, and the repository
should stay small.
