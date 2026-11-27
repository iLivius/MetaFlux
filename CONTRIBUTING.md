# Contributing to MetaFlux

Thank you for taking the time. MetaFlux is maintained by one person next to a day job
in a research institute, and it is used by biologists who read the code before they
trust it. The rules below follow from that: every change must be easy to read, easy to
check, and honest about what it does to a finished run. Questions that are not bug
reports go to the [Discussions](https://github.com/iLivius/MetaFlux/discussions); see
also [SUPPORT.md](SUPPORT.md).

## Reporting a bug

Open an issue with the **Bug report** form. It asks for what is needed to reproduce the
problem without a round of questions:

- the MetaFlux version (the banner the run prints, or `git describe --tags`), the mode,
  and for amplicon runs the marker, method and, for 16S, the SILVA release;
- your run config (the YAML you pass with `--configfile`; shorten paths if you like,
  but keep the keys);
- the command you ran, and what MetaFlux printed: the `[MetaFlux]` start-up lines and
  the error, or the whole dry-run output if the dry run is what fails;
- the log of the failing job, `logs/<rule>/<sample>.log` in the output folder;
- what you expected to happen instead.

Please check the [troubleshooting page](https://iLivius.github.io/MetaFlux/troubleshooting/)
first: the failures worth knowing about are described there with their fix.

## Proposing a change

- **Say what you have in mind before writing much code.** Open an issue (the **Feature
  request** form) or a discussion. Many things that look like a missing feature are
  deliberate (the docs usually say why), and a short exchange first saves both sides a
  rewrite.
- **One change per pull request**, branched from the development branch of the next
  version (`dev/vX.Y.Z`) if there is one, otherwise from `main`, with a title and
  description that say what changes for a user. Typos and small documentation fixes do
  not need an issue first.
- **The maintainer reviews and merges every pull request**, including those of people
  with write access to the repository. Open yours from a fork and keep the box that
  lets maintainers edit your branch ticked (GitHub ticks it by default), so that a
  small fix can be made on your branch directly instead of costing another round of
  review. With write access you can open it from a branch of this repository instead.
- **Say what runs again.** Snakemake keeps a record of every finished step in
  `.snakemake/` of the folder it runs from: the text of the step's command (comments
  included), its parameters, the list of its input files, and its conda environment
  file with that file's pin. A step that calls a script (`script:`) is judged by the
  script file's date, not its text: it runs again whenever the script is newer than
  the step's outputs. Changing any of them, even a comment inside a shell block, a
  script or an environment file, makes that step and everything after it run again in
  the finished output folders of every user who updates with `git pull` (a fresh clone
  reruns every script step anyway). That is often right (a correctness fix must reach
  old folders), but it has to be said in the pull request, and a change that only
  reorders or reformats a rule should be avoided for that reason.
- **No new dependency without a discussion first.** Every tool must be installable from
  bioconda or conda-forge, pinned where the version changes results, and listed in one
  of the environment files under `workflow/envs/`. After editing an environment file,
  refresh its pin file (`workflow/envs/<tool>.linux-64.pin.txt`; the installation page
  says how): as long as the old pin exists, Snakemake installs the old builds from it
  and never sees the edit. A tool used only by the tests or the checks (pytest, mkdocs)
  is not a pipeline dependency and is pinned in `.github/workflows/ci.yml` instead.
- **One pin is rebuilt from another, not solved: `dada2_custom.linux-64.pin.txt`.**
  `workflow/envs/dada2_custom.yaml` is `dada2.yaml` plus `r-mgcv` and `r-scam`, for error
  functions of the user's own, and its pin must hold every build of
  `dada2.linux-64.pin.txt`, so that a model learned in it rests on the same dada2, R and
  maths libraries as one learned in `dada2.yaml`. A fresh solve does not keep them: a
  dry-run solve of the yaml, when this pin was made, changed 97 of the 185 builds of
  dada2's pin, OpenBLAS and the GCC runtime among them. So whenever `dada2.yaml` or its
  pin changes, refresh dada2's pin first, make the same edit in `dada2_custom.yaml`, and
  rebuild its pin from dada2's, in a new folder outside the repository:

  ```bash
  prefix=/path/to/a/new/folder
  # 1. Build the environment of dada2's pin exactly: an explicit list is installed as
  #    it is, without solving.
  conda create --yes --no-default-packages -p "$prefix" --file workflow/envs/dada2.linux-64.pin.txt
  # 2. Install the two packages with every installed package frozen. conda updates
  #    openssl and ca-certificates on every install (its aggressive_update_packages
  #    setting), even with --freeze-installed, so name both with the version and build
  #    they have in dada2's pin.
  conda install --yes --freeze-installed --override-channels -c conda-forge -c bioconda \
      --strict-channel-priority -p "$prefix" r-mgcv r-scam \
      "openssl=<version>=<build>" "ca-certificates=<version>=<build>"
  # 3. Export it in the format of every pin file.
  conda list --explicit --md5 -p "$prefix" > workflow/envs/dada2_custom.linux-64.pin.txt
  ```

  The new pin must hold every line of dada2's pin, and the added lines must be the two
  packages and what they bring (nlme comes with mgcv):
  `comm -23 <(sort workflow/envs/dada2.linux-64.pin.txt) <(sort workflow/envs/dada2_custom.linux-64.pin.txt)`
  prints nothing, and `comm -13` with the same two files lists the added lines.
  `tests/test_error_model_settings.py` checks both, and that the yaml is `dada2.yaml`
  plus the two packages. The versions of mgcv and scam named under "R packages" in
  `docs/amplicon/error-model.md` are those of this pin; update them with it.
- **Run the checks** (next section) before you open the pull request, and add a test
  or a dry-run config when your change can be covered by one.

## Running the checks

The same checks run on GitHub for every push and pull request
(`.github/workflows/ci.yml`). Locally, from the repository root:

```bash
# 1. The unit tests of the pure-Python parts (start-up checks, reference check,
#    OTU-table builder, read-tracking parsers, taxon lookup, k-mer evidence gate,
#    MIDORI2 conversion and COI screens, orientation table and final-cycle decision
#    of primer trimming,
#    settings of DADA2's error learning, the choice of its error function, the list
#    of error functions and the vendored DADA2 function, the files and the MultiQC
#    section of the error-model report, and the check before denoising).
#    Needs pytest and PyYAML; the OTU-table tests also need biom-format and numpy
#    (the packages of workflow/envs/biom.yaml) and are skipped without them. The R
#    checks (the experimental error_model: custom, and the file order, file pairing and
#    texts of the error-model scripts) run only when METAFLUX_RSCRIPT names an Rscript
#    (see tests/helpers.py); the CI skips them.
pytest tests

# 2. The dry-run matrix: every marker, method, SILVA release and shotgun profiler
#    on the tiny data under .test/, plus the configs that must stop at start-up.
#    Needs Snakemake on PATH, nothing else; no conda environment is built.
bash .test/run_dry_runs.sh

# 3. The documentation. --strict turns every warning into an error, so a page left
#    out of the nav or a link to a renamed heading fails the build.
pip install mkdocs-material==9.7.7 mkdocs-glightbox==0.5.3
mkdocs build --strict

# 4. Snakemake's own linter, informational: it flags a few long-standing style points
#    that are kept on purpose, so do not try to silence them; read its output for
#    the rules you touched.
snakemake --lint --configfile .test/config/16S_sintax_silva144.yaml
snakemake --lint --configfile .test/config/shotgun_both.yaml
```

`tests/helpers.py` explains how the tests reach code that lives in rule files and in
scripts that expect Snakemake's `snakemake` object; `.test/README.md` explains the dry
runs, the stand-in databases and how to add a config. A real run on the test datasets
(`docs/about/test-datasets.md`) is still the final check for anything that changes
results; the CI cannot run tools.

Before a release, make one more real run: amplicon mode with a single sample. Copy the
read pair of one sample of a test dataset into a folder of its own, point
`input.fastq_dir` at it and `output.out_dir` at a new, empty folder; any marker will
do, since the step being checked is DADA2's. Never point `output.out_dir` at the
folder of a finished run: Snakemake would run that sample's steps again there and then
rebuild the run's tables from that one sample. It worked when the run finishes and
`stats/read_tracking.txt` has one row, named after the sample, with a number in every
column (an `NA` there means a step lost the sample's name). None of the checks above
can see this case: the unit tests cover the Python parts only, and the dry runs plan
the jobs without running them, so no R script runs. Yet one sample is exactly where
DADA2 behaves differently: given one sample, `derepFastq()`, `dada()` and
`mergePairs()` return that sample's result on its own instead of a list with one
result per sample, and `dada_seqtab` has to put it back into a list.

## How the code is written

MetaFlux's readers are biologists with some coding ability, and every line is read and
understood by a human before it goes in. The style rules are there to keep that
possible:

- **Comment every rule, function, script and block**, in plain language, about the
  biology and the data flow, not the syntax: what the step is for, what it takes in
  and where that comes from, what it does to it, what it writes, and what reads that
  next. A comment that turns out to be obvious costs nothing; a missing one costs an
  hour six months later. Write for a colleague who knows genomics, not software
  patterns; define a technical term in half a sentence the first time you use it.
- **Prefer the obvious solution.** Explicit intermediate variables with meaningful
  names, plain loops, and one step per line beat chained one-liners and clever
  abstractions. If a less obvious approach is genuinely better, use it and say why in
  a comment.
- **Avoid deep nesting and dense functional constructs.** Someone must be able to
  follow a script top to bottom and know what it does without running it.
- **Stop early with a message that says what to do.** A wrong config value should stop
  the run at start-up with the key, the value and the fix named, as the existing
  `sys.exit("[MetaFlux] ...")` checks in `workflow/rules/shared/00_common.smk` do, not
  surface as a crash three hours in.
- **Scripts are plain Python or R** that read the `snakemake` object and nothing else;
  MetaFlux scripts never import from each other (each runs in its own conda
  environment). Where two scripts share a function, they carry a copy and say so, and
  a test keeps the copies identical. One exception: the list of the error functions
  DADA2 can learn its error model with, `workflow/scripts/amplicon/60_error_functions.R`,
  is sourced by both `60d_learn_errors.R` and `60e_error_model_report.R` (both run in
  `dada2.yaml`, or in `dada2_custom.yaml`, which holds the same builds plus mgcv and scam,
  when a run names an error function of the user's own), so that the learning and the
  report always mean the same thing by a function's name, and both rules list it as an
  input. Its header says what a change to
  it costs (every amplicon folder learns its error models again) and how to add a
  function. Another exception: `workflow/scripts/amplicon/30a_final_cycle.py` is a small
  command-line program, not a Snakemake script. Rule `trim_primers` runs it from its
  shell command, once per read direction, to decide whether Cutadapt removes the last
  sequencing cycle, and lists it as an input; its header says what a change to it costs
  (primer trimming runs again in every amplicon folder).
- **Vendored code is never edited.** `workflow/scripts/amplicon/vendor/` holds one
  function copied unchanged from DADA2 (`makeBinnedQualErrfun`, under DADA2's license,
  `LICENSE.dada2`). A test checks its bytes against the md5 in its header, so a hand
  edit or an editor that strips trailing spaces fails the checks. A newer version of
  the function goes in as a new copy, with its source and md5 in the header.
- **Public text stays general and timeless.** The docs, the comments, the commit
  messages and the changelog name versions, tools and data releases, but no dates of
  our own work, commit hashes, issue numbers or private discussions. Facts come with
  how they were measured.
- **Commit messages are plain sentences**: a short first line saying what changes for
  the user, then why, and what runs again in a finished folder if anything does.

## Documentation

The site under `docs/` is built with MkDocs (`mkdocs serve` previews it). Every page
must be listed in the `nav` of `mkdocs.yml`, every link must resolve (the strict build
checks links to headings too), and a change in behaviour must reach every page that
describes it: the per-marker pages, the configuration reference, the output reference
and the troubleshooting page overlap on purpose. The template `config/config.yaml` is
documentation as well: its comments are what most users read first.

## License

MetaFlux is released under the [MIT License](LICENSE). By contributing you agree that
your contribution is released under the same license.
