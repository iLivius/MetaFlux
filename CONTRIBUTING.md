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
- **One change per pull request**, branched from `main`, with a title and description
  that say what changes for a user. Typos and small documentation fixes do not need an
  issue first.
- **Say what runs again.** Snakemake records the text of every rule's command, its
  parameters, its script and its conda environment file. Changing any of them, even a
  comment inside a shell block or an environment file, makes every finished output
  folder run that step and everything after it again. That is often right (a
  correctness fix must reach old folders), but it has to be said in the pull request,
  and a change that only reorders or reformats a rule should be avoided for that reason.
- **No new dependency without a discussion first.** Every tool must be installable from
  bioconda or conda-forge, pinned where the version changes results, and listed in one
  of the environment files under `workflow/envs/`. A tool used only by the tests or the
  checks (pytest, mkdocs) is not a pipeline dependency and is pinned in
  `.github/workflows/ci.yml` instead.
- **Run the checks** (next section) before you open the pull request, and add a test
  or a dry-run config when your change can be covered by one.

## Running the checks

The same checks run on GitHub for every push and pull request
(`.github/workflows/ci.yml`). Locally, from the repository root:

```bash
# 1. The unit tests of the pure-Python parts (start-up checks, reference check,
#    OTU-table builder, read-tracking parsers, taxon lookup). Needs pytest and
#    PyYAML; the OTU-table tests also need biom-format and numpy (the packages of
#    workflow/envs/biom.yaml) and are skipped without them.
pytest tests

# 2. The dry-run matrix: every marker, method, SILVA release and shotgun profiler
#    on the tiny data under .test/, plus the configs that must stop at start-up.
#    Needs Snakemake on PATH, nothing else; no conda environment is built.
bash .test/run_dry_runs.sh

# 3. The documentation. --strict turns every warning into an error, so a page left
#    out of the nav or a link to a renamed heading fails the build.
pip install mkdocs-material==9.7.7
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
  a test keeps the copies identical.
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
