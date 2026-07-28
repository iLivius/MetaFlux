# Configuration

Everything MetaFlux does is decided by one YAML file. There are no command-line
switches for analysis parameters — the run is fully described by its config, which
makes it easy to archive next to the results and re-run later.

The shipped template is `config/config.yaml`. It is a template: the `/path/to/...`
entries are placeholders that must be replaced before anything runs.

## How the file is loaded

The Snakefile hard-codes `configfile: "config/config.yaml"`, so that file is
**always** read as a base layer. A file passed with `--configfile` is merged on top
of it, and `--config key=value` on the command line overrides both. A common pattern
is to copy the template to `config/config.local.yaml` (git-ignored) and run with
`--configfile config/config.local.yaml`, keeping the shipped template untouched.

```bash
snakemake --sdm conda --cores 16 --configfile config/config.local.yaml --dry-run
```

!!! warning "Why the taxon filter lists ship commented out"

    Because the template is always loaded underneath the run config, any concrete
    value left in it silently applies to runs that do not set their own. That is
    harmless for a thread count and destructive for
    `amplicon.taxonomy.filter.keep`: a 16S `keep: [k__Bacteria, k__Archaea]`
    leaking into an 18S run would discard every eukaryotic ASV. The `keep` and
    `discard` lines are therefore commented out in the template, so unset they
    resolve to `[]` — no filtering, which fails safe.

Paths may be absolute or relative. **Relative paths resolve against the directory
snakemake is invoked from**, not against the location of the config file. Absolute
paths are safest, especially on a cluster.

Each section below is tagged with the mode that consumes it: `[shared]`,
`[amplicon]`, or `[shotgun]`. Keys belonging to the mode that is not running are
ignored entirely.

## Which file a setting actually comes from

Four kinds of file get mentioned across this site, and only three of them decide
anything. It is worth being explicit about which is which, because they are easy to
confuse once you are reading the documentation and the config side by side.

| File | Where | Read by | Can it change a run? |
|---|---|---|---|
| These documentation pages | `docs/**.md` | people | **No** |
| The config template | `config/config.yaml` | Snakemake | Yes — always, as the base layer |
| A run config | any `.yaml` given to `--configfile` | Snakemake | Yes — merged on top |
| A marker pack | `workflow/markers/<type>.yaml` | Snakemake | Yes — but in its own namespace |

!!! warning "Documentation pages set nothing"

    Everything on this site is prose rendered from Markdown, including the per-marker
    pages. **Nothing under `docs/` is ever read by Snakemake.** The YAML blocks shown
    on these pages are illustrations of what to put in your own config — they are not
    in force anywhere, and they cannot override or be overridden by anything.

    So the question "if two pages disagree, which one wins?" has no winner: neither
    does. A disagreement between two documentation pages is simply a documentation
    bug. The file that decides is `config/config.yaml` plus whatever you merge on top
    of it.

### The config layers

Three layers, each beating the one before it:

1. **`config/config.yaml`** — the Snakefile hard-codes `configfile: "config/config.yaml"`.
   That path is resolved against the directory snakemake **runs in**, not against the
   Snakefile, so it loads on every run launched from the repository root — the supported
   way to invoke MetaFlux. Started from somewhere else, Snakemake either aborts
   (`Workflow defines configfile config/config.yaml but it is not present`) or, if you
   passed `--configfile`, silently uses that file *alone* with no template underneath.
2. **`--configfile my_run.yaml`** — merged on top, key by key, at any depth.
3. **`--config key=value`** — command line, beats both. Only **top-level** keys can be
   named this way: `--config mode=shotgun` works, `--config amplicon.dada2.merge.min_overlap=99`
   is rejected outright. Nested settings belong in a `--configfile`.

The merge in step 2 is nested rather than wholesale, which is what makes short run
configs practical. A run config naming two settings changes exactly those two:

```yaml
# my_run.yaml — an overlay that changes two settings
amplicon:
  dada2:
    merge:
      min_overlap: 10
output:
  out_dir: /data/my_run_out
```

Running with `--configfile my_run.yaml` gives `min_overlap: 10` and the new output
directory, while `max_mismatch`, `just_concatenate`, every `trunc_len` key and
everything else still come from the template. You never have to restate a block to
change one value inside it.

This is an **overlay, not a complete run config**: it inherits the template's
`/path/to/...` placeholders for `input.fastq_dir` and `amplicon.primers`, which still
have to be replaced before anything runs.

!!! warning "The merge recurses only through blocks"

    Key-by-key merging happens as long as what you write at each level is itself a
    block. Write a key with a list, a scalar, or **no value at all**, and it replaces
    that whole subtree in the template instead of merging into it.

    The dangerous form is the last one, because it looks like an unfinished thought
    rather than an instruction:

    ```yaml
    amplicon:
      trunc_len:        # <- no value: YAML reads this as null
    ```

    That silently deletes all eight `trunc_len` keys the template provides, and the run
    fails later with a `KeyError` from deep inside the workflow rather than a config
    error. Same for `max_ee: [1]`, which replaces `[2, 5]` outright rather than merging
    element-wise. If you do not intend to change a block, leave it out entirely.

!!! tip "The flip side, and why `keep`/`discard` ship commented out"

    Because the template is *always* underneath, a concrete value left in it silently
    applies to any run that does not set its own. Harmless for a thread count; not
    harmless for `amplicon.taxonomy.filter.keep`, where a 16S `[k__Bacteria, k__Archaea]`
    left in the template would discard every ASV of an 18S run. That is why those two
    lines ship commented out — see the warning above.

### Marker packs live beside the config, not inside it

A marker pack holds *facts about the marker*: which database the classifier uses, which
reference the amplicon-length probe measures, whether a region extractor exists, how many
slots the marker's rank model has. A run config holds *choices about this
experiment*: primers, expected length, truncation, taxon filter, paths.

The two are kept in genuinely separate namespaces. A pack is read into `MARKERS[type]`
and never merged into `config`, so a setting cannot silently arrive from both places.
This is enforced by construction rather than convention: across all five shipped packs,
**no pack key appears anywhere under `amplicon:` in the config template**. The only name
the two share at all is `references`, handled below.

If you genuinely need a different pack — another database release, a different rank
model — put a file of the same name in `config/markers/`. That is a whole-file
replacement, not a per-key merge: your file supersedes the shipped pack entirely, so it
must declare everything the pack needs, not just the part you wanted to change.

### Where a pack and the config do interact

Two directions, and they are not symmetrical.

**A config entry overriding a pack** happens in exactly one place: database locations.
Every pack declares a default filename and download URL for each database it owns, and
the run config only has to name one when you want to **relocate** it — to point at a
copy already on a shared filesystem rather than let MetaFlux download its own:

```yaml
references:
  silva:
    train: /shared/dbs/silva_nr99_v138.2_toGenus_trainset.fa.gz
```

A **non-empty** `references.*` entry wins; the pack default is used otherwise. The test
is truthiness, so an entry left `null` or blank counts as absent and falls back to the
pack silently rather than erroring. The key that relocates each database is listed in
[Markers](../amplicon/markers/index.md).

**A pack overriding the config** is the other direction, and it is the one case where
the workflow does not do what your config literally said. It is worth knowing about,
because if you copy a config between markers you will meet it.

It happens when a setting asks for something that does not exist *for the marker you
chose* — not something unwise, something genuinely unavailable. Two cases ship today,
and they behave differently on purpose:

| You asked for | On these markers | What happens |
|---|---|---|
| `extraction.enabled: true` | 18S, gyrB, rpoB | **warning, the run continues** with extraction off |
| `taxonomy.method: sintax` | gyrB, rpoB | **error, the run stops** |

#### Region extraction

Metaxa2 cuts the ribosomal gene out of a 16S ASV; ITSx cuts ITS1 or ITS2 out of a
fungal one. gyrB and rpoB are protein-coding genes, so neither tool applies to them,
and the 18S pack declares no extractor either. All three say `extractor: none`, and the
workflow builds no extraction step at all for them.

If your config still says `extraction.enabled: true`, the pipeline would go looking for
`seqs_extracted.fasta` — a file no rule produces — and Snakemake would stop with a
missing-input error naming a filename you never asked for. Instead MetaFlux turns the
setting off itself and tells you on stderr:

```text
[MetaFlux] warning: marker 18S has no target-region extractor;
ignoring amplicon.extraction.enabled: true (no extraction step will run)
```

Note this one is a **warning** — the run proceeds. If you adapted a 16S config for an
18S run, extraction quietly not happening is the correct outcome, not a fault to chase.

#### SINTAX classification

VSEARCH's `--sintax` needs a database in its own format. SILVA, UNITE and PR2 all have
one, so 16S, ITS and 18S can be classified either way. gyrB's DD7RZ8 trainset and
rpoB's FROGS build ship an RDP-style training set only, so those two packs leave
`taxonomy_sintax_db` empty. Asking for sintax there stops the run at once, naming the
alternative:

```text
[MetaFlux] marker gyrB has no SINTAX reference database, so
amplicon.taxonomy.method: sintax is not available for it. Use method: rdp.
```

This one errors rather than falling back, and the asymmetry with the case above is
deliberate. Switching off a step that could not have run costs you nothing you could
otherwise have had. Silently classifying with a different algorithm than the one you
asked for would change your results without telling you — so it refuses instead.

Neither case is really the pack "winning an argument". Nothing is being weighed against
anything: the config asked for a tool that does not exist for that marker, and the pack
is simply where that fact happens to be written down.

## The four lines that define a run

Most of the file is tuning that can stay at its defaults. Four settings describe
what the analysis actually is:

| Setting | Choices |
|---|---|
| `mode` | `amplicon` or `shotgun` |
| `amplicon.type` | `16S`, `ITS`, `18S`, `gyrB`, `rpoB` (amplicon mode only) |
| `input.fastq_dir` | where the paired FASTQs are |
| `output.out_dir` | where results go |

---

## Mode selection `[shared]`

```yaml
mode: amplicon          # amplicon | shotgun
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `mode` | Chooses which set of rules is built into the workflow graph | `amplicon` | Case-insensitive. Anything other than `amplicon` or `shotgun` stops the run at parse time. |

The mode is chosen here, never auto-detected from the data. Amplicon mode requires
primer FASTAs; shotgun mode supports fetching reads from an NCBI BioProject and
removing an arbitrary host genome. See
[choosing a mode](../getting-started/choosing-a-mode.md).

!!! warning "The mode must match the data"

    Both modes accept `{sample}_R1/_R2` filenames, so pointing amplicon rules at
    shotgun reads does not necessarily error — it just runs the wrong analysis.
    The mode can be flipped for a single run without editing the file:
    `--config mode=shotgun`.

---

## Input `[shared]`

```yaml
input:
  fastq_dir: /path/to/fastq_dir
  bioproject: null        # e.g. PRJNA603575  — shotgun only
  accession_list: null    # e.g. [SRR11605259, SRR11605260]
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `fastq_dir` | Directory holding the paired FASTQs for every sample | `/path/to/fastq_dir` (placeholder) | Created if it does not exist, because BioProject mode downloads into it. Sample names come from globbing this directory. |
| `bioproject` | NCBI BioProject accession to fetch from SRA instead of using local files | `null` | Shotgun only. `null` means use whatever is already in `fastq_dir`. |
| `accession_list` | Restricts a BioProject fetch to named runs | `null` | A YAML list of run accessions. Useful when a BioProject mixes WGS with amplicon or RNA-Seq runs. |

**File naming.** Amplicon mode requires `{sample}_R1.fastq.gz` /
`{sample}_R2.fastq.gz` and stops with a clear message if it finds SRA-style
`{sample}_1/_2` naming instead. Shotgun mode accepts either convention and detects
which one is in use. The reads themselves must be `.fastq.gz` or the uncompressed
`.fastq`: `fq` and `fq.gz` pass the startup check but are not resolved by the
rules afterwards. The extension check only looks at the R1 files — a mixture
there is rejected rather than half-processed, while a mismatched R2 extension is
not caught, so the whole directory is best kept consistent.

Sample names are taken from the part of the filename before `_R1`/`_1`. Names
containing any of `* # @ % ^ / ! ? & : ; | < >` or a space are rejected at parse
time, since Snakemake treats those characters specially. Underscores are fine.

**BioProject mode.** When `bioproject` is set, MetaFlux queries NCBI e-utils,
keeps only runs whose library layout is `PAIRED`, and caches the run table as
`{bioproject}_runinfo.csv` inside `fastq_dir` so later invocations work offline.
The `download_sra` rule then fetches each run with `prefetch` plus `fasterq-dump`
and renames SRA's `_1`/`_2` output to `_R1`/`_R2`.

---

## Output `[shared]`

```yaml
output:
  out_dir: /path/to/output_dir
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `out_dir` | Root directory for every result, log and benchmark of the run | `/path/to/output_dir` (placeholder) | Logs land in `out_dir/logs/`, benchmark timings in `out_dir/benchmarks/`. |

The directory layout inside `out_dir` differs by mode and is described in
[output files](output.md).

---

## References

```yaml
references:
  refdb_root: refdb
  phix:
    fasta: refdb/phix/phix.fna
    fetch_url: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz
  kraken_db: /path/to/kraken2_db
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `refdb_root` | Cache directory for the marker reference databases MetaFlux downloads itself | `refdb` | Amplicon only. Shared across runs and across output directories, so a database is fetched once and reused. |
| `phix.fasta` | Local path of the PhiX genome | `refdb/phix/phix.fna` | Amplicon mode also builds a bowtie2 index next to it (same basename, `.bt2` files). |
| `phix.fetch_url` | Where the PhiX genome is downloaded from if missing | NCBI FTP, `GCF_000819615.1_ViralProj14015` | Fetched once by the `fetch_phix` rule and gunzipped. |
| `kraken_db` | Directory of a prebuilt Kraken2/Bracken index | `/path/to/kraken2_db` (placeholder) | Shotgun only, user-supplied — MetaFlux never builds or downloads it. The run stops immediately if the path does not exist. |

**PhiX is used by both modes, but differently.** In amplicon mode the FASTA above is
fetched and indexed, and bowtie2 maps reads against it. In shotgun mode BBDuk uses
the `phix174_ill.ref.fa.gz` reference bundled inside the BBTools conda package, so
neither the download nor the index is needed there.

### Marker databases are declared by the marker packs, not here

SILVA, UNITE, PR2, SILVA-Eukaryotic, DD7RZ8 and FROGS are **not** listed in the
config. Each marker pack (`workflow/markers/<type>.yaml`) declares its own database
filenames and download URLs, and MetaFlux fetches (and where needed converts) them
on first use into `refdb_root`. An entry under `references` is only needed to
**relocate** a database — for instance to point at a copy already on a shared
filesystem:

```yaml
references:
  silva:
    train: /shared/dbs/silva_nr99_v138.2_toGenus_trainset.fa.gz
```

When such a key is present it wins; otherwise the path is `refdb_root` joined with
the pack's declared filename. The full set of relocatable keys:

| Config key | Marker | Default location under `refdb_root` |
|---|---|---|
| `silva.train` | 16S | `silva/silva_nr99_v138.2_toGenus_trainset.fa.gz` |
| `silva.species` | 16S | `silva/silva_v138.2_assignSpecies.fa.gz` |
| `silva.sintax` | 16S | `silva/silva_nr99_v138.2_sintax.fa.gz` (derived from the trainset, not downloaded) |
| `unite.fasta` | ITS | `unite/sh_general_release_dynamic_s_all_19.02.2025.fasta` |
| `unite.sintax` | ITS | `unite/utax_reference_dataset_all_19.02.2025.fasta.gz` |
| `unite.uchime_its1` | ITS | `unite/uchime_its1.fasta` |
| `unite.uchime_its2` | ITS | `unite/uchime_its2.fasta` |
| `silva_euk.fasta` | 18S | `silva_euk/silva_132.18s.dada2.fa.gz` |
| `pr2.dada2` | 18S | `pr2/pr2_SSU_dada2.fasta.gz` |
| `pr2.utax` | 18S | `pr2/pr2_SSU_UTAX.fasta.gz` |
| `gyrb.dada2` | gyrB | `gyrb/train_set_gyrB_v6.fa.gz` |
| `rpob.frogs` | rpoB | `rpob/rpob_refseq_cc_20240707_dada2.fa.gz` (the converted DADA2 trainset) |

A reference that has neither a config entry nor a pack default stops the run with a
message naming both places it could be fixed.

---

## Amplicon parameters `[amplicon]`

Ignored entirely when `mode: shotgun`. The one exception: if `amplicon.primers.fwd`
or `.rev` is filled in during a shotgun run, MetaFlux prints a warning that they are
being ignored — so a leftover value is visible rather than silently misleading.

### Marker, region and seed

```yaml
amplicon:
  type: 16S             # 16S | ITS | 18S | gyrB | rpoB
  its_region: ITS2      # ITS1 | ITS2 (ignored unless type == ITS)
  seed: 42
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `type` | Selects the marker pack: which databases are used, whether a region extractor runs, how the amplicon length is measured, and which taxonomic ranks the output carries | `16S` | Matched case-insensitively, so `gyrb`, `GYRB` and `gyrB` all resolve to the shipped `gyrB` pack. An unknown value lists the available markers and stops. |
| `its_region` | Which ITS subregion was amplified | `ITS2` | Only read when `type: ITS`. Selects the UNITE UCHIME probe reference and the ITSx region kept after extraction. |
| `seed` | Fixes the random number generator for the stochastic steps of the amplicon path | `42` | Passed to `dada_seqtab` and `assign_taxonomy`. |

Each marker has its own page with the databases it uses and its quirks:
[16S](../amplicon/markers/16S.md), [ITS](../amplicon/markers/its.md),
[18S](../amplicon/markers/18S.md), [gyrB](../amplicon/markers/gyrb.md),
[rpoB](../amplicon/markers/rpob.md).

!!! note "How far the seed goes"

    DADA2's `learnErrors(randomize = TRUE)` and `assignTaxonomy`'s bootstrap
    confidence are fully reproducible with this seed at any thread count —
    verified by running the same data twice at 8 threads and getting identical
    output. VSEARCH's `--sintax` is not: it races several threads on one random
    number stream, so a handful of per-rank confidence values can still drift
    between runs. For byte-identical `sintax` output, set
    `resources.threads.assign_taxonomy: 1`, which is slower but single-streamed.

### PhiX removal

```yaml
  decontamination:
    remove_phix: true
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `remove_phix` | Maps reads against the PhiX genome with bowtie2 and keeps the pairs that do not align | `true` | Set `false` to skip. The `2.no_phix/` stage, its read count and the `nophix` Falco QC stage then do not exist, and Cutadapt reads the raw FASTQs directly. |

PhiX is spiked into most Illumina runs as a low-diversity balancing control; those
reads are not part of the sample and would be carried into ASV inference.

### Primers

```yaml
  primers:
    fwd: /path/to/fwd_primer.fasta
    rev: /path/to/rev_primer.fasta
    orientation: fixed      # fixed | mixed
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `fwd` | FASTA with the forward primer sequence | placeholder | Must exist — amplicon mode stops at parse time if either file is missing. |
| `rev` | FASTA with the reverse primer sequence | placeholder | Reverse complements of both are generated by the `revcomp_primers` rule; nothing needs to be supplied pre-complemented. |
| `orientation` | Whether R1 always carries the forward primer (`fixed`) or reads may come in either orientation (`mixed`) | `fixed` | `mixed` runs a second Cutadapt pass with the primers swapped, reverse-complements its output so both passes agree, and concatenates the two. Anything other than `fixed`/`mixed` stops the run. |

The primer pair also drives the probe cache. The cache key is a short hash of both
primer files, so changing a primer invalidates the cached amplicon-length
measurement automatically, and reusing the same primers across runs reuses it.

### Amplicon length

```yaml
  expected_length: auto     # int | [min, max] | "auto"
  probe_length_stat:
    16S: p95
    18S: p95
    rpoB: p95
  min_overlap: 12
```

`expected_length` is the length of the amplicon the primers produce. It is not used
to filter anything directly — it is the number that decides whether the chosen R1
and R2 truncation points still leave the pair able to merge.

| Key | What it does | Default | Notes |
|---|---|---|---|
| `expected_length` | Amplicon size used for the truncation overlap check | `auto` | `auto` measures it from the reference database; an integer states it directly; a `[min, max]` pair states a window and the **maximum** is used for the overlap check. |
| `probe_length_stat.<key>` | Which statistic of the measured length distribution `auto` should take | `p95` for `16S`, `18S` and `rpoB` | One of `min`, `q1`, `median`, `q3`, `p95`, `p99`, `max`. |
| `min_overlap` | Bases of forward/reverse overlap the truncation picker must preserve | `12` | A design constraint only. The overlap actually enforced when pairs are merged is `dada2.merge.min_overlap`. |

**What `auto` does per marker.** Markers whose reference is a full-length gene are
measured by in-silico PCR: the primers are matched against the reference and the
lengths of the surviving products form the distribution. Markers whose reference is
already cut to amplicon length are measured directly.

| Marker | How `auto` measures the amplicon | Probe reference | Statistic key read |
|---|---|---|---|
| 16S | In-silico PCR | SILVA trainset | `16S` |
| ITS | Direct read of pre-extracted subregion lengths | UNITE UCHIME ITS1 or ITS2 | not used |
| 18S | In-silico PCR | SILVA-Eukaryotic 18S | `18S` |
| gyrB | Direct read — DD7RZ8 already ships trimmed amplicons | DD7RZ8 | `16S` (deliberate reuse) |
| rpoB | In-silico PCR — FROGS ships full-length genes | FROGS RefSeq rpoB | `rpoB` |

!!! warning "Keep all three `probe_length_stat` keys"

    For every marker except ITS, MetaFlux checks at parse time that the key its
    pack names is present in `probe_length_stat`, whatever `expected_length` is set
    to, and exits with the missing key name if it is not. gyrB reads the `16S`
    key — that is a deliberate reuse, not a typo. Deleting any of the three
    entries breaks the markers that point at it.

!!! warning "18S V9 must be sized manually"

    V4 primer pairs work with `auto`: the measured window is wide (roughly
    380–595 bp of core amplicon, depending on the pair). V9 (1391F/EukBr) does not,
    because EukBr sits at the 3′ terminus and most references are truncated before
    it — only 24–41% are recovered in silico, so `auto` would size the window from
    a biased minority. Set it manually, for example `expected_length: [115, 150]`.

An integer must be between 50 and 10000 bp; a `[min, max]` pair must have exactly
two positive numbers with `min <= max`. Out-of-range or malformed values are
rejected with an explicit message rather than quietly used.

Details of the measurement and how it feeds truncation are on
[amplicon length and truncation](../amplicon/length-and-truncation.md).

### Primer trimming (Cutadapt)

```yaml
  cutadapt:
    max_error_rate: 0.2
    min_length: 50
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `max_error_rate` | Cutadapt's `-e` — the fraction of the matched primer length that may differ | `0.2` | Cutadapt reads IUPAC codes in the primer natively, so a degenerate primer costs nothing from this budget; it covers sequencing errors in the first cycles and genuine primer/template mismatch, which is common with universal primers. Cutadapt's own default is 0.1 — 0.2 is deliberately more permissive because the 5′ pass runs with `--discard-untrimmed`, so an unmatched primer costs the whole pair. Check the per-sample match rates in `stats/cutadapt/` rather than assuming. The same value is reused by the in-silico PCR probe. |
| `min_length` | Reads shorter than this after trimming are discarded, together with their mate | `50` | Cutadapt's `--minimum-length`, with `--pair-filter=any`. |

No quality trimming happens at this stage. DADA2's error model is fitted from the
whole run's quality scores, and pre-trimming by quality would remove the
low-quality bases it needs in order to learn the quality-to-error relationship.

### Truncation

```yaml
  trunc_len:
    mode: auto                  # auto | manual
    q_threshold: 20
    q_floor: 15
    resolve_policy: raise_trunc # raise_trunc | relax_q | error
    subsample_files: 12
    manual_r1:                  # e.g. 260
    manual_r2:                  # e.g. 240
```

`truncLen` is one length per read direction; everything past it is discarded. It is
pulled two ways: cutting shorter removes the low-quality 3′ tail the error model
dislikes, but cutting too short stops R1 and R2 from reaching each other, and a pair
that cannot merge is lost entirely.

| Key | What it does | Default | Notes |
|---|---|---|---|
| `mode` | `auto` derives the cuts from observed quality; `manual` uses the two numbers below | `auto` | Manual mode skips the quality analysis completely. |
| `q_threshold` | Quality below which `auto` cuts | `20` | The picker reads the lower quartile (Q1) of the per-base quality at each cycle from Falco's report, takes the median across samples, and cuts just before the first position where that value falls below the threshold. |
| `q_floor` | Lowest threshold `relax_q` may descend to | `15` | Only used by the `relax_q` policy. |
| `resolve_policy` | What to do when the quality-based cuts leave too little overlap | `raise_trunc` | See below. |
| `subsample_files` | Carried in the template but not read by any rule in v2.1.0 | `12` | The picker uses every sample's Falco report, not a subsample. |
| `manual_r1` / `manual_r2` | Explicit R1 and R2 truncation lengths | unset | Required when `mode: manual` — the run fails if either is missing. Ignored for ITS. |
| `min_read_coverage_pct` | `auto` only: a cut is never placed past the length this percentage of a sample's reads still reach | `95` | Read from Falco's "Sequence Length Distribution". See the note below — this is the guard against `auto` discarding an entire library. |

!!! note "What `min_read_coverage_pct` is protecting against"

    Falco reports a quality bin for every position up to the **longest** read in a
    file, however few reads reach it. `filterAndTrim` **discards** any read shorter
    than `truncLen` rather than padding it, so a cut placed where only a handful of
    reads reach does not trim a few bases — it throws away the library. On a run whose
    quality never drops below `q_threshold`, nothing else stops the cut landing there.

    Before quality is looked at, each sample's bins are therefore cut off at the length
    `min_read_coverage_pct` of its reads still reach. On the 18S
    [test dataset](../about/test-datasets.md) this moves the R1 ceiling from 288 bp to
    271 bp, which is the difference between retaining 9 read pairs out of 908,768 and
    retaining 90.9% of the library.

    Mind the direction: the value is a guarantee about reads, so raising it asks for a
    position more reads reach and therefore gives a **shorter** cut with fewer reads
    discarded. Lowering it gives a longer cut and discards more. Between about 80 and 99
    it barely moves, because read lengths cluster. `100` is a cliff rather than an off
    switch — it demands every read reach the cut, collapsing the ceiling onto the single
    shortest read in the file. No value disables the guard.
    `logs/pick_trunclen.log` prints each sample's longest read next to the length the
    requested percentage reaches, and emits a `NOTE` when the two differ by 10 bp or
    more.

The overlap constraint is `truncR1 + truncR2 >= expected_length + min_overlap`.
When the quality-based cuts break it:

| Policy | Behaviour |
|---|---|
| `raise_trunc` | Extends both cuts past the quality drop to recover the missing bases, splitting the deficit between R1 and R2 in proportion to how much read length each has left. Errors if even the combined headroom is not enough. |
| `relax_q` | Lowers the quality threshold by 1 at a time, down to `q_floor`, re-cutting until the overlap holds. Errors if the floor is reached without success. |
| `error` | Stops immediately and reports the shortfall in bp. |

The decision, including which policy resolved it and how much slack is left, is
written to `stats/trunclen.json`.

!!! note "ITS ignores truncLen"

    ITS amplicon length varies genuinely between taxa, so a fixed cut would throw
    away real short amplicons. For ITS the picker still computes and records the
    quality-based cuts, but `dada_filter` overrides `truncLen` to `c(0, 0)`. The
    3′ tail is not left untrimmed: `dada2.trunc_q` still trims each read
    adaptively, with `dada2.max_ee` and the per-read minimum length as the quality
    gate. This override also beats `mode: manual`, so `manual_r1`/`manual_r2` have
    no effect for ITS.

!!! warning "The `auto` picker can overshoot on very clean runs"

    If quality never drops below `q_threshold`, the cut is placed at the last
    quality bin — which can sit above the bulk of the trimmed read lengths and
    drop every read at the filtering step. See
    [troubleshooting](../troubleshooting.md).

### ASV length filter

```yaml
  length_filter:
    mode: auto          # auto | manual
    window_margin: 50
    range: null         # [min, max]
```

This filter runs on the ASVs, **after** target extraction, so the lengths it sees are
directly comparable to the probe distribution the window was built from.

| Key | What it does | Default | Notes |
|---|---|---|---|
| `mode` | `auto` builds the window from the probe; `manual` requires an explicit `range` | `auto` | |
| `window_margin` | Bases added below and above the probe window | `50` | The auto window is `[probe_q1 − margin, probe_p95 + margin]`. |
| `range` | Explicit `[min, max]` in bp | `null` | When set it wins, whatever `mode` says. |

When `mode: auto` runs without a probe (because `expected_length` was given
manually), the window falls back to the resolved expected length ±15%.

The filter writes `stats/dada2/asv_length_stats.json` — pre-extraction,
filter-source and kept distributions side by side — plus a histogram with the window
drawn on it. If the window removes every ASV the run stops rather than writing an
empty table.

!!! warning "ITS auto needs an auto expected length"

    `type: ITS` with `length_filter.mode: auto` also requires
    `expected_length: 'auto'`, since ITS never resolves an expected length through
    the truncation picker and the fallback would have nothing to size from.
    MetaFlux catches this combination at parse time, before the expensive DADA2
    run, and the message offers two fixes: set `expected_length: 'auto'`, or
    switch to `length_filter.mode: manual` with an explicit `range`.

### Target-region extraction

```yaml
  extraction:
    enabled: true
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `enabled` | Runs the marker's region extractor over the ASVs, trimming each to the target region and dropping ASVs the extractor cannot place | `true` | Metaxa2 for 16S, ITSx for ITS. 18S, gyrB and rpoB have no extractor. |

Markers with no extractor (`18S`, `gyrB`, `rpoB`) should be run with
`enabled: false`. Left `true`, MetaFlux prints a warning and forces it off rather
than asking the workflow graph for a file no rule produces. The ASV length filter
still runs either way.

### Taxonomy

```yaml
  taxonomy:
    method: sintax        # rdp | sintax
    min_boot: 80          # rdp only
    sintax_cutoff: 0.8    # sintax only
    try_rc: true          # rdp only
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `method` | Which classifier assigns taxonomy | `sintax` | `rdp` is DADA2's `assignTaxonomy` (plus `addSpecies` for 16S); `sintax` is `vsearch --sintax`, which uses roughly 2–3 GB regardless of ASV count. Both paths apply the same taxon filter and write the same three tables. |
| `min_boot` | RDP bootstrap confidence a rank must reach to be reported | `80` | Read only by the `rdp` path. Range 0–100. |
| `sintax_cutoff` | SINTAX confidence a rank must reach to be reported | `0.8` | Read only by the `sintax` path. Range 0–1; 0.8 is roughly equivalent to `min_boot: 80`. |
| `try_rc` | Also tries the reverse complement of each ASV during classification | `true` | Read only by the `rdp` path. The `sintax` path always searches both strands. |

!!! warning "gyrB and rpoB are rdp-only"

    The DD7RZ8 and FROGS releases ship a DADA2 trainset with no SINTAX build, so
    those packs declare no SINTAX reference. `method: sintax` for either marker is
    refused at startup with a message naming the marker, rather than failing later
    on a missing path.

### Taxon filter

```yaml
    filter:
      enabled: true
      keep:    #[k__Bacteria, k__Archaea]
      discard: #[o__Chloroplast, f__Mitochondria]
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `enabled` | Master switch | `true` | `false` keeps every ASV regardless of the lists. |
| `keep` | Tokens defining what survives | unset → `[]` | An ASV is kept only if one of its ranks matches a token. An empty list does nothing. |
| `discard` | Tokens pruning what `keep` left | unset → `[]` | An ASV is removed if any of its ranks matches a token. Runs after `keep`. |

Each token is matched against a whole segment of the taxonomy string
(`k__Bacteria`, `o__Chloroplast`, `g__Ralstonia`), so it targets exactly one rank
and never matches loosely somewhere else in the lineage. There is no built-in
default: whatever is listed is exactly what is applied.

| Marker | Suggested `keep` | Suggested `discard` |
|---|---|---|
| 16S | `[k__Bacteria, k__Archaea]` | `[o__Chloroplast, f__Mitochondria]` |
| ITS | `[k__Fungi]` | `[]` |
| 18S | `[d__Eukaryota]` (PR2 ranks; `k__Eukaryota` with SILVA-Euk) — note this also drops ASVs left unclassified at Domain, which on some 18S datasets is a large fraction; see [18S](../amplicon/markers/18S.md) before applying it | `[]` |
| gyrB | no recommended default — the DD7RZ8 paralog tag is a separate, structural filter | |
| rpoB | no recommended default — FROGS RefSeq has no equivalent contamination to remove | |

Contaminants a blank extraction control reveals go in `discard`, e.g.
`discard: [o__Chloroplast, f__Mitochondria, g__Ralstonia, g__Bradyrhizobium]`. A
bare string is accepted where a list is expected — `keep: k__Fungi` is treated as
`keep: [k__Fungi]` — so a single token cannot be split into characters by accident.
The full treatment is on
[keeping and discarding taxa](../amplicon/taxon-filter.md).

!!! note "The old regex keys are refused"

    Configs written before the token lists used `include_pattern` /
    `exclude_pattern` regexes. Those keys now stop the run with instructions,
    rather than being ignored and leaving the filter silently disabled.

### DADA2 — tier 1, operational parameters

```yaml
  dada2:
    pool: false             # false | "pseudo" | true
    max_ee: [2, 5]
    trunc_q: 2
    filter:
      min_len: 50
      min_len_stat: null
      max_len: null
    learn_errors:
      nbases: 1e8
      max_consist: 999
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `pool` | How samples are pooled during ASV inference | `false` | `false` treats each sample independently (fastest); `"pseudo"` runs a second pass informed by ASVs seen across samples; `true` pools all samples (most sensitive to rare variants, most expensive). Strings are read case-insensitively. |
| `max_ee` | `filterAndTrim`'s `maxEE` — the expected-error ceiling per read, as `[R1, R2]` | `[2, 5]` | R2 is looser because reverse reads degrade faster on Illumina chemistry. |
| `trunc_q` | `filterAndTrim`'s `truncQ` — trims each read at the first base at or below this quality | `2` | The DADA2 default. This is the adaptive 3′ trim, distinct from the fixed `truncLen` cut. It is applied **before** `truncLen`; see the note below. |
| `filter.min_len` | Per-read length floor after trimming | `20` | DADA2's own `filterAndTrim` default. Used for every marker; for ITS it is the fallback when `min_len_stat` is unset. |
| `filter.min_len_stat` | Optionally derives `min_len` from the probe distribution instead | `null` | **ITS only, opt-in.** Left `null`, every marker uses `filter.min_len`. Setting a statistic (e.g. `q1`) puts the floor near 175 bp, above ~24% of the UNITE reference — see [ITS](../amplicon/markers/its.md). The other markers never take this path: their probe measures the full amplicon, which exceeds read length. |
| `filter.max_len` | Per-read length ceiling | `null` | `null` means no upper limit. |

!!! note "`truncQ` runs before `truncLen`, and the leftovers are discarded rather than kept short"

    This is DADA2's own ordering inside `filterAndTrim`, not something MetaFlux changes,
    but it produces a failure that is easy to misdiagnose. `truncQ` trims each read at its
    first low-quality base; `truncLen` is then applied, and any read that `truncQ` has
    already shortened **below** `truncLen` is **dropped**, not passed through at its
    shorter length.

    So a run that loses an unexpected share of reads at `dada_filter` may be losing them
    to `truncQ`, while the number that looks wrong is `truncLen`. Before lowering
    `truncLen`, check the filtering columns in `stats/read_tracking.txt` and consider
    whether a lower `trunc_q` is the real cause — particularly on a run whose 3′ quality
    degrades early.
| `learn_errors.nbases` | Bases sampled to train the error model | `1e8` | The DADA2 default. Reads are drawn in random order, seeded by `amplicon.seed`. |
| `learn_errors.max_consist` | Maximum self-consistency iterations of the error model | `999` | Effectively unlimited; the model normally converges long before. |

### DADA2 — tier 2, advanced tuning

```yaml
    dada:
      omega_a: 1e-40
    merge:
      min_overlap: 12
      max_mismatch: 0
      just_concatenate: false
      trim_overhang: true
    chimera:
      method: consensus
      min_fold_parent: 4
      allow_one_off: false
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `dada.omega_a` | `OMEGA_A` — how strong the evidence must be before a sequence is called a new ASV rather than an error of an existing one | `1e-40` | The DADA2 default. A lower value is stricter; a higher one yields more ASVs. |
| `merge.min_overlap` | `mergePairs` minimum overlap, in bp, enforced at runtime | `12` | DADA2's default. This is the one that actually decides whether a pair merges; `amplicon.min_overlap` is the design constraint used earlier when choosing truncation lengths. |
| `merge.max_mismatch` | Mismatches tolerated inside the overlap | `0` | DADA2's default, and deliberate: after denoising, a mismatch in the overlap is evidence the pair should not merge rather than a sequencing error. |

!!! note "What these two cost in practice"

    Both were previously set looser than DADA2 ships them (10 and 2). Moving them to
    the package defaults was measured by running each of the five
    [test datasets](../about/test-datasets.md) twice, once with each pair of values:

    | Marker | Reads merged | ASVs |
    |---|---|---|
    | 16S | −1.4% | 273 → 211 |
    | ITS | **−33.5%** | 29 → 12 |
    | 18S | −0.6% | 1,508 → 1,268 |
    | gyrB | −2.7% | 256 → 241 |
    | rpoB | −0.3% | 48 → 44 |

    Read counts through filtering and denoising were identical in every case; the whole
    difference appears at the merge step. ASV counts drop much further than reads
    because pairs that merge only with mismatches tolerated tend to yield distinct,
    low-abundance sequences — so the loss is concentrated in rare ASVs. Loosening
    them again recovers those ASVs, at the cost of accepting overlap disagreements
    that DADA2 treats as evidence against merging.

    **ITS is the exception worth knowing about.** ITS2 length varies widely between
    fungi, so with 2×300 reads the longer amplicons overlap by only a little — exactly
    the pairs a 10 → 12 bp requirement drops, and their short overlap sits in the
    low-quality 3′ tails where `max_mismatch: 0` then removes more. A third of the
    merged reads went. On the mock that cost no genera, but for ITS on a more variable
    community `min_overlap: 10` with `max_mismatch: 2` is a defensible setting.
| `merge.just_concatenate` | Joins R1 and R2 with a run of Ns instead of merging on overlap | `false` | A safety valve for ITS amplicons too long to overlap at all. The resulting sequences are not real contigs and the length filter must be set accordingly. |
| `merge.trim_overhang` | Trims bases that read past the far end of the amplicon | `true` | Safe for 16S; necessary for short ITS amplicons where the read runs off the end into primer sequence. |
| `chimera.method` | `removeBimeraDenovo` method | `consensus` | Chimeras are called per sample and the verdicts pooled. |
| `chimera.min_fold_parent` | How much more abundant a parent must be than its suspected chimera | `4` | DADA2's `minFoldParentOverAbundance`, raised above the package default. Higher is **more permissive**: parents must be more abundant to qualify, so fewer sequences are called chimeric and more ASVs survive. See the note below. |
| `chimera.allow_one_off` | Also flags sequences one mismatch away from a perfect chimera | `false` | A relaxed mode; more aggressive, and more likely to remove real variants. |

!!! tip "Where to look before changing the chimera settings"

    `min_fold_parent: 4` is a deliberate departure from the DADA2 default, in the
    direction of removing fewer sequences. The reasoning, and the range of values the
    DADA2 authors and community have discussed for pooled and amplicon-specific data
    (typically 4 to 8), is set out in the upstream discussions — worth reading before
    changing it, since the right value depends on the marker and on how samples were
    pooled:

    - [dada2 #2030 — Questions about chimera filtering](https://github.com/benjjneb/dada2/issues/2030)
    - [dada2 #1715 — High loss of reads following chimera removal](https://github.com/benjjneb/dada2/issues/1715)
    - [dada2 #1368 — Remove chimeras, consensus vs. pooled](https://github.com/benjjneb/dada2/issues/1368)

    A large chimeric fraction is more often a symptom than a parameter problem: check
    that primers were fully removed first, since residual primer sequence is a common
    cause of apparent chimeras.

---

## Shotgun parameters `[shotgun]`

Ignored entirely when `mode: amplicon`.

### Decontamination

```yaml
shotgun:
  decontamination:
    remove_phix: true
    host_genomes: []
    host_min_id: 0.95
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `remove_phix` | BBDuk k-mer matching against the bundled `phix174_ill.ref.fa.gz`, at `k=31 hdist=1` | `true` | No download or index build is involved. Keep it on unless the lab is known to have skipped the spike-in. |
| `host_genomes` | YAML list of local paths and/or URLs to host FASTA(.gz) references | `[]` | `[]` skips host removal. Multiple entries are concatenated and indexed once, so one BBMap index serves every sample. |
| `host_min_id` | BBMap `minid` — the alignment identity above which a read is called host-derived | `0.95` | Reads that do not map are kept. |

Host removal is alignment-based rather than k-mer based, which is more sensitive to
host reads that diverge from the reference assembly — the approach the BBTools
authors recommend.

!!! warning "Use a masked reference for human"

    An unmasked human genome removes real microbial reads that happen to be
    homologous to ribosomal or low-complexity human sequence. The masked reference
    on Zenodo at `https://zenodo.org/records/4116107` (Handley 2020,
    `human_virus_masked.fasta.gz`, 889 MB) is a masked hg19 — ribosomal RNA,
    plant, animal and fungal homology and low-entropy sequence removed — with a
    further layer of viral masking added. The
    first set of masks is what stops real microbial reads being discarded as host;
    the viral layer additionally keeps human regions homologous to viruses out of
    the filter, which is the right behaviour when the viral fraction is of
    interest. Because the file is roughly 900 MB, download it once to a local path
    and point `host_genomes` at it rather than re-fetching the URL on every fresh
    run.

!!! warning "BBDuk and BBMap both need correctly paired files"

    Both tools read R1 and R2 in lockstep and require the two files to hold the
    same number of records in the same order. Read length is not the issue: ragged
    reads, and mates of one pair with different lengths, are processed normally. A
    crash reading `java.lang.AssertionError: List size mismatch: N vs M` (BBDuk) or
    `There appear to be different numbers of reads in the paired input files`
    (BBMap) means the pairing is broken — usually because reads were filtered
    before deposition and orphaned mates were left in one file. Repair the FASTQs
    with the tool shipped in the same conda package rather than skipping the step.

    ```bash
    repair.sh in1=R1.fastq.gz in2=R2.fastq.gz \
        out1=fixed_R1.fastq.gz out2=fixed_R2.fastq.gz outs=singletons.fastq.gz
    ```

    See [troubleshooting](../troubleshooting.md) and
    [decontamination](../shotgun/decontamination.md).

### Adapter and quality trimming

```yaml
  fastp:
    min_read_length: 100
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `min_read_length` | fastp's `--length_required` — reads shorter than this after trimming are dropped | `100` | Comfortable for 2×150 bp reads; use around 70 for HiSeq 2×100. |

fastp runs with `--detect_adapter_for_pe`, plus `--cut_front` and `--cut_right`
sliding-window trimming from both ends. Its JSON report is read later by Bracken to
pick a read-length model, and by the read-tracking table.

### Kraken2

```yaml
  kraken:
    confidence: 0.15
    hit_groups: 3
    memory_mapping: true
    threads: null
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `confidence` | Fraction of a read's k-mers that must fall inside the clade of the taxon it is assigned to | `0.15` | Higher is stricter. Kraken2 does not simply reject a read that falls short — it moves the label up the taxonomy until the score clears the threshold, and calls the read unclassified only if even root fails. Raising it therefore mostly pushes assignments to coarser ranks, leaving more for Bracken to redistribute, and only secondarily increases the unclassified fraction. Wright et al. (2023) report two optima on simulated data: 0.60 by mean F1, and 0.15 by L1 distance and by recall/precision over classified reads. The default of 0.15 favours accurate relative abundances; raise it towards 0.60 to favour F1. |
| `hit_groups` | `--minimum-hit-groups` — distinct k-mer hit groups a read needs to keep its classification | `3` | |
| `memory_mapping` | Pages the database from disk on demand instead of loading it into RAM | `true` | See below. |
| `threads` | Pins Kraken2's thread count | `null` | `null` falls back to `resources.threads.kraken2`. |

Kraken2 is also run with `--use-names`, `--paired`, `--gzip-compressed`,
`--report-minimizer-data` and `--classified-out`, so classified read pairs are kept
alongside the report.

!!! note "Memory mapping decides how much RAM a Kraken2 job asks for"

    Kraken2's memory request is computed at parse time, not taken from the
    `resources.mem_mb.kraken2` entry, which is not read.

    - **`memory_mapping: true`** — the kernel pages the database on demand and
      shares those pages between concurrent Kraken2 processes through the OS file
      cache, so running several samples in parallel does not multiply RAM use.
      Each job is budgeted a fixed 20000 MB of private workspace.
    - **`memory_mapping: false`** — every process loads the whole hash privately.
      The budget becomes the size of `hash.k2d` plus 10% headroom, so it
      self-adjusts to the database in use (a 16 GB PlusPF, a ~100 GB PlusPF, a
      200+ GB core_nt). If `hash.k2d` cannot be read yet — while the database is
      still downloading, say — a conservative 110000 MB is used instead. (With
      memory mapping on the fixed 20000 MB applies either way, since the index
      size is never consulted.)

### Bracken

```yaml
  bracken:
    tax_lev: S
    threshold: 10
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `tax_lev` | Rank at which abundance is re-estimated: `D`, `P`, `C`, `O`, `F`, `G` or `S` | `S` | Drop to `G` for a shallow database, where species-level calls are not supported by the reference content. |
| `threshold` | Minimum Kraken2 read count a taxon needs at `tax_lev` to enter Bracken's estimate | `10` | Applied to Kraken2's own count for the taxon, before Bracken re-estimates anything. Taxa below it are removed from the model and their reads discarded rather than redistributed (Bracken reports these as "reads discarded"). It is a floor on direct evidence, not on the final estimate: a species whose reads Kraken2 mostly left at genus level — common where the database holds several close relatives — can be dropped despite being abundant. Lower it if that is a concern, and expect a longer tail of one-read taxa. |

Bracken redistributes the reads Kraken2 assigned to internal nodes down to the
target rank. Bracken databases ship k-mer distributions for several read lengths
(50, 75, 100, 150, 200, 250, 300); MetaFlux reads the post-trimming mean read length
out of the fastp JSON and picks the closest one, because applying a 150 bp
distribution to 100 bp reads biases the redistribution.

### OTU-table taxon filter

```yaml
  taxonomy_filter:
    enabled: false
    keep:    []
    discard: []
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `enabled` | Master switch | `false` | Off by default: metagenomic sample types vary too much for a universal default to be safe. |
| `keep` | Tokens defining what survives | `[]` | e.g. `[k__Bacteria, k__Archaea]` for prokaryotes only. |
| `discard` | Tokens pruning what `keep` left | `[]` | e.g. `[g__Homo]` as an explicit host safety net. |

This is the shotgun counterpart of `amplicon.taxonomy.filter` and uses the same
rank-aware token matching: `keep` first, then `discard`, empty list means that
direction does nothing. It is applied to the final OTU table after Bracken, so it
changes which OTUs appear in the table and never how reads were classified.

### Per-taxon read extraction

```yaml
  extract_taxa: []
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `extract_taxa` | NCBI scientific names whose reads should be pulled out into their own FASTQ files | `[]` | `[]` skips extraction entirely. Example: `[Bacteria, Archaea]`. |

Each name is looked up in the Bracken report at whatever rank it appears, so a broad
group (`Bacteria`) and a single species (`Escherichia coli`) work equally well.
KrakenTools then extracts every read under that taxid, including its children.
Output files are named after a slug of the taxon (lower-case, spaces to
underscores). A taxon absent from a sample yields empty files rather than a failure.
See [per-taxon read extraction](../shotgun/read-extraction.md).

---

## Resources `[shared]`

```yaml
resources:
  threads_default: 4
  mem_mb_default: 2000
  threads:
    kraken2: 16
    # … one entry per rule
  mem_mb:
    kraken2: 20000
    # … one entry per rule
```

Per-rule CPU and RAM live in one block with global fallbacks. Rules look their
values up by name; a name that is absent falls back to `threads_default` or
`mem_mb_default`. These values feed straight into Snakemake's resource accounting,
so they are what a cluster executor turns into job requests.

| Key | What it does | Default |
|---|---|---|
| `threads_default` | Threads for any rule with no `threads` entry | `4` |
| `mem_mb_default` | Memory in MB for any rule with no `mem_mb` entry | `2000` |

### Threads

| Rule name | Mode | Default | Notes |
|---|---|---|---|
| `multiqc` | shared | 2 | |
| `falco` | amplicon | 2 | |
| `bowtie2` | amplicon | 8 | Covers both `build_phix_index` and the per-sample `rm_phix` mapping. |
| `cutadapt` | amplicon | 4 | Covers both `trim_primers` and the in-silico PCR probe. |
| `pick_trunclen` | amplicon | 4 | |
| `dada_quality_plots` | amplicon | 4 | |
| `dada_filter` | amplicon | 4 | |
| `dada_seqtab` | amplicon | 8 | |
| `dada_length_filter` | amplicon | 2 | |
| `target_extract` | amplicon | 8 | Metaxa2 or ITSx. |
| `assign_taxonomy` | amplicon | 8 | Set to 1 for byte-reproducible `sintax` output. |
| `aggregate_read_counts` | shared | 2 | Present in the config template but not read — the rule declares no `threads`; see the note below. |
| `decontam_phix` | shotgun | 6 | BBDuk scales poorly past 4–6 worker threads on a 5 kb reference; run more samples in parallel instead. |
| `build_host_index` | shotgun | 16 | |
| `decontam_host` | shotgun | 16 | |
| `fastp` | shotgun | 16 | Read by the `trim_adapters` rule. fastp scales poorly past about 16 threads. |
| `kraken2` | shotgun | 16 | Overridden by `shotgun.kraken.threads` when that is not `null`. |
| `bracken` | shotgun | 4 | |
| `kraken_biom` | shotgun | 2 | |
| `extract_taxon_reads` | shotgun | 4 | |
| `compress_extracted` | shotgun | 6 | pigz parallelism for the extracted read files. |
| `download_sra` | shotgun | 6 | `fasterq-dump` parallelism; past about 8 the bottleneck is disk, not CPU. |

### Memory

| Rule name | Mode | Default (MB) | Notes |
|---|---|---|---|
| `multiqc` | shared | 2000 | |
| `bowtie2` | amplicon | 4000 | |
| `pick_trunclen` | amplicon | 4000 | |
| `dada_quality_plots` | amplicon | 4000 | |
| `dada_filter` | amplicon | 2000 | |
| `dada_seqtab` | amplicon | 16000 | The heaviest amplicon step: error learning, denoising, merging and chimera removal all happen here across every sample at once. |
| `dada_length_filter` | amplicon | 4000 | |
| `target_extract` | amplicon | 8000 | |
| `assign_taxonomy` | amplicon | 16000 | The `rdp` path is the demanding one; `sintax` stays around 2–3 GB whatever the ASV count. |
| `aggregate_read_counts` | shared | 2000 | Present in the config template but not read — the rule declares no `mem_mb`; see the note below. |
| `decontam_phix` | shotgun | 4000 | |
| `build_host_index` | shotgun | 30000 | BBMap index build for the host reference; a masked human genome needs roughly 24–28 GB. |
| `decontam_host` | shotgun | 30000 | Mapping against that index. |
| `fastp` | shotgun | 8000 | Read by the `trim_adapters` rule. |
| `kraken2` | shotgun | 20000 | Not read — the real figure is computed at parse time from `memory_mapping` and the actual `hash.k2d` size. |
| `bracken` | shotgun | 4000 | |
| `kraken_biom` | shotgun | 2000 | |
| `extract_taxon_reads` | shotgun | 2000 | |

!!! note "Names here are lookup labels, not always rule names"

    Most entries match a rule name exactly. A few are lookup labels rather than
    rule names: `bowtie2` covers `build_phix_index` and `rm_phix`, `cutadapt`
    covers `trim_primers` and `amplicon_probe`, and `fastp` is the label used by
    the `trim_adapters` rule — there is no rule called `fastp`, and a
    `trim_adapters` key here is ignored. Rules that declare no `threads` or
    `mem_mb` at all — `link_reads`, `count_reads_*`, `revcomp_primers`, the
    `fetch_*` and `convert_*` reference rules, `samples_manifest` and
    `aggregate_read_counts` — take Snakemake's own defaults and ignore any entry
    of the same name. Adding a key for a rule that does not look it up has no
    effect. Conversely, `finalize_otu_table` does look its resources up but has
    no entry in the template, so it runs on `threads_default` and
    `mem_mb_default`.
