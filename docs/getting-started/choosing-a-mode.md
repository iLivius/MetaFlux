# Choosing a mode

MetaFlux runs one of two pipelines. The choice is not a matter of taste — it follows
from what was in the tube when the library was made.

- **`amplicon`** — a marker gene was PCR-amplified with a specific primer pair, and
  those amplicons were sequenced. 16S rRNA for bacteria and archaea, ITS for fungi,
  18S rRNA for eukaryotes, or the protein-coding markers gyrB and rpoB.
- **`shotgun`** — total DNA was fragmented and sequenced without any target-specific
  amplification, so the reads come from everywhere in every genome present.

If the library prep involved primers targeting a gene, it is amplicon data. If it
did not, it is shotgun data.

## What each mode does with the reads

| | `amplicon` | `shotgun` |
|---|---|---|
| What was sequenced | One marker region, PCR-amplified | Whole DNA, no amplification |
| Core method | DADA2 exact ASV inference | Kraken2 classification + Bracken re-estimation |
| Reference | Marker database, fetched automatically (SILVA, UNITE, PR2, DD7RZ8, FROGS) | A Kraken2 index, supplied by the user |
| Primers | Required — `amplicon.primers.fwd` / `.rev` | Not used |
| Main result | `6.taxonomy/asv_table.txt` — ASVs × samples, with taxonomy | `03.abundance/otu_table.tsv` — taxa × samples, plus the BIOM version |
| Resolution | Exact sequence variants of one gene | Whatever the index supports, down to species with Bracken |
| Input naming | `{sample}_R1` / `{sample}_R2` only | `{sample}_R1/_R2` or SRA-style `{sample}_1/_2` |
| Can fetch its own reads | No | Yes — from an NCBI BioProject accession |

Both modes read the same `input` and `output` blocks, honour the same per-rule
`resources` settings, and finish by writing a MultiQC report and a per-sample
read-tracking table under `out_dir`.

## The switch

One key, at the top of `config/config.yaml`:

```yaml
mode: amplicon          # amplicon | shotgun
```

That key decides which set of rules is loaded at all: the Snakefile includes
`workflow/rules/amplicon/*.smk` or `workflow/rules/shotgun/*.smk` accordingly, plus
the shared rules in both cases. It can also be overridden for a single run without
touching the file:

```bash
snakemake --sdm conda --cores 16 --config mode=shotgun
```

## The other mode's settings are ignored

`config/config.yaml` carries both halves at all times, and each section is labelled
`[shared]`, `[amplicon]` or `[shotgun]`. The half that does not belong to the
selected mode is simply never read: an `amplicon:` block full of primer paths and
DADA2 parameters has no effect on a shotgun run, and a `shotgun:` block full of
Kraken2 and Bracken parameters — plus the `references.kraken_db` path, which sits
in the shared `references:` block — has no effect on an amplicon run.

There is no need to comment anything out. One config file can hold a working setup
for both and be switched with `mode` alone.

The one place this is not silent is primers. Shotgun mode prints a note if
`amplicon.primers.fwd` or `.rev` carries any value at all, then continues:

```text
[MetaFlux] warning: mode=shotgun but amplicon.primers fields are set — they will be ignored
```

The shipped template fills both fields with placeholder paths, and
`config/config.yaml` is always loaded as a base layer under the `--configfile` given
on the command line, so the note appears on most shotgun runs even when the run
config never mentions primers. It can be ignored, or silenced by setting both to
`null` in the run config.

!!! warning "The mode must match the data"

    Nothing is auto-detected. Both modes accept files named `{sample}_R1` /
    `{sample}_R2`, so pointing shotgun mode at amplicon reads — or the reverse —
    does not necessarily fail. It can run to the end and produce a table that looks
    ordinary and means nothing. The mode is a statement about how the library was
    made, and it has to be true.

## What to prepare before running

**Amplicon** needs, on top of the shared input and output paths:

- `amplicon.type` — which marker (`16S`, `ITS`, `18S`, `gyrB`, `rpoB`). This selects
  the reference databases, whether a region extractor runs, and which taxonomy
  methods are available.
- `amplicon.primers.fwd` and `.rev` — the primer sequences as FASTA files. Both must
  exist or the run stops immediately. They are used twice: Cutadapt trims them off
  the reads, and — when `expected_length` is left at `auto` (the default), for markers
  whose reference is a full-length gene — the same pair is run in silico against that
  reference to work out how long the amplicon should be.
- `amplicon.taxonomy.filter` — the keep/discard lists appropriate to the marker.
  There is no built-in default, so an unset filter keeps everything.

Reference databases come down automatically on first use; nothing else to install.
Per-marker detail is in the [markers section](../amplicon/markers/index.md).

**Shotgun** needs:

- `references.kraken_db` — the directory of a Kraken2 index. It is not fetched
  automatically and the run stops at once if the path does not exist. See
  [reference databases](../shotgun/databases.md).
- Optionally `shotgun.decontamination.host_genomes` — one or more host references
  for read removal, e.g. a masked human genome. See
  [decontamination](../shotgun/decontamination.md).
- Optionally `input.bioproject` instead of local FASTQs — MetaFlux resolves the
  PAIRED runs of an NCBI BioProject accession and downloads them into `fastq_dir`.

## Next

- [Installation](installation.md) — databases, environments, what is fetched and what is not.
- [Quick start](quick-start.md) — a first run, start to finish.
- [Amplicon overview](../amplicon/overview.md) and [shotgun overview](../shotgun/overview.md) — what each pipeline does, step by step.
