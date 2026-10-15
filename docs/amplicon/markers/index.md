# Amplicon markers

MetaFlux runs one marker per amplicon run, chosen with `amplicon.type` in the config:
`16S`, `ITS`, `18S`, `gyrB` or `rpoB`. Everything that differs between markers — which
reference the amplicon-length probe measures, which database the classifier uses, whether
a target-region extractor exists, how many taxonomic ranks the output carries — is
declared in a small YAML file called a *marker pack*, one per marker, under
`workflow/markers/`. The rest of the amplicon pipeline is the same code for all five.

Marker packs hold facts about the marker, never run settings. Primers, expected lengths,
truncation values and taxon filters live in the run config; they are the things that change
from one experiment to the next. A pack is read into its own namespace and is never merged
into `config`, so no key has two homes.

!!! note

    A pack in `config/markers/` shadows a shipped pack of the same name. Adding a marker
    is therefore a data-only change — a new pack file, no code edit — as long as its
    databases are plain downloads: the fetch rule is generated from the `url` the pack
    gives each reference. A database that arrives as an archive, or whose headers have
    to be reshaped before a classifier will read them, still needs a rule of its own in
    `workflow/rules/amplicon/10_refdb.smk`. That is how UNITE, the UNITE UCHIME
    subregions, the SILVA 144 training set (whose kingdom names lose SILVA's `--other`
    suffix), the SILVA SINTAX builds, the 18S PR2 SINTAX file (whose domain moves to the
    letter `d:`) and the FROGS rpoB build are handled.

    A pack in `config/markers/` replaces the shipped pack completely, so after an
    upgrade copy the shipped pack again and re-apply your change: every pack must
    declare `reference_values` (see [the reference check](#the-reference-check-reference_values)),
    and a copy made for an older MetaFlux, which has none, stops the run at start-up,
    whichever marker the run uses.

## The five markers side by side

| Marker | What it is, what it resolves | `probe_mode` | Why that mode | Taxonomy reference | SINTAX build? | Extractor | Ranks |
|---|---|---|---|---|---|---|---|
| **16S** | Bacterial/archaeal SSU rRNA; the standard prokaryote survey marker, genus-level for most taxa | `pcr` | SILVA ships full-length SSU genes, so the amplicon window has to be cut out with the primers | SILVA 144 toGenus trainset, all seven SILVA ranks with the kingdom names written without `--other` (`silva_train`, built from `silva_train_raw`; SILVA 138.2 selectable), plus SILVA species assignment (`silva_species`) for DADA2 `addSpecies` | Yes — `silva_sintax`, derived from the trainset by `convert_silva_sintax`. `rdp` and `sintax` both available | Metaxa2 | 8 on SILVA 144 (Domain, Kingdom … Species), 7 on SILVA 138.2 |
| **ITS** | Fungal internal transcribed spacer, ITS1 or ITS2 (`amplicon.its_region`); the fungal barcode, species-level for well-represented lineages | `direct` | The UNITE UCHIME release is already the ITS1/ITS2 subregion plus a 15 bp anchor of flanking rRNA at each end, so lengths are read straight off it | UNITE general release `sh_general_release_dynamic_s_all_19.02.2025` (`unite_fasta`) | Yes — UNITE UTAX release (`unite_sintax`). `rdp` and `sintax` both available | ITSx | 7 |
| **18S** | Eukaryotic SSU rRNA; protists and other micro-eukaryotes | `pcr` | The probe substrate (SILVA-Euk v132) is full-length SSU | PR2 v5.1.1 SSU DADA2 file (`pr2_dada2`) | Yes — `pr2_sintax`, PR2 v5.1.1's UTAX file (`pr2_utax`) with its domain moved to the SINTAX letter `d:`, built by `convert_pr2_sintax`. `rdp` and `sintax` both available | none | 9 under `rdp`, 8 under `sintax` |
| **gyrB** | DNA gyrase subunit B; protein-coding, single-copy, higher intra-genus resolution than 16S | `direct` | The DD7RZ8 v6 trainset is already primer-trimmed in-silico amplicons (~247 bp), not whole genes | DD7RZ8 v6 trainset, GTDB r226 (`gyrb_dada2`) — the same file is the probe substrate | **No.** `rdp` only | none | 7 (first rank is a paralog tag, not a taxon) |
| **rpoB** | RNA polymerase subunit B; protein-coding, single-copy in almost all bacteria, and assigns far more ASVs to species than a 16S V3–V4 amplicon does. Not a strain marker — the reference stops at species, so the output does too | `pcr` | The FROGS RefSeq build ships full-length rpoB genes (~4 kb), so in-silico PCR recovers the amplicon window, ~387 bp once the two primers are trimmed off (~435 bp including them) | FROGS RefSeq rpoB build 20240707, converted to a DADA2 trainset (`rpob_frogs`) — the same file is the probe substrate | **No.** `rdp` only | none | 7 |

Names in backticks are the pack's internal reference symbols. They belong to the
`workflow/markers/*.yaml` packs, which use them to name the probe substrate
(`probe_ref`) and the taxonomy databases (`taxonomy_refdb`, `taxonomy_species_db`,
`taxonomy_sintax_db`); the workflow turns each symbol into a path when it is parsed.
They are not run-config keys. Relocating a database in the run config uses a separate
two-level key under `references:` — one level for the database family, one for the
file — as in the config template:

```yaml
references:
  silva:
    species: /shared/dbs/silva_v144_assignSpecies.fa.gz
```

A top-level `silva_train:` instead does nothing: the entry is ignored and the pack
default is used. The key that relocates each symbol is in the table below.

## Which files get downloaded, and where

Each pack carries the download URL and the default on-disk location of every database it
owns. Paths below are relative to `references.refdb_root` (default `refdb`). Files marked
*converted* or *extracted* are not stored in the form they were downloaded in.

| Marker | Symbol | File under `refdb/` | Role | Run-config key that relocates it |
|---|---|---|---|---|
| 16S | `silva_train_raw` | `silva/silva_nr99_v144_toGenus_trainset.fa.gz` | probe substrate; SILVA 144 training set as published (seven ranks) | `silva.train_raw` |
| 16S | `silva_train` | `silva/silva_nr99_v144_toGenus_7rank_other-verbatim_trainset.fa.gz` | `rdp` trainset — *converted* from `silva_train_raw` (the same seven ranks, kingdom names without `--other`) | `silva.train` |
| 16S | `silva_species` | `silva/silva_v144_assignSpecies.fa.gz` | DADA2 `addSpecies` (16S only) | `silva.species` |
| 16S | `silva_sintax` | `silva/silva_nr99_v144_7rank_other-verbatim_sintax.fa.gz` | `sintax` database — *converted* from `silva_train` | `silva.sintax` |
| ITS | `unite_fasta` | `unite/sh_general_release_dynamic_s_all_19.02.2025.fasta` | `rdp` trainset — *extracted* from a `.tgz` | `unite.fasta` |
| ITS | `unite_sintax` | `unite/utax_reference_dataset_all_19.02.2025.fasta.gz` | `sintax` database | `unite.sintax` |
| ITS | `unite_uchime_ITS1` | `unite/uchime_its1.fasta` | probe substrate when `its_region: ITS1` — *extracted* from a `.zip` | `unite.uchime_its1` |
| ITS | `unite_uchime_ITS2` | `unite/uchime_its2.fasta` | probe substrate when `its_region: ITS2` — *extracted* from a `.zip` | `unite.uchime_its2` |
| 18S | `silva_euk` | `silva_euk/silva_132.18s.dada2.fa.gz` | probe substrate only | `silva_euk.fasta` |
| 18S | `pr2_dada2` | `pr2/pr2_SSU_dada2.fasta.gz` | `rdp` trainset | `pr2.dada2` |
| 18S | `pr2_utax` | `pr2/pr2_SSU_UTAX.fasta.gz` | PR2's UTAX file as published; the input of `pr2_sintax` | `pr2.utax` |
| 18S | `pr2_sintax` | `pr2/pr2_SSU_UTAX_domain-first.fasta.gz` | `sintax` database — *converted* from `pr2_utax` (the domain moved to the letter `d:`) | `pr2.sintax` |
| gyrB | `gyrb_dada2` | `gyrb/train_set_gyrB_v6.fa.gz` | probe substrate + `rdp` trainset | `gyrb.dada2` |
| rpoB | `rpob_frogs` | `rpob/rpob_refseq_cc_20240707_domain-d_dada2.fa.gz` | probe substrate + `rdp` trainset — *converted* from the FROGS archive (the domain written `d__`) | `rpob.frogs` |

The 16S rows are those of the default release, SILVA 144; `other-verbatim` becomes
`other-stripped` with `amplicon.taxonomy.silva_other: strip`. With
`references.silva.release: "138.2"` the 16S files are
`silva/silva_nr99_v138.2_toGenus_trainset.fa.gz`, `silva/silva_v138.2_assignSpecies.fa.gz`
and `silva/silva_nr99_v138.2_sintax.fa.gz`, and there is no `silva_train_raw`. See
[16S rRNA](16S.md#reference-databases). A converted file is written where its
run-config key points, and rebuilt there when its input or its conversion script is
newer (see
[reference files MetaFlux builds itself](../../reference/configuration.md#reference-files-metaflux-builds-itself)).

Everything here is fetched on first use and cached, so a second run reuses it. Apart from
`references.silva.release`, which picks the SILVA release for 16S, the only reason to put
a `references:` entry in the run config is to point MetaFlux at a copy that already
exists somewhere else, for example a shared directory on a cluster.

## What `probe_mode` means

Before DADA2 runs, MetaFlux measures how long the amplicon is *expected* to be, and uses
that distribution to size the truncation lengths and the ASV length filter. That
measurement is the amplicon-length probe (rule `amplicon_probe`), and it only runs when
`amplicon.expected_length` is set to `auto`. `probe_mode` decides how the measurement is
taken.

**`pcr`** — in-silico PCR. The probe takes the supplied primer FASTAs and runs two
cutadapt passes against a full-length reference: first keep the reference sequences
in which the forward primer matches, trimming it and everything upstream, then, of
those, keep the ones in which the reverse-complemented reverse primer matches, trimming
it and everything downstream. Neither match is anchored — the references are whole
genes, so the primers bind inside them.
What survives is the set of amplicons those primers would produce from that reference, and
their lengths are the distribution. This is the mode for 16S, 18S and rpoB, because SILVA,
SILVA-Euk and the FROGS rpoB build all contain whole genes. For rpoB a third pass leaves
out of the statistics the products that still hold a second forward-primer site (see
[A primer that matches twice](../length-and-truncation.md#a-primer-that-matches-twice)).

**`direct`** — no PCR simulation at all. The reference is *already* amplicon-length, so
cutadapt is never called and the sequence lengths are read off the file as they are. This
is the mode for ITS, whose UNITE UCHIME reference is the pre-extracted ITS1 or ITS2
subregion, and for gyrB, whose DD7RZ8 trainset was itself built by cutting in-silico
amplicons with the Barret 2015 primers.

Running the wrong mode fails in a specific way worth knowing: in-silico PCR against a
pre-trimmed reference finds no primer binding sites left to cut against, so it recovers
nothing. For ITS the problem is biological rather than technical — ITS1F and ITS4 bind in
the flanking rRNA genes, and those binding sites are not part of a UNITE sequence (the
UCHIME references keep only a 15 bp anchor of flanking rRNA at each end).

Both modes end at the same place: one JSON of length statistics
(`n_amplicons`, `min`, `q1`, `median`, `q3`, `p95`, `p99`, `max`, `mean`, `stdev`, plus
`excluded_taxa` and `n_excluded`, and for rpoB `n_excluded_second_site`) cached
in the pipeline's own `refdb/cache/`, next to the workflow. That directory is fixed at
the repository root and does not follow `references.refdb_root`, so the cache stays put
across runs and output directories even when the databases themselves live on a shared
path. The file is keyed by marker, reference tag (for 16S including the lineages left
out of the statistics, for rpoB the third pass) and a hash of the two primer sequences. Change a primer and the probe re-runs; leave it alone and it is reused. Which
single statistic becomes `expected_length` comes from `amplicon.probe_length_stat` — see
[amplicon length and truncation](../length-and-truncation.md).

!!! note

    `probe_mode` is a pack field, not a config knob. It is validated when the workflow is
    parsed: anything other than `pcr` or `direct` stops the run immediately with
    `has an invalid probe_mode … must be 'pcr' or 'direct'`.

## What "no extractor" means

For 16S and ITS, MetaFlux can isolate the target region within each ASV before taxonomy is
assigned — Metaxa2 for 16S, ITSx for ITS. This is the `target_extract` rule, switched on
with `amplicon.extraction.enabled`.

18S, gyrB and rpoB have `extractor: none`. There is no equivalent tool for gyrB or rpoB,
and for 18S the pack deliberately does not use Metaxa2 even though it could handle SSU.
Those markers define no `target_extract` rule at all, so if the config still asks for
extraction the workflow would otherwise go looking for a file nothing produces. Instead
the setting is overridden at parse time and the reason is printed:

```text
[MetaFlux] warning: marker 18S has no target-region extractor; ignoring amplicon.extraction.enabled: true (no extraction step will run)
```

Nothing else changes. The ASV length filter still runs, and it is what keeps ASVs inside
the window derived from the probe distribution. Setting `extraction.enabled: false`
explicitly for these three markers has exactly the same effect and avoids the warning.

## Taxonomy methods per marker

`amplicon.taxonomy.method` chooses between `rdp` (DADA2's naive Bayesian classifier, with
`min_boot` as the bootstrap cutoff) and `sintax` (VSEARCH, with `sintax_cutoff`). A marker
can only offer `sintax` if a SINTAX-formatted build of its reference exists. For gyrB and
rpoB none does — DD7RZ8 and the FROGS release both ship a DADA2 trainset only — so the
pack leaves `taxonomy_sintax_db` empty and the run stops at parse time rather than failing
somewhere deep in the taxonomy rule:

```text
[MetaFlux] marker gyrB has no SINTAX reference database, so amplicon.taxonomy.method: sintax is not available for it. Use method: rdp.
```

DADA2's `addSpecies` step, which adds a species call on top of the genus-level
assignment, applies to 16S only — it is the one marker whose pack names a
`taxonomy_species_db`. The other four have no second species step because their single
reference already carries a species rank: PR2 for 18S, DD7RZ8 for gyrB, the converted
FROGS build for rpoB, and UNITE's native `s__Genus_species` for ITS.

## Rank models

Each pack declares the rank columns its classifier emits and the prefixes used to render
them (`d__`, `k__`, `p__`, …). The ranks are those of the marker's reference, nothing
added and nothing removed: `d__` is always the domain slot and `k__` always the kingdom
slot, and a table has a kingdom only where its reference defines one. So the markers
differ:

| Marker | Ranks |
|---|---|
| 16S, SILVA 144 | Domain, Kingdom, Phylum, Class, Order, Family, Genus, Species |
| 16S, SILVA 138.2 | Domain, Phylum, Class, Order, Family, Genus, Species |
| ITS | Kingdom, Phylum, Class, Order, Family, Genus, Species (UNITE has no domain) |
| 18S, `rdp` | Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species |
| 18S, `sintax` | Domain, Supergroup, Division_Subdivision, Class, Order, Family, Genus, Species |
| gyrB | Gene (the paralog tag), Phylum, Class, Order, Family, Genus, Species |
| rpoB | Domain, Phylum, Class, Order, Family, Genus, Species |

The prefix of each rank, and the shotgun tables next to these, are in
[taxonomy strings in the final tables](../../reference/output.md#taxonomy-strings-in-the-final-tables).

16S declares its ranks per SILVA release, because SILVA 144 added the kingdom that 138.2
does not have. 18S declares two models, one per taxonomy method: PR2's DADA2 file has
Division and Subdivision, its UTAX file one field that joins the two. Prefixes are keyed
to what a rank *means* rather than to the letter code in the reference, so a `d__`,
`sg__` or `c__` … `s__` token works under either method; a Division token does not (see
[18S rRNA](18S.md#pr2-ships-two-rank-depths)).

Two prefix styles exist. Under `bare` (16S, 18S) the reference values carry no prefix and
the builder adds them. Under `embedded` (ITS, gyrB, rpoB) the values already arrive with
their prefixes attached (`k__Fungi` in UNITE, `d__Bacteria` in the converted FROGS file)
and are emitted verbatim. Either way the taxonomy strings in the final table look the
same, which is what the keep/discard filter matches against — see
[keeping and discarding taxa](../taxon-filter.md).

gyrB is worth one extra sentence. Its first rank is not a taxon: DD7RZ8 tags every
reference as `gyrB` (the real gene) or `other` (a co-amplified paralog such as *parE*).
Keeping that tag as rank 1 means paralogs can be dropped with the ordinary taxon filter
(`discard: [other]`) instead of a mechanism of their own.

## The reference check (`reference_values`)

Before any ASV is classified, rule `check_taxonomy_reference` reads the exact file the
classifier is about to use, a path pinned in the run config included, and checks every
header against the pack:

- the rank shape: in a DADA2 training set (`rdp`) no lineage may be deeper than the
  pack's ranks (for 16S without the species, which `addSpecies` fills from a file of
  its own), and the deepest lineages must reach them; a shallower lineage must end with
  `;`, as SILVA's do, so that `assignTaxonomy` pads it with empty ranks; the names carry
  the pack's prefixes where the reference embeds them; a SINTAX file must use the pack's
  letters, in the pack's top-down order;
- the names in the ranks the pack has rules for, its `reference_values`.

A file that fails stops the run before any classification. The message names the file,
the rank, the wrong names with their counts, and example headers. The usual causes are a
run config that pins a wrong file (another SILVA release's, say) or a reference that
changed upstream; the fix is the right file, or a new release block in the pack. When the
file passes, the check writes `stats/taxonomy_reference_check.json`, a record of what it
saw, and the taxonomy step can start. The check takes seconds, and it starts as soon as
its reference is there, long before the denoising ends.

A rule names a rank of the pack and lists, per kind, the names allowed there:

- `one_of`: the names allowed in this rank;
- `ends_with`: allowed endings (a name passes when it is in `one_of` or ends with one of
  these);
- `not_one_of`: names that must never appear in this rank.

The shipped packs declare:

| Marker | `reference_values` |
|---|---|
| 16S, both releases | Domain: `one_of [Bacteria, Archaea, Eukaryota]` |
| 16S, SILVA 144 | Kingdom: `ends_with [ati]`, `one_of [Incertae Sedis]` |
| ITS | Kingdom: `not_one_of [Bacteria, Archaea, Eukaryota]` |
| 18S | Domain: `one_of` PR2's eight domain values (Eukaryota, Bacteria, Archaea and the organelle tags `Eukaryota:plas`, `:mito`, `:nucl`, `:apic`, `:chro`); Supergroup: `not_one_of` the same eight |
| gyrB | Gene: `one_of [gyrB, other]` |
| rpoB | Domain: `one_of [Bacteria]` |

With the rank shape, these rules catch a reference whose ranks have shifted: SILVA's own
144 download, whose kingdoms carry `--other`; a SILVA 138.2 file used as 144, with its
phyla where the kingdoms should be; PR2's own UTAX file, which holds the domain under the
letter `k:` and the supergroup under `d:`; the rpoB file of an earlier MetaFlux, which
wrote the domain as `k__`. The same rules check the keep/discard tokens when the run
starts (see
[when MetaFlux refuses a token at start-up](../taxon-filter.md#when-metaflux-refuses-a-token-at-start-up)).

When you write a pack of your own:

- Every pack declares `reference_values`, at its top level or in every release block. A
  pack without it stops the run at start-up, whichever marker the run uses. A pack for a
  marker of your own may write `reference_values: {}` to have the rank shape checked
  and no names.
- A rule applies only where its rank exists. For 18S, a rule for Division (or
  Subdivision) is tested under `method: rdp`, which reads the DADA2 file, and left out
  under `sintax`, with one line at start-up saying so; a rule for Division_Subdivision
  works the other way round. A rank that is in neither rank model stops the run.
- A 16S rule for Species stops the run at start-up: `addSpecies` fills that rank from a
  file of its own, after the classification, so neither file the check reads has a
  species to test.
- A mistake in a rule (a misspelt rank or kind, an empty list, a name YAML does not read
  as text) stops the run at start-up too, under either method.

## Per-marker pages

- [16S rRNA](16S.md) — bacteria and archaea, SILVA, Metaxa2 extraction, species assignment
- [ITS](its.md) — fungi, UNITE, ITSx extraction, ITS1 vs ITS2
- [18S rRNA](18S.md) — eukaryotes, probed on SILVA-Euk and classified on PR2, 9 ranks under `rdp` and 8 under `sintax`
- [gyrB](gyrb.md) — DD7RZ8, direct probe, paralog tag, `rdp` only
- [rpoB](rpob.md) — FROGS RefSeq, in-silico PCR probe, `rdp` only

For how the probe distribution becomes truncation lengths and an ASV length window, see
[amplicon length and truncation](../length-and-truncation.md). For the whole amplicon
path end to end, see the [amplicon overview](../overview.md).
