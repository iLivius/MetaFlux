# ITS (fungi)

The fungal internal transcribed spacer sits between the rRNA genes, and MetaFlux
targets one of its two halves at a time — ITS1 or ITS2 — classified against
[UNITE](https://unite.ut.ee/). The marker facts live in the pack
`workflow/markers/ITS.yaml`; the run config chooses which subregion was amplified and
supplies the primers.

ITS is handled differently from every other marker in MetaFlux. Its length varies a lot
between fungal taxa, so the fixed-length read truncation used for 16S would remove real,
short amplicons. The differences are in the direct-mode probe, the `truncLen` override
and ITSx extraction, described in the sections below.

## Setting up a run

Set `type: ITS`, choose the subregion with `its_region`, and supply the primer FASTAs
for the pair actually used. Everything else below is the shipped default, apart
from `keep` and `discard`: those ship commented out, and `[k__Fungi]` is the
recommended ITS starting point, not something already in effect.

```yaml
amplicon:
  type: ITS
  its_region: ITS2                # ITS1 | ITS2

  primers:
    fwd: /path/to/fwd_primer.fasta
    rev: /path/to/rev_primer.fasta
    orientation: fixed            # fixed | mixed

  expected_length: auto           # lengths read from the UNITE UCHIME subregion

  length_filter:
    mode: auto
    window_margin: 50

  extraction:
    enabled: true                 # ITSx

  taxonomy:
    method: sintax                # rdp | sintax
    min_boot: 80                  # rdp only
    sintax_cutoff: 0.8            # sintax only
    filter:
      enabled: true
      keep:    [k__Fungi]         # recommended, ships commented out
      discard: []                 # recommended, ships commented out

  dada2:
    trunc_q: 2                    # adaptive 3' trim — this is what controls ITS tails
    max_ee: [2, 5]
    filter:
      min_len: 20                 # DADA2's own default, used for every marker
      min_len_stat: null          # ITS-only opt-in; null = use min_len above
```

!!! warning

    `length_filter.mode: auto` and a manual `expected_length` cannot be combined for
    ITS. The ITS branch never resolves an expected length (see below), so the automatic
    length window has nothing to size itself from unless the probe has run. MetaFlux
    checks this combination at start-up and exits with a message; either keep
    `expected_length: auto`, or set `length_filter.mode: manual` with an explicit
    `length_filter.range`.

## What differs from the other markers

| Step | ITS behaviour | Why |
|---|---|---|
| Length probe | `probe_mode: direct` | the UNITE UCHIME release is already trimmed back to (near) the subregion, so lengths are measured, not PCR'd |
| `truncLen` | forced to `c(0, 0)` | a fixed cut would discard genuinely short fungal amplicons |
| 3′ quality | handled by `truncQ` + `maxEE` | adaptive per read instead of one length for all reads |
| `minLen` | `filter.min_len` (20 bp), as for every marker | a probe-derived floor is available but off by default — see below |
| Extraction | ITSx | isolates ITS1 or ITS2 from the flanking rRNA |
| Rank prefixes | `prefix_style: embedded` | UNITE values already carry `k__`, `p__`, … |
| Species step | none | UNITE lineages already end in a species name |

## Reference databases

Four files are declared by the pack, of which a given run uses three. As with every
marker, they are fetched on first use into `references.refdb_root` (default `refdb`)
and cached there.

| Pack symbol | Path under `refdb_root` | Used for | How it arrives |
|---|---|---|---|
| `unite_fasta` | `unite/sh_general_release_dynamic_s_all_19.02.2025.fasta` | RDP training set (`method: rdp`) | `.tgz` unpacked by rule `fetch_unite` |
| `unite_sintax` | `unite/utax_reference_dataset_all_19.02.2025.fasta.gz` | SINTAX database (`method: sintax`) | plain download (rule `fetch_unite_sintax`) |
| `unite_uchime_ITS1` | `unite/uchime_its1.fasta` | length probe when `its_region: ITS1` | `.zip` unpacked by rule `fetch_uchime` |
| `unite_uchime_ITS2` | `unite/uchime_its2.fasta` | length probe when `its_region: ITS2` | same `.zip`, same rule |

The two UCHIME files come out of a single archive, so both are written even though a
run probes only one. They are used **only** to measure amplicon lengths — the
taxonomy always comes from the general release or the UTAX file, never from these.

An entry under `references:` in the run config is only needed to relocate a file, for
example to a shared filesystem:

```yaml
references:
  unite:
    fasta:       /shared/dbs/sh_general_release_dynamic_s_all_19.02.2025.fasta
    sintax:      /shared/dbs/utax_reference_dataset_all_19.02.2025.fasta.gz
    uchime_its1: /shared/dbs/uchime_its1.fasta
    uchime_its2: /shared/dbs/uchime_its2.fasta
```

## Amplicon length probe (`probe_mode: direct`)

For 16S the amplicon window is found by an in-silico PCR against full-length genes.
That approach does not work for ITS: the common ITS primers (ITS1F, ITS4 and
relatives) bind in the flanking rRNA genes, and those binding sites are not part of a
UNITE sequence — there would be nothing for the primers to match.

The UNITE UCHIME release solves it from the other side. Its sequences are already cut
back to the ITS1 or ITS2 subregion plus a 15 bp anchor of flanking rRNA at each end, so
the `amplicon_probe` rule simply reads their lengths straight off the FASTA, without
running Cutadapt at all, and gzip-copies the file into the cache so the rule's declared
outputs exist either way:

```text
refdb/cache/probe_ITS_unite_uchime_ITS2_<primer_hash>.json
refdb/cache/probe_ITS_unite_uchime_ITS2_<primer_hash>.amplicons.fa.gz
```

`<primer_hash>` is the first 12 hex characters of a SHA-256 over the two primer files,
and the reference tag follows `its_region`, so switching ITS2 → ITS1 produces a
separate cache entry rather than reusing the wrong distribution.

The JSON holds the same statistics as any other probe — `n_amplicons`, `min`, `q1`,
`median`, `q3`, `p95`, `p99`, `max`, `mean`, `stdev` — but for ITS they feed two
consumers only:

- **`dada_filter`** can take the per-read length floor `minLen` from
  `probe[dada2.filter.min_len_stat]`, but this is **off by default**
  (`min_len_stat: null`) and the config floor `dada2.filter.min_len` — 20 bp, DADA2's
  own default — is used instead. Only ITS can use a probe-derived floor at all: the
  extracted subregion sits below read length, whereas for 16S the probe measures the
  full amplicon and would exceed the read length and drop everything. Why it is not
  the default is set out under [the length floor](#min_len-and-the-short-tail).
- **`dada_length_filter`** builds the ASV window as `[min, p95 + window_margin]`,
  `window_margin` defaulting to 50 bp: the floor is the shortest region in the
  reference, without margin (see [the length floor of the window](#region-extraction-with-itsx)).
  With `extraction.enabled: false` it is `[q1 − window_margin, p95 + window_margin]`, as
  for the other markers.

Nothing about `truncLen` reads the probe for ITS, and `amplicon.min_overlap` and
`probe_length_stat` are not applied to an ITS run at all.

## Why `truncLen` is switched off

In DADA2, `truncLen` does two things at once: it cuts every read to a fixed length,
**and** it discards any read shorter than that length. For a marker of near-constant
size that is exactly what is wanted. ITS is not such a marker — its length varies
between fungal groups, and a genuinely short ITS amplicon is a real biological
observation, not a truncated read. A fixed cut would throw those taxa away and bias
the community toward whatever happens to be long.

So for ITS the `dada_filter` rule sets `truncLen = c(0, 0)`, meaning "no fixed cut",
before calling `filterAndTrim`. The override is unconditional: it also wins over
`trunc_len.mode: manual`, so `manual_r1` and `manual_r2` have no effect on an ITS run.

The 3′ low-quality tail is still removed — just adaptively, per read:

| Setting | Default | What it does on an ITS run |
|---|---|---|
| `dada2.trunc_q` | `2` | truncates each read at its first base with quality at or below this value |
| `dada2.max_ee` | `[2, 5]` | drops a read whose expected errors exceed the limit, R1 and R2 separately |
| `dada2.filter.min_len` / `min_len_stat` | `20` / `null` | per-read length floor; the probe-derived option is off by default |

If ITS R2 quality is poor, the lever is therefore a stricter `trunc_q` or a tighter
`max_ee` — not a fixed truncation length.

The one fixed cut ITS reads get is the extra last sequencing cycle of a 2 × 151,
2 × 251 or 2 × 301 run, removed at primer trimming (`cutadapt.trim_final_cycle: auto`).
Its calls are often wrong the same way in most reads of a sequence, and without
truncation the wrong base reaches the overlap of the two mates, where `mergePairs`
accepts no mismatch.
See [the extra last cycle](../overview.md#the-extra-last-cycle) and
[the validation](#validation).

### `min_len` and the short tail

!!! warning "A probe-derived `min_len` reintroduces the floor that `truncLen: c(0,0)` removed"

    `truncLen` is switched off for ITS precisely so that short fungal amplicons are not
    cut away. `min_len` is a per-read floor that can undo that, because a read shorter
    than the floor is discarded whatever `truncLen` is doing. MetaFlux therefore leaves
    `min_len_stat: null`, and every marker uses `filter.min_len` — 20 bp, DADA2's own
    default.

    Setting `min_len_stat` turns on a probe-derived floor for ITS. The UNITE UCHIME
    reference dataset (release 16.10.2022, the file this workflow downloads) has these
    lengths:

    | Subregion | n | min | q1 | median | q3 | p95 | max |
    |---|--:|--:|--:|--:|--:|--:|--:|
    | ITS1 | 160,204 | **17** | 174 | 200 | 230 | 293 | 2,427 |
    | ITS2 | 154,918 | **38** | 175 | 204 | 233 | 293 | 4,461 |

    `min_len_stat: q1` puts the floor near **175 bp**, above roughly **24%** of the
    reference distribution in both subregions — so reads from genuinely short ITS
    amplicons are discarded, which is exactly the bias the `truncLen` exemption exists
    to avoid. `min_len_stat: min` avoids that but is not uniform across subregions
    either: 38 bp is a reasonable floor for ITS2, while ITS1 bottoms out at 17 bp, below
    DADA2's own default and low enough to admit uninformative fragments.

    Leaving it `null` sidesteps the choice: a flat 20 bp floor removes nothing of
    biological interest from either subregion, and length filtering of the assembled
    ASVs still happens later in `dada_length_filter`, where the probe distribution is
    applied to whole amplicons rather than to individual reads.

The `pick_trunclen` rule still runs and still writes `stats/trunclen.json`, but what
lands in it depends on `trunc_len.mode`. Under `auto` the rule computes its usual
quality-based cuts, skips the merge-overlap constraint, and records
`"expected_length": null` together with
`"expected_length_source": "not_applicable:ITS_truncLen_overridden"`. Under `manual` it
returns before the quality analysis runs at all and writes the configured `r1` and `r2`
plus `"resolved_via": "manual"`. Either way the
numbers are QC information, not instructions: `dada_filter` overrides them to `c(0, 0)`.

!!! tip

    Pairs that do not overlap after adaptive trimming are lost at `mergePairs`. Two
    config keys exist for short or awkward ITS amplicons: `dada2.merge.trim_overhang`
    (default `true`) removes bases that run past the far end of the amplicon, and
    `dada2.merge.just_concatenate` (default `false`) joins R1 and R2 with a run of Ns
    instead of merging them on overlap.

## Long amplicons that do not merge are lost

ITS length varies biologically, and on a fixed read length the long tail of that
distribution cannot be spanned. Where R1 and R2 no longer reach each other, `mergePairs`
discards the pair. Nothing in the workflow rescues those reads, so an ITS run is
systematically biased against its longest amplicons.

This is a limitation of short-read ITS metabarcoding, not of the workflow. Two points
limit its effect:

- The bias is **consistent across samples in a study**. Every sample loses the same
  length classes for the same reason, so comparisons of composition between samples —
  which is what these designs are built to support — remain valid.
- Exhaustive recovery of everything the marker could theoretically retrieve is not
  available in any case: reference databases are themselves incomplete, so an ASV
  recovered but unassignable adds little.

!!! warning "`just_concatenate` is all-or-nothing, and is not a rescue for the long tail"

    `dada2.merge.just_concatenate: true` does not merge the pairs that overlap and
    concatenate the rest: DADA2 joins **every** pair with a spacer of Ns and never
    attempts an overlap merge, even where the pair would have merged cleanly.

    This has several consequences. The overlap-based error correction that
    merging provides is gone for every ASV; the joined length is R1 + R2 + spacer rather
    than the amplicon length, so the [ASV length filter](../length-and-truncation.md)
    is no longer measuring what it thinks it is; the N run and the unsequenced middle
    degrade taxonomic assignment; and ASVs are no longer comparable between runs
    sequenced at different read lengths.

    Mixing the two treatments within one dataset would also be its own bias — some ASVs
    error-corrected across an overlap, others not — so the parameter deliberately does
    not offer it. Leave it `false` unless there is a specific reason, and treat the lost
    long amplicons as a known and reportable property of the assay.

## Region extraction with ITSx

With `extraction.enabled: true` (the default), the `target_extract` rule runs
[ITSx](https://microbiology.se/software/itsx/) over the ASVs to cut away the flanking
rRNA (SSU, 5.8S, LSU) and keep the spacer itself. This is what brings ASV lengths onto
roughly the same scale as the probe distribution, which was measured on the UNITE
UCHIME ITS1/ITS2 files.

The match is approximate: those UCHIME files were themselves produced with
`ITSx --anchor 15`, so every reference sequence keeps about 15 bp of
flanking rRNA at each end — SSU and 5.8S for ITS1, 5.8S and LSU for ITS2. MetaFlux runs
ITSx over the ASVs with the default anchor of 0, so the reference distribution sits
roughly 30 bp above the ASVs it is being compared with (*Candida haemulonii* in the test
set: a 93 bp reference record, a 63 bp region in the ASV). At the top of the window the
50 bp `window_margin` absorbs most of that.

The floor of the window is the shortest region in the reference, the probe's `min`,
without margin. A short region is a short ITS, not a fragment: ITSx reports a region only
where it finds the 5.8S boundary, and a region whose SSU or LSU boundary it misses keeps
some rRNA, so it comes out longer, not shorter.

| | ITS2 | ITS1 |
|---|---|---|
| Floor (probe `min`) | 38 bp | 17 bp |
| Shortest reference records | 4 records of lichen fungi, 38–41 bp | 62 records, 17–19 bp |
| Reference sequences below the floor, anchors subtracted | 0.15 % | 0.21 % |

The minimum comes from a few unusual records, so the floor guards against fragments
rather than marking a biological limit. The floor of the other markers,
`q1 − window_margin`, would be 125 bp for ITS2: on the ITS test set it removed
*Candida haemulonii* (63 bp), *Metschnikowia pulcherrima* (70 bp) and
*Dipodascus australiensis* (78 bp), 64 % of the reads after chimera removal.

| | |
|---|---|
| Input | `05.dada2/seqs.fasta` and `05.dada2/seqtab_head_names.txt`, from `dada_seqtab` |
| Command | `ITSx --only_full F -t all` |
| Region kept | the file matching `amplicon.its_region` |
| Output | `05.dada2/seqs_extracted.fasta`, `05.dada2/seqtab_extracted_head_names.txt`, `05.dada2/itsx_extraction.summary.txt`, `05.dada2/itsx_collapse_map.tsv` |
| Raw tool output | kept under `05.dada2/_itsx_tmp/` for inspection |
| Consumed by | `dada_length_filter` |

`--only_full F` keeps partial detections, and it is mandatory here: a targeted
ITS1-or-ITS2 amplicon never contains the complete ITS region flanked by both SSU and
LSU, so requiring a full detection would reject everything. `-t all` scans all
eukaryote profiles.

**Collapsing duplicates.** Two distinct DADA2 ASVs can differ only in the flanking
5.8S or LSU bases — common with multi-primer or multiplexed designs. Once ITSx trims
those flanks away the two sequences are identical, and treating them as separate
markers would split one taxon's abundance in half. MetaFlux collapses such ASVs into
one representative (the lowest-numbered ASV of the group, i.e. the most abundant, since
DADA2 numbers ASVs by decreasing abundance) and sums their per-sample counts.
`itsx_collapse_map.tsv` records the outcome for every representative — one row each,
singletons included. One row from the ITS test set
([PRJNA862334](../../about/test-datasets.md)):

```text
representative    merged_asvs         n_merged   extracted_length
ASV_3             ASV_3,ASV_18        2          232
```

ASVs from which ITSx cannot extract the requested region are dropped from both the
FASTA and the count table. If the region file comes back empty, the rule stops with an
error: that usually means the amplicon does not cover the subregion named in
`its_region`, or the primers used do not actually target that subregion.

## Taxonomy

`amplicon.taxonomy.method` selects the classifier. Both paths read the
length-filtered ASVs (`05.dada2/seqs_lenfilt.fasta` and the matching
`seqtab_lenfilt_head_names.txt`), apply the same taxon filter, and write the same
four tables into `06.taxonomy/`: the three filtered tables and `asv_table_unfiltered.txt`,
whose `metaxa2_origin` column is `NA` for ITS.

| | `rdp` | `sintax` |
|---|---|---|
| Tool | DADA2 `assignTaxonomy` | `vsearch --sintax` |
| Reference | `unite_fasta` | `unite_sintax` |
| Confidence setting | `min_boot: 80` | `sintax_cutoff: 0.8` |
| Reverse complement | `try_rc: true` | `--strand both`, always |
| Species step | none — see below | none |
| Conda environment | `workflow/envs/taxonomy.yaml` | `workflow/envs/vsearch.yaml` |

There is no `addSpecies` step for ITS: the pack sets `taxonomy_species_db: null`
because UNITE lineages already run all the way to species. 16S needs the extra step
only because its SILVA training set stops at genus.

The RDP path classifies the *unique* ASV sequences and expands the result back
afterwards, so no sequence is classified twice. On an ITS run that step finds
nothing to do: the ITSx collapse above has already merged the ASVs that extraction
left identical, so `seqs_extracted.fasta` — and the `seqs_lenfilt.fasta` built from
it — carries one entry per distinct sequence. The de-duplication matters for markers
whose extractor collapses nothing, 16S with Metaxa2 being the case in point.

### Embedded rank prefixes

UNITE values arrive already carrying their rank prefix — `k__Fungi`, `p__Basidiomycota`
and so on. The ITS pack therefore sets `prefix_style: embedded`, which tells the
taxonomy-string builder to emit the values verbatim rather than adding a second set of
prefixes on top. The SINTAX path reaches the same place from the other direction: it
strips any prefix the database embedded, then re-adds the pack's own, so the rank
prefixes come out identical whichever method ran.

UNITE has no domain rank: its top rank is the kingdom, so an ITS string starts at `k__`
and has no `d__` slot. The UTAX file (`sintax`) stores that kingdom under the SINTAX
letter `d:` (`tax=d:Fungi,p:Ascomycota,…`); MetaFlux reads it as the kingdom, writes it
`k__`, and never turns it into a domain. The kingdom slot holds whatever UNITE puts at
that rank: the kingdoms (`k__Fungi`, `k__Viridiplantae`, `k__Metazoa`), protist groups
(`k__Alveolata`, `k__Stramenopila`, `k__Rhizaria` …) and placeholders (`k__Protista`,
`k__unidentified`, `k__Eukaryota_kgd_Incertae_sedis`).

The prefixes belong to the taxonomy strings. In `taxon_seq_table.txt`, where each rank
has its own column, the cells hold the bare names — `Fungi`, `Basidiomycota`, …,
`Agaricus_bisporus` — under either method.

Two differences from a 16S string are visible side by side. The ITS string starts at
the kingdom, where the 16S one starts at the domain (`d__`; on SILVA 144 the kingdom
`k__` follows it). And UNITE writes its binomial with an underscore, which MetaFlux keeps
as it is, where the SILVA/16S path builds `s__Genus species` with a space:

```text
ITS   k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Agaricaceae;g__Agaricus;s__Agaricus_bisporus
16S   d__Bacteria;k__Bacillati;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis
```

Filter tokens must be written the way the lineage writes them, underscore included.

!!! note "Species names under `rdp`"

    The UTAX file's `s:Agaricus_bisporus` is carried through whole by `sintax`. DADA2's
    `assignTaxonomy` instead recognises the UNITE general release from its `FU|refs` /
    `FU|reps` headers and cuts the genus out of every species field, so on its own it
    would return `;s__bisporus`. MetaFlux puts the classifier's genus back in front of
    the epithet before anything else reads the result, so under `rdp` the string also
    ends `;s__Agaricus_bisporus`, and a species token such as `s__Agaricus_bisporus`
    works the same with either method, in `keep`, in `discard` and in the
    [best-hit check](../taxon-filter.md#the-best-hit-check-discard_min_support) of `discard`.

    Varieties, subspecies and forms are the exception, about 1 % of UNITE's named
    species. Of those DADA2 keeps only the part from `var.`, `subsp.` or `f.` onwards,
    so the epithet in front of it is lost and the name cannot be rebuilt: under `rdp`,
    `s__Diaporthe_phaseolorum_var._meridionalis` is written `s__var._meridionalis`,
    while `sintax` writes the full name. A token for such a name matches only under
    `sintax`.

### `rdp` and `sintax` on the same ASVs

The two methods write the same form, but they are different classifiers reading two
different UNITE files, so they do not give the same answers:

- **Depth.** On the ITS test set ([PRJNA862334](../../about/test-datasets.md), a yeast
  mock, 12 ASVs, UNITE release 19.02.2025), `sintax` at `sintax_cutoff: 0.8` named the
  species of 8 ASVs and the genus of 2 more. `rdp` at `min_boot: 80` named the species
  of 1 and stopped at phylum or kingdom for 8; the other 3 got no kingdom at all, so
  `keep: [k__Fungi]` removed them (in 10 of 11 identical `rdp` runs; ASV_16, below, is
  the exception).
- **Placeholders.** Where UNITE has no name for a rank, the general release (`rdp`)
  carries an *Incertae sedis* placeholder such as `c__Ascomycota_cls_Incertae_sedis`,
  while the UTAX file (`sintax`) leaves the rank out: the `sintax` string skips it and
  its cell is `NA`. The 3,701 reference sequences of unknown kingdom read
  `k__Eukaryota_kgd_Incertae_sedis` in the general release and `k__unidentified` under
  `sintax`.

## Keeping and discarding taxa

There is no built-in default — the keep and discard lists in the shipped config are
commented out on purpose, so a run that forgets them filters nothing rather than
silently applying another marker's tokens. The recommended ITS starting point:

```yaml
taxonomy:
  filter:
    enabled: true
    keep:    [k__Fungi]
    discard: []
```

`keep: [k__Fungi]` keeps ASVs whose kingdom is Fungi and drops everything else,
including ASVs left unclassified at kingdom level. ITSx scans all eukaryote profiles,
so non-fungal eukaryotic ITS can reach this stage. Nothing is recommended
for `discard` by default; lab or kit contaminants revealed by a negative control can
be added there as genus tokens.

With `method: rdp`, DADA2 breaks ties between equally good references at random, and
this can change an ASV's whole lineage between identical runs. On the ITS test set,
ASV_16 had an orchid as best hit in 8 of 11 identical runs and *Saccharomyces* in the
other 3. Its kingdom support reached `min_boot: 80` in only one of those three, so
`keep: [k__Fungi]` kept the ASV in 1 run and removed it in the other 10. `sintax` gives
the same answer every time.

Each token is matched against a whole rank segment of the taxonomy string, so
`g__Fusarium` matches a genus named exactly that and never a substring elsewhere in
the lineage. A 16S keep list left in the config stops the run at start-up, before any
step: ITS has no domain rank, so `d__Bacteria` can never match, and `k__Bacteria` puts a
domain where ITS has its kingdom. If the filter removes every ASV after the
classification, MetaFlux stops with an error rather than writing an empty table. Full
behaviour on [Keeping and discarding taxa](../taxon-filter.md).

## What an ITS run produces

The final tables are in `06.taxonomy/`:

| File | Rows | Contents |
|---|---|---|
| `asv_table.txt` | `ASV_#` IDs | per-sample counts plus a `taxonomy` column |
| `asv_table_seqs.txt` | ASV sequences | the same table, keyed by sequence |
| `taxon_seq_table.txt` | `ASV_#` IDs | one column per rank (Kingdom…Species) plus the sequence; the cells hold the bare names (`Fungi`, `Hanseniaspora_uvarum`) |
| `asv_table_unfiltered.txt` | `ASV_#` IDs, every ASV that reached the filter | the same counts and `taxonomy`, plus `metaxa2_origin` (`NA` for ITS) and `filter_status`: `kept`, or the filter step that removed the ASV |

Two diagnostics are worth opening on every ITS run.
`stats/dada2/asv_length_hist.png` plots ASV lengths before extraction, after
extraction, and after the length filter — for ITS the first two differ substantially,
since ITSx removes the conserved flanks, and the gap between them is a direct measure
of how much flanking sequence the amplicon carried. `05.dada2/itsx_collapse_map.tsv`
shows which ASVs were merged after extraction. The full output layout is in
[Output files](../../reference/output.md).

## Validation

The six samples of the [ITS test set](../../about/test-datasets.md) (yeast mocks A, B, C
and F of PRJNA862334, MiSeq, reads of 301 bp) were run with the defaults and `keep: [k__Fungi]`,
with the two changes of v2.7.0 switched on one at a time. Each ASV was matched to the
ITS2 sequences of the species in the mocks, so detection does not depend on the
classifier.

| | v2.6.1 | last cycle removed | short-ITS floor | both (v2.7.0) |
|---|--:|--:|--:|--:|
| ASVs | 12 | 14 | 16 | 18 |
| Retention, % of raw read pairs | 16.0 | 33.7 | 44.5 | 62.2 |
| Expected ITS2 types detected, of 26 | 12 | 19 | 19 | 26 |
| Expected genera detected, of 17 | 9 | 10 | 16 | 17 |
| *S. paradoxus*, of 4 samples | 0 | 4 | 0 | 4 |
| *S. uvarum* group, of 2 samples | 0 | 2 | 0 | 2 |
| Short-ITS2 yeasts, of 7 sample–species pairs | 0 | 0 | 7 | 7 |

Counts are summed over the samples: a mock of four species sequenced in two samples
counts 8 expected types. A type is a species, or a group of species with the same ITS2
(the *S. uvarum* group). *Retention* is the share of the raw read pairs that reach the
final table. The short-ITS2 yeasts are *C. haemulonii*, *M. pulcherrima* and
*D. australiensis*.

- **The two changes recover different species.** Removing the last cycle lets the pairs
  of *S. paradoxus* and of the *S. uvarum* group merge; the floor keeps the short-ITS2
  yeasts. Together they detect every expected type in every sample.
- **`max_mismatch: 1` adds little.** With both changes it adds 2 ASVs, 2.8 points of
  retention and no expected type, so the default stays at 0.
- **Class-only ASVs.** SINTAX stops at class (Saccharomycetes) for 5 of the 18 ASVs,
  32 % of the reads: four *Saccharomyces* ASVs, whose species share a near-identical ITS2
  in this amplicon, and one variant of *M. pulcherrima*.
- **Species outside the mock lists.** *Rhodotorula diobovata* (mock A) and *Yamadazyma
  mexicana* (mock B) are found in v2.6.1 and v2.7.0 alike.

Detection is complete; the proportions are not. The next table compares, per sample, the
read share of the short-ITS2 yeasts with their share in the mock, and gives the genus
distance: the Bray–Curtis dissimilarity between the expected and the observed genus
shares, from 0 (identical) to 1 (nothing in common), with genera assigned by sequence.

| Mock (samples) | Short-ITS2 yeasts | Expected | Observed | Genus distance, v2.6.1 | v2.7.0 |
|---|---|--:|--:|--:|--:|
| A (1) | *C. haemulonii*, *M. pulcherrima* | 40 % | 72.3 % | 0.40 | 0.32 |
| B (1) | *D. australiensis* | 20 % | 44.8 % | 0.20 | 0.29 |
| C (2) | *C. haemulonii*, *M. pulcherrima* | 10 % | 47.6 %, 48.1 % | 0.13, 0.15 | 0.38, 0.38 |
| F (2) | none | 0 % | 0.6 %, 0.5 % | 0.15, 0.90 | 0.16, 0.17 |

The short-ITS2 yeasts get two to five times their share. Their reads add to those of the
other species instead of replacing them (*S. cerevisiae* in mock C: 14,638 and 12,086
reads in v2.6.1, 14,814 and 12,245 in v2.7.0), so the excess is already in the reads, as
expected from PCR, which favours short templates. The distance therefore rises in mocks B
and C although every expected type is found. It falls in mock A, where the short-ITS2
yeasts make up 40 % of the community, and in mock F, where the *Saccharomyces* pairs now
merge. Read shares of ITS amplicons this different in length are not abundances.
