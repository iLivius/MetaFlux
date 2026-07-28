# ITS (fungi)

The fungal internal transcribed spacer sits between the rRNA genes, and MetaFlux
targets one of its two halves at a time — ITS1 or ITS2 — classified against
[UNITE](https://unite.ut.ee/). The marker facts live in the pack
`workflow/markers/ITS.yaml`; the run config chooses which subregion was amplified and
supplies the primers.

ITS behaves differently from every other marker in MetaFlux, and the differences are
biological rather than cosmetic. ITS length varies a lot between fungal taxa, so the
fixed-length read truncation that suits 16S would silently delete real, short
amplicons. The three sections that follow — the direct-mode probe, the `truncLen`
override, and ITSx extraction — are where that plays out.

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
    checks this combination at startup and exits with a message; either keep
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
- **`dada_length_filter`** builds the ASV window as
  `[q1 − window_margin, p95 + window_margin]`, `window_margin` defaulting to 50 bp.

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

### `min_len` and the short tail

!!! warning "A probe-derived `min_len` reintroduces the floor that `truncLen: c(0,0)` removed"

    `truncLen` is switched off for ITS precisely so that short fungal amplicons are not
    cut away. `min_len` is a per-read floor that can undo that, because a read shorter
    than the floor is discarded whatever `truncLen` is doing. MetaFlux therefore leaves
    `min_len_stat: null`, and every marker uses `filter.min_len` — 20 bp, DADA2's own
    default.

    Setting `min_len_stat` turns on a probe-derived floor for ITS, and the numbers are
    worth seeing before doing so. Measured directly from the UNITE UCHIME reference
    dataset (release 16.10.2022, the file this workflow downloads):

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

## Long amplicons are lost, and that is accepted

ITS length varies biologically, and on a fixed read length the long tail of that
distribution cannot be spanned. Where R1 and R2 no longer reach each other, `mergePairs`
discards the pair. Nothing in the workflow rescues those reads, so an ITS run is
systematically biased against its longest amplicons.

This is a real limitation of short-read ITS metabarcoding rather than a defect in the
workflow, and MetaFlux does not try to hide it. Two things make it tolerable:

- The bias is **consistent across samples in a study**. Every sample loses the same
  length classes for the same reason, so comparisons of composition between samples —
  which is what these designs are built to support — remain valid.
- Exhaustive recovery of everything the marker could theoretically retrieve is not
  available in any case: reference databases are themselves incomplete, so an ASV
  recovered but unassignable adds little.

!!! warning "`just_concatenate` is all-or-nothing, and is not a rescue for the long tail"

    It is tempting to read `dada2.merge.just_concatenate: true` as "merge what overlaps,
    concatenate the rest". It does not work that way. With this set, DADA2 joins **every**
    pair with a spacer of Ns and never attempts an overlap merge, even where the pair
    would have merged cleanly.

    That has consequences beyond the ITS page. The overlap-based error correction that
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

Only roughly, though, and it is worth knowing why. Those UCHIME files were themselves
produced with `ITSx --anchor 15`, so every reference sequence keeps about 15 bp of
flanking rRNA at each end — SSU and 5.8S for ITS1, 5.8S and LSU for ITS2. MetaFlux runs
ITSx over the ASVs with the default anchor of 0, so the reference distribution sits
roughly 30 bp above the ASVs it is being compared with. The 50 bp `window_margin`
absorbs most of that, and what is left shows up only at the bottom of the window. On
the ITS2 file the workflow downloads, q1 is 175 bp, so the floor lands at
`175 − 50 = 125` bp; that excludes about 5% of the reference sequences once their
anchors are subtracted, where a matched comparison would exclude about 2%. The
sequences in that band are the shortest fungal ITS2, so the low end of the window is
the part worth checking on a dataset expected to contain them.

| | |
|---|---|
| Input | `5.dada2/seqs.fasta` and `5.dada2/seqtab_head_names.txt`, from `dada_seqtab` |
| Command | `ITSx --only_full F -t all` |
| Region kept | the file matching `amplicon.its_region` |
| Output | `5.dada2/seqs_extracted.fasta`, `5.dada2/seqtab_extracted_head_names.txt`, `5.dada2/itsx_extraction.summary.txt`, `5.dada2/itsx_collapse_map.tsv` |
| Raw tool output | kept under `5.dada2/_itsx_tmp/` for inspection |
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
singletons included:

```text
representative    merged_asvs         n_merged   extracted_length
ASV_7             ASV_7,ASV_142       2          231
```

ASVs from which ITSx cannot extract the requested region are dropped from both the
FASTA and the count table. If the region file comes back empty, the rule stops with an
error: that usually means the amplicon does not cover the subregion named in
`its_region`, or the primers do not match the assumption.

## Taxonomy

`amplicon.taxonomy.method` selects the classifier. Both paths read the
length-filtered ASVs (`5.dada2/seqs_lenfilt.fasta` and the matching
`seqtab_lenfilt_head_names.txt`), apply the same taxon filter, and write the same
three tables into `6.taxonomy/`.

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

The visible consequence is in the species slot. UNITE writes its binomial with an
underscore and MetaFlux keeps it as it is, where the SILVA/16S path builds
`s__Genus species` with a space:

```text
ITS   k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Agaricaceae;g__Agaricus;s__Agaricus_bisporus
16S   k__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis
```

Filter tokens must be written the way the lineage writes them, underscore included.

!!! warning "The species slot is the one place the two methods disagree"

    The lineage above is the `sintax` form: the UTAX file's `s:Agaricus_bisporus` is
    carried through whole, so the string ends `;s__Agaricus_bisporus`. Under `rdp` it
    does not. DADA2's `assignTaxonomy` recognises the UNITE general release from its
    `FU|refs` / `FU|reps` headers and strips the genus out of the species field itself,
    so the same ASV comes back as `;s__bisporus`. Species-level `keep`/`discard` tokens
    therefore have to match the method actually being run: an `s__Agaricus_bisporus`
    token matches nothing on an `rdp` run — silently harmless in `discard`, and harmless
    in `keep` too unless it is the only keep token, in which case the filter empties the
    table and the run stops. Ranks above species are unaffected.

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

Each token is matched against a whole rank segment of the taxonomy string, so
`g__Fusarium` matches a genus named exactly that and never a substring elsewhere in
the lineage. If the filter removes every ASV, MetaFlux stops with an error rather than
writing an empty table — with ITS that almost always means a 16S keep list
(`[k__Bacteria, k__Archaea]`) was left in the config. Full behaviour on
[Keeping and discarding taxa](../taxon-filter.md).

## What an ITS run produces

The final tables are in `6.taxonomy/`:

| File | Rows | Contents |
|---|---|---|
| `asv_table.txt` | `ASV_#` IDs | per-sample counts plus a `taxonomy` column |
| `asv_table_seqs.txt` | ASV sequences | the same table, keyed by sequence |
| `taxon_seq_table.txt` | `ASV_#` IDs | one column per rank (Kingdom…Species) plus the sequence |

Two diagnostics are worth opening on every ITS run.
`stats/dada2/asv_length_hist.png` plots ASV lengths before extraction, after
extraction, and after the length filter — for ITS the first two differ substantially,
since ITSx removes the conserved flanks, and the gap between them is a direct measure
of how much flanking sequence the amplicon carried. `5.dada2/itsx_collapse_map.tsv`
shows which ASVs were merged after extraction. The full output layout is in
[Output files](../../reference/output.md).
