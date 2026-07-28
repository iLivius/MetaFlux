# gyrB (bacterial DNA gyrase subunit B)

`gyrB` encodes the B subunit of DNA gyrase, a protein-coding single-copy gene in bacteria.
Because it is protein-coding, it accumulates substitutions faster than the 16S rRNA gene,
so it separates organisms inside a genus that a 16S V3–V4 amplicon leaves merged. The cost
is a reference database much smaller than SILVA, and a marker that is bacterial rather than
universal.

MetaFlux runs gyrB through the same amplicon path as every other marker — PhiX removal,
primer trimming, DADA2 denoising, length filtering, taxonomy. Three things set it apart
from the rRNA markers, all of them consequences of its reference database
(`workflow/markers/gyrB.yaml` declares them): the reference is already amplicon-length,
only the RDP classifier is available, and the reference carries a built-in paralog tag.
The first needs no attention at all; the other two need one config line each.

## At a glance

| | |
|---|---|
| `amplicon.type` | `gyrB` (also accepted as `gyrb` or `GYRB` — matched case-insensitively) |
| Reference database | DD7RZ8 v6, the INRAE gyrB DADA2 trainset built on GTDB r226 |
| Role of that database | both the taxonomy reference and the amplicon-length probe substrate |
| Probe mode | `direct` — lengths read straight off the reference |
| Classifier | RDP only (`amplicon.taxonomy.method: rdp`) |
| Separate species step | none — the trainset itself carries species |
| Region extractor | none (no gyrB equivalent of Metaxa2 or ITSx) |
| Rank model | `Gene;Phylum;Class;Order;Family;Genus;Species` |
| Primers | none are shipped; Barret 2015 F64 / R353 is the pair the reference was built with |

## Setting up a run

Set the marker type, point at the primer FASTAs, and choose the RDP classifier. A minimal
gyrB block looks like this:

```yaml
mode: amplicon

amplicon:
  type: gyrB
  primers:
    fwd: /path/to/gyrB_F64.fasta
    rev: /path/to/gyrB_R353.fasta
    orientation: fixed          # fixed | mixed
  expected_length: auto
  probe_length_stat:
    16S: p95                    # gyrB reads this key — see "Amplicon length" below
  extraction:
    enabled: false              # no gyrB region extractor
  taxonomy:
    method: rdp                 # gyrB has no SINTAX build
    min_boot: 80
    try_rc: true
    filter:
      enabled: true
      keep:    []
      discard: []               # [other] drops paralogs — see "The paralog tag"
```

Primers are never supplied by MetaFlux: `amplicon.primers.fwd` and `amplicon.primers.rev`
must both point at existing FASTA files, and parsing stops immediately if either is
missing. The Barret 2015 F64 / R353 pair is the natural choice, because DD7RZ8's reference
amplicons were cut with those same primers, but any gyrB pair works — the length window is
measured from the reference rather than assumed.

!!! warning "The shipped config template defaults to `method: sintax`"

    `config/config.yaml` ships with `amplicon.taxonomy.method: sintax`, which is fine for
    16S/ITS/18S but not available for gyrB. A run config copied from the template without
    changing that line stops at parse time, before any compute is spent:

    ```text
    [MetaFlux] marker gyrB has no SINTAX reference database, so
    amplicon.taxonomy.method: sintax is not available for it. Use method: rdp.
    ```

## The reference database

DD7RZ8 v6 is a gyrB trainset in DADA2 format, built at INRAE on GTDB release r226. It is
downloaded on first use and cached, so later runs reuse it.

| | |
|---|---|
| Cached at | `refdb/gyrb/train_set_gyrB_v6.fa.gz` |
| Source | Recherche Data Gouv (INRAE Dataverse) file-access endpoint, `https://entrepot.recherche.data.gouv.fr/api/access/datafile/760120` |
| Fetch rule | `fetch_gyrb_dada2`, generated from the URL in the marker pack |
| Log | `<out_dir>/logs/refdb/fetch_gyrb_dada2.log` |

The download is a single plain `.fa.gz` file, so no rule was written by hand for it:
MetaFlux generates one fetch rule per plain-download reference declared by the active
marker pack. The URL is a Dataverse API endpoint that redirects to the object store, which
`wget` follows.

The cache root is `references.refdb_root` (default `refdb`). An entry under `references`
is only needed to relocate a database that already exists elsewhere on disk — for gyrB
that entry is:

```yaml
references:
  gyrb:
    dada2: /shared/dbs/train_set_gyrB_v6.fa.gz
```

Every reference sequence is one primer-trimmed in-silico amplicon, not a full-length gene,
and the header is a GTDB lineage whose ranks already carry their prefixes, preceded by a
bare gene tag:

```text
>{gene};p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Genus_species
```

Because the prefixes are already in the reference, the marker pack sets
`prefix_style: embedded` and MetaFlux emits each rank value verbatim instead of prepending
its own `p__`/`c__`/… — the species slot stays a GTDB binomial such as
`s__Actinotalea_bogoriensis`.

## Amplicon length: measured directly, not by in-silico PCR

For 16S, 18S and rpoB the amplicon-length probe runs an in-silico PCR: two passes of
cutadapt against a full-length reference, keeping the sequences that carry the forward
primer at the 5′ end and the reverse complement of the reverse primer at the 3′ end. That
is impossible against DD7RZ8, whose sequences have already had those primer sites cut off
upstream — the PCR would find nothing to match and return almost no amplicons.

So the gyrB pack sets `probe_mode: direct`: cutadapt is not run at all, and the length
distribution is read straight from the reference sequences. Measured across the shipped
DD7RZ8 v6 file, that distribution is tight — median 247 bp, minimum 234 bp, maximum
270 bp — which is exactly what a pre-trimmed amplicon reference should look like.

The probe only runs when `amplicon.expected_length: auto`. It writes its statistics to a
JSON file whose name encodes the marker, the reference tag and a 12-character hash of both
primer FASTAs:

```text
refdb/cache/probe_gyrB_gyrb_dd7rz8_v6_<primer_hash>.json
```

Changing primers changes the hash, so the probe reruns; keeping the same primers across
runs reuses the cached result. That JSON then feeds two downstream steps: `pick_trunclen`,
which uses the resolved expected length to enforce the forward/reverse overlap constraint
on `truncLen`, and the ASV length filter, which builds its `auto` window from the probe's
q1 and p95 plus `length_filter.window_margin`. See
[Amplicon length and truncation](../length-and-truncation.md) for both.

!!! warning "gyrB reads `probe_length_stat.16S`, not a key of its own"

    Which statistic of the probe distribution becomes `expected_length` is chosen by
    `amplicon.probe_length_stat.<key>`, and the gyrB pack deliberately reuses the `16S`
    key rather than defining a separate one. If the `16S` entry is missing, parsing
    stops with:

    ```text
    [MetaFlux] amplicon.probe_length_stat is missing the '16S' key that marker gyrB
    needs. Add it, e.g.
        probe_length_stat:
          16S: p95
    ```

    The key is checked at parse time for every non-ITS marker, whether or not
    `expected_length` is `auto`. A manual `expected_length` (a single integer, or a
    `[min, max]` pair) skips the probe itself, but the `16S` entry must still be present.

## Taxonomy: RDP only

The DD7RZ8 release ships a DADA2 trainset and nothing else — there is no SINTAX-formatted
build of it — so the marker pack leaves `taxonomy_sintax_db` empty and MetaFlux refuses
`method: sintax` at parse time rather than failing later inside the classifier.

Classification is DADA2's `assignTaxonomy` (the RDP naive Bayesian classifier) against the
trainset, with `taxonomy.min_boot` as the bootstrap confidence minimum and
`taxonomy.try_rc` deciding whether reverse-complemented sequences are also tried. There is
no separate `addSpecies` step: unlike 16S, where species come from a second SILVA file,
the gyrB trainset is itself labelled to species, so the classifier reaches the species rank
directly. `amplicon.seed` fixes the bootstrap RNG, so RDP results are reproducible run to
run.

The classifier writes seven rank columns, taken from the marker pack:

| Slot | Rank column | Prefix in the taxonomy string |
|---|---|---|
| 1 | `Gene` | none — a bare `gyrB` or `other` tag |
| 2 | `Phylum` | `p__` |
| 3 | `Class` | `c__` |
| 4 | `Order` | `o__` |
| 5 | `Family` | `f__` |
| 6 | `Genus` | `g__` |
| 7 | `Species` | `s__` |

## The paralog tag

DNA gyrase subunit B has a close relative, topoisomerase IV subunit B (*parE*), similar
enough that gyrB primers co-amplify it. DD7RZ8 handles this by labelling every reference
sequence in its first rank slot: `gyrB` for the true gene, `other` for a co-amplified
paralog such as *parE*.

MetaFlux keeps that tag as an ordinary taxonomy rank, which means paralogs can be removed
with the same keep/discard filter used to drop chloroplast and mitochondrial ASVs from a
16S run — no extra machinery, no separate step:

```yaml
amplicon:
  taxonomy:
    filter:
      enabled: true
      discard: [other]     # not the same as keep: [gyrB] — see below
```

Filter tokens are matched against whole segments of the taxonomy string, so `other` matches
the gene tag and only the gene tag — never a substring somewhere inside a species name.
Because the tag carries no rank prefix, the token is written bare, unlike the prefixed
tokens (`k__Bacteria`, `o__Chloroplast`) used for the rRNA markers. See
[Keeping and discarding taxa](../taxon-filter.md).

`discard: [other]` drops the paralogs, and `keep: [gyrB]` is not an equivalent way of
writing it. An ASV the classifier could not tag at all has an empty taxonomy string, so it
carries no segment for either list to match: `discard` leaves it in, `keep` throws it out.
In the shipped gyrB test run that is 5 of 256 ASVs — 2% lost silently by the `keep` form.

!!! tip "Leave the filter empty on the first run"

    How much `other` to expect depends on the primers and the samples, so there is no
    universal right answer and MetaFlux applies no default. Running once with empty
    keep/discard lists makes the paralog fraction visible in the output — count the ASVs
    whose taxonomy string starts with `other`, or read the `Gene` column of
    `6.taxonomy/taxon_seq_table.txt` — and then decide whether to discard them.

    The fraction is not random noise: co-amplification of *parE* is taxon-dependent, so
    it concentrates in particular clades rather than spreading evenly across the table.
    A run dominated by the affected groups will show a much larger `other` fraction than
    one that is not, which is worth checking against the composition before assuming a
    primer problem.

If a keep/discard combination removes every ASV, `assign_taxonomy` stops with an error
instead of writing a header-only table, since that is nearly always a marker/config
mismatch.

!!! warning "GTDB splits genera, and the filter matches whole segments exactly"

    DD7RZ8 v6 carries GTDB r226 taxonomy, and GTDB routinely divides a genus that is
    polyphyletic under NCBI into several suffixed genera — `g__Pseudomonas_A`,
    `g__Pseudomonas_B`, and so on. In the shipped reference, **734 of 5,207 genera**
    carry such a suffix, and *Pseudomonas* alone appears as seven separate names:

    ```text
    g__Pseudomonas   g__Pseudomonas_B   g__Pseudomonas_E   g__Pseudomonas_M
    g__Pseudomonas_O   g__Pseudomonas_P   g__Pseudomonas_T
    ```

    The keep/discard filter matches whole `;`-delimited segments exactly — deliberately,
    so that `g__Bacillus` can never match a substring elsewhere in a lineage. The
    consequence here is that `discard: [g__Pseudomonas]` removes **only** the unsuffixed
    genus and leaves the other six untouched, and `keep: [g__Pseudomonas]` keeps only
    that one.

    Check the names actually present before writing a filter — the `Genus` column of
    `6.taxonomy/taxon_seq_table.txt` shows them — and list every variant that applies.
    This affects gyrB only: rpoB's FROGS reference carries NCBI taxonomy, where the
    genus is not split.

## No region extractor

MetaFlux has no gyrB equivalent of Metaxa2 (16S) or ITSx (ITS), so there is nothing to trim
ASVs down to a target region. `amplicon.extraction.enabled: true` is therefore forced off
for gyrB with a warning on stderr:

```text
[MetaFlux] warning: marker gyrB has no target-region extractor; ignoring
amplicon.extraction.enabled: true (no extraction step will run)
```

Setting `enabled: false` explicitly keeps the run log clean. The ASV length filter still
runs either way, so off-target amplicons of the wrong size are still removed.

## What the output looks like

The final tables are the same three files as for every other marker, under
`6.taxonomy/` in the output directory:

| File | Rows | Notes |
|---|---|---|
| `asv_table.txt` | ASV IDs | per-sample counts plus one `taxonomy` column |
| `asv_table_seqs.txt` | ASV sequences | the same table keyed by sequence instead of ID |
| `taxon_seq_table.txt` | ASV IDs | one column per rank (`Gene` … `Species`) plus the sequence |

A gyrB taxonomy string keeps the gene tag in front of the GTDB lineage:

```text
gyrB;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas_aeruginosa
```

Ranks the classifier could not assign above `min_boot` are simply absent from the string,
so a shorter string means a less confident assignment rather than a missing rank.

## Related pages

- [rpoB](rpob.md) — the other protein-coding marker, which differs in almost every detail
  listed above except being RDP-only
- [Amplicon length and truncation](../length-and-truncation.md)
- [Keeping and discarding taxa](../taxon-filter.md)
- [Configuration](../../reference/configuration.md)

## References

**20.** Briand, M., Rué, O. & Barret, M. (2025). gyrB database for taxonomic assignment formatted for DADA2 (train_set_gyrB_v6). *Recherche Data Gouv (INRAE Dataverse)*. https://doi.org/10.57745/DD7RZ8

**22.** Parks, D. H., et al. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research*. (Genome Taxonomy Database, release r226)

**31.** Barret, M., et al. (2015). Emergence shapes the structure of the seed microbiota. *Applied and Environmental Microbiology*. (gyrB primers F64 / R353; also used to build the DD7RZ8 reference amplicons)

The numbering follows the full list on [Citation and references](../../about/citation.md).
