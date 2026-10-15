# rpoB (bacterial RNA polymerase subunit B)

`rpoB` encodes the beta subunit of bacterial RNA polymerase. Like gyrB it is a
protein-coding marker, single-copy in almost all bacteria, and it assigns a far higher
proportion of ASVs to species than a 16S V3–V4 amplicon does. Ogier et al. (2019) put it
forward, together with the degenerate `Univ_rpoB_deg` primer pair, as a marker for
describing bacterial community diversity by amplicon sequencing, and report 60–80% of rpoB
OTUs assigned at species level against about 20% for 16S. It does not resolve strains: the
reference stops at species, so neither does the output.

The MetaFlux path for rpoB is the ordinary amplicon path — PhiX removal, primer trimming,
DADA2 denoising, length filtering, taxonomy. What is marker-specific lives in
`workflow/markers/rpoB.yaml`: the reference is a set of full-length genes rather than
pre-cut amplicons, it needs a one-time reformatting before DADA2 can read it, and only the
RDP classifier is available.

## At a glance

| | |
|---|---|
| `amplicon.type` | `rpoB` (also accepted as `rpob` or `RPOB` — matched case-insensitively) |
| Reference database | FROGS RefSeq rpoB, build 20240707 (complete + chromosome NCBI RefSeq genomes) |
| Role of that database | both the taxonomy reference and the amplicon-length probe substrate |
| Reference content | full-length rpoB genes, roughly 4 kb, about 47,000 sequences |
| Probe mode | `pcr` — in-silico PCR recovers the amplicon window |
| Preparation | download + unpack, then a one-time header conversion to DADA2 format |
| Classifier | RDP only (`amplicon.taxonomy.method: rdp`) |
| Separate species step | none — the trainset itself carries species |
| Region extractor | none (no rpoB equivalent of Metaxa2 or ITSx) |
| Rank model | `Domain;Phylum;Class;Order;Family;Genus;Species` (no kingdom) |
| Primers | none are shipped; Ogier 2019 `Univ_rpoB_deg` F/R is the usual pair |

## Setting up a run

Set the marker type, point at the primer FASTAs, and choose the RDP classifier:

```yaml
mode: amplicon

amplicon:
  type: rpoB
  primers:
    fwd: /path/to/Univ_rpoB_deg_F.fasta
    rev: /path/to/Univ_rpoB_deg_R.fasta
    orientation: fixed          # fixed | mixed
  expected_length: auto
  probe_length_stat:
    rpoB: p95                   # required for rpoB, auto or not
  extraction:
    enabled: false              # no rpoB region extractor
  taxonomy:
    method: rdp                 # rpoB has no SINTAX build
    min_boot: 80
    try_rc: true
    filter:
      enabled: true
      keep:    []
      discard: []               # no recommended default for rpoB
```

Primers are never supplied by MetaFlux: `amplicon.primers.fwd` and `amplicon.primers.rev`
must both point at existing FASTA files, and parsing stops immediately if either is
missing. `Univ_rpoB_deg` is degenerate, which is what makes it broad enough for community
work. The same primer FASTAs are used both to trim the reads and to run the in-silico PCR
described below, and both steps allow the mismatch rate set by
`amplicon.cutadapt.max_error_rate` (default `0.2`).

!!! warning "The shipped config template defaults to `method: sintax`"

    `config/config.yaml` ships with `amplicon.taxonomy.method: sintax`. That is fine for
    16S/ITS/18S but not available for rpoB, so a run config copied from the template
    without changing that line stops at parse time, before any compute is spent:

    ```text
    [MetaFlux] marker rpoB has no SINTAX reference database, so
    amplicon.taxonomy.method: sintax is not available for it. Use method: rdp.
    ```

## The reference database

The reference is the FROGS RefSeq rpoB databank built on 2024-07-07 from complete and
chromosome-level NCBI RefSeq genomes — the rpoB genes extracted from those genomes, about
47,000 sequences, each a full-length gene of roughly 4 kb.

| | |
|---|---|
| Downloaded from | `https://web-genobioinfo.toulouse.inrae.fr/frogs_databanks/assignation/rpoB/rpoB_bacteria_NCBI_refseq_genome_complete_and_chromosome_20240707.tar.gz` |
| Fetch rule | `fetch_rpob_archive` → `<out_dir>/logs/refdb/fetch_rpob_archive.log` |
| Conversion rule | `convert_rpob_to_dada2` → `<out_dir>/logs/refdb/convert_rpob_to_dada2.log` |
| Cached at | `refdb/rpob/rpob_refseq_cc_20240707_domain-d_dada2.fa.gz` |

The cache root is `references.refdb_root` (default `refdb`). An entry under `references`
is only needed to relocate a database that already exists elsewhere on disk — for rpoB that
entry names the **converted** file:

```yaml
references:
  rpob:
    frogs: /shared/dbs/rpob_refseq_cc_20240707_domain-d_dada2.fa.gz
```

The converted file is one MetaFlux builds itself, so a path given here is where
MetaFlux writes it, and it rebuilds it there whenever the conversion script is newer
than the file (see
[reference files MetaFlux builds itself](../../reference/configuration.md#reference-files-metaflux-builds-itself)).

Both rules run once, the first time an rpoB run needs the database, and are skipped
afterwards because their output is already on disk. "domain-d" in the file name marks
the file that writes the domain as `d__` (below). The converted file of MetaFlux v2.4.0
and earlier, `rpob_refseq_cc_20240707_dada2.fa.gz`, wrote it as `k__` and is no longer
read, so the first rpoB run of a newer MetaFlux downloads the FROGS archive again
(173 MB) and converts it once.

### Why a conversion step exists

FROGS distributes the databank as a `.tar.gz` bundling three things: the FASTA, a BLAST
index, and an RDP-classifier tree. `fetch_rpob_archive` downloads and unpacks it, keeps
only the FASTA (the BLAST index files share the `.fasta` stem, so the rule picks the one
true FASTA and skips the `.properties` file), and gzips it to a temporary
`refdb/rpob/rpob_frogs_raw.fasta.gz`. Snakemake deletes that intermediate once the
conversion has consumed it.

The FASTA is shaped for FROGS' own affiliation tool, not for DADA2. Every header is an
accession followed by a lineage under a constant `Root`, in which each node carries a
`[id: N]` database tag:

```text
>WP_095092576.1 Root;k__Bacteria [id: 1];p__Bacillota [id: 2];...;s__Staphylococcus_simiae [id: 7]
```

DADA2's `assignTaxonomy` instead wants the header to *be* the lineage — semicolon
separated, no accession, no per-node tags. `convert_rpob_to_dada2` runs
`workflow/scripts/refdb/frogs_rpob_to_dada2.py`, which makes exactly three edits per
header and leaves sequence lines untouched:

1. drop the leading `>ACCESSION Root;`, so the lineage starts at `k__Bacteria`;
2. strip every ` [id: N]` tag;
3. write the first rank as a domain: `k__Bacteria` becomes `d__Bacteria`.

The third edit is needed because FROGS writes the domain under the `k__` prefix. This
release carries NCBI's taxonomy as it was when FROGS built it (July 2024), before NCBI
added kingdoms for prokaryotes, so every lineage has a domain and no kingdom. In
MetaFlux's tables `d__` is always the domain and `k__` always the kingdom, so the rpoB
pack names its first rank Domain, with the prefix `d__`, and the converted file says the
same. The result is the file DADA2 trains on:

```text
>d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus_simiae
```

Because sequence lines pass through verbatim, nothing that faces the ASVs is altered — the
conversion only rewrites labels.

!!! note "The converter refuses to run on an unexpected header shape"

    All three edits are anchored to the exact shape of the 2024-07-07 release, and the
    URL in the marker pack pins that dated filename rather than a "latest" pointer. If a
    future release changed the header format, an edit that silently did nothing would
    ship corrupted taxonomy strings — leftover accession or `[id:` text mistaken for a
    rank value, or a kingdom written as a domain — with no error at all. So each edit is
    checked: a header that does not start with `>ACCESSION Root;`, that still contains
    `[id:` afterwards, or whose first rank is then anything but `k__Bacteria` or
    `k__Archaea` (a kingdom such as `k__Pseudomonadati`, or a rank already written
    `d__`), stops the run and names the offending line. The log records how many headers
    were rewritten. The [reference check](index.md#the-reference-check-reference_values)
    then confirms that every lineage of the converted file starts with `d__Bacteria`.

## Amplicon length: recovered by in-silico PCR

The rpoB reference is full-length, so reading its sequence lengths would return the length
of the gene (about 4 kb), not of the amplicon. The marker pack therefore sets
`probe_mode: pcr`, the same path used for 16S and 18S: two passes of cutadapt against the
reference, using the supplied primers. For rpoB a third pass follows.

| Pass | What it does |
|---|---|
| 1 | `-g file:fwd` with `--discard-untrimmed` — keep references carrying the forward primer at the 5′ end, and trim it off |
| 2 | `-a file:rev_rc` with `--discard-untrimmed` — of those, keep the ones carrying the reverse complement of the reverse primer at the 3′ end, and trim it off |
| 3 | the forward primer again, required at its full length, with `--discard-trimmed` — leave out of the statistics the bodies that still hold a second forward-primer site |

The reverse primer is reverse-complemented first (with `seqtk seq -r`) so it can be
searched for at the 3′ end. What survives both passes is the set of in-silico amplicon
bodies; their lengths become the probe distribution. With the `Univ_rpoB_deg` primers the
recovered bodies are around 387 bp — the amplicon minus the two primers, which cutadapt
trims off in the two passes above.

Pass 3 is set in the rpoB marker pack (`probe_exclude_second_site: true`). In 17% of the
FROGS genes that both primers amplify, the forward primer also matches, with as many
mismatches as the probe allows, a stretch about 920 bp upstream of its real site, and
cutadapt keeps that earlier match. Without pass 3 the probe gives a median of 389 bp but a
p95 of 1,332 bp, and p95 is the statistic `expected_length: auto` uses; with it, q1 is
386 bp, the median 386 bp and p95 395 bp. The left-out bodies stay in the amplicons file.
See [A primer that matches twice](../length-and-truncation.md#a-primer-that-matches-twice).

This is the opposite arrangement to [gyrB](gyrb.md), whose reference ships already
primer-trimmed and is therefore measured directly. The difference is a property of the two
databases, not of the two genes.

!!! tip "If the probe finds nothing"

    When no reference sequence carries both primers, the probe stops rather than reporting
    an empty distribution:

    ```text
    [amplicon_probe] ERROR: no reference sequences had both primers. Check primer
    sequences against the reference DB orientation and primer error tolerance.
    ```

    That points at the primer FASTAs — the wrong orientation for this reference, or a
    mismatch tolerance too tight for the primers in use. The probe log is
    `<out_dir>/logs/amplicon_probe.log`.

The probe only runs when `amplicon.expected_length: auto`. It writes its statistics, and
the recovered amplicons, to files whose names encode the marker, the reference tag (which
ends in `_excl-second-site` for pass 3), and a 12-character hash of both primer FASTAs:

```text
refdb/cache/probe_rpoB_rpob_refseq_cc_20240707_excl-second-site_<primer_hash>.json
refdb/cache/probe_rpoB_rpob_refseq_cc_20240707_excl-second-site_<primer_hash>.amplicons.fa.gz
```

Changing primers changes the hash, so the probe reruns; keeping the same primers across
runs reuses the cached result. The JSON then feeds two downstream steps: `pick_trunclen`,
which uses the resolved expected length to enforce the forward/reverse overlap constraint
on `truncLen`, and the ASV length filter, which builds its `auto` window from the probe's
q1 and p95 plus `length_filter.window_margin`. See
[Amplicon length and truncation](../length-and-truncation.md) for both.

!!! warning "`probe_length_stat.rpoB` must always be present"

    Which statistic of the probe distribution becomes `expected_length` is chosen by
    `amplicon.probe_length_stat.rpoB` (the shipped template sets `p95`). The key is
    checked at parse time whether or not `expected_length` is `auto`; with no such entry,
    parsing stops with:

    ```text
    [MetaFlux] amplicon.probe_length_stat is missing the 'rpoB' key that marker rpoB
    needs. Add it, e.g.
        probe_length_stat:
          rpoB: p95
    ```

    A manual `expected_length` (a single integer, or a `[min, max]` pair) skips the probe
    itself, but the entry must still be there.

## Taxonomy: RDP only

The FROGS package targets FROGS' own classifier; MetaFlux uses the converted FASTA with
DADA2's `assignTaxonomy` instead. There is no SINTAX-formatted build of this reference, so
the marker pack leaves `taxonomy_sintax_db` empty and MetaFlux refuses `method: sintax` at
parse time rather than failing later inside the classifier.

`taxonomy.min_boot` sets the bootstrap confidence minimum and `taxonomy.try_rc` decides
whether reverse-complemented sequences are also tried. There is no separate `addSpecies`
step: unlike 16S, where species come from a second SILVA file, the converted rpoB trainset
is already labelled to species, so the classifier reaches the species rank directly.
`amplicon.seed` fixes the bootstrap RNG, so RDP results are reproducible run to run,
apart from DADA2's random tie-breaking between equally good reference lineages, which can
move a bootstrap value between identical runs (see
[Troubleshooting](../../troubleshooting.md)).

The classifier writes seven rank columns, taken from the marker pack. The converted
reference already carries `d__`/`p__`/… prefixes, so the pack sets
`prefix_style: embedded` and MetaFlux emits each value verbatim rather than prepending
its own:

| Slot | Rank column | Prefix in the taxonomy string |
|---|---|---|
| 1 | `Domain` | `d__` |
| 2 | `Phylum` | `p__` |
| 3 | `Class` | `c__` |
| 4 | `Order` | `o__` |
| 5 | `Family` | `f__` |
| 6 | `Genus` | `g__` |
| 7 | `Species` | `s__` |

There is no kingdom slot, so a `k__` filter token stops the run at start-up, and so does
a domain other than `d__Bacteria` (the reference holds bacteria only).

Names follow NCBI's taxonomy as it was when FROGS built this release (July 2024): the
phyla carry their current names, `p__Bacillota` and `p__Pseudomonadota`, and NCBI had not
yet added the prokaryote kingdoms. The species slot is an underscore-joined binomial,
`s__Staphylococcus_simiae`. That matters when writing filter tokens: they must match the
reference's spelling, not the older Firmicutes/Proteobacteria names.

These are NCBI names, not GTDB ones. Unlike the GTDB-based gyrB reference, rpoB genus names
carry no GTDB suffixes: it is `g__Pseudomonas`, never `g__Pseudomonas_E`. The two markers
sit on different taxonomic backbones — gyrB on GTDB r226, rpoB on NCBI RefSeq — which
circumscribe genera and name species differently, so labels are not interchangeable between
them.

## No paralog tag, and no default taxon filter

gyrB's reference carries a gene tag that separates true gyrB from co-amplified paralogs
such as *parE*. rpoB has no close paralog of that kind, so the FROGS RefSeq set contains
rpoB and nothing else, there is no equivalent tag to filter on, and MetaFlux recommends no
default `keep`/`discard` list for this marker.

The flip side is worth knowing: because every reference sequence is rpoB, an off-target or
paralogous amplicon cannot show up as a distinct label — it is assigned to whichever rpoB
lineage it most resembles. gyrB's tag makes that visible; for rpoB the checks are the ASV
length window, the bootstrap confidences, and the fraction of ASVs left unassigned.

The filter is still available, and the usual reason to reach for it on an rpoB run is
dropping lab or kit contaminants that a negative control revealed:

```yaml
amplicon:
  taxonomy:
    filter:
      enabled: true
      keep:    []
      discard: [g__Ralstonia, g__Bradyrhizobium]
```

Tokens are matched against whole segments of the taxonomy string, so `g__Ralstonia` hits a
genus named exactly that and never a substring elsewhere in the lineage. Use the exact
rank-prefixed name as it appears in the `taxonomy` column of `asv_table.txt`. See
[Keeping and discarding taxa](../taxon-filter.md).

!!! note "An empty list is not the same as a disabled filter"

    An empty `keep` or `discard` list makes that direction a no-op, which is the
    recommended starting point for rpoB. `enabled: false` switches the whole filter off.
    Either way every ASV survives — but if a filter ever removes *all* of them,
    `assign_taxonomy` stops with an error instead of writing a header-only table, since
    that is nearly always a marker/config mismatch.

## No region extractor

MetaFlux has no rpoB equivalent of Metaxa2 (16S) or ITSx (ITS), so there is nothing to trim
ASVs down to a target region. `amplicon.extraction.enabled: true` is therefore forced off
for rpoB with a warning on stderr:

```text
[MetaFlux] warning: marker rpoB has no target-region extractor; ignoring
amplicon.extraction.enabled: true (no extraction step will run)
```

Setting `enabled: false` explicitly keeps the run log clean. The ASV length filter still
runs either way, so off-target amplicons of the wrong size are still removed.

## What the output looks like

The final tables are the same four files as for every other marker, under
`06.taxonomy/` in the output directory:

| File | Rows | Notes |
|---|---|---|
| `asv_table.txt` | ASV IDs | per-sample counts plus one `taxonomy` column |
| `asv_table_seqs.txt` | ASV sequences | the same table keyed by sequence instead of ID |
| `taxon_seq_table.txt` | ASV IDs | one column per rank (`Domain` … `Species`) plus the sequence; the cells hold the bare names (`Pseudomonadota`, `Xenorhabdus_nematophila`), while the `taxonomy` strings keep the prefixes |
| `asv_table_unfiltered.txt` | ASV IDs, every ASV that reached the filter | the counts and `taxonomy` before the taxon filter, plus `metaxa2_origin` (`NA` for rpoB) and `filter_status` |

An rpoB taxonomy string reads like a 16S one on SILVA 138.2 (a domain and no kingdom),
with NCBI names:

```text
d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas_fluorescens
```

Ranks the classifier could not assign above `min_boot` are simply absent from the string,
so a shorter string means a less confident assignment rather than a missing rank.

## Related pages

- [gyrB](gyrb.md) — the other protein-coding marker, pre-trimmed reference and a paralog
  tag, but likewise RDP-only
- [Amplicon length and truncation](../length-and-truncation.md)
- [Keeping and discarding taxa](../taxon-filter.md)
- [Configuration](../../reference/configuration.md)

## References

**21.** FROGS rpoB reference databank (2024). Bacterial rpoB genes from complete and chromosome NCBI RefSeq genomes, build 20240707. *INRAE Toulouse*. https://web-genobioinfo.toulouse.inrae.fr/frogs_databanks/assignation/rpoB/

**32.** Ogier, J.-C., et al. (2019). rpoB, a promising marker for analyzing the diversity of bacterial communities by amplicon sequencing. *BMC Microbiology*. (rpoB primers Univ_rpoB_deg)

The numbering follows the full list on [Citation and references](../../about/citation.md).
