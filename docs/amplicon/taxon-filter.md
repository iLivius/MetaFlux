# Keeping and discarding taxa

Every amplicon run ends with a set of ASVs, each carrying a taxonomy string. Some of
those ASVs are not what the study is about: chloroplast and mitochondrial reads
co-amplified by 16S primers, plant or host sequences in an ITS run, a *Ralstonia* that
came in with the extraction kit, a *parE* paralog picked up by gyrB primers. The
taxonomy filter is where they are removed.

It is a plain keep/discard list of taxon names, set in
`amplicon.taxonomy.filter`. There is no hidden default and no per-marker magic in the
code — whatever is listed in the config is exactly what gets applied.

## Where the filter sits in the run

The filter is the last thing that happens inside the `assign_taxonomy` rule, after
classification and before any table is written.

| Step | What happens |
|------|--------------|
| Input | `5.dada2/seqs_lenfilt.fasta` (ASV sequences after the length filter, and after Metaxa2/ITSx extraction where that runs) plus the matching count table `5.dada2/seqtab_lenfilt_head_names.txt` |
| Classify | `rdp` path — DADA2 `assignTaxonomy` (+ `addSpecies` for 16S); `sintax` path — `vsearch --sintax` |
| Build lineage | ranks joined into one taxonomy string per ASV |
| **Filter** | **`keep`, then `discard`, applied to that string** |
| Output | `6.taxonomy/asv_table.txt`, `asv_table_seqs.txt`, `taxon_seq_table.txt` |

The two classifier paths are separate scripts — `workflow/scripts/amplicon/80a_assign_taxonomy.R`
for `rdp` and `workflow/scripts/amplicon/80c_parse_sintax.py` for `sintax` — but they
implement the same filter with the same behaviour, so the choice of classifier does not
change how the lists are written.

The filter removes whole ASVs. A discarded ASV disappears from all three output tables
at once, counts included, and nothing downstream re-adds it. The cost in reads is
visible in `stats/read_tracking.txt`, whose `post_taxonomy_filter` column is summed from
the filtered `asv_table.txt`: comparing it with `post_length_filter` shows how many
reads the filter removed, per sample.

!!! note

    Read classification itself is untouched. The filter only decides which of the
    classified ASVs are written out; it never changes the taxonomy an ASV was given.

## How matching works

Each ASV's taxonomy string is a lineage with ranks joined by `;` (no space), each rank
carrying its own prefix:

```text
k__…;p__…;c__…;o__…;f__…;g__…;s__…
```

A token in `keep` or `discard` is compared against **whole rank segments** of that
string — the pieces between the semicolons — not against the string as a whole. A token
matches when it is exactly equal to one of those segments.

That is what makes the lists safe to use at any rank. `g__Bacillus` matches an ASV whose
genus is exactly *Bacillus*, and nothing else: not a differently suffixed genus such as
`g__Bacillus_A`, and not a species or family segment that happens to contain the same
letters. Substring accidents cannot happen.

An ASV passes the filter as follows:

- **`keep` runs first and defines what survives.** An ASV is kept only if at least one
  of its ranks matches a `keep` token. `keep: [k__Bacteria, k__Archaea]` therefore keeps
  prokaryotes and drops both eukaryotic hits and anything left unclassified at kingdom
  level — an ASV with an empty lineage matches nothing, so it does not survive a
  non-empty `keep`.
- **`discard` then prunes what remains.** An ASV is removed if any of its ranks matches
  a `discard` token. `discard: [o__Chloroplast, f__Mitochondria]` removes the usual 16S
  plastid and host carry-over.
- Because `discard` runs second, it wins: a taxon named in both lists is removed.

## The three keys

```yaml
amplicon:
  taxonomy:
    filter:
      enabled: true
      keep:    [k__Bacteria, k__Archaea]
      discard: [o__Chloroplast, f__Mitochondria]
```

| Key | Effect |
|-----|--------|
| `enabled` | Master switch. `false` turns the filter off entirely and every classified ASV is written out, whatever `keep` and `discard` contain. |
| `keep` | Tokens that define what survives. An empty list (or the key left unset) means this direction does nothing — no ASV is dropped for failing to match. |
| `discard` | Tokens that are then removed from the survivors. An empty list means this direction does nothing. |

Both lists accept a bare scalar as well as a YAML list, so `keep: k__Fungi` and
`keep: [k__Fungi]` are equivalent.

Two consequences worth spelling out, because they are easy to conflate:

- `enabled: true` with both lists empty is a working filter that removes nothing. That
  is the fail-safe state, not an error.
- `enabled: false` with populated lists ignores those lists completely. It is the quick
  way to test what the filter was actually costing.

### Why the shipped config has the values commented out

`config/config.yaml` is always loaded as a base layer underneath the `--configfile`
given on the command line. Any concrete value written there would silently apply to
every run that does not override it — 16S tokens leaking into an 18S run, for instance,
where `keep: [k__Bacteria, k__Archaea]` would discard every eukaryotic ASV. So the
template ships the recommended values commented out:

```yaml
      keep:    #[k__Bacteria, k__Archaea]
      discard: #[o__Chloroplast, f__Mitochondria]
```

Left unset they resolve to `[]`, which filters nothing. Copy the line for the marker in
use into the run config and uncomment it there.

## Recommended starting points per marker

These are starting points, not defaults — nothing in the code applies them.

| Marker | `keep` | `discard` |
|--------|--------|-----------|
| 16S | `[k__Bacteria, k__Archaea]` | `[o__Chloroplast, f__Mitochondria]` |
| ITS | `[k__Fungi]` | `[]` |
| 18S | `[d__Eukaryota]` — note this also drops ASVs left unclassified at Domain, which on some 18S datasets is a large fraction; see [18S rRNA](markers/18S.md) before applying it | `[]` |
| gyrB | *no recommended default* — see below | |
| rpoB | *no recommended default* — see below | |

**18S** uses `d__Eukaryota`, not `k__Eukaryota`: the taxonomy reference is PR2, whose
top rank is Domain and is emitted with the `d__` prefix. The config template notes that
a SILVA-Eukaryotic taxonomy reference would carry `k__` instead.

**gyrB** has no recommended keep/discard, but it does have a filter worth applying, of a
different kind. The DD7RZ8 reference tags every entry's first rank as `gyrB` (the true
gene) or `other` (a co-amplified paralog such as *parE*). That tag rides along as an
ordinary taxonomy rank, so the same machinery drops paralogs with no extra code:

```yaml
discard: [other]     # drops paralogs, keeps ASVs the classifier left unassigned
```

`keep: [gyrB]` is stricter, not equivalent — it also drops every ASV that came back
unclassified, since an empty lineage matches no keep token.

Note the tag has no rank prefix — it is emitted verbatim as `other` or `gyrB`. How much
`other` to expect depends on the primers and the samples, so a sensible first run leaves
the lists empty, looks at the paralog fraction in the output, and decides after. See
[gyrB](markers/gyrb.md).

**rpoB** has nothing structurally equivalent. The FROGS RefSeq reference contains only
rpoB, so there is no paralog or organellar class to filter at this stage. Leave the
lists empty unless a particular run shows a reason not to. See [rpoB](markers/rpob.md).

## Writing tokens that actually match

Tokens are compared literally, so they have to be spelled exactly as the reference
database spells them. The reliable way to get that right is to copy the name from the
`taxonomy` column of `6.taxonomy/asv_table.txt` from a first, unfiltered run.

The rank prefixes for the Linnaean seven ranks (16S, ITS, gyrB from phylum down, rpoB):

| Rank | Prefix |
|------|--------|
| Kingdom | `k__` |
| Phylum | `p__` |
| Class | `c__` |
| Order | `o__` |
| Family | `f__` |
| Genus | `g__` |
| Species | `s__` |

18S is the exception, because PR2 has more ranks: Domain `d__`, Supergroup `sg__`,
Division `dv__`, Subdivision `sbd__`, then Class, Order, Family, Genus, Species with the
usual prefixes. The `sintax` path for 18S merges Division and Subdivision into a single
`dv__` field, so a subdivision-level token only exists on the `rdp` path.

Species-level tokens need extra care, because the exact form follows the reference and
not the pipeline: SILVA (16S) via `rdp` yields a binomial with a space, e.g.
`s__Bacillus subtilis`, while UNITE (ITS), the GTDB-derived gyrB reference and the
NCBI RefSeq-derived rpoB reference use an underscore, e.g. `s__Staphylococcus_simiae`.
The `sintax` path emits whatever label its database carries. Copy, do not retype.

!!! tip

    To drop lab or kit contaminants that a negative control revealed, add them to
    `discard` at the rank that identifies them cleanly — usually genus:

    ```yaml
    discard: [o__Chloroplast, f__Mitochondria, g__Ralstonia, g__Bradyrhizobium]
    ```

## When the filter empties the table

If the filter removes 100% of ASVs, the run stops with an error instead of writing a
header-only table. In practice that outcome almost always means the lists do not match
the marker — 16S keep tokens on an ITS run, or `k__Eukaryota` on an 18S run classified
against PR2, where the top rank is `d__`. A silent empty table would be discovered much
later, in R, with no clue as to why.

The error names the marker and the number of ASVs lost (numbers below are illustrative):

```text
the keep/discard filter removed ALL 1843 ASVs (amp_type=ITS). The keep/discard lists
likely do not match this marker: for ITS use keep: [k__Fungi]; for 16S use keep:
[k__Bacteria, k__Archaea]. Fix amplicon.taxonomy.filter, or set
amplicon.taxonomy.filter.enabled: false.
```

The check only fires when the filter is enabled and there were ASVs to begin with, so a
run that produced no ASVs at all fails earlier and for its own reasons.

What to check, in order:

1. The rank prefix. Does the marker's reference use `k__` or `d__` at the top?
2. The marker itself — is `amplicon.type` what the tokens were written for?
3. The spelling, against the `taxonomy` column of a run made with `enabled: false`.

## Seeing what the filter did

Both classifier paths log their counts to `logs/assign_taxonomy.log`, three lines that
say exactly how many ASVs each direction cost:

```text
[assign_taxonomy] Keep filter (k__Bacteria,k__Archaea): 1620 / 1843 ASVs retained
[assign_taxonomy] Discard filter (o__Chloroplast,f__Mitochondria): 47 ASV(s) removed
[assign_taxonomy] Final ASV count: 1573
```

The `sintax` path writes the same three lines under an `[assign_taxonomy/sintax]` tag.
A `keep` line that retains far fewer ASVs than expected is worth following up before
trusting the table: often the reference simply failed to classify those ASVs at the rank
the token names, rather than the ASVs being genuinely off-target.

!!! warning "Retired configuration keys"

    Earlier versions used regular expressions, `include_pattern` and `exclude_pattern`,
    which carried per-marker defaults. Those keys are gone. A config that still contains
    either of them stops the run at startup with a message showing the keep/discard
    equivalent, rather than running with the filter silently disabled.

## The shotgun counterpart

Shotgun runs have their own filter, `shotgun.taxonomy_filter`, with the same three keys
and exactly the same whole-segment token matching: `keep` first, then `discard`, empty
list means no-op. It is applied to the final OTU table after Bracken, so it changes
which OTUs appear in `03.abundance/otu_table.tsv` and `otu_table.biom`, never how reads
were classified.

Two differences from the amplicon side:

- It ships **off** (`enabled: false`). Metagenomic sample types vary too much for a
  universal starting point; common uses are `keep: [k__Bacteria, k__Archaea]` for
  prokaryotes only, or `discard: [g__Homo]` as an explicit host safety net.
- If it removes every OTU, it logs a warning and writes an empty table rather than
  failing. The amplicon filter is stricter because an amplicon run has exactly one
  marker to get right, and getting it wrong is the common cause of an empty result.

See the [shotgun overview](../shotgun/overview.md) and the
[configuration reference](../reference/configuration.md).
