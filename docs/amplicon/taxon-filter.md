# Keeping and discarding taxa

Every amplicon run ends with a set of ASVs, each carrying a taxonomy string. Some of
those ASVs are not what the study is about: chloroplast and mitochondrial reads
co-amplified by 16S primers, plant or host sequences in an ITS run, a *Ralstonia* that
came in with the extraction kit, a *parE* paralog picked up by gyrB primers. The
taxonomy filter is where they are removed.

It is a plain keep/discard list of taxon names, set in
`amplicon.taxonomy.filter`. The keep and discard lists have no hidden default and no
per-marker magic in the code — whatever is listed in the config is exactly what gets
applied. For 16S there is
one more step that is not a list: ASVs that Metaxa2 calls mitochondria or chloroplast
are removed too ([the organelle filter](#the-organelle-filter-organelles), on by
default).

## Where the filter sits in the run

The filter is the last thing that happens inside the `assign_taxonomy` rule, after
classification and before any table is written.

| Step | What happens |
|------|--------------|
| Input | `05.dada2/seqs_lenfilt.fasta` (ASV sequences after the length filter, and after Metaxa2/ITSx extraction where that runs) plus the matching count table `05.dada2/seqtab_lenfilt_head_names.txt` |
| Classify | `rdp` path — DADA2 `assignTaxonomy` (+ `addSpecies` for 16S); `sintax` path — `vsearch --sintax` |
| Build lineage | ranks joined into one taxonomy string per ASV |
| **Filter** | **`keep`, then `discard`, applied to that string; then `discard`'s best-hit check; for 16S, then Metaxa2's organelle call** |
| Output | `06.taxonomy/asv_table.txt`, `asv_table_seqs.txt`, `taxon_seq_table.txt` (filtered), and `asv_table_unfiltered.txt` (every ASV, with what the filter did to it) |

The two classifier paths are separate scripts — `workflow/scripts/amplicon/80a_assign_taxonomy.R`
for `rdp` and `workflow/scripts/amplicon/80c_parse_sintax.py` for `sintax` — but they
implement the same filter with the same behaviour, so the choice of classifier does not
change how the lists are written.

The filter removes whole ASVs. A discarded ASV disappears from all three output tables
at once, counts included, and nothing downstream re-adds it. The one file that still
lists it is [`asv_table_unfiltered.txt`](#the-table-with-everything-asv_table_unfilteredtxt),
which exists for exactly that. The cost in reads is
visible in `stats/read_tracking.txt`, whose `post_taxonomy_filter` column is summed from
the filtered `asv_table.txt`: comparing it with `post_length_filter` shows how many
reads the filter removed, per sample.

!!! note

    Read classification itself is untouched. The filter only decides which of the
    classified ASVs are written out; it never changes the taxonomy an ASV was given.

## How matching works

Each ASV's taxonomy string is a lineage with ranks joined by `;` (no space), each rank
carrying its own prefix. Which ranks a string has follows the marker's reference: a 16S
string on SILVA 144, for example, reads

```text
d__…;k__…;p__…;c__…;o__…;f__…;g__…;s__…
```

while on SILVA 138.2 it has no `k__` slot and on ITS no `d__` slot. The slots of every
marker are listed under
[taxonomy strings in the final tables](../reference/output.md#taxonomy-strings-in-the-final-tables).

A token in `keep` or `discard` is compared against **whole rank segments** of that
string — the pieces between the semicolons — not against the string as a whole. A token
matches when it is exactly equal to one of those segments.

That is what makes the lists safe to use at any rank. `g__Bacillus` matches an ASV whose
genus is exactly *Bacillus*, and nothing else: not a differently suffixed genus such as
`g__Bacillus_A`, and not a species or family segment that happens to contain the same
letters. Substring accidents cannot happen.

One exception is deliberate. SILVA 144 marks some clades `X--other` (a clade that does
not contain the type of taxon *X* but is still considered part of it), and the filter
ignores that suffix on both sides: `g__Bradyrhizobium` also matches
`g__Bradyrhizobium--other`. Labels that SILVA 144 merges with `+`
(`g__Enterobacter--other+Klebsiella--other`) are not split; to match one, use the whole
label as the token. See [16S rRNA](markers/16S.md#the-other-suffix-amplicontaxonomysilva_other).

An ASV passes the filter as follows:

- **`keep` runs first and defines what survives.** An ASV is kept only if at least one
  of its ranks matches a `keep` token. `keep: [d__Bacteria, d__Archaea]` therefore keeps
  prokaryotes and drops both eukaryotic hits and anything left unclassified at domain
  level — an ASV with an empty lineage matches nothing, so it does not survive a
  non-empty `keep`.
- **`discard` then prunes what remains.** An ASV is removed if any of its ranks matches
  a `discard` token. `discard: [o__Chloroplast, f__Mitochondria]` removes the usual 16S
  plastid and host carry-over.
- **`discard` also checks the best hit** (`discard_min_support`, below): an ASV is
  removed when the classifier's best hit carries a `discard` token at that token's rank
  with enough support, even if that rank fell below the confidence cutoff and is
  missing from the reported lineage.
- **For 16S, Metaxa2's organelle call comes last** (`organelles`, below): an ASV that
  Metaxa2 called mitochondria or chloroplast is removed even when no token matched it.
- Because `discard` runs second, it wins: a taxon named in both lists is removed.

## The five keys

```yaml
amplicon:
  taxonomy:
    filter:
      enabled: true
      keep:    [d__Bacteria, d__Archaea]
      discard: [o__Chloroplast, f__Mitochondria]
      discard_min_support: 50
      organelles: metaxa2
```

| Key | Effect |
|-----|--------|
| `enabled` | Master switch. `false` turns the filter off entirely and every classified ASV is written out, whatever `keep` and `discard` contain. |
| `keep` | Tokens that define what survives. An empty list (or the key left unset) means this direction does nothing — no ASV is dropped for failing to match. |
| `discard` | Tokens that are then removed from the survivors. An empty list means this direction does nothing. |
| `discard_min_support` | Best-hit check for `discard`, in percent (0–100). Default `50`; `null` switches it off. See below. |
| `organelles` | 16S only: `metaxa2` (default) also removes the ASVs Metaxa2 calls mitochondria or chloroplast; `off` does not. Does nothing for the other markers. See below. |

Both lists accept a bare scalar as well as a YAML list, so `keep: k__Fungi` and
`keep: [k__Fungi]` are equivalent.

Two consequences worth spelling out, because they are easy to conflate:

- `enabled: true` with both lists empty is a working filter that removes nothing, except
  on 16S, where the organelle filter still removes the ASVs Metaxa2 calls mitochondria or
  chloroplast (`organelles: off` stops that too). That is the fail-safe state, not an
  error.
- `enabled: false` with populated lists ignores those lists completely, and switches
  the organelle filter off as well. It is the quick way to test what the filter was
  actually costing. The tokens are still checked when the run starts, so that they are
  right on the run that switches the filter back on (see
  [when MetaFlux refuses a token at start-up](#when-metaflux-refuses-a-token-at-start-up)).

### Why the shipped config has the values commented out

`config/config.yaml` is always loaded as a base layer underneath the `--configfile`
given on the command line. Any concrete value written there would silently apply to
every run that does not override it — 16S tokens leaking into an 18S run, for instance,
where `keep: [d__Bacteria, d__Archaea]` would discard every eukaryotic ASV. So the
template ships the recommended values commented out:

```yaml
      keep:    #[d__Bacteria, d__Archaea]
      discard: #[o__Chloroplast, f__Mitochondria]
```

Left unset they resolve to `[]`, which filters nothing. Copy the line for the marker in
use into the run config and uncomment it there. `discard_min_support` and `organelles`
ship set. The first does nothing until `discard` lists a token; the second acts only on
16S with extraction on, and does nothing for any other marker.

### The best-hit check (`discard_min_support`)

The taxonomy each ASV is reported with stops at the last rank whose confidence reaches
`min_boot` (rdp) or `sintax_cutoff` (sintax). Below that rank the classifier still has a
best guess; it is simply not confident enough to report it. For most taxa that is the
right behaviour. For organelle 16S it lets contaminants through: a mitochondrial ASV
whose best hit is *Mitochondria* at 60% support is reported as `p__Pseudomonadota`, so
`f__Mitochondria` never sees it. Plant and algal mitochondria are diverse, SILVA 144
holds many of them, and the classifiers often stay below 80% at order and family level
for them.

So, for every ASV that survived `keep` and `discard`, each `discard` token is also
compared with the best hit's name at the token's own rank (`o__` = Order, `f__` =
Family, and so on), and the ASV is removed when the names match and the support at that
rank is at least `discard_min_support`:

| Method | Support used |
|---|---|
| `rdp` | the bootstrap support of that rank (0–100), from DADA2 `assignTaxonomy` |
| `sintax` | the confidence of that rank in the full best-hit lineage VSEARCH reports (0–1, compared ×100) |

The `keep` list is not affected, and `s__` tokens are only ever compared with the
reported lineage (the species from `addSpecies` has no support value). Each ASV removed
this way gets its own log line:

```text
[assign_taxonomy]   ASV_19: best hit f__Mitochondria at 70% support (reported as d__Bacteria;k__Pseudomonadati;p__Pseudomonadota;c__Alphaproteobacteria)
[assign_taxonomy] Discard filter, best-hit check (support >= 50%): 5 more ASV(s) removed
```

The default of 50 means the best hit needs at least half the support at that rank (50
of 100 bootstrap draws, or a SINTAX confidence of 0.5). Set it
higher to be more conservative, or to `null` to get the filter of MetaFlux v2.4.0 and
earlier.

### The organelle filter (`organelles`)

For 16S, Metaxa2 has already looked at every ASV before taxonomy is assigned. It is the
region extractor (rule `target_extract`), and after cutting out the rRNA it aligns each
ASV against its own reference set and sorts it into one origin: bacteria, archaea,
eukaryota, chloroplast, mitochondria or uncertain. With `organelles: metaxa2`, the
default, the filter also removes every ASV that Metaxa2 called mitochondria or
chloroplast.

This is for organelle 16S that no token can reach. A plant mitochondrion whose
classifier best hit is a bacterium is reported as a bacterium, so neither
`f__Mitochondria` nor the best-hit check has anything to match.

How much is left for the organelle filter depends on the classifier and the reference.
Plant mitochondrial 16S is often closest to a mitochondrion in SILVA, yet `rdp` can
report it as *Enterobacterales* on either release, and `sintax` on SILVA 138.2 may stop
at phylum, class or order with a bacterial best hit below that. `sintax` on SILVA 144
more often has *Mitochondria* as the best hit, so there `f__Mitochondria` and the
best-hit check already remove most of it. Metaxa2, which aligns each ASV against its own
reference set, organelles included, calls such ASVs mitochondrial, so with the organelle
filter on, the choice of classifier and release matters much less for which ASVs are
kept.

The tokens stay useful next to it: now and then the classifier calls an ASV
*Mitochondria* that Metaxa2 calls bacterial. The two signals overlap a lot but not
completely, so the filter uses both.

What the filter takes from Metaxa2, and what it does not:

- Only the organelle calls. Metaxa2's other calls (eukaryota, uncertain, the domains)
  remove nothing; which domains survive is still decided by `keep`.
- Only Metaxa2's final answer, from its per-origin FASTA files under `05.dada2/metaxa2/`,
  collected per ASV in `05.dada2/metaxa2_origin.tsv`. Not the preliminary letter in the
  extraction headers, which calls far too many sequences mitochondrial (see
  [16S rRNA](markers/16S.md#metaxa2-reports-origin-twice-and-only-the-second-answer-is-the-real-one)).
  That final call is stable from run to run, unlike Metaxa2's own taxonomy file,
  which the filter does not read.
- It runs after `keep`, `discard` and the best-hit check, so it only removes what those
  three let through.

It can only act where Metaxa2 ran: 16S with `extraction.enabled: true`, and
`filter.enabled: true`. For ITS, 18S, gyrB and rpoB, or with extraction off, the setting
is inactive and the run says so once at start-up, instead of failing:

```text
[MetaFlux] organelle filter: inactive — marker ITS does not run Metaxa2, so there is no organelle call to use
```

To switch it off, set `organelles: off` (`null` works too). The three filtered tables are
then exactly what the tokens and the best-hit check alone give. Any other value stops the
run at start-up.

The log lists what Metaxa2 found, how much of it the earlier steps had already removed,
and every ASV that only Metaxa2 caught (numbers illustrative):

```text
[assign_taxonomy/sintax] Organelle filter (Metaxa2 origin): 29 ASV(s), 2984 reads, called mitochondria (29) or chloroplast (0)
[assign_taxonomy/sintax]   already removed above by the tokens or the best-hit check: 26 ASV(s), 2977 reads (keep 0, discard 20, best-hit check 6)
[assign_taxonomy/sintax]   removed only by the organelle filter: 3 ASV(s), 7 reads
[assign_taxonomy/sintax]     ASV_199: mitochondria by Metaxa2, 3 reads (reported as d__Bacteria;k__Pseudomonadati;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales--other)
```

### The table with everything: `asv_table_unfiltered.txt`

Removed ASVs are not lost. Every run writes `06.taxonomy/asv_table_unfiltered.txt`: the
same layout as `asv_table.txt`, for every ASV that reached the filter, with two more
columns.

| Column | Values |
|---|---|
| `metaxa2_origin` | Metaxa2's call: `mitochondria`, `chloroplast`, `bacteria`, `archaea`, `eukaryota`, `uncertain`. `NA` where Metaxa2 did not run (other markers, or extraction off) or made no call. Filled whenever Metaxa2 ran, also with `organelles: off`. |
| `filter_status` | `kept`, `removed_keep`, `removed_discard`, `removed_best_hit` or `removed_organelle` |

An ASV that more than one step would remove is labelled with the first one, in the order
the filter runs them: keep, discard, best-hit check, organelle filter. So
`removed_organelle` means that only Metaxa2 caught it, while `metaxa2_origin` shows
Metaxa2's opinion on every row, whatever the status. The `kept` rows are `asv_table.txt`
exactly, and the read counts of all rows add up to the `post_length_filter` column of
`stats/read_tracking.txt`. The tree and the other tables are built from the kept rows
only.

This is the table for questions about the removed ASVs themselves — for example, which
primer pair amplifies fewer organelles. The organelle share of each sample, in R:

```r
tab    <- read.table("06.taxonomy/asv_table_unfiltered.txt", header = TRUE, sep = "\t", row.names = 1)
counts <- tab[, !(names(tab) %in% c("taxonomy", "metaxa2_origin", "filter_status"))]
colSums(counts[tab$metaxa2_origin %in% c("mitochondria", "chloroplast"), ]) / colSums(counts)
```

This counts what Metaxa2 calls an organelle. To count what the tokens removed as well, add
`| tab$filter_status %in% c("removed_discard", "removed_best_hit")` to the row selection —
as long as `discard` lists only organelle tokens.

## Recommended starting points per marker

These are starting points, not defaults — nothing in the code applies them.

| Marker | `keep` | `discard` |
|--------|--------|-----------|
| 16S | `[d__Bacteria, d__Archaea]` | `[o__Chloroplast, f__Mitochondria]` |
| ITS | `[k__Fungi]` | `[]` |
| 18S | `[d__Eukaryota]` — note this also drops PR2's organelle-tagged domains (`d__Eukaryota:plas`, `:mito`, `:nucl`, `:apic`, `:chro`) and the few ASVs left unclassified at Domain; see [18S rRNA](markers/18S.md) | `[]` |
| gyrB | *no recommended default* — see below | |
| rpoB | *no recommended default* — see below | |

**16S and 18S** start with the domain, written `d__` under both methods; ITS starts with
the kingdom, `k__`, because UNITE has no domain rank. For 18S the taxonomy reference is
PR2: SILVA-Eukaryotic only sizes the amplicon (the length probe) and never names an ASV,
so its rank labels play no part in the filter.

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
`taxonomy` column of `06.taxonomy/asv_table_unfiltered.txt`, which lists every ASV of a
run, the removed ones included.

The rank prefixes, and the markers whose strings have each rank:

| Rank | Prefix | Markers |
|------|--------|---------|
| Domain | `d__` | 16S, 18S, rpoB |
| Kingdom | `k__` | 16S on SILVA 144, ITS |
| Phylum | `p__` | 16S, ITS, gyrB, rpoB |
| Class | `c__` | all five |
| Order | `o__` | all five |
| Family | `f__` | all five |
| Genus | `g__` | all five |
| Species | `s__` | all five; on 16S under `rdp` only |

gyrB's first slot is its gene tag, `gyrB` or `other`, with no prefix. 18S has PR2's own
ranks between the domain and the class: Supergroup `sg__`, Division `dv__`, Subdivision
`sbd__`. The `sintax` path for 18S merges Division and Subdivision into a single `dv__`
field, so a subdivision-level token only exists on the `rdp` path.

Species-level tokens need extra care, because the exact form follows the reference and
not the pipeline: SILVA (16S) via `rdp` yields a binomial with a space, e.g.
`s__Bacillus subtilis`, while UNITE (ITS, under both methods), PR2 (18S), the
GTDB-derived gyrB reference and the NCBI RefSeq-derived rpoB reference use an
underscore, e.g. `s__Staphylococcus_simiae`. Copy, do not retype.

!!! tip

    To drop lab or kit contaminants that a negative control revealed, add them to
    `discard` at the rank that identifies them cleanly — usually genus:

    ```yaml
    discard: [o__Chloroplast, f__Mitochondria, g__Ralstonia, g__Bradyrhizobium]
    ```

## When MetaFlux refuses a token at start-up

A token that no taxonomy string of the run can contain would remove every ASV in
`keep`, or nothing in `discard`, and the run would only find out after the denoising and
the classification. So the tokens are checked when the workflow is read, before any step
runs, against the slots of the strings this run writes: the marker's ranks, for the
taxonomy method and, for 16S, the SILVA release in use. One stop lists every problem and,
where there is one, the token to write instead. The check runs also with
`enabled: false`, so a wrong token is caught when it is written, not on the run that
switches the filter on.

A token stops the run when:

| Problem | Examples |
|---|---|
| its prefix is one this run's strings never have | `k__…` on 16S with SILVA 138.2, or on rpoB; `d__…` on ITS; `sbd__…` on 18S under `sintax`; `s__…` on 16S under `sintax` |
| it has no prefix (gyrB's gene tag excepted), or no name after the prefix | `Bacteria`, `d__` |
| it puts a domain under `k__`, also where the strings have a `k__` slot | `k__Bacteria`, which stops with "write d__Bacteria" |
| its name can never stand in that slot of the marker's reference | `d__Fungi` on 16S; `k__Fungi` or `k__Bacteroidota` on 16S with SILVA 144; `sg__Eukaryota` on 18S; `Other` on gyrB; `d__Archaea` on rpoB |
| it is not text: YAML reads a bare `yes`, `no`, `null` or a number as something else | `keep: [yes]` |

The last-but-one row uses the marker pack's `reference_values`, the rules the reference
file itself has to pass before any classification (see
[the reference check](markers/index.md#the-reference-check-reference_values)). A trailing
`--other` is dropped from a name before the check, as the filter itself drops it. The
rules cover the top ranks only, and an `ends_with` rule checks only the ending: on SILVA
144 any name ending in `ati` passes as a kingdom, a misspelt `k__Bacilati` included.

A typical stop, for a 16S config written for MetaFlux v2.4.0 or earlier:

```text
[MetaFlux] these amplicon.taxonomy.filter tokens can never match a taxonomy string of this run (16S, method sintax, release 144), so the run stops before any step:
  - amplicon.taxonomy.filter.keep: 'k__Bacteria' puts a domain under k__, which holds kingdoms in this table. Bacteria is a domain: write d__Bacteria.
  - amplicon.taxonomy.filter.keep: 'k__Archaea' puts a domain under k__, which holds kingdoms in this table. Archaea is a domain: write d__Archaea.
A token must equal one whole slot of the taxonomy string, prefix included. Since v2.5.0 the domain is always written d__, and k__ is the kingdom where the reference has one.
```

Where a prefix belongs to the marker but not to this run's taxonomy method, the message
says why: an `s__` token on 16S under `sintax` gets "16S under method: sintax writes no
species", an `sbd__` token on 18S under `sintax` "Under method sintax it is part of
Division_Subdivision, written dv__".

What the check cannot see, because the names themselves are not checked there:

- a misspelt name in a slot without a rule — a genus, a family, a species — passes and
  then matches nothing;
- a 16S species written with an underscore: 16S writes `s__Escherichia coli`, with a
  space, so `s__Escherichia_coli` passes under `rdp` and matches nothing;
- on 18S, a `dv__` name in the other method's form: under `sintax` PR2 joins Division and
  Subdivision with "-" (`dv__Opisthokonta-Fungi`), so `dv__Opisthokonta` matches nothing
  there, and under `rdp` the joined form matches nothing (see
  [18S rRNA](markers/18S.md#pr2-ships-two-rank-depths)).

## When the filter empties the table

If the filter removes 100% of ASVs, the run stops with an error instead of writing a
header-only table. Tokens that can never match stop the run at start-up already (above),
so at this point the tokens fit the marker's strings but not the run: a token the marker
has but the study does not want (an 18S `d__Eukaryota` left in a 16S config, where SILVA
holds a few eukaryote records), a misspelt genus or family, or a taxon the samples simply
do not hold. A silent empty table would be discovered much later, in R, with no clue as
to why.

The error names the marker and the number of ASVs lost (numbers below are illustrative):

```text
the keep/discard filter removed ALL 1843 ASVs (amp_type=ITS). The keep/discard lists likely do not match this marker: for 16S use keep: [d__Bacteria, d__Archaea]; for 18S keep: [d__Eukaryota]; for ITS keep: [k__Fungi]. Fix amplicon.taxonomy.filter, or set amplicon.taxonomy.filter.enabled: false.
```

If the organelle filter removed some of them, the message ends with how many, and that
`organelles: off` keeps those: a 16S library that is nearly all organelle 16S (leaf or
root tissue, say) can empty the table with perfectly good lists.

The check only fires when the filter is enabled and there were ASVs to begin with, so a
run that produced no ASVs at all fails earlier and for its own reasons.

What to check, in order:

1. The marker itself — is `amplicon.type` what the tokens were written for?
2. The spelling, against the `taxonomy` column of a run made with `enabled: false`.
3. Whether the samples hold the taxon at all, in `06.taxonomy/asv_table_unfiltered.txt`
   of that run.

!!! note "Run-to-run reproducibility of the filtered set"
    With a fixed `amplicon.seed`, the confidence values the filter reads come out the
    same at any thread count (`sintax` from v2.5.0 on), and so does the filtered ASV
    set. With `rdp`, DADA2's random tie-breaking can still move a bootstrap value
    between identical runs, and now and then flip an ASV's whole lineage, so that a
    token keeps or removes that ASV from one run to the next (see
    [Troubleshooting](../troubleshooting.md#amplicon-asv-counts-or-taxonomy-calls-differ-between-identical-re-runs)).

## Seeing what the filter did

Both classifier paths log their counts to `logs/assign_taxonomy.log`, one line per
step that say exactly how many ASVs each direction cost:

```text
[assign_taxonomy] Keep filter (d__Bacteria,d__Archaea): 1620 / 1843 ASVs retained
[assign_taxonomy] Discard filter (o__Chloroplast,f__Mitochondria): 47 ASV(s) removed
[assign_taxonomy] Discard filter, best-hit check (support >= 50%): 3 more ASV(s) removed
[assign_taxonomy] Organelle filter (Metaxa2 origin): 61 ASV(s), 8210 reads, called mitochondria (40) or chloroplast (21)
[assign_taxonomy]   already removed above by the tokens or the best-hit check: 48 ASV(s), 7930 reads (keep 0, discard 45, best-hit check 3)
[assign_taxonomy]   removed only by the organelle filter: 13 ASV(s), 280 reads
[assign_taxonomy] Final ASV count: 1557
```

The organelle lines appear on 16S runs only (followed by one line per ASV that only
Metaxa2 caught, see [above](#the-organelle-filter-organelles)), and
`06.taxonomy/asv_table_unfiltered.txt` has the same information per ASV, in its
`filter_status` column. The `sintax` path writes the same lines under an
`[assign_taxonomy/sintax]` tag, except that it gives a best hit's support as a SINTAX
confidence (`best hit f__Mitochondria at 0.68 confidence`) and the threshold as
`support >= 50.0%`.
A `keep` line that retains far fewer ASVs than expected is worth following up before
trusting the table: often the reference simply failed to classify those ASVs at the rank
the token names, rather than the ASVs being genuinely off-target.

!!! warning "Retired configuration keys"

    Earlier versions used regular expressions, `include_pattern` and `exclude_pattern`,
    which carried per-marker defaults. Those keys are gone. A config that still contains
    either of them stops the run at startup with a message showing the keep/discard
    equivalent, rather than running with the filter silently disabled.

## The shotgun counterpart

Shotgun runs have their own filter, `shotgun.taxonomy_filter`, with three of the same
keys (`enabled`, `keep`, `discard`; there is no best-hit check and no organelle filter)
and exactly the same whole-segment token matching: `keep` first, then `discard`, empty
list means no-op. It is applied to the final OTU table after Bracken, so it changes
which OTUs appear in `03.abundance/otu_table.tsv` and `otu_table.biom`, never how reads
were classified.

One difference from the amplicon side: it ships **off** (`enabled: false`). Metagenomic
sample types vary too much for a universal starting point; common uses are
`keep: [d__Bacteria, d__Archaea]` for prokaryotes only, or `discard: [g__Homo]` as an
explicit host safety net.

The OTU table always has eight slots, `d__ k__ p__ c__ o__ f__ g__ s__`, on any Kraken2
database. A virus has no domain: its `d__` slot holds its realm (`d__Duplodnaviria`) or
stays empty, so `keep: [d__Bacteria, d__Archaea]` also drops every virus, and so does
`[d__Bacteria, d__Archaea, d__Eukaryota]`, together with anything else without a domain
(plasmids). To keep viruses, list their realms. A virus that NCBI places in no realm has
an empty `d__`, so no realm token matches it: a `keep` list of realms drops it, and a
`discard` list of realms leaves it in. Name it by a lower rank it does have: its family
where there is one, otherwise its genus or species. Realm names follow the NCBI taxonomy
the database was built with (the realm called Monodnaviria in older builds is
Floreoviria in newer ones).

When Kraken2 runs, these tokens are checked at start-up as well, against the eight slots
and against the database itself, whose domains, kingdoms and realms MetaFlux reads from
its `taxo.k2d` file:

- the run stops on a token with a wrong or missing prefix, a domain under `k__`
  (`k__Bacteria`), `Viruses` under any prefix (no slot holds it; the message lists the
  database's realms instead), a realm under `k__` (`k__Duplodnaviria`, write
  `d__Duplodnaviria`) and, with the filter on, a `keep` list none of whose tokens can
  match a row of this database;
- it only warns about a well-formed token that this database gives no row: any `k__`
  token on a GTDB database, which has no kingdom rank, or a domain, realm or kingdom
  the database does not have. One config can serve several databases, so this alone
  does not stop the run.

With `profiler: sylph` the tokens are not checked, because the filter acts on
`otu_table.tsv` only and sylph's tables are written unfiltered.

Like the amplicon filter, it stops the run when it removes every OTU, instead of
writing an empty table. The error names the tokens and lists the `d__` values the table
holds, which is usually enough to spot a token that does not fit the table. A table
that is empty before the filter (no taxon at `tax_lev` in any sample) is still written,
with a warning in the log.

See the [shotgun overview](../shotgun/overview.md) and the
[configuration reference](../reference/configuration.md).
