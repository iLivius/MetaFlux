# Troubleshooting

This page covers these problems: an automatic truncation length that throws away
every read, a PhiX-removal step that stops on
the read pairs, a shotgun OTU table that comes out empty at a rank other than
species, a Bracken rank that MetaFlux refuses, a shotgun table step that stops on the
codes of the Kraken2 reports, a keep/discard token refused at start-up, a shotgun
taxon filter that removes every taxon, small differences between two runs on the same
data, a DADA2 error model that never converges, a run that stops because the reads
carry too few quality values or the samples carry different ones, `error_model: binned`
refused on continuous quality values, straight segments in the error plot of binned
data, the experimental
`custom_error_function` refused at start-up, `learn_errors` stopping in your own error
function, the environment of that function failing to build on Apple Silicon Macs or
Linux on ARM, an outdated marker pack in
`config/markers/`, a reference file that fails the
reference check, a reference download that fails its md5 or SHA-256 check, and YAML
mistakes in the config. Each entry gives the symptom, the reason behind it, and what
to change.

## Amplicon: automatic `truncLen` can drop all reads on very clean runs

!!! warning "Symptom"

    The run gets through primer trimming and QC, then stops at the DADA2
    filtering step (rule `dada_filter`) with a message that begins:

    ```text
    [dada_filter] STOP: filterAndTrim kept none of the <N> read pairs of sample <sample>, so it wrote no filtered files. The settings that remove reads: truncLen R1=<cut>, R2=<cut> (a read shorter than truncLen is dropped), ...
    ```

**Cause.** With `amplicon.trunc_len.mode: auto`, the `pick_trunclen` rule reads
the per-base quality tables that Falco wrote for the primer-trimmed reads
(`stats/falco/{sample}_R{1,2}_stripped/fastqc_data.txt`) and cuts just before the
first position whose per-base lower quartile (Q1, the 25th percentile of quality
at that cycle) falls below `amplicon.trunc_len.q_threshold` — 20 by default — so
the cut keeps the last cycle still at or above it. If quality never
drops that far, and clean MiSeq runs often stay high to the last cycle,
`truncLen` lands on the ceiling instead.

Falco reports a quality bin for every position up to the **longest** read in a file,
however few reads reach it — so a ceiling taken from the quality table alone describes
the longest reads, not typical ones. On a library where a few reads run to full length
and most stop a little short, a cut placed there discards nearly everything.

MetaFlux guards against this by default: `amplicon.trunc_len.min_read_coverage_pct`
(95) caps each sample at the length that percentage of its reads still reach, before
quality is considered at all. The guard is on by default, so if this error appears
anyway, the cause is one of the two cases below.

After primer trimming most reads end a little before the full read length, and
DADA2's `filterAndTrim` **discards any read shorter than `truncLen`** rather than
padding it. So a cut placed at the very end of the read removes almost everything,
and DADA2 has nothing left to work with.

!!! example "What it looks like on a real library without the guard"

    The 18S test dataset (PRJDB12231; see [test datasets](about/test-datasets.md))
    shows it. Its R1 quality never falls below Q20, so without the guard (MetaFlux
    v2.1.0 and earlier) `auto` puts the R1 cut on the ceiling and returns
    `truncLen = 288 / 239`. After primer trimming the reads themselves are:

    | | min | p05 | median | max |
    |---|--:|--:|--:|--:|
    | R1 | 78 | 271 | 281 | 289 |
    | R2 | 54 | 282 | 282 | 298 |

    The ceiling comes out at 288 because the strictest sample sets it: its longest R1
    read is 288 bp, while the other three reach 289. Most R1 reads (70%) stop at
    exactly 281. Since `filterAndTrim` discards anything shorter than `truncLen`, only
    9 of 908,768 pairs are long enough for the cut, and without the cap the run
    finished with 3 ASVs.

    Retention against the R1 cut, measured on this library:

    | truncLen R1 | 288 | 281 | 275 | 270 |
    |---|--:|--:|--:|--:|
    | pairs retained | 0.0% | 70.3% | 92.4% | 100.0% |

    Moving the R1 cut from 288 to 270 bp raises retention from 0% to 100%.

    With the 95% read-coverage cap the R1 ceiling comes in at 271 instead of 288, and
    `auto` returns 271 / 239 on this library with no manual intervention. It carries
    90.9% of the raw reads through to the final table, practically the same as a
    hand-picked 270 / 239.

**Why it can still happen with the guard on.** Two cases:

1. **`min_read_coverage_pct` was lowered.** The value is a guarantee about reads, so
   a *lower* percentage permits a *longer* cut — that is the direction that reintroduces
   this failure, and it is the first thing to check. Raising it does the opposite, and
   `100` is not an off switch: it collapses the ceiling onto the single shortest read in
   the file, which fails the overlap constraint instead.
2. **Read lengths are ragged enough that even the 95% length is unrepresentative.**
   Heavily pre-trimmed input does this, as do markers whose amplicons genuinely vary
   in length. `logs/pick_trunclen.log` shows it. It lists each sample's read-length
   limit (the length 95% of its reads reach, not counting reads shorter than half the
   expected amplicon length). For R1 and for R2 it then names the sample that sets the
   ceiling, gives that sample's longest read next to its limit, and prints a `NOTE`
   when the two are 10 bp or more apart. A very large gap there means the ceiling is being set
   from a distribution with a long tail, and a hand-picked value will beat it.

**Fix — set the truncation length by hand.**

1. Look at the read lengths *after* primer trimming: the
   `QC | primer-trimmed reads | {sample}_R{1,2}` entries under "Sequence Length
   Distribution" in the MultiQC report, or the per-sample Falco reports under
   `stats/falco/{sample}_R{1,2}_stripped/`.
2. Switch to `mode: manual` and set `manual_r1` / `manual_r2` a little below the
   length most reads still reach, so `filterAndTrim` stops discarding them.
3. Keep the pairs long enough to still merge:
   `truncLen_R1 + truncLen_R2 ≥ amplicon_length + min_overlap`
   (`amplicon.min_overlap` is 12 bp by default).

```yaml
amplicon:
  trunc_len:
    mode: manual
    manual_r1: 270
    manual_r2: 239
```

!!! tip "Pick the R1 and R2 cuts separately"

    Only the read whose quality stayed flat is stuck on the ceiling. In the example
    above R2 *did* degrade, so `auto` found a quality-based cut for it at 239. Only R1
    needed a new value; R2 kept the cut that `auto` found.

!!! tip "Retrying is cheap"

    Snakemake picks up where the run left off: PhiX removal, primer trimming and
    QC are already done and stay done, so only DADA2 and the steps after it run
    again. Trying a second value costs minutes, not a full re-run.

!!! note "ITS is not affected"

    For `type: ITS`, `dada_filter` overrides `truncLen` to `c(0, 0)`
    unconditionally — ITS amplicons vary in length, so a fixed cut would throw
    away genuinely short sequences. The 3′ tail is still cleaned, by `truncQ`
    and `maxEE` instead. Because the override is unconditional, `manual_r1` /
    `manual_r2` have no effect on an ITS run either.

Background on how the cut is chosen, and on the overlap constraint, is in
[Amplicon length and truncation](amplicon/length-and-truncation.md).

## Shotgun: BBDuk PhiX removal stops on the read pairs

!!! warning "Symptom"

    The shotgun run aborts at rule `decontam_phix` with one of these two messages:

    ```text
    java.lang.AssertionError:
    There appear to be different numbers of reads in the paired input files.
    The pairing may have been corrupted by an upstream process.  It may be fixable by running repair.sh.
    ```

    ```text
    java.lang.AssertionError: List size mismatch: 600 vs 601
    ```

    The first means the pairing **is** broken. The second does not: it comes from a
    newer bbmap than the one MetaFlux pins (see the end of this entry).

**Cause.** BBDuk reads the R1 and R2 files in lockstep, one record from each at
a time, and requires both files to hold the same number of records in the same
order. Read length has nothing to do with it — ragged lengths, and mates of the
same pair ending up different lengths, are processed normally. "There appear to be
different numbers of reads" means R1 and R2 contain different numbers of reads,
which usually happens when reads were filtered before deposition and orphaned mates
were left behind in one file. SRA and BioProject submissions are the common
case. MetaFlux's own fetch can produce it as well: rule `download_sra` calls
`fasterq-dump --split-files`, which writes a record to `_1` only for any spot
holding a single read (`--split-3` is the option that keeps the two files in
step).

The BBMap host-removal step (rule `decontam_host`) fails on the same input, for
the same reason and with the same message: BBMap reads the pair in lockstep too and
makes the same requirement. The message names the fix outright.

**Fix — repair the pairing.** `repair.sh` ships in the same `bbmap` conda
package MetaFlux already installs for these rules, so no extra tool is needed:

```bash
repair.sh in1=R1.fastq.gz in2=R2.fastq.gz \
  out1=fixed_R1.fastq.gz out2=fixed_R2.fastq.gz \
  outs=singletons.fastq.gz
```

It re-syncs the two files on read name and writes the leftover single-ended
reads to `outs=`. Point `input.fastq_dir` at the repaired files and run again.
Where the raw, untrimmed reads are still available they are the better starting
point: they are properly paired, and MetaFlux trims with fastp anyway, so
pre-trimmed input buys nothing.

!!! warning "`remove_phix: false` is not a fix here"

    Switching PhiX removal off only moves the failure to BBMap, which stops on
    the same files. Switching host removal off as well moves it to fastp, which
    stops on unequal record counts too. Repair the files instead.

**`List size mismatch` is a different problem.** MetaFlux pins bbmap 39.52, whose
BBDuk and BBMap never print this message. It comes from the rewritten BBDuk that
`bbduk.sh` runs in bbmap 39.76 to 40.02: its paired-read reader loses track of
mates and stops with this message, or hangs, on correctly paired reads, above all on
trimmed reads whose mates differ in length. It prints the same message on broken
pairs too, so it says nothing about the pairing, and `repair.sh` does not help. The
first lines of `logs/decontam_phix/{sample}.log` name the version that ran
(`Version 39.52` for the pinned one).

**Fix — run with MetaFlux's environment.** Start the run with `--sdm conda`, so that
Snakemake installs and uses the pinned environment (`workflow/envs/bbtools.yaml`),
and keep the pin if you edit that file. MetaFlux before v2.5.0 did not pin bbmap, so
its environment could get one of these releases.

More on host references and masking is in
[Decontamination](shotgun/decontamination.md).

## Shotgun: empty or nearly empty OTU table with `tax_lev` other than `S`

!!! warning "Symptom"

    `shotgun.bracken.tax_lev` is set to `G`, or to any rank other than `S`. The run
    finishes without an error, but `03.abundance/otu_table.tsv` has no rows, or only a
    handful, while Bracken's own `{sample}_output.txt` files list plenty of taxa.

**Cause.** Fixed in v2.5.0. Up to v2.4.0, rule `kraken_biom` ran kraken-biom with its
default settings, which assume a species-level report: it takes the clade count for
species, and from genus up to order only the reads placed directly on the taxon.
Plain Bracken 3.1 writes 0 in that column for almost every taxon above species (a bug
in Bracken itself, fixed in the 3.1p1 build that MetaFlux pins from v2.5.0 on), so
those taxa dropped out. A genus survived only when Kraken2 had placed no reads on any
of its species, and at class, phylum or domain level, which lie outside that range,
nothing survived at all. Bracken's `{sample}_output.txt` files were never affected, and
neither were species-level runs.

**Fix.** Update to v2.5.0 or later and run the same command again, from the same
MetaFlux folder. From v2.5.0 on, MetaFlux builds the table with rule
`finalize_otu_table` straight from the Bracken reports, from their clade counts, and
pins Bracken 3.1p1. In a folder made with an earlier version every step runs again,
Kraken2 included, because the output folders have new names (see
[What happens on a re-run](reference/running.md#what-happens-on-a-re-run)), so nothing
has to be forced.

A config with `tax_lev: D` stops at start-up from v2.5.0 on; see the next entry.

## Shotgun: the run refuses `tax_lev: D` or `K`

!!! warning "Symptom"

    A run with Kraken2 (`profiler: kraken2` or `both`) stops before any step starts,
    with one of these messages (wrapped here):

    ```text
    [MetaFlux] shotgun.bracken.tax_lev: D is not allowed. Viruses have no domain in
    NCBI's taxonomy (since 2025), and at this rank Bracken adds their reads, with a
    share of the reads left at the root, to the domain the Kraken2 report lists last
    (usually Archaea). Use P or a lower rank. To compare domains, add up the rows of
    otu_table.tsv by their d__ slot; that sum counts the reads Bracken placed at the
    table's rank, not every read of the domain.
    ```

    ```text
    [MetaFlux] shotgun.bracken.tax_lev: K is not allowed. Most protists have no
    kingdom, and a GTDB database has no kingdom rank at all. At this rank Bracken adds
    the reads of a taxon without a kingdom to the kingdom the Kraken2 report lists just
    before it, and on a GTDB database it stops after Kraken2. Use P or a lower rank;
    the k__ slot of otu_table.tsv gives the kingdom of each row.
    ```

**Cause.** Bracken reads a Kraken2 report from top to bottom and adds the reads of
every row below the target rank to the last row of that rank it has read. For a taxon
that has the target rank, that last row is the taxon's own ancestor, because the report
lists every branch right after its parent. A taxon without the rank is counted towards
whichever row of that rank comes before it. At domain and kingdom level whole groups of
taxa lack the rank:

- **Domain.** Viruses have no domain in NCBI's taxonomy since 2025, and their branch
  comes after all cellular organisms. Their reads, and the share of the root's reads
  that Bracken attributes to them, go to the domain listed last. Kraken2 lists sister
  rows by read count, largest first, so this is the domain with the fewest reads,
  usually Archaea. In a saltmarsh rhizosphere sample (SRR6394731, BioProject
  PRJNA422917) classified with the Kraken2 index `k2_NCBI_reference_20251007`, Archaea
  had 2,612 reads at phylum level and 35,636 at domain level, from only 71 viral reads
  and the root reads attributed to them.
- **Kingdom.** Most protists have no kingdom in NCBI's taxonomy, and GTDB has no
  kingdom rank at all. The protists' reads go to the kingdom listed just before them:
  in the same rhizosphere sample, Fungi had 7,957 reads in its phyla but 12,001 at kingdom
  level. On a GTDB database Bracken finds no kingdom and stops with "no reads found",
  after Kraken2 has run.

**Fix.** Set `tax_lev` to `S`, the default. It is the only rank where this cannot
happen, because every row below a species row is a strain of that species. At `P` to
`G` it hits only the few taxa that lack that rank, but one row can still gain several
times its own reads: in the rhizosphere sample, protists without a phylum took
Oomycota from 56 to 1,095 reads at `P`, and bacteria without an order took Nevskiales
from 1,393 to 6,237 at `O`. At these ranks, compare a surprising row with the
species-level table.

Domain and kingdom totals can be read from any of these tables: add up the rows of
`otu_table.tsv` by their `d__` or `k__` slot. That sum counts the reads Bracken placed
at the table's rank, not every read of the domain or kingdom: reads that Kraken2 left
above that rank, and that Bracken could not push down, are in no row (see
[Bracken abundance re-estimation](shotgun/overview.md#6-bracken-abundance-re-estimation)).
With `profiler: sylph` Bracken does not run, and `tax_lev` is not checked.

## Shotgun: the table step stops on a report code

!!! warning "Symptom"

    The run stops at rule `finalize_otu_table` with a message like this one (wrapped
    here):

    ```text
    [finalize_otu_table] SRR123_report.txt: Eukaryota (taxid 2759) is a domain in the
    Kraken2 database set as references.kraken_db, but the report codes it R2, not D.
    Kraken2 releases before 2.17 write reports like this (they do not know the rank name
    'domain'), and so do reports made with another database. Run Kraken2 again on
    references.kraken_db with the release MetaFlux pins (workflow/envs/kraken2.yaml), for
    example with -R kraken2.
    ```

    Its relatives say that a taxon "is coded D, but it is not a domain", that a kingdom is
    not coded `K` (or a `K` row is not a kingdom), or that one taxon "has two lineages".
    No table is written.

**Cause.** The OTU table takes its domain and kingdom slots from the rank codes of the
Bracken reports, and MetaFlux checks every `D` and `K` code against the taxonomy of the
database named in `references.kraken_db`, its `taxo.k2d` file. The two disagree when the
Kraken2 reports were written by a Kraken2 older than 2.17, which codes NCBI's domains
`R1` or `R2` (a run started without `--sdm conda` uses whatever Kraken2 is on the
`PATH`), or when they were made with another database than the one the config names
now. Reports of two databases can also give one taxon two lineages. A table built from
such reports would have an empty `d__` slot in every row, or slots filled from another
taxonomy, with no other sign that anything is wrong.

**Fix.** Run Kraken2 again with MetaFlux's pinned environment (`--sdm conda`), against
the database the config names; `-R kraken2` forces it when Snakemake sees nothing to
redo. If instead the reports are right and the config names the wrong database, point
`references.kraken_db` back at the database they were made with.

## Start-up: the run refuses a keep/discard token

!!! warning "Symptom"

    The run stops before any step starts, with a message like this one (an amplicon
    run; in a shotgun run the message names `shotgun.taxonomy_filter` and the Kraken2
    database):

    ```text
    [MetaFlux] these amplicon.taxonomy.filter tokens can never match a taxonomy string of this run (16S, method sintax, release 144), so the run stops before any step:
      - amplicon.taxonomy.filter.keep: 'k__Bacteria' puts a domain under k__, which holds kingdoms in this table. Bacteria is a domain: write d__Bacteria.
      - amplicon.taxonomy.filter.keep: 'k__Archaea' puts a domain under k__, which holds kingdoms in this table. Archaea is a domain: write d__Archaea.
    A token must equal one whole slot of the taxonomy string, prefix included. Since v2.5.0 the domain is always written d__, and k__ is the kingdom where the reference has one.
    ```

**Cause.** The taxon filter compares each token with whole slots of the taxonomy
strings, so a token that no string of the run can contain would remove every ASV or
taxon in `keep`, or nothing in `discard`, and the run would only find out at its end.
MetaFlux therefore checks the tokens when it starts, also with the filter switched off.
The usual cause is a config written for MetaFlux v2.4.0 or earlier, which wrote the
domain under `k__` (`k__Bacteria`); since v2.5.0 `d__` is always the domain and `k__`
always the kingdom. The message names the other causes too: a prefix the marker, method
or SILVA release never writes (`k__` on SILVA 138.2 or rpoB, `d__` on ITS, `s__` on 16S
under `sintax`), a name that cannot stand in that slot of the reference (`k__Fungi` on
16S), and in a shotgun run `Viruses`, or a virus realm under `k__`.

**Fix.** Write the token the message names (`d__Bacteria` for `k__Bacteria`) and run the
same command again; nothing has run yet, so nothing is lost. The rules, and what the
check cannot see, are under
[when MetaFlux refuses a token at start-up](amplicon/taxon-filter.md#when-metaflux-refuses-a-token-at-start-up).

## Shotgun: the taxon filter removed every taxon

!!! warning "Symptom"

    The run stops at rule `finalize_otu_table` with a message like this one (wrapped
    here; the numbers depend on the run):

    ```text
    [finalize_otu_table] shotgun.taxonomy_filter removed every taxon of the table
    (111 taxa; keep: [p__Firmicutes], discard: []). The d__ values in this table are:
    d__Archaea (1 taxon), d__Bacteria (107 taxa), d__Eukaryota (3 taxa). A token must
    equal one whole slot of the lineage, prefix included (d__ k__ p__ c__ o__ f__ g__
    s__). Check the tokens, or set shotgun.taxonomy_filter.enabled: false.
    ```

    No `otu_table.tsv` or `otu_table.biom` is written.

**Cause.** `shotgun.taxonomy_filter` compares each token with whole slots of the
lineage, so a token with a name the table does not have matches nothing. A `keep` list
in which no token matches removes every row, and so does a `discard` list that covers
every row. An empty table would look like a finished run, so MetaFlux stops instead, as
the amplicon filter does. Tokens that can never match any table (a wrong prefix, a
domain under `k__`) stop the run at start-up already (see
[the entry above](#start-up-the-run-refuses-a-keepdiscard-token)), but names below the
kingdom are not checked there, so a misspelt name, or one the database does not use,
gets this far: NCBI's taxonomy names the phylum `p__Bacillota`, not `p__Firmicutes`.

**Fix.** Compare the tokens with the `d__` values the message lists, or with the
`taxonomy` column of a run made with `shotgun.taxonomy_filter.enabled: false`, correct
them and run the same command again. Kraken2 and Bracken do not run again; the table is
built again from the Bracken reports already in `03.abundance/`. To keep every taxon,
switch the filter off.

## Amplicon: ASV counts or taxonomy calls differ between identical re-runs

Three steps in the amplicon path draw on a random-number generator, so two runs
on the exact same input can differ slightly:

| Step | Rule | What the random draw is for |
|------|------|-----------------------------|
| DADA2 `learnErrors(randomize = TRUE)` | `learn_errors` | Subsamples reads to build the error model |
| DADA2 `assignTaxonomy` (`method: rdp`) | `assign_taxonomy` | Bootstrap confidence |
| VSEARCH `--sintax` (`method: sintax`) | `assign_taxonomy` | Bootstrap confidence, seeded from a random source by default |

Without a fixed seed each of these draws a different sample every run, so a
handful of ASVs or taxonomy calls sitting near a decision threshold can flip
between two runs. This is expected DADA2/VSEARCH behaviour — not a bug, and not
corrupted data.

Before v2.4.0 the order of the reads also differed between runs, and the seed had no
effect on it: bowtie2 removed phiX on several threads and wrote the surviving pairs in
whatever order they finished, so every run
handed DADA2 the same reads in a different sequence (pairing was never broken, and
cutadapt and `filterAndTrim` only drop reads, they never reorder). DADA2 numbers ASVs by
decreasing abundance and keeps the input order for ties, so two ASVs with exactly the
same total count could swap IDs between runs. **Fixed in v2.4.0**: bowtie2 now runs with
`--reorder`, so every run hands DADA2 the reads in the same order and gets the same ASV
sequences and counts back. If you are reading output from an earlier version, compare
runs by sequence rather than by ASV ID. See
[Reproducibility](amplicon/phylogeny.md#reproducibility) for what this meant for the
optional phylogeny stage, which is sensitive to the order of the alignment rows.

`amplicon.seed` (42 by default) fixes the generator for all three, and the thread count
then no longer changes the result:

- **`learnErrors`: reproducible at any thread count.**
- **`--sintax`: reproducible at any thread count, from v2.5.0 on.** MetaFlux v2.5.0 runs
  VSEARCH 2.32, which gives every ASV its own random stream. Older VSEARCH releases
  shared one stream between all threads, so up to v2.4.0 the `sintax` results were
  identical only on one thread (`resources.threads.assign_taxonomy: 1`). A new VSEARCH
  release can also draw differently, as 2.32 did, so compared with a run made with an
  older MetaFlux version many confidence values move slightly, and a lineage near the
  cutoff can gain or lose a trailing rank, even with the same seed. The raw
  `_sintax_raw.tsv` lists the ASVs in the order the threads finish them, so sort it
  before comparing two runs; the tables built from it do not depend on that order.
  All of this needs a positive `amplicon.seed`: VSEARCH reads `0` as "pick a random
  seed".
- **`assignTaxonomy` (rdp): reproducible at any thread count, with one exception.**
  DADA2's `assignTaxonomy` breaks exact ties between equally good reference lineages at
  random, outside the seed, so a few bootstrap values can differ between identical runs;
  `sintax`, the default, is not affected. Now and then such a tie flips an ASV's whole
  lineage between two identical runs, even its kingdom: on the ITS test set
  ([PRJNA862334](about/test-datasets.md)), ASV_19 was classified down to the genus
  *Saccharomyces* in 2 of 11 identical runs. In the other 9 its kingdom support stayed
  below `min_boot: 80`, so it got no lineage at all and `keep: [k__Fungi]` removed it.

A few raw working files of the extraction step differ between identical runs as well,
without any effect on the tables: the summaries of Metaxa2 and ITSx carry their start
and end times, Metaxa2's own classification (`05.dada2/metaxa2/seqs.taxonomy.txt`) can
change in a few rows, and ITSx writes `05.dada2/_itsx_tmp/seqs.graph` in the order its
searches finish. Leave them out when you compare two runs (see
[Output files](reference/output.md#inside-05dada2)).

## Amplicon: DADA2's error model never converges

!!! warning "Symptom"

    The run stops at rule `check_error_model`, after the error models were learned and
    before any denoising. The terminal shows MetaFlux's message, which starts with

    ```text
    [check_error_model] STOP: the DADA2 error model of the R2 reads did not converge (R2:7ebe3bed),
    so MetaFlux stops before denoising. What to check, and how to go on with this model, follows;
    ```

    followed by Snakemake's own lines, which end with

    ```text
    Error in rule check_error_model:
        ...
        log: /data/run1/logs/check_error_model.log (check log file(s) for error details)
    ...
    WorkflowError:
    At least one job did not complete successfully.
    ```

    The whole message is in `logs/check_error_model.log` in the output directory: how
    the learning ended, what to check, and the line that accepts this model. A stop
    inside `learn_errors` itself, with a message that starts with
    `[learn_errors] STOP:` (for example when the error function fails), has its reason at
    the end of `logs/learn_errors.log`.

    Before the stop, rule `learn_errors` can run for hours. Its log gets one
    `selfConsist step` line per round, each followed by a MetaFlux line with the change
    of the model in that round; `tail -f logs/learn_errors.log` follows it. Of 194
    learning runs on real data, 186 converged within 8 rounds and the slowest took 101.
    On a COI library of 1.6 million pairs learned as one sample, a round took about 4
    minutes, so the cap of 200 rounds would take about 12 hours per read direction there.

**How the learning ends.** DADA2's `learnErrors` fits one error model per read
direction: for every quality score, how often each base is read as another base. It
repeats two steps, denoising the reads with the current error rates and estimating the
rates again from the result, until a round gives back a table that already went into a
round. That is DADA2's own stop, and MetaFlux calls such a model *converged*; the table
can repeat the one that went into the same round, or come back after a cycle of a few
rounds, when a few reads keep switching between two ASVs. MetaFlux ends the loop in two
more ways, and counts both as not converged:

- **The cap**, `amplicon.dada2.learn_errors.max_consist` (200 rounds by default).
- **The flat-rounds rule**, `max_flat_rounds` (20 by default; 0 switches it off). On a
  slow library the table can keep changing by tiny amounts for a hundred rounds or
  more. MetaFlux ends the learning when the model has stopped improving: the change in
  a round is less than 1/1000 of the change in round 1, and the lowest change so far
  fell by less than 10 % over the last 20 rounds. DADA2's author has said that ending
  the loop where the change stops decreasing, and using that model, makes no biological
  difference. The rule was calibrated by replaying four learning runs round by round,
  all with the loess fit: on a 101-round run it stopped at round 42, where denoising
  with the stopped model changed 0.02 % of the reads, and over every value of
  `max_flat_rounds` tested the largest change was 0.04 % of the reads. Six learning runs
  replayed with the binned fit agree: the rule ended one 95-round run at round 49, where
  denoising with the stopped model changed 0.007 % of the reads, and over every value
  tested the largest change was 0.02 %.

A model ended in either way stops the run, because the numbers alone cannot tell a
slow but healthy library from one with faulty reads: a COI library with reads of both
strands in one file (below) changed by about 0.003 per round from about round 4 on, like a
healthy slow one, and the flat-rounds rule ended it at round 32. DADA2 on its own does
not stop: at its limit (10 rounds by default) it prints `Self-consistency loop
terminated before convergence.` and carries on with the last table, as QIIME 2 and
nf-core/ampliseq do. MetaFlux stops instead, so that no ASV table is made from such a
model before someone has looked at it.

**Cause.** The loop keeps moving when the reads do not share one error profile, or when
it has little to settle on:

- **Reads of both strands in one file.** When some reads are reverse complements of the
  others, the sequencer's errors appear mirrored in part of the file (a true A read as G
  in the forward reads is a T read as C in the reversed ones), and the poor bases of the
  reversed reads sit at the start of the read instead of the end, where no truncation
  removes them.
- **Reads that still carry primer bases.** A primer with degenerate positions leaves the
  reads of one template differing from each other at several positions of high quality,
  which to DADA2 look either like real variants or like errors at a quality score that
  should hardly have any.
- **A folder trimmed under `orientation: mixed` by MetaFlux v2.5.0 or earlier**, which
  has both problems at once: the old second pass reverse-complemented the reversed pairs
  instead of swapping their mates, and the first pass also kept, through a chance match
  of three primer bases, pairs that were the other way round with most of a primer still
  on.
- **Binned quality scores.** NovaSeq, NextSeq and other recent Illumina instruments
  write only a few quality values, and DADA2's usual fit expects a smooth range of them.
  The slowest run seen, 101 rounds, was such data (three values: 9, 20 and 34), learned
  with that usual fit. With `error_model: auto`, the default, MetaFlux learns binned data
  with DADA2's binned fit instead; with `error_model: loess` they keep the usual fit.
  Whether the binned fit settles faster on the slow libraries is not measured yet.
- **Very diverse samples with few repeated reads**, such as soil or sediment. DADA2
  tells an error from a real variant by how often a sequence recurs; when most reads are
  unique, the loop settles slowly. MetaFlux learns from whole samples in a random order
  fixed by `amplicon.seed`, and cannot be steered to the better-duplicated ones, so with
  clean reads such a model is a candidate for acceptance (below).

Measured on a genuinely mixed COI library (two libraries of the same insect mock in
opposite orientations, 1.6 million pairs as one sample) trimmed the old way: of the R1
reads DADA2 learned from, half were reverse-complemented, with their quality running
backwards, and 7 % still began with 17 bases of the primer; the R1 model had not
converged after 171 rounds and 10.5 hours when the run was stopped. The same pairs
trimmed with v2.6.0 converged in 6 and 5 rounds. `pick_trunclen` had shown the problem
first: with the poor last bases of the reversed reads now at the start of the read, the
lower quartile of the quality is below the Q20 threshold at the very first cycle (18 in
R1; 16 and 18 in the first two cycles of R2) and at or above it afterwards, and
`pick_trunclen` cuts before the first cycle that fails, so the quality-based cut lands at
0 (`Primary cuts at Q1>=20.0: R1=0, R2=0` and `WARNING: very short primary cuts` in
`logs/pick_trunclen.log`) and the truncation lengths come from `raise_trunc` alone
(270 / 160 here).

**What a stop leaves on disk.** The two models (`05.dada2/error_model_R{1,2}.rds`), the
record of the learning (`stats/dada2/error_model.json`) and the report on the models (the
figures and tables `stats/dada2/error_model_*`, and DADA2's own plots
`stats/dada2/filtered_read_R{1,2}_error_plot.png`) are kept: they are written before the
check. Started with `--keep-going`, MetaFlux also builds the MultiQC report, whose last
module, "DADA2 error model", shows the same report. The stop deletes nothing else either:
if an earlier run had finished in this folder, its tables in `05.dada2/`, `06.taxonomy/`
and `stats/read_tracking.txt` are still there, made with the earlier model, until the run
goes on.

**Fix.** Check the reads first. In `03.stripped/{sample}_R1.fastq.gz` every read should
start with the first bases after the forward primer. One line shows the most common
starts of a sample's reads, to compare with the primer sequences (with your output folder
in place of `/data/run1` and a sample name in place of `<sample>`):

```bash
zcat /data/run1/03.stripped/<sample>_R1.fastq.gz | head -4000 | awk 'NR%4==2' | cut -c1-20 | sort | uniq -c | sort -rn | head
```

A start that matches the primer itself, or the last bases of the other primer, was not
trimmed cleanly. In `stats/falco/<sample>_R1_stripped/fastqc_report.html` (and `_R2_`) the
per-base quality should fall along the read; a rise from low values at the start is the
sign of reversed reads, and flat lines at 2 to 4 values mean binned quality scores. A
`very short primary cuts` warning in `logs/pick_trunclen.log` can mean the same. With
`orientation: mixed`, `stats/orientation/{sample}.tsv` says how many pairs each pass kept;
a folder trimmed by MetaFlux v2.5.0 or earlier has no such file, and a run with v2.6.0 or
later trims it again, with everything after it (see
[the mixed-orientation second pass](amplicon/overview.md#the-mixed-orientation-second-pass)).
If the reads were trimmed before MetaFlux saw them, or the two files of a sample are
swapped, correct the input: a library in which every pair is swapped is run with its two
files swapped and `orientation: fixed`, a library with pairs in both orientations with
`mixed`. Corrected reads give a new model, which is learned again.

If you lowered `max_consist`, for example to the 10 rounds of QIIME 2, raise it back
first (or delete the line): some libraries need 20 to 100 rounds, and their models are
fine; a model that converges within the new limit needs no acceptance. More rounds rarely
help once 200 were not enough. A run config copied from a template up to v2.6.1 says
`max_consist: 999`, and MetaFlux prints a note at start-up when it sees it: with 999 a
model that keeps changing runs for days before the check, unless the flat-rounds rule
ends it first. Delete the line to use 200.

**Before you accept a model.** Rule out the reads first: a model fitted to a mixture
describes the mixture, however good it looks, and accepting it carries the problem into
every table. Then make the two checks that DADA2's FAQ gives for a model that stopped
before converging:

- **The fit figures**, `stats/dada2/error_model_fit_R1.png` and `error_model_fit_R2.png`.
  Each panel is one kind of misread: `A2C` is a true A read as C. Along the bottom is the
  consensus quality score (the mean quality of identical reads), up the side the rate, on
  a log scale. The grey points are the rates DADA2 counted in the reads: filled at the
  quality values the reads carry, open at the averaged qualities of identical reads,
  consensus qualities no read carries. The solid line is the model; the dashed lines are the other
  error functions refitted on the same counts; the dotted line is what the quality scores
  promise on paper (Q30: one error in 1,000 bases). A model worth keeping has a solid line
  that follows the filled points in every panel and falls as the quality score rises. Be
  wary of a line that misses most of the filled points, or that rises with the quality
  score over a wide range of scores. A small upturn at the very top is common: on the
  models of the eight MiSeq test sets, 288 of the 342 rises of the loess fit lay among
  the five highest quality scores (34 to 38). `rising_steps` in
  `stats/dada2/error_model_fit.tsv` counts the rises, and the same file gives how far the
  line typically lies from the counted rates. DADA2's
  own pictures of the models, `filtered_read_R1_error_plot.png` and `_R2_`, show the same
  points and line in DADA2's style (red: the promised rate).
- **The convergence**, `stats/dada2/error_model_convergence.png`, with the numbers behind
  it in `stats/dada2/error_model.json`, per read direction: `change_per_round` (how much
  the model moved in each round, summed over all its rates) and `distance_to_final` (how
  far each round's starting table is from the final one, DADA2's `checkConvergence`).
  Close to convergence, the last values are several orders of magnitude below the first
  two.

If both checks pass, DADA2's FAQ says to use the model ("in almost all cases, just
proceed with the learned error rates"). If either fails, look at the reads again rather
than allowing more rounds. How the figures and numbers are made, with the values of the
models MetaFlux was tested on to compare with, is described under
[the report](amplicon/error-model.md#the-report).

**Accepting a model.** The stop message ends with the line that accepts the model,
built from its code:

```yaml
accept_unconverged: ["R2:7ebe3bed"]
```

It goes in the run config of this output folder (the file you give with `--configfile`,
not `config/config.yaml`), inside the `amplicon:` → `dada2:` → `learn_errors:` block it
already has. Do not paste a second `amplicon:` block: YAML keeps only the last of two
blocks with the same name, and your other amplicon settings would be lost (MetaFlux stops
at start-up when it finds one). In a run config that has none of these blocks yet:

```yaml
amplicon:
  dada2:
    learn_errors:
      accept_unconverged: ["R2:7ebe3bed"]
```

Then start MetaFlux again with the same command. The model is not learned again:
`check_error_model` writes its verdict, `stats/dada2/error_model_check.txt` (for example
`R1:0a55be94 converged; R2:7ebe3bed accepted (not converged: flat-rounds rule, round
42)`), and `dada_seqtab` and the steps after it run. Changing the list itself makes no
step run again: the check reads it when it runs. So a code changed or removed in a
finished folder takes effect only when the check runs again (after a new model, after
the report is drawn again, or with `--forcerun check_error_model`): a model whose code is
then not listed stops the run again, and only the verdict file is removed, while with its
code listed the run goes on with `dada_seqtab` and every step after it. Keep the line in
that run config: when a
later run of the folder learns the models again (after a MetaFlux update or from a fresh
clone, for example) from the same reads with the same settings and the same DADA2, it gets
the same model and goes on, while other reads, settings or software, and possibly another
type of computer, give another model with another code, and the run stops again for you
to look at it. The check also runs again whenever the error-model report is drawn again
(after a MetaFlux update that changes the report, for example), and it then needs the
line as well. Never put the line in a config shared by several runs. Together, the verdict file and the run config record
which model the tables were made with.

## Amplicon: the run stops because the reads carry fewer than three quality values

!!! warning "Symptom"

    The run stops at rule `choose_error_model`, after the filtering and before the error
    models are learned, with a message that starts like one of these:

    ```text
    [choose_error_model] STOP: the filtered R1 reads carry fewer than three distinct quality values above Q2 (values seen in the first 100,000 reads of each voting sample: 3, 30). DADA2 cannot learn an error model without quality information, whatever error_model says. ...
    [choose_error_model] STOP: sample S7 (100,000 R1 reads scanned) carries fewer than three distinct quality values above Q2 (3, 30), unlike the other voting samples (32 values, 7-38). ...
    [choose_error_model] STOP: no voting sample carries three or more distinct quality values above Q2 in its R1 reads: S1 (100,000 reads scanned: 3, 30), S2 (100,000 reads scanned: 3, 40). ...
    ```

**Cause.** DADA2 learns its error model from the quality scores of the reads: for each
score, how often a base is misread. A sequencer writes at least three distinct values
above Q2 (instruments with binned quality scores three to seven, a MiSeq some thirty);
values at or below Q2 mark no-calls on Illumina instruments, and the default
`trunc_q: 2` removes them. One or two values mean that the values were rewritten on the
way: ENA serves some runs with Q30 at every base (and Q3 for rejected bases), and SRA
Lite files carry made-up values too. DADA2 cannot learn from such reads: its usual fit
fails with `Error rates could not be estimated`, and the binned fit has at most one
straight line to draw, through two made-up values. One such
sample among good ones would make the one error model of the run wrong for every sample,
so it stops the run as well. MetaFlux counts the values in the first 100,000 filtered
reads of each sample (rule `dada_filter`); only samples with at least 1,000 filtered
reads vote, so a near-empty blank never causes this stop. A raised `trunc_q` removes
values too: `filterAndTrim` cuts each read at its first base at or below that value, and
the message then says so first.

**Fix.** Download the reads from NCBI with sra-tools (`prefetch`, then `fasterq-dump`),
which gives the sequencer's own values (see "Are the quality values real?" in
[checking a public dataset](about/test-datasets.md#check-a-public-dataset-before-adopting-it)),
or leave the sample out of the input folder. If `trunc_q` was raised, set it back to 2.
`stats/dada2/quality_values.tsv` in the output directory lists the values of every sample
and read direction, and the Falco report of each sample's filtered reads
(`stats/falco/<sample>_R1_filtered/fastqc_report.html`, and `_R2_`) shows them along the
read. The setting `error_model` cannot get past this stop: no error function can learn
from reads without quality information. The stop comes before the error models are
learned, so the run makes no new MultiQC report: MultiQC waits for the report on the error
models. In a folder that ran before, whatever lies on disk from the steps after this one
(the error models and their report, the tables of `05.dada2/`, the MultiQC report) comes
from the earlier run. The files named here are on disk.

## Amplicon: the samples of a run carry different quality values

!!! warning "Symptom"

    The run stops at rule `choose_error_model`, before the error models are learned,
    with a message that starts like one of these:

    ```text
    [choose_error_model] STOP: the voting samples do not share one kind of quality values in the R1 reads: 12 samples carry binned values (9, 20, 34) and 3 carry continuous ones (32 values, 7-38); ...
    [choose_error_model] STOP: the voting samples carry two different sets of binned quality values in the R1 reads: most carry 11, 25, 37, but S4, S5 also carry 9. ...
    [choose_error_model] STOP: the filtered R1 reads carry 15 distinct quality values above Q2 (20, 21, ..., 34). MetaFlux reads 3 to 10 values as binned quality scores ...
    [choose_error_model] STOP: the filtered R1 reads carry 10 distinct quality values above Q2 (29, 30, 31, 32, 33, 34, 35, 36, 37, 38), and some of them are neighbours (29 and 30). ...
    ```

**Cause.** DADA2 learns one error model per read direction for all samples of a run, and
`error_model: auto` (the default) picks its error function from the quality values the
reads carry: the loess fit for continuous values (20 or more distinct values above Q2, as
a MiSeq writes them) and the binned fit for binned values (3 to 10, as NovaSeq, NextSeq,
MiSeq i100 and other recent instruments write them). Samples from two instruments, or
from one instrument under two versions of its control software (two binned tables, such
as 11, 25, 37 and 9, 25, 37), cannot share one model fitted to both, so `auto` does not
choose for them. Between 11 and 19 distinct values is seen on no instrument: a mix of
files, values rewritten by a sequencing facility, or MiSeq reads under a raised
`trunc_q` (only the values above it remain) are the likely causes. Neighbouring values
among 3 to 10, such as 29 and 30, are seen on no instrument either: every instrument
table leaves gaps between its values, while MiSeq reads under a `trunc_q` of 28 to 30 keep
only the 8 to 10 values above it, all neighbours, and reads trimmed at a quality
threshold before they reached MetaFlux look the same; the message names `trunc_q` when it
is above 2. Only samples with at least 1,000 filtered reads vote; a smaller sample with
other values gets a note and is ignored.

**Fix.** Run the groups in separate output folders, one input folder per instrument.
If `trunc_q` was raised, set it back to 2, the default, so that the choice sees every
value the reads carry. If you know what the data are and accept one model for all
samples, set the error function in the run config:

```yaml
error_model: loess
```

The line goes under `amplicon:` → `dada2:` → `learn_errors:`, inside the `amplicon:`
block the run config already has; a second `amplicon:` block would replace the first,
and the other amplicon settings would be lost (MetaFlux stops at start-up when it finds
one). `loess` serves any mix with one smooth curve (it
is the fit MetaFlux used for every run up to v2.6.1); the stop then becomes a warning.
`binned` is not advised for two binned tables: `learnErrors` learns from a few whole
samples (until `nbases` bases are in), and the values of a group it does not draw would
be fitted on the averaged qualities of the other group. The per-sample values are in the
message and in `stats/dada2/quality_values.tsv`. As with every stop before the error
models are learned, the run makes no new MultiQC report (in a folder that ran before,
the one on disk, like the other results after this step, comes from the earlier run).

## Amplicon: `error_model: binned` stops on continuous quality values

!!! warning "Symptom"

    The run stops at rule `choose_error_model` with:

    ```text
    [choose_error_model] STOP: amplicon.dada2.learn_errors.error_model is binned, but the filtered R1 reads carry continuous quality values (32 distinct values above Q2, 7-38), as a MiSeq writes them. The binned fit would join the error rates of 32 values with straight lines instead of fitting a curve. Set error_model to auto or loess, or delete the line from your run config.
    ```

    or, when the samples of the run carry different kinds of quality values, with a text
    that starts `... error_model is binned, but together the voting samples' filtered R1
    reads carry 27 distinct quality values above Q2 (9-39), which MetaFlux reads as
    continuous, although not every sample carries that many on its own`.

**Cause.** The binned fit draws straight lines between the error rates counted at the
quality values the reads carry. On binned data those are three to six values, where a
curve has too few points to rest on. On continuous data it would join some thirty
counted rates one by one and follow their noise: on the models of the eight MiSeq test
sets, its rates rose with the quality score in a quarter to two fifths of the steps
between neighbouring scores (92 to 142 of 348 to 372), against at most one in nine for
the loess fit (2 to 40).
A run config copied from a run on binned data is the usual cause. The second text usually
comes from binned and continuous samples in one run: a sample with continuous values
carries 20 or more values on its own, so together the samples count as continuous.

**Fix.** Set `error_model: auto` or `error_model: loess` in the run config, or delete
the line to get the default. For samples of two instruments, run the groups in separate
output folders, or set `error_model: loess` (see
[the samples of a run carry different quality values](#amplicon-the-samples-of-a-run-carry-different-quality-values)).

## Amplicon: the error plot shows straight segments between a few quality values

!!! note "What it looks like"

    Not an error. On a run with binned quality scores, the line of the model in the fit
    figures (`stats/dada2/error_model_fit_R1.png` and `_R2_`, solid) and in DADA2's error
    plots (`stats/dada2/filtered_read_R1_error_plot.png` and `_R2_`, black) is made of a
    few pieces that meet at the quality values the reads carry, bent on the log axis,
    with points at many more quality scores than the reads carry. The title names the fit
    and its values, for example
    `R2 error model: binned fit on 9, 20, 34 - converged (round 6 repeated its starting table)`.

**Cause.** This is how the binned fit works. It joins the error rates counted at the
values the reads carry (9, 20 and 34 in the example) with straight lines; on the log
scale of the plot, a straight line between two rates is drawn bowed upward. The points
between those values are consensus qualities: DADA2 counts the errors of each unique
sequence at the mean quality of its identical reads, so averaged scores such as 27
appear although no read carries them. Between the two upper values they come from
sequences seen in many reads, mostly read without error, and most of their bases lie
below the line. How many of the misreads they hold depends on the data: on one NextSeq
2000 run, with two samples learned together for R1 and one for R2, they held 8 % of the
R1 and 32 % of the R2 bases but only 1.2 % and 2.3 % of the counted misreads; one sample
of the same run, learned alone, gave 7.6 % and 35 % of the bases and 15 % and 13 % of
the misreads. Between the two lowest values they rest on a tenth of a percent of the
bases or less on the runs checked, yet they can hold a good part of the misreads (6 %
and 15 % in that single sample) and lie far above the line. The share of your run is
`share_substitutions_at_averaged_columns` in `stats/dada2/error_model_report.json`, and
the notes in `stats/dada2/error_model_fit.tsv` give it too.

**What would be a problem.** A line that misses the points at the quality values the
reads carry (the filled points of the fit figures), or error rates that rise as the
quality score rises. A line above the nominal line (dotted in the fit figures, red in
DADA2's plots) is normal for every fit: the models MetaFlux was tested on sat 1.2 to 25
times above it at their most common quality value (1.2 to 2.0 on the Zymo D6331 mock and
on one NextSeq 2000 sample learned alone, with either fit; 7 to 25 on the others), on
healthy and faulty libraries alike. The fit figures also show the loess fit refitted on the same counts
(dashed). On the NextSeq 2000 run checked, its curve dipped by up to four orders of
magnitude between the two upper values, where only averaged qualities hold bases, and
rose again towards the top value (a dip below the panel is cut off at its edge); the
binned fit draws a straight line there. The page
[DADA2 error model](amplicon/error-model.md#the-report) explains the figures and gives
the numbers of the models MetaFlux was tested on.

## Amplicon: the run refuses `custom_error_function` at start-up

!!! warning "Symptom"

    The run stops before any step, with a message about the experimental setting
    `amplicon.dada2.learn_errors.custom_error_function` (your own error function), such
    as:

    ```text
    [MetaFlux] amplicon.dada2.learn_errors.custom_error_function.file: /home/me/MetaFlux/functions/my_errfun.R is not an existing file (written in the config as 'functions/my_errfun.R'; a relative path is read from the folder snakemake is started in, here /home/me/MetaFlux)
    [MetaFlux] amplicon.dada2.learn_errors.custom_error_function.name must be the name of an R function as you would type it in R, without quotes: ...; got 'my-errfun'
    ```

**Cause and fix.** MetaFlux checks the block before any step, with a message that names
the setting: a missing file would otherwise stop Snakemake at once with a bare "Missing
input files" for rule `learn_errors` (or `error_model_report`), and a bad name would
surface only after hours of filtering, as an R error in `learn_errors`. Each message
names the key:

| The message says | Cause | Fix |
|---|---|---|
| `error_model is custom, so ... custom_error_function.file must give the R file ...; it is empty` | `error_model: custom` without a file | Set `file`, or set `error_model` back to `auto`, `loess` or `binned`. |
| `error_model is custom, so ... custom_error_function.name must give the name of your error function ...` | `error_model: custom` without a name | Set `name` to the function's name in the file. |
| `custom_error_function sets a file but no name` (or `a name but no file`) | Half of the block, under another `error_model` | Set both, and the report draws the function next to the others for comparison; or leave both empty. |
| `custom_error_function must be a block with two keys` | A path written directly after `custom_error_function:` | Write `file:` and `name:` on the two lines below it, indented. |
| `custom_error_function.file must be text` (or `.name`) | A number or a list | Write the path, or the name, as text. |
| `custom_error_function.file: ... is not an existing file` | A typo, a folder, or a relative path read from another folder than intended | Check the path the message prints. A relative path is read from the folder `snakemake` is started in, not from the folder of the config file; an absolute path avoids the question. |
| `custom_error_function.file: ... cannot be read` | The file's permissions | Make the file readable for the user who runs Snakemake (`chmod u+r`). |
| `custom_error_function.name must be the name of an R function` | A dash or a space in the name, an accented letter, a first character that is a digit or an underscore, or a word R reserves (`function`, `TRUE`, `NA` ...) | Write the name exactly as the file defines the function, for example `my_errfun`. |

Whether the file defines a function of that name, and what the function returns, can only
be checked in R, when `learn_errors` runs (next entry). See
[your own error function](amplicon/error-model.md#your-own-error-function-experimental).

## Amplicon: `learn_errors` stops in your own error function

!!! warning "Symptom"

    With `error_model: custom` (experimental), the run stops at rule `learn_errors`,
    with a message on the terminal and at the end of `logs/learn_errors.log` that starts
    like one of these:

    ```text
    [learn_errors] STOP: MetaFlux could not load your error function (amplicon.dada2.learn_errors.custom_error_function, experimental): /data/error_functions/my_errfun.R does not define my_errfn (it defines: my_errfun). Check the name in the config against the file.
    [learn_errors] STOP: the R1 error function failed in round 0: the custom error function my_errfun (/data/error_functions/my_errfun.R) returned a rate of exactly 0 in 41 cell(s), the first A2G at quality 0. ...
    ```

**Cause.** MetaFlux reads the file, fetches the function by its name and checks every table
the function returns before DADA2 uses it (see
[what MetaFlux checks](amplicon/error-model.md#what-metaflux-checks)). Plain DADA2 stops
on some of these problems with a short message of its own:
`Error matrix must have 16 rows (A2A, A2C, ...).` for a wrong number of rows,
`All error matrix entries must be <=1.` for a value above 1,
`missing value where TRUE/FALSE needed` for a missing value, and
`Error matrix is NULL.` for a data frame or, after
`Error rates could not be estimated`, for an R error inside the function. It only warns
about a rate of 0, and goes on without a word with other columns, other names or rows
that do not sum to 1. MetaFlux checks every table and stops with the reason. "Round 0"
is DADA2's start-up pass, the first call of the function.

**Fix.** The reason follows the first colon:

| The reason says | Cause | Fix |
|---|---|---|
| `reading ..., the file of the custom error function, failed with this R error: ...` | The file cannot be read as R code: a syntax error, an error in a line outside the function, or `library()` of a package the environment of your function (`dada2_custom.yaml`) does not provide | Read R's message; `source()` the file in an R session to see the error. The packages the environment provides are listed under [R packages](amplicon/error-model.md#r-packages). |
| `... does not define ... (it defines: ...)` | The name in the config differs from the function's name in the file | Correct `name`; the message lists what the file defines. |
| `... is not a function but an object of class ...` | The name names data, not a function | Set `name` to the function. |
| `the custom error function ... stopped with this R error: ...` | An R error inside the function | Read R's message. The counts of any saved model have the shape the function receives: `my_errfun(readRDS("05.dada2/error_model_R1.rds")$trans)` tries the function in R on counts of the same shape (those of a model saved by an earlier run; the counts of the failing call are not saved). |
| `the file of the custom error function, ..., cannot be read` | The file was moved, deleted or made unreadable after start-up | Check the file at the path the message prints. |
| `returned an object of class ..., not a numeric matrix` | A data frame or a list | Return a numeric matrix (`as.matrix()`). |
| `returned a matrix of type ..., not a numeric matrix` | A matrix of text or of TRUE/FALSE values | Return the rates as numbers. |
| `returned a matrix of ... rows and ... columns; it must have the shape of the counts` | Columns dropped or added | Return one row per row and one column per column of `trans`. |
| `returned a table whose rows are named ...` (or `columns`) | Names other than those of `trans`, or in another order | Use `rownames(trans)` and `colnames(trans)`, or no names: MetaFlux then copies them. |
| `returned missing or infinite values` | Often 0/0 at quality scores without bases, or the logarithm of 0 | Fill the scores without bases from their neighbours, as DADA2's functions do. |
| `returned values below 0 or above 1` | A rate that is not a probability | Keep every rate between 0 and 1. |
| `returned a rate of exactly 0` | Counted rates used as they are, where no misread was counted, or a number too small for the computer | Keep every misread rate at or above 1e-7, as DADA2's own functions do; a rate of 0 can turn sequences with that misread into ASVs of their own. A rate far below the counted one acts almost the same and is not stopped, only warned about when it is below 1e-7. |
| `returned rows ... that do not sum to 1` | The read-correctly rows (`A2A` and so on) not computed from the misread rates | Set each to 1 minus the three misread rates of its base, at every quality score. |

As with every stop inside `learn_errors`, the run makes no new MultiQC report; the Falco
reports of the filtered reads in `stats/falco/` are on disk. In a folder that ran before
(such as a copy made to try a new version of the function), only the models and the
record are removed: the error-model report, the verdict of the check
(`stats/dada2/error_model_check.txt`), the tables of `05.dada2/` and every result after
them, and the MultiQC report on disk still come from the earlier models until the
learning succeeds again. A warning of the function does not stop the run: it is written
to `logs/learn_errors.log` and to the record.

## Amplicon: the environment of your own error function cannot be built

!!! warning "Symptom"

    On an Apple Silicon Mac or on Linux on ARM, a run that names an error function of
    your own (`custom_error_function`, experimental) stops before any job starts, while
    Snakemake creates its conda environments:

    ```text
    Could not create conda environment from .../workflow/envs/dada2_custom.yaml:
    ...
    PackagesNotFoundInChannelsError: The following packages are not available from current channels:

      - r-scam
    ```

**Cause.** The steps that call your function run in `workflow/envs/dada2_custom.yaml`: the
DADA2 environment plus the R packages mgcv and scam (see
[R packages](amplicon/error-model.md#r-packages)). conda-forge has no r-scam build for
these machines (`osx-arm64`, `linux-aarch64`). Older conda versions call the error
`PackagesNotFoundError`. A run that names no function of its own never builds this
environment.

**Fix.** If your function does not use scam, delete the line `- r-scam` from
`workflow/envs/dada2_custom.yaml` in the MetaFlux folder you run from, and start the run
again: mgcv exists for these machines. A function that needs scam cannot run there until
conda-forge builds it for them; run it on Linux (x86-64) or an Intel Mac, where the file
can be built as shipped.

## Amplicon: an outdated marker pack in `config/markers/` stops the run

!!! warning "Symptom"

    Every amplicon run stops before any step starts, whichever marker it uses, with:

    ```text
    [MetaFlux] config/markers/18S.yaml was written for an older MetaFlux: it has no reference_values. A pack in config/markers/ replaces the shipped one completely, so copy workflow/markers/18S.yaml again and re-apply your change.
    ```

**Cause.** A pack in `config/markers/` replaces the shipped pack of the same name
completely, so an old copy would bring back everything the shipped pack has changed
since: an 18S pack from MetaFlux v2.4.0 would classify against PR2's own letter order,
which hides the domain of many ASVs. Since v2.5.0 every pack declares
`reference_values`, the rules its reference file is checked against before
classification, so a pack without them was written for an older MetaFlux. Every pack in
`config/markers/` is checked on every amplicon run, which is why an outdated 18S pack
stops a 16S run too. A pack for a marker of your own without `reference_values` gets a
message of its own, which says what to add.

**Fix.** Copy the shipped pack (`workflow/markers/<marker>.yaml`) to `config/markers/`
again and re-apply your change, or delete the copy if the shipped pack now does what
you need. For a marker of your own, add `reference_values` with rules for its top
ranks, or `reference_values: {}` to have only the rank shape checked (see
[the reference check](amplicon/markers/index.md#the-reference-check-reference_values)).
An 18S pack that does declare `reference_values` but keeps PR2's original SINTAX
letters (`k, d, p, …`) stops under `method: sintax` with a message about VSEARCH's
fixed letter order; the fix is the same.

## Amplicon: the reference check stops the run

!!! warning "Symptom"

    The run stops at rule `check_taxonomy_reference`, before any classification, with
    a message like this one (shortened):

    ```text
    [check_taxonomy_reference] /shared/dbs/silva_nr99_v144_toGenus_7rank_other-verbatim_trainset.fa.gz does not hold the ranks declared by the 16S marker pack, release 138.2, method rdp (Domain, Phylum, Class, Order, Family, Genus, Species):
      - a lineage deeper than the 6 ranks the rank model expects in this file (the Species rank comes from addSpecies): 473,550 header(s)
          e.g. >Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;
      - the deepest lineage has 7 ranks, but the rank model expects 6 in this file (the Species rank comes from addSpecies)
      The classifier would put names on the wrong ranks, so the run stops before any classification. Either the run config points at a wrong file (references:, and for 16S the SILVA release chosen with references.silva.release), or the reference changed upstream and its marker pack needs a new release block (see the reference_values notes in workflow/markers/16S.yaml).
    ```

**Cause.** A classifier never checks what a rank means: if a reference has one rank more
or less than the marker pack expects, or two ranks swapped, every name lands in the
wrong column and nothing errors. So before any classification MetaFlux checks the exact
file the classifier is about to read against the rank set and the names its pack
declares (see
[the reference check](amplicon/markers/index.md#the-reference-check-reference_values)).
The message names the file, what does not fit, and example headers. Almost always a path
under `references:` points at the wrong file: in the example, a SILVA 144 training set
(seven ranks, with the kingdom) pinned as `silva.train` for a run on SILVA 138.2 (six
ranks). Less often the reference itself has changed upstream.

**Fix.** Point the `references:` entry at the file that belongs to the run's marker and
SILVA release, or remove it, so that MetaFlux uses, and if needed builds, its own (see
[reference files MetaFlux builds itself](reference/configuration.md#reference-files-metaflux-builds-itself)).
For 16S, keep `references.silva.release` and the pinned SILVA files together. If the
reference itself has changed (a new release with another rank), the marker pack needs a
new release block. `logs/refdb/check_taxonomy_reference.log` in the output directory
lists every problem the check found.

## Amplicon: a reference download stops with "md5 mismatch"

!!! warning "Symptom"

    A `fetch_*` rule fails, and its log (`logs/refdb/fetch_<name>.log` in the output
    directory) ends with one of these lines (the second comes from the COI downloads,
    whose MIDORI2 FASTA and bacterial reference files are checked against SHA-256 sums):

    ```text
    md5 mismatch for <file>: expected …, got …
    sha256 mismatch for <file>: expected …, got …
    ```

**Cause.** The downloaded file is not the one the marker pack expects: the download was
cut off, a proxy returned an error page, or the file behind the address has changed.
MetaFlux deletes the file, so it is never used.

**Fix.** Run the same command again; a broken download usually works the second time.
If the same checksum comes back every time, the file on the server has changed. Report
it as an [issue on GitHub](https://github.com/iLivius/MetaFlux/issues) rather than
editing the md5 or SHA-256 in the marker pack.

## Config (YAML) errors

Snakemake reports any config problem with one generic message
(*"Config file is not valid JSON or YAML…"*). That almost always means a YAML
indentation or list-syntax issue, or a stray tab — not a problem with the
parameter values themselves. Validate the file directly to get the exact line:

```bash
python -c "import yaml; yaml.safe_load(open('config/config.yaml'))"
```

A *"expected block end, but found block sequence start"* message points to a
misplaced or over-indented `- ` list item — commonly in `host_genomes`.

A key written twice in the same block is no YAML error: YAML keeps the last of the two
without a word. A second `amplicon:` block pasted at the end of a run config, for
example to add a line a MetaFlux message asks for, would replace the first one whole,
and the marker, the primers and every other amplicon setting of the first would fall
back to the template. MetaFlux therefore reads each file given with `--configfile` once
more at start-up and stops on such a repeat, with the file, the key and both lines:

```text
[MetaFlux] /data/run1.yaml: the key 'amplicon' is written twice in the same block (lines 18 and 29). YAML keeps only the last of the two, without a word, ...
```

Move the settings of the second block into the first, each once, and delete the second.

The full key-by-key reference is in [Configuration](reference/configuration.md).
