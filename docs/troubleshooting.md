# Troubleshooting

Eight situations come up often enough to be worth writing down: an automatic
truncation length that throws away every read, a PhiX-removal step that stops on
the read pairs, a shotgun OTU table that comes out empty at a rank other than
species, a Bracken rank that MetaFlux refuses, a shotgun taxon filter that removes
every taxon, small differences between two runs on the same data, a reference
download that fails its md5 check, and YAML mistakes in the config. Each entry gives
the symptom, the reason behind it, and what to change.

## Amplicon: automatic `truncLen` can drop all reads on very clean runs

!!! warning "Symptom"

    The run gets through primer trimming and QC, then stops at the DADA2
    filtering step (rule `dada_filter`) with:

    ```text
    Error: No reads passed the filter ...
    ```

**Cause.** With `amplicon.trunc_len.mode: auto`, the `pick_trunclen` rule reads
the per-base quality tables that Falco wrote for the primer-trimmed reads
(`stats/falco/{sample}_R{1,2}_stripped/fastqc_data.txt`) and cuts just before the
first position whose per-base lower quartile (Q1, the 25th percentile of quality
at that cycle) falls below `amplicon.trunc_len.q_threshold` — 20 by default — so
the cut keeps the last cycle still at or above it. If quality never
drops that far, and clean MiSeq runs often stay high to the last cycle,
`truncLen` lands on the ceiling instead.

Where that ceiling sits is the whole question. Falco reports a quality bin for every
position up to the **longest** read in a file, however few reads reach it — so a
ceiling taken from the quality table alone describes the longest reads, not typical
ones. On a library where a few reads run to full length and most stop a little short,
a cut placed there discards nearly everything.

MetaFlux guards against this by default: `amplicon.trunc_len.min_read_coverage_pct`
(95) caps each sample at the length that percentage of its reads still reach, before
quality is considered at all. **If you are hitting this error on a default config, the
guard is already on, so the cause is one of the two edge cases below, not a missing
guard** — check them first.

That is the trap. After primer trimming most reads end a little earlier than the
full read length, and DADA2's `filterAndTrim` **discards any read shorter than
`truncLen`** rather than padding it. So a cut placed at the very end of the read
removes almost everything, and DADA2 has nothing left to work with.

!!! example "What it looks like on a real library without the guard"

    The 18S set in [test datasets](about/test-datasets.md) shows it. Its R1 quality
    never falls below Q20, so without the guard `auto` puts the R1 cut on the ceiling
    and returns `truncLen = 288 / 239`. After primer trimming the reads themselves are:

    | | min | p05 | median | max |
    |---|--:|--:|--:|--:|
    | R1 | 78 | 271 | 281 | 289 |
    | R2 | 54 | 282 | 282 | 298 |

    The ceiling comes out at 288 because one sample contains R1 reads that long. Half
    the reads stop at 281. Since `filterAndTrim` discards anything shorter than
    `truncLen`, the cut keeps 9 of 908,768 pairs and DADA2 finishes with 3 ASVs.

    Retention against the R1 cut, measured on this library:

    | truncLen R1 | 288 | 281 | 275 | 270 |
    |---|--:|--:|--:|--:|
    | pairs retained | 0.0% | 70.3% | 92.4% | 100.0% |

    Eighteen base pairs separate "loses everything" from "loses nothing".

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
   in length. `logs/pick_trunclen.log` tells you directly — it prints each sample's
   longest read alongside the length 95% of its reads reach, and emits a `NOTE` when
   the gap is 10 bp or more. A very large gap there means the ceiling is being set
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
    above R2 *did* degrade, so `auto` found a genuine cut for it at 239 and that value
    was kept unchanged; only R1 needed moving. Overriding both to the same round
    number would have thrown away 30 bp of good R2 sequence for nothing.

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

**Fix.** Update to v2.5.0 and run the same command again, from the same MetaFlux
folder. Snakemake keeps a record of how each output was made (in `.snakemake/`) and
sees that rule `kraken_biom`, rule `bracken` and the Bracken environment have changed.
It runs `bracken` again for each sample, on the Kraken2 reports already in the output
directory, then `kraken_biom`, `finalize_otu_table` and MultiQC (and per-taxon read
extraction, when it is on). Kraken2 does not run again. A dry run (`-n`) lists these
jobs before anything runs. If `kraken_biom` is not among them, for example because the
run used `--rerun-triggers mtime` or the `.snakemake/` folder was deleted, force it:

```bash
snakemake --sdm conda --cores 4 --configfile <your config> -R kraken_biom
```

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

**Cause.** Bracken reads a Kraken2 report from top to bottom and counts every row
below the target rank towards the last row of that rank it has read. For a taxon that
has the target rank, that row is its own ancestor, because the report lists every
branch right after its parent. A taxon without the rank is counted towards whichever
row of that rank comes before it. At domain and kingdom level whole groups of taxa
lack the rank:

- **Domain.** Viruses have no domain in NCBI's taxonomy since 2025, and their branch
  comes after all cellular organisms. Their reads, and the share of the root's reads
  that Bracken attributes to them, go to the domain listed last. Kraken2 lists sister
  rows by read count, largest first, so this is the domain with the fewest reads,
  usually Archaea. In a soil sample, Archaea had 2,612 reads at phylum level and 35,636
  at domain level, from only 71 viral reads and the root reads attributed to them. In a
  sample of the ZymoBIOMICS Gut Microbiome Standard (D6331), Archaea went from 6,415
  reads at phylum level to 14,028 at domain level.
- **Kingdom.** Most protists have no kingdom in NCBI's taxonomy, and GTDB has no
  kingdom rank at all. The protists' reads go to the kingdom listed just before them:
  in the same soil sample, Fungi had 7,957 reads in its phyla but 12,001 at kingdom
  level. On a GTDB database Bracken finds no kingdom and stops with "no reads found",
  after Kraken2 has run.

**Fix.** Set `tax_lev` to `S` (the default) or another rank from `P` to `G`. Domain
and kingdom totals can be read from any of these tables: add up the rows of
`otu_table.tsv` by their `d__` or `k__` slot. That sum counts the reads Bracken placed
at the table's rank, not every read of the domain or kingdom: reads that Kraken2 left
above that rank, and that Bracken could not push down, are in no row (see
[Bracken abundance re-estimation](shotgun/overview.md#6-bracken-abundance-re-estimation)).
With `profiler: sylph` Bracken does not run, and `tax_lev` is not checked.

## Shotgun: the taxon filter removed every taxon

!!! warning "Symptom"

    The run stops at rule `finalize_otu_table` with a message like this one (wrapped
    here; the numbers depend on the run):

    ```text
    [finalize_otu_table] shotgun.taxonomy_filter removed every taxon of the table
    (111 taxa; keep: [k__Bacteria], discard: []). The d__ values in this table are:
    d__Archaea (1 taxon), d__Bacteria (107 taxa), d__Eukaryota (3 taxa). A token must
    equal one whole slot of the lineage, prefix included (d__ k__ p__ c__ o__ f__ g__
    s__). Check the tokens, or set shotgun.taxonomy_filter.enabled: false.
    ```

    No `otu_table.tsv` or `otu_table.biom` is written.

**Cause.** `shotgun.taxonomy_filter` compares each token with whole slots of the
lineage, so a token with the wrong prefix or a misspelt name matches nothing. A `keep`
list in which no token matches removes every row, and so does a `discard` list that
covers every row. An empty table would look like a finished run, so MetaFlux stops
instead, as the amplicon filter does. A common cause is a token written for another
table: in MetaFlux's shotgun tables `d__` holds the domain (`d__Bacteria`) and `k__` the
kingdom (`k__Pseudomonadati`), so `k__Bacteria` matches nothing.

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
| DADA2 `learnErrors(randomize = TRUE)` | `dada_seqtab` | Subsamples reads to build the error model |
| DADA2 `assignTaxonomy` (`method: rdp`) | `assign_taxonomy` | Bootstrap confidence |
| VSEARCH `--sintax` (`method: sintax`) | `assign_taxonomy` | Bootstrap confidence, seeded from a random source by default |

Without a fixed seed each of these draws a different sample every run, so a
handful of ASVs or taxonomy calls sitting near a decision threshold can flip
between two runs. This is expected DADA2/VSEARCH behaviour — not a bug, and not
corrupted data.

One source of run-to-run difference used to have nothing to do with a random draw, and
the seed did not touch it: the order of the reads. Bowtie2 removes phiX on several
threads and wrote the surviving pairs in whatever order they finished, so every run
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
  `sintax`, the default, is not affected.

## Amplicon: a reference download stops with "md5 mismatch"

!!! warning "Symptom"

    A `fetch_*` rule fails, and its log (`logs/refdb/fetch_<name>.log` in the output
    directory) ends with:

    ```text
    md5 mismatch for <file>: expected …, got …
    ```

**Cause.** The downloaded file is not the one the marker pack expects: the download was
cut off, a proxy returned an error page, or the file behind the address has changed.
MetaFlux deletes the file, so it is never used.

**Fix.** Run the same command again; a broken download usually works the second time.
If the same checksum comes back every time, the file on the server has changed. Report
it as an [issue on GitHub](https://github.com/iLivius/MetaFlux/issues) rather than
editing the md5 in the marker pack.

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

The full key-by-key reference is in [Configuration](reference/configuration.md).
