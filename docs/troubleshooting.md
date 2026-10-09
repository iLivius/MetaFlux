# Troubleshooting

Four situations come up often enough to be worth writing down: an automatic
truncation length that throws away every read, a PhiX-removal step that stops on
mis-paired FASTQs, small differences between two runs on the same data, and
YAML mistakes in the config. Each entry gives the symptom, the reason behind it,
and what to change.

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

!!! example "What it looked like on a real library, before the guard existed"

    The 18S set in [test datasets](about/test-datasets.md) is where this was found.
    Its R1 quality never falls below Q20, so `auto` put the R1 cut on the ceiling and
    returned `truncLen = 288 / 239`. After primer trimming the reads themselves were:

    | | min | p05 | median | max |
    |---|--:|--:|--:|--:|
    | R1 | 78 | 271 | 281 | 289 |
    | R2 | 54 | 282 | 282 | 298 |

    The ceiling came out at 288 because one sample contained R1 reads that long. Half
    the reads stop at 281. Since `filterAndTrim` discards anything shorter than
    `truncLen`, the cut kept 9 of 908,768 pairs and DADA2 finished with 3 ASVs.

    Retention against the R1 cut, measured on this library:

    | truncLen R1 | 288 | 281 | 275 | 270 |
    |---|--:|--:|--:|--:|
    | pairs retained | 0.0% | 70.3% | 92.4% | 100.0% |

    Eighteen base pairs separate "loses everything" from "loses nothing".

    With the 95% read-coverage cap the R1 ceiling comes in at 271 instead of 288, and
    `auto` now returns 271 / 239 on this library with no manual intervention. It carries
    90.9% of its raw reads through to the final table: 828,018 reads, against 828,170
    for the 270 / 239 that previously had to be worked out by hand — a difference of
    0.02%.

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

## Shotgun: BBDuk PhiX removal fails on mis-paired (pre-filtered) reads

!!! warning "Symptom"

    The shotgun run aborts at rule `decontam_phix` with:

    ```text
    java.lang.AssertionError: List size mismatch: 200 vs 195
    ```

    The two FASTQs hold different numbers of records: the pairing **is** broken.

**Cause.** BBDuk reads the R1 and R2 files in lockstep, one record from each at
a time, and requires both files to hold the same number of records in the same
order. Read length has nothing to do with it — ragged lengths, and mates of the
same pair ending up different lengths, are processed normally. `List size
mismatch: N vs M` means R1 and R2 contain different numbers of reads, which
usually happens when reads were filtered before deposition and orphaned mates
were left behind in one file. SRA and BioProject submissions are the common
case. MetaFlux's own fetch can produce it as well: rule `download_sra` calls
`fasterq-dump --split-files`, which writes a record to `_1` only for any spot
holding a single read (`--split-3` is the option that keeps the two files in
step).

The BBMap host-removal step (rule `decontam_host`) fails on the same input, for
the same reason: BBMap reads the pair in lockstep too and makes the same
requirement. BBMap's message names the fix outright — *"There appear to be
different numbers of reads in the paired input files. The pairing may have been
corrupted by an upstream process. It may be fixable by running repair.sh."*

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

More on host references and masking is in
[Decontamination](shotgun/decontamination.md).

## Shotgun: empty or nearly empty OTU table with `tax_lev` other than `S`

!!! warning "Symptom"

    `shotgun.bracken.tax_lev` is set to `G`, or to any rank other than `S`. The run
    finishes without an error, but `03.abundance/otu_table.tsv` has no rows, or only a
    handful, while Bracken's own `{sample}_output.txt` files list plenty of taxa.

**Cause.** Fixed in v2.5.0. Up to v2.4.0, rule `kraken_biom` ran kraken-biom with its
default settings, which assume a species-level report: it takes the clade count for
species, and for every coarser rank only the reads placed directly on the taxon.
Bracken 3.1 writes 0 in that column for almost every taxon above species (an upstream
bug, [Bracken #315](https://github.com/jenniferlu717/Bracken/issues/315)), so those
taxa dropped out. A genus survived only when Kraken2 had placed no reads on any of its
species, and at class, phylum or domain level nothing survived at all. On the
benchmark data the genus tables kept between 0.06 % and 9 % of the reads: 2 to 114
genera, where Bracken had found 117 to 687. Bracken's `{sample}_output.txt` files were
never affected, and neither were species-level runs.

**Fix.** Update to v2.5.0 and run the same command again, from the same MetaFlux
folder. Snakemake keeps a record of how each output was made (in `.snakemake/`), sees
that rule `kraken_biom` has changed, and rebuilds the table from the Bracken files
already in the output directory: only `kraken_biom` and `finalize_otu_table` run,
Kraken2 and Bracken do not. A dry run (`-n`) shows those two jobs before anything
runs. If they do not appear, for example because the run used
`--rerun-triggers mtime` or the `.snakemake/` folder was deleted, force them:

```bash
snakemake --sdm conda --cores 4 --configfile <your config> -R kraken_biom
```

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
- **`assignTaxonomy` (rdp): reproducible at any thread count, with one exception.**
  DADA2's `assignTaxonomy` breaks exact ties between equally good reference genera at
  random, outside the seed, so a few bootstrap values can differ between identical runs;
  `sintax`, the default, is not affected.

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
