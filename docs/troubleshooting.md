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
guard is already on and something else is going on** — check the two situations at the
bottom of this entry first.

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
    `auto` now returns 271 / 239 on this library with no manual intervention — within
    0.02% of the 270 / 239 that had to be worked out by hand. It carries 90.9% of its
    raw reads through to the final table.

**Why it can still happen with the guard on.** Two cases:

1. **`min_read_coverage_pct` was raised.** At `100` the guard is off and the old
   behaviour is back exactly; anything above about 99 is close enough to be suspect.
   This is the first thing to check, since it is the only way to reintroduce the
   failure through configuration alone.
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

`amplicon.seed` (42 by default) fixes the generator for all three, but only two
of them become fully reproducible:

- **`learnErrors` and `assignTaxonomy` (rdp): fully reproducible at any thread
  count.** Verified by running the amplicon test set twice at 8 threads with the
  same seed — `seqs.fasta`, ASV counts and rdp-path taxonomy came back
  byte-identical.
- **`--sintax`: not.** VSEARCH runs its bootstrap confidence step across threads
  that share one random-number stream, so which thread consumes which draw still
  depends on scheduling, seed or no seed. Verified the same way: with
  `amplicon.seed` fixed, 8-threaded `assign_taxonomy` re-runs still shifted a
  handful of per-rank confidence values (a genus call at 0.99 in one run, 1.00 in
  the next), while single-threaded re-runs came back identical.

!!! tip "If byte-identical sintax output matters more than speed"

    Set the taxonomy step to one thread:

    ```yaml
    resources:
      threads:
        assign_taxonomy: 1
    ```

    Otherwise, treat sub-percent confidence drift near `sintax_cutoff` as normal
    run-to-run noise — the same behaviour as before the seed was added.

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
