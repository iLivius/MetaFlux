# Rationale

Over many years of amplicon-sequencing and metagenomics work, a set of methods
and parameter choices accumulates into an in-house pipeline. MetaFlux
re-implements that methodology as a more modern, modular Snakemake workflow.

The re-write was not a straight port. It added five things the previous
implementation never had, all aimed at the same target: decisions that used to be
made by eye, or hard-coded for one marker and one sequencing run, are now derived
from the data itself.

## 1. Data-driven read truncation

Choosing where to cut Illumina reads before denoising is the classic by-eye step:
open the quality profile, pick a cycle where the plot "looks bad", hope the pairs
still overlap. MetaFlux can derive it instead.

With `amplicon.trunc_len.mode: auto`, the `pick_trunclen` rule reads the per-base
quality tables Falco produced for the primer-trimmed reads. For each cycle it
takes the lower quartile (Q1, the 25th percentile of quality at that position),
aggregates it across samples as a median, and cuts just before the first
position where that value drops below `q_threshold` (20 by default). The cut
then has to pass an overlap check, because a read pair that no longer overlaps
cannot be merged:

```text
truncLen_R1 + truncLen_R2  ≥  amplicon_length + min_overlap
```

`amplicon.min_overlap` (12 bp) is the design constraint used here; the runtime
threshold DADA2's `mergePairs` actually applies is a separate key,
`amplicon.dada2.merge.min_overlap` (12 bp). When the quality-based cut fails the
overlap check, `resolve_policy` decides what happens:

| Policy | What it does |
|--------|--------------|
| `raise_trunc` *(default)* | Extends the cuts past the quality drop to recover the missing bases, taking each base from whichever read has the better quality at the position it would gain |
| `relax_q` | Lowers the quality threshold one step at a time, down to `q_floor` (15), re-cutting until the overlap holds |
| `error` | Aborts and reports the overlap shortfall in bp, leaving the choice with the analyst |

The chosen values are written to `stats/trunclen.json` and consumed by
`dada_filter`. Setting `mode: manual` bypasses the analysis and uses `manual_r1`
/ `manual_r2` directly. ITS is the exception throughout: the picker still runs
and records its numbers, but because ITS amplicons vary in length `dada_filter`
overrides `truncLen` to `c(0, 0)` and leaves the 3′ tail to `truncQ` and `maxEE`.

See [Amplicon length and truncation](../amplicon/length-and-truncation.md) for
the details, and [Troubleshooting](../troubleshooting.md) for the one case where
`auto` misfires.

## 2. Agnostic amplicon-length inference

The measured amplicon length drives the overlap check on the truncation cuts, the
ASV length filter window, and — for ITS only — DADA2's per-read minimum length.
Hard-coding it means
re-deriving a number by hand for every new primer pair. MetaFlux measures it
instead, from the reference database, using the primers supplied in the config.

The `amplicon_probe` rule works in one of two ways depending on the marker:

| Probe mode | Markers | What it does | Reference used |
|------------|---------|--------------|----------------|
| PCR | 16S, 18S, rpoB | Two-pass in-silico PCR with Cutadapt: keep reference sequences carrying the forward primer at the 5′ end and the reverse complement of the reverse primer at the 3′ end; the surviving trimmed sequences are the predicted products | SILVA, SILVA-Euk, FROGS |
| Direct | ITS, gyrB | Measure lengths straight off a reference that is already amplicon-length — no primer binding sites needed | UNITE UCHIME pre-extracted ITS1/ITS2, DD7RZ8 pre-trimmed amplicons |

Direct mode exists because in-silico PCR is unreliable for these two, for two
different reasons. The common ITS primers (ITS1F, ITS4) bind in flanking rRNA
regions that UNITE sequences do not contain, so a PCR-style probe would recover
almost nothing. gyrB's DD7RZ8 reference has the opposite problem: it is not a
whole-gene database at all — its sequences were already cut down to the
in-silico amplicon (with the Barret 2015 primer sites removed in the process)
when the trainset itself was built. An in-silico PCR against it would be
searching for primer-binding sites that are no longer there, so it would also
recover almost nothing, just for the reverse reason: too little sequence
outside the amplicon for ITS, none at all for gyrB.

Both modes produce the same thing: a length distribution summarised as JSON.
The consumers read it differently. `probe_length_stat` picks which
statistic of that distribution becomes `expected_length` (`min`, `q1`, `median`,
`q3`, `p95`, `p99`, `max`; p95 by default), and that number is used only in the
overlap check on the truncation cuts described above. Unless `length_filter.range`
is set explicitly — which overrides everything below — the ASV length filter
ignores that choice whenever the probe has run and builds its window from the
distribution itself, keeping ASVs inside
`[probe q1 − window_margin, probe p95 + window_margin]`, with
`length_filter.window_margin` at 50 bp. Changing `probe_length_stat` therefore
moves the overlap constraint, not the ASV window.

The probe result is cached under `refdb/cache/`. The filename encodes the marker
type, the reference, and a hash of both primer sequences, so the probe re-runs
only when one of those three actually changes.

!!! note "Automatic is not always right"

    Where the reference cannot support the measurement, `auto` is the wrong
    choice. A probe that recovers nothing at all does stop the run with an error,
    but one that recovers a biased minority looks exactly like a good one, so the
    marker pages call those cases out instead. The 18S V9 primers (1391F/EukBr) are the documented
    case: EukBr sits at the 3′ terminus and most references are truncated before
    it, so only a biased minority of references yield a product. V9 runs need a
    manual `expected_length` — see [18S rRNA](../amplicon/markers/18S.md).

## 3. Orientation-aware primer trimming

Cutadapt strips the 5′ primers and any read-through 3′ adapter before denoising.
No quality trimming happens at this stage, on purpose: DADA2 builds its error
model from per-cycle quality profiles and needs them intact.

Some libraries come back with reads in mixed orientation: part of the R1 file
actually holds reverse reads, and part of R2 holds forward reads. Trimming those
with the forward primer alone simply loses them. Setting
`amplicon.primers.orientation: mixed` adds a second Cutadapt pass with the
primers swapped; its output is reverse-complemented and concatenated with the
first pass, so the swapped reads rejoin the main set instead of being discarded.
With `orientation: fixed` only the first pass runs. The reverse-complemented
primer sequences both passes need are built once by the `revcomp_primers` rule.

The two knobs that matter here are `amplicon.cutadapt.max_error_rate` (0.2,
generous enough for degenerate primer positions) and
`amplicon.cutadapt.min_length` (50 bp, a floor that drops reads left as fragments
after trimming).

## 4. Fine-grained resource control

Every rule's CPU and RAM are set in one `resources` block, keyed by rule name,
with `threads_default` (4) and `mem_mb_default` (2000) covering anything not
listed:

```yaml
resources:
  threads_default: 4
  mem_mb_default: 2000
  threads:
    kraken2: 16
    decontam_host: 16
    assign_taxonomy: 8
  mem_mb:
    kraken2: 20000
    assign_taxonomy: 16000
```

Rules read these through the `threads_for()` and `mem_mb_for()` helpers defined
in `workflow/rules/shared/00_common.smk`, and the values feed straight into
Snakemake's cluster executor. The same workflow therefore runs on a laptop and on
an HPC queue unchanged — only the numbers move. One value is computed rather than
read: Kraken2's memory request is overridden when the workflow is parsed, from
`shotgun.kraken.memory_mapping` and the real size of the database's `hash.k2d`
file, because that is the one requirement that depends on which database was
supplied.

## 5. One pipeline, two modes

A single `mode` key chooses between the amplicon and shotgun paths. They are not
two workflows in one repository: read QC and its MultiQC report, read tracking,
reference handling and the resource settings are shared, and only the processing
core differs.

The practical effect is that a lab runs one installation, one config format and
one command for both kinds of data, and the read-tracking and QC outputs look the
same whichever mode produced them. See
[Choosing a mode](../getting-started/choosing-a-mode.md).

## Reproducibility

The five capabilities above are about getting the analysis right. The
infrastructure around them is about being able to get the same answer again.

- **One conda environment per step.** Every rule that calls an external tool
  points at its own environment file under `workflow/envs/`, and `--sdm conda`
  (software deployment method) has Snakemake build those environments on first
  run and reuse them afterwards. Steps do not share one large environment, so
  the dependencies of one tool cannot force a change in another's. Where a
  specific version matters it is pinned in the file — `r-rcppparallel=5.1.9` in
  the DADA2 environments, for instance, because newer builds break DADA2 at load
  time. A few rules have no environment of their own and run in Snakemake's: those
  that only link files, count reads, fetch or concatenate references, or crunch
  numbers in Python, using standard system utilities whose version does not change
  the result.
- **Resumable.** Snakemake tracks which outputs already exist and are current, so
  an interrupted run continues rather than restarting; `--rerun-incomplete`
  re-does any job that was cut off mid-write.
- **One version-controlled config.** Every parameter lives in a single YAML file.
  Committing that file alongside the results records exactly what was run — the
  mode, the marker, the primers, the reference databases, the thresholds.
- **Fixed random seed.** `amplicon.seed` pins the random-number generator for the
  stochastic steps in the amplicon path. One caveat applies to multithreaded
  sintax; it is written up in [Troubleshooting](../troubleshooting.md).
