# DADA2 error model

When a sequence differs from a more abundant one by only a few bases, DADA2 checks whether sequencing errors could account for the reads carrying those differences. A handful of reads may be consistent with errors. If the sequence occurs more often than the error model predicts, DADA2 can infer a separate amplicon sequence variant (ASV). Otherwise, it assigns those reads to the more abundant sequence.

The *error model* describes how often each base is read as another at a given quality score, such as an A being called as G at Q30. These rates can vary between instruments and runs, and between R1 and R2. In the MetaFlux test sets, the estimated error rate at the most common quality score was 1.2 to 25 times the rate implied by that score. See [Goodness of fit](#goodness-of-fit).

DADA2 learns these rates from the run’s own reads before denoising. If the model underestimates errors, DADA2 can call spurious ASVs. If it overestimates them, DADA2 can merge genuine variants with more abundant sequences.

This page explains how MetaFlux fits a separate model for R1 and R2, chooses a fitting function from the observed quality scores and reports the results for each run. It also covers what to do if a model fails to converge and MetaFlux stops before denoising, and how to try a custom error function (experimental).

The terms used throughout:

| Term | Meaning |
|---|---|
| transition | a true base read as a base: `A2G` is a true A read as G, `A2A` a true A read correctly; 16 in all, 12 of them misreads |
| transition counts | how often each transition was seen at each quality score, counted by comparing reads with the ASV they were assigned to |
| error function | the step that turns the transition counts into a table of error rates, one per transition and quality score (`loess`, `binned`, or `scam`) |
| round | one pass of denoising with the current rates and counting; learning repeats rounds until the rates stop changing |
| converged | learning ended because a round gave back a table already seen |
| continuous, binned | quality values written as nearly every score from 2 to 38 (MiSeq), or as a few fixed values such as 9, 20 and 34 (recent instruments) |
| model code | a short code that names one model, such as `R2:7ebe3bed`, used to accept a model that did not converge |

Six rules of [step 6 of the amplicon path](overview.md#6-asv-inference) take part:

| Order | Rule | What it does | Writes |
|---|---|---|---|
| 1 | `dada_filter` | Filters each sample's reads and counts the quality values in the first 100,000 filtered reads of each file. | `04.filtered/`, `stats/dada2/{sample}.filter_stats.json`, `stats/dada2/{sample}.quality_values.json` |
| 2 | `choose_error_model` | Picks the error function for R1 and for R2 from those values. | `stats/dada2/error_model_choice.json`, `stats/dada2/quality_values.tsv` |
| 3 | `learn_errors` | Learns the two error models and records how their learning ended. | `05.dada2/error_model_R{1,2}.rds`, `stats/dada2/error_model.json` |
| 4 | `error_model_report` | Writes the report on the models: tables, figures and a section of the MultiQC report. | `stats/dada2/error_model_*`, `stats/dada2/filtered_read_R{1,2}_error_plot.*`, `stats/dada2/multiqc/` |
| 5 | `check_error_model` | Lets the run go on, or stops it before denoising when a model did not converge and was not accepted. | `stats/dada2/error_model_check.txt` |
| 6 | `dada_seqtab` | Denoises every sample with the two models. | the tables of `05.dada2/` |

Up to v2.6.1, `dada_seqtab` learned the models itself, with DADA2's default error
function, and used them however their learning had ended.

## What `learnErrors` fits

DADA2's `learnErrors` fits one error model per read direction for the whole run, from
the filtered reads of all samples: one for R1 and one for R2, because the second read
of an Illumina run is usually the noisier one. A model is a table with one row per
*transition*, a true base read as A, C, G or T (16 rows; the 12 that are misreads are
named like `A2C`, a true A read as C), and one column per quality score. Each cell holds
the probability of that transition at that quality score.

How `learnErrors` gets there:

1. **It draws reads.** It reads whole samples in a random order until
   `amplicon.dada2.learn_errors.nbases` bases are in (default `1e8`, DADA2's own
   default); the last sample read is taken whole, so the model can rest on more bases
   than that. The order is drawn with R's random numbers from `amplicon.seed`, and R2
   is learned after R1 in the same R session, so its draw continues from R1's. Before
   the draw MetaFlux sorts the files byte by byte, so the language setting (locale) of
   the shell that runs MetaFlux does not change their order. The same files, seed and
   software give the same samples, and the same models.
2. **It starts from the worst case:** every error rate at the highest possible value.
3. **It repeats two steps** on the drawn reads, a loop DADA2 calls *self-consistency*:
   it denoises them with the current error rates, then counts, at each quality score,
   how often each true base was read as each base, by comparing every read with the
   ASV it was assigned to. These are the *transition counts*. The *error function* then
   turns the counts into the error rates of the next round.
4. **It stops** when a round gives back a table that already went into a round, or at
   a limit (see [How learning ends](#how-learning-ends)).

The quality score of a count is a *consensus* quality: DADA2 gives each unique
sequence the mean quality of its identical reads at each position, rounded to a whole
number, and counts the sequence's transitions there. That is why the bottom axis of
every error plot reads "consensus quality score", and why reads with binned quality
scores, which carry only a few values, give counts at many more. On a 16S run made on a
NextSeq 2000, whose reads carry the values 9, 20 and 34, the counts had bases at every
score from 9 to 34. The scores in between, which no read carries, come from sequences
seen in several reads, mostly reads identical to their ASV: they held 8.2 % of the R1
and 32.4 % of the R2 bases, but only 1.4 % and 2.8 % of the counted misreads.

## Continuous and binned quality scores

The kind of quality values the reads carry decides which error function fits the
models (see [Choosing the error function](#choosing-the-error-function)), so MetaFlux
reads it from the reads before learning.

A MiSeq writes nearly every quality value from 2 to 38: filtered MiSeq reads carry 28 to
32 distinct values above Q2 (counted on 56 files of the six test sets). Recent
Illumina instruments write only a few values, called *bins*: 9, 20 and 34 on a
NextSeq 1000/2000 with the 600-cycle kit, for example. The box "Binned quality scores"
under [choosing the cut from the quality profile](length-and-truncation.md#choosing-the-cut-from-the-quality-profile)
lists them per instrument.

The values depend on the instrument, its control software and even its reagent kit, so
MetaFlux never takes them from that table: it reads them from your reads. Rule
`dada_filter` counts the values in the first 100,000 reads of each filtered file and
writes them, with the bases at each, to `stats/dada2/{sample}.quality_values.json`. Every
value a file carries showed up within its first 1,000 reads, and the count takes about
7 seconds per file on 2×300 reads. Rule `choose_error_model` then pools the samples:

- **Only samples with at least 1,000 filtered reads vote.** A blank or a failed sample
  with a few dozen reads can lack values by chance. It is listed, with its values, in
  `stats/dada2/quality_values.tsv`, but it decides nothing. If no sample reaches 1,000
  filtered reads, every sample votes, with a warning.
- **Q0 to Q2 do not count** towards the kind of quality values. They are the no-call
  values of every Illumina table, and the default `trunc_q: 2` cuts each read at its
  first base at or below Q2, so they rarely reach the filtered reads. When they do, they
  stay among the values the binned fit may be anchored on (below).

Three places show the values. The report draws the bases at each value
(`stats/dada2/error_model_quality_values.png`). The Falco report of each sample's
filtered reads (`stats/falco/<sample>_R1_filtered/fastqc_report.html`, and `_R2_`) draws
the quality along the read as box plots, which on binned reads sit flat on the common
values. The MultiQC report shows only the mean quality at each position of the read, and
a mean of a few values can fall anywhere between them, so it hides the bins.

## Choosing the error function

In every round of learning, the error function turns the transition counts into error
rates. The counted rates cannot be used as they are: at quality scores that hold few
bases they rest on a handful of misreads, and at scores without bases there is no rate at
all. The function therefore fits, for each kind of misread, a curve of the error rate
against the quality score, and reads the rates off the curve. MetaFlux offers two:

- **`loess`**, DADA2's default (`loessErrfun`). For each kind of misread it fits a
  smooth curve of the error rate (on a log scale) against the quality score, weighted
  by the bases at each score. It is made for continuous quality values, and it is the
  fit of every MetaFlux run up to v2.6.1.
- **`binned`**, DADA2's `makeBinnedQualErrfun`. It joins the error rates counted at the
  quality values the reads carry with straight lines; DADA2's author recommends it for
  binned quality values. It first appeared in DADA2 1.36.0, and MetaFlux, which runs
  DADA2 1.34.0, ships an unchanged copy of it (see
  [Citation and references](../about/citation.md#license)).

A third value, `custom`, learns both models with an R function of your own; it is
experimental and described under
[Your own error function](#your-own-error-function-experimental).

The setting `amplicon.dada2.learn_errors.error_model` picks between them, for R1 and R2
separately, from the number of distinct values above Q2 that the voting samples carry
together:

| Distinct values above Q2 | Read as | `auto` | `loess` | `binned` | `scam` |
|---|---|---|---|---|---|
| fewer than 3 | no quality information | stop | stop | stop | stop |
| 3 to 10, no two of them neighbours | binned | binned fit | loess fit, with a warning | binned fit | SCAM fit |
| 3 to 10 with neighbours (such as 29 and 30), or 11 to 19 | undetermined | stop | loess fit, with a warning | binned fit, with a warning | SCAM fit |
| 20 or more | continuous | loess fit | loess fit | stop | SCAM fit |

No instrument table has fewer than three values above Q2, or 11 to 19 of them, or two
neighbouring values: binned values always leave gaps (11, 25, 37 on a NovaSeq 6000). One
or two values mean that the values were rewritten on the way; 11 to 19 usually mean a mix
of instruments, values rewritten by a sequencing facility, or a raised `trunc_q`, which
removes the lower values. Neighbours among 3 to 10 values are continuous values with the
lower ones gone: `trunc_q` 28 to 30 leaves MiSeq reads with only the 8 to 10 values above
it, which the binned fit would otherwise take without a word. Reads trimmed at a quality
threshold before they reached MetaFlux look the same. Neighbours that only appear between
samples, each of which carries a binned table of its own (a NovaSeq X writes 40 on some
versions of its control software and 41 on others), are two binned tables, not an
undetermined kind (see the stops below).

DADA2 learns one model per read direction for all samples of a run, so all voting
samples must carry the same kind of values. Under `auto`, MetaFlux does not choose for
a run in which they do not, and says why. The binned fit is anchored on the values the
voting samples carry, Q2 and below included, that have bases in the reads `learnErrors`
drew. A value inside the range at which one of the four nucleotides has no bases is left
out, and `logs/learn_errors.log` says so; with fewer than two anchors left, the learning
stops with the reason.

**The stops.** The messages name the samples and the quality values concerned and, where
it helps, the Falco reports of the samples' filtered reads; the troubleshooting page has
an entry for each kind.

| Message (it starts with `[choose_error_model] STOP:`) | Setting | Usual cause | How to go on |
|---|---|---|---|
| `the filtered R1 reads carry fewer than three distinct quality values above Q2` | any | Values rewritten on the way (reads downloaded from ENA, SRA Lite files), or a raised `trunc_q` ([entry](../troubleshooting.md#amplicon-the-run-stops-because-the-reads-carry-fewer-than-three-quality-values)). | Fetch the reads from NCBI with sra-tools; set `trunc_q` back to 2. |
| `sample S7 (...) carries fewer than three distinct quality values above Q2` | any | One sample with rewritten values among good ones. | Fetch that sample again, or leave it out. |
| `no voting sample carries three or more distinct quality values above Q2` | any | Every sample's values rewritten, with different made-up values (two per sample, three or more together). | Fetch the reads from NCBI with sra-tools. |
| `the filtered R1 reads carry 15 distinct quality values above Q2` | `auto` | A mix of instruments, values rewritten by a facility, or MiSeq reads under a raised `trunc_q` ([entry](../troubleshooting.md#amplicon-the-samples-of-a-run-carry-different-quality-values)). | If you know what the data are, set `error_model: loess` or `binned`. |
| `the filtered R1 reads carry 10 distinct quality values above Q2 (29, 30, ...), and some of them are neighbours` | `auto` | MiSeq reads under a raised `trunc_q` (the message then names it), or reads trimmed at a quality threshold before MetaFlux ([entry](../troubleshooting.md#amplicon-the-samples-of-a-run-carry-different-quality-values)). | Set `trunc_q` back to 2; or, if you know what the data are, set `error_model: loess` or `binned`. |
| `the voting samples do not share one kind of quality values` | `auto` | Binned and continuous samples in one run. | Run the groups in separate output folders, or set `error_model: loess`. |
| `the voting samples carry two different sets of binned quality values` | `auto` | Two instruments, or two versions of the control software. | Separate output folders, or `error_model: loess`; `binned` is not advised here. |
| `error_model is binned, but the filtered R1 reads carry continuous quality values` | `binned` | A run config copied from a run on binned data ([entry](../troubleshooting.md#amplicon-error_model-binned-stops-on-continuous-quality-values)). | Set `error_model: auto` or `loess`, or delete the line. |
| `error_model is binned, but together the voting samples' filtered R1 reads carry ...` | `binned` | Binned and continuous samples in one run: together their values count as continuous. | Run the groups in separate output folders, or set `error_model: loess`. |

Set in the run config, `error_model: loess` turns the three stops that only `auto` makes
into warnings, and the loess fit is then used for every sample, as up to v2.6.1.
`error_model: binned` does the same for undetermined values and for two sets of binned values,
but it is not advised for the second: `learnErrors` learns from a few whole samples, and
the values of a group it does not draw would be fitted on the averaged qualities of the
other group. On binned and continuous samples together, `binned` stops as it does on
continuous values.

**The warnings and notes**, in `logs/choose_error_model.log` and on the terminal:

- `loess` set on binned values: DADA2's author recommends the binned fit for them; the
  run goes on with the loess fit, as up to v2.6.1.
- A sample with fewer than 1,000 filtered reads carries values the voting samples do
  not: its extra values are neither used for the choice nor as anchors. If
  `learnErrors` draws that sample under the binned fit, a value outside the anchors
  stops the learning, and `logs/learn_errors.log` gives the reason.
- No sample has 1,000 filtered reads: every sample votes.
- R1 and R2 carry different kinds of values, which is unusual for one sequencing run:
  check the input files.

Where a message asks you to set the error function, add the line under `amplicon:` →
`dada2:` → `learn_errors:` inside the `amplicon:` block your run config already has. A
second `amplicon:` block would replace the first, because YAML keeps only the last of two
blocks with the same name, and your other amplicon settings would be lost; MetaFlux stops
at start-up when a run config writes a key twice in the same block, and names both lines.

**The default** is `auto`. Under `auto`, binned data get the binned fit, so their ASV
tables differ from those of v2.6.1 on purpose; `error_model: loess` keeps the fit of
v2.6.1. Continuous data get the loess fit, as before. The choice follows the kind of
quality values, not how closely each function fits the counts; the reason is under
[Evidence and limits](#evidence-and-limits).

## How learning ends

The loop ends in one of three ways, which `stats/dada2/error_model.json` records as
`stopped_by`. In each round MetaFlux measures the *change*: how far the table that came
out of the round is from the table that went in, summed over all its cells.

1. **DADA2's own stop (`dada2_repeat`).** A round gives back a table that already went
   into a round. MetaFlux calls such a model *converged*. Usually it is the table that
   went into the same round: both models of the 16S test set stopped that way, in round
   five. On slow libraries a few reads can keep switching between two ASVs, so the loop
   goes round in a circle until a table comes back from several rounds earlier. The
   slow but healthy runs seen ended in such cycles, of 11, 26 and 49 rounds: the R2
   model of the BF3/BR2 library of the COI insect mock after 32 rounds, and the R1 and
   R2 models of a NextSeq 2000 16S run after 38 and 101 rounds.
2. **The flat-rounds rule (`flat_rounds`),** `max_flat_rounds` (default `20`; `0`
   switches it off). MetaFlux ends the learning when the model has stopped improving:
   the change in a round is less than 1/1000 of the change in round 1, and the lowest
   change so far fell by less than 10 % over the last 20 rounds. MetaFlux then hands
   DADA2 back the table that went into that round. DADA2 takes it for a repeat, reports
   convergence after that round in its usual line, and keeps that table as the model;
   MetaFlux's own line below it says that the model did not converge.
3. **The cap (`cap`),** `max_consist` (default `200`). DADA2 prints `Self-consistency
   loop terminated before convergence.` and keeps the last table. DADA2 prints that line
   whenever the last allowed round is reached, even when the table repeated in that very
   round; MetaFlux's line then says that it did converge.

A model ended by the flat-rounds rule or the cap counts as *not converged*, and the run
stops before denoising until its code is accepted
([next section](#the-check-and-acceptance-by-code)).
DADA2 on its own does not stop: at its limit it carries on with the last table, as
QIIME 2 and nf-core/ampliseq do.

![The convergence figure of the COI insect mock library SRR8082166, two panels on a log scale. Left, R1: the change in the round falls from about 46 in round 1 to about 0.0002 by round 5, and round 6 gives back the table of round 5. Right, R2: the change falls below the dashed line at a thousandth of round 1 by round 3, then stays between about 0.0001 and 0.015 for thirty rounds, until round 32 gives back the table of round 22, marked by a grey band.](../assets/error_model_convergence_coi_mock.png)

*The convergence figure of the BF3/BR2 library of the COI insect mock (SRR8082166, one
sample of 653,378 filtered read pairs, MiSeq 2×300). Black: the change in each round;
grey, dashed: how far each round's starting table is from the final model (DADA2's
`checkConvergence`). R1 ended in round 6 by giving back the table of round 5 (a cycle of
2 rounds), so the distance of round 5 is 0: its line breaks there, and a triangle at the
bottom marks the 0. R2 changed by less than a thousandth of round 1's change from round 3 on, and
took until round 32 to give back a table seen before, that of round 22 (grey band: a
cycle of 11 rounds). The flat-rounds rule did not act: the lowest change kept reaching
new lows more than 10 % lower (in rounds 13, 14 and 21), so no stretch of 20 rounds
passed without one. Both models converged, and the run went on.*

**`max_consist` is an upper limit, not a speed setting.** A model that converges stops at
the same round whatever the limit. Of 194 learning runs on real data, 186 converged within
8 rounds and the slowest took 101. A round took 2 to 5 seconds on the 16S test set and 38
minutes on two samples of 1.6 million NextSeq 2000 reads learned together, so on a large
library the cap alone can take days to reach; the flat-rounds rule is there for that
case. Up to v2.6.1 the template said `999`; MetaFlux prints a note at start-up when a run
config still says so. A lower limit, such as the 10 rounds of DADA2's and QIIME 2's
defaults, saves time only on the libraries that need more rounds, and those are the runs
it stops.

**How the flat-rounds rule was set.** The rule was calibrated by replaying four learning
runs round by round, all with the loess fit: R1 (two samples, 38 rounds) and R2 (one
sample, 101 rounds) of one binned NextSeq 2000 16S run, a MiSeq COI library (32 rounds),
and a COI library with reads of both strands in one file, which never converged. With
`max_flat_rounds: 20` the rule stopped the 101-round run at round 42, where denoising that
sample with the stopped model changed 0.02 % of its reads, and the mixed library at round
32; it let the other two end by DADA2's own repeat. Over every value tested, 5 to 50, the
largest change was 0.04 % of the reads.

The same replays were then run with the binned fit, on six learning runs: R1 and R2 of the
same 16S run, on its full filtered reads and on a subsample of 1e8 bases, and R1 and R2 of
the Zymo D6331 standard (NextSeq 2000). Five converged by DADA2's own repeat in 4 to 18
rounds. R1 of the full 16S reads needed 95 rounds (85 h); the rule would have stopped it at
round 49, 34 h earlier, where denoising with the stopped model changed 0.007 % of the
reads. Over every value tested, 5 to 50, the largest change was 0.02 %. No replay reached
the cap of 200 rounds.

**Why a model ended early counts as not converged.** The numbers cannot tell a slow but
healthy library from one with faulty reads. The COI library with reads of both strands in
one file changed by about 0.003 per round from round 4 on, as a slow healthy library
does, and the goodness of fit of its model lay within the range of the healthy ones
([reference values](#goodness-of-fit)). Only the reads showed the fault. So the rule saves
the time of the loop, and the check after it asks for a look at the reads.

**Following a long learning.** `tail -f logs/learn_errors.log` in the output folder shows
one MetaFlux line per round, after DADA2's `selfConsist step` line, with the change, the
lowest change so far and the seconds the round took. On the 16S test set (six samples,
MiSeq 2×300), the R1 part of the log reads:

```text
[learn_errors] R1: loess fit (loessErrfun, DADA2's default), chosen by auto: 32 distinct quality values above Q2 (7-38)
[learn_errors] Learning error rates (R1, nbases=1e+08)...
52975572 total bases in 217113 reads from 6 samples will be used for learning the error rates.
Initializing error rates to maximum possible estimate.
selfConsist step 1 ......
[learn_errors] R1 round 1: change 27.8; 3 s
   selfConsist step 2
[learn_errors] R1 round 2: change 0.543 (round 1: 27.8; lowest so far 0.543); 5 s
   selfConsist step 3
[learn_errors] R1 round 3: change 0.0404 (round 1: 27.8; lowest so far 0.0404); 5 s
   selfConsist step 4
[learn_errors] R1 round 4: change 0.0141 (round 1: 27.8; lowest so far 0.0141); 5 s
   selfConsist step 5
[learn_errors] R1 round 5: change 0 (round 1: 27.8; lowest so far 0); 5 s
Convergence after  5  rounds.
```

and, after R2, how each learning ended:

```text
[learn_errors] R1: converged: round 5 gave back the table that went into it.
[learn_errors] R1: change per round: 27.8, 0.543, 0.0404, 0.0141, 0
[learn_errors] R1: distance of each round's starting table to the final one: 28, 0.554, 0.0363, 0.0141, 0
[learn_errors] R1: model code R1:d73e63a3
```

When the error function itself fails, the job stops with a line that starts with
`[learn_errors] STOP:` and gives the reason, which DADA2 would hide behind `Error rates
could not be estimated`. On a cluster this line goes to the scheduler's log of the job;
it is always at the end of `logs/learn_errors.log`.

## The check and acceptance by code

A model that did not converge is not used without a look at the reads: the run stops
before denoising, and goes on once the model is accepted by its code.

When both models converged, `check_error_model` writes a one-line verdict,
`stats/dada2/error_model_check.txt`, and the run goes on. When a model was ended by the
flat-rounds rule or the cap, the run stops here, before any denoising, unless the model's
code is listed in `amplicon.dada2.learn_errors.accept_unconverged`. The check runs on the
machine Snakemake runs on, also under a cluster profile, so its message reaches your
terminal; the same text is in `logs/check_error_model.log`. It comes in this order:

1. which model did not converge, with its code;
2. how its learning ended, with the change in rounds 1 and 2 and in the last round;
3. what to check in the reads first: the warnings of `pick_trunclen`, the Falco report of
   the primer-trimmed reads, the start of the reads (for primer bases left on them), and
   the quality values;
4. which figures and tables of the report to look at;
5. the line that accepts the model, and where it goes.

The troubleshooting entry
[DADA2's error model never converges](../troubleshooting.md#amplicon-dada2s-error-model-never-converges)
walks through these checks.

**Model codes.** A code such as `R2:7ebe3bed` is the read direction and the first 8
characters of an md5 sum of the model's numbers: every rate of the final table, written
with 17 significant digits. It names one model. The same reads, settings and software
give the same model, and so the same code.

**Accepting a model.** The stop message ends with the line to add, built from the codes
of every model of the run that did not converge:

```yaml
accept_unconverged: ["R2:7ebe3bed"]
```

It goes in the run config of this output folder (the file you give with `--configfile`,
not `config/config.yaml`), under `amplicon:` → `dada2:` → `learn_errors:`, inside the
`amplicon:` block the run config already has. Do not paste a second `amplicon:` block:
YAML keeps only the last of two blocks with the same name, and your other amplicon
settings would be lost (MetaFlux stops at start-up when it finds one). In a run config
that has none of these blocks yet:

```yaml
amplicon:
  dada2:
    learn_errors:
      accept_unconverged: ["R2:7ebe3bed"]
```

Then start MetaFlux again with the same command. The model is not learned again: the
check writes its verdict, for example
`R1:0a55be94 converged; R2:7ebe3bed accepted (not converged: flat-rounds rule, round 42)`,
with a warning on the terminal, and `dada_seqtab` and the steps after it run. The check
reads the list when it runs, not as a parameter of a step, so adding or removing a code
makes no step run again. A code changed or removed in a finished folder therefore takes
effect only when the check runs again: after a new model, after the report is drawn
again, or when you ask for it with `--forcerun check_error_model`. A model whose code is
then not listed stops the run again, and only the verdict file is removed; with its code
listed, the run goes on, and `dada_seqtab` and every step after it run again. Keep the
line in that run config: when a later run of the folder
learns the models again (after a MetaFlux update, for example) from the same reads with
the same settings and the same DADA2, it gets the same model and goes on, while other
reads, settings or software, and possibly another type of computer (below), give another
model with another code, and the run stops again for a look. The check also runs again,
with `dada_seqtab` and the steps after it, whenever the report is drawn again (after a
MetaFlux update that changes the report, for example), and it then needs the line as
well. Never put the line in a config shared by several runs. A listed code that
names no model of the run gets a note and changes nothing.

**Where the decision is recorded.** In the verdict file and in the run config. In a
finished run, `"converged": false` in `stats/dada2/error_model.json`, or "Converged: no"
in the MultiQC summary, therefore means a model accepted by its code; in a run that
stopped at the check, the verdict file is missing. The report and the MultiQC section are
made before the check, so they never say "accepted".

**Another type of computer.** The loess fit runs Fortran code that calls the BLAS
linear-algebra library, and the OpenBLAS that conda installs picks its routines by
processor type. Learning the same reads with the same software on another type of
computer may therefore give a model that differs in its last digits, with another code;
the run then stops again for a look. This has not been tested yet.

**What a stop leaves on disk.** The two models, the record of their learning and the
whole report are written before the check, so they are there. Started with
`--keep-going`, MetaFlux also builds the MultiQC report, with its "DADA2 error model"
section. Nothing is deleted: if an earlier run had finished in this folder, its tables in
`05.dada2/` and `06.taxonomy/` stay, made with the earlier model, until the run goes on.

## The report

The report shows how each model was learned and how well it describes the counts it was
learned from. Open it when the run stops at the check, or to look over a run that did not
stop; for a first look, [What to expect](#goodness-of-fit) and the warning signs after
the reference table say what a healthy model looks like.

`error_model_report` writes the report on the two models in every amplicon run, whether
they converged or not, as tables and figures in `stats/dada2/` and as a section of the
MultiQC report. It is the evidence for accepting a model that did not converge, and the
same evidence for judging one that did. For each read direction it takes the transition
counts the model was fitted to and refits every error function MetaFlux offers on them.
This takes seconds and learns nothing new, so the refit of the function used gives back
the model exactly (the report checks that, as `refit_reproduces_model`). The report states
what this comparison can show, in the same words in `error_model_report.json`, in the
MultiQC section and here:

!!! note "What the comparison can and cannot show"

    The comparison refits only the smoothing step: each error function MetaFlux offers is
    applied to the transition counts that the model used was fitted to (for a model ended
    by the flat-rounds rule, those of the round before the last). It shows how closely
    each function describes those counts, not which function would have given the better
    ASVs: learning with another function would have produced other counts, and a more
    flexible function always comes closer to the counts it is fitted to. The numbers
    compare functions within this run only; they cannot judge a run, nor tell faulty reads
    from good ones.

### The figures

Each figure is written as a PNG, for the MultiQC report and for quick viewing, and as a
PDF.

**The fit figures**, `error_model_fit_R1.png` and `error_model_fit_R2.png`: one panel per
kind of misread, the consensus quality score along the bottom, the error rate up the side
on a log scale.

- **Grey points** are the rates counted in the reads `learnErrors` drew, sized by the
  bases at that score. They are *filled* at the quality values the reads carry and
  *open* at the averaged qualities of identical reads, which only binned data produce.
  Small *triangles at the bottom* mark scores with bases but no misread of that kind: on
  a log scale a zero cannot be drawn, and leaving it out would make the counted rates
  look higher than they are.
- **The solid line** is the model used, in the colour of its error function. **Dashed
  lines** are the other error functions refitted on the same counts. The rate axis spans
  the counted rates, the model and the nominal line; a dashed line that runs outside it
  (a refit can drop by orders of magnitude where the counts are thin) is cut off at the
  edge of the panel, and the caption says so.
- **The black dotted line** is the nominal rate of the quality score, the rate the score
  promises on paper (Q30: one error in 1,000 bases, a third of it for each of the three
  wrong bases).
- **The ticks along the bottom** mark the quality values found in the filtered reads; on
  binned data, values the binned fit left out get grey ticks.
- **The title** names the error function, with the anchors of a binned fit, and how the
  learning ended, for example
  `R2 error model: loess fit - converged (cycle of 11 rounds, rounds 22-32)` or
  `R2 error model: loess fit - NOT converged: flat-rounds rule, round 42 (code R2:7ebe3bed)`.

![The fit figure of the R1 model of the 16S test set: twelve panels, one per kind of misread, with the error rate on a log scale against the consensus quality score from 7 to 38. Grey points, larger at high quality scores where most bases are, fall with the quality score; the solid blue loess curve follows them smoothly; the dashed orange binned refit joins the points one by one and zigzags; a dotted black line marks the nominal rate of each score.](../assets/error_model_fit_16S_R1.png)

*The fit figure of the R1 model of the 16S test set (six samples, MiSeq 2×300, 32
distinct quality values). The loess fit, the model used (solid blue), follows the counted
rates and falls as the quality score rises; at the high scores that hold most bases it
sits above the nominal rate (dotted), as every model tested does. The binned fit refitted
on the same counts (dashed orange) joins the counted rates point to point and follows
their noise, which is why `error_model: binned` stops on continuous quality values; in
the C2T panel it drops out of the panel at Q7, where no C was misread as T, and is cut
off at the panel's edge.*

**The convergence figure**, `error_model_convergence.png` (shown under
[How learning ends](#how-learning-ends)): one panel per read direction, the round along
the bottom, on a log scale up the side the change in each round (black) and how far each
round's starting table is from the final model (grey, dashed; DADA2's
`checkConvergence`). Both are sums of rate differences, so they share one axis.

- A dashed line marks 1/1000 of the change in round 1, the depth the flat-rounds rule
  needs.
- A value of 0 cannot be drawn on a log scale, so its line breaks there and a marker at
  the bottom shows it: a triangle for DADA2's repeat; for the flat-rounds rule a square,
  because there the last distance is 0 only since MetaFlux handed back the table of that
  round.
- A grey band marks DADA2's repeat (labelled with the round repeated and the length of
  the cycle) or the last `max_flat_rounds` rounds before a flat-rounds stop.
- A dashed vertical line marks a flat-rounds stop, a solid one the cap; when the cap was
  not reached, the panel's title says so.

**The quality-values figure**, `error_model_quality_values.png`: the bases at each quality
value in the first 100,000 filtered reads of each voting sample, pooled, per read
direction, with the number of values above Q2 and the function chosen in the panel's
title (for example `R1: 32 values above Q2 (7-38) -> loess (auto)`). Values at or below Q2
are grey; with ten values or fewer, each is labelled with its share of the bases.

**DADA2's own plots**, `filtered_read_R1_error_plot.png` and `_R2_`: the picture DADA2's
`plotErrors` draws of each model (points: counted rates; black line: the model; red line:
the nominal rate), as MetaFlux has drawn it since v1.0.0, now with the same titles as the
fit figures. For a model ended by the flat-rounds rule, the points are the counts the
kept table was fitted to, those of the round before the last, as in the fit figures.

A figure the report cannot draw is replaced by a picture that states the error, and the
problem is listed in the report's warnings, which the check repeats on the terminal. The
report fails only when the record or a model cannot be read, so that a problem with a
figure never blocks the denoising.

### The tables

- `error_model_fit.tsv`: one row per read direction and error function, with the measures
  of goodness of fit defined below and a `note` where a row's numbers mislead.
- `error_model_rates.tsv`: one row per read direction, kind of misread and quality score:
  whether the score is a value the reads carry (`read_value`), an averaged one (`averaged`)
  or has no bases (`no_data`), the bases and misreads counted, the counted rate, the
  nominal rate, the rate of the model and the rate of each refit. The numbers behind the
  fit figures.
- `error_model_rounds.tsv`: one row per read direction and round, with the change, the
  lowest change so far, the distance to the final model, 1/1000 of round 1's change,
  whether the round lies in the window of the flat-rounds rule, and the event of the round
  (the round repeated, DADA2's repeat, the flat-rounds stop or the cap). The numbers behind
  the convergence figure.
- `error_model_report.json`: the report in one file, with the statement above, how each
  learning ended, the model codes and the report's warnings.

### Goodness of fit

For each read direction and error function, `error_model_fit.tsv` gives these measures.
Write N for the transition counts, E for the fitted table, and n for the bases of one true
nucleotide at one quality score (the sum of its four transitions).

| Column | What it measures |
|---|---|
| `loglik_per_base` | The log-likelihood of the counts under the fitted rates, Σ N·log E over the cells with counts, divided by the bases. Closer to 0 is closer. |
| `deviance_per_million_bases` | The multinomial deviance, 2·(Σ N·log(N/n) − Σ N·log E), per million bases: 0 when the fitted rates equal the counted ones. Given in total and split into the columns at the quality values the reads carry (`_at_read_values`) and the averaged columns (`_at_averaged_columns`); the two parts add up to the total. |
| `typical_fold_difference` | How many times the fitted rate is above or below the counted rate, typically: 10 to the power of the mean distance between the two on the log10 scale, over the misread cells with at least one misread, weighted by the bases n. `1.00` means the curve sits on the counted rates. Given in total, at the read values and at the averaged columns. |
| `rising_steps` of `steps_compared` | Places where the fitted rate rises by more than 10 % from one quality score to the next, over the 12 kinds of misread, between the lowest and the highest score with bases. A model should give a higher quality score fewer errors. `substitutions_with_rise` counts the kinds of misread with at least one rise. |

No column ranks the functions, and no row is marked best: each function wins on what it
is built for. The binned fit passes through the counted rates at the quality values the
reads carry, by construction, so it is always closest there; on continuous data it
therefore joins some thirty counted rates point to point and rises with quality in about
a third of the steps. The loess fit smooths across all columns, so it sits a little away
from the counted rates; on binned data it rises with quality in a sixth to a quarter of
the steps. MetaFlux chooses the function by the kind of quality values (DADA2's author
recommends the binned fit for binned values) and by whether the rates fall as quality
rises, not by closeness to the counts.

**Reference values.** The models the report was checked on, with the function used (or,
for the first NextSeq 2000 models, the function `auto` now picks) marked ★. The MiSeq
models (the eight MiSeq test sets, a single 16S sample, and a COI library with a fault)
and the first NextSeq 2000 models were learned with the loess fit, so their binned rows
are refits on the counts of a loess model. One sample of the same NextSeq 2000 run,
learned alone, and the Zymo D6331 mock under `auto` were learned with the binned fit, and
their loess rows are the refits; the Zymo mock was also learned with `error_model: loess`
(last two models), whose binned rows are again refits. The third column divides the
misread rate of the model (all three wrong bases together) by the nominal rate, at the
quality score that holds the most bases. The second column gives the share of the
counted bases that sit at quality values the reads carry; the rest sit at averaged
qualities, which only binned data produce. The last two columns are the measures defined
in the table above. Where only one typical difference is given, the model has no averaged
columns.

| Model | Bases at the reads' values | Model ÷ nominal rate at the commonest score | Function | Rising steps | Typical difference: all / read values / averaged |
|---|---|---|---|---|---|
| 16S test set, six samples, R1 (MiSeq) | 100 % | 7.4 (Q37) | loess ★ | 20 of 372 | 1.55× |
| | | | binned | 118 of 372 | 1.00× |
| 16S test set, six samples, R2 | 100 % | 8.3 (Q38) | loess ★ | 6 of 372 | 1.22× |
| | | | binned | 117 of 372 | 1.00× |
| 16S test set, one sample learned alone, R1 | 100 % | 21.3 (Q38) | loess ★ | 30 of 372 | 1.49× |
| | | | binned | 142 of 372 | 1.00× |
| ITS test set, six samples, R1 (MiSeq) | 100 % | 7.6 (Q38) | loess ★ | 17 of 372 | 1.47× |
| | | | binned | 142 of 372 | 1.00× |
| ITS test set, six samples, R2 | 100 % | 8.1 (Q38) | loess ★ | 5 of 372 | 1.27× |
| | | | binned | 109 of 372 | 1.00× |
| 18S test set, four samples, R1 (MiSeq) | 100 % | 8.7 (Q38) | loess ★ | 12 of 372 | 1.18× |
| | | | binned | 136 of 372 | 1.00× |
| 18S test set, four samples, R2 | 100 % | 8.3 (Q38) | loess ★ | 2 of 372 | 1.11× |
| | | | binned | 92 of 372 | 1.00× |
| gyrB test set, five samples, R1 (MiSeq) | 100 % | 7.9 (Q38) | loess ★ | 27 of 348 | 1.45× / 1.45× / 3.16× |
| | | | binned | 108 of 348 | 1.00× / 1.00× / 1.84× |
| gyrB test set, five samples, R2 | 100 % | 11.6 (Q38) | loess ★ | 14 of 372 | 1.35× |
| | | | binned | 110 of 372 | 1.00× |
| rpoB test set, five samples, R1 (MiSeq) | 100 % | 16.8 (Q38) | loess ★ | 40 of 372 | 1.41× |
| | | | binned | 113 of 372 | 1.00× |
| rpoB test set, five samples, R2 | 100 % | 8.6 (Q38) | loess ★ | 33 of 372 | 1.60× |
| | | | binned | 130 of 372 | 1.00× |
| COI insect mock, BF3/BR2 library, R1 (MiSeq) | 100 % | 11.2 (Q37) | loess ★ | 26 of 372 | 1.42× |
| | | | binned | 112 of 372 | 1.00× |
| COI insect mock, BF3/BR2 library, R2 | 100 % | 14.2 (Q38) | loess ★ | 19 of 372 | 1.15× |
| | | | binned | 131 of 372 | 1.00× |
| COI insect mock, Leray-XT library, R1 (MiSeq) | 100 % | 11.7 (Q38) | loess ★ | 34 of 372 | 1.35× |
| | | | binned | 132 of 372 | 1.00× |
| COI insect mock, Leray-XT library, R2 | 100 % | 13.5 (Q38) | loess ★ | 15 of 372 | 1.25× |
| | | | binned | 116 of 372 | 1.00× |
| COI insect mock, fwhF2/fwhR2n library, R1 (MiSeq) | 100 % | 7.2 (Q37) | loess ★ | 36 of 372 | 1.78× |
| | | | binned | 134 of 372 | 1.00× |
| COI insect mock, fwhF2/fwhR2n library, R2 | 100 % | 11.0 (Q38) | loess ★ | 37 of 360 | 1.44× / 1.44× / 5.39× |
| | | | binned | 117 of 360 | 1.00× / 1.00× / 2.30× |
| COI insect mock, reads of both strands in one file, R1 (stopped after 30 rounds) | 100 % | 25.4 (Q38) | loess ★ | 32 of 372 | 1.35× |
| | | | binned | 123 of 372 | 1.00× |
| 16S, NextSeq 2000 (9, 20, 34), R1 | 91.8 % | 24.7 (Q34) | loess | 77 of 300 | 1.40× / 1.13× / 14.4× |
| | | | binned ★ | 0 of 300 | 1.46× / 1.00× / 98.4× |
| 16S, NextSeq 2000 (9, 20, 34), R2 | 67.6 % | 20.1 (Q34) | loess | 68 of 300 | 3.26× / 2.07× / 9.90× |
| | | | binned ★ | 0 of 300 | 3.95× / 1.00× / 113× |
| 16S, NextSeq 2000 (9, 20, 34), one sample learned alone with the binned fit, R1 | 92.4 % | 1.7 (Q34) | loess | 76 of 300 | 1.10× / 1.04× / 2.39× |
| | | | binned ★ | 0 of 300 | 1.14× / 1.00× / 6.46× |
| 16S, NextSeq 2000 (9, 20, 34), one sample learned alone with the binned fit, R2 | 64.7 % | 2.0 (Q34) | loess | 75 of 300 | 2.06× / 1.40× / 4.20× |
| | | | binned ★ | 0 of 300 | 3.41× / 1.00× / 33.5× |
| 16S V3-V4, Zymo D6331 gut mock, NextSeq 2000 (9, 20, 34), one sample, binned fit chosen by `auto`, R1 | 88.2 % | 1.2 (Q34) | loess | 52 of 300 | 1.18× / 1.09× / 3.16× |
| | | | binned ★ | 0 of 300 | 1.20× / 1.00× / 10.7× |
| Zymo D6331 gut mock, binned fit chosen by `auto`, R2 | 65.1 % | 1.8 (Q34) | loess | 67 of 300 | 1.96× / 1.54× / 4.26× |
| | | | binned ★ | 0 of 300 | 1.92× / 1.00× / 16.1× |
| Zymo D6331 gut mock, `error_model: loess`, R1 | 88.2 % | 1.2 (Q34) | loess ★ | 54 of 300 | 1.17× / 1.09× / 3.08× |
| | | | binned | 0 of 300 | 1.16× / 1.00× / 8.88× |
| Zymo D6331 gut mock, `error_model: loess`, R2 | 65.1 % | 1.3 (Q34) | loess ★ | 66 of 300 | 1.88× / 1.56× / 3.94× |
| | | | binned | 0 of 300 | 1.74× / 1.00× / 16.4× |

On binned data, most of the bases of the averaged columns lie below both curves: they
come mostly from reads identical to their ASV. How many of the misreads these columns
hold, and so how many open points lie above the curves, depends on the data. In the
first NextSeq 2000 models of the table (two samples learned together for R1, one sample
for R2), they held 8.2 % and 32.4 % of the bases but only 1.4 % and 2.8 % of the counted
misreads, and both fits sat many times above their counted rates, the binned fit more,
since it draws straight lines between the reads' values. In the sample of the same run
learned alone, they held 7.6 % and 35.3 % of the bases and 20.8 % and 28.4 % of the
misreads: most of their bases still lie below the binned curve, but most open points lie
above it (163 of 212 cells with misreads in R1, 173 of 242 in R2), among them 60 of the
70 and 92 of the 93 below Q20, which together hold less than a tenth of a percent of the
bases. A
loess model learned from the same sample gave the same picture (17 % and 27 % of the
misreads at the averaged columns), so the data cause it, not the fit. On the Zymo D6331
mock, the averaged columns held 11.8 % and 34.9 % of the bases and 16.5 % and 23.0 % of
the misreads under the binned model (16.4 % and 22.3 % under the loess model).
`share_substitutions_at_averaged_columns` in `error_model_report.json`, and the notes in
`error_model_fit.tsv`, give the shares of your run.

**What to expect.** Under the loess fit, a smooth curve through the filled points of
each panel that falls as the quality score rises. In several panels it turns up a little
at the top: on the models of the eight MiSeq test sets above, 289 of the 342 rising steps
of the loess fit lay among the five highest quality scores (34 to 38). Under the binned fit, straight
pieces that meet at the quality values the reads carry (drawn bowed on the log scale),
passing through the filled points, with the open points of the averaged columns on both
sides: the large ones, which hold most of their bases, mostly below; small ones often
above, at the low scores most of them.

**Warning signs.** A solid line that misses most of the filled points, or one that rises
with the quality score over a wide range of scores, or clearly more often than on the
models above. A line above the nominal line is not one: at the quality score that holds
the most bases, every model tested gave 1.2 to 25 times the nominal rate (1.2 to 2.0 on
the Zymo D6331 mock and on the NextSeq 2000 sample learned alone, with either fit; 7 to
25 on the others), healthy and faulty libraries alike.

### In the MultiQC report

The last module of `multiqc/multiqc_report.html`, "DADA2 error model", shows the report in
six sections: a summary (per read direction the error function and how it was chosen, the
quality values, the rounds and the cap, how the learning ended, whether it converged, the
change in rounds 1, 2 and the last, the bases learned from, the model code and the
report's warnings), the quality values, the convergence figure, the fit figures of R1 and
R2, and a goodness-of-fit table with the rising steps and the typical differences. Its
description is the statement above. MultiQC waits for the report but not for the check,
so with `--keep-going` a run that the check stops still gets its MultiQC report. A stop
before the report, in `choose_error_model` or inside `learn_errors`, leaves no MultiQC
report; the Falco reports of each sample in `stats/falco/` are there all the same. On a
slow library the MultiQC report therefore appears only once the learning has ended.

## Evidence and limits

- **The flat-rounds rule changes very little.** On three healthy libraries learned with
  the loess fit, a model stopped by the rule changed at most 0.04 % of the reads compared
  with the model the loop ended with (0.02 % at the default of 20 rounds). With the
  binned fit the rule acted on one of six replayed learning runs and changed at most
  0.02 % of the reads (0.007 % at 20 rounds); no replay reached the cap of 200 rounds.
- **The plots and numbers cannot tell faulty reads from a harmless slow loop.** The COI
  library with reads of both strands in one file flattened at a change of about 0.003 per
  round from round 4 on, like a healthy slow library; only the reads showed the fault.
- **A converged model is not necessarily a correct one.** A model fitted to a mixture
  describes the mixture, however quickly it converged.
- **The goodness-of-fit numbers compare the functions within one run.** They cannot rank
  candidate functions by the ASVs they would give; that needs data of known composition.
- **The closest fit does not give the best ASVs.** Eleven error functions (`loess`,
  `binned`, and variants of the loess fit proposed by DADA2 users and by nf-core/ampliseq,
  and a monotone GAM) each learned the models of three mock data sets with known
  sequences: the ZymoBIOMICS D6331 standard on binned quality values (NextSeq 2000, one
  sample), and five rpoB mock libraries on MiSeq quality values, as sequenced and
  re-binned in silico. On no set did the function closest to the counts give the fewest
  spurious ASVs. On the MiSeq rpoB data, a binned fit with every quality score as a bin
  had a deviance of 581, against 39,945 for the loess fit, and gave one spurious ASV more.
  On the re-binned rpoB data, the loess fit came closest, and the binned fit gave 54 ASVs
  instead of 59, with 0.4 percentage points more reads in true ASVs. Every function
  recovered the same true sequences; the variants weighted by log10 of the base count gave
  the most spurious ASVs on all three sets. Only one sample of binned data from an
  instrument was tested.
- **Splitting the learning off changes no result where the loess fit converges.** The new
  steps were run outside Snakemake, with the default seed and `nbases`, on two test sets.
  On the 16S test set, learning and denoising the filtered reads of v2.5.0 gave its four
  DADA2 tables byte for byte. On the BF3/BR2 library of the COI insect mock, the steps
  from the filtering on, starting from the primer-trimmed reads of v2.6.0, gave the same
  filtered reads, the R1 and R2 models converged in rounds 6 and 32 as in v2.6.0, and the
  four DADA2 tables were byte for byte those of v2.6.0. (The DADA2 steps of v2.6.1 are
  those of v2.6.0, which gave the tables of v2.5.0 on the 16S test set.) Full runs with
  Snakemake of the eight MiSeq test sets (16S, ITS, 18S, gyrB, rpoB and three COI
  libraries) then gave, byte for byte, the tables of v2.6.1: in `05.dada2/`, `seqs.fasta`,
  `seqtab_head_seqs.txt`, `seqtab_head_names.txt` and `read.counts`, the tables of the
  extraction and of the length filter; the taxonomy tables of `06.taxonomy/` (all but the
  raw SINTAX output); `stats/read_tracking.txt`; and the alignment and the tree of the
  16S run. The raw working files of Metaxa2 (`05.dada2/metaxa2/seqs.*`), of ITSx
  (`05.dada2/_itsx_tmp/`, and `05.dada2/itsx_extraction.summary.txt`, a copy of ITSx's
  summary) and of SINTAX (`06.taxonomy/_sintax_raw.tsv`) can differ between two identical
  runs of any version: they hold start times, Metaxa2's own classification varies, and
  ITSx and VSEARCH write their results in the order their searches finish. No table
  depends on what varies in them.

## Comparing error functions on your own data

The report's refits show how each function describes the counts of one model, not which
function gives the better ASVs. To see that on your own data, learn the models again with
the other function on a copy of the output folder, and compare the results:

1. **Copy the output folder with its file dates** (the copy needs the disk space of the
   folder):

    ```bash
    cp -a /data/run1 /data/run1_binned
    ```

2. **Make a run config for the copy**: a copy of the run config with `output.out_dir`
   pointing at `/data/run1_binned`, and the other function under `amplicon:` → `dada2:` →
   `learn_errors:`, inside the `amplicon:` block it already has:

    ```yaml
    error_model: binned
    ```

3. **Run it with `--forcerun choose_error_model`**, after a dry run with `-n`:

    ```bash
    snakemake --sdm conda --cores 16 --configfile run1_binned.yaml --forcerun choose_error_model
    ```

    Snakemake keeps its records of finished steps in the MetaFlux folder, under the paths
    of their outputs, so it has none for the copy, and a changed setting alone is not
    noticed there. `--forcerun choose_error_model` makes the choice again, and the models
    are learned again from it, with everything after them; the steps before are kept.

4. **Compare** the ASV tables of the two folders (`06.taxonomy/asv_table.txt`), or, on a
   mock community, how many of its members each run recovers.

A model learned with the other function can stop at the check like any other; accept it by
its code in the copy's run config.

## Your own error function (experimental)

MetaFlux can learn the two error models with an R function of your own instead of `loess`
or `binned`, or draw such a function in the report next to them. This setting is meant
for developing and testing error functions on real runs. It is experimental: its settings,
checks and outputs can change in a later release.

### Settings

```yaml
amplicon:
  dada2:
    learn_errors:
      error_model: custom
      custom_error_function:
        file: /data/error_functions/my_errfun.R   # the R file that defines the function
        name: my_errfun                           # the function's name in that file
```

- `file` is an absolute path, or one relative to the folder `snakemake` is started in,
  like every path of the config (not relative to the config file).
- With `error_model: custom`, both models (R1 and R2) are learned with the function.
- With `auto`, `loess` or `binned` and a function named in `custom_error_function`, the
  models are learned as usual, and the report refits the function on the same counts as
  the others, for comparison only.
- With `file` and `name` both empty, the default, there is no function of your own.

As with every setting MetaFlux asks for, the lines go inside the `amplicon:` block your
run config already has; a second `amplicon:` block would replace the first, and stops the
run at start-up.

Under `custom`, `choose_error_model` still reads the quality values of the filtered reads,
and records them for the report. Its stops about the data itself still apply, because no
error function can learn a rate per quality score from made-up values: fewer than three
values above Q2, or one sample without them. The stops that exist only because `auto`
chooses between `loess` and `binned` (mixed instruments, two binned tables, undetermined
values) become warnings, as for any function set in the config.

### What the function receives and returns

The contract is DADA2's own `errorEstimationFunction`, so the same function also works in
plain DADA2:

```r
source("my_errfun.R")
err <- learnErrors(filtered_files, errorEstimationFunction = my_errfun)
```

MetaFlux asks more of the table than DADA2 does (see [What MetaFlux checks](#what-metaflux-checks)),
so a function plain DADA2 accepts can still be stopped here.

- **Input:** one argument, `trans`, the transition counts of one round of learning. It is
  a matrix of whole numbers (an integer matrix in DADA2 1.34.0) with 16 rows named `A2A`,
  `A2C`, `A2G`, `A2T`, `C2A`, `C2C`, `C2G`, `C2T`, `G2A`, `G2C`, `G2G`, `G2T`, `T2A`,
  `T2C`, `T2G`, `T2T`, in this order (the true base, `2`, the base it was read as), and
  one column per consensus quality score, named `"0"`, `"1"`, `"2"` and so on: at least 41
  columns (scores 0 to 40), more when the reads carry higher scores. A cell counts the
  bases of that true base read as that base at that score, in the reads `learnErrors`
  drew.
- **Output:** a numeric matrix of the same shape. A cell is the probability that the true
  base of the row is read as the second base of the row at that quality score. Every value
  is above 0 and at most 1, and the four rows of each true base (`A2A`, `A2C`, `A2G` and
  `A2T` for A, and so on) sum to 1 at every score. Keep every misread rate between 1e-7
  and 0.25, as DADA2's own functions do (`loessErrfun` clamps its estimates to that
  range): MetaFlux stops only on a rate of exactly 0 and warns on a rate below 1e-7, but
  a rate far below what the counts show acts almost like a 0 (see *Why a rate of 0 stops
  the run* below). Row and column names can be left out: MetaFlux then copies those of
  `trans`. Attributes (labels R can attach to an object, such as a class, or a fitted
  model kept with the table) other than the row and column names are dropped before
  DADA2 sees the table: DADA2 compares the tables of two rounds with `identical()`, which
  also compares attributes, so one that changes from round to round would keep it from
  ever seeing a repeat.
- **The same counts, the same table:** DADA2 ends the learning only when the function
  returns a table identical, bit for bit, to a table that already went into a round. The
  function must therefore return the same table whenever it gets the same counts: no
  random numbers without a seed of its own (see [Reproducibility](#reproducibility)), and
  nothing that depends on the time or on earlier calls.
- **Calls:** once in DADA2's start-up pass (round 0) and once per round, for R1 and then
  for R2; the report calls it once more per read direction, for its refit. In round 0
  DADA2 denoises with every error rate set to 1 and one cluster per sample, so the
  round-0 counts count every difference from the most abundant sequence, real variants
  included, as a misread. The function must still return a valid table for them.

The counts of any model MetaFlux has saved have the shape your function receives (the
counts of the model's last round), so a function can be tried in R before a run:

```r
trans <- readRDS("/data/run1/05.dada2/error_model_R1.rds")$trans
rates <- my_errfun(trans)
```

### A minimal example

`.test/custom_error_functions/loess_wrapper.R` hands the counts to DADA2's default
function and returns its table unchanged:

```r
loess_wrapper <- function(trans) {
  dada2::loessErrfun(trans)
}
```

With this function, `error_model: custom` gave, on the six samples of the 16S test set
(default seed and `nbases`), byte for byte the error models (`.rds`), the model codes and
the four DADA2 tables of `error_model: loess`. A function of your own replaces the body
with its own fit.

The file can also hold helper functions and constants. Keep the function and everything
it needs in this one file: MetaFlux tracks only this file, so editing another file it
reads (a helper read with `source()`, a fitted object read with `readRDS()`) does not make
the models be learned again, and that file is not in the md5 sum of the record. A relative
path inside the file is read from the folder `snakemake` is started in, not from the
file's own folder.

MetaFlux reads the file into an R environment of its own (a separate set of named objects
in the R session, not a conda environment). What the file defines with `<-` stays there
and cannot replace MetaFlux's own variables or functions. Two things still reach
MetaFlux's script: `<<-` or `assign(..., envir = globalenv())`, which write into it, and
variables your function uses without defining them, which R then looks up among
MetaFlux's own. The learning step and the report define different ones, so such a
function can work in one and fail in the other, and it would not work in plain DADA2.
Define at the top of the file everything the function keeps or uses (for example
`cache <- NULL`: a `<<-` inside the function then changes that variable, not one of
MetaFlux's), use only `trans` and what the file defines, and never `assign()` into the
global environment. The file sees the packages loaded in R and those it loads itself with
`library()`.

### Two more examples

Two more files in `.test/custom_error_functions/` fit the error rates themselves. Each
keeps everything it needs in the one file, as described above, and also works in plain
DADA2.

**`binned_from_counts.R`** joins the error rates counted at the quality scores of the
counts with straight lines, with DADA2's binned fit, `makeBinnedQualErrfun`. Its bins are
the scores at which all four true bases have bases, plus the lowest and the highest score
with any bases, because `makeBinnedQualErrfun` stops when the counts hold bases outside
its bins:

```r
binned_from_counts <- function(trans) {
  quality_scores <- as.numeric(colnames(trans))
  # ... has_all_four: TRUE at the scores at which A, C, G and T all have bases
  bins <- quality_scores[has_all_four]
  scores_with_bases <- quality_scores[colSums(trans) > 0]
  bins <- sort(unique(c(min(scores_with_bases), bins, max(scores_with_bases))))
  makeBinnedQualErrfun(bins)(trans)
}
```

`makeBinnedQualErrfun` first appeared in DADA2 1.36.0, and MetaFlux pins DADA2 1.34.0, so
the file carries an unchanged copy of it, under DADA2's license; MetaFlux's own copy is
read by MetaFlux's scripts, not by a file of yours. The bins come from the counts, whose
scores are consensus qualities (see [What `learnErrors` fits](#what-learnerrors-fits)). On
MiSeq reads they are the 30 or so scores the reads carry, and the fit follows the counted
rates without smoothing; on the five rpoB mock libraries, the file learned, bit for bit,
the models of the binned fit with every quality score as a bin under
[Evidence and limits](#evidence-and-limits). On binned reads the bins also include the
scores in between, at which few misreads were counted (25 bins from 9 to 34 on the R1
counts of the NextSeq 2000 test set), so the result is not that of `error_model: binned`.
On counts with bases only at 9, 20 and 34, the file returns the table of the built-in fit,
bit for bit.

**`monotone_gam.R`** fits, for each kind of misread, the error rate against the quality
score with a GAM of the R package scam whose curve may only fall as the score rises
(`bs = "mpd"`, a monotone decreasing curve). The fit is binomial, with the bases at each
score as weights:

```r
counted <- data.frame(quality = quality_scores[has_bases],
                      rate    = misreads[has_bases] / bases[has_bases],
                      weight  = bases[has_bases])
curve_model <- scam(rate ~ s(quality, bs = "mpd", k = k), family = binomial(link = "logit"),
                    data = counted, weights = weight)
```

Around the fit, the file does what DADA2's own functions do: it holds the curve flat below
the lowest and above the highest score with bases, where no read gives a rate; it keeps
every misread rate between 1e-7 and 0.25; and it sets each read-correctly row (`A2A` and
so on) to 1 minus the three misread rates of its base. When fewer than six scores hold
bases of a true base, or scam fails, it falls back to a straight line on the log-odds scale
(a binomial GLM, used only when the line does not rise), then to one pooled rate. scam's
default curve has 10 pieces (`k = 10`) and stops with "Model has more coefficients than
data" when fewer than 10 scores hold bases of a true base; the file sets `k` from the
scores available.
On the counts of the test sets, scam fitted every misread: the counts of the NextSeq 2000
test set, whose reads carry three values, held bases at 25 or 26 scores of each true base.
The file needs scam, which `dada2_custom.yaml` provides (see [R packages](#r-packages)). On
the three mock data sets under [Evidence and limits](#evidence-and-limits), it learned, bit
for bit, the models of the monotone GAM there.

Both files returned tables that passed every check below, and the same table when called
again, on the counts of the last round of the eight MiSeq test sets, of the ZymoBIOMICS
D6331 sample and of the rpoB mock libraries re-binned in silico (20 tables). In plain DADA2 (`learnErrors` with `nbases` 1e7, on the
ZymoBIOMICS D6331 sample, two rpoB mock libraries and 50,000 reads of a COI library), both
learned models that converged in 4 or 5 rounds.

Error functions posted by DADA2 users, for example in DADA2's issues, follow the same
contract, DADA2's `errorEstimationFunction` (one argument, the counts; a table of the same
shape back), and work here unchanged, as long as the R packages they load are in
`dada2_custom.yaml` and their tables pass the checks below.

### What MetaFlux checks

At start-up, before any job: with `error_model: custom`, `file` and `name` must both be
set; with any other `error_model`, both or neither; the file must exist and be readable;
the name must be one R can give a function without quotes (the letters a-z and A-Z, with
no accented letters, digits, dots and underscores, starting with a letter or with a dot
not followed by a digit, and no word R reserves, such as `function` or `TRUE`).

When `learn_errors` runs, MetaFlux reads the file and fetches the function by its name,
then checks every table the function returns, before DADA2 uses it:

1. a numeric matrix with as many rows and columns as `trans`;
2. row and column names, when given, equal to those of `trans`, in their order;
3. no missing (`NA`, `NaN`) or infinite value;
4. every value from 0 to 1;
5. no value of exactly 0;
6. the four rows of each true base summing to 1, within 1e-6, at every quality score.

A rate below 1e-7, the lowest DADA2's own functions return, does not stop the learning:
it gives a warning (see *Why a rate of 0 stops the run* below). The first problem of the
list stops the learning, with a message that names the function, its file and the
problem, such as:

```text
[learn_errors] STOP: the R1 error function failed in round 0: the custom error function my_errfun (/data/error_functions/my_errfun.R) returned rows A2A, A2C, A2G, A2T that do not sum to 1 at 41 quality score(s), the first at quality 0 (sum 1.25). ...
```

An R error inside the function stops the learning in the same way, with R's own message;
DADA2 itself would report only "Error rates could not be estimated". The messages, with
their causes, are listed in [Troubleshooting](../troubleshooting.md#amplicon-learn_errors-stops-in-your-own-error-function).
A warning of the function, or MetaFlux's warning about rates below 1e-7, does not stop
the learning: each distinct text is written once per read direction, with the round in
which it first appeared, to `logs/learn_errors.log` and to the `warnings` of its read
direction in `stats/dada2/error_model.json`.

**Why a rate of 0 stops the run.** DADA2 accepts a rate of 0 with a warning ("Zero in
error matrix."), but reads it as a misread that cannot happen: a sequence that carries
that misread at that quality score can no longer be explained as a misread of a more
abundant one, and if it is seen in more than one read it can become an ASV of its own.
On the R1 reads of one sample of the 16S test set (28,648 reads, denoised with DADA2
1.34.0), the loess model gave 138 ASVs; the same model with `A2G` set to 0 at every
quality score gave 226, with `A2G` at 0 at Q37 alone 147, and with all 12 misreads at 0
at Q37 279. Under the loess model, 1,504 distinct sequences differed from their ASV by an
A-to-G misread; 120 of them were seen in more than one read, and the 0 added 88 ASVs.
DADA2's own error functions never return a 0: both keep every misread rate between 1e-7
and 0.25 and set each read-correctly row (`A2A` and so on) to 1 minus the three misread
rates of its base. A function of your own can do the same.

A rate that is not 0 but far below the counted one is not stopped, yet acts almost the
same. On the same reads, `A2G` at 1e-300 at every quality score gave the 226 ASVs of a 0,
at 1e-30 186 and at 1e-12 159; even 1e-7, the floor of DADA2's own functions, gave 146,
because the model's own `A2G` rates are all above 3e-4. MetaFlux therefore warns when a
table holds a rate below 1e-7, without stopping: a fitted curve can dip there at the
highest quality scores, or at scores without bases, where no read carries the rate.

Everything after the fit works as for `loess` and `binned`: the early stop, the cap, the
check before denoising and acceptance by model code.

### Reproducibility

- **The same function, reads and settings give the same models and model codes.** MetaFlux
  sets R's seed (`amplicon.seed`) once before learning, and `learnErrors` draws with R's
  random numbers the order in which it reads the samples, R2's draw continuing from R1's.
  MetaFlux saves R's random-number state before each call of your function and puts it
  back after it, so the samples drawn are the same whatever the function does with random
  numbers. Without this, a function that set a seed of its own changed the samples the R2
  model was learned from: with `nbases` 1e7 on the 16S test set, R2 was learned from one
  sample instead of four. R's kind of random numbers is put back too, after the file is
  read: a file can switch it (code written for parallel work often starts with
  `RNGkind("L'Ecuyer-CMRG")`), and the seed would then start another series of numbers.
  Without this, such a file made R1 be learned from two samples instead of four, and R2
  from three.
- **A function that draws random numbers must set its own seed** (for example `set.seed(1)`
  at its start), so that it returns the same table for the same counts. The report refits
  the function outside the learning; without a seed of its own, the refit returns another
  table, `refit_reproduces_model` is false in `error_model_fit.tsv`, and the report warns.
- **Editing the file makes steps run again.** Under `error_model: custom` the file is an
  input of `learn_errors`: an edit makes the models be learned again, with the report, the
  check, `dada_seqtab` and every step after it. A function named for comparison only is an
  input of the report alone: an edit draws the report again and, because the check waits
  for the report, the check, `dada_seqtab` and the steps after it run again too, with the
  same models. Writing the same path differently (`./my_errfun.R` instead of
  `my_errfun.R`, or as an absolute path) counts as a change too, because the path as
  written is a parameter of the steps that read the file. To try several versions of a
  function on a large data set, work on a copy of the output folder (see
  [Comparing error functions on your own data](#comparing-error-functions-on-your-own-data)).
  When the learning stops in a folder that ran before, only the models and the record are
  removed: the report, the verdict of the check, the tables of `05.dada2/` and every
  result after them, and the MultiQC report on disk still come from the earlier models
  until the learning succeeds again.
- **The record says which version of the function made a model.** Under `error_model:
  custom`, `stats/dada2/error_model.json` holds `custom_error_function`: the function's
  name, the file as written in the config and as resolved, the md5 sum of the file, and
  the R packages loaded while the file was read or the function ran, with their versions
  (`packages_loaded`; packages MetaFlux had loaded already are not listed). The log line
  of each read direction names the function, the file and the first 8 characters of the
  md5 sum.

### R packages

Your function runs in a conda environment of its own, `workflow/envs/dada2_custom.yaml`:
the DADA2 environment of the other DADA2 steps (`dada2.yaml`) with two R packages added
for error functions built on generalised additive models (GAMs: a smooth curve, here of
error rate against quality score, whose shape is learned from the counts):

- **mgcv** 1.9-4, which fits GAMs (`mgcv::gam`), with nlme, which mgcv needs;
- **scam** 1.2-22, which fits GAMs whose curve may only fall or only rise, for example
  error rates that never go up as the quality score goes up (`bs = "mpd"`, a monotone
  decreasing curve).

Next to these, the environment provides what `dada2.yaml` provides: dada2 1.34.0 and the
packages it loads (among them Biostrings, ShortRead and RcppParallel), ggplot2, jsonlite,
R's base packages (stats, splines, methods and the others) and MASS, Matrix and lattice.
A package it does not provide cannot be loaded, and the learning stops with R's message
that there is no package of that name. The versions above are those of the pin file used
on Linux (x86-64); on other machines conda picks them when it builds the environment.

**Used only for your own function.** `learn_errors` runs in it under `error_model:
custom`, and `error_model_report` whenever a function is named, also as a comparison;
every other step, and every run that names no function of its own, runs in `dada2.yaml`
and never builds it. On Linux (x86-64) it is installed from
`workflow/envs/dada2_custom.linux-64.pin.txt`, which holds every build of
`dada2.linux-64.pin.txt` plus mgcv, nlme and scam. A model learned in it therefore rests
on the same dada2, R and linear-algebra library as one learned in `dada2.yaml`. On the
16S test set, `error_model: custom` with the example function above gave in it, byte for
byte, the models of `error_model: loess` learned in `dada2.yaml`, and the report drawn in
it had the same tables and figures as one drawn in `dada2.yaml`. On machines without pin
files (macOS, Linux on ARM), conda solves the two environments separately and may pick
other builds of the packages they share. With a function named only for comparison, the
models are learned in `dada2.yaml` and refitted in `dada2_custom.yaml`, so the report may
then warn that the refit of the function used does not give back the model
(`refit_reproduces_model` false); the models themselves are not affected.

**Apple Silicon Macs and Linux on ARM.** conda-forge has no r-scam build for these
machines (`osx-arm64`, `linux-aarch64`), so the environment cannot be built there as
shipped: a run that names a function of its own stops before any job starts, when
Snakemake creates the environment, with `Could not create conda environment from
.../dada2_custom.yaml` and conda's message that `r-scam` is not available from the current
channels. If your function does not use scam, delete the line `- r-scam` from
`workflow/envs/dada2_custom.yaml` in the MetaFlux folder you run from; without it the
environment can be built on both (mgcv exists there). A function that needs scam cannot
run on these machines until conda-forge builds it for them. On Intel Macs and on Linux
(x86-64) the file can be built as shipped. This is also why the two packages are not in
`dada2.yaml`: there, r-scam would stop every amplicon run on these machines.

### In the report

- With `error_model: custom`, the solid line of the fit figures is your function, labelled
  `custom (my_errfun): model used`, and the plot titles name it, for example
  `R1 error model: custom (my_errfun) - converged (round 5 repeated its starting table)`.
  Its refit must give back the model, as for the other functions (`refit_reproduces_model`).
- With any other `error_model` and a function named, the fit figures draw it as a dashed
  curve, `custom (my_errfun): refitted on the same counts`.
- `error_model_fit.tsv` gets a row `custom` per read direction, with the note
  `experimental: your own error function my_errfun (<file as configured>, md5 <first 8
  characters>)`, and `error_model_rates.tsv` a column `rate_custom`.
- In the MultiQC report, the summary gets a column "Custom error function" with the
  function's name, the first 8 characters of the md5 sum of its file and whether the model
  was learned with it, and the goodness-of-fit table a row `R1 | custom (my_errfun)` per
  read direction.
- A comparison refit that fails (an R error, or a table that fails one of the checks
  above) stops nothing and is not a warning of the report. The figure caption says
  `custom (my_errfun): refit failed: <reason> (no curve)`, `refit_status` in
  `error_model_fit.tsv` gives the reason, and the goodness-of-fit table of MultiQC says
  `refit failed (see error_model_fit.tsv)`. Warnings the function gives during a refit
  are not shown. A file that cannot be read when the report runs is a warning of the
  report, and so is a refit of the function the models were learned with that does not
  give back the model.

The statement under [The report](#the-report) holds for your function too: how closely a
function describes the counts of one model says nothing about the ASVs it would give.

## Files

[Output files](../reference/output.md) describes every file of these steps and what reads
it: the two models in `05.dada2/` ([inside 05.dada2](../reference/output.md#inside-05dada2)),
the record, the choice, the quality values, the report and the verdict in `stats/dada2/`
([stats](../reference/output.md#stats)), the MultiQC section
([multiqc](../reference/output.md#multiqc)), and the logs `logs/choose_error_model.log`,
`logs/learn_errors.log`, `logs/error_model_report.log` and `logs/check_error_model.log`
([logs](../reference/output.md#logs-benchmarks-and-_aux)). The settings are in the
[configuration reference](../reference/configuration.md#dada2-tier-1-operational-parameters)
and in the comments of `config/config.yaml`.
