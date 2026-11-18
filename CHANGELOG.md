# Changelog

What changed in each MetaFlux release, newest first. Each entry says what changes
for a user, what runs again in an output folder made with the previous release, and
what to edit in an existing run config. The documentation site has the details; the
release notes on GitHub repeat the longer entries. Every release is archived on
Zenodo; the concept DOI in the README always resolves to the newest one.

## 2.7.0 — [To be settled before release: the date]

### DADA2's error model: learned in its own step, checked before denoising, shown in the report

DADA2 denoises with an error model per read direction, learned from the reads in rounds
until a round gives back a table it has seen before. Up to v2.6.1 `dada_seqtab` learned
the two models and denoised with them in one step, and went on with a model whatever way
its learning had ended. Now:

- **Learning is a step of its own, `learn_errors`.** It writes the two models as DADA2
  returned them (`05.dada2/error_model_R1.rds` and `_R2.rds`), which `dada_seqtab` reads,
  and a record of how their learning ended (`stats/dada2/error_model.json`): the rounds,
  the change of the model in each round, how the loop ended, and a code for each model,
  such as `R2:7ebe3bed` (the read direction and the start of an md5 sum of the model's
  numbers). R1 is learned before R2 with the same seed, as before, so the same reads give
  the same models.
- **`amplicon.dada2.learn_errors.max_consist` is 200** (it was 999): an upper limit on
  the rounds, not a speed setting. Of 194 learning runs on real data, 186 converged
  within 8 rounds and the slowest took 101.
- **An early stop, `max_flat_rounds: 20`.** On some libraries the model keeps changing by
  tiny amounts for a hundred rounds or more. MetaFlux now ends the learning when the
  change in a round is less than 1/1000 of the change in round 1 and the lowest change so
  far fell by less than 10 % over the last 20 rounds (`0` switches the rule off).
  Calibrated by replaying four learning runs round by round, all with the loess fit: it
  ended a 101-round run at round 42, where denoising that sample with the stopped model
  changed 0.02 % of its reads, and a COI library with reads of both strands in one file at
  round 32; over every value tested (5 to 50) the largest change was 0.04 % of the reads.
  Six learning runs replayed with the binned fit (next section) agree: the rule ended one
  95-round run at round 49, where the stopped model changed 0.007 % of the reads (at most
  0.02 % over every value tested), and no run reached the cap.
- **A check before denoising, `check_error_model`.** A model ended by the cap or by the
  early stop has not converged, and the numbers alone cannot tell a slow but healthy
  library from one with faulty reads (a library with reads of both strands in one file
  looked just like a healthy one). The run now stops before any denoising, prints what to
  check in the reads and in the report, and the line that accepts this model by its code:
  `accept_unconverged: ["R2:7ebe3bed"]` under `amplicon.dada2.learn_errors`. With the line
  in the run config the run goes on without learning again, and a later model with
  another code stops it again. The verdict is written to
  `stats/dada2/error_model_check.txt`. A stop deletes nothing: the tables of an earlier
  run in the folder stay until the run goes on.
- **A report on the models, in every amplicon run, `error_model_report`.** Before the
  check, so that it exists when the check stops the run: for each read direction a fit
  figure (the error rates counted at each quality score for the 12 kinds of misread, the
  model, the other error functions refitted on the same counts, and the nominal rate of
  the quality score), a convergence figure (the change in each round, with DADA2's stop,
  the early stop and the cap marked), the quality values of the filtered reads, DADA2's
  own error plots as before, and four tables with the goodness of fit of each error
  function, the rates and the rounds (`stats/dada2/error_model_*`). The same report forms
  a new MultiQC module, "DADA2 error model". The comparison refits only the smoothing step
  on the counts behind the model used, so it shows how closely each function describes
  those counts, not which function would give the better ASVs.
- **The documentation** has a new page, "DADA2 error model", and the troubleshooting entry
  "DADA2's error model never converges" is rewritten for the check.

For maintainers: the denoising script is now `60g_dada_seqtab.R` and the ASV length filter
`60h_dada_length_filter.py`; 60c to 60f are the new steps.

### Binned quality scores: the error function is chosen from the reads

NovaSeq, NextSeq, MiSeq i100 and other recent Illumina instruments write only a few
quality values (for example 9, 20 and 34), and DADA2's author recommends a binned error
function for them, `makeBinnedQualErrfun`, instead of the loess fit DADA2 uses by default.

- **`amplicon.dada2.learn_errors.error_model: auto | loess | binned`.** `dada_filter` now
  counts the quality values in the first 100,000 filtered reads of each file (a new file
  per sample, `stats/dada2/{sample}.quality_values.json`; `filter_stats.json` is
  unchanged), and the new step `choose_error_model` pools the samples with at least 1,000
  filtered reads. Under `auto`, 3 to 10 distinct values above Q2 with gaps between them
  give the binned fit, anchored on the values the reads carry, and 20 or more the loess
  fit, for R1 and R2 separately. The values are read from the reads, never from an
  instrument table. The choice goes to `stats/dada2/error_model_choice.json`, every
  sample's values to `stats/dada2/quality_values.tsv`.
- **The run stops before learning** when the reads carry fewer than three values above
  Q2 (as files re-encoded by ENA do), when one sample does, and, under `auto`, when
  samples of different instruments are mixed, or when 11 to 19 values, or 3 to 10 that
  include neighbours, leave the kind open. Each message says how to go on;
  `error_model: loess` turns the last two into warnings.
  `binned` set on continuous values stops; `loess` set on binned values warns and goes on.
- **Neighbouring quality values are not binned.** No instrument table has two
  neighbouring values (NovaSeq 6000: 11, 25, 37), but MiSeq reads filtered with a
  `trunc_q` of 28 to 30 keep only the 8 to 10 values above it, all neighbours; counted
  as 3 to 10 values, they would get the binned fit without a word. Under `auto` such
  values stop the run, with a message that names `trunc_q`; with a function set in the
  config they give a warning. Neighbours that only arise between samples of two binned
  tables (a NovaSeq X writes 40 or, on some software versions, 41) count as two binned
  tables.
- **`makeBinnedQualErrfun` ships with MetaFlux**, copied unchanged from DADA2 1.38.0 into
  `workflow/scripts/amplicon/vendor/` under DADA2's license, because the function first
  appeared in DADA2 1.36.0. MetaFlux stays on DADA2 1.34.0: `dada2.yaml` and its pin do not
  change, and nothing changes for continuous data, which keep the loess fit.
- **Experimental: your own error function.** `error_model: custom`, with the R file and
  the function's name under `custom_error_function`, learns both models with an R
  function of your own, written as DADA2's `errorEstimationFunction` (it takes the table
  of transition counts and returns error rates of the same shape, so it also works in
  plain DADA2); named under any other `error_model`, the function is only refitted and
  drawn in the report, next to loess and binned. The block is checked at start-up, and
  every table the function returns is checked before DADA2 uses it (shape, names, no
  missing value, every rate above 0 and at most 1, the rows of each base summing to 1; a
  rate below 1e-7 gives a warning); the record keeps the function's name, its file, the
  file's md5 sum and the R packages it loaded. With a wrapper of DADA2's `loessErrfun`,
  the 16S test set gave the error models, model codes and DADA2 tables of
  `error_model: loess` byte for byte.
- **Two example error functions** in `.test/custom_error_functions/`, shown under "Your
  own error function" in the documentation: `binned_from_counts.R`, DADA2's binned fit
  with the bins taken from the counts (it carries an unchanged copy of
  `makeBinnedQualErrfun`, under DADA2's license), and `monotone_gam.R`, a monotone
  decreasing GAM of the R package scam for each kind of misread. Both also work in plain
  DADA2.
- **A conda environment for your own error function, `dada2_custom.yaml`:** the DADA2
  environment plus the R packages mgcv and scam, for error functions built on GAMs
  (smooth curves fitted to the counts; scam can force a curve that never rises).
  `learn_errors` runs in it under `error_model: custom`, and `error_model_report` whenever
  a function is named; every other step, and every run that names no function of its
  own, stays in `dada2.yaml` and never builds it. On Linux its pin file holds every build
  of dada2's pin plus mgcv, nlme and scam, so a model learned in it rests on the same
  dada2, R and maths libraries. conda-forge has no r-scam for Apple Silicon Macs or Linux
  on ARM: there the file cannot be built as shipped, and a function that does not use
  scam runs once the line `- r-scam` is deleted from it. That is why the two packages
  are not in `dada2.yaml`. For maintainers: this pin is rebuilt from dada2's pin, not
  solved (`CONTRIBUTING.md` says how).
- **The default of `error_model` is `auto`:** binned data get the binned fit, so their ASV
  tables change on purpose; `error_model: loess` keeps the fit of v2.6.1. On three mock data
  sets, eleven error functions recovered the same true sequences, and the function closest
  to the counts never gave the fewest spurious ASVs, so `auto` chooses by the kind of
  quality values, not by fit.

### Falco on the filtered reads

Falco now also profiles the reads DADA2's `filterAndTrim` keeps, a fourth stage,
`filtered` (`stats/falco/{sample}_R{1,2}_filtered/`), shown in MultiQC as
`QC | quality-filtered reads`: the reads the error models are learned from. MultiQC shows
the mean quality along the read, which hides binned values; each sample's own Falco page
draws them as box plots on their common values, and the error-model report lists them
all. One more Falco job per sample and read direction, seconds each.

### The samples the error models learn from: the same in every language setting; R1 paired with R2 by name

- **`learn_errors` sorts the filtered files byte by byte.** `learnErrors` learns from
  whole samples drawn in random order from the sorted list, so the order decides which
  samples the model is learned from. Up to v2.6.1 the list was sorted alphabetically in
  the language setting (locale) of the shell, and for sample names that mix capitals and
  small letters, or `-`, `_` and `.`, that order differs between `LANG`/`LC_ALL=C`
  (common on clusters and in containers) and a UTF-8 setting: learning a folder again
  from another shell could give another model, and now, with models checked by their
  code, a second stop after a code had been accepted. For such names a run can draw other
  samples than v2.6.1 did in a UTF-8 shell; names that sort alike in every setting, as
  those of the test sets do, draw the same samples as before.
- **`dada_seqtab` pairs each R1 file with the R2 file of the same name.** It paired the
  two sorted lists by position, which under `LANG`/`LC_ALL=C` swapped the R2 files of
  samples named like `Soil` and `Soil_R2`: the job then stopped with DADA2's error
  "Non-corresponding derep-class and dada-class objects.", or, when the two samples had
  the same number of filtered reads, merged each sample's R1 reads with the other's R2
  reads without an error. The row order of the tables is unchanged; it still follows the
  language setting, as up to v2.6.1, and so does the numbering of ASVs with equal read
  counts.

### The 16S length probe: products that start at a look-alike of the forward primer

The probe that sizes the amplicon (`expected_length: auto`) cuts the products out of SILVA
with Cutadapt, which keeps the leftmost place where the forward primer matches within the
error rate. The 17-mer 341F (`CCTACGGGNGGCWGCAG`, the V3-V4 primer of Illumina's 16S
protocol) also matches, with 2 mismatches, a stretch in V2 of *Pseudomonas*,
*Acinetobacter*, *Stenotrophomonas* and their relatives, so on SILVA 144 6.3 % of the
341F/805R products started there, 583 bp long instead of 427 bp, and put p95 at 582 bp.
Every 2 x 300 run with these primers stopped at `pick_trunclen` on the defaults, although
its reads were long enough.

- **The 16S marker pack now sets `probe_exclude_second_site: true`**, the third pass the
  rpoB pack already had: products that still hold a second whole forward-primer site are
  left out of the length statistics (they stay in the amplicons FASTA), for both SILVA
  releases. On SILVA 144, 341F/805R: p95 582 -> 429 bp; 799F/1175R (V5-V7, the 16S test
  set): 0.69 % of the products, p95 386 -> 385 bp; 515F/806R (V4): 0.15 %, p95 255 bp
  either way. On SILVA 138.2, 799F/1175R: 0.61 %, p95 384 bp either way, q1 375 -> 374 bp
  (the lower edge of the `auto` ASV length window moves from 325 to 324 bp). The probe
  JSON gains `n_excluded_second_site`, and the cache name ends in `_excl-second-site`.
  For 341F/805R the `auto` ASV length window becomes [353, 479] bp instead of [353, 632].
- **The stop at the overlap check names the numbers and the settings.** Under each
  `resolve_policy` the message now says how far the reads reach after primer removal
  (R1 + R2), what merging needs (expected amplicon length + `amplicon.min_overlap`), where
  the expected length came from (`amplicon.expected_length`, or the probe statistic with
  its exact key, `amplicon.probe_length_stat.16S` for 16S and gyrB), and what to check.
  It used to give only the deficit and two headrooms.

### Automatic truncLen: what the cut is based on, short off-target products, and a second method

`pick_trunclen` picks one `truncLen` per read direction in `trunc_len.mode: auto`. Three
changes:

- **The log says what the cut is based on.** It names the quality values of the
  primer-trimmed reads (read from the first 10,000 reads of each file) and their kind,
  binned or continuous, by the same rule as the choice of the error function; what
  `q_threshold` means on them (on 9, 20, 34 data, the default 20 cuts where at least a
  quarter of the bases carry 9); and whether each cut came from quality or from the
  read-length limit. It warns when `q_threshold` equals one of the binned values, as the
  default does on NextSeq 1000/2000 600-cycle data, or lies above the highest. No cut
  changes.
- **A short off-target product no longer sets the read-length limit.** The length that
  `min_read_coverage_pct` of a sample's reads reach is now counted over the reads at least
  half the expected amplicon length long; shorter reads come from a primer dimer or an
  unspecific product, and the ASV length window would remove them anyway. Up to v2.6.1 one
  sample with more than 5 % of such reads could cap the whole run too short to merge, and
  the run stopped at the overlap check: on a 78-sample NextSeq 16S soil study, 13 samples
  carried 5 to 40 % of their reads in a 150-165 bp product. The study now gets 272/280.
  The log lists every sample's limit and the share of its reads left out. A sample with
  more than half of its reads below the floor keeps the old limit, with a warning. ITS is
  unchanged. On the eight MiSeq test libraries every `truncLen` is the same as before.
- **A second method, `trunc_len.method: simulate`** (the default stays `quality`). It
  replays `filterAndTrim` with the run's own `trunc_q`, maxN and `max_ee` on 20,000 read
  pairs per sample spread over each whole file, finds for each read the longest cut at
  which it still passes, and takes the R1/R2 cuts that keep the most pairs while leaving
  `expected_length + min_overlap + simulate_extra_overlap` (new key, 20 bp) bases; of the
  cut pairs within 0.5 points of the most, the longest. It reads the same quality values
  as the filter, so binning does not matter to it, and it sees reads with an N at a cycle
  that varies from tile to tile, which no quality summary shows. The quality rule's cut is
  logged next to it, with the share of the sampled pairs each keeps. `trunclen.json` keeps
  its fields (`resolved_via: simulate`). ITS stays with the quality rule. On the test
  sets it changed the cuts of every library but ITS; the section "The simulate method"
  of "Amplicon length and truncation" gives the cuts, the reads kept and the ASVs on
  mock communities.
- **The documentation** ("Amplicon length and truncation" and the configuration reference)
  describes all three.

### Primer trimming: the extra last sequencing cycle is removed from untruncated reads

An Illumina run sequences one cycle more than the read length it is set up for (2 × 151,
2 × 251, 2 × 301): the extra cycle serves the phasing correction of the cycle before it,
and its own calls are poor. On the MiSeq test sets the last base of R2 agrees with the
base the other mate reads at the same place of the template in 60 to 86 % of the pairs,
against 96 to 98 % one cycle earlier; in the ITS2 yeast mocks, where the template carries
a T, R2's cycle 301 reads C in 79 % of the pairs. ITS reads are not truncated, so that
base reached DADA2, survived denoising where most reads of a sequence carry the same
wrong call, and made `mergePairs`, which accepts no mismatch in the overlap, reject the
pairs of that sequence.

- **`amplicon.cutadapt.trim_final_cycle: auto | true | false`** (new key, default `auto`).
  `trim_primers` removes the last base of every read before it searches for the primers
  (Cutadapt's `--cut -1` and `-U -1`, in the 5′ call of both passes). `auto` does so for
  each sample and read direction that DADA2 does not truncate (ITS; the other markers only
  with `trunc_len.mode: manual` and a cut of 0) and whose longest read, among the first
  10,000, is 76, 101, 151, 251 or 301 bp; `true` always; `false` never, as up to v2.6.1.
  The decision, with the longest read and the reason, opens
  `logs/trim_primers/{sample}.log`. It is made by a new program,
  `workflow/scripts/amplicon/30a_final_cycle.py`, an input of the rule.
- **ITS2 yeast mocks** (six samples): with the removal alone, retention rises from 16.0
  to 33.7 % of the raw read pairs, and *S. paradoxus* (4 samples) and the *S. uvarum*
  group (2 samples), whose pairs did not merge, are found. With the new window floor for
  ITS (next section) as well: 18 ASVs instead of 12, retention 62.2 %, and every expected
  ITS2 type found in every sample (26 of 26, against 12). `max_mismatch: 1` on top adds
  2 ASVs, 2.8 points of retention and no expected type, so its default stays 0.
- **Why `auto` leaves truncated reads alone.** Forced on the other MiSeq test sets, the
  removal changed no final table of 16S and gyrB (the cut removes the base anyway) and
  cost reads where the cut sits at the read-length limit: on 18S the R1 cut moved from
  271 to 274 bp, and the final table lost 1.2 % of its reads and 38 of its 1,262 ASVs,
  with merging at 98 % either way; on a copy of the rpoB set with binned quality scores,
  one of the 19 strains (42 reads).
- **Primer trimming of the ITS2 test library** (20,000 pairs of SRR20648181) keeps the
  same pairs with and without the cut, and no trimmed read holds cycle 301. 45 % of the
  reads, those that run into the far primer, are unchanged and 50 to 55 % are one base
  shorter; 0.2 % of the R1 and 4.6 % of the R2 reads end a few bases earlier or later,
  because the chance matches of a few bases of the far primer at the very end of a read
  now fall on other reads. Where the base is kept, the trimmed reads are those of v2.6.1,
  byte for byte.
- **The documentation** (the amplicon overview, "The extra last cycle"; the ITS page and
  its new "Validation" section; the test datasets; the configuration and output
  references) describes it.

### The ASV length window for ITS

- **With ITSx extraction, the floor of the `auto` window is the shortest ITS region in
  the probe reference** (38 bp for ITS2, 17 bp for ITS1), without `window_margin`. ITSx
  reports a region only where it finds the 5.8S boundary, and a boundary it misses on
  the other side leaves rRNA on the region, so a short region is a short ITS, not a
  fragment. Up to v2.6.1 the floor was q1 minus the margin, as for the other markers
  (125 bp for ITS2), and it removed the yeasts with the shortest ITS2 of the ITS2 test
  set (63-78 bp), 64 % of the reads after chimera removal. The `source` of the window in
  `asv_length_stats.json` reads `auto:ITS_probe_min_p95+50`. The ceiling, the other
  markers and ITS without extraction are unchanged.

### COI: version 2 of the set of bacterial COI homologues

- **COI runs use `bacterial_cox1_gtdb-r232.v2`** (Zenodo DOI 10.5281/zenodo.24189850)
  instead of `gtdb-r232.v1`: deeper searches against MIDORI2, records 80-90 % identical
  examined in lineages that read TGA as tryptophan, and a chimera rule for animal COI
  records that carry a stretch of an endosymbiont gene. 242,696 whole genes and 31,709
  fragments; its MIDORI2 list examines 500 GB273 records and removes 236 (187 distinct
  sequences), against 186 and 170 in version 1. The release README lists the changes.
- **The MIDORI2 conversion removes records by sequence.** The list now gives the SHA-1 of
  each record's sequence, and `prepare_midori2_coi_sintax` leaves out every record that
  carries a removed sequence, under whatever name MIDORI2 files it. On GB273 it leaves out
  236 of 3,213,682 records, as the build of the set does. A list from before version 2
  (without `sequence_sha1`) stops the run. A list made from the same MIDORI2 file whose
  removed sequences are not all in it stops it too; a list made from another MIDORI2 file
  counts the removed sequences the file lacks in the log and the report, and the run goes
  on.
- **On the three libraries of the COI insect mock**, the specimens recovered and the reads
  of the kept ASVs within 97 % of a specimen are unchanged. The off-target guard moves 71
  ASVs instead of 69: two *Tisiphia* ASVs (5 and 228 reads), which version 1 left among
  the animal ASVs as *Mycetophila lunata* because MIDORI2 still held a record carrying a
  *Tisiphia* stretch under that name (KM679400.1). SINTAX names none of the 71 at species
  level; with version 1 it called 41 of them the stink bug *Picromerus bidens*, after a
  *Wolbachia* record filed under that name (KR036430.1), which version 2 removes.
- **What runs again in a COI output folder:** the three files of the set are downloaded
  into `refdb/coi/bacterial_cox1/gtdb-r232.v2/` and the converted MIDORI2 file is written
  again. Every step that reads that file runs again: the length probe (`probe_ref` of the
  COI pack), and with it, where `expected_length: auto`, `pick_trunclen` and every step
  after it; the reference check, the classification and the COI screens. The files of
  version 1 stay in `refdb/` and can be deleted.

### The messages of the COI steps say bacterial COI homologues

COX1 is the name of the mitochondrial gene. The genes of the bacterial and archaeal
reference set that the COI steps use are its homologues (*coxA* or *ctaD*, *cyoB*,
*qoxB* and *ccoN*, among others), together with the related nitric oxide reductases.
The messages and comments of the COI steps, the comments of `config/config.yaml` and of
the COI marker pack included, now call them bacterial COI homologues, as the manual and
the Zenodo record of the set do. Four messages read differently:

- **The report of the MIDORI2 conversion** (`…_metaflux_sintax.report.tsv`, next to the
  converted reference): the row `records left out: bacterial COX1 under a eukaryote
  name` is now `records left out: bacterial COI homologue under a eukaryote name`, with
  the same count.
- **The log of the conversion** (`logs/refdb/prepare_midori2_coi_sintax.log`):
  `N left out as bacterial COX1` is now `N left out as bacterial COI homologues`.
- **The stop of the taxonomy step** when the COI screens remove or move every ASV:
  `(N by the frame check, M as bacterial COX1)` is now `(N by the frame check, M as
  bacterial COI homologues)`.
- **The stop of `fetch_bacterial_cox1`** when COI.yaml names no URL or SHA-256 for the
  set, written to its log: `for the bacterial COX1 set yet` is now `for the set of
  bacterial COI homologues yet`.

File names, rule names, settings, output columns and the values of `filter_status` keep
the short form `cox1` (for example `bacterial_cox1_asvs.tsv`, `post_bacterial_cox1` and
`moved_bacterial_cox1`), so no run config, and no script that reads the tables of an
output folder, needs a change; only a script that looks up the row of the conversion
report by its name needs the new wording above. No count, threshold or decision changes.

### Raw reads named `.fq` or `.fq.gz`, and the R2 file of every sample

- **Raw reads with the extension `.fq.gz` or `.fq` work**, in both modes; some sequencing
  providers deliver reads named so. Up to v2.6.1 they passed the start-up check, but every
  step that reads the raw reads looked for `.fastq.gz` or `.fastq` only: an amplicon run
  stopped on a missing input, and a shotgun run under `_R1`/`_R2` naming planned an SRA
  download for each local sample. The steps now build the path from the extension found
  at start-up, and the raw read count of amplicon mode (`stats/01.raw_reads.counts`) takes
  the sample name from file names with any of the four extensions.
- **Every sample's R2 file is checked at start-up**: it must exist next to the R1 file,
  with the same extension, and the message names the files not found. Up to v2.6.1 only
  the R1 files were looked at, and a missing R2 file showed up later. In a shotgun run with
  local files under `_R1`/`_R2` naming it made Snakemake plan an SRA download for that
  sample (rule `download_sra` writes `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz`),
  and before a job runs, Snakemake deletes those of its outputs that exist: the sample's
  own R1 file. On a test sample with the code of this release before the change, the R1
  file was deleted and the download then failed. A folder with gzipped R1 and uncompressed
  R2 files (or the reverse), which earlier versions happened to read, now stops at
  start-up too.
- **Tested with Snakemake**: the 16S test set read as `.fq.gz` and as uncompressed `.fq`,
  and the shotgun Zymo D6331 library read as `SRR38976715_1.fq.gz` / `_2.fq.gz`, gave the
  result tables of the same runs on `.fastq.gz` files, byte for byte.

### Metaxa2 and ITSx with a long output path

- **The extraction (`target_extract`) works whatever the length of `out_dir`.** Metaxa2
  (16S) and ITSx (ITS) name their temporary folder after the output prefix MetaFlux gives
  them, with every `/` turned into `-`, and create it as one folder name in the folder they
  are started from. MetaFlux gave the full path, so with an `out_dir` longer than 210
  bytes (once made absolute) the name passed the 255 bytes Linux file systems allow for
  one name. Both tools then exited as if they had succeeded: Metaxa2 extracted nothing and
  the run stopped later in `dada_length_filter` with a bare `KeyError: 'min'`; under ITS,
  `target_extract` stopped with a message that blamed the primers or `amplicon.its_region`.
  Both tools now run inside their output folder (`05.dada2/metaxa2/`, `05.dada2/_itsx_tmp/`)
  with the short prefix `seqs`, so the name no longer depends on `out_dir`. The output
  files are those of before; the line `Running in <folder>: ...` of
  `logs/target_extract.log` names the folder.
- **Metaxa2's temporary folder is no longer made in the MetaFlux folder** (the folder
  Snakemake runs from) but in `05.dada2/metaxa2/`; both tools delete it when they end.
- **Tested with Snakemake**: the 16S and ITS test sets with an `out_dir` of 300 bytes gave
  the result tables of the same runs with a short `out_dir`, byte for byte.

### A manual truncation with the automatic length window, and two clearer stops

- **A manual truncation works with a length of your own and the automatic ASV length
  window.** With `trunc_len.mode: manual`, `expected_length` set to a number or a
  `[min, max]` range, and `length_filter.mode: auto` (the template's), `dada_length_filter`
  crashed with a bare `KeyError: 'expected_length'` after the whole DADA2 run. With no
  probe, its window is the expected length -15 % / +15 %, read from `stats/trunclen.json`,
  and in manual mode `pick_trunclen` did not write it there; the start-up check caught
  only the ITS form of this combination. `pick_trunclen` now records the expected length in
  both modes, resolved the same way (a number as it is, a range as its max, `auto` as the
  probe's statistic), so `trunclen.json` gains `expected_length` and
  `expected_length_source` in manual mode. A length it cannot resolve (for example
  `expected_length: null`) is recorded as `null`, with a note in its log, since the manual
  cuts do not need it; the length filter then stops with one `ERROR` line that names the
  settings to change.
- **`dada_length_filter` stops with an `ERROR` line instead of a Python traceback** when its
  settings leave it nothing to size the window from: the case above, and a
  `length_filter.mode` that is neither `auto` nor given a `range`.
- **The length probe says when all its products came from lineages it leaves out.** When
  every in-silico product of a probe that leaves lineages out of its statistics (16S:
  chloroplast, mitochondrial and eukaryote sequences, the pack's `probe_exclude_taxa`)
  came from those lineages, the probe stopped with "no reference sequences had both
  primers", which sends the user to check the primers' orientation. The message now gives
  the number of products and the lineages they came from. The other two texts of this
  stop are unchanged.
- **Tested.** With Snakemake, the 16S test set with a manual truncation (243/154, the cuts
  `auto` chooses on it), `expected_length: 385` and the template's length filter stopped
  with the `KeyError` before the change and now finishes, with the window [327, 443]. The
  probe, run with its cutadapt and seqtk on a synthetic reference of chloroplast records,
  gave the old text before the change and the new one after; on half bacterial records it
  measures as before.

### The amplicon tables: a stable format, and a named first column

- **The format of the four tables in `06.taxonomy/` is kept from this release on**, for
  the scripts and tools that read them, and a change to it will be announced here: R's
  `write.table` conventions (tab-separated, text in double quotes, an unquoted `NA` for a
  missing value) with the first column named in the header (below), the `taxonomy`
  column right after the sample columns, the rank prefixes each marker pack lists, the
  whole species name in the `s__` slot (genus included), and no confidence values. The
  ASV IDs (`ASV_1`, `ASV_2`, ...) are not part of it: they number the ASVs of one run by
  read count, and the same sequence can get another number in another run or version;
  match ASVs across runs by their sequence, the first column of `asv_table_seqs.txt`. See
  "Output files" in the manual.
- **Every table names its first column**: `asv_id` in `asv_table.txt`,
  `taxon_seq_table.txt` and `asv_table_unfiltered.txt`, `sequence` in
  `asv_table_seqs.txt`, under both taxonomy methods. The header started with an empty
  cell, as R's `write.table(col.names = NA)` writes it above a column of row names. In
  `taxon_seq_table.txt`, whose other cells are all text, readers that guess whether the
  first line is a header (`data.table::fread`, for one) took that line for a data row.
  The rows are unchanged. `read.table(..., header = TRUE, row.names = 1)` in R and
  `pd.read_csv(..., sep = "\t", index_col = 0)` in pandas read the tables as before;
  read as an ordinary column, the first column is now called `asv_id` or `sequence`
  instead of `X` (R's `read.table`), `V1` (`fread`) or `Unnamed: 0` (pandas). The two
  tables COI runs add to `06.taxonomy/`, `bacterial_cox1_asvs.tsv` and
  `coi_asv_flags.tsv`, name theirs `asv_id` too (`coi_screens.tsv` already did).
- **Tested with Snakemake**: `assign_taxonomy` run again with the changed code in copies of
  finished runs of the 16S test set (`method: sintax`, 199 ASVs) and of the gyrB test set
  (`method: rdp`, 237 ASVs) wrote the four tables row for row as before, with only the
  first header cell changed; so did a run of the COI test library amplified with BF3/BR2
  (28 ASVs moved by the off-target guard), for all six of its tables.

### 16S: ASVs that Metaxa2 trims to one sequence are merged

- **The 16S extraction (`target_extract`, Metaxa2) merges ASVs that come out of it with
  the same sequence**, as the ITS extraction always did after ITSx: the lowest-numbered
  ASV of the group (the most abundant) stands for it, with the group's per-sample counts
  added up, and keeps its own call in `metaxa2_origin.tsv`. Metaxa2 seldom trims a 16S
  ASV, but it can: on a V4 run of 144 field samples it cut 31 bases off the 3′ end of five
  ASVs. Two ASVs that differed only in the bases it cuts away would reach the taxonomy
  with one sequence, and then `method: rdp` stopped at `asv_table_seqs.txt`, on R's
  "duplicate 'row.names' are not allowed", and `method: sintax` wrote the sequence on two
  rows of that table. No such pair turned up in the runs checked: in the 16S test runs
  of this release Metaxa2 trimmed no ASV, and the five ASVs of the V4 run kept sequences
  of their own.
- **A new file, `05.dada2/metaxa2_collapse_map.tsv`**, in the layout of
  `itsx_collapse_map.tsv`: one row per ASV left after merging, with the IDs merged into
  it, their number and the length of the extracted sequence. Where nothing is merged, the
  FASTA and the count table of the extraction are written as before, byte for byte.
- **The taxonomy stops when one sequence labels two ASVs**, under both methods, with a
  message that names them, instead of writing such a table. After the merge this should
  not happen.
- **Tested with Snakemake**, in copies of a finished run of the 16S test set (231 ASVs
  extracted by Metaxa2):
  - `target_extract` run again with the changed code wrote the extracted FASTA, its count
    table and `metaxa2_origin.tsv` byte for byte as before, and a collapse map of 231
    groups of one;
  - in a copy whose `ASV_7` had been given the sequence of `ASV_5`, so that Metaxa2
    returned the two identical, the step merged `ASV_7` into `ASV_5`, adding up their
    counts, and the taxonomy wrote one row per sequence under both methods, with
    `stats/read_tracking.txt` as before up to the taxon filter (reads move, none are
    lost). With the code before the change, the same copy made `method: rdp` stop on
    "duplicate 'row.names' are not allowed" and `method: sintax` write one sequence on
    two rows;
  - run on the output of the code before the change (two ASVs, one sequence), the new
    taxonomy step stopped with its message under both methods, and none of the four
    tables was left.

### The SRA run info of each sample (shotgun, BioProject mode)

- **A new file, `sra_runinfo.tsv`, next to `samples.tsv`.** In BioProject mode MetaFlux
  has always fetched NCBI's run info of the BioProject at start-up, to pick the paired
  runs, and kept it as `{bioproject}_runinfo.csv` in `fastq_dir`, but nothing of it
  reached the output folder. The new step `sra_runinfo` writes the rows of the run's
  samples there: the column `sample`, then every column of the run info, named and filled
  as NCBI writes them (run, BioSample, sample name, library, instrument, centre and more),
  one row per sample in the order of `samples.tsv`. It is the start of a sample metadata
  table. It holds no BioSample attributes (habitat, host, place, collection date): NCBI's
  run info does not carry them. Runs on local reads are unchanged.
- **Tested with Snakemake**: a BioProject-mode config for three runs of PRJNA562911 (their
  run info already cached in `fastq_dir`, so nothing was fetched) planned the step, which
  wrote one row per run, with the 47 columns of NCBI's run info after `sample`.

### Documentation

After 2.6.1, the online manual was updated to call the set bacterial COI homologues: the
page "Bacterial COX1 reference set (COI)" became "Bacterial COI homologues (reference
set)", at the same address, and the other pages use the same name. A note on that page
says what the name covers and that file names, settings and output columns keep the
short form `cox1`. The taxon filter page now quotes the stop of the COI screens in its
new wording.

### Checks at start-up

- `max_consist` must be a whole number of at least 1 (also written `200.0` or `"200"`).
  A fraction, which DADA2 cut down without a word, `Inf` or text, which failed at the end
  of the first round, and `true`/`false` now stop the run with one message. A note at
  start-up when it is 999, the value of every template up to v2.6.1.
- `max_flat_rounds` must be a whole number of at least 0; a note when it is not smaller
  than `max_consist`, because the early stop then can never act.
- `error_model` must be `auto`, `loess`, `binned` or the experimental `custom`, in any
  capitals.
- `trim_final_cycle` must be `auto`, `true` or `false`, in any capitals.
- `custom_error_function` (experimental): with `error_model: custom`, `file` and `name`
  must both be set; with any other value, both or neither. The file must exist and be
  readable, and the name must be one R can give a function.
- `accept_unconverged` must hold model codes such as `"R2:0a55be94"`, as a list or a
  single text; `true` or `false` stops the run with a message that says what it takes.
- A `learn_errors:` or `dada2:` block left empty in a run config (all its lines commented
  out, as after commenting out `accept_unconverged`) stops the run with a message that
  names the block. YAML reads such a block as empty, and it replaced the template's whole
  block, so MetaFlux stopped with a bare Python error.
- A key written twice in the same block of a file given with `--configfile` stops the
  run, with the file, the key and both lines. YAML keeps only the last of the two without
  a word, so a second `amplicon:` block pasted at the end of a run config, for example to
  add a line a message asks for, replaced the first one whole, and the marker, the primers
  and the other amplicon settings of the first fell back to the template.
- `dada_filter` stops with a message that names the filtering settings when no read pair
  of a sample passes the filter. It used to end with Snakemake's bare "missing output
  files" error.
- Every sample needs its R2 (or `_2`) file, with the extension of its R1 file ("Raw reads
  named `.fq` or `.fq.gz`, and the R2 file of every sample" above).
- `mode`, `amplicon.type`, `amplicon.its_region`, `amplicon.primers.orientation` and
  `amplicon.taxonomy.method` left empty, or given a number or `true`/`false`, stop the run
  with their usual message, which shows the value. They ended with a bare Python
  `AttributeError`, and an empty `its_region` did so for every marker, although only ITS
  uses it; for the others it is now ignored.
- In amplicon mode, sample names that look like a number, a missing value or true/false
  (`01`, `002`, `1e5`, `NA`, `true`) stop the run, with the same test as in shotgun runs
  with sylph. MetaFlux's rdp taxonomy step reads the sequence table with
  `data.table::fread`, which reads such names back as numbers: `002` came out as sample 2,
  and `01` and `1` became one name. R and pandas do the same when they read the tables
  afterwards. Any name with a letter, such as an SRA run accession, passes.

### What runs again in a folder made with v2.6.1

The changed wording of the COI steps ("The messages of the COI steps say bacterial COI
homologues" above) changes comments and messages only; the steps it touches (the read
tracking, and the taxonomy under `method: sintax`) run again with the other changes below,
with the same results.

- **After a `git pull`** in the MetaFlux folder that made the output folder, in every
  amplicon folder: `dada_filter` for every sample (it has a new output, the quality values
  of the filtered reads), Falco on the filtered reads, `choose_error_model`,
  `learn_errors`, `error_model_report` and `check_error_model` (new), and `dada_seqtab`
  with every step after it: the extraction where it runs, the length filter, the COI
  screens, the taxonomy, the read tracking and the phylogeny. MultiQC runs again in both
  modes, because its settings changed. Also in every amplicon folder: `pick_trunclen` (its
  script changed: the stop text of the overlap check, the new log lines and the length
  floor; it also reads the primer-trimmed reads now) and `dada_quality_plots` (a comment
  in its script changed), and, where `expected_length: auto`, the length probe, once per
  marker, reference and primer pair (a comment in its script changed; its result is cached
  in `refdb/cache/` and shared by the folders that use it). Outside 16S the probe gives
  the same result, and `pick_trunclen` the same cuts on the same primer-trimmed reads.
  **16S folders with
  `expected_length: auto`** compute the probe anew under its new name
  (`..._excl-second-site_...`), once per primer pair and SILVA release, and `pick_trunclen`
  and every step after it run on the new statistics; where p95 or q1 changed, the cuts and
  the ASV length window can change, and with them the tables (see "Results" below). Also
  in every amplicon folder: `trim_primers` for every sample, because its new program
  (`30a_final_cycle.py`, the decision on the last cycle) is an input newer than the
  trimmed reads, and with it Falco and the read count of the primer-trimmed reads. The
  reads it writes change only where the last cycle is removed, under `auto` on ITS. Also
  the raw read count (`count_reads_raw`), whose command now also reads `.fq` and `.fq.gz`
  names, with the same counts. Kept: the read links, the PhiX removal, Falco of the raw and
  PhiX-filtered reads and the reference files. In shotgun mode only MultiQC runs again,
  and, in BioProject mode, the new step `sra_runinfo` writes its file. (The named first
  header cells come with the taxonomy, and the merge after Metaxa2 with the extraction,
  which both run again in every amplicon folder anyway; in 16S folders the extraction
  also writes its new `metaxa2_collapse_map.tsv`.)
- **From a fresh clone, a release archive or a second worktree**: every step that runs a
  script, with the steps after it, as at every release ("What runs again in a folder made
  with v2.5.0", under 2.6.0, says why), and `trim_primers` with every step after it, for
  the same reason as above: its new input is newer than the trimmed reads.
- **`--rerun-triggers`.** The record, the report and the verdict of the new steps are
  files a run asks for, and they are missing in a folder made with v2.6.1, so learning,
  the report, the check and `dada_seqtab` with every step after it run whatever
  `--rerun-triggers` says. So does `dada_filter`: the choice of the error function reads
  the quality values it now writes to a file of its own per sample, which such a folder
  lacks, also under a `--rerun-triggers` set without `code`, such as the
  `mtime params input software-env` that earlier release notes suggested for skipping
  reruns. The 16S probe is a missing file too (it has a new name), so it is computed, and
  `pick_trunclen` runs again on it, whatever `--rerun-triggers` says; `pick_trunclen` also
  has new parameters (the key of `amplicon.probe_length_stat` its stop text names, the
  method and the filter settings it replays) and new inputs (the primer-trimmed reads).
  `trim_primers` runs again under any set that includes `mtime`, because its new input
  is newer than the trimmed reads, and it has a new parameter (`trim_final_cycle`).
  Skipping the reruns saves nothing in this release.
- **Results.** Outside ITS, under `auto`, the primer-trimmed and the filtered reads are
  those of v2.6.1 (see the end of this item). Where neither the cap nor the early
  stop acts and the loess fit is used, the models and the tables are those of v2.6.1:
  full runs with Snakemake of the eight MiSeq test sets (16S, ITS, 18S, gyrB, rpoB and
  three COI libraries) gave, byte for byte, the tables of v2.6.1: in `05.dada2/`,
  `seqs.fasta`, `seqtab_head_seqs.txt`, `seqtab_head_names.txt`, `read.counts` and the
  tables of the extraction and of the length filter; the taxonomy tables of
  `06.taxonomy/` (all but the raw SINTAX output); `stats/read_tracking.txt`; and the
  alignment and the tree of the 16S run. The raw working files
  of Metaxa2 (`05.dada2/metaxa2/seqs.*`), ITSx (`05.dada2/_itsx_tmp/`, and the copy of
  ITSx's summary, `05.dada2/itsx_extraction.summary.txt`) and SINTAX
  (`06.taxonomy/_sintax_raw.tsv`) can differ between two identical runs of any version
  (start times, Metaxa2's own classification, the order in which ITSx and VSEARCH write
  their results); no table depends on what varies in them. A model that the cap or the
  early stop ends now stops the run for acceptance instead of being used. Binned data
  under `auto` get the binned fit, so different tables, on purpose. These comparisons
  were made before the removal of the last cycle and the new ITS window ("Primer
  trimming: the extra last sequencing cycle is removed from untruncated reads" and "The
  ASV length window for ITS" above), which change the ITS tables on purpose. Repeated on
  the final code: 18S, gyrB, rpoB and the three COI libraries gave the tables of v2.6.1,
  byte for byte, and 16S differs only by the probe change below. The copies of the rpoB
  and gyrB sets re-binned in silico, the Zymo D6331 mock and the NIST mix A gave, byte for
  byte, the tables of the same code without the two changes.
- **Results of the 16S probe change.** The comparison above was made before it. Run again
  with Snakemake, the 16S test set (799F/1175R, SILVA 144) gets p95 385 bp instead of
  386, so `raise_trunc` stops one base earlier: `truncLen` 243/154 instead of 244/154
  (overlap slack 0 in both). The filter keeps 0.55 % more read pairs (218,313 instead of
  217,113), and the run ends with
  237 ASVs after merging instead of 234, 199 in the final table instead of 197, and
  194,158 final reads instead of 193,097 (+0.55 %); the error models get new codes, and
  every table from `05.dada2/` on, the taxonomy and the tree differ accordingly. The ASV
  length window is [324, 435] bp instead of [324, 436]. The Zymo D6331 library
  (341F/805R, NextSeq 2000) now passes `pick_trunclen` with the template's `p95` (429 bp):
  with the same cuts as with `q3` (303/259), its error models, model codes, ASVs and read
  counts are those of the run with `q3`; its ASV length window is [353, 479] bp instead of
  [353, 632].
- **Results of the ITS window change.** The comparison above was also made before it.
  Replayed on the inputs of the length filter of the v2.6.1 ITS run, the window is
  [38, 343] bp instead of [125, 343], and the filter keeps all 16 ASVs that ITSx
  extracted instead of 12, with 205,569 reads instead of 74,077. The four ASVs it adds
  are the short-ITS2 yeasts (*Candida haemulonii*, *Metschnikowia pulcherrima* and a
  *Dipodascus*); the tables of the length filter, the taxonomy tables and the read
  tracking of the ITS run differ accordingly. Replayed the same way, the seven other
  test sets give the outputs of the length filter of v2.6.1 byte for byte. [To be
  confirmed before release: a full ITS run on the final code, and the ITS results of
  "Test datasets" and of the ITS page brought in line with it.]
- **Results of the changes for `.fq` names and long output paths.** Full runs with
  Snakemake of the eight MiSeq test sets and the Zymo D6331 mock above, and of two shotgun
  libraries (Zymo D6331 and PRJNA1477390), gave, byte for byte, the result tables of the
  same runs before these two changes.
- **Results of the named first column and of the merge after Metaxa2.** The comparisons
  above were made before both. The four tables of `06.taxonomy/` now differ from those of
  v2.6.1 in their first header cell only ("The amplicon tables" above). The merge changes
  nothing where no two ASVs come out of Metaxa2 identical, which was the case in every
  16S run compared above.

### What to edit in a config written for v2.6.1

- `max_consist: 999`, in a run config copied from an older template: change it to 200 or
  delete the line. With 999 a model that keeps changing runs for days before the check,
  unless the early stop ends it first; MetaFlux prints a note.
- A lower `max_consist`, such as the 10 rounds of QIIME 2: the run now stops where it
  used to go on with a model that had not converged. Raise it, or delete the line.
- A cluster profile that gives `dada_seqtab` its own resources (a longer run time, for
  example) needs the same for `learn_errors`, which now does the long part of that work.
  In the run config, `learn_errors` takes the `dada_seqtab` values of `resources` unless
  it is named there. `error_model_report` is a short job that uses one thread: unless it
  is named there, it takes 1 thread and the memory of `dada_quality_plots` (4000 MB in
  the template).
- To keep the fit of v2.6.1 on binned data: `error_model: loess`.
- **New stops that data finished with v2.6.1 can meet.** Under `auto` only: samples from
  two instruments or two sets of binned values, and 11 to 19 quality values or 3 to 10
  that include neighbours (MiSeq reads under a `trunc_q` that leaves fewer than 20 of
  their values);
  `error_model: loess` turns these into warnings and runs as v2.6.1 did. Under any
  setting: a sample with rewritten quality values, reads with fewer than three quality
  values above Q2, and a sample of which the filter keeps no read.
- An amplicon run whose sample names look like numbers (`01_R1.fastq.gz`) now stops at
  start-up: rename the read files (`S01_R1.fastq.gz`, for example) and run them into a new
  output folder.
- A script of your own that reads the tables of `06.taxonomy/` and finds their first
  column by the stand-in name its reader gave the empty header cell (`X` in R's
  `read.table`, `V1` in `data.table::fread`, `Unnamed: 0` in pandas): the column is now
  called `asv_id` (`sequence` in `asv_table_seqs.txt`). Tables read with `row.names = 1`
  in R or `index_col = 0` in pandas come out as before.
- A 16S marker pack of your own in `config/markers/`: add `probe_exclude_second_site: true`
  at its top level (see `workflow/markers/16S.yaml`). Without it the probe keeps the
  statistics of v2.6.1, and 341F/805R runs still stop at `pick_trunclen`.
- A `probe_length_stat` of `q3`, or a fixed `expected_length`, set only to get a 341F/805R
  run past `pick_trunclen`, is no longer needed: the default `p95` now gives 429 bp on
  SILVA 144.
- A run that stopped at the overlap check because of a short off-target product and was
  then finished with `trunc_len.mode: manual` can go back to `auto`.
- Nothing to add for the new key `cutadapt.trim_final_cycle`: the template supplies
  `auto`. To keep the primer-trimmed ITS reads of v2.6.1, set `trim_final_cycle: false`
  under `amplicon.cutadapt` (primer trimming still runs again once, with the same result).
- `trim_final_cycle: true` with `trunc_len.mode: manual` and a cut at the full length of
  the primer-trimmed reads: the reads are one base shorter, and `filterAndTrim` discards
  every read shorter than its cut. Lower the cut by one base.
- Nothing to add for the new keys `trunc_len.method` and `simulate_extra_overlap`: the
  template supplies them.
- A setting MetaFlux asks for in a message goes inside the `amplicon:` block the run
  config already has; a second `amplicon:` block would replace the first, and the other
  amplicon settings would be lost. MetaFlux now stops at start-up when it finds one (see
  "Checks at start-up").

## 2.6.1 — 2026-11-02

### Settings typed inline with `--config` are read as in a config file

A setting inside a block can be typed on the command line as inline YAML, for example
`--config 'amplicon={phylogeny: {enabled: false}}'`. Snakemake hands every value of such
a block over as text, and two kinds of setting then went wrong without a word:

- **On/off settings** were switched on by `false`, because Python counts any text that
  is not empty as true: the phylogeny and its support values, the taxonomy filters of
  both modes, PhiX removal in both modes, the extraction, the two COI screens, the k-mer
  evidence gate and Kraken2's memory mapping.
- **`amplicon.dada2.max_ee`** reached DADA2 as text, and DADA2 compared each read's
  expected errors with it as text, so a read with 10.5 expected errors passed
  `max_ee: [2, 2]` ("10.5" sorts before "2").

Settings written in a config file were never affected. MetaFlux now reads every value
given with `--config` the way YAML reads it in a config file, before any setting is used.
The "Running MetaFlux" page also said that an inline block replaces the whole block of
the config; Snakemake merges it key by key, as it merges a run config, and the page now
says so.

### The MIDORI2 conversion writes its log

`prepare_midori2_coi_sintax` declared a log file but wrote nothing to it. It now writes
there the files it read and wrote and the main counts, or the message it stopped with.
The converted reference is byte for byte the same.

### Documentation

The pages on COI and on the bacterial COX1 set, the configuration reference and the COI
marker pack said that the MIDORI2 records in the set's list are at least 90 % identical
to a gene of the set. They are at least 90 % identical to a gene that passes checks
5.1–5.3 of its build, and 17 of the 186 match only genes that check 5.4 then left out.
The set's README reports the sensitivity of its counts for twelve thresholds, not for
each of them, and the COI screens search far deeper than VSEARCH's default, not eight
times deeper.

On the home page of the documentation, the metro map now opens in an image viewer on the
page instead of as a bare PNG: a second click shows it at full size, and x or Esc closes
it. The viewer is the mkdocs-glightbox plugin, a new requirement of the documentation
build only (`pip install mkdocs-material==9.7.7 mkdocs-glightbox==0.5.3`). The map and
the README overview now say v2.6.1 and are otherwise unchanged; the map PNGs are 3509 x
1291 px instead of 5848 x 2152, so that at full size its station labels are about 20 px
high.

### What runs again in a folder made with v2.6.0

- **After a `git pull`** in the MetaFlux folder that made the output folder: MultiQC in
  every folder, because its settings file, one of its inputs, now names v2.6.1. In a COI
  folder also the MIDORI2 conversion, whose script changed, and every step that follows
  from the converted reference: the length probe (`expected_length: auto`, the default),
  `pick_trunclen`, DADA2, the COI screens, the taxonomy and the read tracking. The steps
  up to primer trimming are kept. The change to `--config` is read at start-up and makes
  no step run again.
- **From a fresh clone, a release archive or a second worktree**: every step that runs a
  script, with the steps after it, as at every release ("What runs again in a folder made
  with v2.5.0" below says why).
- The steps that run again give the results of v2.6.0, as between any two identical runs:
  the converted MIDORI2 file is byte for byte the same.
- **To skip the reruns**, give `--rerun-triggers mtime params input software-env` to that
  run and to every later run on the same output folder; only MultiQC then runs again.
- A run made with one of the settings above typed inside a block with `--config` gave a
  wrong result: run it again with the same command. A corrected `max_ee` is a changed
  step parameter, which Snakemake notices; a switch typed `false` now leaves its step
  out, and the files that step wrote before stay in the output folder.

### What to edit in a config written for v2.6.0

Nothing.

## 2.6.0 — 2026-10-31

### COI, a new amplicon marker

`amplicon.type: COI` runs the mitochondrial COI gene, the animal barcode, through the
amplicon path and classifies it with SINTAX against MIDORI2 (release GB273, UNIQ NUC_SP
variant). `method: rdp` is refused for COI: a DADA2 classifier trained on MIDORI2 would
need about 225 GB of memory at species level.

- **The reference** is downloaded from MIDORI2's Zenodo record (zip MD5 and FASTA SHA-256
  checked) and converted once: the domain moves to `d:`, the kingdom (Metazoa, Fungi,
  Viridiplantae) is added from MIDORI2's own NCBI taxonomy snapshot, the taxid suffixes
  go, placeholder ranks and interim species names ("sp.", BOLD codes) are left out of the
  headers, and the MIDORI2 records that are bacterial COX1 filed under eukaryote names
  are removed. `keep: [k__Metazoa]` keeps the animals of a run.
- **Two COI screens** run before the taxon filter, both on by default (`amplicon.coi`):
  a reading-frame check removes ASVs with a TAA or TAG in the run's frame (allowing
  for a primer cut a base off at the ASV's start), or a length too far from the run's
  reference length (pseudogenes, nuclear copies of the gene called NUMTs, artefacts);
  an off-target guard moves to `06.taxonomy/bacterial_cox1_asvs.tsv` the ASVs that are
  closer to a set of bacterial COX1 genes (B) than to MIDORI2 (E), *Wolbachia* above
  all. `06.taxonomy/coi_screens.tsv` holds the measurements behind both;
  `asv_table_unfiltered.txt` marks the ASVs `removed_frame` or `moved_bacterial_cox1`;
  the read-tracking table has the columns `post_frame_check` and `post_bacterial_cox1`.
- **Species flags.** `06.taxonomy/coi_asv_flags.tsv` says, per kept ASV, how far a
  species call can be trusted (`low_identity` below 98 %, `shared_sequence` when a
  record of another named species is as close) and flags ASVs less than 80 % identical
  to every MIDORI2 record (`far_from_midori2`).
- **The bacterial COX1 set** (`bacterial_cox1_gtdb-r232.v1`, Zenodo
  10.5281/zenodo.23820760) was built for this release from the 199,923
  species representatives of GTDB R232; how, and how it was tested, is on the new
  manual page "Bacterial COX1 reference set (COI)".

Checked on the insect mock of Elbrecht et al. (2019), one library per primer pair run as
one sample: 361 (BF3/BR2), 355 (Leray-XT) and 362 (fwhF2/fwhR2n) of its 374 specimens
are recovered by an ASV at least 97 % identical to their barcode, the ASVs that match a
specimen hold 99.1–99.3 % of the reads of the final tables, the off-target guard caught
every *Wolbachia* ASV that SINTAX had called a stink bug absent from the mock (41 ASVs,
4,315 reads), and no specimen lost its only ASV to either screen.

### Mixed orientation: reversed pairs swapped back, and kept once

With `amplicon.primers.orientation: mixed`, `trim_primers` runs a second Cutadapt pass
with the primers swapped, for the pairs whose R1 carries the reverse primer. Three things
change under `mixed`. With `orientation: fixed`, the default, Cutadapt gets the same
arguments and writes the same reads as in v2.5.0, and no new file is written.

- **The two mates of a reversed pair are swapped, not reverse-complemented.** Every
  version since v1.0.0 reverse-complemented each mate. That put the reads on the right
  strand but left R1 on the 3′ end of the amplicon and R2 on the 5′ end, so DADA2's
  `mergePairs` (with `trim_overhang: true`, the default) kept only the stretch where the
  two overlap. Now the pair's R2 goes to the R1 file and its R1 to the R2 file, and the
  pair reaches DADA2 like the pairs of the first pass. On the BF3/BR2 COI insect mock
  (SRR8082166, in which every R1 carries the reverse primer), none of the merged reads
  came within 9 bp of the full 418 bp with the old handling, against 99.4 % with the
  mates swapped, in the DADA2 test that led to the fix; the v2.6.0 run of the same
  library gives 99.5 %.
- **Both 5′ calls need 10 primer bases** (`--overlap 10`, under `mixed` only). Both
  passes read every pair, and Cutadapt counts a 5′ primer as found when a read starts
  with the primer's last three bases. Where the start of each primer resembles the end
  of the other, a pair passed both and was written twice, once with its mates the wrong
  way round and most of a primer still on: 13–17 % of the pairs in the two BF3/BR2
  libraries of the COI mock and 0.14 % in its fwhF2/fwhR2n library (none in its
  Leray-XT library), 0.3 % in the rpoB and 0.04 % in the gyrB test set, 2 pairs in the
  18S test set, none in the 16S and ITS test sets. With 10 bases, 1 pair in
  SRR8082166 and 32 in the gyrB test set are still written twice. The cost is the pairs
  whose primer is found from only 3 to 9 of its bases at the start of a mate, which are
  now dropped: 0.006–0.15 % of the pairs in the 16S, ITS, 18S, gyrB and rpoB test
  sets, and at most 0.06 % in the four COI mock libraries. In those pairs the mate cut
  short is hardly ever a cleanly trimmed amplicon read.
- **A new file, `stats/orientation/{sample}.tsv`**, under `mixed` only: the read pairs
  that went into primer trimming, the pairs the first pass kept as sequenced and the
  pairs the second pass kept and swapped, each with its share of the pairs that went
  in, from Cutadapt's reports of the two 5′ calls. Nothing in the pipeline reads it. A
  folder trimmed with `mixed` by v2.5.0 lacks it, and a run asks for it, so primer
  trimming and every step after it run again in such a folder, however v2.6.0 was
  obtained: that is what brings the two changes above to old folders.

For a user of `mixed`: reversed pairs now merge to the full amplicon and are hardly
ever written twice, so the tables of a `mixed` run change; re-run such folders, and
MetaFlux processes them again from primer trimming on (see below). The Cutadapt
section of the MultiQC report shows the first pass only. A library in which every pair
is swapped is better run with its two input files swapped and `fixed`: Cutadapt then
reads every pair once, and its report covers all the reads.

In a truly mixed library the R1 file now holds read 1 of the pairs that came the right
way round and read 2 of the reversed ones, and DADA2 fits its R1 and R2 error models to
that blend. On a genuinely mixed library (two BF3/BR2 libraries of the COI insect mock
in opposite orientations, 1.6 million pairs as one sample) the blend recovered 363 of
the 374 specimens against 364 for the two libraries run on their own, with 600 kept
ASVs shared of 654 and 664 and tables almost identical (Bray-Curtis dissimilarity 0.011
on specimen proportions, on a 0-to-1 scale); it kept 7 % fewer reads, lost at DADA2's
quality filter because the truncation lengths are chosen on the blended quality
profile, and that is an upper bound, since the two libraries also come from two
sequencing runs. With the old handling the same library gave no result: DADA2's error
model had not converged after 171 rounds and 10.5 hours, on reads of which half were
reverse-complemented and 7 % still carried most of a primer (see the troubleshooting
page).

### Amplicon runs with a single sample

A run with one sample no longer stops at `dada_seqtab`. Given one sample, DADA2's
`derepFastq()`, `dada()` and `mergePairs()` return that sample's result on its own
instead of a list with one result per sample, and the read tracking of `dada_seqtab`
stopped on it (had it got further, the sequence tables would have named the sample
"1"). `dada_seqtab` now puts a lone result back into a one-element list named after the
sample. With one sample, `pool` has nothing to pool: DADA2 denoises it on its own.

- **Two or more samples: nothing changes.** On the 16S test set the sequence tables,
  the FASTA, the read counts and the error-plot PNGs came out byte for byte the same.
- **Release check.** `CONTRIBUTING.md` adds a real run with a single sample to the
  release checks, since the unit tests and the dry runs never run R.
- **Checked on five one-sample runs** (16S with SINTAX and with `rdp` plus the species
  step, ITS, gyrB under `rdp`, and the shotgun test set with `profiler: both`): all
  finished. The shotgun sample came out identical to its column in the four-sample
  run; the amplicon tables differ from the same sample's column in a multi-sample run
  only through what a run decides from all of its samples (the truncation lengths,
  DADA2's error model and its chimera vote, and the classifier's random draws, which
  depend on which ASVs the run has), as expected.

### k-mer evidence floor `auto`, sized to each sample's depth

`shotgun.kmer_evidence.min_distinct_minimizers: auto` gives each sample its own floor,
worked out after Kraken2 from the sample's own report: (classified read pairs −
*Homo sapiens* read pairs) / 1,000, at most 6,500. A species stays when its distinct
minimizers reach the unrounded floor; `min_reads` and the 1 % safety rail work as
before. The human read pairs are left out because PlusPF holds the human genome and
host removal is off by default (the row is recognised by taxid 9606 and its name, so
nothing is left out on a GTDB-based index, where 9606 is another taxon). The cap is
there because MetaFlux's benchmark stops at 6.5 million classified pairs.
`logs/kmer_evidence/{sample}.log` gives the classified pairs, the human pairs left out
and the floor; the reasons of the evidence table and the plot give the floor rounded
up, marked `(auto)`. The evidence table keeps its columns.

- **Opt-in; 300 stays the default.** On MetaFlux's benchmark `auto` raised species F1
  more than 300 did and lowered no sample's F1, and on one CAMI II sample run at four
  real depths (0.76 to 6.31 million classified pairs) it removed no true species; but
  three quarters of the benchmark gain came from one real mock community, and it
  changes the table far more than 300 does. On the four test-dataset samples it keeps
  25–61 species against 66–104, and 2.7–15.1 % of the read pairs leave the table,
  against at most 0.09 % at 300, most of them from one order that is probably really
  there (the 16S rRNA fragments in the same reads point that way) but has no close
  reference genome. That read loss is why `auto` stays opt-in. It also assumes paired
  2×150 bp reads. See "Scaling the floor with depth" on the k-mer evidence page before
  switching to it.
- **With a whole number nothing changes.** The rule gets the same number as before. On
  the four test-dataset reports and the MBARC-26 report, with the gate off at floors 0
  and 300 and on at 300 and 1,000, the evidence table, the gated report, the plot and
  the log are byte for byte those of v2.5.0.
- **The setting is checked more strictly at start-up.** A whole number (also written
  `300.0` or `"300"`) or `auto`, in any mix of capitals and with spaces around it, is
  accepted. A fraction (v2.5.0 cut `300.5` to 300 without a word) and `true`/`false`
  (read as 1 and 0) now stop the run, with one message that says what is allowed;
  other text and an empty value, which stopped v2.5.0 with a Python traceback, get the
  same message. A negative number stops the run as before. A config with such a value
  needs a whole number or `auto`.
- **A correction to the v2.5.0 k-mer evidence page.** Its section on scaling the floor
  with depth said that the read pairs of the species a depth-scaled floor removes move
  to relatives (10–33 % of the read pairs at a floor of 1,673 on the test dataset, and
  "13 to 42 % end up under a different species" at 3,000). Measured with `auto` on the
  same four samples (floors of 1,673 to 3,458), 2.7–15.1 % of the read pairs leave the
  table instead, because no species of their genus, family or order is left to take
  them; at 300 the loss is at most 0.09 %. The pipeline's behaviour does not change. A
  table made with a floor sized by hand after that page should be compared with the
  same run at the gate off before it is relied on.

### The Zymo mock in the documentation

The documentation now quotes the ZymoBIOMICS D6331 mock from a published run, SRA run
SRR29240428 (Martí et al. 2025), subsampled to 5.30 million read pairs after fastp, the
depth of the run used before. The earlier run, SRR38976715, is not linked to a
publication and is no longer used. Every D6331 figure was measured again on the new run
with the same method. Two change visibly: on the mock, the best Bracken alpha is now 5.7
to 21 times the default (was 4.5 to 7 times), and the abundance error against Zymo's
composition (L1) is about 19.7 percentage points (was 40.9). The pipeline does not
change, and nothing runs again.

### What runs again in a folder made with v2.5.0

What runs again depends on how you get v2.6.0, because Snakemake goes by two kinds of
evidence. It keeps a record of every finished step in `.snakemake/`, inside the MetaFlux
folder it is started from: the step's command, its parameters, the list of its input
files and its software environment. A step whose record no longer matches runs again. A
step that runs a script is judged by the script file's date instead of its text: it
runs again whenever the script is newer than the step's outputs, whatever the script
says. Every step after a step that runs again runs again too.

**Updating the MetaFlux folder that made the output folder (`git pull`).** Git gives a
new date only to the files it changes, and Snakemake has the records.

- Amplicon mode, every marker, both orientation settings: `trim_primers` and every step
  after it run again, because its command changed (with `mixed`, also because its new
  file, `stats/orientation/{sample}.tsv`, is missing): falco on the trimmed reads, the
  stripped-read counts, `pick_trunclen`, DADA2, the extraction where it runs, the
  taxonomy, the read tracking, the phylogeny and MultiQC. The raw-read links, counts
  and QC, the PhiX removal, the primer files, the reference databases and the length
  probe are kept. (The scripts of `dada_seqtab`, `assign_taxonomy` under `sintax` and
  `aggregate_read_counts` changed too; those steps come after `trim_primers` anyway.)
- Shotgun mode with Kraken2 (`profiler: kraken2` or `both`): `kmer_evidence` runs again
  for every sample, because its script changed, and with the gate on
  (`kmer_evidence.enabled: true`) so do `bracken` and `finalize_otu_table`. Kraken2 does
  not run again.
- Every folder, in both modes and with every profiler: MultiQC runs again, because the
  release changes the version its report shows, in `workflow/multiqc_config.yaml`,
  which is one of MultiQC's inputs.

**Starting from a fresh clone, a release archive or a second worktree.** A worktree is
another checkout of the same repository, made with `git worktree`. Every file of such a
copy is newer than the old outputs (a clone or a worktree dates every file when it is
made, a release archive at the time of the release commit), and it has no records of
the old folder. So every step that runs a script runs again, with everything after it,
whether v2.6.0 changed that script or not: in amplicon mode the quality plots of the
trimmed reads, `pick_trunclen` and everything after it (fetching the reference databases
and computing the length probe again unless `refdb/` was carried over), in shotgun mode
the k-mer evidence (and, with the gate on, Bracken), the sylph tables, the OTU table,
the read tracking and MultiQC. A step without a script has nothing to be compared with,
and a changed setting goes unnoticed, so in a folder made with `fixed` `trim_primers`
does not run again. In a folder made with `mixed` it does, for every sample and with
every step after it, because the orientation table it now writes is missing there.

**Results.** With `orientation: fixed` the steps that run again give the results of
v2.5.0, except that under `method: rdp` a few bootstrap values, and now and then a
lineage, can move, as between any two identical runs (DADA2's `assignTaxonomy` breaks
exact ties at random; see the troubleshooting page). Each change was checked on its
own:

- `fixed`: `trim_primers` writes the same reads as v2.5.0 (30 samples: the five test
  sets and four COI mock libraries).
- `mixed`, on the 16S test set with every pair swapped: the taxonomy, the ASV names and
  the tree are those of the v2.5.0 run of the original files, with 7 of 193,097 final
  reads fewer (the pairs the 10 bases cost).
- Two or more samples: `dada_seqtab` writes the same tables (16S test set).
- The COI code: the 18S run of 120 samples made with v2.5.0 gave the four tables of
  `06.taxonomy/` and `stats/read_tracking.txt` byte for byte the same.
- The default floor of 300: the k-mer evidence files are byte for byte those of v2.5.0.

All the changes together were checked too: the seven test-dataset runs made for v2.5.0
(the 16S, ITS, 18S, gyrB and rpoB test sets with `fixed`, gyrB and rpoB under `rdp`, and
the shotgun test set with `profiler: both` and with `profiler: sylph`) were made again
with v2.6.0 into new folders. Every step that processes the reads ran as in v2.5.0, and
every ASV and taxonomy table, FASTA file, tree and read count came out byte for byte as
with v2.5.0, the `rdp` taxonomy included; every read file held the same reads (the
compressed classified-read files of the shotgun run differ only in the time stamp of
their gzip header). Four kinds of file differ, as they can between any two runs, and
none changes a result: Cutadapt's reports in `stats/cutadapt/` name the temporary folder
of `trim_primers`, which gets a random name; the raw SINTAX output
(`06.taxonomy/_sintax_raw.tsv`) and ITSx's chart of profile hits
(`05.dada2/_itsx_tmp/seqs.graph`) hold the same lines in another order; and Metaxa2's
own taxonomy file (`05.dada2/metaxa2/seqs.taxonomy.txt`), which nothing reads, quotes
for one ASV of the 16S run another of four equally good hits (91.95 % identity over 373
bases instead of 91.91 % over 371).

With `orientation: mixed` the results change on purpose, and nothing needs to be added
to the command: in every folder made with `mixed`, `trim_primers` and every step after
it run again, after a `git pull` as from a fresh clone, a release archive or a new
worktree, because v2.5.0 did not write the orientation table.

**Skipping the reruns.** `--rerun-triggers mtime params input software-env` leaves out
the check of commands and script dates. With it only MultiQC runs again, except in a
folder made with `orientation: mixed` (see below): after a `git pull` the release has
changed MultiQC's settings file, `workflow/multiqc_config.yaml`, and in a fresh clone
that file is newer than the old report; Snakemake checks the input files of a step
whatever `--rerun-triggers` says. The option holds only for the run it is given to:
Snakemake keeps the old records, so a later run without it finds the same changes.
Give it to every later run on that output folder, or let the steps run once; while it
is given, a later change of a command or a script, a future fix included, goes
unnoticed too. In a folder made with `orientation: mixed` it does not keep the
fix away: `trim_primers` and every step after it run all the same, because the
orientation table is missing, and a missing file the run asks for is made whatever
`--rerun-triggers` says. From a fresh clone with this option, a changed setting is not
noticed either, since there is no record to compare it with: switching such a folder to
`min_distinct_minimizers: auto` then needs `--forcerun kmer_evidence`. The "Running
MetaFlux" page of the documentation explains these checks ("What happens on a
re-run").

### What to edit in a config written for v2.5.0

- Nothing, for a run that stays as it was. Every new key has a default: `amplicon.coi`
  (both screens on; used by COI runs only, though `length_tolerance` is checked in every
  amplicon run), `amplicon.probe_length_stat.COI` (needed by COI runs only) and the
  `resources` rows for `coi_screens`. A folder made with `orientation: mixed` needs no
  edit either, only the re-run described above.
- `shotgun.kmer_evidence.min_distinct_minimizers` must hold a whole number or `auto`.
  A fraction, `true`/`false` or other text, which v2.5.0 took without a word or stopped
  on with a traceback, now stops the run with the message that says what is allowed.
- A COI run needs `amplicon.type: COI` and `amplicon.taxonomy.method: sintax` (`rdp` is
  refused); `keep: [k__Metazoa]` keeps the animals. The COI page has a complete example.
  The entry `COI: p95` under `amplicon.probe_length_stat` and the `resources` rows for
  `coi_screens` come from the template, which is read underneath the run config in every
  run launched from the repository root; only a `config/config.yaml` kept from v2.5.0
  lacks them. Add them from the v2.6.0 template: without `COI: p95` the run stops at
  start-up and names it, and without the rows `coi_screens` falls back to
  `threads_default` and `mem_mb_default` (4 threads and 2,000 MB as shipped), well below
  the 8.5–8.6 GB that VSEARCH holds, which matters on a cluster or with
  `--resources mem_mb=...`.

## 2.5.0 — 2026-10-16

### Taxonomy strings: `d__` is the domain, `k__` the kingdom

SILVA 144 and NCBI now carry both a domain and a kingdom for prokaryotes. MetaFlux
used to write the domain under `k__`. From this release on, `d__` is always the domain
slot and `k__` always the kingdom slot, and a table carries a kingdom only where its
reference defines one.

| Output | v2.4.0 | v2.5.0 |
|---|---|---|
| 16S, default SILVA (138.2, now 144) | `k__Bacteria;p__Pseudomonadota;…` | `d__Bacteria;k__Pseudomonadati;p__Pseudomonadota;…` |
| 16S, SILVA 138.2; rpoB | `k__Bacteria;p__…` | `d__Bacteria;p__…` |
| ITS; gyrB | unchanged | unchanged |
| 18S, `sintax` | empty for an ASV with an uncertain supergroup | at least `d__Eukaryota` |
| Kraken2, NCBI database | `k__Bacteria;p__…;s__Escherichia coli` | `d__Bacteria;k__Pseudomonadati;p__…;s__Escherichia coli` |
| Kraken2, a virus | `k__;p__Uroviricota;…;s__P37` | `d__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;…;s__Peduovirus P37` |
| Kraken2, GTDB database | `k__Archaea;p__…` | `d__Archaea;k__;p__…` |
| sylph (new) | — | `d__Bacteria;p__…` |

- **16S classifies on SILVA 144 by default**, on all of its ranks, kingdom included
  and classified with the same cutoff as every other rank; the kingdom is written
  without SILVA's `--other` suffix (`amplicon.taxonomy.silva_other` says what to do
  with the other `--other` names). SILVA 138.2 stays available with
  `references: {silva: {release: "138.2"}}`. The SILVA downloads of both releases are
  checked against the md5 that SILVA or Zenodo publishes for them; a mismatch stops
  the run.
- **18S under `sintax`: the domain is classified first.** VSEARCH reads the SINTAX
  letters as a fixed top-down order, and PR2's own file puts the domain under the
  second letter, so an ASV lost its domain whenever its supergroup was uncertain.
  MetaFlux now classifies against a copy with the two letters swapped, so the domain
  is no longer lost: on the 18S test set, 158 of the 159 ASVs that had an empty
  taxonomy now carry `d__Eukaryota`. 18S `sintax` runs made with earlier versions
  have this problem; re-run them.
- **The shotgun OTU table is built from the Bracken reports**, without kraken-biom:
  eight lineage slots on every Kraken2 database (`d__` to `s__`), NCBI's kingdoms kept,
  a virus's realm in `d__`, species names complete (`Peduovirus P37`, not `P37`), and a
  slot never filled from a neighbouring branch of the tree. At `tax_lev: S` the counts
  are the same Bracken estimates as before. At `tax_lev` `P` to `G` the v2.4.0 table was
  empty or nearly empty, without an error; it now holds one row per taxon at that rank,
  with Bracken's re-estimated reads. In `otu_table.biom` the sample `Id` is the plain
  sample name.
- **`rdp` and `sintax` tables have one format** (`""` as the first header cell, `NA`
  for a rank the classifier did not reach), ITS tables made with `rdp` write the
  species as `Genus_epithet`, and the per-rank cells of `taxon_seq_table.txt` hold
  bare names for every marker.
- The 16S `rdp` `Species` column still holds the epithet alone, as `addSpecies` gives
  it; the documentation shows the R line that joins it to the genus.

### sylph, an alternative shotgun profiler

`shotgun.profiler: kraken2 | sylph | both` chooses between Kraken2 + Bracken (the
default, unchanged), sylph + sylph-tax (genomes detected from k-mer sketches,
abundance from genome coverage, `03.abundance/sylph_table.tsv`), or both side by side
on the same cleaned reads. sylph databases and their sylph-tax lineage tables are
user-supplied (`references.sylph_databases`, `references.sylph_taxonomy_dir`) and are
checked at start-up: the taxonomy must match the database's release, GlobDB r226 is
refused (754 of its genomes have no domain), and a GTDB and a GlobDB database are
refused in the same run (GlobDB already holds the GTDB genomes). Read extraction and
the k-mer evidence gate work on Kraken2's output and are switched off under `sylph`.
The sylph tables are written unfiltered. To add sylph to a finished Kraken2 run, use a
new output folder; MetaFlux stops otherwise and says why.

### Checks at start-up, before any job runs

- **Filter tokens.** Every `keep` and `discard` token, in both modes and also with the
  filter switched off, is checked against the slots the run will write and, for
  shotgun, against the Kraken2 database's domains, kingdoms and realms. A token that
  can never match stops the run and the message names the one to write
  (`k__Bacteria` → `d__Bacteria`); `Viruses` and a realm under `k__` are refused; a
  name this database lacks only warns.
- **Reference files.** Each taxonomy reference is checked, header by header, against
  the rank set and the name rules (`reference_values`) its marker pack declares,
  before any classification (`stats/taxonomy_reference_check.json` records what was
  seen). A file of another release, SILVA's raw download, PR2's own letter order or
  an old converted rpoB file stop the run with the file, the slot and example values
  named.
- **Marker packs.** A pack in `config/markers/` written for an older MetaFlux stops
  every amplicon run; a pack for a marker of your own needs `reference_values`
  (`{}` checks the rank shape only). A pack whose SINTAX letters break VSEARCH's
  order stops too.
- **Bracken rank.** `shotgun.bracken.tax_lev` must be one of `P`, `C`, `O`, `F`, `G`
  and `S`. `D` and `K` are refused (Bracken files the reads of taxa without the rank
  under the wrong taxon), and so is any other value, such as a lower-case letter, which
  used to fail in Bracken only after Kraken2 had run. `P`, `C`, `O` and `F` give a
  warning. `tax_lev` is read only when Kraken2 runs.
- `discard_min_support: true` is refused (YAML reads it as a yes/no value, not a
  number).

### Taxon filter

- `amplicon.taxonomy.filter.discard_min_support` (default 50): an ASV is also removed
  when the classifier's best hit carries a `discard` token at that rank with at
  least this support, even when the rank is below the cutoff and missing from the
  reported lineage; it catches mitochondrial 16S reported only to phylum or class.
  `null` gives the previous behaviour.
- `amplicon.taxonomy.filter.organelles: metaxa2` (16S with extraction): the ASVs that
  Metaxa2 calls mitochondria or chloroplast are removed as well, which catches plant
  mitochondria the classifier reports as bacteria. `06.taxonomy/asv_table_unfiltered.txt`
  lists every ASV with its filter status, organelles included.
- The shotgun filter stops the run when it removes every taxon of a table that had
  some, instead of writing an empty table, and lists the `d__` values the table holds.

### Output folders

Two-digit folder numbers in both modes: the amplicon folders `1.reads` … `7.phylogeny`
become `01.reads` … `07.phylogeny`, and the read-count files in `stats/` take the same
two digits (`01.raw_reads.counts` …). Shotgun mode already used two digits; its
`02.classification/` now holds one folder per profiler, `kraken2/` (the Kraken2 files,
until now directly in `02.classification/`, and the k-mer evidence files of
`02b.evidence/`) and `sylph/`. `01.preprocessing`, `03.abundance`,
`04.extracted_reads`, `stats/`, `logs/`, `benchmarks/` and `multiqc/` keep their names.

### Length probe, read extraction, read tracking

- **rpoB length probe:** in-silico products that still hold a second full-length
  forward-primer site are left out of the statistics (in many rpoB genes the primer
  also fits a look-alike stretch upstream, which put p95 at 1,332 bp; without them
  p95 is 395 bp, close to the 387 bp amplicon). The probe cache name records this.
- **16S length probe on SILVA 138.2:** chloroplast, mitochondrial and eukaryote
  references are left out, as on SILVA 144.
- **Read extraction** (`shotgun.extract_taxa`): a name must be exactly the name of one
  row of the sample's Kraken2 report, at any rank, so `Prevotella` no longer picks
  `Prevotellaceae`; two-word names and sample names with `_` work; a genus and a
  subgenus of the same name count as the genus; the reads are taken through
  Kraken2's own report, so strain-level reads are included; the per-job log keeps
  only KrakenTools' summary lines instead of hundreds of megabytes of progress
  counters.
- **Read tracking:** `nophix` counts exactly the read pairs BBDuk kept, so `raw` minus
  `nophix` is the pairs it removed (a pair with PhiX in one mate only used to count as
  half a pair); MultiQC labels BBDuk's statistics "PhiX removal".

### k-mer evidence floor, `allow_one_off` description

- **`shotgun.kmer_evidence.min_distinct_minimizers` defaults to 300** instead of 333.
  The published floor of 1,000 unique k-mers converts to about 333 distinct
  minimizers, and the default rounds that down because the conversion is approximate.
  The gate is off by default, so the default OTU tables do not change; the evidence
  table and the gated report written on every Kraken2 run now use 300. With the gate
  on, species with 300 to 332 distinct minimizers now pass. The k-mer evidence page
  describes a floor scaled to sequencing depth.
- **`amplicon.dada2.chimera.allow_one_off`:** the comment in `config/config.yaml` and
  the configuration page had the two settings the wrong way round. `false`, DADA2's
  default, flags exact bimeras only; `true` also flags ASVs one mismatch or indel away
  from an exact bimera, so it removes more. The pipeline's behaviour does not change,
  and runs at the default need no re-run. A config that set `true` to relax the check,
  as the old text advised, removed more ASVs instead: set it back to `false` and
  re-run.

### Pinned tool versions

- **Pinned conda environments on Linux.** Next to every `workflow/envs/<tool>.yaml` sits
  a `<tool>.linux-64.pin.txt` with the exact package builds MetaFlux was tested with. On
  Linux (x86-64) Snakemake installs each environment from that file, without solving,
  so two installations of the same release run the same software down to the build;
  other platforms solve from the yaml as before. The installation page says how to
  refresh a pin after editing a yaml.
- **Kraken2 2.17.1.** Releases before 2.17 write NCBI's 2025 domains as `R1`/`R2`
  instead of `D`; the table step stops on such reports.
- **Bracken 3.1p1,** which fixes Bracken's own count column above species level.
- **bbmap 39.52.** BBDuk in bbmap 39.76 to 40.02 loses track of read pairs on
  correctly paired reads ("List size mismatch") or hangs. BBDuk and BBMap also write
  the reads they keep in input order (`ordered=t`), so a run is reproducible.
- **VSEARCH 2.32.0,** which gives every ASV its own random stream, so `sintax`
  results are the same at any thread count. Compared with an older VSEARCH many
  confidence values move slightly, and a lineage near the cutoff can gain or lose a
  trailing rank.

### What runs again in a folder made with v2.4.0

Every step: the output folders have new names, and the pin files count as a software
change for every step that runs in a conda environment, which also catches the steps
whose files kept their names (the raw-read QC in `stats/`, the primer reverse
complements, the reference files already on disk and, in BioProject mode, the SRA
downloads; only `samples.tsv` is kept). A run started once with `--rerun-triggers
mtime params input code` keeps those files. The old folders and count files are left
as they are and can be deleted. Reference files: the SILVA 144 files
are downloaded and built once; the 18S PR2 SINTAX file is built once; the rpoB
reference is built again under a new name (the FROGS archive is downloaded again);
the SILVA 138.2 SINTAX file is rebuilt under its old name, with the same content, so
a `refdb` folder shared with a v2.4.0 checkout keeps working for both. Results that
change on purpose are those listed above; besides them, an `rdp` re-run can move a
few calls at `min_boot`, as between two identical runs. The documentation's
"Running MetaFlux" page has the full note.

### What to edit in a config written for v2.4.0

- Filter tokens that write a domain under `k__`: `k__Bacteria` → `d__Bacteria`,
  `k__Archaea` → `d__Archaea`, `k__Eukaryota` → `d__Eukaryota`. The run stops at
  start-up and names the token to write. ITS keeps `k__Fungi`.
- To stay on SILVA 138.2 for 16S, add `references: {silva: {release: "138.2"}}`;
  without it the run uses SILVA 144 and classifies the kingdom.
- `shotgun.bracken.tax_lev: D` is refused; use `S` (the default) or another rank.
- A marker pack of your own in `config/markers/`: copy the shipped pack again and
  re-apply your change; a pack for a new marker needs `reference_values`.
- `discard_min_support: true` is refused; write a number from 0 to 100, or `null`.
- Under `resources`, the OTU-table rule is now `finalize_otu_table`; a `kraken_biom`
  line is ignored and can go.
- A config that sets `shotgun.kmer_evidence.min_distinct_minimizers: 333` keeps that
  floor; write `300`, or delete the line, for the new default.
- New keys with safe defaults, so an old config parses unchanged: `shotgun.profiler`
  (`kraken2`), `references.sylph_databases` and `sylph_taxonomy_dir`,
  `amplicon.taxonomy.silva_other` (`verbatim`), `filter.discard_min_support` (50),
  `filter.organelles` (`metaxa2`), and `resources` rows for the sylph rules.

### Documentation and repository

New sylph page; the marker, taxon-filter, output, configuration, running and
troubleshooting pages rewritten for the domain and kingdom slots, the start-up checks
and the new folders; the docs deployment runs a strict build before publishing. Every
example names the dataset and the database build it comes from; the test-datasets
page lists the shotgun test dataset and the benchmark data, with numbers from runs
made with this release's defaults; the shotgun figures are recomputed on this
release's runs, and the benchmark scorer was checked against OPAL. The threshold page
and the config comment no longer advise a smaller `shotgun.bracken.threshold_alpha`
for defined mixtures: on the ZymoBIOMICS D6331 mock the best value is 4.5 to 7 times
the default of 5×10⁻⁵ (the page gave 2.4×10⁻⁵ before); the default is unchanged. The
repository
gains a pin file for every conda environment (`workflow/envs/*.linux-64.pin.txt`), a
test suite (`tests/`), a dry-run matrix with tiny test data (`.test/`), a CI
workflow that runs both with a strict docs build on every push, contribution
guidelines, a code of conduct, support notes, issue forms, this changelog and a
`CITATION.cff`.

## 2.4.0 — 2026-09-21

- **Optional 16S phylogeny stage** (`amplicon.phylogeny`, off by default; stage 85,
  `7.phylogeny/`): the retained ASVs are aligned with MAFFT and a maximum-likelihood
  tree is built with IQ-TREE 3 (default), FastTree or RAxML-NG, exported unrooted as
  Newick with every ASV of the table as a tip, with provenance and QC files. Pruning,
  rooting and diversity statistics stay in R; the docs carry a tested example.
- **The amplicon path is reproducible from the raw reads to the tree:** bowtie2 runs
  with `--reorder`, the phylogeny input is sorted by sequence, and the tree search
  runs on one thread by default (`resources.threads.phylo_tree`). Under `sintax` on
  more than one thread VSEARCH still drifted near the cutoff (fixed in 2.5.0 by the
  VSEARCH pin).
- New dependencies, created only when the phylogeny stage is enabled: mafft, iqtree,
  fasttree, raxml-ng, dendropy.
- Docs: *Phylogeny (16S)* and *Phylogeny validation* pages, *Choosing a mode*
  rewritten, acknowledgements.
- Upgrading: configs parse unchanged; a leftover `phylogeny.root` key is rejected
  with an explanation. A re-run does not reproduce a v2.3.1 run byte for byte: ASVs
  tied in abundance may be numbered differently, so compare runs by sequence.

## 2.3.1 — 2026-08-19

- The PhiX reference is fetched over HTTPS instead of FTP, which many networks block.
  The file is byte-for-byte the same, so nothing has to be re-run. A personal run
  config keeps the old URL: change `ftp://` to `https://` in `references.phix.fetch_url`.

## 2.3.0 — 2026-08-19

- Shotgun: a species-level k-mer evidence gate (`shotgun.kmer_evidence`, off by
  default) that removes taxa without enough distinct k-mers before Bracken, and a
  depth-scaled Bracken read threshold (`shotgun.bracken.threshold: auto`).
- The documentation site is deployed to GitHub Pages automatically.
- `sintax`: a malformed VSEARCH output line stops the run instead of being guessed at.

## 2.2.0 — 2026-08-19

- The documentation site (MkDocs); the README is slimmed to a pointer at it.
- Automatic `truncLen` is never placed above the bulk of the reads
  (`min_read_coverage_pct`), and `raise_trunc` allocates by base quality.
- `addSpecies` gets the same reverse-complement chance as `assignTaxonomy`.
- The DADA2 filter and merge parameters move to the package defaults.
- Comments and docstrings expanded across the whole pipeline; fixes from a
  full-codebase audit; the tagline reads "short-read" instead of "Illumina".

## 2.1.0 — 2026-07-24

- Two protein-coding markers: **gyrB** (INRAE DD7RZ8 v6 trainset; the reference's gene
  tag doubles as a paralog filter, `discard: [other]`) and **rpoB** (FROGS RefSeq build,
  converted to a DADA2 trainset). Both `rdp` only, both verified on their published
  mock communities.
- Reference-preparation scripts grouped under `workflow/scripts/refdb/`; `amplicon.type`
  matches marker packs case-insensitively.

## 2.0.0 — 2026-07-23

- **Breaking:** `amplicon.taxonomy.filter`'s `include_pattern` / `exclude_pattern` are
  retired in favour of rank-aware `keep` / `discard` token lists, and no marker gets an
  implicit default filter any more.
- The 18S rRNA marker (PR2 taxonomy, SILVA-Euk length probe).
- Marker packs (`workflow/markers/*.yaml`): marker facts and reference URLs move out
  of the code into data that a user can shadow.
- `amplicon.seed` fixes the random number generators of DADA2 and VSEARCH.
- Automatic `truncLen` fixed for markers other than 16S.

## 1.0.0 — 2026-06-26

- First release: amplicon mode (16S and ITS, DADA2) and shotgun mode (Kraken2 +
  Bracken) in one Snakemake workflow, with one config file, one command and one
  MultiQC report.
