# Changelog

What changed in each MetaFlux release, newest first. Each entry says what changes
for a user, what runs again in an output folder made with the previous release, and
what to edit in an existing run config. The documentation site has the details; the
release notes on GitHub repeat the longer entries. Every release is archived on
Zenodo; the concept DOI in the README always resolves to the newest one.

## 2.6.1 — [To be settled before release: the date]

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
