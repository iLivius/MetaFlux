# Test datasets

Every marker MetaFlux supports is exercised on a small published dataset before a
release. The runs are deliberately small — four to six samples each — so the whole
amplicon path can be re-run in minutes after a change and the output compared against
the previous result.

All the data on this page are public and reachable from the accessions and references
given, so any
result on it can be reproduced from scratch with the primers, references and databases
listed. Exact counts of ASVs and of reads in the final table also depend on the reference
release and the taxon-filter settings, and both have changed between MetaFlux versions,
so treat the counts on this page as indicative: a re-run can differ somewhat. Within each
comparison below, both runs used the same reference and taxon-filter settings, so the
difference between them holds even where the absolute counts have since moved.

| Marker | Dataset | Samples | Runs used | Primers | Reference |
|---|---|--:|---|---|---|
| 16S (V5–V7) | [PRJNA305879](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA305879) — plant-associated bacterial communities, MiSeq 2×300 | 6 of 27 | SRR3032384, SRR3032389, SRR3032393, SRR3032396, SRR3032402, SRR3032405 | 799F / 1175R | SILVA 144 (138.2 selectable) |
| ITS (ITS2) | [PRJNA862334](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA862334) — yeast mock communities, MiSeq 2×300 | 6 of 36 | SRR20648162, SRR20648163, SRR20648167, SRR20648181, SRR20648185, SRR20648187 | ITS3 / ITS4 | UNITE |
| 18S (V4) | [PRJDB12231](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB12231) — aquatic mesocosms ± calcined dolomite, MiSeq 2×300 | 4 | DRR316761–DRR316764 | TAReuk454FWD1 / TAReukREV3 | PR2 (SILVA-Euk probe) |
| gyrB | [PRJNA459479](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA459479) — 15-strain mock communities | 5 of 40 | SRR7127623–SRR7127625, SRR7127639, SRR7127640 | F64 / R353 | DD7RZ8 v6 |
| rpoB | [PRJEB81812](https://www.ebi.ac.uk/ena/browser/view/PRJEB81812) — 19-strain synthetic mock | 5 of 96 | ERR13956464–ERR13956468 | Univ_rpoB_deg | FROGS RefSeq |

The three mock-community sets (ITS, gyrB, rpoB) have known composition, so they check
taxonomy as well as mechanics. The ITS run (UNITE release 19.02.2025, `method: sintax`)
recovers four of the expected yeast genera — *Hanseniaspora*, *Rhodotorula*,
*Debaryomyces*, *Yamadazyma* — but stops at class (Saccharomycetes) for two ASVs that
between them carry 43% of all reads, one of them alone accounting for 42%. That is the
largest component of the run. BLAST against NCBI nt and a direct search against the
UNITE reference used here both place the two sequences at 99% identity or more to
*Saccharomyces*, most likely *S. cerevisiae*. The reference contains this species, but
VSEARCH's `--sintax` bootstrap gives no consensus at order, family or genus for these
sequences, probably because several *Saccharomyces* species have near-identical ITS2 in
this amplicon. On this mock it is therefore a classification gap for what is probably
its most abundant species. See
[the merge-parameter section below](#effect-of-moving-the-merge-parameters-to-dada2s-defaults)
for how this same gap shows up under the other `dada2.merge` setting too, as a different
set of ASVs carrying a similar fraction of reads.

The 16S and 18S sets are environmental rather than defined, so they test the mechanics
of their path — probe measurement, truncation, merging, rank handling in SILVA and PR2 —
and are judged on whether the output is coherent rather than against a known answer.
For 18S (PR2 v5.1.1, `method: sintax`) that check is that 99.97% of reads classify as
Eukaryota before the `keep: [d__Eukaryota]` filter, splitting into TSAR (53%), Obazoa
(28%) and Archaeplastida (14%), and that the observed median ASV length of 382 bp lies
within a base of the 381 bp median predicted by in-silico PCR against SILVA-Euk 18S v132,
the probe reference. Defined 18S mock communities are scarce; the ground-truth work is
carried by the other three sets.

## Shotgun test dataset

Shotgun mode is tested on four metagenomes from BioProject
[PRJNA1477390](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1477390): microbial consortia
enriched on naproxen from wastewater treatment plants, sequenced paired-end on an Illumina
NextSeq 500. The BioProject also holds metatranscriptomes; only the four metagenome runs
are used.

| Runs | Read pairs per run | Kraken2 index | sylph databases |
|---|--:|---|---|
| SRR39106521, SRR39106522, SRR39106523, SRR39106524 | 4.6–5.9 million | PlusPF `k2_pluspf_20260626` | `gtdb-r232-c200-dbv1` and `fungi-refseq-2025-10-11-c200-dbv1` |

The cultures have no known composition, so the run checks the mechanics of shotgun mode
rather than its accuracy. The real-data examples on the
[k-mer evidence](../shotgun/kmer-evidence.md) and [sylph](../shotgun/sylph.md) pages come
from these four runs and from the D6331 run listed under [Benchmark data](#benchmark-data).

## Benchmark data

Accuracy figures for shotgun mode come from communities of known composition:

- the ZymoBIOMICS Gut Microbiome Standard (D6331), SRA run
  [SRR38976715](https://www.ncbi.nlm.nih.gov/sra/SRR38976715) (BioProject PRJNA1474435),
  Illumina NovaSeq 6000, 8.1 million read pairs, scored against the composition in
  Zymo's instruction manual;
- CAMI II short-read samples with their gold-standard profiles
  ([Meyer et al. 2022](https://doi.org/10.1038/s41592-022-01431-4)): marine samples 0–9,
  plant-associated (rhizosphere) samples 0–20 and strain-madness samples 0–9.

Zymo's manual lists the standard's *Fusobacterium* as *F. nucleatum* 2/1/50A, but sylph
(GTDB r232), MetaPhlAn 4 and Kraken2 (PlusPF `k2_pluspf_20260626` and GTDB r226) put
nearly all of its reads in *F. animalis*, formerly *F. nucleatum* subsp. *animalis*, so
*F. animalis* is scored as the true member and calls of *F. nucleatum*, *F. polymorphum*
and *F. vincentii* count as false positives.

The mock also contains a little human DNA, a contaminant of the sample: 933 read pairs,
0.02 % of the classified pairs with PlusPF `k2_pluspf_20260626`. Kraken2 reports it as
*Homo sapiens* because PlusPF contains the human genome; MetaPhlAn cannot report it,
because its database holds microbial markers only. It is not scored. MetaFlux's
[host removal](../shotgun/decontamination.md#host-removal) is the step that removes such
reads before classification.

## Check a public dataset before adopting it

Assembling the five amplicon sets meant rejecting other candidates. The checks below are
the ones that caught them. A dataset can fail in other ways too, so check its run table
and a sample of its reads before building a run on it.

| Check | How to run it | What it caught |
|---|---|---|
| **Do the reads still carry their primers?** | On a few thousand reads from more than one sample, check that R1 starts with the forward primer and R2 with the reverse primer, allowing for degenerate bases. The commonest 22-mer at the 5′ end of R1 shows what is there instead: if it changes with the organisms in each sample, the primer is gone. | [PRJNA1138067](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1138067), ITS2 mocks of two fungi in different proportions. The archived reads carry no primer: each read starts at the first base after its primer site, and the commonest 5′ 22-mers of R1 are the 5.8S sequences of the two fungi (75% and 21% of reads in SRR29900459). MetaFlux trims primers with `--discard-untrimmed`, so a library like this loses every read at that step. |
| **Is it really paired?** | In the SRA run table, check that `spots_with_mates` is not 0 and that `avgLength` matches the stated read length and layout: a 2×300 paired run reads about 600, not about 300. | A 16S V5–V7 deposit labelled `PAIRED`. The forward and reverse reads of each sample were deposited as two separate runs with one read per spot (`spots_with_mates` 0, `avgLength` at most 300), so each run downloads as one mate only. fasterq-dump writes a single FASTQ for such a run and reports the second read of every spot as zero-length. |
| **Is it the platform and layout you expect?** | Read the platform and layout fields of the ENA or SRA run table, not a paper or a summary of the study. | [PRJEB23971](https://www.ebi.ac.uk/ena/browser/view/PRJEB23971), 18S V4: 454 GS FLX Titanium, single-end. [PRJEB24936](https://www.ebi.ac.uk/ena/browser/view/PRJEB24936), rpoB and 16S mock communities on Illumina MiSeq, deposited as single-end reads: the rpoB mock runs average 431–434 bases per read, the length of the whole amplicon, so the pairs had already been merged. |

## Effect of moving the merge parameters to DADA2's defaults

MetaFlux uses DADA2's own defaults for `dada2.merge.min_overlap` and `max_mismatch`
(12 and 0). To measure what they cost against the looser 10 and 2, every marker was run
on its set with each pair of values and nothing else changed. These runs were made when
the defaults changed, in v2.2.0, with 16S classified against SILVA 138.2 and truncation
set by hand for gyrB (210 / 170) and rpoB (250 / 200). With the v2.5.0 defaults the 16S,
gyrB and rpoB counts differ. Read counts are summed across the samples of each set; ASV
counts are rows of the final `06.taxonomy/asv_table.txt`, for 18S counted before the
`keep: [d__Eukaryota]` filter.

| Marker | Reads merged, at 10/2 | at 12/0 | change | ASVs at 10/2 | at 12/0 | change |
|---|--:|--:|--:|--:|--:|--:|
| 16S | 214,301 | 211,268 | −1.4% | 273 | 211 | −22.7% |
| ITS | 350,726 | 233,114 | **−33.5%** | 29 | 12 | −58.6% |
| 18S | 848,352 | 842,971 | −0.6% | 1,508 | 1,268 | −15.9% |
| gyrB | 275,075 | 267,686 | −2.7% | 256 | 241 | −5.9% |
| rpoB | 139,113 | 138,734 | −0.3% | 48 | 44 | −8.3% |

The table shows three things.

**Every read counted through `filterAndTrim` and denoising is identical** in both runs
of every marker. The change appears only at the merge step, which is where these two
parameters act. The `min_len` change that accompanied them (50 → 20, DADA2's default)
altered nothing either, because reads are far longer than either floor by the time the
filter sees them.

**ASV counts fall much further than read counts.** For 16S a 1.4% read loss removes 23%
of ASVs. That is expected: a pair that merges only because two mismatches were tolerated
in the overlap tends to form a sequence of its own with very few reads, so the discarded
pairs are concentrated in rare ASVs. DADA2 sets `maxMismatch = 0` for this reason: after
denoising has corrected sequencing errors, a mismatch left in the overlap suggests that
the two mates do not belong together.

**ITS loses the most reads, for a biological reason.** ITS2 length varies widely between
fungal taxa, so with 2×300 reads the longer amplicons in a community overlap by only a
few bases. Raising the minimum overlap from 10 to 12 bp removes those pairs, and because
their short overlap lies in the low-quality 3′ ends of both reads, where mismatches are
most likely, `maxMismatch: 0` removes many of the others. Together the two changes
removed a third of the merged reads.

On this mock it cost no taxa: all four expected genera are still recovered from the 12
remaining ASVs. The class-only fraction described above (almost certainly *Saccharomyces*)
does not come from this parameter change either. In a run with the looser 10/2 setting
(made with v2.2.0) it is 7 separate ASVs rather than 2 (the looser overlap tolerance
lets more marginal read variants through, so the same underlying signal is split across
more low-abundance ASVs), carrying 42.5% of that run's final reads, against 43.3% in the
v2.5.0 run.
The exact ASVs differ, but the fraction stuck at class level is essentially unchanged —
this classification gap is a property of the reference and classifier for this organism,
not of the merge parameters. A defined community of only four confidently-named genera is
a weak test of that distinction, though — on a richer or more length-variable community
the same read loss could take real taxa down with it, not just leave them unresolved. If
an ITS run merges far fewer reads than expected, these are the two parameters to look at
first, and `min_overlap: 10` with `max_mismatch: 2` is a defensible choice for ITS.

## What the two truncation paths are tested on

Between them the 16S and 18S sets exercise both ways `trunc_len.mode: auto` can end up,
which is the main reason they are kept as they are.

### 18S: the set that needs the read-coverage cap

The 18S library is clean enough on R1 that its aggregated lower quartile never falls
below `q_threshold`, so the picker has no quality signal to cut on and takes the ceiling
instead. Before the read-coverage cap was added in v2.2.0, that ceiling came from the
longest reads present: 288 bp for R1, set by the sample whose longest R1 read is 288 bp
(the other three reach 289), while 70% of R1 reads stop at 281. Since `filterAndTrim`
discards anything shorter than `truncLen`, 9 of 908,768 pairs were long enough for the
288 / 239 cut, and the run ended with 3 ASVs. A hand-picked `mode: manual` at
270 / 239 avoids that; `amplicon.trunc_len.min_read_coverage_pct` makes it unnecessary
by capping each sample at the length 95% of its reads still reach before quality is
consulted. On this library it brings the R1 ceiling to 271, and `auto` returns
**271 / 239** unaided, carrying 90.9% of raw reads through to the ASV table — 828,018
read pairs reach the length filter, against 828,170 for the hand-picked 270 / 239, a
difference of 0.02%. (With the recommended `keep: [d__Eukaryota]`, a further 214 pairs
leave the final table, in six ASVs that have no domain or the nucleomorph label
`d__Eukaryota:nucl`.)

The 18S test run therefore uses `mode: auto` on purpose, as the test that the
read-coverage cap works. Only R1 was affected: R2's quality does fall, the picker places
a quality-based cut at 239, and that value is the same in all the runs above, so R1 and
R2 should be assessed separately rather than both set by hand.

The full diagnosis is in
[Troubleshooting](../troubleshooting.md#amplicon-automatic-trunclen-can-drop-all-reads-on-very-clean-runs).

### 16S: the overlap-constrained path

The 16S set is the opposite case and runs on `mode: auto` unchanged. Its quality does
fall away towards the 3′ end, so the picker finds quality-based cuts, but they come out
too short to leave the 12 bp of overlap the amplicon needs (386 bp, the p95 of
the in-silico products on SILVA 144; 384 bp on SILVA 138.2). That triggers
`resolve_policy: raise_trunc`, which extends the cuts back out until the pair can
merge, and the run settles at 244 / 154 with the overlap requirement met exactly
(242 / 154 on SILVA 138.2).

This makes the set the working test of that resolution path, which is otherwise easy to
leave unexercised: most libraries either merge comfortably or fail outright.
