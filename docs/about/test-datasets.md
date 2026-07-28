# Test datasets

Every marker MetaFlux supports is exercised on a small published dataset before a
release. The runs are deliberately small — four to six samples each — so the whole
amplicon path can be re-run in minutes after a change and the output compared against
the previous result.

All five sets are public and reachable from the accessions below, so any result on this
page can be reproduced from scratch. The run configs live in `config/` and are not
tracked by git, since they carry local paths.

| Marker | Dataset | Samples | Primers | Reference |
|---|---|--:|---|---|
| 16S (V5–V7) | [PRJNA305879](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA305879) — plant-associated bacterial communities, MiSeq 2×300 | 6 | 799F / 1175R | SILVA v138.2 |
| ITS (ITS2) | [PRJNA862334](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA862334) — yeast mock communities, MiSeq 2×300 | 6 | ITS3 / ITS4 | UNITE |
| 18S (V4) | [PRJDB12231](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB12231) — aquatic mesocosms ± calcined dolomite, MiSeq 2×300 | 4 | TAReuk454FWD1 / TAReukREV3 | PR2 (SILVA-Euk probe) |
| gyrB | [PRJNA459479](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA459479) — 15-strain mock communities | 5 | F64 / R353 | DD7RZ8 v6 |
| rpoB | [PRJEB81812](https://www.ebi.ac.uk/ena/browser/view/PRJEB81812) — 19-strain synthetic mock | 5 | Univ_rpoB_deg | FROGS RefSeq |

The three mock-community sets (ITS, gyrB, rpoB) have known composition, so they check
taxonomy as well as mechanics. The ITS run recovers the expected yeast genera —
*Hanseniaspora*, *Rhodotorula*, *Debaryomyces*, *Yamadazyma* — with the remainder
resolved to Saccharomycetes.

The 16S and 18S sets are environmental rather than defined, so they test the mechanics
of their path — probe measurement, truncation, merging, rank handling in SILVA and PR2 —
and are judged on whether the output is coherent rather than against a known answer.
For 18S that check is that 98.5% of reads classify as Eukaryota, splitting into TSAR
(53%), Obazoa (28%) and Archaeplastida (14%), and that the observed ASV median length
of 381 bp lands on the 381 bp median predicted by in-silico PCR against PR2. Defined
18S mock communities are scarce; the ground-truth work is carried by the other three
sets.

## Check a public dataset before adopting it

Assembling these five sets meant rejecting three others, each for a different reason,
and in every case the submission metadata looked fine. It is worth running the same
three checks on any BioProject before building a run around it — they cost a few
minutes and each of them caught something.

**Do the reads still carry their primers?** MetaFlux trims primers with cutadapt under
`--discard-untrimmed`, so a library deposited *after* primer removal loses every read
at that step. [PRJNA1138067](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1138067), an
otherwise ideal two-species ITS2 mock, is deposited primer-trimmed and was rejected.

To tell without trusting the metadata, count the most common 22-mer at the 5′ end of R1
across a few thousand reads. A primer is the same in every sample; biological sequence
tracks the community. In that rejected dataset the 5′ end split 76% / 20% between two
mock samples — the ratio of the two species, not a primer. The adopted 18S set passes
the same test cleanly: TAReuk454FWD1 opens 98.5% of R1 and TAReukREV3 98.1% of R2.

**Is it genuinely paired?** [PRJNA326072](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA326072)
is labelled PAIRED but was archived R1-only. The tell is in the run metadata before
downloading anything: every run reports `avgLength` around 290, where a real 2×300
paired run reports around 600. Converting it produces one full file and one file of
zero-length reads. Its sister deposit PRJNA305879, from the same paper, is properly
paired — 27 of 27 runs at `avgLength` 600–601 — and is the set now used.

**Is it the platform you think?** PRJEB23971 was proposed as an 18S candidate on the
strength of a description reading "Illumina MiSeq". It is 454 GS FLX Titanium, and
single-end. Read the `platform` and `library_layout` fields from the ENA or SRA run
table rather than the free-text study description.

## Effect of moving the merge parameters to DADA2's defaults

`dada2.merge.min_overlap` and `max_mismatch` were moved to DADA2's own defaults
(12 and 0, from 10 and 2). Every marker was then run twice on the set above, once with
each pair of values and nothing else changed, to measure what that costs. Read counts
are summed across the samples of each set; ASV counts are rows of the final
`6.taxonomy/asv_table.txt`.

| Marker | Reads merged, at 10/2 | at 12/0 | change | ASVs at 10/2 | at 12/0 | change |
|---|--:|--:|--:|--:|--:|--:|
| 16S | 214,301 | 211,268 | −1.4% | 273 | 211 | −22.7% |
| ITS | 350,726 | 233,114 | **−33.5%** | 29 | 12 | −58.6% |
| 18S | 848,352 | 842,971 | −0.6% | 1,508 | 1,268 | −15.9% |
| gyrB | 275,075 | 267,686 | −2.7% | 256 | 241 | −5.9% |
| rpoB | 139,113 | 138,734 | −0.3% | 48 | 44 | −8.3% |

Three things are worth reading out of this.

**Every read counted through `filterAndTrim` and denoising is identical** in both runs
of every marker. The change appears only at the merge step, which is where these two
parameters act. The `min_len` change that accompanied them (50 → 20, DADA2's default)
altered nothing either, because reads are far longer than either floor by the time the
filter sees them.

**ASV counts fall much further than read counts.** For 16S a 1.4% read loss removes 23%
of ASVs. That is the expected shape rather than a problem: a pair that merges only
because two mismatches were tolerated in the overlap tends to produce a sequence of its
own, carrying very few reads. So the discarded pairs are concentrated in rare ASVs.
DADA2's `maxMismatch = 0` is deliberate for exactly this reason — once denoising has
already corrected sequencing error, a mismatch left in the overlap is evidence the two
mates do not belong together, rather than noise to be tolerated.

**ITS is the marker that actually pays for this**, and the reason is biological. ITS2
length varies a great deal between fungal taxa, so with 2×300 reads the longer
amplicons in a community overlap by only a little. Those marginal pairs are precisely
the ones a 10 → 12 bp requirement removes, and the short overlap they do have sits in
the low-quality 3′ tails of both reads, where a mismatch is most likely — so
`maxMismatch: 0` removes many of the rest. A third of the merged reads went with them.

On this mock it cost no taxa: all four expected genera are still recovered from the 12
remaining ASVs, along with the two Saccharomycetes sequences that carry the unresolved
remainder. A defined community of four genera is a soft test of that, though — on a
richer or more length-variable community the same loss could take real taxa with it. If
an ITS run merges far fewer reads than expected, these are the two parameters to look at
first, and `min_overlap: 10` with `max_mismatch: 2` is a defensible choice for ITS.

## What the two truncation paths are tested on

Between them the 16S and 18S sets exercise both ways `trunc_len.mode: auto` can end up,
which is the main reason they are kept as they are.

### 18S: the set that produced the read-coverage cap

The 18S library is clean enough on R1 that its aggregated lower quartile never falls
below `q_threshold`, so the picker has no quality signal to cut on and takes the ceiling
instead. That ceiling used to come from the longest read present — 288 bp, from the one
sample holding R1 reads that long, while half the reads stop at 281. Since
`filterAndTrim` discards anything shorter than `truncLen`, the cut kept 9 of 908,768
pairs and the run came out with 3 ASVs.

The working fix at the time was `mode: manual` at 270 / 239, worked out by measuring the
read-length distribution by hand. That this was necessary at all was the argument for
`amplicon.trunc_len.min_read_coverage_pct`, which now caps each sample at the length 95%
of its reads still reach before quality is consulted. On this library it brings the R1
ceiling to 271, and `auto` returns **271 / 239** unaided — within 0.02% of the
hand-picked values, carrying 90.9% of raw reads through to the final table.

The config is therefore back on `mode: auto`, deliberately: it is the working test that
the guard does its job. Note also that only R1 was ever the problem. R2's quality did
degrade, the picker found an honest cut at 239, and that number is unchanged throughout —
a reminder to treat the two directions separately rather than overriding both.

The full diagnosis is in
[Troubleshooting](../troubleshooting.md#amplicon-automatic-trunclen-can-drop-all-reads-on-very-clean-runs).

### 16S: the overlap-constrained path

The 16S set is the opposite case and runs on `mode: auto` unchanged. Its quality does
fall away towards the 3′ end, so the picker finds honest quality-based cuts — but they
come out too short to leave the 12 bp of overlap the 384 bp amplicon needs. That
triggers `resolve_policy: raise_trunc`, which extends the cuts back out until the pair
can merge, and the run settles at 242 / 154 with the overlap requirement met exactly.

This makes the set the working test of that resolution path, which is otherwise easy to
leave unexercised: most libraries either merge comfortably or fail outright.
