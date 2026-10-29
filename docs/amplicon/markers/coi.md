# COI (mitochondrial cytochrome c oxidase subunit I)

COI is the animal barcode: the 5′ part of the mitochondrial gene for subunit I of
cytochrome c oxidase, about 650 bp in the classic barcode, of which metabarcoding primers
amplify 200–420 bp. It separates most animal species, and it is the standard marker for
surveys of arthropods, other invertebrates and fish, from bulk samples, Malaise traps or
environmental DNA.

The MetaFlux path for COI is the ordinary amplicon path — PhiX removal, primer trimming,
DADA2 denoising, length filtering, taxonomy — with three things that are specific to COI:

- the reference is **MIDORI2**, converted once into the rank model MetaFlux uses;
- the taxonomy is assigned with **SINTAX only**;
- two **COI screens** run before the taxon filter: a reading-frame check against
  pseudogenes and nuclear copies of the gene (NUMTs), and an off-target guard that moves
  bacterial COX1 sequences, *Wolbachia* above all, out of the main tables.

## At a glance

| | |
|---|---|
| `amplicon.type` | `COI` |
| Reference database | MIDORI2 release GB273 (Leray et al. 2022), variant UNIQ NUC_SP, SINTAX file; Zenodo version DOI 10.5281/zenodo.22685701 |
| Role of that database | taxonomy reference, amplicon-length probe substrate, and identity E of the off-target guard |
| Second reference | the bacterial COX1 reference set `bacterial_cox1_gtdb-r232.v1`, built from GTDB R232: identity B of the off-target guard, and the list of MIDORI2 records at least 90 % identical to a candidate gene of the set (186 in GB273), each marked `removed` (170, bacterial COX1 under eukaryote names) or `listed` (16, kept) |
| Probe mode | `pcr` — MIDORI2 holds whole barcodes and whole genes |
| Preparation | download and unzip, then a one-time conversion of the headers (below) |
| Classifier | SINTAX only (`amplicon.taxonomy.method: sintax`) |
| Region extractor | none |
| Rank model | `Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species` |
| COI screens | reading-frame check and off-target guard, both on by default (`amplicon.coi`) |
| Primers | none are shipped; validated with BF3/BR2, Leray-XT and fwhF2/fwhR2n (below) |
| Memory and time | SINTAX against MIDORI2 holds 8.4 GB at the peak and takes 1 to 1.5 minutes; the COI screens up to 8.6 GB and about a minute (measured on 8 threads; the resources table of the [configuration reference](../../reference/configuration.md) has the details) |

## Setting up a run

```yaml
mode: amplicon

amplicon:
  type: COI
  primers:
    fwd: /path/to/BF3.fasta           # CCHGAYATRGCHTTYCCHCG
    rev: /path/to/BR2.fasta           # TCDGGRTGNCCRAARAAYCA
    orientation: fixed                # fixed | mixed
  expected_length: auto
  probe_length_stat:
    COI: p95
  extraction:
    enabled: false                    # no COI region extractor
  taxonomy:
    method: sintax                    # COI has no rdp trainset
    filter:
      enabled: true
      keep:    [k__Metazoa]           # animal surveys; leave empty to see everything
      discard: []
  coi:
    frame_check: true
    length_tolerance: 9
    offtarget_guard: true
```

`method: rdp` stops the run at start-up: DADA2's `assignTaxonomy` would need about
225 GB of memory to classify against MIDORI2 NUC_SP down to species (about 64 GB for the
smaller NUC variant), so MetaFlux ships no trainset for COI.

!!! warning "Check which read carries which primer"

    In some COI libraries every R1 starts with the reverse primer; the BF3/BR2 and
    fwhF2/fwhR2n runs of the insect mock used for validation (below) are two of them.
    With `orientation: fixed` such a library loses most of its reads at primer trimming,
    but not all. Cutadapt counts a 5′ primer as found when a read starts with the
    primer's last three bases, and the end of each BF3/BR2 primer resembles the start of
    the other, so a `fixed` run on the BF3/BR2 files as distributed keeps about 13 % of
    the pairs, the wrong way round and with most of a primer still on (0.14 % with
    fwhF2/fwhR2n).

    For a library in which every pair is swapped, swap the two files of each sample
    (R1 ↔ R2) and keep `orientation: fixed`. That is the recommended way: Cutadapt then
    reads every pair once, and its report in MultiQC describes all the reads.
    `orientation: mixed` is for libraries that really are mixed, with pairs in both
    orientations, and from v2.6.0 on it handles them: it swaps the two mates of each
    reversed pair back, so that DADA2 merges them to the full amplicon. On the BF3/BR2
    library run both ways, `mixed` on the files as distributed and `fixed` on the swapped
    files gave the same 595 kept ASVs and the same 361 specimens; `mixed` dropped 509
    (0.06 %) of the 910,014 pairs that `fixed` keeps after primer trimming, those whose
    primer Cutadapt found from fewer than 10 bases. See
    [the mixed-orientation second pass](../overview.md#the-mixed-orientation-second-pass).

## The reference database: MIDORI2

MIDORI2 (Leray et al. 2022) collects the mitochondrial protein-coding and rRNA genes of
eukaryotes from GenBank, release by release, with taxonomy from NCBI. MetaFlux uses its COI
file for SINTAX in the UNIQ NUC_SP variant: every unique sequence of each species,
including the records whose species name is incomplete ("sp.", "cf.", BOLD codes). Those
records add about 1.35 million barcodes, mostly of specimens that were barcoded but not
named, and they cover lineages that the named records do not.

| | |
|---|---|
| Downloaded from | `https://zenodo.org/records/22685701/files/MIDORI2_GB273_SINTAX_sp.zip` (584 MB; member `SINTAX_sp/uniq/MIDORI2_UNIQ_NUC_SP_GB273_CO1_SINTAX.fasta.gz`), and `MIDORI2_GB273_KRAKEN2_taxonomy.zip` (69 MB) from the same record |
| Checked | MD5 of each zip as Zenodo lists it; SHA-256 of the extracted FASTA as MIDORI2's `release_sha256.txt` lists it |
| Fetch rules | `fetch_midori2_coi`, `fetch_midori2_taxonomy`, `fetch_bacterial_cox1` |
| Conversion rule | `prepare_midori2_coi_sintax`. It writes no log: its counts go to the report next to the converted file (next row), and when it stops, its message is printed to stderr with Snakemake's own messages |
| Cached at | `refdb/midori2/GB273/midori2_gb273_uniq_nuc_sp_co1_metaflux_sintax.fasta.gz`, with a report of every count the conversion changed next to it (`…_metaflux_sintax.report.tsv`) |

### What the conversion changes, and why

A MIDORI2 SINTAX header reads

```text
>KU746593.1.<1.>658;tax=k:Eukaryota_2759,p:Arthropoda_6656,c:Insecta_50557,o:Diptera_7147,f:Sarcophagidae_7381,g:Sarcophaga_7384,s:Sarcophaga munronis_1806722;
```

and after conversion

```text
>KU746593.1.<1.>658;tax=d:Eukaryota,k:Metazoa,p:Arthropoda,c:Insecta,o:Diptera,f:Sarcophagidae,g:Sarcophaga,s:Sarcophaga_munronis;
```

- **The domain moves to `d:`.** MIDORI2 writes Eukaryota under the letter `k:`, but
  VSEARCH reads the letters in a fixed order (`d` above `k` above `p`) and MetaFlux writes
  `d__` for a domain and `k__` for a kingdom.
- **A kingdom is added** from MIDORI2's own NCBI taxonomy snapshot: the first ancestor of
  rank "kingdom" of the record's species (Metazoa, Fungi or Viridiplantae). Protists have
  no kingdom in NCBI, so their records get no `k:` field. With the kingdom in place,
  `keep: [k__Metazoa]` keeps the animals of a run.
- **The NCBI taxid is cut off** every name (`Arthropoda_6656` → `Arthropoda`), so that
  names and filter tokens read as names.
- **Placeholder ranks are left out.** Where NCBI has no name at a rank, MIDORI2 fills the
  slot with a name built from the ranks below (`genus_Lepismatidae sp. …`). SINTAX would
  report such a name as a taxon. In the NUC_SP variant, half of the records added over NUC
  carry a placeholder genus (287,117 placeholder genera against 58,627 real ones in GB273).
- **Interim species names are left out of the species rank** ("Sarcophaga sp. BIOUG…",
  "cf.", "aff."), so that a species column holds only formal binomials. The record keeps
  its genus and higher ranks.
- **Bacterial COX1 filed under eukaryote names is removed.** MIDORI2 is built from the
  eukaryotic divisions of GenBank, so a bacterial sequence submitted under an animal's name
  enters it. The list that comes with the bacterial COX1 reference set marks 170 such records
  in GB273 for removal, mostly *Wolbachia*, *Tisiphia* and other symbionts filed as insects.
  How they were found is described on [Bacterial COX1 reference set (COI)](bacterial-cox1.md).

The conversion stops, rather than guessing, on a header without `;tax=`, a first field
that is not `k:Eukaryota_<taxid>`, a taxid missing from the taxonomy snapshot, or a
removal list that names records absent from the FASTA (a list made for another MIDORI2
release).

### The reference check

`reference_values` in the COI pack allows only Eukaryota in the domain rank and only
Metazoa, Fungi and Viridiplantae in the kingdom rank. A phylum in the kingdom rank would
mean that MIDORI2's own file is in use instead of the converted one; the run then stops
before any classification.

## The COI screens

COI is a protein-coding gene, and two kinds of sequence reach a COI ASV table that are not
the mitochondrial COI of a eukaryote: pseudogenes and NUMTs, and bacterial COX1. Rule
`coi_screens` measures both for every ASV and writes `06.taxonomy/coi_screens.tsv`; the
taxonomy step applies the decisions **before** keep/discard, so a pseudogene never reaches
the token lists, and a bacterial ASV that SINTAX called an animal is not kept by
`k__Metazoa`. `06.taxonomy/asv_table_unfiltered.txt` keeps every ASV, with the step that
took it out (`removed_frame`, `moved_bacterial_cox1`), and the read-tracking table has a
column for each screen (`post_frame_check`, `post_bacterial_cox1`).

With a screen switched off, the rule still measures, and `coi_screens.tsv` says what the
screen would have done. Its `verdict` column applies only the screens that are on, the
frame check first (an ASV that fails both is `frame`); for a screen that is off, read
`frame_ok` and `length_ok` (frame check) or compare B with E (off-target guard).

### Reading-frame check (`amplicon.coi.frame_check`)

Every ASV starts right after the forward primer, so the ASVs of a run share one codon
phase. MetaFlux takes as the run's reading frame the frame that has no TAA or TAG, the
stop codons of every animal mitochondrial code, in most reads. The reference length is
the most frequent length, counted by reads, among the ASVs open in that frame. An ASV
passes when

- it has no TAA or TAG in the run's frame, or in the frame that an offset of its start
  gives it, and
- its length is within `length_tolerance` bases (default 9, three codons) of the reference
  length.

The offset matters with primers that match some templates imperfectly. Cutadapt can then
cut a base before or after the end of the primer, the same way for every read of that
template, and the ASV starts one base early or late: it is one base longer or shorter, and
open in the neighbouring frame. On the Leray-XT run of the insect mock (below), whose
primers carry inosines, this happened to 9 % of the reads, among them the only ASVs of 29
specimens.
The frame an offset of *d* bases gives is the run's frame plus *d*, so an ASV one base
short is tested in the frame before the run's. A frameshift inside an ASV, as in a
pseudogene or an indel artefact, leaves stop codons in both frames, one on each side of the
shift.

TGA is not used, because it codes for tryptophan in animal and fungal mitochondria, and
neither are AGA and AGG, stops in the vertebrate code only, so the check works the same for
invertebrates and vertebrates. On MIDORI2 GB273 inserts of the three validated primer
pairs, the frame starting at the second base of the insert is free of TAA and TAG in
99.3–100 % of genuine records (per group).

A NUMT that kept its reading frame passes; the check catches the ones with a stop codon or
a frameshift. Tools that judge NUMTs by abundance and phylogeny, such as metaMATE
(Andújar et al. 2021), are not part of MetaFlux.

### Off-target guard (`amplicon.coi.offtarget_guard`)

Most bacteria that respire oxygen carry COX1 homologues, and common COI primer pairs
amplify some of them, *Wolbachia* and other Rickettsiales of insects above all. MIDORI2
holds eukaryotes only, so such an ASV can only be assigned to a eukaryote, usually an
animal, and often with high SINTAX confidence.

The guard compares every ASV with the cleaned MIDORI2 (identity E) and with the bacterial
COX1 reference set, whole genes and fragments (identity B). An ASV with **B > E** is moved
to `06.taxonomy/bacterial_cox1_asvs.tsv`, with its counts, what SINTAX had called it, B,
E and its closest bacterial gene. Nothing is deleted, so a *Wolbachia* infection stays
visible.

- **Search.** VSEARCH `--usearch_global --id 0.5 --strand plus --maxaccepts 256
  --maxrejects 2048`, eight times VSEARCH's default search depth: with distant sequences
  the best alignment is often not among the first targets VSEARCH tries. In the
  development tests this depth recognised the same share of bacterial amplicons as an
  exhaustive search (94.8 %), and a shallower one (20 accepts, 200 rejects) 2.2 points
  fewer.
- **Identity.** Recounted from the aligned rows over A/C/G/T columns: identical A/C/G/T
  columns divided by A/C/G/T columns plus internal gap columns. VSEARCH's own identity
  counts N and other ambiguity codes as identical, and about 9 % of MIDORI2 records
  contain them. A hit counts only when the columns cover at least 80 % of the ASV.
- **Why a comparison and not an identity floor.** Bacterial amplicons and animal sequences
  whose relatives are missing from MIDORI2 fall in the same range of identities to
  MIDORI2: a floor that catches 95 % of the bacterial amplicons also removes 36–87 % of the
  animal sequences whose family is missing. In the final test of the bacterial set, 936
  bacterial amplicons from genomes not used to build its rules were each compared with their
  own species left out of the set, as for a bacterium whose species GTDB lacks: the
  comparison recognised 98.4 % of them (99.9 % for BF3/BR2, 88.2 % for Leray-XT, 93.5 % for
  fwhF2/fwhR2n), and 96.0 % with the whole genus left out. It moved 0.08 % or less of 10,008
  eukaryotic inserts whose species, genus, family or order was left out of MIDORI2. How the
  set was built and tested is on [Bacterial COX1 reference set (COI)](bacterial-cox1.md).
- **What it cannot catch.** A bacterial lineage that is missing from the set but has a
  record in MIDORI2 under an animal's name, near-identical to the ASV, stays with the
  animal: E is then higher than any B. Step 6 of the set removes such records only for
  the lineages the set contains.

### Species flags

`06.taxonomy/coi_asv_flags.tsv` has one row per ASV of `asv_table.txt`: E, the best
MIDORI2 record, the species of every record tied with it at that identity (joined by a
vertical bar; `no_species` there stands for a tied record without a formal species name),
and two flags that add what the SINTAX confidence cannot show.

| `species_flag` | Meaning |
|---|---|
| `ok` | SINTAX called a species, the best MIDORI2 hit is at least 98 % identical, and no record of another named species is as close |
| `low_identity` | the best MIDORI2 hit is less than 98 % identical: the species, or a close relative, may be missing from MIDORI2 |
| `shared_sequence` | a record of another named species is as close as the best one (sister species with shared haplotypes, synonyms, mislabelled records): this amplicon cannot tell them apart |
| `no_species` | SINTAX called no species |

`far_from_midori2` is 1 when E is below 80 %: an ASV of a group that MIDORI2 covers
poorly, or a bacterium that the bacterial set does not contain. It is a flag only.

Why flags and not a filter: in a leave-out test with MetaFlux's SINTAX settings on 418-bp
insect sequences, a species present in MIDORI2 was called for 91.4 % of the queries, and
correctly in 96 % of those calls; with the query's species removed from MIDORI2, 24.2 %
still received a (wrong) species at confidence 0.8 or more. A 98 % identity limit removed
only 35 % of those wrong calls, because many of them come from near-identical records of
other species, and it also removed 1.6 % of the right ones.

## What the run produces

| File | Contents |
|---|---|
| `refdb/cache/probe_COI_midori2_gb273_uniq_nuc_sp_co1_<hash>.json` | Probed amplicon-length distribution for the primer pair, from MIDORI2 |
| `06.taxonomy/asv_table.txt` | Counts per sample and the MIDORI2 taxonomy string, without the ASVs the screens took out |
| `06.taxonomy/asv_table_unfiltered.txt` | Every ASV that reached the filter, with `filter_status` (`removed_frame`, `moved_bacterial_cox1` and the statuses of the taxon filter) |
| `06.taxonomy/coi_screens.tsv` | Per ASV: its reads summed over all samples, length, open frames, frame and length verdicts, B and the closest bacterial gene, E, the closest MIDORI2 record and the species of the records tied with it, and the verdict of the screens |
| `06.taxonomy/bacterial_cox1_asvs.tsv` | The ASVs the off-target guard moved: counts, SINTAX call, B, E, closest bacterial gene and its GTDB taxonomy, sequence |
| `06.taxonomy/coi_asv_flags.tsv` | Per kept ASV: E, best MIDORI2 record, the species of the records tied with it, `species_flag`, `far_from_midori2` |
| `stats/read_tracking.txt` | With the columns `post_frame_check` and `post_bacterial_cox1` |
| `logs/coi_screens.log` | The run's reading frame and reference length, and the counts of each screen |

## Validation

The three primer-pair runs of the insect mock of Elbrecht et al. (2019) (bulk abdomens of
374 arthropod specimens; SRR8082166 BF3/BR2, SRR8082178 Leray-XT, SRR8082148
fwhF2/fwhR2n) were each run as one sample, with the defaults above and
`keep: [k__Metazoa]`, and their ASVs were compared with the 374 specimen barcodes. In the
BF3/BR2 and fwhF2/fwhR2n libraries
every R1 carries the reverse primer, so their two files were swapped and the runs kept
`orientation: fixed`, as the warning above recommends; the Leray-XT library, whose R1
carries the forward primer, ran with `fixed` as distributed.

| | BF3/BR2 | Leray-XT | fwhF2/fwhR2n |
|---|---|---|---|
| reads after the length filter (ASVs) | 638,878 (804) | 897,924 (767) | 1,199,927 (888) |
| reference length, run frame | 418 bp, frame 1 | 313 bp, frame 1 | 205 bp, frame 1 |
| removed by the frame check | 177 ASVs, 3,716 reads (0.58 %) | 95, 2,160 (0.24 %) | 64, 1,540 (0.13 %) |
| moved by the off-target guard | 27 ASVs, 3,283 reads (0.51 %) | 26, 2,353 (0.26 %) | 16, 895 (0.07 %) |
| removed by `keep: [k__Metazoa]` | 5 ASVs, 299 reads (0.05 %) | 5, 107 (0.01 %) | 7, 80 (0.01 %) |
| kept | 595 ASVs, 631,580 reads | 641, 893,304 | 801, 1,197,412 |
| reads of kept ASVs within 97 % of a specimen | 99.2 % | 99.3 % | 99.1 % |
| specimens recovered at ≥ 97 % (identical) | 361 (357) | 355 (352) | 362 (362) |

- **Nothing a specimen needed was removed.** No specimen that was missed has an ASV within
  97 % among those the two screens took out, so no specimen lost its only ASV to a
  screen. The few removed ASVs that match a specimen (14, 3 and 2) belong to specimens
  that keep another ASV: in BF3/BR2, 11 of them are length variants and 3 are full-length
  sequences with a stop codon in every frame; the 3 and 2 of the other runs pass the
  length check (one Leray-XT ASV is two bases longer than the reference) and carry a stop
  codon in every frame.
- **The offset of the frame check.** In the Leray-XT run, 9 % of the reads, among them
  the only ASVs of 29 specimens, belong to ASVs open in frame 0, most of them one base
  short: their first base went with the primer, which cutadapt cut one base too far for
  these templates. In the fwhF2/fwhR2n run the same happened to 3.5 % of the reads (14
  specimens). Without the offset rule the frame check would have removed them.
- **The off-target guard.** It moved 69 ASVs over the three runs, most closest to
  *Wolbachia* (57), the others to *Rickettsia*, *Mesenet*, *Tisiphia* and one other
  bacterium. SINTAX had called 41 of them, with 4,315 reads, the stink bug *Picromerus
  bidens*, which is not in the mock, at species level. All 41 are closest to a
  *Wolbachia* gene in the bacterial set, and their best MIDORI2 hit is a *Wolbachia*-like
  record filed under that name (KR036430.1, no TGA among five tryptophan codons);
  `keep: [k__Metazoa]` would have kept them.
- **Species flags.** Of the kept ASVs, 53 % (BF3/BR2), 46 % (Leray-XT) and 38 %
  (fwhF2/fwhR2n) are `ok`, and `low_identity` and `shared_sequence` together flag 9–11 %;
  the others have no species call, more of them the shorter the amplicon. None is far from
  MIDORI2.
- **Counts near the detection limit move.** When the same library is run again with
  another setting (split into two samples or run as one, other truncation lengths,
  another error model), the specimens it gains or loses are recovered by a single ASV of
  3 to 10 reads, at the limit where DADA2 tells a rare variant from an error. The
  specimen counts above therefore shift by one or two between such runs, while the shares
  of reads do not. Specimens that one library recovers and another misses are a
  different matter: they are missed with thousands of reads elsewhere (in the BF3/BR2
  library SRR8082166 no kept ASV comes near the barcode of *Camponotus pennsylvanicus*,
  while the Leray-XT library gives it 5,295 reads and a second BF3/BR2 library of the
  same mock, SRR8082198, 633), so the PCR that made the library, not DADA2's detection
  limit, decides them.

## References

- Andújar, C., et al. (2021). Validated removal of nuclear pseudogenes and sequencing
  artefacts from mitochondrial metabarcode data. *Molecular Ecology Resources*
  21:1772–1787. doi:10.1111/1755-0998.13337
- Elbrecht, V., et al. (2019). Validation of COI metabarcoding primers for terrestrial
  arthropods. *PeerJ* 7:e7745. doi:10.7717/peerj.7745
- Leray, M., Knowlton, N., Machida, R. J. (2022). MIDORI2: a collection of quality
  controlled, preformatted, and regularly updated reference databases for taxonomic
  assignment of eukaryotic mitochondrial sequences. *Environmental DNA* 4:894–907.
  doi:10.1002/edn3.303

The full list, including the tools used at every step, is on
[Citation and references](../../about/citation.md).
