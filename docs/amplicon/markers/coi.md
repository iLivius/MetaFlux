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
| Second reference | the bacterial COX1 reference set `bacterial_cox1_gtdb-r232.v1`, built from GTDB R232: identity B of the off-target guard, and the list of MIDORI2 records that are bacterial COX1 under eukaryote names |
| Probe mode | `pcr` — MIDORI2 holds whole barcodes and whole genes |
| Preparation | download and unzip, then a one-time conversion of the headers (below) |
| Classifier | SINTAX only (`amplicon.taxonomy.method: sintax`) |
| Region extractor | none |
| Rank model | `Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species` |
| COI screens | reading-frame check and off-target guard, both on by default (`amplicon.coi`) |
| Primers | none are shipped; validated with BF3/BR2, Leray-XT and fwhF2/fwhR2n (below) |
| Memory and time | SINTAX against MIDORI2 holds about 8.3 GB (8.5 minutes for 717 ASVs on 8 threads); the COI screens about 8.4 GB (about a minute) |

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

    In some COI libraries every R1 starts with the reverse primer; the BF3/BR2 run of the
    insect mock used for validation (below) is one. With `orientation: fixed` such a
    library loses all its reads at primer trimming. Swap the two files of each sample
    (R1 ↔ R2) and keep `orientation: fixed`. [To be settled before release:
    `orientation: mixed` reverse-complements both mates of a swapped pair instead of
    swapping them, and for an amplicon longer than the overlap of the two reads, such as
    the 418-bp BF3/BR2 insert, DADA2 then keeps only that overlap.]

## The reference database: MIDORI2

MIDORI2 (Leray et al. 2022) collects the mitochondrial protein-coding and rRNA genes of
eukaryotes from GenBank, release by release, with taxonomy from NCBI. MetaFlux uses its COI
file for SINTAX in the UNIQ NUC_SP variant: every unique sequence of each species,
including the records whose species name is incomplete ("sp.", "cf.", BOLD codes). Those
records add about 1.35 million barcodes, mostly of specimens that were barcoded but not
named, and they cover lineages that the named records do not.

| | |
|---|---|
| Downloaded from | `https://zenodo.org/records/22685701/files/MIDORI2_GB273_SINTAX_sp.zip` (member `SINTAX_sp/uniq/MIDORI2_UNIQ_NUC_SP_GB273_CO1_SINTAX.fasta.gz`), and `MIDORI2_GB273_KRAKEN2_taxonomy.zip` from the same record |
| Checked | MD5 of each zip as Zenodo lists it; SHA-256 of the extracted FASTA as MIDORI2's `release_sha256.txt` lists it |
| Fetch rules | `fetch_midori2_coi`, `fetch_midori2_taxonomy`, `fetch_bacterial_cox1` |
| Conversion rule | `prepare_midori2_coi_sintax` → `<out_dir>/logs/refdb/prepare_midori2_coi_sintax.log` |
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
  enters it. The bacterial COX1 reference set lists 170 such records in GB273, mostly
  *Wolbachia*, *Tisiphia* and other symbionts filed as insects. How they were found is
  described on [Bacterial COX1 reference set (COI)](bacterial-cox1.md).

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
screen would have done.

### Reading-frame check (`amplicon.coi.frame_check`)

Every ASV starts right after the forward primer, so all ASVs of a run share one codon
phase. MetaFlux takes as the run's reading frame the frame free of TAA and TAG, the stop
codons of every animal mitochondrial code, in most reads, and as the reference length the
most frequent length, by reads, of the ASVs open in that frame. An ASV passes when

- it has no TAA or TAG in the run's frame, and
- its length differs from the reference length by whole codons only, by at most
  `length_tolerance` bases (default 9, three codons).

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
  animal sequences whose family is missing. In the final test of the bacterial set, the
  comparison recognised 98.4 % of 936 bacterial amplicons from genomes not used to build its
  rules (99.9 % for BF3/BR2, 88.2 % for Leray-XT, 93.5 % for fwhF2/fwhR2n), 96.0 % when the
  amplicon's whole genus was left out of the set, and moved 0.08 % or less of 10,008
  eukaryotic inserts whose species, genus, family or order was left out of MIDORI2. How the
  set was built and tested is on [Bacterial COX1 reference set (COI)](bacterial-cox1.md).
- **What it cannot catch.** A bacterial lineage that is missing from the set but has a
  record in MIDORI2 under an animal's name, near-identical to the ASV, stays with the
  animal: E is then higher than any B. Step 6 of the set removes such records only for
  the lineages the set contains.

### Species flags

`06.taxonomy/coi_asv_flags.tsv` has one row per ASV of `asv_table.txt`: E, the best
MIDORI2 record and its species, and two flags that add what the SINTAX confidence cannot
show.

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
| `06.taxonomy/coi_screens.tsv` | Per ASV: length, open frames, frame and length verdicts, B and the closest bacterial gene, E and the closest MIDORI2 record and species, and the verdict of the screens |
| `06.taxonomy/bacterial_cox1_asvs.tsv` | The ASVs the off-target guard moved: counts, SINTAX call, B, E, closest bacterial gene and its GTDB taxonomy, sequence |
| `06.taxonomy/coi_asv_flags.tsv` | Per kept ASV: E, best MIDORI2 record and species, `species_flag`, `far_from_midori2` |
| `stats/read_tracking.txt` | With the columns `post_frame_check` and `post_bacterial_cox1` |
| `logs/coi_screens.log` | The run's reading frame and reference length, and the counts of each screen |

## Validation

The BF3/BR2 run of the insect mock of Elbrecht et al. (2019) (SRR8082166: bulk abdomens of
374 arthropod BINs, 918,986 read pairs) was run with the defaults above and
`keep: [k__Metazoa]`, and its ASVs were compared with the 374 specimen barcodes. The library
was split into two pseudo-samples of alternate read pairs, because MetaFlux stops on a run
with a single sample (`dada_seqtab`); the R1 and R2 files were swapped (see the warning
above).

- **Reading frame.** 99.5 % of the reads after the length filter are open in frame 1 (the
  second base after the primer), and the most frequent length among them is 418 bp, the
  BF3/BR2 insert.
- **Frame check.** It removed 165 of the 717 ASVs, with 3,791 reads (0.6 %), and none of
  them was the only ASV of a mock specimen: the 14 that match a specimen at 97 % or more are
  length or frame variants of a kept ASV.
- **Off-target guard.** It moved 22 ASVs, with 3,073 reads (0.5 %), to
  `bacterial_cox1_asvs.tsv`: 18 closest to *Wolbachia*, 2 to *Mesenet*, 1 to *Rickettsia*,
  1 to another bacterium. SINTAX had called 12 of them, with 1,879 reads, the stink bug
  *Picromerus bidens*, which is not in the mock, at species level: MIDORI2 holds a
  *Wolbachia*-like record under that name (KR036430.1, no TGA among five tryptophan codons),
  and `keep: [k__Metazoa]` would have kept them.
- **What is left.** 526 ASVs with 631,678 reads. 423 of them, with 99.3 % of the reads, are
  at least 97 % identical to a specimen barcode, and 359 of the 374 specimens are recovered
  at 97 % or more (354 of them identically). The 15 others have no ASV within 97 %.
- **Species flags.** Of the 526 kept ASVs, 294 are `ok`, 192 `no_species`, 23
  `shared_sequence` and 17 `low_identity`; none is far from MIDORI2.

The Leray-XT and fwhF2/fwhR2n runs of the same mock (SRR8082178, SRR8082148) have not been
run through MetaFlux yet; for those pairs the off-target guard is so far tested in silico
only (see [Bacterial COX1 reference set (COI)](bacterial-cox1.md#step-7-final-test)).

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
