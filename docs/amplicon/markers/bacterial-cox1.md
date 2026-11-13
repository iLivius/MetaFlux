# Bacterial COI homologues (reference set)

COI primers do not amplify only animals. Most bacteria that respire oxygen carry genes
homologous to the mitochondrial cytochrome c oxidase subunit I, and common COI primer pairs
amplify some of them, those of *Wolbachia* above all. MIDORI2, the COI reference of MetaFlux,
contains only records filed under eukaryotes, so an amplicon of bacterial origin cannot be
recognised as such: it can only be assigned to a eukaryote, usually an animal.

MetaFlux therefore uses a second reference: a set of bacterial and archaeal subunit-I genes
taken from the species representatives of the Genome Taxonomy Database, GTDB (Parks et al.
2026). For every ASV, MetaFlux compares its best identity to this set (B) with its best
identity to MIDORI2 (E). An ASV with B > E is reported in a separate table of bacterial COI
homologues. The same set is used once, when MIDORI2 is prepared, to remove MIDORI2 records
that are themselves bacterial COI homologues filed under eukaryotic names.

!!! note "Names"
    COI, or COX1, names the mitochondrial gene. The genes of this set are its bacterial and
    archaeal homologues: subunit I of the heme-copper oxidases (genes *coxA* or *ctaD*, *cyoB*,
    *qoxB* and *ccoN*, among others) and of the related nitric oxide reductases, which the build
    does not separate (step 3). This manual calls them bacterial COI homologues; "bacterial"
    includes the archaea. File names, settings and output columns keep the short form `cox1`
    (for example `bacterial_cox1_gtdb-r232.v2` and `bacterial_cox1_asvs.tsv`).

!!! note "You do not build this set"
    MetaFlux downloads the released set and checks its SHA-256 sum. This page explains how
    the set was built, why each step is done as it is, how other tools deal with the same
    problem, and what the set cannot do.

All rules on this page were set on a **development sample** of 3,785 species representatives
of GTDB R226, the release before R232: all 485 species of the order Rickettsiales (56 of them
*Wolbachia*), which holds *Wolbachia*, *Rickettsia* and other symbionts of arthropods, plus
3,000 other bacteria and 300 archaea drawn per phylum (at least three per phylum, the rest in
proportion to the number of species in the phylum). Shares given for "all bacteria" are
weighted back to the number of GTDB species in each group. Some rules of step 6 and check 5.4
were set after looking at the R232 genes and the MIDORI2 GB273 records. The set
was then tested on data not used to set the rules (step 7).

## At a glance

| | |
|---|---|
| Release | `bacterial_cox1_gtdb-r232.v2` (Zenodo DOI: [10.5281/zenodo.24189850](https://doi.org/10.5281/zenodo.24189850); all versions: [10.5281/zenodo.23820759](https://doi.org/10.5281/zenodo.23820759)) |
| Source | GTDB R232: 199,923 species representatives, 189,801 bacterial and 10,122 archaeal |
| Genes | GTDB's own gene predictions, made with Prodigal 2.6.3 (Hyatt et al. 2010), file `gtdb_proteins_nt_reps_r232.tar.gz`; the build runs no gene prediction of its own |
| Content | subunit-I genes of the heme-copper oxidase superfamily, as DNA: whole genes in one file, fragments of at least 200 amino acids in a second file |
| Size | 242,696 whole genes (121 MB compressed) and 31,709 fragments (11 MB) from 123,358 genomes; counts per step in the release README |
| MIDORI2 list | 500 MIDORI2 GB273 records, one row per record, keyed by the SHA-1 of the sequence: 236 records (187 sequences) `removed`, 264 `listed` |
| Taxonomy in the headers | GTDB R232, in SINTAX format (`d:`, `p:`, `c:`, `o:`, `f:`, `g:`, `s:`) |
| Used by MetaFlux | the B-versus-E comparison for every ASV, and the cleaning of MIDORI2 (MetaFlux 2.7.0 and later) |
| Licence | CC BY-SA 4.0, as GTDB |
| Tools for the build | HMMER 3.4 (Eddy 2011) and VSEARCH 2.32.0 (Rognes et al. 2016) from bioconda; Python 3.12.15, numpy 1.26.4 and pigz 2.8 from conda-forge; cutadapt 4.9 for the final test |

## How MetaFlux uses the set

For a COI run, MetaFlux downloads the release and MIDORI2, removes from MIDORI2 every record
whose sequence step 6 marks as removed, and computes for every ASV its best identity to the
cleaned MIDORI2 (E) and to the set (B), as described in step 5, over at least 80 % of the ASV.

| ASV | What MetaFlux does |
|---|---|
| B > E | moved from the main ASV and taxonomy tables to `bacterial_cox1_asvs.tsv`, with its read counts, B, E and its best bacterial match; nothing is deleted, so that, for example, a *Wolbachia* infection stays visible |
| all other ASVs | E and the best MIDORI2 match are reported; ASVs with E below 80 % are flagged |

In the read-tracking table, `post_bacterial_cox1` holds the reads left after the move; the
moved reads are `post_frame_check` minus `post_bacterial_cox1`. The 80 % flag marks ASVs far
from everything in MIDORI2: in the development tests it marked 64–90 % of bacterial amplicons,
depending on the primer pair, and at most 5 % of genuine animal sequences whose genus was
missing from MIDORI2. Removing flagged ASVs would cost too many genuine eukaryotes (see "What
MetaFlux does, and why"), so they are only flagged.

The list of MIDORI2 records belongs to MIDORI2 GB273. Each row carries the SHA-1 of the
record's sequence (column `sequence_sha1`: upper case, line breaks removed) and the SHA-256 of
the MIDORI2 file the list was made from (column `midori2_file_sha256`). MetaFlux removes every
record whose sequence is marked removed, under whatever name MIDORI2 files it; on GB273 this
leaves out 236 of the 3,213,682 records (187 distinct sequences), the records the build
removed. The MIDORI2 conversion stops on a list without `sequence_sha1`,
and on a list made from the same MIDORI2 file whose removed sequences are not all in it. With
a list made from another MIDORI2 file, the missing removed sequences are only counted and the
records the new file adds are not examined, so a new MIDORI2 release needs a new list, made
with step 6 of the build scripts. Settings are on the [COI marker page](coi.md).

## Overview of the build

```text
0  source: GTDB R232 gene file (DNA, all species representatives) and metadata
   (taxonomy, CheckM2, GC content, protein counts), checked against GTDB's MD5 sums
      │
1  per genome: number of genes = GTDB's protein count; translate with genetic code 11,
   or code 4 where the genome's genes show TGA read as tryptophan
      │
2  HMM search of all proteins with the Pfam COX1 model (PF00115, gathering threshold);
   COX2, COX3 and cytochrome b models searched as well, to describe each contig
      │
3  per hit, from residues at fixed alignment positions: is it an oxidase subunit I
   (six metal-binding histidines), and is it family A (subtype A1 or A2)?
   type label (aa3, bo3, aa3-600, cbb3) from TIGRFAM and NCBIfam models where they apply
      │
4  whole gene or fragment: start and stop codon present, ≥ 70 % of the model covered,
   400–1,100 amino acids; otherwise a fragment if ≥ 200 amino acids
      │
5  contamination checks (a failure removes the gene):
     5.1 reading frame   5.2 GC content   5.3 match to mitochondrial COI
     5.4 closer to MIDORI2 (cleaned in step 6) than to genes of the same GTDB family
   flags only: cross-phylum cluster; mitochondrion-like contig; duplicate copy
      │
6  cleaning of MIDORI2, once per sequence: records close to a gene are removed when
   the evidence says "bacterial" and nothing says "mitochondrial"; also isolated
   records and chimeras; others are listed; repeated with 5.4 until nothing changes
      │
7  final test on data not used to set the rules, and regression checks
      │
8  release: genes and fragments (FASTA), candidate table, MIDORI2 list, README,
   scripts, SHA-256 sums  ──►  Zenodo
```

## Why a bacterial reference is needed

### Bacterial COI homologues in COI libraries

Subunit I of the heme-copper oxygen reductases is the catalytic subunit of the terminal
oxidases of aerobic respiration. The enzymes form three families, A, B and C. Family A
contains the mitochondrial enzyme and the aa3-type oxidases of many bacteria, and is split
into subtypes A1 and A2 (Pereira et al. 2001). Heme-copper oxidases are widely distributed in
Bacteria and Archaea (Sousa et al. 2012), and many genomes carry several of them. The
bacterial genes are known as *ctaD* or *coxA* (aa3-type cytochrome c oxidase), *cyoB*
(bo3-type quinol oxidase), *qoxB* (aa3-600 quinol oxidase) and *ccoN* (cbb3-type oxidase,
family C). COI is not a general marker for bacteria; *Wolbachia* is the exception: *coxA* is
one of the five loci of its multilocus sequence typing scheme (Baldo et al. 2006), and Nowak
et al. (2026) use insect COI barcode data to survey *Wolbachia* strains.

COI primers bind the bacterial genes because the primer sites lie in conserved parts of the
protein. The table gives the share of genomes with a family-A (aa3-like) gene in which each
primer matches perfectly.

| Primer | Nowak et al. (2026, preprint): 59,318 RefSeq genomes, 36.8 % with an aa3-like gene | Development sample (weighted to GTDB) |
|---|---|---|
| jgHCO2198 (reverse, Leray pairs) | 79.8 % | 78.0 % |
| BR2 (reverse) | 70.6 % | 64.6 % |
| BF3 (forward) | 1.9 % | 1.4 % |

The reverse primers bind widely; the forward primer decides which bacteria amplify. A product
forms only when both primers bind, so an in-silico PCR on the development sample counted a
genome as amplified when both primers bind its DNA, facing each other, with the expected
insert (418, 313 and 205 bp for the three pairs, ± 30 bp), at most one mismatch per primer,
none in the last two bases at the 3′ end, and no insertions or deletions. The tests of B > E
below use every product with up to three mismatches, to include more distant amplicons.

| Primer pair | *Wolbachia* (56 genomes) | other Rickettsiales (429 genomes) | all bacteria (weighted to GTDB) |
|---|---|---|---|
| BF3/BR2 | 93 % | 25 % | 3.2 % |
| Leray-XT | 70 % | 15 % | 0.2 % |
| fwhF2/fwhR2n | 0 % | 0.5 % | 0 % |

Each value is the share of genomes in the group that give a PCR product. BF3/BR2, for example,
gives a product in 52 of the 56 *Wolbachia* genomes, but in only about 3 of every 100
bacterial species. How much of a library is bacterial depends on the sample:

| Sample | Bacterial share | Source |
|---|---|---|
| water eDNA, UK estuarine and coastal sites | 45.8 % of Leray-XT reads, 55 % of the reads of a fish-targeted, less degenerate pair, had a bacterial best hit in NCBI nt | Collins et al. 2019 |
| single insects, 7,901 BF3/BR2 libraries | *Wolbachia* COI in 35.1 % of libraries, median 1.8 % of the reads where present | Nowak et al. 2026 |
| insect specimens, 1.09 million Sanger barcode traces | *Wolbachia* in 0.16 % | Smith et al. 2012 |

| Safeguard of a eukaryote-only reference | Why it fails for bacterial sequences |
|---|---|
| top rank | MIDORI2 has Eukaryota in every record, so SINTAX (Edgar 2016) reports it with support 1.00 for any sequence |
| Metazoa filter | with a kingdom rank from MIDORI2's taxonomy files, SINTAX assigned 42–69 % of the bacterial in-silico PCR amplicons, and 94–100 % of the *Wolbachia* amplicons, to Metazoa with confidence ≥ 0.8 |
| pseudogene filter | it rejects TAA or TAG in frame (the invertebrate mitochondrial stops), which a gene of a table-11 bacterium, as almost all are, does not contain |
| clean reference | MIDORI2 is built from the eukaryotic divisions of GenBank (Leray et al. 2022), so bacterial sequences submitted under a eukaryote's name enter it: with the development genes alone, 132 records of GB273 were ≥ 90 % identical to a bacterial gene |

With a reference made only of eukaryotic records, none of these safeguards recognises a
bacterial sequence.

### How other tools handle it

COI tools handle non-target sequences in two ways: bacteria in the reference, or a
eukaryote-oriented reference with similarity limits or taxonomic filters. The table gives the
approach and the bacterial content or limits of each tool.

| Tool | Approach | Bacterial content, limits and handling |
|---|---|---|
| CO1 Classifier (Porter & Hajibabaei 2018), used with MetaWorks (Porter & Hajibabaei 2022) | bacteria in the reference | "bacterial outgroup sequences were added", selection not documented; version 5.1.0: 262 bacterial of 2,216,547 records, 229 of them *Wolbachia* barcode-region records; RDP classifier (Wang et al. 2007), since 5.0.0 also SINTAX and BLAST; with pseudogene filtering on, MetaWorks keeps only ESVs (its term for ASVs) assigned to the target group (Arthropoda in the shipped configuration) |
| COInr (Meglécz 2023) | bacteria in the reference | all cellular organisms, origin taken at face value; release of 3 May 2026: 9,775 bacterial sequences from 970 genera, 9,454 from BOLD; VTAM and mkLTG point to files of the release of 6 May 2022, not counted here |
| MJOLNIR (db_COI_MBPK, Wangensteen et al. 2018; DUFA) | bacteria in the reference | some prokaryotes, mostly Rickettsiales; MJOLNIR3 removes by default the MOTUs (clusters of similar sequences) assigned to prokaryotes, the original MJOLNIR does not |
| DARN (Zafeiropoulos et al. 2021) | placement on a reference tree | 463 bacterial and 21 archaeal consensus sequences, built from 3,917 BOLD and 9,154 genome-derived (PF00115) bacterial sequences and 217 archaeal ones; in the authors' mock test BLAST assigned all bacterial sequences to eukaryotes, placement recognised all |
| study-specific MIDORI (Nowak et al. 2026) | bacteria in the reference | alphaproteobacterial sequences added to MIDORI release GB239 |
| BOLD | contaminant check | up to version 3, uploads were checked against contaminants including *Wolbachia* (Smith et al. 2012); no description found for version 5 |
| VTAM (González et al. 2023) | identity steps | lowest group holding ≥ 90 % of the BLAST hits, from 100 % down to 70 %; groups defined at ≤ 80 % "not very robust" |
| eKOI (González-Miguéns et al. 2025) | identity limit | ASVs below 84 % ignored in the re-analysis of 15 datasets, to be applied with caution |
| Anacapa BLCA (Curd et al. 2019) | identity limit | hits of at least 80 % |
| BOLDigger3 3.0.3 (Buchner & Leese 2020) | rank limited by identity | 97, 95, 90, 85 and 75 % from species to class; 50 % for phylum |
| apscale-blast 2.0.6 | rank limited by identity | 98, 95, 90, 87 and 85 % from species to class |
| nf-core/ampliseq 3.0.0 (Straub et al. 2020) | identity limit and filters | MIDORI2 CO1 for DADA2 and a VSEARCH consensus classifier (90 %, full query coverage); optional codon filter; name-based exclusion of mitochondria and chloroplasts |
| BGE Bioscan workflow (BF3/BR2) | taxonomic filter | SINTAX against BOLDistilled; keeps ASVs assigned to Animalia, and to a class, at confidence ≥ 0.8 |

With a eukaryote-oriented reference and an identity limit or a taxonomic filter (VTAM to the
BGE Bioscan workflow in the table), a bacterial ASV below the lowest limit stays unassigned;
above it, it is assigned to a eukaryotic group at the rank its identity allows; where only the
target group is kept, it is removed unless assigned to that group. No published measurement of
how often each case occurs for bacterial ASVs was found, nor a study comparing the approaches.

### What MetaFlux does, and why

For every ASV, MetaFlux measures two identities: B, to the closest gene of the bacterial set,
and E, to the closest record of MIDORI2. The ASV counts as bacterial when B > E; there is no
fixed identity limit.

The rule was chosen on test sequences whose origin is known, and on the hard case of real data:
sequences whose close relatives are missing from the references.

- **Bacterial test sequences:** the COI amplicons the three primer pairs would produce from the
  development genomes (in-silico PCR). The genomes were split into two halves by species, and
  each half served as the bacterial reference for the other, so no amplicon could find its own
  species; the test was repeated with the split by genus.
- **Eukaryotic test sequences:** real MIDORI2 COI records, each compared with MIDORI2 after
  every record of its own species was removed, then of its genus, family and order.

The final test (step 7) repeats these measurements on genomes and records not used to set the
rule.

| Primer pair | Identity limit needed to catch 95 % of bacterial amplicons | Genuine animal sequences lost at that limit (family missing from MIDORI2) | Bacterial amplicons recognised by B > E |
|---|---|---|---|
| BF3/BR2 | 83 % | 36 % | 96 % |
| Leray-XT | 89 % | 72 % | 91 % |
| fwhF2/fwhR2n | 95 % | 87 % | 85 % |

The first three columns show that no identity limit separates bacteria from eukaryotes:
bacterial amplicons and eukaryotes whose relatives are missing fall in the same range of
identities to MIDORI2. The last column shows the comparison B > E, which recognised 94 % of the
bacterial amplicons pooled and 97 % of the *Wolbachia* amplicons. Of the genuine animal
sequences, it moved none when their species or genus was missing from MIDORI2, 0.04 % when
their family and 0.14 % when their order was missing. Two simpler variants were weaker or added
nothing:

- A one-sided screen removing any ASV at least 78 % identical to a bacterial amplicon, the
  lowest limit that removed no eukaryote (highest eukaryotic identity 77.6 %, once two genomes
  with eukaryotic contamination were left out), recognised 71 % of the bacterial amplicons.
- Bacterial sequences added to the SINTAX reference pulled none of 4,560 eukaryotic test
  sequences into Bacteria, but added at most 0.6 percentage points to the share recognised by
  B > E. MetaFlux does not add them.

## The build, step by step

### Step 0. Source of the bacterial sequences

| Source | For | Against |
|---|---|---|
| bacterial records labelled COI in GenBank or BOLD (COInr, CO1 Classifier) | ready to use | labels from submitters, some wrong (two BOLD "Alphaproteobacteria" are fly COI: 94–97 % to fly records, GC 31–37 %); few lineages (*Wolbachia* 229 of 262 in the CO1 Classifier; Actinomycetes and Alphaproteobacteria three quarters of COInr's) |
| in-silico PCR products cut from genomes | close to what the primers amplify | tied to the primer pairs of the build |
| whole genes predicted in genomes (used in part by DARN) | every subunit-I gene, independent of primers, with the genome's taxonomy | MAGs cut genes at contig ends and can contain foreign contigs |

Each row is a source considered, with its advantage and its drawback. MetaFlux uses the
third, from GTDB: one taxonomy, a representative per species, and GTDB's own gene predictions,
used as they are. Only species representatives are used:
with the species left out of the bacterial reference, the development tests still recognised
94 % of the bacterial amplicons (78 % with the genus left out), so the other 700,000 genomes of
the release would multiply the size of the set with little gain. Step 4 deals with cut genes
in MAGs (metagenome-assembled genomes), step 5 with contamination.

### Step 1. Integrity and translation

The gene file and the metadata are checked against GTDB's MD5 sums, and the number of genes of
each genome must equal its protein count in the metadata (a genome that fails is left out and
listed). The build translates the DNA genes itself rather than downloading GTDB's protein file
(another 123 GB), so protein and DNA always come from the same record. The genetic code must
then be worked out per genome: GTDB predicts some lineages, mainly Mycoplasmatales (GTDB phylum
Bacillota_I), with table 4 (TGA codes for tryptophan), but the table it used is not in the
release files, and the metadata give NCBI's table only, none for 92,894 of the 189,801
bacterial representatives.

Under table 11 no predicted gene can contain TGA in frame, so a genome is translated with
table 4 when more than 1 % of its genes with a start and a stop codon (`partial=00`) contain an
in-frame TGA before the stop codon. In R232 this assigned 918 of the 199,923 genomes to table
4; for 137 genomes the inferred table differs from NCBI's. "Mollicutes use table 4" would have
been wrong for 11 development genomes, and the coding-density rule of GTDB-Tk disagreed with
GTDB's table for 2 of 3,785. None of the 15 development genomes predicted with table 4 carried
a candidate gene.

### Step 2. Finding subunit-I genes

All proteins are searched with HMMER `hmmsearch` and the Pfam model PF00115 (Pfam 38.2;
Paysan-Lafosse et al. 2025) at its gathering threshold (34.8 bits), the score above which Pfam's
curators count a protein as a family member; PF00115 covers subunit I of all three oxidase
families and the nitric oxide reductases. COX2, COX3 and cytochrome b models are searched only
to describe each contig (step 5). The table gives the alternatives tested and not used.

| Alternative | Why not used for selection |
|---|---|
| family-A motif `W.[FW]GHPEV[YF]....P` | finds only subtype A1; missed the 476 A2-type genes of the development sample, among them every cyanobacterial aa3 oxidase |
| TIGRFAM models per oxidase type, at their trusted cutoffs | too strict: left 758 whole family-A genes without a label, holding 143 of the 1,239 in-silico PCR products in kept genes; used only for labels (step 3) |
| RefSeq gene annotations | depend on the pipeline: NCBI names its copy of PF00115 "cbb3-type cytochrome c oxidase subunit I", so family-B genes carry that name |

### Step 3. Oxidase or not, and which family

A protein can match PF00115 without being an oxidase subunit I: in the development sample,
249 of the 4,304 genes that passed all other checks (5.8 %) failed the oxidase rule below (193
with two or more ligand histidines replaced, 56 below 100 bits), and no in-silico PCR product
lay in them. Each hit is aligned to PF00115 together with human COX1 (GenBank HQ713447.1),
which gives the alignment columns their usual residue numbers (identical to bovine numbering
at these positions).

| Bovine position | Role | Used for |
|---|---|---|
| His61, His378 | ligands of heme a | oxidase |
| His376 | ligand of heme a3 | oxidase |
| His240, His290, His291 | ligands of CuB | oxidase |
| Tyr244 | cross-linked to His240 (helix VI) | family A |
| Glu242 | proton channel, motif `XGHPEV` | subtype A1 |
| Tyr238, Ser239 | motif `YSHPXV`, replacing the A1 glutamate | subtype A2 |

An oxidase scores at least 100 bits and has at most one histidine replaced (a gap is not a
replacement), because one alignment can misplace a residue in a divergent sequence. Family A
and its subtypes are defined by the residues in the table (Pereira et al. 2001). The type label
comes from TIGRFAM and NCBIfam models of NCBI hmm_PGAP 20.0 (Li et al. 2021; parameters table).

Families B and C (C carries its cross-linked tyrosine on helix VII) and the nitric oxide
reductases (no cross-linked tyrosine) remain "family not assigned" unless a TIGRFAM model names
them (Sousa et al. 2011). All oxidase genes are kept: family A held every development PCR
product, but other primers or bacteria may amplify families B or C; a family-A-only set would
be about 40 % as large; and the other genes cannot raise B above E, because a eukaryotic COI
amplicon is far less similar to them than to MIDORI2.

### Step 4. Whole genes and fragments

A gene counts as whole if it has a start and a stop codon (`partial=00`), covers at least
70 % of the PF00115 model (all domain hits together, so that a broken alignment is measured as
one), and is 400–1,100 amino acids long. Development whole genes had a median of 555 amino
acids; the range includes the subunit I–III fusions of some haloarchaea, *Thermus* and
*Deinococcus* (517–901 amino acids) and two haloarchaeal genes of 1,007 and 1,065; the upper
limit catches predictions that run on into a neighbouring gene.

A gene that is not whole is kept as a fragment, in a separate file, if it has at least 200
amino acids; in R232 this includes 51 oxidase genes longer than 1,100 amino acids. The primers
amplify only 200–420 bp of a gene of about 1,600 bp, so a fragment covering that region works
as well as a whole gene; in development, 5.9 % of the PCR products lay only in fragments,
mostly from MAGs of uncultured lineages. Fragments are kept apart because they show only part
of the residues of step 3. Shorter pieces carry too little sequence for step 5; the shortest
fragment holding a PCR product had 253 amino acids. MetaFlux uses both files to compute B.

Whole genes and fragments held 1,239 of the 1,259 development PCR products (98.4 %); of the
other 20, five came from contigs of eukaryotic origin, five lay across a gene split in two,
and ten lay in genes broken by a frameshift or an internal stop codon.

### Step 5. Contamination checks

GTDB admits a genome only if CheckM (Parks et al. 2015) and CheckM2 (Chklovski et al. 2023)
estimate completeness ≥ 50 %, contamination < 5 % and a quality score ≥ 50 (Parks et al.
2026); one foreign contig is below what such estimates detect, so every candidate gene is
tested. Contig GC and length are estimated from the contig's genes (same GC verdict as the
contig sequence for all 5,060 development candidates). The table gives the four checks.

| Check | The gene is removed when | Why |
|---|---|---|
| 5.1 reading frame | its length is not a multiple of 3, or it has an internal stop under the table of step 1 | |
| 5.2 GC content | the GC of the gene or of its contig differs from the genome's by more than 0.15 | largest difference of a genuine development gene: 0.132; a thraustochytrid-like COX1 in a planctomycete MAG (gene 0.39, genome 0.66) fails |
| 5.3 match to mitochondrial COI | it is at least 90 % identical, over 80 % of the shorter sequence, to a certainly mitochondrial MIDORI2 record | genuine bacterial genes reached at most 85.4 %; a copepod COI in a bacterial MAG is 99.5 % identical to a mitochondrial record |
| 5.4 closer to eukaryotes | E ≥ 80 % and E > B_other (see below) | a eukaryotic sequence binned into a MAG is closer to MIDORI2 than to any bacterial gene |

- **Certainly mitochondrial** records have no reading frame open under the bacterial code
  (each frame has at least two TAA, TAG or TGA), at least two in-frame TGA, and TGA at least as
  often as TGG, which a table-11 bacterial gene cannot have. A looser definition admitted two
  bacterial records filed as a fish and a planthopper, which then rejected 19 genuine genes.
  2,904,469 of the 3,213,682 records of MIDORI2 GB273 (UNIQ NUC_SP, incomplete names such as
  "sp." included) qualify; chimeric records (step 6) do not count. Coverage is measured on the
  shorter sequence because many records are barcodes of about 650 bp. In R232, check 5.3
  removed no gene.
- **Check 5.4** runs after step 6, on MIDORI2 without the records removed there, so that a
  bacterial record filed under an animal's name cannot make a genuine gene look eukaryotic.
  E counts every record that is not removed, as in MetaFlux, over hits of at least 300
  columns; B_other is the best identity to a partner gene (below). A failing gene is removed if
  it is a fragment or its GC differs from the genome's by more than 0.08; otherwise it is
  marked "review" and also left out, so that no manual decision enters the build. In R232, 100
  genes had an E of at least 80 %: check 5.4 removed 12 and left out 9 for review and 6 under
  the partner rule (`left_out:5.4_partner_rule`).
- **In the development sample** check 5.4 failed three candidates: the thraustochytrid gene of
  check 5.2 (80.1 % to *Aurantiochytrium*), a 55-amino-acid fungal piece in a *Rickettsia* MAG
  (89.3 % to *Phialocephala*), and a 274-amino-acid *Polaromonas* fragment probably removed
  wrongly, through MF319222.1, a record filed as the clam *Corbicula fluminea* that is probably
  bacterial itself (GC 0.64, one TGA).

An animal or fungal COI binned into a bacterial MAG is mostly removed before any check: its
TGA codons are stops under table 11, so at most short reading frames below the 200-amino-acid
floor are predicted (none at a copepod COI, 55 amino acids at a fungal COX1 in development).

**How identity is computed.** Identities in steps 5 and 6 and in MetaFlux are identical
A/C/G/T columns divided by all aligned columns of the VSEARCH alignment, end gaps left out,
internal gaps counted as differences, and columns with N or another ambiguity code left out.
VSEARCH's own identity counts N and compatible ambiguity codes as matches; about 9 % of MIDORI2
records contain ambiguity codes, and in development VSEARCH's value was inflated for 8–15 % of
the bacterial amplicons (one archaeal amplicon appeared 100 % identical to an insect record,
52.8 % over A/C/G/T). Because incompatible codes count as mismatches, the searched copies of
MIDORI2 have every ambiguity code replaced by N (53,877 records carry codes other than N), so a
search cut-off cannot lose a hit; the final test and MetaFlux use the original sequences.

Three observations are recorded as flags and do not remove a gene:

| Flag | Meaning | In R232 |
|---|---|---|
| cross-phylum cluster | clusters at 95 % with genes of another phylum, as one contaminant binned into MAGs of several lineages would; B is not affected, but the reported bacterial lineage can be wrong | 522 kept genes |
| mitochondrion-like contig | contig shorter than 20 kb that also carries cytochrome b and COX2 or COX3; bacterial bc1–oxidase operons trigger it too (37 development genes, all bacterial) | 2,247 genes |
| duplicate copy | two whole genes of the same genome at least 99 % identical, as in a duplicated contig | 373 kept genes |

#### Searches against MIDORI2

VSEARCH tries long records first, so a search that keeps only the first few hits of each
gene misses close short records. The build therefore searches in both directions, and again
wherever a search stopped at its limit:

| Search | What it covers |
|---|---|
| reverse search | every MIDORI2 record against the genes that passed check 5.2 (16 accepted, 128 rejected candidates) |
| deeper searches | records with 0–15 hits again with no accept limit (1,024 rejected); records that filled 16 slots again with no limit |
| records rich in N | the 1,762 records with more than 20 % non-A/C/G/T, searched on their own; 2 had a hit of at least 77 % |
| genes near a decision | genes within 3 points of a decision, again with no limit against MIDORI2 and the other genes; 186 in the first round, no outcome changed, no new gene in the second |
| merging | each gene–record pair decided on its best alignment; 88,644 of the 142,580 pairs were found only by the new searches |

On 400 GB273 records with 1–15 hits in the first search, the searches still missed 12 of
1,056 pairs at ≥ 78 %, 2 of 445 at ≥ 80 %, and none of 425 at ≥ 80 % over at least 300
columns; no gene's best identity changed. Every search and its settings are in the README.

#### Length floor of the 80 % level

A MIDORI2 hit counts at the 80 % level only with at least 300 counted columns (A/C/G/T columns
and internal gaps): for E, for the certainly mitochondrial identity used as gene evidence, and
for records examined below 90 % in step 6, but not for check 5.3 or records examined at 90 %.
In R232, 740 genes had a best hit of at least 80 % that was shorter.

| Outcome without the floor | Outcome with the floor | Genes |
|---|---|---|
| left out for review (5.4) | pass | 128 |
| removed (5.4) | pass | 24 |
| left out by the partner rule (5.4) | pass | 4 |

Each row counts genes whose outcome the floor changes, with the same B_other and the same
cleaned MIDORI2.

#### Partner rule

B_other counts only a gene of another genome of the same GTDB family (the name exactly,
suffixes included) that passes check 5.4 itself, so that genes cannot vouch for each other
(three genes close to the mitochondrial COX1 of the green alga *Ostreococcus tauri* would
otherwise pass as each other's partners). Because a failing gene can take away
another gene's partner, step 6 and check 5.4 are repeated until nothing changes (2 passes in
R232). A gene that fails only because of this rule is left out as `left_out:5.4_partner_rule`:
in R232 those three (`GB_GCA_002726005.1|PBRS01000001.1_23`,
`GB_GCA_963999645.1|CAXEUY010000032.1_4`, `GB_GCA_002700775.1|PAEL01000127.1_2`) and `GB_GCA_043746065.1|JAYTRS010000187.1_6`,
`GB_GCA_043776615.1|JAYVRC010000321.1_1` and `GB_GCA_913031655.1|CAJWFY010000134.1_4`.

### Step 6. Cleaning MIDORI2

The genes that pass checks 5.1–5.3 are compared with every record of MIDORI2 GB273 (UNIQ
NUC_SP, which contains every record of the NUC variant). Genuine COI and its bacterial
homologues are never very close (in development, no eukaryotic sequence was more than 77.6 %
identical to a bacterial gene), so a record at least 90 % identical to a gene is either a
bacterial sequence under a eukaryote's name, or the gene is a eukaryotic sequence binned into a
bacterial genome. Records 80–90 % identical are examined too, under the conditions below. Two
independent lines of evidence say which side is misfiled:

| Evidence | Bacterial | Mitochondrial | None |
|---|---|---|---|
| the gene's relatives | closest partner gene (step 5) ≥ 80 % and at least as close as the closest certainly mitochondrial record | a certainly mitochondrial record ≥ 80 % (over ≥ 300 columns) and closer than any partner gene | neither |
| the record's tryptophan codons | likelihood ratio ≥ 20 for a bacterial gene | likelihood ratio ≥ 20 for genuine COI of its lineage | fewer than 5 tryptophan codons, no lineage share, or a lineage whose mitochondria read TGA as a stop |

Each row gives when a line of evidence says bacterial, mitochondrial or nothing. In a
bacterial gene every tryptophan is TGG, while mitochondria that read TGA as tryptophan use it
at very different rates (94 % of tryptophan codons in insect records, 21 % in mushrooms, 15 % in
the liver-fluke family Fasciolidae), so the codons are compared with the TGA share of the
record's own genus or family (parameters table); 20 is the level Kass and Raftery (1995) call
strong. Taxa without a mitochondrial code take that of the nearest ancestor (5,565 records;
5,562 then read TGA as tryptophan). TGA is a stop in the mitochondria of plants, brown algae,
oomycetes, diatoms, chrysophytes, and green algae such as *Micromonas*.

**Decision.** A record is removed when the gene's evidence is bacterial, from a gene that passes
check 5.4, and neither the gene nor the codons say mitochondrial; otherwise it is listed, and
check 5.4 judges the gene with the record still in MIDORI2. Codons never remove a record on
their own, because a missing TGA is no evidence in lineages that rarely use TGA; they can only
prevent a removal. Two lines are needed because a single bacterial record left in MIDORI2 can
make check 5.4 reject many genes of its bacterium: MF319222.1 (one TGA among six tryptophan
codons, 97.9 % identical to a *Paenacidovorax* gene) made 22 genes of Burkholderiales fail.

| Records | Matched gene | Evidence | Decision |
|---|---|---|---|
| six of *Rhizoctonia*, *Ganoderma* and *Fomes* (mushrooms), 91–98 % | MAG of a soil archaeon (Nitrososphaeraceae); in fact the first exon of a fungal COX1 read into the open reading frame of its group I intron (two LAGLIDADG motifs, Pfam PF00961); the closest record, KC352446.1, is that intron frame alone, filed as COI | codons: none (lineages use TGA for 1–5 % of tryptophan codons); gene: no relative in another genome | listed; check 5.4 removes the gene |
| carrion beetle *Calosilpha brunnicollis*, one TGA among seven tryptophan codons, 96.8 % | *Aquirickettsiella*, a bacterial pathogen of arthropods | gene: same-family relative, bacterial; codons: the genus uses TGA for almost every tryptophan codon, they agree | removed |
| three *Dinobryon pediforme* (chrysophyte), 99.3–99.6 % | MAG of the Paracaedibacterales, alphaproteobacteria that include endosymbionts of protists | codons: none (TGA is a stop); gene: same-family relative, bacterial | removed |
| *Fasciola hepatica* (liver fluke), 96.0 % | *Neorickettsia findlayensis*, an endosymbiont of trematodes, with relatives in two other *Neorickettsia* genomes (89–90 %) | codons: none (rare TGA); gene: bacterial | removed |

Each row is an example from the R232 build: the records, the gene they match, the two lines of
evidence and the decision. The mushroom row shows a gene removed instead of the records.

In R232, 446 sequences (498 records) were examined: 141 at ≥ 90 % identity to a gene that
passed check 5.3 and 305 at 80–90 %. 187 sequences (236 records) were removed and 261
sequences (264 records) listed.

| Decision | Reason | Sequences |
|---|---|---|
| removed | gene looks bacterial | 154 |
| removed | gene looks bacterial; codons agree | 26 |
| removed | no MIDORI2 relative (isolation line) | 4 |
| removed | chimera (window scan) | 3 |
| listed | gene looks mitochondrial | 239 |
| listed | no evidence either way | 15 |
| listed | codons look mitochondrial | 7 |

Counts are per sequence and include the chimeric records found by the window scan (below). Of
the removed sequences, 123 were examined at ≥ 90 % and 61 at 80–90 %; the other 3 are the
chimeras. Of the 956 genes behind these decisions, 932 look bacterial, 6 mitochondrial and 18
give no evidence. MIDORI2 itself already excludes 2,215 GenBank records found mislabelled at
the order, class or phylum level (Leray et al. 2019, 2022); bacterial sequences submitted under
the host's name are the separate problem this step addresses.

#### Records 80–90 % identical, and the lineage guard

A bacterial record 80–90 % identical to the gene of its own bacterium can make that gene fail
check 5.4. Records 80–90 % identical, over at least 300 columns, to a gene that passed check
5.3 are therefore examined when the record's lineage reads TGA as tryptophan. Genuine COI can
reach 80 % identity to a bacterial gene in some lineages, and where TGA is rare the codons
cannot prevent a wrong removal; the **lineage guard** therefore skips records whose lineage
uses TGA for less than 0.05 of its tryptophan codons (the family's share, otherwise the
genus's; a record whose share is unknown is examined). In R232, 421 sequences were candidates
at 80–90 %; the guard left out 116 and 305 were examined.

| Family of the record | Sequences left out by the guard |
|---|---|
| Pleurotaceae (oyster mushrooms) | 42 |
| Boletaceae (boletes) | 34 |
| Peronosporaceae (oomycetes) | 24 |
| Bacillariaceae (diatoms) | 10 |
| Naviculaceae (diatoms) | 3 |
| Ceratobasidiaceae (fungi) | 2 |
| Sarcophytidae (soft corals) | 1 |

All the families left out use TGA rarely or not at all. On the release, records kept in
MIDORI2 and genes that pass with bacterial evidence form 695 pairs at ≥ 80 % over at least 300
columns, between 27 records and 73 genes: 692 pairs with raphidophyte records
(Raphidophyceae) and 2 with a green-algal record (Picocystophyceae), both of mitochondrial
code 1 and so not examined, and 1 with a sponge record (Demospongiae), examined and listed.

#### Isolation line

A bacterial record without TGA, 80–90 % identical to a gene with no partner gene, gives no
evidence on either side; left in MIDORI2 it would make its gene fail check 5.4, and MetaFlux
would report the bacterium under the record's eukaryotic name. A genuine COI record usually has
relatives in MIDORI2. An examined record is therefore removed when no other MIDORI2 record
reaches 75 % identity to it, confirmed by a search with no limits, and it is at least 80 %
identical to a released gene, unless its codons or one of its genes look mitochondrial (full
rule in the parameters table).

In R232, 12 sequences had no MIDORI2 relative of at least 75 %. Four were removed only through
this line: AF120671.1 (filed as the shipworm *Bankia carinata*), HQ198594.1 (the amphipod
*Gammarus komareki*), MF743154.1 (a water flea, Daphniidae) and MK308055.1 (the comb jelly
*Beroe ovata*). The genes they matched are among those the regression checks expect to be
kept (step 7).

#### Chimeric records

Some animal COI records in MIDORI2 carry a stretch of the COI homologue of an endosymbiont such
as *Wolbachia* or *Tisiphia*, while the rest reads TGA as tryptophan; such a record can pass as
certainly mitochondrial.

- **Chimera screen.** Every gene–record alignment is scanned in windows of 300 columns; a
  certainly mitochondrial record with one window at ≥ 90 % and another below 80 % is chimeric,
  as is one ≥ 90 % identical to a record removed in step 6. Chimeric records are not used as
  certainly mitochondrial evidence.
- **Chimera rule.** A record with such an alignment is removed when its best window is ≥ 95 %
  identical to a released gene that is not on a mitochondrion-like contig, the whole record is
  below 90 % identical to that gene, and the window carries no in-frame TGA while the rest of
  the record does (full rule in the parameters table). Left in MIDORI2, such a record hides the endosymbiont:
  its ASV matches the stretch as closely as the bacterial gene and stays among the animal
  ASVs. Other chimeric certainly mitochondrial records stay in MIDORI2 and are listed.

In R232 the window scan flagged EF518436.1 (filed as the ant *Pheidole truncula*, *Wolbachia*
stretch) and KM679400.1 (*Mycetophila lunata*, *Tisiphia* stretch). The chimera rule removed
both and GU446981.1 (*Clastoptera proteus*, a stretch 96.0 % identical to a *Rickettsia*
gene). The scan sees only alignments above the search cut-off (75 % for the forward search
against certainly mitochondrial records, 77 % for the others), so a chimera with a smaller
bacterial stretch is caught only through a removed record ≥ 90 % identical to it; on GB273 the
two chimeras were reported just above the cut-off (EF518436.1 by 0.0–0.4 points, KM679400.1 by
1.3–2.5 points).

#### Removal by sequence

Step 6 decides once per sequence, identified by the SHA-1 of the sequence in upper case, so
identical records filed under different names are removed or listed together. Of the 2,997,664
distinct sequences of GB273, 100,255 are carried by two or more records (316,273 records); 52
rows of the list are identical copies, listed with the decision of their sequence. MetaFlux
removes by sequence in the same way, so the MIDORI2 it prepares equals the cleaned MIDORI2 of
the build: 3,213,446 records.

#### Built-in stops

The build stops instead of releasing contradictory decisions, for example a gene whose E is
below a listed record paired with it, or a removed record without a released gene at its
examination level; the five stops are in the parameters table. No stop was reached.

#### Sensitivity of step 6 and check 5.4

| Setting | Value | Moved to | Genes | Records |
|---|---|---|---|---|
| 5.3 mitochondrial match (%) | 90 | 88 | 3 | 136 |
| | | 92 | 0 (lower bound) | 0 |
| 5.4 closer to eukaryotes (%) | 80 | 78 | ≥ 18 | 0 |
| | | 82 | 9 | 0 |
| 5.4 GC limit | 0.08 | 0.06 | 1 | 0 |
| | | 0.10 | 1 | 0 |
| records examined from (%, any lineage) | 90 | 88 | 0 | ≥ 6 |
| | | 92 | 3 | 2 |
| records examined from, below 90 (%) | 80 | 78 | 0 | ≥ 1,284 |
| | | 82 | 4 | 121 |
| records below 90 % examined | yes | no | 34 | 315 |
| lineage guard (TGA share) | 0.05 | 0 | 0 | 122 |
| | | 0.10 | 4 | 7 |
| length floor of the 80 % level (columns) | 300 | 200 | 4 | 25 |
| | | 400 | 0 | 173 |
| gene-evidence identity (%) | 80 | 78 | 0 | 0 |
| | | 82 | 4 | 5 |
| codon error rate | 0.01 | 0.003 or 0.03 | 0 | 0 |
| likelihood-ratio threshold | 20 | 10 or 50 | 0 | 0 |
| tryptophan codons per lineage | 200 | 100 | 0 | 1 |
| | | 500 | 0 | 41 |
| tryptophan codons per record | 5 | 3 | 0 | 1 |
| | | 10 | 0 | 0 |
| partner rank | family | order | 0 | 0 |
| | | any genome | 6 | 0 |
| isolation line | on | off | 4 | 6 |
| isolation limit (%) | 75 | 73 | 1 | 1 |
| | | 77 | 0 | 0 |

Each row moves one setting and runs the rules of step 6 and check 5.4 again on the hits of the
build; "Genes" counts genes whose outcome changes, "Records" MIDORI2 records whose decision
changes. A count marked ≥ or "lower bound" may be higher, because the searches do not fully
cover the moved value. No codon setting changes a gene's outcome; the largest changes come from
the examination below 90 % (34 genes and 315 records when switched off, at least 1,284 records
at 78 %).

### Step 7. Final test

The final test asks how the set handles bacteria and eukaryotes that the references do not
contain, on genomes and MIDORI2 records not used to set any rule. It searches at the depth
MetaFlux uses (256 accepted and 2,048 rejected
candidates, `--id 0.5`, whole genes and fragments searched separately), comparing B and E
unrounded. The rules were fixed before the test, except the chimera rule of step 6, after
which step 6 and the final test were run again.

#### Bacterial amplicons

- **Test genomes.** 3,476 R232 species representatives not in the development sample, chosen
  with the development design and a new seed: all 176 such Rickettsiales (4 of them
  *Wolbachia*), 3,000 other bacteria drawn per phylum, and 300 archaea.
- **In-silico PCR.** As in development, but with up to three mismatches per primer, on purpose
  permissive to get more and more distant amplicons. It runs on GTDB's predicted genes; in
  development this missed only products in genes broken by a frameshift or internal stop (10
  of 1,259, 0.8 %). Products with an intact reading frame under the bacterial code are the
  queries: 936 from 738 genomes (803 BF3/BR2, 102 Leray-XT, 31 fwhF2/fwhR2n), 181 of them with
  at most one mismatch per primer. 928 of the 936 come from subtype A1.

| Genomes | BF3/BR2 | Leray-XT | fwhF2/fwhR2n |
|---|---|---|---|
| *Wolbachia* (4) | 4 (4) | 4 (4) | 3 (0) |
| other Rickettsiales (172) | 72 (45) | 73 (31) | 25 (0) |
| other bacteria (3,000) | 610 (86) | 17 (3) | 3 (0) |
| archaea (300) | 7 (1) | 0 | 0 |

Each cell is the number of genomes giving a product with up to three mismatches per primer; in
brackets, with at most one. With at most one mismatch, BF3/BR2 amplifies 86 of the 3,000 other
bacteria (2.9 %), close to the 3.2 % of the development table.

Every test genome is in the set, so its amplicon would find itself at 100 %. B is therefore
computed without the genes of the amplicon's own species and, in a second run, its genus, and a
record removed in step 6 whose closest gene belongs to that species (or genus) is put back for
E. Measured is the share of amplicons with B > E, which MetaFlux would recognise as bacterial.

#### Eukaryotic sequences

The test draws MIDORI2 GB273 records not used during development, up to 250 taxa per group
(arthropods, chordates, other animals, other eukaryotes), at four levels: 1,000 queries at
species, 1,000 at genus, 858 at family and 478 at order level. The 3,336 queries come from
3,245 distinct records, because 87 records serve at more than one level. Each query is cut to
the inserts of the three primer pairs (10,008 inserts). For each
insert, every record of its own species (genus, family or order) is left out of MIDORI2, and B
is computed against the whole set. Measured is the share with B > E: genuine eukaryotes that
MetaFlux would wrongly move to the bacterial table.

#### Known contaminants

The three contaminants found during development must still be rejected; unlike the other
tests, this one uses development data.

#### Release criteria and results

The criteria, fixed before the test, are a recognised share, with the species left out, of at
least 90 % for BF3/BR2, 85 % for Leray-XT and 80 % for fwhF2/fwhR2n, and at most 0.5 % of
eukaryotic sequences moved at each of the four levels. All criteria were met.

| Bacterial amplicons recognised (B > E) | species left out | genus left out | development (species left out) |
|---|---|---|---|
| all (936) | 98.4 % | 95.7 % | 94 % |
| BF3/BR2 (803) | 99.9 % | 98.8 % | 96 % |
| Leray-XT (102) | 87.3 % | 77.5 % | 91 % |
| fwhF2/fwhR2n (31) | 96.8 % | 77.4 % | 85 % |
| *Wolbachia* (11) | 72.7 % | 0 % | 97 % |
| other Rickettsiales (178) | 93.3 % | 85.4 % | |
| at most 1 mismatch per primer (181) | 94.5 % | 87.8 % | |

Each value is the share of bacterial amplicons with B > E, with the amplicon's species or genus
left out of the set; the last column gives the share measured during development.

- **Eukaryotic sequences moved to the bacterial side:** 0.03 % with the species left out (1
  insert of 3,000), 0 % with the genus, 0.08 % with the family (2 of 2,574) and 0 % with the
  order left out (1,434 inserts). All three are far from both
  references (E 58.5–71.2 %, B 59.9–73.2 %).
- **Known contaminants:** all three stay out of the set.
- **The *Wolbachia* figures rest on 11 amplicons from four genomes.** The three missed come
  from *Wolbachia* sp964277405, whose own *coxA* is in MIDORI2 under the pseudoscorpion
  *Neobisium carcinoides*; step 6 removes that record, but the test puts it back, and with the
  species left out the closest bacterial gene is less close than the record. With the genus
  left out, every *Wolbachia* record of MIDORI2 comes back and no *Wolbachia* amplicon wins, as
  in development (0.7 %): a *Wolbachia* lineage missing from the set is recognised only if
  MIDORI2 holds no record of it under an animal's name.

#### Regression checks

| Check | Result |
|---|---|
| eleven gene–record pairs the searches must find, among them KF697361.1 with `GB_GCA_004211595.1|SEDO01000025.1_7` (92.2 % over 578 columns) and pairs only the deeper searches find | all present, within 0.5 points of the expected identity |
| the 299 chrysophyte and synurophyte records the build does not remove, each against the cleaned MIDORI2 without itself and against the release | none with B > E; 1 with B and E below 70 % (no assignment possible) |
| eukaryotic test records of the final test | none removed, by record ID or by sequence |
| 37 records where genuine COI is close to bacterial genes (27 of them *Chattonella*, a raphidophyte) | all stay in MIDORI2; the 72 genes at 80.2–82.2 % to them (69 Pelagibacteraceae) all pass |
| 16 genes named in advance as genuine, among them those of the isolation line | all kept |

Each row is a check run with the final test, and its result; all passed. Any rule changed after
this test makes a new version.

### Step 8. Release files

| File | Content |
|---|---|
| `bacterial_cox1_gtdb-r232.v2.fasta.gz` | whole genes that pass all checks; DNA; SINTAX headers |
| `bacterial_cox1_gtdb-r232.v2_fragments.fasta.gz` | fragments of at least 200 amino acids that pass all checks |
| `bacterial_cox1_gtdb-r232.v2_metadata.tsv.gz` | every candidate gene: genome, lineage, coordinates, scores, residues, family, type, the outcome of every check, E, B_other and the searches that found them, flags and verdict |
| `bacterial_cox1_gtdb-r232.v2_midori2_gb273_records.tsv` | the MIDORI2 GB273 records examined in step 6 and the chimeric records, one row per record: the matching gene, the decision (removed or listed) and its reason, the SHA-1 of the sequence (`sequence_sha1`) and the SHA-256 of the MIDORI2 file (`midori2_file_sha256`) |
| `README.md` | inputs and their checksums, tool and model versions, every VSEARCH search with its settings, the rules of step 6, counts per step, built-in stops, sensitivity of the counts to the thresholds, results of the final test and the regression checks |
| `scripts.zip`, `environment.yaml`, `environment_explicit.txt`, `hmm_models.md5` | build scripts as they were run, with unit tests; the conda environment, also with its exact packages; MD5 sums of the nine HMM models; to rebuild the release from the same inputs |
| `SHA256SUMS` | checksums of all files |

Each FASTA header carries the gene ID, the GTDB accession, the family and type label and the
GTDB taxonomy, for example:

```text
>RS_GCF_000008385.1|NC_006833.1_429;family=A1;type=aa3_caa3;tax=d:Bacteria,p:Pseudomonadota,c:Alphaproteobacteria,o:Rickettsiales,f:Anaplasmataceae,g:Wolbachia,s:Wolbachia_sp000008385
```

The release is named after the GTDB release and a build number (`gtdb-r232.v2`). A later GTDB
release, or any change of a rule, produces a new version under the same Zenodo concept DOI
(10.5281/zenodo.23820759). MetaFlux pins the version DOI and the SHA-256 of the files it
downloads. The `.gz` files are written without a time stamp, so a rebuild with the same content
gives the same SHA-256. The steps up to check 5.2 reuse the outputs of the earlier build of
the set, made with the same scripts and inputs. If you use
the set, cite the release by its version DOI, together with GTDB (Parks et al. 2026) and
MIDORI2 (Leray et al. 2022).

### Parameters

The table lists the inputs and every threshold of the set, step by step.

| Step | Parameter | Value |
|---|---|---|
| 0 | GTDB gene file | `gtdb_proteins_nt_reps_r232.tar.gz`, 187 GB (174 GiB), MD5 `bd0117f737be1fea6819365f2375c048` |
| 0 | MIDORI2 | GB273, `MIDORI2_UNIQ_NUC_SP_GB273_CO1_SINTAX.fasta.gz`, version DOI 10.5281/zenodo.22685701 |
| 0 | HMMs | Pfam 38.2: PF00115.26, PF00116.27, PF00510.24, PF00033.25; NCBI hmm_PGAP 20.0: TIGR02891.1, TIGR02843.1, TIGR02882.1, TIGR00780.1, NF050116.1 |
| 1 | genetic code | table 11; table 4 if > 1 % of the genes with `partial=00` contain an in-frame TGA before the stop codon |
| 2 | HMM search | `hmmsearch --cut_ga` (PF00115 gathering threshold 34.8 bits) |
| 3 | oxidase | PF00115 ≥ 100 bits and at most 1 of the 6 ligand histidines replaced |
| 3 | family A | Tyr244 with Glu242 (A1) or with Tyr238 and Ser239 (A2) |
| 3 | type label | best model at its trusted cutoff (`--cut_tc`) |
| 4 | whole gene | `partial=00`, model coverage ≥ 0.70 (all domain hits together), 400–1,100 amino acids |
| 4 | fragment | not whole, ≥ 200 amino acids |
| 5.1 | reading frame | length a multiple of 3, no internal stop codon |
| 5.2 | GC content | gene and contig within 0.15 of the genome |
| 5.3 | mitochondrial match | < 90 % identity over ≥ 80 % of the shorter sequence to certainly mitochondrial MIDORI2 records |
| 5.3 | chimeric records | not counted as certainly mitochondrial: a window of 300 columns ≥ 90 % and another < 80 % to a gene, or ≥ 90 % to a removed record |
| 5.4 | closer to eukaryotes | run after step 6, on the cleaned MIDORI2 (every record not removed counts); E ≥ 80 % over ≥ 300 columns and above B_other: removed if fragment or GC difference > 0.08, otherwise "review" (left out) |
| 5.4 | partner rule | B_other from genes of another genome of the same GTDB family that pass 5.4 themselves; step 6 and 5.4 repeated until nothing changes; a gene that fails only through this rule: `left_out:5.4_partner_rule` |
| 5 | flags | cross-phylum cluster at 95 % (a gene joins the closest accepted cluster centre); mitochondrion-like contig < 20 kb; duplicate copy ≥ 99 % |
| 5–6 | identity | identical A/C/G/T columns / (A/C/G/T columns + internal gap columns); end gaps and ambiguity codes left out; VSEARCH 2.32.0 `--usearch_global`, coding strand only, settings in the README |
| 5–6 | searches | gene → MIDORI2 (20 accepted, 200 rejected); MIDORI2 → genes (16, 128), again with no accept limit (1,024 rejected) for records with 0–15 hits and with no limit for records with 16; records with > 20 % non-A/C/G/T searched on their own; genes within 3 points of a decision searched with no limit (at most 3 rounds); best alignment per gene–record pair |
| 5–6 | search copies of MIDORI2 | every ambiguity code replaced by N |
| 6 | MIDORI2 record examined | ≥ 90 % to a gene that passes 5.1–5.3; or 80–90 % over ≥ 300 columns when the record's lineage reads TGA as tryptophan and does not use TGA for < 0.05 of its tryptophan codons (family's share, otherwise genus's; examined when unknown) |
| 6 | unit of decision | the sequence (SHA-1 of the sequence in upper case); every record of a removed sequence is removed |
| 6 | gene's evidence | bacterial: closest partner gene (5.4) ≥ 80 % and at least as close as the closest certainly mitochondrial record; mitochondrial: certainly mitochondrial record ≥ 80 % over ≥ 300 columns and closer; otherwise none |
| 6 | record's evidence | likelihood ratio of its TGA/TGG counts, bacterial gene (TGA error rate 0.01 per tryptophan codon) vs genuine COI (TGA share of its genus, otherwise its family, each with ≥ 200 tryptophan codons); ≥ 20 bacterial, ≤ 1/20 mitochondrial; none with fewer than 5 tryptophan codons, or where TGA is a mitochondrial stop (code inherited from the nearest ancestor when the taxid has none) |
| 6 | MIDORI2 record removed | gene's evidence bacterial, from a gene that passes 5.4, and neither the gene nor the codons mitochondrial; never on codons alone; otherwise listed and kept |
| 6 | isolation line | removed: closest other MIDORI2 record < 75 % (identical copies aside; relatives over ≥ 300 columns or 80 % of a shorter record; confirmed by a search with no limits), ≥ 80 % over ≥ 300 columns to a gene that passes 5.4 and is not on a mitochondrion-like contig, codons and genes not mitochondrial |
| 6 | chimera rule | removed: chimeric alignment, best window of 300 columns ≥ 95 % to a released gene not on a mitochondrion-like contig, whole alignment < 90 %, in-frame TGA only outside the window (at least 1) |
| 6 | built-in stops | E below a listed paired record; listed record ≥ 90 % to a released gene with B_other < E; removed record without a released gene at its examination level; genes near a decision left unsearched; gene passing 5.3 only through a late chimera |
| 7 | release criteria | species left out: ≥ 90 % (BF3/BR2), ≥ 85 % (Leray-XT), ≥ 80 % (fwhF2/fwhR2n) recognised; eukaryotes: ≤ 0.5 % moved at each level; searches at MetaFlux's depth (`--id 0.5`, 256 accepted, 2,048 rejected) |
| README | sensitivity | verdict changes when each value moves: GC limits ± 0.02, fragment floor ± 50 amino acids, identity limits ± 2 points, model coverage ± 0.10, oxidase score ± 20 bits; step 6: error rate 0.003 and 0.03, likelihood-ratio threshold 10 and 50, 100 and 500 tryptophan codons per lineage, 3 and 10 per record, length floor 200 and 400 columns, lineage guard 0 and 0.10, isolation limit 73 and 77 %, partner rank order and any genome; records below 90 % and the isolation line switched off |

## Copy number

| Near-complete genomes (CheckM2 completeness ≥ 90 %, contamination ≤ 5 %) | Development sample (weighted to GTDB) | R232 (120,230 bacterial genomes) |
|---|---|---|
| family A: carriers with one copy | 65 % (35 % with two to five) | 64 % of 71,675 carriers, up to nine |
| family A in *Wolbachia* | all 51 carriers have one copy (99 % of other Rickettsiales) | 55 of 56 genomes have one copy, one has none |
| all subunit-I genes: carriers with two or more | 69 %, up to nine | 67 % of 79,032 carriers, up to twelve |

Family A is mostly single-copy; all subunit-I genes together are mostly multi-copy, because
many bacteria carry several oxidase types, and in 90 % of the multi-copy development genomes
all copies are less than 80 % identical (distinct genes, not assembly artefacts). All copies
that pass are kept: MetaFlux uses only the best match of each ASV, so extra copies give a
bacterial ASV more chances to find a close relative without biasing the result.

## Limitations, open issues and possible improvements

| Limitation | What it means |
|---|---|
| coverage of bacterial diversity | the set contains what GTDB contains: with the genus left out, 95.7 % of the final-test amplicons were recognised (98.4 % with the species left out; BF3/BR2 98.8 %, Leray-XT 77.5 %, fwhF2/fwhR2n 77.4 %), and no *Wolbachia* amplicon, because MIDORI2 holds *Wolbachia* *coxA* under animal names; lineages without cultured or assembled members are not represented |
| one genome per species | an amplicon from a divergent strain should be recognised at least as often as one whose whole species is missing (94 % in development); not tested directly |
| rules set on the data they apply to | the residue rule of step 3, the 200-amino-acid floor, check 5.4 with its GC limit of 0.08 and the definition of certainly mitochondrial records were set on the development sample; the evidence rule of step 6 and the rules set later were fixed after looking at the R232 genes and GB273 records, so no independent data have tested them; the README reports the changes when each of twenty settings moves |
| contaminants that can pass | check 5.3 relies on TGA read as tryptophan, so contamination from plants, brown algae, oomycetes, diatoms or thraustochytrids meets only checks 5.2 and 5.4: a contaminant with a GC close to the genome's and no close MIDORI2 relative passes, and so does one in several MAGs of the same GTDB family, each copy being the other's partner |
| bacterial records left in MIDORI2 | records of lineages GTDB lacks, below the examination levels of step 6, left out by the lineage guard, or with unclear evidence stay in MIDORI2 and in E (42 records whose codons look bacterial remain); such a record can make check 5.4 reject a genuine gene, as with the *Polaromonas* fragment in development |
| genes left out | the 9 genes marked "review" and the 6 left out by the partner rule may include genuine genes, leaving their lineages fewer genes on the bacterial side |
| eukaryotes moved to the bacterial table | in development, 0.04 % of animal sequences whose family, and 0.14 % of those whose order, was missing from MIDORI2; in groups MIDORI2 covers poorly, check that table before concluding that a group is absent |
| family labels | one alignment to the Pfam model; divergent sequences can be misaligned, and nitric oxide reductases are not separated from families B and C; the classifier of Sousa et al. (2011) is no longer available at the address given in the paper |
| genetic codes | the table GTDB used is inferred; lineages with code 25 (Gracilibacteria, Absconditabacteria; TGA read as glycine) are predicted under table 4 or 11: under 11 some genes are truncated or missing, under 4 TGA is translated as tryptophan |
| the MIDORI2 list is a lower bound | short records with few tryptophan codons, records whose gene has no partner and whose codons say nothing (unless the isolation line applies), and chimeras with a small bacterial stretch can be missed; mitochondrial sequences that are not COI, such as the intron reading frame of step 6, stay |
| three primer pairs, shares of taxa | recognition rates were measured for BF3/BR2, Leray-XT and fwhF2/fwhR2n only, as shares of genomes or sequences, not reads; a few abundant ASVs, often *Wolbachia*, can hold most bacterial reads |
| what the set does not predict | whether an amplicon is bacterial, not how much bacterial DNA a sample yields; that depends on primer mismatches, annealing temperature and target DNA, and off-target amplification is usually much larger in eDNA than in bulk animal tissue (Elbrecht et al. 2019; Collins et al. 2019) |

Each row names a limitation of the set or of its tests and what follows from it for a reader
of MetaFlux results. Possible improvements: placement on a reference tree for ASVs far from
both references, as DARN does; protein-level comparison for very divergent ASVs; a table of the
bacterial lineages each common primer pair can amplify; more genomes per species; genomes
outside GTDB, with genes predicted by Pyrodigal (Larralde 2022), tested to give the same
predictions as the current Prodigal code (v2.6.3+31b300a) and faster; a curated classifier for
the oxidase families; sharing the MIDORI2 list with the MIDORI2 maintainers.

## Related pages

- [COI](coi.md)
- [Keeping and discarding taxa](../taxon-filter.md)
- [Output files](../../reference/output.md)
- [Citation and references](../../about/citation.md)

## References

References cited on this page. The tools MetaFlux runs, and the references of the COI
marker itself, are on [Citation and references](../../about/citation.md).

- Baldo, L., et al. (2006). Multilocus sequence typing system for the endosymbiont *Wolbachia
  pipientis*. *Applied and Environmental Microbiology* 72:7098–7110.
  doi:10.1128/AEM.00731-06
- Buchner, D., Leese, F. (2020). BOLDigger – a Python package to identify and organise
  sequences with the Barcode of Life Data systems. *Metabarcoding and Metagenomics* 4:e53535.
  doi:10.3897/mbmg.4.53535
- Chklovski, A., et al. (2023). CheckM2: a rapid, scalable and accurate tool for assessing
  microbial genome quality using machine learning. *Nature Methods* 20:1203–1212.
  doi:10.1038/s41592-023-01940-w
- Collins, R. A., et al. (2019). Non-specific amplification compromises environmental DNA
  metabarcoding with COI. *Methods in Ecology and Evolution* 10:1985–2001.
  doi:10.1111/2041-210X.13276
- Curd, E. E., et al. (2019). Anacapa Toolkit: an environmental DNA toolkit for processing
  multilocus metabarcode datasets. *Methods in Ecology and Evolution* 10:1469–1475.
  doi:10.1111/2041-210X.13214
- Eddy, S. R. (2011). Accelerated profile HMM searches. *PLoS Computational Biology*
  7:e1002195. doi:10.1371/journal.pcbi.1002195
- Edgar, R. C. (2016). SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS
  sequences. *bioRxiv*. doi:10.1101/074161
- Elbrecht, V., et al. (2019). Validation of COI metabarcoding primers for terrestrial
  arthropods. *PeerJ* 7:e7745. doi:10.7717/peerj.7745
- González, A., et al. (2023). VTAM: a robust pipeline for validating metabarcoding data using
  controls. *Computational and Structural Biotechnology Journal* 21:1151–1156.
  doi:10.1016/j.csbj.2023.01.034
- González-Miguéns, R., et al. (2025). A novel taxonomic database for eukaryotic mitochondrial
  cytochrome oxidase subunit I gene (eKOI), with a focus on protists diversity. *Database*
  2025:baaf057. doi:10.1093/database/baaf057
- Hyatt, D., et al. (2010). Prodigal: prokaryotic gene recognition and translation initiation
  site identification. *BMC Bioinformatics* 11:119. doi:10.1186/1471-2105-11-119
- Kass, R. E., Raftery, A. E. (1995). Bayes factors. *Journal of the American Statistical
  Association* 90:773–795. doi:10.1080/01621459.1995.10476572
- Larralde, M. (2022). Pyrodigal: Python bindings and interface to Prodigal, an efficient
  method for gene prediction in prokaryotes. *Journal of Open Source Software* 7(72):4296.
  doi:10.21105/joss.04296
- Leray, M., et al. (2019). GenBank is a reliable resource for 21st century biodiversity
  research. *PNAS* 116:22651–22656. doi:10.1073/pnas.1911714116
- Leray, M., Knowlton, N., Machida, R. J. (2022). MIDORI2: a collection of quality controlled,
  preformatted, and regularly updated reference databases for taxonomic assignment of
  eukaryotic mitochondrial sequences. *Environmental DNA* 4:894–907. doi:10.1002/edn3.303
- Li, W., et al. (2021). RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach
  with protein family model curation. *Nucleic Acids Research* 49(D1):D1020–D1028.
  doi:10.1093/nar/gkaa1105
- Meglécz, E. (2023). COInr and mkCOInr: building and customizing a nonredundant barcoding
  reference database from BOLD and NCBI using a semi-automated pipeline. *Molecular Ecology
  Resources* 23:933–945. doi:10.1111/1755-0998.13756
- Nowak, K. H., Buczek, M., Marszałek, M., Prus-Frankowska, M., Valdivia, C., Deng, J.,
  Shropshire, J. D., Łukasik, P. (2026). Insect COI barcoding data as an untapped resource for
  surveying *Wolbachia* symbioses. *bioRxiv* (preprint, not peer reviewed).
  doi:10.64898/2026.06.24.734267
- Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., Tyson, G. W. (2015). CheckM:
  assessing the quality of microbial genomes recovered from isolates, single cells, and
  metagenomes. *Genome Research* 25:1043–1055. doi:10.1101/gr.186072.114
- Parks, D. H., Chaumeil, P.-A., Mussig, A. J., Rinke, C., Chuvochina, M., Hugenholtz, P.
  (2026). GTDB release 10: a complete and systematic taxonomy for 715 230 bacterial and 17 245
  archaeal genomes. *Nucleic Acids Research* 54(D1):D743–D754. doi:10.1093/nar/gkaf1040
  (covers releases up to R10-RS226; R232 is release 11)
- Paysan-Lafosse, T., et al. (2025). The Pfam protein families database: embracing AI/ML.
  *Nucleic Acids Research* 53(D1):D523–D534. doi:10.1093/nar/gkae997
- Pereira, M. M., Santana, M., Teixeira, M. (2001). A novel scenario for the evolution of
  haem-copper oxygen reductases. *Biochimica et Biophysica Acta* 1505:185–208.
  doi:10.1016/S0005-2728(01)00169-4
- Porter, T. M., Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode
  classification. *Scientific Reports* 8:4226. doi:10.1038/s41598-018-22505-4
- Porter, T. M., Hajibabaei, M. (2022). MetaWorks: a flexible, scalable bioinformatic pipeline
  for high-throughput multi-marker biodiversity assessments. *PLoS ONE* 17:e0274260.
  doi:10.1371/journal.pone.0274260
- Rognes, T., et al. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ*
  4:e2584. doi:10.7717/peerj.2584
- Smith, M. A., et al. (2012). *Wolbachia* and DNA barcoding insects: patterns, potential, and
  problems. *PLoS ONE* 7:e36514. doi:10.1371/journal.pone.0036514
- Sousa, F. L., et al. (2011). A bioinformatics classifier and database for heme-copper oxygen
  reductases. *PLoS ONE* 6:e19117. doi:10.1371/journal.pone.0019117
- Sousa, F. L., et al. (2012). The superfamily of heme-copper oxygen reductases: types and
  evolutionary considerations. *Biochimica et Biophysica Acta* 1817(4):629–637.
  doi:10.1016/j.bbabio.2011.09.020
- Straub, D., et al. (2020). Interpretations of environmental microbial community studies are
  biased by the selected 16S rRNA (gene) amplicon sequencing pipeline. *Frontiers in
  Microbiology* 11:550420. doi:10.3389/fmicb.2020.550420
- Wang, Q., et al. (2007). Naïve Bayesian classifier for rapid assignment of rRNA sequences
  into the new bacterial taxonomy. *Applied and Environmental Microbiology* 73:5261–5267.
  doi:10.1128/AEM.00062-07
- Wangensteen, O. S., et al. (2018). DNA metabarcoding of littoral hard-bottom communities:
  high diversity and database gaps revealed by two molecular markers. *PeerJ* 6:e4705.
  doi:10.7717/peerj.4705
- Zafeiropoulos, H., et al. (2021). The Dark mAtteR iNvestigator (DARN) tool: getting to know
  the known unknowns in COI amplicon data. *Metabarcoding and Metagenomics* 5:e69657.
  doi:10.3897/mbmg.5.69657
- Web sources: BGE BioscanFlow COI chapter
  (https://bioscanflow.readthedocs.io/en/latest/6bioinfo_COI.html); MJOLNIR3
  (https://github.com/adriantich/MJOLNIR3); BOLDigger3
  (https://github.com/DominikBuchner/BOLDigger3); apscale-blast
  (https://github.com/TillMacher/apscale_blast); Pfam models via InterPro
  (https://www.ebi.ac.uk/interpro/); NCBI hmm_PGAP models (https://ftp.ncbi.nlm.nih.gov/hmm/).
