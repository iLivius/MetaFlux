# Bacterial COX1 reference set (COI)

COI primers do not amplify only animals. Most bacteria that respire oxygen carry genes
homologous to the mitochondrial cytochrome c oxidase subunit I, and common COI primer pairs
amplify some of them, those of *Wolbachia* above all. MIDORI2, the COI reference of MetaFlux,
contains only records filed under eukaryotes, so an amplicon of bacterial origin cannot be
recognised as such: it can only be assigned to a eukaryote, usually an animal.

MetaFlux therefore uses a second reference: a set of bacterial and archaeal subunit-I genes
taken from the species representatives of the Genome Taxonomy Database, GTDB (Parks et al.
2026). For every ASV, MetaFlux compares its best identity to this set (B) with its best
identity to MIDORI2 (E). An ASV with B > E is reported in a separate table of bacterial COX1
sequences. The same set is used once, when MIDORI2 is prepared, to remove MIDORI2 records
that are themselves bacterial COX1 filed under eukaryotic names.

!!! note "You do not build this set"
    MetaFlux downloads the released set and checks its SHA-256 sum. This page explains how
    the set was built, why each step is done as it is, how other tools deal with the same
    problem, and what the set cannot do.

All rules on this page were set on a **development sample** of 3,785 species representatives
of GTDB R226, the release before R232. The sample contains all 485 species of the order
Rickettsiales (56 of them *Wolbachia*), because this order holds *Wolbachia*, *Rickettsia* and
other symbionts of arthropods, plus 3,000 other bacteria and 300 archaea drawn per phylum (at
least three per phylum, the rest in proportion to the number of species in the phylum).
Because Rickettsiales are over-represented, shares given for "all bacteria" are weighted back
to the number of GTDB species in each group. The rules were fixed before the set was built
from GTDB R232, except the evidence rule of step 6 and its use in check 5.4, which were
revised during the build. The set was then tested once on data not used during development
(step 7).

## At a glance

| | |
|---|---|
| Release | `bacterial_cox1_gtdb-r232.v1` (Zenodo DOI: [To be settled before release]) |
| Source | GTDB R232: 199,923 species representatives, 189,801 bacterial and 10,122 archaeal |
| Genes | GTDB's own gene predictions, made with Prodigal 2.6.3 (Hyatt et al. 2010), file `gtdb_proteins_nt_reps_r232.tar.gz`; the build runs no gene prediction of its own |
| Content | subunit-I genes of the heme-copper oxidase superfamily, as DNA: whole genes in one file, fragments of at least 200 amino acids in a second file |
| Size | 242,698 whole genes (121 MB compressed) and 31,710 fragments (11 MB) from 123,356 genomes; counts per step in the release README |
| Taxonomy in the headers | GTDB R232, in SINTAX format (`d:`, `p:`, `c:`, `o:`, `f:`, `g:`, `s:`) |
| Used by MetaFlux | the B-versus-E comparison for every ASV, and the cleaning of MIDORI2 |
| Licence | CC BY-SA 4.0, as GTDB |
| Tools for the build | HMMER 3.4 (Eddy 2011) and VSEARCH 2.32.0 (Rognes et al. 2016) from bioconda, Python 3 from conda-forge |

## How MetaFlux uses the set

For a COI run, MetaFlux downloads the release and MIDORI2, removes from MIDORI2 the records
that step 6 marks as bacterial, and computes for every ASV its best identity to the cleaned
MIDORI2 (E) and to the set (B). Both identities are computed as described in step 5 and must
cover at least 80 % of the ASV.

- **B > E.** The ASV is moved from the main ASV and taxonomy tables to
  `bacterial_cox1_asvs.tsv`, with its read counts, B, E and its best bacterial
  match. In the read-tracking table, the column `post_bacterial_cox1` holds the reads left
  after the move; the moved reads are `post_frame_check` minus `post_bacterial_cox1`.
  Nothing is deleted, so that, for example, a *Wolbachia* infection stays visible.
- **All other ASVs.** E and the best MIDORI2 match are reported for every ASV, and ASVs
  with E below 80 % are flagged.

The 80 % flag marks ASVs that are far from everything in MIDORI2. In the development tests,
at most 5 % of genuine animal sequences whose genus was missing from MIDORI2 fell below 80 %,
and the flag alone marked 64–90 % of bacterial amplicons, depending on the primer pair. A
flagged ASV may therefore be a eukaryote of a group that MIDORI2 covers poorly, or a bacterium
that the set does not contain. Removing flagged ASVs would cost too many genuine eukaryotes
(see "What MetaFlux does, and why"), so they are only flagged.

The list of bacterial records belongs to MIDORI2 GB273. MetaFlux stops when the list names
records that are missing from the MIDORI2 file it prepares, as a list made for another
release would; a new MIDORI2 release needs a new list, made with step 6 of the build scripts.
Settings are described on the [COI marker page](coi.md).

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
     5.1 reading frame: DNA length a multiple of 3, no internal stop codon
     5.2 GC content: gene and contig within 0.15 of the genome
     5.3 mitochondrial match: < 90 % identical to MIDORI2 records that are certainly
         mitochondrial
     5.4 closer to eukaryotes: not closer to MIDORI2 (cleaned in step 6) than to the
         genes of other genomes
   flags only, the gene is kept: cross-phylum cluster; mitochondrion-like contig;
   duplicate copy
      │
6  cleaning of MIDORI2, with the genes that pass 5.1–5.3: a record ≥ 90 % identical to
   a gene is removed when the gene's relatives or the record's tryptophan codons say
   "bacterial" and neither says "mitochondrial"; other matches are listed
      │
7  final test on data not used to set the rules (results go into the README)
      │
8  release: whole genes and fragments (FASTA), a table of every candidate with the
   outcome of every check, the MIDORI2 list of step 6, README, scripts,
   SHA-256 sums  ──►  Zenodo
```

## Why a bacterial reference is needed

### Bacterial COX1 in COI libraries

Subunit I of the heme-copper oxygen reductases is the catalytic subunit of the terminal
oxidases of aerobic respiration. The enzymes form three families, A, B and C. Family A
contains the mitochondrial enzyme and the aa3-type oxidases of many bacteria, and is split
into subtypes A1 and A2 (Pereira et al. 2001). Heme-copper oxidases are widely distributed in
Bacteria and Archaea (Sousa et al. 2012), and many genomes carry several of them. The
bacterial genes are known as *ctaD* or *coxA* (aa3-type cytochrome c oxidase), *cyoB*
(bo3-type quinol oxidase), *qoxB* (aa3-600 quinol oxidase) and *ccoN* (cbb3-type oxidase,
family C).

COI is not a general marker for bacteria. *Wolbachia* is the exception: *coxA* is one of the
five loci of its multilocus sequence typing scheme (Baldo et al. 2006), and Nowak et al.
(2026) use insect COI barcode data to survey *Wolbachia* strains.

COI primers bind the bacterial genes because the primer sites lie in conserved parts of the
protein. Nowak et al. (2026, preprint) searched 59,318 bacterial RefSeq genomes: 36.8 % carry a
gene annotated as aa3-like *ctaD/coxA*, and in these the reverse primer jgHCO2198 of the Leray
primer pairs matches the gene perfectly in 79.8 % of genomes, BR2 in 70.6 %, but the forward
primer BF3 in only 1.9 %. In the development sample, among the genomes that carry a family-A
gene, the shares are almost the same: 78.0 %, 64.6 % and 1.4 % (weighted to GTDB).

A single primer that binds is not enough: a product forms only when both primers of a pair
bind. An in-silico PCR on the development sample therefore counted a genome as amplified by a
primer pair when both primers bind its DNA, facing each other, and the stretch between them has
the expected length (418, 313 and 205 bp for the three pairs, ± 30 bp). For the table, each
primer could have at most one mismatch, so that it counts the products likely to form in a
real PCR; the last two bases at its 3′ end had to match, and insertions or deletions were not
allowed. The tests of B > E below use every product with up to three mismatches, so that they
include more distant amplicons. The table gives the share of genomes amplified in each group.

Share of genomes that give a PCR product with each primer pair:

| Primer pair | *Wolbachia* (56 genomes) | other Rickettsiales (429 genomes) | all bacteria (weighted to GTDB) |
|---|---|---|---|
| BF3/BR2 | 93 % | 25 % | 3.2 % |
| Leray-XT | 70 % | 15 % | 0.2 % |
| fwhF2/fwhR2n | 0 % | 0.5 % | 0 % |

For example, BF3/BR2 gives a product in 52 of the 56 *Wolbachia* genomes, but in only about 3
of every 100 bacterial species. The last column is an estimate for all bacterial species in
GTDB: the sample holds every Rickettsiales species but only a few genomes of each other phylum,
so each group is counted in proportion to its number of GTDB species. The reverse primers bind
widely; the forward primer decides which bacteria amplify.

How much of a library is bacterial depends on the sample:

- **Water eDNA.** In samples from UK estuarine and coastal sites, 45.8 % of the Leray-XT reads
  had a bacterial best hit in NCBI nt, and so did 55 % of the reads of a fish-targeted, less
  degenerate COI primer pair (Collins et al. 2019).
- **Insects.** The best-documented source is *Wolbachia*. Its COI appeared in 35.1 % of 7,901
  BF3/BR2 libraries of single insects, at a median of 1.8 % of the reads where present (Nowak
  et al. 2026). Smith et al. (2012) found it in 0.16 % of 1.09 million Sanger barcode traces
  of insect specimens.

With a reference made only of eukaryotic records, these sequences cannot be recognised:

- **No informative top rank.** MIDORI2 has the same top rank (Eukaryota) in every record, so
  SINTAX (Edgar 2016) reports it with support 1.00 for any sequence, bacterial or not.
- **A Metazoa filter keeps them.** When the MIDORI2 headers are given a kingdom rank, taken
  from MIDORI2's own taxonomy files, SINTAX assigned 42–69 % of the bacterial amplicons of the
  in-silico PCR to Metazoa with a confidence of at least 0.8, and 94–100 % of the *Wolbachia*
  amplicons.
- **Pseudogene filters pass them.** These filters reject ASVs with a stop codon in frame
  under the invertebrate mitochondrial code, where only TAA and TAG are stops. A gene of a
  bacterium that uses the bacterial code (table 11), as almost all do, has no TAA or TAG
  inside its reading frame, so it passes.
- **MIDORI2 holds bacterial records itself.** MIDORI2 is built from the eukaryotic divisions
  of GenBank (Leray et al. 2022), so a bacterial sequence enters it when it was submitted under
  a eukaryote's name, usually an animal's. With the genes of the development sample alone, 132
  records of MIDORI2 release GB273 (built from GenBank release 273) were at least 90 %
  identical to a bacterial gene (step 6).

### How other tools handle it

COI tools handle non-target sequences in two ways.

**Bacteria in the reference.**

- **CO1 Classifier.** The CO1 Classifier (Porter & Hajibabaei 2018), the COI reference used
  with MetaWorks (Porter & Hajibabaei 2022), states that "bacterial outgroup sequences were
  added"; how they were chosen is not documented.
    - Version 5.1.0 holds 2,216,547 records, 262 of them bacterial. 229 of these are
      *Wolbachia* COI barcode-region records from GenBank.
    - It is distributed for the RDP classifier (Wang et al. 2007), and since version 5.0.0 also
      in SINTAX and BLAST formats.
    - With pseudogene filtering switched on, the MetaWorks COI workflow keeps only ESVs (exact
      sequence variants, MetaWorks' term for ASVs) assigned to the target group (Arthropoda in
      the shipped configuration), so bacterial ESVs are dropped together with the other
      non-target ESVs.
- **COInr.** COInr (Meglécz 2023) includes all cellular organisms, and its documentation
  states that the taxonomic origin of a sequence is taken at face value from the source
  database.
    - Its release of 3 May 2026 contains 9,775 bacterial sequences from 970 genera, 9,454 of
      them from BOLD.
    - VTAM and mkLTG point their users to ready-made files of the release of 6 May 2022, whose
      bacterial content was not counted for this page.
- **MJOLNIR.** The COI references of MJOLNIR contain some prokaryotes, mostly Rickettsiales.
  These are db_COI_MBPK (Wangensteen et al. 2018), built by in-silico PCR on EMBL and from
  BOLD, and the later DUFA, built from BOLD, NCBI and local sequences. In its final step,
  MJOLNIR3 removes by default the MOTUs (molecular operational taxonomic units, clusters of
  similar sequences) assigned to prokaryotes; in the original MJOLNIR this step is off by default.
- **DARN.** DARN (Zafeiropoulos et al. 2021) places COI sequences on a reference tree that
  includes 463 bacterial and 21 archaeal consensus sequences.
    - The bacterial ones were built from 3,917 bacterial COI sequences from BOLD, and from
      9,154 bacterial sequences whose DNA was extracted from genomes for the proteins that
      Pfam assigns to COX1 (PF00115). 217 archaeal sequences were obtained in the same way.
    - In the authors' test with mock data, a BLAST search against the consensus sequences
      assigned all bacterial sequences to eukaryotes, while placement on the tree recognised all
      of them.
- **A study-specific MIDORI.** Nowak et al. (2026) added alphaproteobacterial sequences to
  MIDORI release GB239 for their study.
- **BOLD.** Up to version 3, BOLD checked uploaded sequences against a library of
  contaminants that includes *Wolbachia* and flagged matching records (Smith et al. 2012). No
  description of this check was found for BOLD version 5.

**A eukaryote-oriented reference with similarity limits or taxonomic filters.**

- **VTAM.** VTAM (González et al. 2023) assigns the lowest taxonomic group that contains at
  least 90 % of the BLAST hits, trying identity steps from 100 % down to 70 %. Its
  documentation notes that when this group "could only be defined at 80% or lower identity
  level, the results are not very robust".
- **eKOI.** When re-analysing 15 published COI datasets against eKOI and MIDORI2, the eKOI
  authors (González-Miguéns et al. 2025) ignored ASVs below 84 % similarity to the
  references, a value they recommend applying with caution.
- **Anacapa.** The BLCA classifier of the Anacapa Toolkit (Curd et al. 2019) considers hits of
  at least 80 % identity.
- **BOLDigger3 and apscale-blast.** BOLDigger3 (version 3.0.3; Buchner & Leese 2020) and
  apscale-blast (version 2.0.6) limit the assigned rank by identity.
    - BOLDigger3: 97, 95, 90, 85 and 75 % from species to class, and 50 % for phylum.
    - apscale-blast: 98, 95, 90, 87 and 85 % from species to class.
- **nf-core/ampliseq.** nf-core/ampliseq 3.0.0 (Straub et al. 2020) offers MIDORI2 CO1 for its
  DADA2 classifier and for a VSEARCH consensus classifier (90 % identity, full query
  coverage). It adds an optional codon filter and a name-based exclusion that by default
  removes mitochondria and chloroplasts.
- **BGE Bioscan.** The BGE Bioscan workflow for BF3/BR2 classifies ASVs with SINTAX against
  BOLDistilled and keeps those assigned to Animalia, and to a class, with a confidence of at
  least 0.8.

With the tools of the second group, what happens to a bacterial ASV depends on its identity
to the nearest eukaryotic record. Below the lowest limit it stays unassigned; above that
limit it is assigned to a eukaryotic group, at the rank that its identity allows; and where
only the target group is kept, it is removed unless it was assigned to that group. No
published measurement of how often each case occurs for bacterial COX1 was found, nor a study
that compares the two approaches. The next section gives the measurements made for MetaFlux.

### What MetaFlux does, and why

MetaFlux compares each ASV with both references. The choice rests on three results of the
development tests, which used the following design.

- **Bacterial queries** were the in-silico PCR products of the development sample, split into
  two halves by species (and, in a second run, by genus). Each half served once as queries
  and once as the bacterial reference for the other half, so that no query could find its
  own species (or genus). MIDORI2 was cleaned with the reference half only.
- **Eukaryotic queries** were genuine COI sequences from MIDORI2, compared after every record
  of their own species, genus, family or order had been left out of MIDORI2, to mimic
  organisms that the reference lacks.

The bacterial reference in these tests was made of in-silico PCR products, not of whole
genes. The final test (step 7) repeats the measurements with the R232 gene set.

**First, no identity limit separates bacteria from eukaryotes.** Bacterial amplicons and
eukaryotic sequences whose relatives are missing from the reference fall in the same range
of identities to MIDORI2:

| Primer pair | Identity limit needed to catch 95 % of bacterial amplicons | Genuine animal sequences lost at that limit (family missing from MIDORI2) |
|---|---|---|
| BF3/BR2 | 83 % | 36 % |
| Leray-XT | 89 % | 72 % |
| fwhF2/fwhR2n | 95 % | 87 % |

**Second, the comparison B > E recognised 94 % of the bacterial amplicons**, pooled over the
three primer pairs: BF3/BR2 96 %, Leray-XT 91 %, fwhF2/fwhR2n 85 %, and *Wolbachia* 97 %. Of the
genuine animal sequences, it moved none to the bacterial side when their species or genus
was missing from MIDORI2, 0.04 % when their family was missing and 0.14 % when their order was
missing.

**Third, simpler variants are weaker or add nothing.**

- A one-sided screen that removes any ASV at least 78 % identical to a bacterial amplicon
  recognised 71 % of the bacterial amplicons. The limit of 78 % lies just above the highest
  identity of any genuine eukaryotic test sequence to the bacterial amplicons (77.6 %, once
  two genomes with eukaryotic contamination were left out of the bacterial reference), so it
  is the lowest limit that removed no eukaryote.
- Bacterial sequences added to the SINTAX reference pulled none of 4,560 eukaryotic test
  sequences into Bacteria, but added at most 0.6 percentage points to the share recognised by
  B > E. MetaFlux does not add them.

## The build, step by step

### Step 0. Source of the bacterial sequences

Three sources were considered.

1. **Bacterial records labelled COI in GenBank or BOLD,** as used by COInr and the CO1
   Classifier.
    - They are ready to use, but their labels come from the submitters, and a few are wrong:
      two BOLD records labelled Alphaproteobacteria are COI of flies (94–97 % identical to fly
      records, GC content 31–37 %).
    - Most records also come from a few lineages. *Wolbachia* makes up 229 of the 262
      bacterial records of the CO1 Classifier, and Actinomycetes and Alphaproteobacteria
      three quarters of those in COInr.
2. **In-silico PCR products cut from genomes.** They are close to what the primers amplify,
   but tied to the primer pairs used for the build; any other primer pair would need a new
   build.
3. **Whole genes predicted in genomes.** Every subunit-I gene of every genome, independent of
   primers, with the taxonomy of the genome. DARN used this source for part of its bacterial
   sequences.

MetaFlux uses the third source, taken from GTDB, for three reasons:

- all GTDB genomes are classified in a single taxonomy;
- each species has a representative genome;
- GTDB distributes its own gene predictions for every representative.

The build uses these predictions as they are and does not predict genes again, so every gene
in the set can be traced to GTDB's annotation of the release.

Only species representatives are used. Recognition does not need the query's own species: in
the development tests, with the species left out of the bacterial reference, 94 % of the
bacterial amplicons were still recognised, and 78 % with the genus left out. One genome per
species is therefore enough, and the other 700,000 genomes of the release would multiply the
size of the set with little gain.

The cost of taking genes from genomes is that many representatives are metagenome-assembled
genomes (MAGs). MAGs often cut genes at contig ends, and they can contain contigs of other
organisms, eukaryotes included. Step 4 deals with cut genes, step 5 with contamination.

### Step 1. Integrity and translation

The gene file and the metadata are checked against GTDB's published MD5 sums. For each
genome, the number of genes in the file must equal the protein count in GTDB's metadata; a
genome that fails is left out and listed.

The build reads the DNA genes and translates them itself, rather than also downloading
GTDB's protein file (another 123 GB). Protein and DNA then always come from the same record.
The price is that the genetic code has to be worked out for each genome:

- **The default.** Each gene is translated with the bacterial genetic code (NCBI table 11).
- **Table 4.** For some lineages, mainly Mycoplasmatales (GTDB phylum Bacillota_I; NCBI phylum
  Mycoplasmatota), GTDB predicts genes with table 4, in which TGA codes for tryptophan.
- **The table is not in the files.** The table GTDB used is not stored in the release files.
  The metadata give only NCBI's table, and for 92,894 of the 189,801 bacterial species
  representatives they give none.
- **The rule.** The table is therefore read back from the genes. Under table 11 TGA is a
  stop codon, so no gene predicted with table 11 can contain TGA inside its reading frame;
  under table 4, genes that encode tryptophan with TGA do. If more than 1 % of a genome's
  genes with both a start and a stop codon (`partial=00`) contain an in-frame TGA before
  their stop codon, the genome is translated with table 4.
- **Recording.** The inferred table is recorded for every gene and compared with NCBI's table
  where NCBI gives one; disagreements are listed in the README.

Two simpler rules were tested on the development sample and not used:

- **A rule by lineage** ("Mollicutes use table 4") would have been wrong for 11 genomes: some
  Acholeplasmatales use table 11, and some table-4 genomes lie outside the Mollicutes.
- **The coding-density rule of GTDB-Tk** disagreed with the table GTDB had used for 2 of the
  3,785 genomes.

The lineages concerned rarely carry respiratory oxidases. In the development sample, none of
the 15 genomes predicted with table 4 (10 Mycoplasmatales, 3 BD1-5, *Anadelfobacter* and
*Stammera*) carried a candidate gene. The README reports how many R232 genomes the 1 % rule
assigned to table 4.

### Step 2. Finding subunit-I genes

All proteins are searched with HMMER `hmmsearch` and the Pfam model PF00115 ("Cytochrome C and
Quinol oxidase polypeptide I", Pfam release 38.2; Paysan-Lafosse et al. 2025). The search uses
the model's gathering threshold, the score above which Pfam's curators count a protein as a
member of the family (34.8 bits for PF00115). PF00115 covers subunit I of all three oxidase
families and of the related nitric oxide reductases.

The Pfam models for COX2 (PF00116), COX3 (PF00510) and cytochrome b (PF00033) are searched
at the same time. Their hits are used only to describe the contig on which a gene lies
(step 5).

Two alternatives were tested on the development sample and not used for selection:

- **A sequence motif of family A** (`W.[FW]GHPEV[YF]....P`). It finds only subtype A1, and
  missed the 476 A2-type genes of the sample, among them every cyanobacterial aa3 oxidase.
- **The TIGRFAM models for the individual oxidase types.** Their trusted cutoffs are the
  scores above which the curators consider every hit a true member of that oxidase type. They
  are set high so that proteins can be named safely, and are too strict for selection: they
  left 758 whole family-A genes without a label, and these held 143 of the 1,239 in-silico
  PCR products that lay in kept genes. The TIGRFAM models are therefore used only for labels
  (step 3).

Gene annotations in RefSeq were not used either. They depend on the annotation pipeline:
NCBI attaches the name "cbb3-type cytochrome c oxidase subunit I" to its copy of the generic
PF00115 model, so genes recognised only by that model carry this name, family-B genes
included.

### Step 3. Oxidase or not, and which family

A protein can match PF00115 above the gathering threshold without being an oxidase subunit I.
In the development sample, 249 of the 4,304 genes that passed all other checks (5.8 %) failed
the oxidase rule below. 193 of them had two or more of the six ligand histidines replaced,
most often the CuB pair His290 and His291, and 56 scored below 100 bits. Most were 425–500
amino acids long, with a median score of 55 bits. Such genes are kept out of the set and listed
in the metadata table. None of the in-silico PCR products of the development sample lay in
them.

Each hit is aligned to PF00115 together with human COX1 (GenBank HQ713447.1). Human COX1 is
included only to give the alignment columns their usual residue numbers; at the positions
below, human and bovine numbering are identical.

| Bovine position | Role | Used for |
|---|---|---|
| His61, His378 | ligands of heme a | oxidase |
| His376 | ligand of heme a3 | oxidase |
| His240, His290, His291 | ligands of CuB | oxidase |
| Tyr244 | cross-linked to His240 (helix VI) | family A |
| Glu242 | proton channel, motif `XGHPEV` | subtype A1 |
| Tyr238, Ser239 | motif `YSHPXV`, replacing the A1 glutamate | subtype A2 |

- **Oxidase.** A hit counts as an oxidase subunit I if it scores at least 100 bits against
  PF00115 and at most one of the six histidines is replaced. One replacement is allowed, and a
  gap in the alignment is not counted as a replacement, because a single alignment to the
  Pfam model can misplace a residue in a divergent sequence.
- **Family A.** A gene is assigned to family A if it carries Tyr244 together with Glu242
  (subtype A1) or with Tyr238 and Ser239 (subtype A2); these are the residues that define the
  subtypes (Pereira et al. 2001).
- **Type label.** The label (aa3 or caa3, bo3, aa3-600, cbb3) comes from the best model above
  its trusted cutoff among TIGR02891, TIGR02843, TIGR02882 and TIGR00780 (TIGRFAMs) and
  NF050116.1 (NCBIfam, aa3 oxidases of Actinomycetota), all from NCBI release hmm_PGAP 20.0
  (Li et al. 2021).

Family B has no HMM in Pfam or TIGRFAM, and family C carries its cross-linked tyrosine on
helix VII instead of helix VI (Sousa et al. 2011). Genes of both families therefore remain
"family not assigned" unless a TIGRFAM model names them. Nitric oxide reductases share the six
histidine positions but lack the cross-linked tyrosine (Sousa et al. 2011), and under these
rules they also end up as "family not assigned". They do no harm on the bacterial side of the
comparison: a eukaryotic COI amplicon is far less similar to them than to MIDORI2, so they
cannot raise B above E.

All oxidase subunit-I genes are kept, not only family A:

- **Why.** In the development sample, family A held every in-silico PCR product of the three
  primer pairs. Other primer pairs, or bacteria absent from the sample, may amplify genes of
  families B or C, and a set limited to family A would miss them.
- **The cost.** A set of family A only would be about 40 % as large.
- **No risk.** The other genes cannot cause a false B > E, for the reason given above.

### Step 4. Whole genes and fragments

In MAGs, genes are often cut at contig ends. A gene counts as whole if it meets three
conditions:

- it has a start and a stop codon (`partial=00` in GTDB's prediction);
- it covers at least 70 % of the PF00115 model, counted over all its domain hits together, so
  that a protein whose alignment to the model is broken in two, for example by an insertion,
  is measured as one;
- it is 400–1,100 amino acids long. Whole genes in the development sample had a median
  length of 555 amino acids. The range includes the fusions of subunits I and III found in
  some haloarchaea, *Thermus* and *Deinococcus* (517–901 amino acids; 43 genes in the
  development sample) and two haloarchaeal genes of 1,007 and 1,065 amino acids. No candidate
  was longer; the upper limit is there to catch gene predictions that run on into a
  neighbouring gene.

A gene that is not whole is kept as a fragment if it has at least 200 amino acids.

- **Why fragments are kept.** The primers amplify only 200–420 bp of a gene of about
  1,600 bp, so a fragment that covers that region works as well as a whole gene in the
  comparison. In the development sample, 5.9 % of the in-silico PCR products lay only in
  fragments, mostly from MAGs of uncultured lineages.
- **Why they are kept apart.** A fragment shows only part of the residues of step 3, so the
  evidence that it is an oxidase is weaker, and a gene cut in two would be counted twice in
  copy numbers.
- **Why 200 amino acids.** Shorter pieces carry too few residues and too little sequence for
  the checks of step 5; the shortest fragment that held an in-silico PCR product had 253
  amino acids.

Fragments go into a separate file. MetaFlux uses both files to compute B.

Whole genes and fragments together held 1,239 of the 1,259 in-silico PCR products of the
development sample (98.4 %). Of the other 20, five came from contigs of eukaryotic origin and
were rightly left out, five lay across a gene split in two and were partly kept, and ten lay
at loci where the gene is broken by a frameshift or an internal stop codon.

### Step 5. Contamination checks

GTDB admits a genome only if CheckM (Parks et al. 2015), which uses marker genes, and CheckM2
(Chklovski et al. 2023), a machine-learning model, both estimate completeness of at least
50 %, contamination below 5 % and a quality score of at least 50 (Parks et al. 2026). Among
the genomes that pass, it prefers type material and high-quality assemblies as species
representatives. A single foreign contig carrying one gene is below what such genome-wide
estimates can detect, so every candidate gene is tested individually.

Two other approaches are possible:

- **Leaving out MAGs** would remove most contamination. But for many uncultured lineages a MAG
  is the only genome, and those lineages hold most of the in-silico PCR products found only in
  fragments (step 4).
- **Classifying each contig taxonomically** would need the genome sequences, another archive
  of about 180 GB that the build does not download.

Without the genome sequences, the GC content and length of a contig are estimated from the
genes on it. In the development sample this gave the same GC verdict as the contig sequence
for all 5,060 candidate genes, and the same contig-length flag for all but one.

Four checks remove a gene.

- **5.1 Reading frame.** The DNA length is a multiple of 3, and there is no internal stop
  codon under the table chosen in step 1.
- **5.2 GC content.** The GC content of the gene, and that of its contig, lies within 0.15 of
  the genome's GC content.
    - The largest difference observed for a genuine gene in the development sample was 0.132.
    - Example: a thraustochytrid-like COX1 in a planctomycete MAG (gene GC 0.39, genome GC 0.66)
      fails this check. Check 5.3 cannot catch it, because thraustochytrids read TGA as a stop,
      so none of their records count as certainly mitochondrial; check 5.4 removes it as well
      (80.1 % identical to *Aurantiochytrium*).
- **5.3 Match to mitochondrial COI.** The gene is less than 90 % identical, over at least 80 %
  of the shorter sequence, to every MIDORI2 record that is certainly mitochondrial.
    - **Why this works.** Animal and fungal mitochondria read TGA as tryptophan, and most of
      them use it often, while a gene of a bacterium that uses table 11 cannot contain TGA in
      frame. Lineages that use TGA rarely, such as mushrooms and liver flukes (step 6), have few
      certainly mitochondrial records, so a COI of theirs binned into a MAG is left to step 6
      and check 5.4.
    - **Why 90 %, over the shorter sequence.** Genuine bacterial genes of the development sample
      reached at most 85.4 % (a *Tisiphia* gene, against a record that is probably chimeric),
      while the copepod COI found in a bacterial MAG (see "Contamination removed before any
      check" below) is 99.5 % identical to a mitochondrial record. Coverage is measured on the
      shorter sequence because many MIDORI2 records are barcodes of about 650 bp, which cannot
      cover 80 % of a whole gene.
    - **The definition.** A record counts as certainly mitochondrial if none of its three
      reading frames is open under the bacterial code (each has at least two TAA, TAG or TGA),
      its reading frame contains at least two TGA codons, and TGA is used at least as often as
      TGG.
    - **Why it is strict.** MIDORI2 contains bacterial records. A looser definition admitted two
      bacterial records filed as a fish and a planthopper, which then caused 19 genuine bacterial
      genes to be rejected.
    - **How many records qualify.** MIDORI2 GB273 is used in its UNIQ NUC_SP variant: every
      unique sequence of each species, with records whose species name is incomplete, such as
      "sp." or "cf.", included. Of its 3,213,682 records, 2,904,469 meet the definition.
- **5.4 Closer to eukaryotes than to bacteria.** The gene's best identity to MIDORI2 is at
  least 80 % and higher than its best identity to the genes of the other GTDB genomes. This
  check is run after step 6, on MIDORI2 without the records removed there, so that a
  bacterial record filed under an animal's name cannot make a genuine bacterial gene look
  eukaryotic.
    - **Records whose codons look bacterial do not count.** For the same reason, a MIDORI2
      record whose own TGA/TGG counts look bacterial by the codon test of step 6 is not used
      as evidence here. Most such records are less than 90 % identical to every gene, so
      step 6 never examined them, but they are close enough to the genes of their own
      bacteria to make those genes look eukaryotic. In R232, 58 records were left out of E
      this way, and 20 genes that they would have made check 5.4 reject passed.
    - **Outcome.** Such a gene is removed if it is a fragment or if its GC content differs from
      the genome's by more than 0.08. Otherwise it is marked "review" and also left out of the
      release, so that no manual decision enters the build.
    - **In the development sample** this check failed three candidates and marked 20 whole genes
      for review:
        - the thraustochytrid gene of check 5.2;
        - a fungal mitochondrial piece in a *Rickettsia* MAG (89.3 % identical to *Phialocephala*,
          just below the limit of check 5.3). At 55 amino acids it also falls below the floor of
          step 4;
        - a 274-amino-acid *Polaromonas* fragment, probably removed wrongly. Its best MIDORI2
          match, MF319222.1, is filed as the clam *Corbicula fluminea* but has a GC content of 0.64
          and one TGA in frame, and is probably a bacterial record itself.

**How identity is computed.** Identities in checks 5.3 and 5.4, in step 6 and in MetaFlux
itself are computed from the VSEARCH alignment:

- identical A/C/G/T columns are divided by all aligned columns;
- end gaps are left out, and internal gaps count as differences;
- columns where either sequence has N or another ambiguity code are left out of both counts.

VSEARCH's own identity counts a column as identical whenever the two symbols can stand for
the same base, so an N matches anything. Its option `--n_mismatch` changes this for N but not
for the other ambiguity codes. About 9 % of MIDORI2 records contain ambiguity codes, some in
long runs, and a sequence aligned over such a run would otherwise appear identical. In the
development tests, VSEARCH's own value was inflated in this way for 8–15 % of the bacterial
amplicons; one archaeal amplicon appeared 100 % identical to an insect record, while its best
identity over A/C/G/T columns was 52.8 %.

Three observations are recorded as flags but do not remove a gene:

- **Cross-phylum cluster.** The gene clusters at 95 % identity with genes from another phylum,
  as the same contaminant binned into MAGs of different lineages would. There were none in the
  development sample. In R232, 524 kept genes carry the flag, in 182 clusters, 154 of them
  made only of MAGs or single-cell genomes; most join a Pseudomonadota gene to one of another
  phylum. Such a gene is still bacterial, so B is not affected, but the bacterial lineage
  MetaFlux reports for an ASV that matches it can be the wrong one.
- **Mitochondrion-like contig.** The gene lies on a contig shorter than 20 kb that also
  carries cytochrome b and COX2 or COX3, as a mitochondrial contig would. Bacterial operons in
  which the bc1 complex lies next to the oxidase genes trigger this flag too (37 genes in the
  development sample, all bacterial; 2,247 in R232), which is why it is not a filter. The one contig of
  eukaryotic origin in the development sample, which carries the thraustochytrid-like gene of
  check 5.2 together with COX2, COX3 and cytochrome b, is 29.5 kb long and is not flagged.
- **Duplicate copy.** Two whole genes of the same genome are at least 99 % identical, as
  happens when a contig is duplicated in an assembly.

**Contamination removed before any check.** An animal or fungal mitochondrial COI binned into
a bacterial MAG contains TGA codons, which the gene prediction reads as stops under table 11.
At most, short open reading frames are predicted at such a locus, and they fall below the
200-amino-acid floor of step 4. In the development sample, no gene was predicted at a copepod
COI in a bacterial MAG, and a fungal COX1 in a *Rickettsia* MAG gave only a 55-amino-acid
open reading frame.

### Step 6. Cleaning MIDORI2

The genes that pass checks 5.1–5.3 are compared with every record of MIDORI2 GB273 (UNIQ
NUC_SP). The NUC_SP variant contains every record of the NUC variant, so the resulting list
applies to both.

**A close match means that one side is misfiled.** Genuine COI and genuine bacterial COX1 are
never this close: in the development tests, no genuine eukaryotic sequence was more than
77.6 % identical to a bacterial gene. A record at least 90 % identical to a gene (identity as
in step 5, over at least 80 % of the shorter sequence) is therefore either a bacterial sequence
filed under a eukaryote's name, or the gene is a eukaryotic sequence binned into a bacterial
genome. Two independent lines of evidence say which.

- **The gene's side: are its relatives bacterial or mitochondrial?** A genuine bacterial gene
  usually has relatives in other GTDB genomes; a mitochondrial COI in a MAG has mitochondrial
  relatives instead.
    - The gene looks bacterial when its closest gene in another genome is at least 80 %
      identical and at least as close as its closest certainly mitochondrial record (step 5.3).
    - It looks mitochondrial when a certainly mitochondrial record is at least 80 % identical
      and closer than any gene of another genome.
    - Otherwise it gives no evidence. Only certainly mitochondrial records are used on this
      side, because they cannot be bacterial: they carry TGA codons in frame.
- **The record's side: do its tryptophan codons fit its lineage?** In a bacterial gene every
  tryptophan is TGG, because TGA is a stop. Most animal and fungal mitochondria read TGA as
  tryptophan, but how often they use it instead of TGG differs widely between lineages. In the
  insect records of MIDORI2, 94 % of tryptophan codons are TGA; in the mushrooms
  (Agaricomycetes), 21 %; in the liver-fluke family Fasciolidae, 15 %. A record without TGA is
  therefore strong evidence in an insect and none in a mushroom.
    - The record's TGA and TGG counts are compared under two explanations. As a bacterial gene,
      a TGA can only be an error, taken to occur at 1 in 100 tryptophan codons. As a genuine
      COI of its lineage, each tryptophan is TGA with the share seen in the lineage's other
      MIDORI2 records: those of its family, or of the lowest higher rank (order, class,
      phylum) whose records hold at least 200 tryptophan codons. Records in conflict with a
      gene are left out of these shares.
    - The ratio of the two likelihoods counts as evidence when it is at least 20 in either
      direction, the level that Kass and Raftery (1995) call strong.
    - In lineages whose mitochondria read TGA as a stop (plants, brown algae, oomycetes,
      diatoms, chrysophytes, and green algae such as *Micromonas*), the codons give no evidence.

**The decision.** A record is removed when at least one line of evidence says bacterial and
none says mitochondrial. Otherwise it is listed and kept, and check 5.4 then judges the gene
with the record still in MIDORI2. Examples from the R232 build:

- **Kept: six mushroom records.** Six records of *Rhizoctonia*, *Ganoderma* and *Fomes* were
  91–98 % identical to a gene in a MAG of a soil archaeon (Nitrososphaeraceae). Their lineages
  use TGA for only 1–5 % of tryptophan codons, so their TGA counts are no evidence. The gene has
  no relative in any other genome. The COX1 model matches only amino acids 12–151 of its 531,
  and the rest carries the two LAGLIDADG motifs of the homing endonucleases encoded in
  group I introns (Pfam PF00961): the gene is the first exon of a fungal COX1 read on into the
  open reading frame of its intron, a common arrangement in fungal mitochondria. The closest
  record, KC352446.1 from a *Rhizoctonia solani* mitochondrial genome, is that intron reading
  frame alone, filed in MIDORI2 as COI. The records stay, and check 5.4 removes the gene.
- **Removed: a beetle record with one TGA.** A record filed as the carrion beetle *Calosilpha
  brunnicollis* has one TGA among seven tryptophan codons and is 96.8 % identical to a gene of
  *Aquirickettsiella*, a bacterial pathogen of arthropods. Its family uses TGA for 93 % of
  tryptophan codons, so one TGA in seven is about 60,000 times likelier for a bacterial gene
  with one sequencing error than for the beetle's own COI.
- **Removed: three chrysophyte records.** Three records of *Dinobryon pediforme* are 99.3–99.6 %
  identical to a gene of a MAG of the Paracaedibacterales, an order of alphaproteobacteria that
  includes endosymbionts of protists. Chrysophytes read TGA as a stop, so the codons say nothing, but the
  gene has relatives in four other genomes of the order (75–86 %).
- **Removed: a liver-fluke record.** A record filed as *Fasciola hepatica* is 96.0 % identical
  to the gene of *Neorickettsia findlayensis*, an endosymbiont of trematodes whose complete
  genome has relatives in two other *Neorickettsia* genomes (89–90 %). Liver flukes rarely use
  TGA, so the codons say nothing either way, and the gene decides.

In R232, 186 records were at least 90 % identical to one of 95 genes that passed check 5.3,
and 170 of them were removed:

- **By lineage.** 137 of the removed records are filed under arthropods, 13 under diatoms,
  6 under chordates, and 14 under seven other groups (gastrotrichs, chrysophytes, red
  algae, annelids, molluscs, a dinoflagellate and a liver fluke). The 16 listed records are
  the six mushroom records, four *Micromonas* records, four nematode records (*Plectus*,
  whose gene looks mitochondrial), and one record each of a synurophyte and a gastrotrich.
- **By evidence.** 128 records were removed on both lines of evidence, 40 on the gene's
  alone and 2 on the codons alone. Of the 40, 13 belong to lineages whose mitochondria read
  TGA as a stop; the other 27 have too few tryptophan codons, or belong to lineages that use
  TGA too rarely, for the codons to tell. 12 records were listed for lack of evidence, 4 because their
  gene looks mitochondrial.
- **Genes.** Of the 95 genes in conflict with a record, 89 look bacterial, 1 mitochondrial,
  and 5 give no evidence (no relative of at least 80 % on either side).

**Why two lines of evidence.** The first version of this step, fixed before the R232 build,
removed a record when it had no TGA in frame and its lineage reads TGA as tryptophan, and kept
every record with at least one TGA. Applied to R232, it went wrong in two ways:

- it removed five of the six mushroom records above, those without a TGA, because a missing
  TGA is no evidence in lineages that rarely use TGA, and the gene they matched was never
  examined;
- it kept bacterial records with a single TGA, such as the beetle record above. Left in
  MIDORI2, such records made check 5.4 reject the genes of the bacteria they come from: a
  single record filed as the clam *Corbicula fluminea* (MF319222.1: one TGA among six
  tryptophan codons, 97.9 % identical to a *Paenacidovorax* gene) made 22 genes of
  Burkholderiales fail.

Under the first rule, check 5.4 left out 102 genes. With the revised rule, and with the
records whose codons look bacterial left out of E (check 5.4), it leaves out 24; the other 78
pass, and no gene fails that passed before. Of the 24, 16 are close to the COI of green
algae, cryptophytes, diatoms, oomycetes, a malawimonad, a fungus or a nematode and are most
likely eukaryotic contaminants; the other 8 are close to records filed as animals (a scallop,
a gastropod, an amphipod, a shipworm) whose codons give no clear evidence, and are probably
genuine bacterial genes left out.

The rule was therefore revised before the final test. Because it was set while looking at the
R232 records, it has not been tested on independent data. The README reports how many
decisions change when each of its settings moves: none for an error rate of 0.003 or 0.03,
a likelihood-ratio threshold of 10 or 50, or 100 or 500 tryptophan codons per lineage; 2
and 1 of the 186 when the 80 % of the gene's side moves to 78 or 82 %. Examining records
from 88 % instead of 90 % would remove 6 of the 14 records between the two values; from 92 %,
15 records would be decided differently or no longer examined.

MIDORI2 already excludes 2,215 GenBank records that an earlier study had found mislabelled at
the order, class or phylum level (Leray et al. 2019, 2022). Bacterial sequences submitted
under the host's name are a separate problem, which this step addresses.

### Step 7. Final test

The rules were fixed before the R232 build, except the evidence rule of step 6 and its use
in check 5.4, which were revised during the build; the final test came after that revision.
It is run once, before release, on genomes and MIDORI2 records not used to set any rule, and
it asks how the set handles bacteria and eukaryotes that the references do not contain.

#### Bacterial amplicons

- **Test genomes.** 3,476 R232 species representatives that are not in the development
  sample, chosen with the development design and a new random seed: all 176 such
  Rickettsiales (4 of them *Wolbachia*), 3,000 other bacteria drawn per phylum, and 300
  archaea. All the Rickettsiales are taken, as in development, because the order holds
  *Wolbachia*, *Rickettsia* and other symbionts of arthropods.
- **In-silico PCR.** The rules of development: both primers bind, facing each other, the
  insert between them is within 30 bp of the expected length, a primer may have up to three
  mismatches but none in its last two bases at the 3′ end, and insertions and deletions are
  not allowed. Three mismatches is permissive, on purpose: the test gets more amplicons, and
  more distant ones. The products with an intact reading frame under the bacterial code are
  the queries. The PCR is run on the genes that GTDB predicted in each genome, not on the
  genome sequence, which the build does not download: an amplicon lies inside the COX1 gene,
  and in development this missed only the products in genes broken by a frameshift or an
  internal stop (10 of 1,259, 0.8 %).
- **Why most genomes give no product.** 738 of the 3,476 genomes gave a product, 936
  products in all: 803 of BF3/BR2, 102 of Leray-XT and 31 of fwhF2/fwhR2n. A genome gives one
  product per primer pair, or more when it carries several COX1 genes.
    - Almost all products come from family-A genes, the family of the mitochondrial enzyme
      (928 of the 936 from subtype A1), and 1,648 of the test genomes (47 %) have no family-A
      gene in the set.
    - Where there is one, the primers decide. BF3/BR2 gives a product in a fifth of the
      other bacteria (610 of 3,000), but in only 86 (2.9 %) with at most one mismatch per
      primer, close to the 3.2 % of the development table. Leray-XT and fwhF2/fwhR2n give
      one in hardly any bacteria outside the Rickettsiales.
    - With at most one mismatch per primer, the products likely to form in a real PCR, 181
      of the 936 remain.

Genomes giving a product, with up to three mismatches per primer (in brackets: at most one):

| Genomes | BF3/BR2 | Leray-XT | fwhF2/fwhR2n |
|---|---|---|---|
| *Wolbachia* (4) | 4 (4) | 4 (4) | 3 (0) |
| other Rickettsiales (172) | 72 (45) | 73 (31) | 25 (0) |
| other bacteria (3,000) | 610 (86) | 17 (3) | 3 (0) |
| archaea (300) | 7 (1) | 0 | 0 |

- **Each amplicon is tested with its own species left out of the set.** Every test genome
  is a GTDB species representative, so its own genes are in the set: its amplicon would find
  itself at 100 % identity and pass without telling anything. A bacterium in a real library
  is rarely the genome that GTDB holds; often its species is missing, sometimes its genus. B
  is therefore computed without the genes of the amplicon's own species, and, in a second
  run, without those of its whole genus, as for a bacterium of a genus absent from GTDB. The
  search keeps up to 3,000 hits per amplicon, and the hits to the genes left out are
  skipped.
- **The same on the MIDORI2 side.** Step 6 removes from MIDORI2 the bacterial records filed
  under eukaryote names, which it finds through the genes of the set. A removed record whose
  closest gene belongs to the amplicon's own species (or genus) is put back for that
  amplicon, as if that species were missing from GTDB. This is the cautious choice: through
  the genes of other species, step 6 might still have removed the record. E is computed
  against MIDORI2 with these records back.
- **Measured** is the share of amplicons with B > E: those that MetaFlux would recognise as
  bacterial.

#### Eukaryotic sequences

- **Test records.** 3,336 MIDORI2 GB273 records not used during development: 1,000 at the
  species level, 1,000 at the genus level, 858 at the family level and 478 at the order
  level, up to 250 taxa per group (arthropods, chordates, other animals, other eukaryotes).
  At the order level, and for the non-animal families, development had already used every
  eligible taxon, so these taxa are tested again with records not used before.
- **Each insert is tested with its own taxon left out of MIDORI2.** Each record is cut to
  the inserts of the three primer pairs (10,008 inserts). For each insert, every record of
  its own species (genus, family or order) is left out of MIDORI2, as for an animal that
  MIDORI2 does not contain, and B is computed against the whole set.
- **Measured** is the share with B > E: genuine eukaryotic sequences that MetaFlux would
  wrongly move to the table of bacterial COX1.

#### Known contaminants

As a check that nothing was lost, the contaminants found during development must still be
rejected. Unlike the other tests, this one uses development data.

#### Release criteria and results

The release is published if the results are close to those of development, by criteria
fixed before the test:

- for the bacterial amplicons with the species left out, a recognised share of at least 90 %
  for BF3/BR2, 85 % for Leray-XT and 80 % for fwhF2/fwhR2n;
- for the eukaryotic sequences, at most 0.5 % moved to the bacterial side at any of the four
  levels.

All criteria were met; the full table is in the README. Each value below is the share of
bacterial amplicons with B > E, with the amplicon's species or genus left out as described
above, next to the same measurement in development (with the species left out).

| Bacterial amplicons recognised (B > E) | species left out | genus left out | development (species left out) |
|---|---|---|---|
| all (936) | 98.4 % | 96.0 % | 94 % |
| BF3/BR2 (803) | 99.9 % | 98.8 % | 96 % |
| Leray-XT (102) | 88.2 % | 80.4 % | 91 % |
| fwhF2/fwhR2n (31) | 93.5 % | 77.4 % | 85 % |
| *Wolbachia* (11) | 72.7 % | 0 % | 97 % |
| other Rickettsiales (178) | 93.3 % | 87.1 % | |
| at most 1 mismatch per primer (181) | 95.6 % | 89.5 % | |

- **Eukaryotic sequences moved to the bacterial side:** 0.03 % with the species left out of
  MIDORI2 (1 insert of 3,000), 0 % with the genus left out, 0.08 % with the family left out
  (2 of 2,574) and 0 % with the order left out. All three are far from both references
  (E 59–71 %, B 60–73 %).
- **Known contaminants:** all three stay out of the set.
- **The *Wolbachia* figures rest on 11 amplicons from four genomes.** The three that were
  missed come from one genome, *Wolbachia* sp964277405. Its own COX1 is in MIDORI2 as a
  record filed under the pseudoscorpion *Neobisium carcinoides*; step 6 removes that record
  in the release, but the test puts it back, as described above. With the species left out,
  the closest bacterial gene is another *Wolbachia* at 90–92 % (B), while the record matches
  at 98–99 % (E), so E wins. The BF3/BR2 amplicon of *Wolbachia* sp051710285, by contrast,
  has B 97 % and E 85 %. With the genus left out, every *Wolbachia* record of MIDORI2 comes
  back and no *Wolbachia* amplicon wins, as in development (0.7 %): a *Wolbachia* lineage
  missing from the set is recognised only if MIDORI2 holds no record of it under an animal's
  name.

Any rule changed after this test makes a new version.

### Step 8. Release files

| File | Content |
|---|---|
| `bacterial_cox1_gtdb-r232.v1.fasta.gz` | whole genes that pass all checks; DNA; SINTAX headers |
| `bacterial_cox1_gtdb-r232.v1_fragments.fasta.gz` | fragments of at least 200 amino acids that pass all checks |
| `bacterial_cox1_gtdb-r232.v1_metadata.tsv.gz` | every candidate gene: genome, lineage, coordinates, scores, residues, family, type, the outcome of every check, flags and verdict |
| `midori2_gb273_bacterial_records.tsv` | MIDORI2 records matching a gene of the set, with the matching gene and the decision (removed or listed) |
| `README.md` | inputs and their checksums, tool and model versions, thresholds, counts per step, results of the final test, sensitivity of the counts to each threshold |
| build scripts and environment file | to rebuild the release from the same inputs |
| `SHA256SUMS` | checksums of all files |

Each FASTA header carries the gene ID, the GTDB accession, the family and type label and the
GTDB taxonomy, for example:

```text
>RS_GCF_000008385.1|NC_006833.1_429;family=A1;type=aa3_caa3;tax=d:Bacteria,p:Pseudomonadota,c:Alphaproteobacteria,o:Rickettsiales,f:Anaplasmataceae,g:Wolbachia,s:Wolbachia_sp000008385
```

The release is named after the GTDB release and a build number (`gtdb-r232.v1`). A later GTDB
release, or any change of a rule, produces a new version under the same Zenodo concept DOI.
MetaFlux pins the version DOI and the SHA-256 of the files it downloads. If you use the set,
cite the release by its version DOI, together with GTDB (Parks et al. 2026) and MIDORI2
(Leray et al. 2022).

### Parameters of version 1

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
| 5.4 | closer to eukaryotes | run after step 6, on the cleaned MIDORI2, without the records whose codons look bacterial (step 6 test); best MIDORI2 identity ≥ 80 % and above the best identity to other genomes: removed if fragment or GC difference > 0.08, otherwise "review" (left out) |
| 5 | flags | cross-phylum cluster at 95 %; mitochondrion-like contig < 20 kb; duplicate copy ≥ 99 % |
| 5–6 | identity | identical A/C/G/T columns / (A/C/G/T columns + internal gap columns); end gaps and ambiguity codes left out; VSEARCH 2.32.0 `--usearch_global`, settings in the README |
| 6 | MIDORI2 record examined | ≥ 90 % to a gene that passes 5.1–5.3 |
| 6 | gene's evidence | bacterial: closest gene of another genome ≥ 80 % and at least as close as the closest certainly mitochondrial record; mitochondrial: certainly mitochondrial record ≥ 80 % and closer; otherwise none |
| 6 | record's evidence | likelihood ratio of its TGA/TGG counts, bacterial gene (TGA error rate 0.01 per tryptophan codon) vs genuine COI (TGA share of its family, or of the lowest higher rank with ≥ 200 tryptophan codons); ≥ 20 bacterial, ≤ 1/20 mitochondrial; none where TGA is a mitochondrial stop |
| 6 | MIDORI2 record removed | at least one line of evidence bacterial and none mitochondrial; otherwise listed and kept |
| 7 | release criteria | species left out: ≥ 90 % (BF3/BR2), ≥ 85 % (Leray-XT), ≥ 80 % (fwhF2/fwhR2n) recognised; eukaryotes: ≤ 0.5 % moved at each level |
| README | sensitivity | verdict changes when each value moves: GC limits ± 0.02, fragment floor ± 50 amino acids, identity limits ± 2 points, model coverage ± 0.10, oxidase score ± 20 bits; step 6: error rate 0.003 and 0.03, likelihood-ratio threshold 10 and 50, 100 and 500 tryptophan codons per lineage |

## Copy number

How many subunit-I genes a genome carries depends on what is counted. In near-complete
genomes of the development sample (CheckM2 completeness ≥ 90 %, contamination ≤ 5 %; shares
weighted to GTDB), with the counts of the 120,230 near-complete bacterial genomes of R232 in
brackets:

- **Family A is mostly single-copy.** 65 % of the genomes that carry it have one copy and 35 %
  have two to five [R232: 64 % of 71,669 carriers have one copy, up to nine]. All 51
  near-complete *Wolbachia* genomes that carry it have one copy, and so do 99 % of the other
  Rickettsiales [R232: 55 of 56 *Wolbachia* genomes have one copy, one has none].
- **All subunit-I genes together are mostly multi-copy.** 69 % of carriers have two or more,
  up to nine, because many bacteria carry several oxidase types [R232: 67 % of 79,030
  carriers, up to twelve]. In 90 % of the multi-copy
  genomes all copies are less than 80 % identical to each other; they are distinct genes, not
  assembly artefacts.

All copies that pass the checks are kept. The primers can bind any of them, and MetaFlux uses
only the best match of each ASV, so extra copies give a bacterial ASV more chances to find a
close relative without biasing the result.

## Limitations, open issues and possible improvements

**Coverage of bacterial diversity.** The set contains what GTDB contains. In the final test,
96.0 % of the amplicons were recognised with their genus left out of the set, against 98.4 %
with only their species left out (78 % and 94 % in the development tests). Most of them were
BF3/BR2 amplicons (803 of 936, 98.8 % recognised with the genus left out); for Leray-XT the
share was 80.4 %, for fwhF2/fwhR2n 77.4 %. None of the 11 *Wolbachia* amplicons was
recognised with its genus left out, because MIDORI2 holds *Wolbachia* COX1 under animal
names: a *Wolbachia* lineage missing from the set is recognised only if MIDORI2 holds no
record of it under an animal's name (step 7). Lineages without cultured or assembled members
are not represented.

**One genome per species.** Variation within a species is represented by one genome. The
development tests removed the whole species and still recognised 94 % of the amplicons, so an
amplicon from a divergent strain should be recognised at least as often; this was not tested
directly.

**Rules set on development data.** These values were chosen after inspection of the
development sample:

- the residue rule of step 3;
- the 200-amino-acid floor for fragments;
- check 5.4, with its GC limit of 0.08;
- the strict definition of certainly mitochondrial records.

They are frozen for version 1 and tested once on data not used to set them. Other values are
defensible; the README reports how many genes change verdict when each of them moves. The
evidence rule of step 6 is the exception: it was revised during the R232 build, after the first
rule had been applied to the R232 records, so no independent data have tested it (step 6).

**Contaminants that can pass.** The mitochondrial-match check (5.3) relies on TGA being read as
tryptophan, so it cannot recognise contamination from plants, brown algae, oomycetes, diatoms
or thraustochytrids. For these, only the GC check (5.2) and the closer-to-eukaryotes check (5.4)
apply. A contaminant with a GC content close to that of the genome and no close relative in
MIDORI2 would pass. Check 5.4 also misses a contaminant that occurs in more than one MAG, for
example in MAGs from the same study, because each copy then finds a close match in another
genome.

**Bacterial records left in MIDORI2.** Step 6 removes a record only when it is at least 90 %
identical to a gene of the set and the evidence points to bacteria. Bacterial records of
lineages that GTDB lacks, or just below 90 %, stay in MIDORI2 and in the E of every ASV. Such
a record can also make check 5.4 reject a genuine bacterial gene, as happened to one
*Polaromonas* fragment in development and probably to 8 genes in R232 (step 6).

**Genes left out for review.** Whole genes marked "review" by check 5.4 are left out of
version 1, although some may be genuine. Their lineages then have fewer genes on the bacterial
side of the comparison.

**Eukaryotes moved to the bacterial table.** In development, the comparison moved 0.04 % of
animal sequences whose family was missing from MIDORI2, and 0.14 % of those whose order was
missing. In groups that MIDORI2 covers poorly, a few genuine eukaryotic ASVs can therefore end
up in the table of bacterial COX1 sequences; check that table before concluding that an
expected group is absent.

**Labels of families B and C and of nitric oxide reductases.** The family labels rest on one
alignment to the Pfam model and on residues read from it. Divergent sequences can be
misaligned, and nitric oxide reductases are not separated from families B and C. A dedicated
classifier for the oxidase families was published (Sousa et al. 2011), but it is no longer
available at the address given in the paper.

**Genetic codes.** The table GTDB used is inferred from the genes, not read from the files.
Lineages with genetic code 25 (Gracilibacteria, and Absconditabacteria, formerly SR1; both read
TGA as glycine) are predicted by GTDB under table 4 or 11. Under table 11 some of their genes
are truncated or missing; under table 4 their TGA codons are translated as tryptophan instead
of glycine.

**The MIDORI2 list is a lower bound.** Short records with few tryptophan codons, chimeric
records (part animal, part bacterial), and records whose gene has no relatives in other genomes
and whose codons say nothing are hard to judge; they are listed, not removed. MIDORI2 also
holds sequences that are mitochondrial but not COI, such as the intron reading frame of step 6;
they are left in place.

**Three primer pairs tested.** The recognition rates were measured for BF3/BR2, Leray-XT and
fwhF2/fwhR2n. The set itself does not depend on primers, but the rates for other primer pairs
are not known.

**Shares count taxa, not reads.** All test results are shares of genomes or sequences, not of
reads. In a real library a few abundant ASVs, often from *Wolbachia*, can hold most of the
bacterial reads, so the share of bacterial reads that is recognised can be higher or lower than
the shares given here.

**What the set does not predict.** The set tells whether an amplicon is bacterial, not how
much bacterial DNA a sample will yield. That depends on primer mismatches, annealing
temperature and how much target DNA the sample holds. Off-target amplification is usually
much larger in eDNA, where animal DNA is scarce, than in bulk samples of animal tissue
(Elbrecht et al. 2019; Collins et al. 2019).

Possible improvements:

- **Tree placement.** Placing ASVs that are far from both references on a reference tree, as
  DARN does, instead of a single identity comparison.
- **Protein-level comparison** for very divergent ASVs, where nucleotide identity loses
  resolution.
- **Which bacteria each primer pair amplifies.** A table of the bacterial lineages that common
  primer pairs can amplify, distributed with the release, for study design.
- **More genomes per species,** for variation within species.
- **Genomes outside GTDB,** with genes predicted by Pyrodigal (Larralde 2022). Pyrodigal is
  tested to give the same predictions as the current Prodigal code (v2.6.3+31b300a, which
  differs from the 2.6.3 release in a few edge cases) and runs faster.
- **A curated classifier for the oxidase families,** to separate family B, family C and the
  nitric oxide reductases.
- **Sharing the MIDORI2 list** with the MIDORI2 maintainers.

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
