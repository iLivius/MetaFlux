# Changelog

What changed in each MetaFlux release, newest first. Each entry says what changes
for a user, what runs again in an output folder made with the previous release, and
what to edit in an existing run config. The documentation site has the details; the
release notes on GitHub repeat the longer entries. Every release is archived on
Zenodo; the concept DOI in the README always resolves to the newest one.

## 2.5.0 — unreleased

### Taxonomy strings: `d__` is the domain, `k__` the kingdom

SILVA 144 and NCBI now carry both a domain and a kingdom for prokaryotes. MetaFlux
used to write the domain under `k__` and kept SILVA 144's kingdom in a side file. From
this release on, `d__` is always the domain slot and `k__` always the kingdom slot, and
a table carries a kingdom only where its reference defines one.

| Output | v2.4.0 | v2.5.0 |
|---|---|---|
| 16S, SILVA 144 | `k__Bacteria;p__Pseudomonadota;…` (kingdom in a side file) | `d__Bacteria;k__Pseudomonadati;p__Pseudomonadota;…` |
| 16S, SILVA 138.2; rpoB | `k__Bacteria;p__…` | `d__Bacteria;p__…` |
| ITS; gyrB | unchanged | unchanged |
| 18S, `sintax` | empty for an ASV with an uncertain supergroup | at least `d__Eukaryota` |
| Kraken2, NCBI database | `k__Bacteria;p__…;s__coli` | `d__Bacteria;k__Pseudomonadati;p__…;s__Escherichia coli` |
| Kraken2, a virus | `k__;p__Uroviricota;…;s__P37` | `d__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;…;s__Peduovirus P37` |
| Kraken2, GTDB database | `k__Archaea;p__…` | `d__Archaea;k__;p__…` |
| sylph (new) | — | `d__Bacteria;p__…` |

- **16S classifies on SILVA 144 by default**, on all of its ranks, kingdom included
  and classified with the same cutoff as every other rank; the kingdom is written
  without SILVA's `--other` suffix (`amplicon.taxonomy.silva_other` says what to do
  with the other `--other` names). SILVA 138.2 stays available with
  `references: {silva: {release: "138.2"}}`. `06.taxonomy/taxonomy_lineage_full.tsv`
  is no longer written.
- **18S under `sintax`: the domain is classified first.** VSEARCH reads the SINTAX
  letters as a fixed top-down order, and PR2's own file puts the domain under the
  second letter, so an ASV lost its domain whenever its supergroup was uncertain.
  MetaFlux now classifies against a copy with the two letters swapped. Only the
  domain label is gained; nothing below it changes. 18S `sintax` runs made with
  earlier versions have this problem; re-run them.
- **The shotgun OTU table is built from the Bracken reports**, without kraken-biom:
  eight lineage slots on every Kraken2 database (`d__` to `s__`), NCBI's kingdoms kept,
  a virus's realm in `d__`, species names complete, and a slot never filled from a
  neighbouring branch of the tree. The counts are the same Bracken estimates as
  before. In `otu_table.biom` the sample `Id` is the plain sample name.
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
the k-mer evidence gate need Kraken2's per-read output and are switched off under
`sylph`. The sylph tables are written unfiltered.

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
- **Bracken rank.** `shotgun.bracken.tax_lev: D` is refused like `K` (Bracken files
  the reads of taxa without the rank under the wrong taxon); `P`, `C`, `O` and `F`
  give a warning. `tax_lev` is read only when Kraken2 runs.
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

Two-digit folder numbers in both modes: amplicon `01.reads` … `07.phylogeny`, and the
read-count files in `stats/` take the same two digits (`01.raw_reads.counts` …);
shotgun `01.preprocessing`, `02.classification/kraken2/` and `02.classification/sylph/`
(one folder per profiler; the k-mer evidence files move in with Kraken2's),
`03.abundance`, `04.extracted_reads`. `stats/`, `logs/`, `benchmarks/` and `multiqc/`
keep their names.

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
- New keys with safe defaults, so an old config parses unchanged: `shotgun.profiler`
  (`kraken2`), `references.sylph_databases` and `sylph_taxonomy_dir`,
  `amplicon.taxonomy.silva_other` (`verbatim`), `filter.discard_min_support` (50),
  `filter.organelles` (`metaxa2`).

### Documentation and repository

New sylph page; the taxonomy, output, configuration, running and troubleshooting
pages rewritten for the domain and kingdom slots, the start-up checks and the new
folders; the docs deployment runs a strict build before publishing. The repository
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

## 2.3.0 — 2026-08-04

- Shotgun: a species-level k-mer evidence gate (`shotgun.kmer_evidence`, off by
  default) that removes taxa without enough distinct k-mers before Bracken, and a
  depth-scaled Bracken read threshold (`shotgun.bracken.threshold: auto`).
- The documentation site is deployed to GitHub Pages automatically.
- `sintax`: a malformed VSEARCH output line stops the run instead of being guessed at.

## 2.2.0 — 2026-07-30

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
