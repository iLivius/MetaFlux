# sylph, the alternative profiler

[sylph](https://github.com/bluenote-1577/sylph) is a second way to turn the cleaned
shotgun reads into a taxon table. It runs instead of Kraken2 + Bracken, or next to them,
when `shotgun.profiler` asks for it. The cleaning in front (PhiX removal, host removal,
fastp) is the same whichever profiler runs.

!!! note "Kraken2 + Bracken stays the default"

    sylph is available from v2.5.0 and is used only when `shotgun.profiler` is set to
    `sylph` or `both`. With the default, `kraken2`, nothing on this page applies.

## What sylph does differently

- **It compares sketches, not reads.** sylph reduces each sample to a sketch, a fixed
  subset of about one in 200 of its 31-bp k-mers, and compares it with sketches of whole
  reference genomes.
- **It reports genomes, at about species level.** A genome is reported when the sample
  holds something at 95 % average nucleotide identity (ANI) or more to it, roughly the
  species boundary. Each row of the table is one reference genome.
- **Abundance comes from coverage.** From how many of a genome's k-mers turn up, and how
  often, sylph estimates the genome's coverage. Coverage times genome size gives the
  share of the sample's sequence that genome explains.
- **No read gets a label.** sylph writes no per-read output, so
  [read extraction](read-extraction.md), which picks reads by the taxon Kraken2 gave each
  one, has nothing to work from. The [k-mer evidence gate](kmer-evidence.md) needs
  Kraken2 too: it filters the Kraken2 report before Bracken, and under `sylph` neither
  runs.
- **Same reads as Kraken2.** sylph sketches exactly the trimmed read pairs Kraken2 would
  classify. With `profiler: both` the two tables come from identical input; they are
  written side by side and never merged.

## Switching it on

```yaml
shotgun:
  profiler: sylph              # kraken2 (default) | sylph | both

references:
  sylph_databases:
    - path: /db/sylph/gtdb-r232-c200-dbv1.syldb
      taxonomy: GTDB_r232
    - path: /db/sylph/fungi-refseq-2025-10-11-c200-dbv1.syldb
      taxonomy: FungiRefSeq-latest
  sylph_taxonomy_dir: /db/sylph_taxonomy
```

With `sylph`, the Kraken2 database is not needed and `references.kraken_db` can keep its
placeholder. With `both`, Kraken2 + Bracken and sylph run on the same reads in one run,
and fastp still runs once per sample.

All databases are profiled in one sylph run, so within a sample an organism present in
two of them is reported once, by its closer genome.

Everything in the block above is checked before the first job starts: the file names, a
taxonomy for every database and the lineage table behind each. When a database file name
shows its release (`gtdb-r232-c200-dbv1.syldb`, or a copy renamed `gtdb_r232_c200.syldb`),
MetaFlux also checks that the taxonomy comes from the same release. That covers a
built-in name and a lineage table given by its path under a built-in table's file name
(`gtdb_r226_metadata.tsv.gz`, say); a table of your own under another name cannot be
matched to a release at this point. After profiling, MetaFlux checks the GTDB release once
more, from the genome paths in sylph's output, whatever the files are called. A mistake
stops the run with a message that says what to change.

## Databases and lineage tables, set up once

MetaFlux does not download sylph databases, just as it does not download a Kraken2
index. The pre-built ones are listed on sylph's
[Pre-built databases](https://sylph-docs.github.io/pre%E2%80%90built-databases/) page.

- **Sketch density.** MetaFlux sketches the reads at c 200, and a database must be built
  at the same density or a coarser one. c 200 builds find low-coverage genomes best; c 1000
  builds need about a fifth of the memory. A c 100 build would make sylph skip every
  sample, so MetaFlux stops on it.
- **File types.** Both `.syldb` and the two-stage `.syl2db` files work.
- **One database per family.** Two GTDB releases, or the c 200 and c 1000 builds of one
  release, share genomes and names, and listing both stops the run. GTDB and GlobDB
  together are allowed, with a warning, but best avoided. GlobDB already holds every
  genome of the GTDB database, with the same lineage, under its accession without the
  version suffix (`GCF_003096595` for GTDB's `GCF_003096595.1`). sylph reports such a
  genome from one database or the other, and which one can change from sample to
  sample, so one genome can end up on two rows, each with zeros where the other copy was
  reported. Genome-level counts, prevalence and diversity are then wrong. The pair also
  roughly doubles the memory without adding species; use GlobDB alone instead.
- **No viral databases.** IMG/VR and UHGV use viral ranks that do not fit the table
  and are refused.

Each database is paired with the sylph-tax taxonomy that names its genomes:

| Database file | `taxonomy:` |
|---|---|
| `gtdb-r232-c200-dbv1.syldb` or `gtdb-r232-c1000-dbv1.syldb` | `GTDB_r232` |
| `gtdb-r226-c200-dbv1.syldb` or `gtdb-r226-c1000-dbv1.syldb` | `GTDB_r226` |
| `fungi-refseq-2025-10-11-c200-dbv1.syldb` | `FungiRefSeq-latest` |
| a database you built | the path of your own `accession<TAB>lineage` table (`.tsv` or `.tsv.gz`) |

The built-in taxonomies read their lineage tables from `references.sylph_taxonomy_dir`.
Fill that folder once, with the sylph-tax version MetaFlux pins in
`workflow/envs/sylph.yaml` (sylph 1.0.0, sylph-tax 1.9.2):

```bash
conda env create -n metaflux-sylph -f workflow/envs/sylph.yaml
conda run -n metaflux-sylph sylph-tax --no-config --taxonomy-dir /db/sylph_taxonomy download
```

This fetches every table sylph-tax knows, about 80 MB, from Zenodo and, for GlobDB r232,
from the GlobDB server. One failed download stops the whole command. The tables can also
be fetched one at a time, for example the two used above:

```bash
wget -P /db/sylph_taxonomy https://zenodo.org/records/19646381/files/gtdb_r232_metadata.tsv.gz
wget -P /db/sylph_taxonomy https://zenodo.org/records/17330476/files/fungi_refseq_2025-10-11_metadata.tsv.gz
```

MetaFlux always calls sylph-tax with `--no-config --taxonomy-dir`, so nothing is read
from or written to `~/.config/sylph-tax`.

## What you get

| File | What it holds |
|---|---|
| `03.abundance/sylph_table.tsv` | the main table: estimated read pairs per genome and sample, and a last row `unassigned` |
| `03.abundance/sylph_table.biom` | the same table in BIOM (HDF5) format |
| `03.abundance/sylph_table.genome_copies.tsv` | each genome's share of genome copies (%) among the detected genomes; same rows, no `unassigned` |
| `03.abundance/sylph_genomes.tsv` | one row per genome: species name, whether it is a placeholder, the FASTA header of the reference genome, the placeholder-genus flag, best coverage and ANI |
| `03.abundance/sylph_profile.tsv` | sylph's own output, unchanged: one row per sample and genome |
| `stats/sylph/sylph_summary.tsv` | one row per sample: read pairs after QC, how much of them the detected genomes explain, genomes found, flagged genomes |

`stats/read_tracking.txt` is written as in every shotgun run; under `profiler: sylph` it
has no Kraken2 columns and ends with `trimmed`. The MultiQC report gets a table section
for sylph. `02c.sylph/` keeps the sketches, so a later run in the same folder with other
databases does not sketch again.

In `sylph_table.tsv`, line 1 is a comment that says what the numbers are and which
databases they come from, and line 2 is the header, with `#OTU ID`, the samples and
`taxonomy`. Each row ID is a genome accession, and the lineage has seven slots, the GTDB
ranks from domain to species:

```text
#OTU ID	sampleA	sampleB	taxonomy
GCF_000146045.2_R64	1520	0	d__Eukaryota;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces cerevisiae
unassigned	0	98211	d__Unassigned;p__;c__;o__;f__;g__;s__
```

The first slot always holds the domain and is always written `d__`: `d__Bacteria`,
`d__Archaea` or `d__Eukaryota`, and `d__Unassigned` on the `unassigned` row. GTDB has no
kingdom rank, so there is no `k__` slot. The FungiRefSeq-latest lineage table writes its
domain as `k__Eukaryota`, and the TaraEukaryoticSMAG table spells it `Eukarya` in some
rows; MetaFlux writes both as `d__Eukaryota`, so every row uses the same seven prefixes.
A lineage whose first rank is not a domain stops the run with a message that names the
lineage table. The lineage of the Kraken2 table is described on the
[output files](../reference/output.md#taxonomy-strings-in-the-final-tables) page.

## Reading the numbers

- **Estimated read pairs, not counted ones.** A value in `sylph_table.tsv` is sylph's
  `Sequence_abundance`, the share of the sample's reads that genome explains, estimated
  from its coverage, times the sample's read pairs after QC, rounded to a whole pair. No
  read was assigned to the genome. Each column adds up to the `trimmed` count of
  `stats/read_tracking.txt`, except by a few pairs per million in a sample sylph explains
  to 100 % (`explained_capped: yes`). There the rounded genome rows can add up to
  slightly more than `trimmed` (`unassigned` is then 0, and `logs/sylph_tables.log` says
  so), or leave a few pairs over in `unassigned`.
- **`unassigned` is a model estimate as well.** MetaFlux runs `sylph profile -u`, which
  estimates how much of each sample the detected genomes do not explain; that remainder
  is the `unassigned` row. It depends on the database, since a fuller database explains
  more. `unassigned` cannot be compared with Kraken2's unclassified read pairs. Leave the
  row out before computing relative abundances of the detected genomes.
- **Read accuracy.** The estimates also depend on sylph's estimate of read accuracy.
  When the sample's genomes are covered too thinly for that estimate (common in shallow
  samples, and in very diverse ones such as soil even when deeply sequenced), sylph
  assumes 99.5 % instead; `read_accuracy` in the summary then reads
  `fixed 99.5% (low depth)`, the depth being how deeply the genomes are covered. That
  assumption moves every genome row of the sample by one common factor, and the
  `unassigned` row with them, which can change the explained share by tens of percentage
  points; which genomes are found, and their genome copies, do not change. Compare
  explained shares, `unassigned` and estimated pairs only between samples with the same
  `read_accuracy`. Otherwise leave out `unassigned` and compare shares among the detected
  genomes, or use `sylph_table.genome_copies.tsv`.
- **Capped samples.** `explained_capped: yes` usually means that sylph's coverage
  estimates added up to at least all of the sample's sequence (sylph does not say by how
  much; in a mock community it can be several percent). sylph then shares 100 % among
  the detected genomes, as it does without `-u`, so every genome's estimated read pairs
  in that sample are scaled down by that unknown amount, and they change when a database
  is added or removed, even for a genome whose coverage stays the same. Nothing is left
  for `unassigned` to estimate. Compare a capped sample with others by shares among the
  detected genomes or by genome copies, not by estimated pairs.
- **Genome copies are a different share.** `sylph_table.genome_copies.tsv` gives each
  genome's share of genome copies, roughly of cells, among the detected genomes only. A
  genome that is large but present in few copies has a high read share and a low copy
  share.
- **Names are kept as GTDB writes them.** That includes placeholder species such as
  `s__JABGUH01 sp018646675`, genera and species split with a suffix (`g__Blautia_A`;
  `s__Escherichia coli_E` is a real, named species), and GTDB's own genus names
  (`s__Methanocatella smithii` for the species long known as *Methanobrevibacter
  smithii*). Comparisons at genus level are less affected than at species level. The
  `name_status` column of `sylph_genomes.tsv` says for every genome whether its name is
  a placeholder.

### The placeholder-genus flag

GTDB sometimes puts a genome in a genus known only by a code (a placeholder genus such as
`g__G047199095`) while the genome's own FASTA header, written by whoever submitted it,
names a Latin species of another genus. In GTDB r232 this is the case for some
*Escherichia coli* assemblies that are new in that release: reads of an *E. coli* strain
can land on such a genome, so
that the strain appears under a placeholder name and not as *Escherichia coli*. MetaFlux
renames and merges nothing.
It flags the genome in `sylph_genomes.tsv` (`GTDB placeholder genus; genome labelled
Escherichia coli`) and counts the flagged genomes of each sample in the summary.

## Reading the tables

The first line of each table is a comment, so skip it:

```r
tab <- read.delim("03.abundance/sylph_table.tsv", skip = 1, row.names = 1,
                  check.names = FALSE, comment.char = "")
```

```python
import pandas as pd
tab = pd.read_csv("03.abundance/sylph_table.tsv", sep="\t", skiprows=1, index_col=0)
```

Do not pass `comment="#"` to pandas: the header line itself starts with `#OTU ID`.

phyloseq reads the BIOM file. Its `parse_taxonomy_greengenes` knows the prefixes `k__` to
`s__` but not `d__`, so the first column comes out named `Rank1`; rename it:

```r
library(phyloseq)
ps <- import_biom("03.abundance/sylph_table.biom",
                  parseFunction = parse_taxonomy_greengenes)
colnames(tax_table(ps))[1] <- "Domain"
```

## Settings that do not apply

Under `profiler: sylph`, `shotgun.extract_taxa` (it needs Kraken2's per-read output) and
`shotgun.kmer_evidence` (it filters the Kraken2 report before Bracken) are switched off
with a warning; under `both` they still act on the Kraken2 side. `shotgun.taxonomy_filter`
and the `shotgun.kraken` and `shotgun.bracken` settings act on the Kraken2 chain only.
The sylph tables are written unfiltered, so each column still adds up to the read pairs
after QC (with the one exception under [Reading the numbers](#reading-the-numbers)).

## Limits

- **Sample names must not look like numbers.** sylph-tax reads sample names with pandas,
  which turns `001` into `1`, merges `01` with `1` and drops `NA`. MetaFlux stops at
  start-up on such names; rename the files so that each name contains a letter (`S001`).
  Names such as `12_3` are fine.
- **A sample with no read pairs left after QC stops the run**, because sylph cannot
  sketch an empty file. Remove it from `fastq_dir`, or profile it with Kraken2.
- **Adding sylph to a finished Kraken2 run** needs a new `out_dir`, with `profiler: sylph`.
  The trimmed reads of the finished run were deleted after use, so `both` in the same
  folder would make them again and run Kraken2 again. MetaFlux stops and says so.
- **Lineages must start with a domain.** Among the built-in taxonomies, only GlobDB r226
  has genomes without one (`d__Unclassified`); a sample in which sylph detects one of
  them stops at `sylph_tables`. GlobDB r232 places every genome in a domain.
- **Memory.** sylph holds every database in memory, in a single `sylph_profile` job for
  all samples. MetaFlux asks for 1.2 times the size of the database files plus 4 GB; set
  `resources.mem_mb.sylph_profile` to a number to override it, for example on a cluster,
  but keep it above about 1.1 times the size of the files, plus some headroom. GTDB r232
  at c 200 with the fungal database is about 27 GB of files. With a few samples, most of
  the job's time goes into reading the databases from disk. Do not size the job from
  `benchmarks/sylph_profile.tsv`: sylph's memory climbs while it reads the databases and
  peaks only in the last seconds of the job, and after its first 30 checks (about the
  first 15 to 20 seconds) Snakemake records memory only every 30 seconds, so the recorded
  `max_rss` usually misses the peak and can be far below it.
- **Memory under `both`.** The `sylph_profile` job can start while Kraken2 jobs are
  running. On a local run, Snakemake keeps their memory apart only when it is given a
  cap, with `--resources mem_mb=...`. Give a cap at least as large as the `sylph_profile`
  figure printed at start-up plus what the Kraken2 jobs declare. A smaller cap does not
  stop the run: Snakemake lowers the job's request to the cap and runs it anyway.
- **Reproducibility.** Sketches are made with one thread
  (`resources.threads.sylph_sketch: 1`), because with more threads they differ slightly
  from run to run; with one thread they are the same every time, also with PhiX or host
  removal on. The thread count of `sylph_profile` changes no number. Genomes with
  exactly the same abundance, common in mock communities, can come out in a different
  order from run to run in `sylph_profile.tsv` and in sylph-tax's files, so sort those
  before comparing two runs; the MetaFlux tables sort their own rows. Each sketch records
  the path of the reads it was made from, so sketches of the same reads in two `out_dir`s
  are not byte-identical.
