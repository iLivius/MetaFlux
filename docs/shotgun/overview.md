# Shotgun mode overview

Shotgun mode profiles whole-metagenome Illumina libraries. There is no PCR and no
primer: every read pair that comes off the sequencer is matched against a k-mer
database, so what limits the result is sequencing depth and the reference database
rather than the resolving power of one marker gene.

The pipeline is [Kraken2](https://github.com/DerrickWood/kraken2) for read
classification and [Bracken](https://github.com/jenniferlu717/Bracken) for abundance
re-estimation, with decontamination and quality trimming in front of them and a
combined OTU table at the end. It is selected with one config key:

```yaml
mode: shotgun          # amplicon | shotgun
```

The mode can also be set on the command line without editing the file:
`--config mode=shotgun`. Which mode suits which data is covered in
[Choosing a mode](../getting-started/choosing-a-mode.md).

!!! warning "`mode` must match the data"

    Both modes accept `{sample}_R1` / `{sample}_R2` filenames, so pointing a shotgun
    run at amplicon data (or the reverse) does not always produce an error — it
    produces the wrong analysis. Check `mode` before every run.

Shotgun mode needs one thing that is not fetched automatically: a Kraken2 database.
See [Reference databases](databases.md).

## The chain

```text
raw R1/R2 ──▶ BBDuk PhiX* ──▶ BBMap host* ──▶ fastp ──▶ Kraken2 ──▶ Bracken ──▶ kraken-biom
                _dephix         _dehost        _trim      report      report     otu_table
                                                            │
                                                            └──▶ KrakenTools per-taxon reads*

  * optional, toggled in config
```

Each stage writes files tagged with the step it came off (`_dephix`, `_dehost`,
`_trim`). When an optional step is switched off, the next step simply reads from
whatever ran before it — the chain closes up on its own, and no placeholder files are
produced.

The rules carry descending Snakemake priorities (`download_sra` 10 down to
`kraken_biom` 1), so when many jobs are ready at once the scheduler works through the
early stages first rather than starting classification while reads are still being
fetched.

## 1. Input — local FASTQs or a BioProject

Reads come from one of two places.

**Local files** in `input.fastq_dir`. Shotgun mode accepts either naming convention
and detects which one is in use:

| Layout | Example | Where it usually comes from |
|---|---|---|
| `{sample}_R1` / `{sample}_R2` | `soil01_R1.fastq.gz` | Illumina/bcl2fastq output, MetaFlux's default |
| `{sample}_1` / `{sample}_2` | `SRR11605259_1.fastq.gz` | raw `fasterq-dump` output |

Reads must carry the `fastq.gz` extension, or plain `fastq` if uncompressed. `fq` and
`fq.gz` are accepted by the startup check but are never resolved afterwards — the
rules look only for `.fastq.gz` and `.fastq` — and in shotgun mode the missing
`.fastq.gz` path matches the `download_sra` rule under `_R1`/`_R2` naming, so the run
tries to fetch each local sample from SRA rather than reporting that the extension is
unsupported; with `_1`/`_2` naming it stops with a missing-input error instead. The check itself sees
only the R1 files: a mixture of extensions among them stops the run at startup, but
the R2 files are never inspected, so uniformity across the directory is not enforced.

**Direct from NCBI SRA.** Setting `input.bioproject` to a BioProject accession makes
MetaFlux resolve that project's runs itself:

```yaml
input:
  fastq_dir: /path/to/fastq_dir      # must be writable — downloads land here
  bioproject: PRJNA603575
  accession_list: null               # e.g. [SRR11605259, SRR11605260]
```

The accession is resolved through NCBI e-utils (`esearch` then `efetch` for the
runinfo table) and the run list is cached in `fastq_dir` as
`{BIOPROJECT}_runinfo.csv`, so a re-run does not query NCBI again. Only runs whose
`LibraryLayout` is `PAIRED` are kept — single-end runs do not fit the rest of the
pipeline. `accession_list` narrows that set further, which is how a project that
mixes WGS with amplicon or RNA-Seq runs is reduced to the metagenomes.

Each run is then fetched with `prefetch` and unpacked with `fasterq-dump`
(`--split-files --skip-technical`), and the two mates are compressed to
`{sample}_R1.fastq.gz` / `{sample}_R2.fastq.gz` in `fastq_dir` — so downstream rules
see the same filenames whether the data was downloaded or already local.

Whatever the source, the sample list is written to `samples.tsv` at the top of the
output directory, with one row per sample: `sample`, `source` (`local` or
`bioproject`), `bioproject`, `r1`, `r2`. It exists so downstream R or Python scripts
can read the manifest instead of re-parsing the config or globbing the FASTQ
directory again.

## 2. PhiX removal (optional)

*Rule `decontam_phix` — off with `shotgun.decontamination.remove_phix: false`.*

PhiX is the control genome spiked into most Illumina runs; those reads are real
sequence but they are not part of the sample.
[BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) matches the reads
against the PhiX reference bundled inside the BBTools package, so nothing is
downloaded or indexed.

Input: the raw FASTQs. Output: `{sample}_dephix_R1/_R2.fastq.gz` (temporary — see
below) plus `{sample}_dephix_stats.txt`, which feeds both MultiQC and the read
tracking table. Consumed by the host step, or by fastp when host removal is off.

!!! warning "BBDuk needs R1 and R2 to be properly paired"

    BBDuk reads the two files in lockstep and requires them to hold the same number of
    records in the same order. Read length does not matter — ragged reads, and mates
    of a pair with different lengths, are processed normally. A crash with
    `AssertionError: List size mismatch: N vs M` means the pairing is broken, usually
    because reads were filtered before deposition and orphaned mates were left in one
    file — common for SRA data. `repair.sh`, from the same BBTools package, puts the
    files back in step. Details are in [Decontamination](decontamination.md).

## 3. Host and contaminant removal (optional)

*Rules `fetch_host_refs`, `build_host_index`, `decontam_host` — off when
`shotgun.decontamination.host_genomes` is an empty list.*

Host DNA dominates many metagenomes (gut, skin, plant tissue) and, left in, it costs
classification time and can be mis-assigned to microbial taxa.
[BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) aligns the reads
against one or more user-supplied genomes and keeps only the reads that do **not**
map.

Input: the PhiX-filtered reads, or the raw reads when PhiX removal is off. Output:
`{sample}_dehost_R1/_R2.fastq.gz` plus `{sample}_dehost_stats.txt` with the surviving
pair count. Consumed by fastp.

Which reference to use — and why a masked human genome matters — is covered in
[Decontamination](decontamination.md).

## 4. Adapter and quality trimming

*Rule `trim_adapters`, always runs.*

[fastp](https://github.com/OpenGene/fastp) removes adapter read-through, trims poor
quality from both ends with a sliding window (`--cut_front --cut_right`), and drops
reads that end up shorter than `shotgun.fastp.min_read_length` (default `100`;
around `70` suits HiSeq 2×100 data). Adapters are detected from the data itself
(`--detect_adapter_for_pe`), so no adapter file is needed.

Input: whichever decontaminated stage ran last, or the raw reads if both are off.
Output: `{sample}_trim_R1/_R2.fastq.gz` (temporary) and the `{sample}_fastp.html` /
`{sample}_fastp.json` reports. The JSON is used three times over: by MultiQC, by the
read tracking table, and by Bracken, which reads the post-trim mean read length out
of it.

## 5. Kraken2 classification

*Rule `kraken2`.*

Kraken2 walks each read pair k-mer by k-mer, looking up the minimizer of each k-mer in
the database — every k-mer sharing a minimizer inherits that minimizer's stored lowest
common ancestor, which is why one database entry can account for a whole run of
apparent k-mer hits. The taxa hit form a small tree, each node weighted by how many
hits it received, and the read is labelled with the leaf of the highest-scoring
root-to-leaf path (the LCA of the leaves if two paths tie). A handful of stray hits to
a distant taxon therefore does not drag the whole read up the taxonomy.

```yaml
shotgun:
  kraken:
    confidence: 0.15      # higher = stricter
    hit_groups: 3         # minimum distinct hit groups per read
    memory_mapping: true
    threads: null         # null = use resources.threads.kraken2
```

`confidence` is the fraction of a read's k-mers that must fall inside the clade of the
taxon the read is assigned to. A read that falls short is not simply rejected: Kraken2
moves the label up the taxonomy until the score clears the threshold, and calls the
read unclassified only when not even root qualifies. Raising it therefore mostly
pushes assignments to coarser ranks — which leaves more work for Bracken in step 6 —
and only secondarily increases the unclassified fraction. A sweep from 0.00 to 1.00 on
simulated samples (Wright et al. 2023) found the best setting depends on what is being
optimised: 0.60 gave the highest mean F1, while 0.15 gave the lowest L1 distance and the
best recall and precision over the reads actually classified. MetaFlux defaults to 0.15,
which favours accurate relative abundances; 0.15 is also the threshold used for CAMI2
benchmarking by Nyström-Persson et al. (2025).

`hit_groups` counts hit *groups*, not k-mer hits. A hit group is a run of overlapping
k-mers that share one minimizer found in the database — one lookup, at one place in
the read — so a read can show dozens of k-mer hits and still amount to a single group.
`hit_groups: 3` requires three such separate hits before a classification stands
(Kraken2's own default is 2), leaving a read unclassified when all its evidence sits
at one locus.

Input: the fastp-trimmed reads and the Kraken2 database. Output, in
`02.classification/`:

| File | Contents |
|---|---|
| `{sample}_report.txt` | the per-taxon summary tree, written with `--report-minimizer-data` |
| `{sample}_output.txt` | one line per read pair with its assignment |
| `{sample}_R1.fastq.gz` / `_R2.fastq.gz` | the reads Kraken2 could classify |

The classified-read FASTQs are kept, not deleted, because per-taxon read extraction
draws from them. The report is consumed by Bracken and by MultiQC; the per-read
output is consumed by read extraction.

!!! note "The report format is not the standard one"

    `--report-minimizer-data` adds two columns to the usual Kraken2 report, which
    shifts the rank and taxid columns. KrakenTools reads the report as it stands,
    but anything else parsing these reports needs to account for the shift.

Memory behaviour is the main practical concern here and is covered in
[Reference databases](databases.md#memory).

## 6. Bracken abundance re-estimation

*Rule `bracken`.*

Kraken2 leaves many reads assigned above species level, at whatever ancestor the
k-mers agreed on. Bracken pushes those reads back down using the database's k-mer
distributions, giving abundances at a chosen rank.

```yaml
shotgun:
  bracken:
    tax_lev: S            # D | P | C | O | F | G | S
    threshold: 10         # minimum read count per taxon
```

Bracken's estimates depend on read length, and the database ships pre-computed
distributions for several lengths (50, 75, 100, 150, 200, 250 and 300 bp). MetaFlux
reads the mean post-trim read length from the sample's fastp JSON and picks the
closest available distribution automatically — using a 150 bp distribution on 100 bp
reads would bias the redistribution.

Input: the Kraken2 report plus the fastp JSON. Output, in `03.abundance/`:
`{sample}_report.txt` (Bracken's re-estimated report) and `{sample}_output.txt`
(the per-taxon abundance table). Consumed by kraken-biom, by MultiQC, and by the
optional read extraction step, which looks up taxids in the Bracken report.

!!! tip "Drop to genus on a shallow database"

    Species-level re-estimation is only as good as the database's species
    representation. With a small or size-capped index, `tax_lev: G` gives a more
    honest table.

!!! warning "The numbers in the OTU table are estimates, and they do not sum to the classified reads"

    Two things about Bracken's output are easy to misread, and both matter when the
    numbers reach a methods section or a normalisation step.

    **They are estimated, not observed.** A value in `otu_table.tsv` is not a count of
    reads that Kraken2 assigned to that taxon. It is Bracken's estimate of how many
    reads *originated* from it, after redistributing reads that Kraken2 had left at
    coarser ranks. Reporting Bracken's estimates as observed assignments overstates what
    was measured.

    **They do not account for every classified read.** Bracken reports only at the
    target rank. Reads classified above it that cannot be pushed down — because the
    evidence does not separate the descendants — are absent from the output altogether.
    So the column total of `otu_table.tsv` is smaller than the `classified` column in
    `stats/read_tracking.txt`, and the two are not meant to reconcile. When a relative
    abundance is computed, be explicit about which denominator was used: the OTU-table
    total, or the classified reads.

## 7. OTU table

*Rules `kraken_biom` and `finalize_otu_table`.*

[kraken-biom](https://github.com/smdabdoub/kraken-biom) merges every sample's Bracken
report into one taxon × sample table in BIOM (HDF5) format. That raw table is then
tidied in a second step, which:

- strips the `_report` suffix that kraken-biom copies from the filenames onto the
  sample IDs;
- rebuilds the species field into a `Genus species` binomial (kraken-biom emits a
  bare epithet, e.g. `s__denitrificans`), so shotgun and amplicon tables carry the
  same kind of lineage string — ranks joined by `;` with no space, each prefixed
  `k__ p__ c__ o__ f__ g__ s__`;
- applies the optional taxon filter.

```yaml
shotgun:
  taxonomy_filter:
    enabled: false
    keep:    []           # e.g. [k__Bacteria, k__Archaea]
    discard: []           # e.g. [g__Homo]
```

The filter uses the same rank-aware token matching as the amplicon side — each token
must match a whole segment of the lineage, never a loose substring — and it runs on
the final table only, after Bracken. It therefore changes which OTUs appear in the
table, never how a read was classified. It is off by default because metagenomic
sample types vary too much for a sensible universal setting. Common uses are keeping
prokaryotes (`keep: [k__Bacteria, k__Archaea]`) and dropping residual host
(`discard: [g__Homo]`). The machinery is described in full under
[Keeping and discarding taxa](../amplicon/taxon-filter.md).

Output: `03.abundance/otu_table.biom` and `03.abundance/otu_table.tsv`. Both are
rendered from the same filtered table, so they always agree. These are the main
deliverable of a shotgun run.

## 8. Per-taxon read extraction (optional)

*Rules `extract_taxon_reads` and `compress_extracted` — only run when
`shotgun.extract_taxa` is non-empty.*

[KrakenTools](https://github.com/jenniferlu717/KrakenTools) pulls the reads assigned
to a named taxon, and everything below it, into per-sample FASTQ files — the starting
point for a targeted assembly or a closer look at one organism. See
[Per-taxon read extraction](read-extraction.md).

## What lands on disk

```text
out_dir/
├── 01.preprocessing/     # dephix / dehost stats, fastp JSON + HTML, host index
├── 02.classification/    # Kraken2 reports, per-read output, classified reads
├── 03.abundance/         # Bracken reports, otu_table.tsv / otu_table.biom
├── 04.extracted_reads/   # KrakenTools per-taxon reads (only if extract_taxa set)
├── stats/                # read_tracking.txt
├── multiqc/              # multiqc_report.html
├── samples.tsv           # sample → source → FASTQ manifest
└── logs/                 # per-rule logs
```

The cleaned reads themselves are marked temporary: Snakemake deletes the `_dephix`,
`_dehost` and `_trim` FASTQs as soon as the rule that needed them has finished, so
`01.preprocessing/` keeps QC and statistics rather than several copies of the
library. What stays there is:

| File | Present when |
|---|---|
| `{sample}_fastp.html` / `_fastp.json` | always |
| `{sample}_dephix_stats.txt` | `remove_phix: true` |
| `{sample}_dehost_stats.txt` | `host_genomes` is set |
| `refs/bbmap_index/` (BBMap host index) | `host_genomes` is set |

The concatenated host reference, `refs/host_combined.fna.gz`, is temporary as well:
it is deleted once the index has been built, which is why a URL in `host_genomes`
is fetched again on every fresh run.

`stats/read_tracking.txt` is the one file worth reading first after a run. One row
per sample, counts in read **pairs** throughout, columns in pipeline order:

```text
sample → raw → [nophix] → [no_host] → trimmed → classified → unclassified
```

The bracketed columns appear only when the corresponding step ran. A large drop at
`no_host` means the sample was host-dominated — but only when PhiX removal ran too.
With `remove_phix: false` the `raw` column is taken from the fastp report, that is,
from the reads *entering* fastp, which is the output of host removal; `raw` and
`no_host` then count exactly the same reads and no host drop can appear however
host-dominated the sample is. In that configuration read the mapping rate out of
`logs/decontam_host/{sample}.log` instead.

A large `unclassified` fraction usually points at the database rather than the data —
an environment poorly covered by the chosen index, not a failed run.

`multiqc/multiqc_report.html` collects the fastp, Kraken2, Bracken and BBTools
statistics into one page.

## Where to go next

- [Reference databases](databases.md) — picking a Kraken2 index and keeping its
  memory use under control.
- [Decontamination](decontamination.md) — PhiX and host removal in detail.
- [Per-taxon read extraction](read-extraction.md) — pulling reads for one taxon.
