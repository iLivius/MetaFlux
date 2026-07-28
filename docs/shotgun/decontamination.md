# Decontamination

Two optional steps sit between the raw reads and quality trimming: PhiX removal with
BBDuk, and host removal with BBMap. Both exist for the same reason — reads that are
not part of the microbial community should not reach the classifier. They cost
classification time, they inflate the denominator of every relative abundance, and
host reads in particular can be assigned to real microbial taxa when a region of the
host genome resembles something in the database.

Both steps are configured in one block:

```yaml
shotgun:
  decontamination:
    remove_phix: true       # BBDuk vs BBTools' bundled PhiX; false = skip
    host_genomes: []        # list of paths/URLs to FASTA(.gz); [] = skip
    host_min_id: 0.95       # BBMap approximate minimum alignment identity
```

The chain is `raw → (PhiX) → (host) → fastp`. Each step reads whatever ran before
it, so switching one off is safe: with PhiX removal disabled, host removal reads the
raw FASTQs; with both disabled, fastp does. Nothing needs to be renamed and no empty
placeholder files are written.

Note that the decontaminated FASTQs are marked temporary — Snakemake deletes them
once the next step has consumed them. The statistics files they produce are kept, and
are what `stats/read_tracking.txt` reports on; of the two, only the BBDuk PhiX file is
in a format MultiQC parses.

## PhiX removal

*Rule `decontam_phix`.*

PhiX is the small bacteriophage genome routinely spiked into Illumina runs as a
sequencing control. Those reads are genuine sequence and will not be classified as
anything real, but removing them early keeps the read counts honest.

BBDuk matches reads by k-mer against `phix174_ill.ref.fa.gz`, which ships inside the
`bbmap` conda package. Nothing is downloaded and no index is built. The settings are
BBTools' own recommendation for PhiX:

```bash
bbduk.sh ref=phix in1=... in2=... out1=... out2=... stats=... k=31 hdist=1
```

`k=31` is the match length and `hdist=1` allows one mismatch within a k-mer, so a
read with a sequencing error still matches.

| | |
|---|---|
| **In** | raw `{sample}_R1/_R2` |
| **Out** | `01.preprocessing/{sample}_dephix_R1/_R2.fastq.gz` (temporary), `{sample}_dephix_stats.txt` |
| **Next** | host removal, or fastp when host removal is off |

The stats file carries BBDuk's `#Total` and `#Matched` counts. MultiQC picks it up
automatically, and the read tracking table derives both the `raw` and `nophix`
columns from it (halving BBDuk's per-read counts to get read pairs).

Leaving this on is the sensible default unless the lab is known to have skipped the
spike-in. In shotgun data PhiX is a small fraction anyway, so switching it off costs
few reads — but it also removes the only independent count of the reads *entering*
the pipeline, which is what the read-tracking table uses for its `raw` column.

## Host removal

*Rules `fetch_host_refs`, `build_host_index`, `decontam_host`.*

Host removal is alignment-based rather than k-mer-based. BBMap maps the reads against
the host reference and MetaFlux keeps the reads that do **not** map. Alignment is
more sensitive than k-mer matching for host reads that differ from the reference
assembly — an individual's variants, a related species standing in for the real host
— which is the approach the BBTools authors recommend (compare their `removehuman.sh`).

```yaml
shotgun:
  decontamination:
    host_genomes:
      - /path/to/human_virus_masked.fasta.gz
      - /path/to/plant_chloroplast.fasta.gz
    host_min_id: 0.95
```

`host_genomes` is a YAML **list**: each reference is its own `- ` item, and each may
be a local path or a URL (`http`, `https` or `ftp`). An empty list skips the whole
step. Multiple entries are concatenated into a single reference and indexed once, so
adding a second host costs one index build, not a second pass over the reads.

The three rules do this in order:

1. **`fetch_host_refs`** resolves every entry — downloading URLs, copying local
   files, gzipping any plain FASTA — and concatenates them into
   `01.preprocessing/refs/host_combined.fna.gz`. Gzip streams can be joined as they
   are, so a 3 GB genome is never re-compressed. A local path that does not exist
   stops the run with a clear error.
2. **`build_host_index`** builds the BBMap index once, under
   `01.preprocessing/refs/bbmap_index/`, and every sample reuses it.
3. **`decontam_host`** maps each sample against that index and writes the unmapped
   reads (BBMap's `outu1`/`outu2`).

| | |
|---|---|
| **In** | PhiX-filtered reads, or raw reads when PhiX removal is off |
| **Out** | `01.preprocessing/{sample}_dehost_R1/_R2.fastq.gz` (temporary), `{sample}_dehost_stats.txt` |
| **Next** | fastp |

BBMap does not print a total/matched table the way BBDuk does, so the rule counts the
surviving R1 records itself and writes a single `surviving_pairs` line. That is what
fills the `no_host` column of `stats/read_tracking.txt`.

### Use a masked human reference

For human host removal, the reference should be a **masked** one. Regions that are
conserved across life — ribosomal RNA genes above all — plus low-complexity repeats
and stretches with clear homology to bacteria, plants or fungi are masked out of the
reference. Without that masking, a real bacterial 16S read can align
convincingly to the human rDNA locus and be thrown away as host, quietly depleting
exactly the organisms the study is about.

The masked human reference on Zenodo
([Handley 2020](https://doi.org/10.5281/zenodo.4116107), 889 MB) is Brian Bushnell's
masked hg19 — ribosomal RNA, plant, animal and fungal homology and low-entropy
sequence already taken out, deleted from this file rather than N-filled — with a
further layer of viral masking added on top. The first set of masks is what stops
real microbial reads being discarded as host. The viral layer additionally keeps
human regions homologous to viruses out of the filter: the right behaviour when the
viral fraction is of interest, and a small loss of host depletion otherwise.

```bash
wget -O /path/to/db/human_virus_masked.fasta.gz \
  "https://zenodo.org/records/4116107/files/human_virus_masked.fasta.gz"
```

!!! tip "Download it once"

    A URL in `host_genomes` is re-fetched by `fetch_host_refs` on every fresh run.
    For a ~900 MB reference used repeatedly, download it once and point
    `host_genomes` at the local file instead.

The same reasoning applies to non-human hosts. An unmasked plant or insect genome
carries the same rRNA and low-complexity traps; if only an unmasked assembly is
available, the surviving read counts in `stats/read_tracking.txt` are worth watching
for an implausibly large drop at the `no_host` column — provided `remove_phix: true`.
With PhiX removal off there is no independent raw count: `raw` is then read from the
fastp report, i.e. from the reads that reached fastp *after* host removal, so it
equals `no_host` and the drop can never appear. In that case judge the host fraction
from BBMap's mapping rate in `logs/decontam_host/{sample}.log`.

### `host_min_id`

`host_min_id` is BBMap's `minid`, and the default is `0.95`. BBMap's own
documentation calls it an *approximate* minimum alignment identity: it tunes how hard
BBMap looks for an alignment rather than imposing a hard cut, so the realised identity
of the reads that get removed scatters around the value. Treat it as a sensitivity
dial, not a threshold.

- **Lower** (e.g. `0.90`) removes more: it catches host reads that diverge from the
  reference, at the cost of pulling in microbial reads that happen to align
  reasonably well.
- **Higher** (e.g. `0.98`) removes less: only near-exact host matches go, so more
  host survives into the table.

`0.95` is a middle setting that suits a reference from the same species as the host.
A reference that is only a relative of the true host argues for a lower value.

A methods section that needs a guaranteed cut-off — "host reads were removed at ≥95%
identity" — cannot get it from `minid`. BBMap's `idfilter` is the explicitly exact
parameter for that, and MetaFlux does not set it.

### Resources

The host steps are the memory-hungry part of preprocessing. The defaults are sized
for a masked human genome, whose BBMap index needs roughly 24–28 GB:

| Rule | `threads` | `mem_mb` |
|---|---:|---:|
| `build_host_index` | 16 | 30000 |
| `decontam_host` | 16 | 30000 |
| `decontam_phix` | 6 | 4000 |

Both BBMap rules pass 85% of their `mem_mb` to the JVM as `-Xmx`, leaving the rest
for the process itself. A smaller host genome can take a smaller budget; a larger one
must be given more, or the index build will run out of heap. BBDuk is capped at 6
threads on purpose — it scales poorly past four to six workers on a 5 kb reference,
so running several samples in parallel is the better use of cores.

## Mis-paired FASTQs stop both steps

BBDuk reads R1 and R2 in lockstep, one record from each at a time, and requires the
two files to hold the same number of records in the same order. Read length is not
part of that requirement: ragged read lengths, including mates of one pair ending up
different lengths, are processed normally.

So a crash at `decontam_phix` reading

```text
java.lang.AssertionError: List size mismatch: 200 vs 195
```

means the two files hold different numbers of reads — the pairing is broken. The
usual cause is a library that was filtered before deposition, leaving orphaned mates
behind in one file. SRA and BioProject submissions are where this turns up.

Turning PhiX removal off does not get round it. BBMap reads the pair in lockstep too
and makes the same requirement, so `decontam_host` stops on the same input, and says
so plainly:

```text
There appear to be different numbers of reads in the paired input files.
The pairing may have been corrupted by an upstream process.
It may be fixable by running repair.sh.
```

With both steps off the run stops at fastp instead, which requires the same thing of
the two files.

**Repair the pairing first.** `repair.sh` comes from the same `bbmap` conda package
MetaFlux already installs for these rules, so nothing new has to be added:

```bash
repair.sh in1=R1.fastq.gz in2=R2.fastq.gz \
  out1=fixed_R1.fastq.gz out2=fixed_R2.fastq.gz \
  outs=singletons.fastq.gz
```

It re-syncs the two files on read name and writes the leftover single-ended reads to
`outs=`. Point `input.fastq_dir` at the repaired files and run again.

Where the raw, untrimmed reads are still available they are the better starting
point — they are properly paired, and MetaFlux trims with fastp anyway, so
pre-trimmed input buys nothing.

## What to check afterwards

- `stats/read_tracking.txt` — the `nophix` and `no_host` columns, per sample. Both
  only appear when the step ran. With `remove_phix: false` the `raw` column is taken
  from the fastp report rather than from BBDuk, so it matches `no_host` exactly.
- `multiqc/multiqc_report.html` — the BBDuk PhiX statistics alongside fastp's. BBMap
  host removal does not appear in MultiQC: `{sample}_dehost_stats.txt` is a single
  `surviving_pairs` count written by MetaFlux, not a format MultiQC parses. BBMap's
  own mapping rate is in the log below.
- `logs/decontam_host/{sample}.log` and `logs/build_host_index.log` — BBMap's own
  output, including the mapping rate, if a sample looks wrong.
