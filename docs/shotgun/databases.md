# Reference databases (shotgun)

Amplicon mode fetches its reference databases on first use and caches them under
`refdb/`. Shotgun mode does not: the Kraken2 database is large, version-sensitive,
and usually shared between projects on a machine, so it is supplied rather than
downloaded by the workflow.

| Reference | Provisioning | Config key |
|---|---|---|
| Kraken2 / Bracken index | **user-supplied**, required | `references.kraken_db` |
| Host genome(s) | **user-supplied**, optional | `shotgun.decontamination.host_genomes` |
| PhiX | bundled inside the BBTools package — nothing to fetch | *(none)* |

PhiX is worth a note: the `references.phix` block in the config is used by amplicon
mode, which downloads the genome and builds a bowtie2 index from it. Shotgun mode
uses BBDuk's own bundled `phix174_ill.ref.fa.gz` instead, so no shotgun rule uses
that block. Keep it anyway: `references.phix.fasta` is read when the workflow is
parsed, in both modes.

## Pointing MetaFlux at the database

```yaml
references:
  kraken_db: /path/to/db/k2_pluspf      # [shotgun] — user-supplied
```

The path is checked when the workflow is parsed, before any job runs. A missing or
mistyped path stops the run immediately with:

```text
[MetaFlux] mode=shotgun requires kraken_db but not found: /path/to/db/k2_pluspf
```

The key points at the *directory*, not at a file. A Kraken2 database directory holds:

| File | Used by | What it is |
|---|---|---|
| `hash.k2d` | Kraken2 | the k-mer → taxon hash — by far the largest file, and the one that sets the memory requirement |
| `opts.k2d` | Kraken2 | the options the index was built with |
| `taxo.k2d` | Kraken2 | the taxonomy tree |
| `*.kmer_distrib` | Bracken | per-read-length k-mer distributions, one file per supported read length |

Bracken is pointed at the same directory, which is why the `.kmer_distrib` files must
be present alongside the index. Pre-built downloads from the catalog below include
them.

## Choosing an index

Pre-built indexes — with sizes, build dates and md5 checksums — are listed in Ben
Langmead's [catalog](https://benlangmead.github.io/aws-indexes/k2). A standard
choice is **PlusPF** — the Standard collection with protozoa and fungi added — at
roughly 100 GB for the full build.

The same collections are also published as size-capped builds — `pluspf_16gb`, for
example — which are reduced to fit a memory ceiling, and much larger ones such as
`core_nt` exist at the other end. All of them work. With `memory_mapping: false`,
MetaFlux reads the index size from disk and sizes its memory request accordingly; with
memory mapping on — the shipped default — it requests a fixed 20 GB of private
workspace and lets the kernel page the index in, so the size on disk is not consulted
at all. See [Memory](#memory).

The trade-off between them is biological, not technical. A capped index holds fewer
k-mers, so more reads go unclassified and species-level calls become less reliable —
on a shallow index it is often more honest to re-estimate at genus level with
`shotgun.bracken.tax_lev: G`.

Two practical rules when choosing:

- The database defines what can be found. Fungi and protozoa are absent from the
  Standard collection, so a soil or gut study that expects them needs PlusPF or
  wider.
- Record the exact dated build. Two runs against different builds of "PlusPF" are
  not comparable, and the build date is what belongs in a methods section.

!!! warning "A missing organism does not show up as unclassified"

    The unclassified fraction is the visible half of database incompleteness, and it is
    the reassuring half. The other half is silent.

    When an organism is absent from the index but a relative is present, its reads do not
    fail to classify. Their k-mers genuinely occur in that relative, so they map to it,
    the lowest common ancestor resolves within that relative's lineage, and the
    confidence score is **high** — because agreement with the assigned clade really is
    high. The result is a confident species-level call for an organism that was never in
    the sample.

    Raising `confidence` does not help here. The score measures how much of a read's
    evidence supports the taxon it was given, and nothing in that calculation can
    distinguish "this read is from *X*" from "this read is from an unsequenced relative
    of *X*". So a well-populated unclassified fraction is not evidence that the database
    was adequate, and a low one is not evidence that the species calls are right.

    The practical defences are to use an index that covers the expected community, to
    re-estimate at genus level where species representation is thin
    (`shotgun.bracken.tax_lev: G`), and to treat species-level calls from a sparse index
    as hypotheses rather than identifications.

## Downloading

A single `wget` stream works but is slow on a ~100 GB file. A multi-connection
downloader is faster and resumable:

```bash
mkdir -p /path/to/db && cd /path/to/db

# substitute the latest dated build from the catalog
aria2c -x16 -s16 https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260226.tar.gz

# verify integrity against the md5 listed in the catalog
md5sum k2_pluspf_20260226.tar.gz

# contains hash.k2d, opts.k2d, taxo.k2d, *.kmer_distrib
tar -xzf k2_pluspf_20260226.tar.gz
```

If the AWS CLI is installed, `aws s3 cp --no-sign-request <url>` is the fastest
option (parallel multipart). Plain `wget -c <url>` also works and resumes a partial
download.

!!! warning "Check the free space before extracting"

    The tarball and the unpacked index sit on disk at the same time, so a full PlusPF
    build needs room for both.

## Memory

This is the part of a shotgun run most likely to go wrong on a shared machine. The
whole hash has to be reachable while reads are being classified, and how it gets
there is set by one key:

```yaml
shotgun:
  kraken:
    memory_mapping: true      # true = mmap the DB; false = load it into RAM
```

Only one of the two strategies applies at runtime.

| | `memory_mapping: true` (default) | `memory_mapping: false` |
|---|---|---|
| How the index is reached | the kernel maps the file and pages it in on demand | each process reads the whole hash into its own memory |
| Concurrent Kraken2 jobs | share the same pages through the operating system's file cache, so parallel jobs do **not** multiply RAM use | each job holds a private full copy |
| Per-process budget MetaFlux requests | fixed 20 GB | `hash.k2d` size + 10% |
| Speed | comparable to an in-RAM run once the whole hash is cached; I/O-bound for the entire run if it never fits | fastest once loaded |

Memory mapping is the safer default, and the only workable one when several samples
are classified at the same time or when the index is larger than the machine's RAM.
Turning it off is worth it when the index comfortably fits in memory and samples are
processed one at a time.

The speed row hides two very different cases. Kraken2's lookups land at random across
the whole hash, so a memory-mapped run only catches up with an in-RAM one once the
entire file sits in the operating system's page cache — that is, when the machine has
appreciably more RAM than the index. If the index does not fit, pages are evicted and
read back continuously and classification stays I/O-bound from start to finish: hours
rather than minutes per sample on a spinning disk or a shared network filesystem.
Memory mapping makes an oversized index *possible*, not fast; where that is the
situation, a size-capped build that does fit is the better trade.

!!! note "How the kraken2 memory request is computed"

    The memory request for the `kraken2` rule is settled when the workflow is
    parsed: with memory mapping on it is a flat 20 GB; with it off, the real size of
    `hash.k2d` on disk plus 10% headroom, so a 16 GB capped index and a 100 GB PlusPF
    are handled without editing anything — or, with memory mapping off and `hash.k2d`
    not readable yet (the database is still downloading, say), a conservative
    ~110 GB. This value is what the rule uses, so the `resources.mem_mb.kraken2`
    entry in the config does not override it.

    On a cluster this matters: the value is what gets requested from the scheduler,
    so an underestimate is a killed job.

The first Kraken2 call of a run reads the whole database off disk and is slow while
the cache is cold; later samples are much faster once the operating system has the
pages in memory — provided, as above, that there is room to keep them there.

Threads for the classification step come from `resources.threads.kraken2` unless
pinned directly:

```yaml
shotgun:
  kraken:
    threads: null       # null = use resources.threads.kraken2 (16 by default)
```

## Bracken settings

Bracken uses the same database directory. Its two settings:

| Key | Default | What it does |
|---|---|---|
| `shotgun.bracken.tax_lev` | `S` | rank the abundances are re-estimated at — `D`, `P`, `C`, `O`, `F`, `G` or `S` |
| `shotgun.bracken.threshold` | `10` | minimum number of reads Kraken2 must already have put in a taxon's clade for it to enter the re-estimation |

`tax_lev: S` is the usual choice and gives a species-level table. Drop it to `G` when
the index is shallow, or when species-level assignments are not credible for the
group in question.

`threshold: 10` is applied to Kraken2's own clade read count for a taxon at the target
rank, before Bracken re-estimates anything. Taxa below it are dropped from the model
altogether and their reads are discarded rather than redistributed — Bracken prints
them in its log as "Total reads discarded".

It is therefore a floor on direct evidence, not on the final estimate, and that has
two consequences worth knowing. A species whose reads Kraken2 mostly parked at genus
level — the usual outcome where the database holds several close relatives — can be
removed even though redistribution would have credited it with thousands of reads.
And because the discarded reads leave the table rather than moving to a surviving
relative, the threshold shrinks the total, not just the tail.

The count is flat rather than proportional, so it also bites harder on shallowly
sequenced samples than on deep ones. Raise it for deep libraries; lower it if either
of the above is a concern, and expect a longer tail of one-read taxa.

The read length Bracken uses is not configurable: it is read per sample from the
fastp report and matched to the nearest distribution the database ships (50, 75, 100,
150, 200, 250 or 300 bp). Mixing libraries with very different read lengths in one
run is fine — each sample is handled on its own.

## Host genomes

The optional host reference is the other user-supplied database. It is configured
under `shotgun.decontamination.host_genomes` and covered in
[Decontamination](decontamination.md).
