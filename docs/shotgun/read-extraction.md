# Per-taxon read extraction

A classification table says which organisms are present and in what proportion. It
does not give the reads back. When the next question is about sequence rather than
abundance — assemble the *Pseudomonas* fraction, look for a resistance gene in the
*Escherichia coli* reads, hand the eukaryotic fraction to a different tool — the
reads belonging to one taxon have to be pulled out of the library.

That is what this optional step does, using
[KrakenTools](https://github.com/jenniferlu717/KrakenTools). It is off unless taxa
are named:

!!! note "Extracted read counts will not match the OTU table"

    The two numbers answer different questions and are expected to differ, usually with
    fewer reads extracted than the OTU table suggests. Extraction returns the reads
    **Kraken2 actually assigned** to that clade. The OTU table holds **Bracken's
    estimates** after redistributing reads from coarser ranks — an estimate of origin,
    not a set of reads that can be handed back. A taxon credited with 50,000 reads in
    `otu_table.tsv` may yield considerably fewer here, and neither figure is wrong.

```yaml
shotgun:
  extract_taxa: []          # empty = skip extraction entirely
```

```yaml
shotgun:
  extract_taxa: [Bacteria, Escherichia coli]
```

With an empty list the extraction rules are never added to the run, so nothing is
computed and no `04.extracted_reads/` directory appears.

## Naming taxa

Entries are **NCBI scientific names**, written as they appear in the classification
report — `Bacteria`, `Fungi`, `Pseudomonas`, `Escherichia coli`.

The name is matched at whatever rank it occurs in the report. There is no rank
setting and no restriction: a domain, a phylum, a genus and a single species are all
named the same way and handled identically. Reads assigned to the named taxon **and
to everything below it** are extracted, so `Bacteria` yields the bacterial fraction of
the sample and `Escherichia coli` that one species, its strains included.

Both the lookup and the walk down the tree (`--include-children`) use the sample's own
Kraken2 report, `02.classification/kraken2/{sample}_report.txt`. It lists every taxon that
Kraken2 assigned at least one read to, strains and other ranks below species included,
so the read pairs extracted for a taxon are exactly its clade count in that report (the
second column). The Bracken settings (`tax_lev`, `threshold`) and the
[k-mer evidence gate](kmer-evidence.md) do not change this report, so they do not
change which reads are extracted.

Up to v2.4.0 MetaFlux used the Bracken report here, which has no rows below `tax_lev`
and none for the taxa Bracken dropped below its threshold. The reads Kraken2 had placed
on those rows were left out, without a warning: for a species whose genomes sit under a
strain taxid, often all of them. In the ZymoBIOMICS Gut Microbiome Standard (D6331), SRA
run SRR38976715 (BioProject PRJNA1474435; see
[Benchmark data](../about/test-datasets.md#benchmark-data)), Kraken2 with PlusPF
`k2_pluspf_20260626` placed all 29,618 *Candida albicans* read pairs on the strain
SC5314 and none on the species itself. An extraction based on the Bracken report would
have found none of them.

!!! warning "The name must match one report row exactly"

    The taxid is found by looking for the row of the sample's Kraken2 report whose
    name, without the indentation that shows its place in the tree, is exactly the
    name given, capital letters included. A partial name therefore finds nothing:
    `coli` does not match *Escherichia coli*, and `Prevotella` finds the genus, not
    the family Prevotellaceae listed above it. Write the full scientific name.

    Some names occur in more than one row of a report. A genus and a subgenus of the
    same name, such as *Plasmodium* and *Leishmania*, are resolved: the subgenus row
    (rank code `G1`) lies inside the genus row (`G`), so MetaFlux extracts the genus,
    whose reads include those of the subgenus, and writes a note in the log. The rule
    is narrow: exactly one of the rows has a plain rank code, and every other row has
    that letter with a number and lies inside it. In any other case, such as a
    placeholder name that GTDB gives to an order, a family and a genus alike, or two
    taxa that share a name, MetaFlux cannot tell which of them is meant, so the run
    stops with a message that lists the rows, with their rank codes and taxids. Remove
    the name from `extract_taxa`, and extract the taxon you mean by hand with
    KrakenTools (the bioconda package `krakentools`), giving it the taxid from the
    message:

    ```bash
    extract_kraken_reads.py -t <taxid> --include-children \
        -k 02.classification/kraken2/<sample>_output.txt \
        -r 02.classification/kraken2/<sample>_report.txt \
        -1 02.classification/kraken2/<sample>_R1.fastq.gz \
        -2 02.classification/kraken2/<sample>_R2.fastq.gz \
        --fastq-output -o <taxon>_R1.fastq -o2 <taxon>_R2.fastq
    ```

    KrakenTools reads the report as it is. Anything else that picks fields out of it
    by position must allow for the two columns that `--report-minimizer-data` adds:
    the taxid and the name sit in fields 7 and 8, not in the usual 5 and 6.

!!! warning "A single species gives an approximate read set"

    Kraken2 assigns each read on its own, from the k-mers it shares with the genomes
    in the database, so the reads extracted for one species are not exactly that
    species' reads:

    - reads of close relatives get mixed in, above all of relatives the database
      lacks: their reads land on the nearest species it holds;
    - part of the species' own reads are missing, often most of them. Kraken2 places
      a read that matches several taxa on the lowest node they share, and extraction
      takes only the reads placed on the named taxon or below it. *Escherichia coli*
      reads that also match *Shigella*, a separate genus in NCBI's taxonomy, are
      placed on the family Enterobacteriaceae. There can be many times more of them
      there than on *E. coli* itself (about 30 times more in the D6331 run
      SRR38976715), while Bracken's estimate in the OTU table gives them back to
      *E. coli*.

    To collect reads for an assembly or a gene search, name the lowest rank that holds
    them: often the genus, but for *E. coli*, whose reads shared with *Shigella* sit on
    the family, Enterobacteriaceae, which also brings in the reads of its other genera.
    The direct-read column of `02.classification/kraken2/{sample}_report.txt` (the third) shows
    on which rows Kraken2 placed them. Use a species-level extraction to check a
    doubtful species call instead: map the extracted reads to a genome of that species
    and look at how evenly they cover it. Reads spread over the whole genome support the
    call; reads piled up on a few regions, such as conserved genes or mobile elements,
    point to a relative or to a false positive. Either way it is a check, not proof.

## What it does

For each sample and each named taxon:

1. The taxon's NCBI taxid is looked up in that sample's Kraken2 report
   (`02.classification/kraken2/{sample}_report.txt`): the row whose name is exactly the
   requested name. The lookup runs with the Python of the KrakenTools environment, so
   it does not depend on the tools installed on the host.
2. `extract_kraken_reads.py` reads the Kraken2 per-read output
   (`02.classification/kraken2/{sample}_output.txt`), which records the assignment of every
   read pair, together with the same Kraken2 report so it can walk the taxonomy below
   the target taxid (`--include-children`).
3. The matching read pairs are written out as FASTQ and then compressed with `pigz` —
   `extract_kraken_reads.py` cannot write gzip itself, so a separate rule does it.

The reads come from the classified-read FASTQs Kraken2 wrote
(`02.classification/kraken2/{sample}_R1.fastq.gz` and `_R2.fastq.gz`), which is why those
files are kept rather than deleted. Reads Kraken2 could not classify are not in that
pool and can never be extracted — extraction only ever pulls out reads that already
carry a Kraken2 assignment.

## Output

```text
out_dir/04.extracted_reads/
├── {sample}_bacteria_R1.fastq.gz
├── {sample}_bacteria_R2.fastq.gz
├── {sample}_escherichia_coli_R1.fastq.gz
└── {sample}_escherichia_coli_R2.fastq.gz
```

The taxon part of the filename is the configured name lowercased with spaces
replaced by underscores, so `Escherichia coli` becomes `escherichia_coli`. Pairing is
preserved: R1 and R2 stay in step and can go straight into an assembler or a mapper.

One pair of files is produced per sample per taxon, so three taxa across forty
samples is 240 files. Naming a broad group such as `Bacteria` in a host-depleted
metagenome writes most of the classified reads out a second time, so check the free
disk space first.

## When a taxon is absent

If no row in a sample's Kraken2 report matches the name, the rule writes **empty**
FASTQ files and records the reason in
`logs/extract_taxon_reads/{sample}_{taxon}.log`:

```text
no taxid found for 'Escherichia coli'
```

Empty output is deliberate: Snakemake needs the declared files to exist for the run
to finish, and an absent taxon in one sample should not fail the whole workflow. An
unexpectedly empty file is a signal to check the log — either Kraken2 assigned no read
of that sample to the taxon, or the name did not match anything in the report.

The OTU table does not decide what can be extracted. A taxon that Bracken dropped for
falling below `shotgun.bracken.threshold` (see
[Reference databases](databases.md#bracken-settings)) has no row in `otu_table.tsv`,
but its reads are in Kraken2's report, and they are extracted.

## Resources

| Rule | `threads` | `mem_mb` | What it does |
|---|---:|---:|---|
| `extract_taxon_reads` | 4 | 2000 | taxid lookup + read extraction |
| `compress_extracted` | 6 | *(not declared)* | `pigz` compression of the extracted FASTQs |

The uncompressed intermediates are marked temporary and are removed once compressed,
so only the `.fastq.gz` files remain.
