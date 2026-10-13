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
to everything below it that appears in the Bracken report** are extracted, so
`Bacteria` yields the bacterial fraction of the sample and `Escherichia coli` that one
species.

That qualification matters. `--include-children` works out the descendant taxids from
the report file it is handed, and MetaFlux hands it the Bracken report, which is
pruned to `shotgun.bracken.tax_lev`. Sub-species and strain rows are therefore not in
it: under the default `tax_lev: S`, reads Kraken2 assigned to a strain taxid are left
behind. That is a real loss on RefSeq-derived indexes, where many genomes are
deposited under strain taxids — *Escherichia coli* str. K-12 substr. MG1655, taxid
511145, for instance — and Kraken2 routinely assigns reads to them. Lineages Bracken
dropped below `shotgun.bracken.threshold` are missing for the same reason.

At `tax_lev: G` or any coarser rank the loss is much larger, because the report then
has no species rows at all. A species name is not found, so the sample gets empty
files. A genus name pulls only the reads Kraken2 placed on the genus itself, not
those placed on its species, so it can miss most of the reads in the genus's clade.

Where those reads matter, run `extract_kraken_reads.py` by hand against
`02.classification/{sample}_report.txt`, which keeps the full tree. KrakenTools reads
that report as it stands. Note the report is written with `--report-minimizer-data`,
which adds two columns, so the taxid and name sit in fields 7 and 8 rather than the
usual 5 and 6 — which matters for anything else that picks fields out of it.

!!! warning "The name must match one report row exactly"

    The taxid is found by looking for the row of the sample's Bracken report whose
    name, without the indentation that shows its place in the tree, is exactly the
    name given, capital letters included. A partial name therefore finds nothing:
    `coli` does not match *Escherichia coli*, and `Prevotella` finds the genus, not
    the family Prevotellaceae listed above it. Write the full scientific name.

    Some names occur in more than one row of a report: a genus and a subgenus of the
    same name (*Leishmania*, for instance), or a placeholder name that GTDB gives to a
    family and to a genus alike. MetaFlux cannot tell which of them is meant, so the
    run stops with a message that lists the rows, with their rank codes and taxids.
    Remove the name from `extract_taxa`, and extract the taxon you mean by hand, as
    described above, with `extract_kraken_reads.py -t <taxid>`.

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
      there than on *E. coli* itself, while Bracken's estimate in the OTU table gives
      them back to *E. coli*.

    To collect reads for an assembly or a gene search, name the genus or a higher
    rank (under the default `tax_lev: S`; coarser settings lose reads, as described
    above). The direct-read column of `02.classification/{sample}_report.txt` (the
    third) shows on which rows Kraken2 placed them. Use a species-level extraction to
    check a doubtful species call instead: map the extracted reads to a genome of that
    species and look at how evenly they cover it. Reads spread over the whole genome
    support the call; reads piled up on a few regions, such as conserved genes or
    mobile elements, point to a relative or to a false positive. Either way it is a
    check, not proof.

## What it does

For each sample and each named taxon:

1. The taxon's NCBI taxid is looked up in that sample's Bracken report
   (`03.abundance/{sample}_report.txt`): the row whose name is exactly the requested
   name. The lookup runs with the Python of the KrakenTools environment, so it does
   not depend on the tools installed on the host.
2. `extract_kraken_reads.py` reads the Kraken2 per-read output
   (`02.classification/{sample}_output.txt`), which records the assignment of every
   read pair, together with the Bracken report so it can walk the taxonomy below the
   target taxid (`--include-children`).
3. The matching read pairs are written out as FASTQ and then compressed with `pigz` —
   `extract_kraken_reads.py` cannot write gzip itself, so a separate rule does it.

The reads come from the classified-read FASTQs Kraken2 wrote
(`02.classification/{sample}_R1.fastq.gz` and `_R2.fastq.gz`), which is why those
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
metagenome effectively writes a second copy of the classified library — worth
checking the disk has room.

## When a taxon is absent

If no row in a sample's Bracken report matches the name, the rule writes **empty**
FASTQ files and records the reason in
`logs/extract_taxon_reads/{sample}_{taxon}.log`:

```text
no taxid found for 'Escherichia coli'
```

Empty output is deliberate: Snakemake needs the declared files to exist for the run
to finish, and an absent taxon in one sample should not fail the whole workflow. An
unexpectedly empty file is a signal to check the log — either the taxon really is not
in that sample, or the name did not match anything in the report.

Since the lookup goes through the Bracken report, only taxa Bracken reported can be
extracted. A taxon that fell below `shotgun.bracken.threshold` will not be found —
see [Reference databases](databases.md#bracken-settings).

## Resources

| Rule | `threads` | `mem_mb` | What it does |
|---|---:|---:|---|
| `extract_taxon_reads` | 4 | 2000 | taxid lookup + read extraction |
| `compress_extracted` | 6 | *(not declared)* | `pigz` compression of the extracted FASTQs |

The uncompressed intermediates are marked temporary and are removed once compressed,
so only the `.fastq.gz` files remain.
