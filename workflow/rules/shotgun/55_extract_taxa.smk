# Per-taxon read extraction via KrakenTools.
#
# wildcards.taxon is the slug (e.g. "bacteria"); we map it back to the real
# scientific name via params.target so we can look up its NCBI taxid in the
# sample's Kraken2 report. The name must be exactly the name of one report row (the
# comment above rule extract_taxon_reads says why). It is matched at WHATEVER
# rank it appears in the report (domain, genus, species, ...) — there is no rank
# restriction, so extract_taxa names a broad group (Bacteria) and a single
# species (Escherichia coli) the same way. --include-children then pulls every
# read under that taxid, not just reads classified exactly at it. The reads of a
# single species are only an approximate set (docs/shotgun/read-extraction.md
# says why and what they are still good for).
#
# Only activated when shotgun.extract_taxa is non-empty (the corresponding
# targets are only added to rule all in that case — see _shotgun_targets()
# in 00_common.smk).
#
# Terminal output: the extracted per-taxon FASTQs are not read by any other
# rule in this pipeline. They exist for whatever the user wants to do next
# with just that taxon's reads — assemble a genome from a single genus's
# reads, hand a species's reads to a dedicated strain-typing tool, and so on.

# The extracted files are named <sample>_<taxon slug>_R1.fastq(.gz), and both parts may
# contain '_' (the slug of "Escherichia coli" is escherichia_coli; a sample may be called
# Ctrl_2). From the file name alone Snakemake cannot tell where the sample name ends, so
# both wildcards are pinned to the values this run really has: the run's sample names
# (SAMPLES: the FASTQs found in fastq_dir, or the runs of input.bioproject) and the slugs
# of shotgun.extract_taxa, made exactly as _shotgun_targets() (00_common.smk) makes them
# for the file names.
_extract_sample_choices = []
for _sample_name in SAMPLES:
    _extract_sample_choices.append(re.escape(_sample_name))
_extract_taxon_choices = []
for _taxon_name in EXTRACT_TAXA:
    _extract_taxon_choices.append(re.escape(_taxon_name.lower().replace(" ", "_")))
EXTRACT_SAMPLE_PATTERN = "|".join(_extract_sample_choices)
EXTRACT_TAXON_PATTERN = "|".join(_extract_taxon_choices)


# Which report the rule reads: Kraken2's own report of the sample,
# 02.classification/{sample}_report.txt (rule kraken2, 35_kraken2.smk), both for the taxid
# lookup below and for KrakenTools (-r). extract_kraken_reads.py works out which taxa lie
# below the requested one (--include-children) only from the rows of the report it is
# handed, and it skips, without a warning, every read placed on a taxon that report
# lacks. Kraken2's report lists every taxon that got at least one read, strains
# included, so the read pairs extracted for a taxon are exactly the clade count of its
# row (column 2). The Bracken report (03.abundance/{sample}_report.txt) would lose
# reads, because it has no rows below the Bracken rank (at tax_lev S, no strain rows)
# and no row for a taxon that fell below Bracken's threshold. Many genomes of a
# database sit under a strain taxid: in the Zymo gut mock all 78,214 Candida albicans
# read pairs are on the strain SC5314, so an extraction through the Bracken report gave
# none of them. Kraken2's report also does not depend on the Bracken settings or on the
# k-mer evidence gate, so changing those never changes which reads are extracted.
#
# The taxid lookup at the start of the shell command below. extract_kraken_reads.py
# needs the NCBI taxid of the taxon, while the config gives its name, so the rule first
# finds that name in the same Kraken2 report. Its rows are tab-separated, and their last
# three columns are always the rank code, the taxid and the name. MetaFlux writes the
# report with --report-minimizer-data, which puts two minimizer columns before them
# (eight columns instead of the standard six), so the lookup counts the columns from the
# end and works on either layout. The name is indented by two spaces per level of the
# taxonomy tree.
# A row counts only when its name, once that indentation is stripped, is exactly the
# requested name. Taking the first row whose name merely CONTAINS it picks the wrong
# taxon whenever a longer name that contains it comes first in the report: Prevotella
# would find the family Prevotellaceae, Campylobacter the phylum Campylobacterota,
# Marinobacter the genus Marinobacterium, and Bacillus cereus the Bacillus cereus
# group. The test compares plain text, not a pattern, so a name with regex characters
# in it (e.g. "Bacteroides sp.") is matched literally.
#   one matching row       its taxid goes to KrakenTools;
#   no matching row        Kraken2 assigned no read of this sample to the taxon (or the
#                          name is misspelt), and the rule writes empty files (the
#                          else-branch of the shell command);
#   several matching rows  one name at two ranks (a genus and a subgenus of the same
#                          name, such as Plasmodium; a GTDB placeholder name used for
#                          a family and a genus) or two taxa that share a name.
#                          MetaFlux cannot tell which one is meant, so the lookup
#                          writes a message listing them to the terminal and to the
#                          log, and exits with an error. Bash strict mode, which
#                          Snakemake always switches on, then ends the job before any
#                          read file is written.
# The lookup runs in the Python that the krakentools env already contains (KrakenTools
# is a set of Python scripts), not in the host's awk, so it does not depend on the
# machine. Its code must contain no curly braces (Snakemake would read them as
# placeholders) and no single quotes (they would end the bash string that holds it).
rule extract_taxon_reads:
    wildcard_constraints:
        sample = EXTRACT_SAMPLE_PATTERN,
        taxon  = EXTRACT_TAXON_PATTERN,
    input:
        # Kraken2's own report, not Bracken's (the comment above the rule says why).
        report = OUT / "02.classification" / "{sample}_report.txt",
        kraken = OUT / "02.classification" / "{sample}_output.txt",
        r1     = OUT / "02.classification" / "{sample}_R1.fastq.gz",
        r2     = OUT / "02.classification" / "{sample}_R2.fastq.gz",
    output:
        r1 = temp(OUT / "04.extracted_reads" / "{sample}_{taxon}_R1.fastq"),
        r2 = temp(OUT / "04.extracted_reads" / "{sample}_{taxon}_R2.fastq"),
    params:
        target = lambda wc: {t.lower().replace(" ", "_"): t for t in EXTRACT_TAXA}[wc.taxon],
    threads: lambda wc: threads_for("extract_taxon_reads")
    resources:
        mem_mb = lambda wc: mem_mb_for("extract_taxon_reads"),
    conda:
        "../../envs/krakentools.yaml"
    log:
        LOGS / "extract_taxon_reads" / "{sample}_{taxon}.log"
    shell:
        r"""
        target="{params.target}"
        # Find the NCBI taxid of the requested name in this sample's Kraken2 report:
        # the row whose name is exactly that name (the comment above the rule says why,
        # and what happens when no row or several rows carry it).
        taxid=$(python -c '
import sys

target = sys.argv[1]
report_path = sys.argv[2]
log_path = sys.argv[3]

# Every report row whose name, without its indentation, is the requested name,
# kept as (rank code, taxid).
matching_rows = []
with open(report_path) as report:
    for line in report:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        # The last three columns are the rank code, the taxid and the name, in the
        # eight-column report Kraken2 writes here and in the standard six-column one.
        rank_code = fields[-3]
        row_taxid = fields[-2]
        row_name = fields[-1].strip()
        if row_name == target:
            matching_rows.append((rank_code, row_taxid))

# Exactly one row: print its taxid, which the shell keeps as the taxid to extract.
if len(matching_rows) == 1:
    only_rank_code, only_taxid = matching_rows[0]
    print(only_taxid)

# Several rows: say which rows they are and what to do instead, then stop.
if len(matching_rows) > 1:
    row_list = []
    for rank_code, row_taxid in matching_rows:
        row_list.append("rank " + rank_code + " taxid " + row_taxid)
    message = (
        "[MetaFlux] extract_taxa: " + repr(target) + " is the name of "
        + str(len(matching_rows)) + " rows in " + report_path
        + " (" + ", ".join(row_list) + "), so MetaFlux cannot tell which taxon to"
        + " extract. Remove the name from shotgun.extract_taxa, and extract the taxon"
        + " you mean by hand with KrakenTools (extract_kraken_reads.py -t TAXID), as"
        + " the read extraction page of the docs describes."
    )
    with open(log_path, "w") as log_file:
        log_file.write(message + "\n")
    sys.stderr.write(message + "\n")
    sys.exit(1)

# No matching row: nothing is printed, so the taxid stays empty.
' "$target" {input.report} {log})
        echo "{wildcards.sample}: taxon='$target' taxid=$taxid" > {log}
        if [ -n "$taxid" ]; then
            extract_kraken_reads.py -k {input.kraken} \
                -1 {input.r1} -2 {input.r2} \
                -t "$taxid" --include-children --fastq-output \
                -o {output.r1} -o2 {output.r2} \
                -r {input.report} >> {log} 2>&1
        else
            # Taxon absent from this sample — emit empty files to satisfy the DAG.
            echo "no taxid found for '$target'" >> {log}
            : > {output.r1}; : > {output.r2}
        fi
        """


# extract_kraken_reads.py doesn't write .gz natively, so we pigz afterwards.
# Kept in the kraken2 env since pigz is already there.
rule compress_extracted:
    wildcard_constraints:
        sample = EXTRACT_SAMPLE_PATTERN,
        taxon  = EXTRACT_TAXON_PATTERN,
    input:
        r1 = OUT / "04.extracted_reads" / "{sample}_{taxon}_R1.fastq",
        r2 = OUT / "04.extracted_reads" / "{sample}_{taxon}_R2.fastq",
    output:
        r1 = OUT / "04.extracted_reads" / "{sample}_{taxon}_R1.fastq.gz",
        r2 = OUT / "04.extracted_reads" / "{sample}_{taxon}_R2.fastq.gz",
    threads: lambda wc: threads_for("compress_extracted")
    conda:
        "../../envs/kraken2.yaml"
    log:
        LOGS / "compress_extracted" / "{sample}_{taxon}.log"
    shell:
        """
        pigz -p {threads} -c {input.r1} > {output.r1} 2> {log}
        pigz -p {threads} -c {input.r2} > {output.r2} 2>> {log}
        """
