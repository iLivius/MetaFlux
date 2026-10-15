# Per-taxon read extraction via KrakenTools.
#
# wildcards.taxon is the slug (e.g. "bacteria"); we map it back to the real
# scientific name via params.target so we can look up its NCBI taxid in the
# sample's Kraken2 report. The name must be exactly the name of one report row, or of
# a genus and a subgenus of the same name, which count as the genus (the comment above
# rule extract_taxon_reads says why). It is matched at WHATEVER
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
# 02.classification/kraken2/{sample}_report.txt (rule kraken2, 35_kraken2.smk), both for
# the taxid lookup below and for KrakenTools (-r). extract_kraken_reads.py works out which
# taxa lie below the requested one (--include-children) only from the rows of the report
# it is handed, and it skips, without a warning, every read placed on a taxon that report
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
#                          an order, a family and a genus) or two taxa that share a
#                          name. One case is resolved, the genus with a subgenus of its
#                          own name: exactly one matching row has a plain rank code
#                          (one letter, G), and every other matching row has that
#                          letter with a number (G1, G2, ...) and lies inside the plain
#                          row's branch of the report (it is listed after it, before
#                          the next row at the plain row's level or above). Kraken2
#                          codes a node without a standard rank by the standard rank
#                          above it plus its distance from it, so G1 is a child of a
#                          genus: Plasmodium is G 5820 with G1 418103 below it, and
#                          Leishmania G 5658 with G1 38568. The plain row is then
#                          taken: --include-children extracts every row of its branch,
#                          the sub-rank rows with their reads among them. The lookup
#                          writes a note on this to the log. In the 331 Kraken2 reports
#                          MetaFlux was tested on (NCBI and GTDB databases), only
#                          Plasmodium and Leishmania are resolved this way, and the
#                          987 other names found on several rows still stop.
#                          Every other case stops: MetaFlux cannot tell which row is
#                          meant (taking the outermost row would silently turn a GTDB
#                          placeholder name into an order or family extraction). The
#                          lookup writes a message listing the rows to the terminal and
#                          to the log, and exits with an error. Bash strict mode, which
#                          Snakemake always switches on, then ends the job before any
#                          read file is written.
# The log of each job (logs/extract_taxon_reads/{sample}_{taxon}.log) gets the lookup's
# note or stop message, the taxid line, and what extract_kraken_reads.py prints, cut to
# its summary lines. KrakenTools 1.2.1 updates a progress counter for every 10,000
# lines of Kraken2's per-read output, every 1,000 reads of each FASTQ and every read it
# writes out. Each update starts with a carriage return, so on a terminal it overwrites
# the previous one, but in a file all of them stay: the log of a Bacteria extraction
# from a gut mock sample (5.1 million read pairs) reached 530 MB. A small filter between
# KrakenTools and the log keeps, for every line, only the text after its last carriage
# return, which is what a terminal shows when the line is done: the report parsed, the
# number of taxa searched, the per-read output lines read, the read IDs saved and found
# in R1 and in R2, the reads written, the files and the start and end times. Any error
# message passes through unchanged. Only the log changes: the reads go from KrakenTools
# straight to the output files. Bash strict mode (pipefail) still ends the job when
# KrakenTools fails.
# The lookup and the filter run in the Python that the krakentools env already contains
# (KrakenTools is a set of Python scripts), not in the host's awk, so they do not depend
# on the machine. Their code must contain no curly braces (Snakemake would read them as
# placeholders), no single quotes (they would end the bash string that holds it) and
# no three double quotes in a row (they would end the shell string of the rule, so
# their functions are explained in comments, not docstrings).
rule extract_taxon_reads:
    wildcard_constraints:
        sample = EXTRACT_SAMPLE_PATTERN,
        taxon  = EXTRACT_TAXON_PATTERN,
    input:
        # Kraken2's own report, not Bracken's (the comment above the rule says why).
        report = OUT / "02.classification" / "kraken2" / "{sample}_report.txt",
        kraken = OUT / "02.classification" / "kraken2" / "{sample}_output.txt",
        r1     = OUT / "02.classification" / "kraken2" / "{sample}_R1.fastq.gz",
        r2     = OUT / "02.classification" / "kraken2" / "{sample}_R2.fastq.gz",
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
        # Start the log empty. The lookup below may add a line to it (the note on a
        # genus taken together with its subgenus, or the message of a stop), and every
        # line after that is added below it.
        : > {log}
        # Find the NCBI taxid of the requested name in this sample's Kraken2 report:
        # the row whose name is exactly that name (the comment above the rule says why,
        # and what happens when no row or several rows carry it).
        taxid=$(python -c '
import sys

target = sys.argv[1]
report_path = sys.argv[2]
log_path = sys.argv[3]

# Read the report once, top to bottom. For every row keep its depth in the
# taxonomy tree (the name is indented by two spaces per level), in report order;
# a row whose name, without that indentation, is the requested name is also kept
# as (row number, depth, rank code, taxid).
row_depths = []
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
        name_field = fields[-1]
        row_name = name_field.strip()
        depth = (len(name_field) - len(name_field.lstrip(" "))) // 2
        row_number = len(row_depths)
        row_depths.append(depth)
        if row_name == target:
            matching_rows.append((row_number, depth, rank_code, row_taxid))


# For a name found on several rows: the row to extract, or None to stop.
# A row is returned only when exactly one of the rows has a plain rank code (one
# letter, such as G) and every other row has that same letter with a number (G1,
# G2, ...) and lies inside the branch of the plain row: it comes after the plain
# row and before the first later row at the depth of the plain row or above.
# That is a genus with a subgenus of the same name (Plasmodium: G 5820 with G1
# 418103 inside it). Extracting the plain row with its children then includes
# the reads of the other rows. (No docstring here: its three double quotes would
# end the shell string of the rule.)
def plain_row_holding_the_others(matching_rows, row_depths):
    plain_rows = []
    numbered_rows = []
    for row in matching_rows:
        if len(row[2]) == 1:
            plain_rows.append(row)
        else:
            numbered_rows.append(row)
    if len(plain_rows) != 1:
        return None
    plain_number, plain_depth, plain_code, plain_taxid = plain_rows[0]
    # Where the branch of the plain row ends: the first later row at its depth or
    # above it, or the end of the report.
    branch_end = len(row_depths)
    for later_number in range(plain_number + 1, len(row_depths)):
        if row_depths[later_number] <= plain_depth:
            branch_end = later_number
            break
    for row_number, depth, rank_code, row_taxid in numbered_rows:
        same_rank = rank_code[0] == plain_code and rank_code[1:].isdigit()
        inside_branch = plain_number < row_number < branch_end
        if not same_rank or not inside_branch:
            return None
    return plain_rows[0]


# Exactly one row: print its taxid, which the shell keeps as the taxid to extract.
if len(matching_rows) == 1:
    print(matching_rows[0][3])

# Several rows: take the plain row when it holds all the others, or stop.
if len(matching_rows) > 1:
    row_list = []
    other_codes = []
    for row_number, depth, rank_code, row_taxid in matching_rows:
        row_list.append("rank " + rank_code + " taxid " + row_taxid)
        if len(rank_code) > 1:
            other_codes.append(rank_code)
    rows_text = ", ".join(row_list)
    taken_row = plain_row_holding_the_others(matching_rows, row_depths)

    # The plain row holds the others: note it in the log, and print its taxid.
    if taken_row is not None:
        taken_code = taken_row[2]
        taken_taxid = taken_row[3]
        note = (
            "[MetaFlux] extract_taxa: " + repr(target) + " is the name of "
            + str(len(matching_rows)) + " rows in " + report_path + " (" + rows_text
            + "). Every other row (rank " + ", ".join(other_codes) + ") lies inside"
            + " the rank " + taken_code + " row, as a subgenus lies inside its genus,"
            + " so MetaFlux extracts the rank " + taken_code + " row (taxid "
            + taken_taxid + "), whose children include the other rows."
        )
        with open(log_path, "a") as log_file:
            log_file.write(note + "\n")
        print(taken_taxid)

    # Any other case: say which rows they are and what to do instead, then stop.
    if taken_row is None:
        message = (
            "[MetaFlux] extract_taxa: " + repr(target) + " is the name of "
            + str(len(matching_rows)) + " rows in " + report_path
            + " (" + rows_text + "), so MetaFlux cannot tell which taxon to"
            + " extract. Remove the name from shotgun.extract_taxa, and extract the taxon"
            + " you mean by hand with KrakenTools (extract_kraken_reads.py -t TAXID), as"
            + " the read extraction page of the docs describes."
        )
        with open(log_path, "a") as log_file:
            log_file.write(message + "\n")
        sys.stderr.write(message + "\n")
        sys.exit(1)

# No matching row: nothing is printed, so the taxid stays empty.
' "$target" {input.report} {log})
        echo "{wildcards.sample}: taxon='$target' taxid=$taxid" >> {log}
        if [ -n "$taxid" ]; then
            # Extract the reads, and pass what KrakenTools prints through a small
            # filter on its way to the log (the comment above the rule says why). The
            # reads go straight to the output files and never through the filter.
            extract_kraken_reads.py -k {input.kraken} \
                -1 {input.r1} -2 {input.r2} \
                -t "$taxid" --include-children --fastq-output \
                -o {output.r1} -o2 {output.r2} \
                -r {input.report} 2>&1 | python -c '
import sys

# Bytes in, bytes out: read as text, Python would turn every carriage return into
# a newline, and the progress counter would become one line per update.
source = sys.stdin.buffer
target = sys.stdout.buffer
carriage_return = b"\r"
newline = b"\n"

# Read the stream in blocks of 1 MiB, so that a progress line of several hundred
# MB is never held in memory whole. For every finished line, write only the text
# after its last carriage return: the final count. The unfinished end of the
# current line is cut the same way after every block, so it stays short.
unfinished = b""
while True:
    block = source.read(1048576)
    if not block:
        break
    pieces = (unfinished + block).split(newline)
    unfinished = pieces.pop()
    for piece in pieces:
        final_text = piece.rsplit(carriage_return, 1)[-1]
        target.write(final_text + newline)
    unfinished = unfinished.rsplit(carriage_return, 1)[-1]

# A last line without a newline (an error message, say) is kept too.
if unfinished:
    target.write(unfinished + newline)
' >> {log}
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
