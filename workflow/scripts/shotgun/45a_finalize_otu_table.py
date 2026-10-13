#!/usr/bin/env python3
"""Build the shotgun OTU table, otu_table.tsv and otu_table.biom, from the Bracken reports.

WHAT THIS SCRIPT MAKES
  One table: one row per taxon at the Bracken rank (shotgun.bracken.tax_lev, S by
  default), one column per sample, holding Bracken's read counts. Each row carries a
  lineage of eight fixed slots:

      d__ k__ p__ c__ o__ f__ g__ s__
      domain, kingdom, phylum, class, order, family, genus, species

  These are the eight rank codes Kraken2 prints in its reports (D K P C O F G S), so
  every Kraken2 table has the same eight slots, whatever the database. A rank that a
  taxon's lineage does not have leaves its slot empty: on a GTDB database, which has
  no kingdom rank, k__ is empty in every row; on an NCBI database k__ is empty for most
  protists, which NCBI places in no kingdom; at tax_lev G, s__ is empty.

WHERE EVERY VALUE COMES FROM
  Input: the Bracken report of each sample, 03.abundance/{sample}_report.txt, written
  by rule bracken (45_abundance.smk). A Bracken report has the layout of a Kraken2
  report: one row per taxon in tree order, the name indented by two spaces for every
  level below the root. Its six columns are: percentage of reads, clade reads (reads
  on the taxon and everything below it), reads on the taxon itself, rank code, taxid,
  name. Bracken has re-estimated the read numbers at its rank.
  - Counts: the clade reads of every row whose rank code is exactly the Bracken rank
    ("S", not "S1"), when above 0. These are Bracken's re-estimated reads, the numbers
    of {sample}_output.txt, and the numbers kraken-biom put in this table up to
    MetaFlux v2.4.0.
  - Lineage: while a report is read from top to bottom, the path from the root to the
    current row is kept. Each slot takes the name of the row on that path whose rank
    code is exactly the slot's letter. Rows with a numbered code (K1, P3, R2, ...) are
    in-between nodes such as "Opisthokonta" or "FCB group": they have no slot and are
    skipped. The path is cut back at every row, so a slot never keeps a name from a
    neighbouring branch of the tree. (kraken-biom did: depending on which report it
    read first, it gave Tepidimonas, which NCBI places in no family, the family of the
    branch listed before it, Alcaligenaceae.)
  - Species: the whole name of the species row, as the database writes it
    ("Escherichia coli", "Peduovirus P37", "uncultured Martelella sp."). kraken-biom
    cut the genus off species names. MetaFlux used to glue it back, which could not
    restore names such as "Peduovirus P37" (left as "P37") and doubled the genus in
    names such as "uncultured Martelella sp.".
  - Domain and kingdom: the rows coded D and K. Which taxa the database itself ranks
    as a domain or a kingdom is read from its taxonomy file (taxo.k2d) when the
    workflow starts (00_common.smk) and handed to this script as params.

VIRUSES: THE REALM GOES IN d__
  Viruses have no domain. ICTV, which names virus taxa, made the realm its top rank in
  place of the domain, and NCBI keeps the realm in the same "domain/realm" column of
  its taxonomy as the domain of cellular organisms. So the d__ slot of a virus holds
  its realm (e.g. d__Duplodnaviria), and its viral kingdom (e.g. Heunggongvirae) sits
  in k__ like any other kingdom. Kraken2 has no rank code for a realm: it prints a
  realm as R2, the code it also gives other unranked nodes at that depth
  ("unclassified viruses" is R2 too). A realm is therefore recognised by its rank in
  taxo.k2d (params.realms), not by its code. A virus without a realm (about 1 in 6
  NCBI virus species) keeps an empty d__, as does anything else without a domain,
  plasmids for example.
  Viruses itself (taxid 10239) is never written as a domain. NCBI ranked it a
  superkingdom before 2025, so on a database built from that older taxonomy Kraken2
  codes it D (params.virus_root then lists it). That D is skipped, and the virus gets
  its realm in d__, exactly as on newer databases.

CHECKS THAT STOP THE RUN
  - A rank code that disagrees with the database. Every taxon the database ranks as a
    domain must be coded D in the reports, every kingdom K, and every row coded D or K
    must be a domain or a kingdom of the database (Viruses on an older database
    excepted, see above). Kraken2 releases before 2.17 code NCBI's 2025 domains R1 or
    R2 instead of D: their reports would leave d__ empty in every row, with no other
    sign of trouble. Reports made with another database than references.kraken_db
    code other taxa D or K.
  - One taxon with two lineages. All reports come from one database, so every sample
    gives a taxon the same lineage; two lineages mean reports of different databases.
  - A filter that removes every taxon. A table that had rows before the filter and has
    none after it would look like a finished run with nothing in it. Almost always the
    tokens do not fit this table: a misspelt name, or a token written for another
    table, such as k__Bacteria where this table writes the domain as d__Bacteria. The
    message names the tokens and lists the d__ values the table holds. The amplicon
    filter stops in the same case. A table that is already empty before the filter (no
    sample has a taxon at the Bracken rank) is still written, with a warning in the log.

THE OPTIONAL FILTER (shotgun.taxonomy_filter)
  A keep or discard token such as "g__Homo" must equal one whole slot of the lineage,
  so it can never match inside a longer name. keep runs first (a row must match at
  least one keep token), then discard (a row matching any discard token is removed).
  It changes only which rows this table shows; the Bracken reports stay as they are.
  If it removes every row, the job stops (see CHECKS THAT STOP THE RUN above).

OUTPUT
  otu_table.tsv and otu_table.biom in 03.abundance/, written from the same rows and
  numbers, so the two never disagree. otu_table.tsv is a final target listed in
  _shotgun_targets() (00_common.smk), and otu_table.biom is written by the same job;
  nothing else in MetaFlux reads them. Rows keep the order kraken-biom gave them: taxa
  in order of first appearance, reading the samples in the order of params.samples.
  Tables made before and after this change can therefore be compared line by line.
"""
import os
import sys

import biom
import numpy
from biom.util import biom_open

# 'snakemake' is the object Snakemake hands to every script: this rule's input, output,
# params and log. (noqa: F821 tells a code checker that the name is not undefined.)
sm = snakemake  # noqa: F821

# ── Settings handed over by rule finalize_otu_table (45_abundance.smk) ────────────
# input.reports is expanded over SAMPLES, so report_paths[i] belongs to sample_names[i].
report_paths = []
for path in sm.input.reports:
    report_paths.append(str(path))
sample_names      = list(sm.params.samples)     # SAMPLES, sorted: the column order
bracken_rank      = str(sm.params.tax_lev)      # the rank of every row, e.g. "S"
# What the Kraken2 database ranks as a domain, a kingdom or a virus realm, and Viruses
# itself on a database built before 2025: {taxid: name} each, read from taxo.k2d.
domain_taxids     = dict(sm.params.domains)
kingdom_taxids    = dict(sm.params.kingdoms)
realm_taxids      = dict(sm.params.realms)
virus_root_taxids = dict(sm.params.virus_root)
filter_enabled    = bool(sm.params.filter_enabled)
filter_keep       = list(sm.params.filter_keep)
filter_discard    = list(sm.params.filter_discard)
out_biom          = str(sm.output.biom)
out_tsv           = str(sm.output.tsv)

# The eight rank codes that get a slot, in lineage order.
SLOT_CODES = ["D", "K", "P", "C", "O", "F", "G", "S"]

log_fh = open(sm.log[0], "w")


def logmsg(message):
    """Write one line to this rule's log, flushed at once so that the log is complete
    even when the script stops on one of its checks."""
    log_fh.write(f"[finalize_otu_table] {message}\n")
    log_fh.flush()


def stop(message):
    """Write the reason to the log, then end the job with an error. Snakemake prints the
    message in the terminal and the run stops; no table is written."""
    logmsg("ERROR: " + message)
    log_fh.close()
    sys.exit(f"[finalize_otu_table] {message}")


def read_report(report_path):
    """Read one Bracken report and return its rows, from top to bottom.

    Each row becomes a tuple (depth, rank_code, taxid, name, clade_reads):
      depth        how many levels below the root the row sits (indentation / 2);
      rank_code    "D", "K", ..., "S", or a numbered code such as "P3";
      taxid        as text, as the report writes it;
      name         without its indentation;
      clade_reads  reads on the taxon and everything below it (column 2).
    Stops on a line without the six columns Bracken writes.
    """
    rows = []
    with open(report_path) as handle:
        line_number = 0
        for line in handle:
            line_number += 1
            line = line.rstrip("\n")
            if line == "":
                continue
            fields = line.split("\t")
            if len(fields) != 6:
                stop(f"{report_path}, line {line_number}: {len(fields)} columns, but a "
                     f"Bracken report has 6.")
            raw_name = fields[5]
            indentation = len(raw_name) - len(raw_name.lstrip(" "))
            depth = indentation // 2
            rank_code = fields[3].strip()
            taxid = fields[4].strip()
            name = raw_name.strip()
            clade_reads = int(fields[1])
            rows.append((depth, rank_code, taxid, name, clade_reads))
    return rows


def check_rank_code(report_name, rank_code, taxid, name):
    """Stop when the code a report gives a taxon disagrees with the rank the Kraken2
    database gives it (see CHECKS THAT STOP THE RUN above). Only D and K are checked,
    because only the domains and kingdoms of the database are read from taxo.k2d; that
    is enough to catch both causes, an old Kraken2 and reports of another database."""
    # Viruses on a database built before 2025: ranked superkingdom, so coded D, and
    # never counted as a domain (see VIRUSES above).
    if taxid in virus_root_taxids:
        return
    rerun_hint = ("Run Kraken2 again on references.kraken_db with the release MetaFlux "
                  "pins (workflow/envs/kraken2.yaml), for example with -R kraken2.")
    if taxid in domain_taxids and rank_code != "D":
        stop(f"{report_name}: {name} (taxid {taxid}) is a domain in the Kraken2 database "
             f"set as references.kraken_db, but the report codes it {rank_code}, not D. "
             f"Kraken2 releases before 2.17 write reports like this (they do not know the "
             f"rank name 'domain'), and so do reports made with another database. "
             f"{rerun_hint}")
    if rank_code == "D" and taxid not in domain_taxids:
        stop(f"{report_name}: {name} (taxid {taxid}) is coded D, but it is not a domain "
             f"in the Kraken2 database set as references.kraken_db. Were these reports "
             f"made with another database? {rerun_hint}")
    if taxid in kingdom_taxids and rank_code != "K":
        stop(f"{report_name}: {name} (taxid {taxid}) is a kingdom in the Kraken2 database "
             f"set as references.kraken_db, but the report codes it {rank_code}, not K. "
             f"Were these reports made with another database? {rerun_hint}")
    if rank_code == "K" and taxid not in kingdom_taxids:
        stop(f"{report_name}: {name} (taxid {taxid}) is coded K, but it is not a kingdom "
             f"in the Kraken2 database set as references.kraken_db. Were these reports "
             f"made with another database? {rerun_hint}")


def lineage_of(path):
    """The eight slots ('d__Bacteria', 'k__Pseudomonadati', ..., 's__Escherichia coli')
    of the taxon at the end of `path`, the list of (rank_code, taxid, name) of every row
    from the root down to it. A rank missing from the path gives an empty slot ('o__')."""
    name_for_code = {}
    for rank_code, taxid, name in path:
        if rank_code not in SLOT_CODES:
            continue      # an in-between node (K1, P3, R2, ...): no slot
        if taxid in virus_root_taxids:
            continue      # Viruses, coded D on an older database: never a domain
        name_for_code[rank_code] = name

    # No domain on the path: a virus, or something else without a domain (a plasmid,
    # for example). The d__ slot then takes the virus realm on the path, if there is one.
    if "D" not in name_for_code:
        for rank_code, taxid, name in path:
            if taxid in realm_taxids:
                name_for_code["D"] = name
                break

    slots = []
    for code in SLOT_CODES:
        slot_name = name_for_code.get(code, "")
        slots.append(code.lower() + "__" + slot_name)
    return slots


def matches_any(lineage_slots, tokens):
    """True when at least one filter token equals one of the lineage's slots."""
    for token in tokens:
        if token in lineage_slots:
            return True
    return False


def passes_filter(lineage):
    """True when a row with this lineage survives keep (it must match one keep token,
    if any are given), then discard (it must match no discard token)."""
    lineage_slots = set(lineage)
    if filter_keep and not matches_any(lineage_slots, filter_keep):
        return False
    if matches_any(lineage_slots, filter_discard):
        return False
    return True


def taxa_count_text(count):
    """A number of table rows in words: '1 taxon', '2 taxa', '1,394 taxa'."""
    if count == 1:
        return "1 taxon"
    return f"{count:,} taxa"


def as_yaml_list(tokens):
    """Filter tokens written as in the run config, '[d__Bacteria, d__Archaea]' or '[]',
    so that the stop message can be compared with the config at a glance."""
    return "[" + ", ".join(tokens) + "]"


def d_values_text(taxids, lineages):
    """The d__ values of these table rows, each with its number of rows, for the stop
    message when the filter removes every taxon. `lineages` maps each taxid to its
    eight slots (lineage_for, filled in block 1 below). Example:
        'd__Bacteria (95 taxa), d__Eukaryota (6 taxa), an empty d__ (2 taxa)'
    Rows without a domain or a virus realm (viruses without a realm, plasmids) have an
    empty d__ slot; it is listed last and named in words, because 'd__' alone is easy
    to miss in a list."""
    rows_per_value = {}
    for taxid in taxids:
        d_slot = lineages[taxid][0]           # the first of the eight slots is d__
        if d_slot not in rows_per_value:
            rows_per_value[d_slot] = 0
        rows_per_value[d_slot] += 1

    pieces = []
    for d_slot in sorted(rows_per_value):
        if d_slot == "d__":
            continue                          # the empty slot, added last below
        pieces.append(f"{d_slot} ({taxa_count_text(rows_per_value[d_slot])})")
    if "d__" in rows_per_value:
        pieces.append(f"an empty d__ ({taxa_count_text(rows_per_value['d__'])})")
    return ", ".join(pieces)


# ── 1. Read every report: the counts and the lineage of every taxon ──────────────
# The reports are read in the order of sample_names, so taxa enter row_order in the
# order of their first appearance, sample by sample: the row order kraken-biom used.
row_order = []          # taxids, in the order the rows will be written
reads_of = {}           # taxid -> {sample: Bracken's reads}
lineage_for = {}        # taxid -> its eight slots
first_report_of = {}    # taxid -> the report its lineage was first read from

for sample, report_path in zip(sample_names, report_paths):
    report_name = os.path.basename(report_path)
    rows = read_report(report_path)
    path = []           # (rank_code, taxid, name) of every row from the root down
    taxa_at_rank = 0
    for depth, rank_code, taxid, name, clade_reads in rows:
        check_rank_code(report_name, rank_code, taxid, name)

        # The rows above this one's depth are its ancestors; the rest of the path
        # belonged to the previous branch of the tree.
        ancestors = path[:depth]
        path = ancestors + [(rank_code, taxid, name)]

        # Only the rows at the Bracken rank that have reads become table rows.
        if rank_code != bracken_rank:
            continue
        if clade_reads <= 0:
            continue
        taxa_at_rank += 1

        lineage = lineage_of(path)
        if taxid not in reads_of:
            row_order.append(taxid)
            reads_of[taxid] = {}
            lineage_for[taxid] = lineage
            first_report_of[taxid] = report_name
        elif lineage != lineage_for[taxid]:
            stop(f"{name} (taxid {taxid}) has two lineages: "
                 f"{';'.join(lineage_for[taxid])} in {first_report_of[taxid]}, and "
                 f"{';'.join(lineage)} in {report_name}. All reports must come from the "
                 f"same Kraken2 database (references.kraken_db).")
        reads_of[taxid][sample] = clade_reads
    logmsg(f"{sample}: {taxa_at_rank} taxa at rank {bracken_rank}")

logmsg(f"{len(row_order)} taxa x {len(sample_names)} samples before the taxon filter")
if not row_order:
    logmsg(f"WARNING: no sample has a taxon at rank {bracken_rank}; writing an empty table")

# ── 2. The optional keep/discard filter (shotgun.taxonomy_filter) ────────────────
# See THE OPTIONAL FILTER above. The tokens come from the run config unchanged.
# A table that had rows and lost every one of them to the filter stops the job here,
# before anything is written (see CHECKS THAT STOP THE RUN above): an empty table
# would look like a finished run. A table that had no rows to begin with has already
# been warned about above, and is written empty.
if filter_enabled and (filter_keep or filter_discard):
    kept_taxids = []
    for taxid in row_order:
        if passes_filter(lineage_for[taxid]):
            kept_taxids.append(taxid)
    logmsg(f"filter ON  (keep: {as_yaml_list(filter_keep)}, "
           f"discard: {as_yaml_list(filter_discard)}): "
           f"{len(kept_taxids)} / {len(row_order)} OTUs retained")
    if row_order and not kept_taxids:
        stop(f"shotgun.taxonomy_filter removed every taxon of the table "
             f"({taxa_count_text(len(row_order))}; keep: {as_yaml_list(filter_keep)}, "
             f"discard: {as_yaml_list(filter_discard)}). The d__ values in this table "
             f"are: {d_values_text(row_order, lineage_for)}. A token must equal one whole "
             f"slot of the lineage, prefix included (d__ k__ p__ c__ o__ f__ g__ s__). "
             f"Check the tokens, or set shotgun.taxonomy_filter.enabled: false.")
else:
    kept_taxids = list(row_order)
    logmsg("filter OFF (shotgun.taxonomy_filter disabled or empty)")

# ── 3. The count matrix, one row per kept taxon and one column per sample ────────
# A taxon that a sample does not have gets 0 there. The counts are stored as
# floating-point numbers, as biom-format stores every BIOM matrix, so the TSV writes
# 726906.0, exactly as the tables made with kraken-biom did.
count_rows = []
for taxid in kept_taxids:
    counts = []
    for sample in sample_names:
        counts.append(float(reads_of[taxid].get(sample, 0)))
    count_rows.append(counts)

# ── 4. otu_table.biom ────────────────────────────────────────────────────────────
# Rows are the taxids, columns the samples. The eight slots are stored as the
# 'taxonomy' of each row, where phyloseq (import_biom) and mia look for the lineage.
# Each sample gets the metadata field 'Id' with its name, as in the tables made with
# kraken-biom (which wrote '<sample>_report' there). An empty table needs its shape
# spelled out: numpy cannot tell the number of columns from an empty list of rows.
if count_rows:
    matrix = numpy.array(count_rows, dtype=float)
else:
    matrix = numpy.zeros((0, len(sample_names)), dtype=float)
observation_metadata = []
for taxid in kept_taxids:
    observation_metadata.append({"taxonomy": lineage_for[taxid]})
sample_metadata = []
for sample in sample_names:
    sample_metadata.append({"Id": sample})
table = biom.Table(matrix, kept_taxids, sample_names,
                   observation_metadata=observation_metadata,
                   sample_metadata=sample_metadata,
                   type="OTU table")
# The second argument of to_hdf5 is a free-text "generated by" tag stored in the file;
# nothing reads it back. The file also records its creation time, so two runs give
# different bytes: compare BIOM files by content.
with biom_open(out_biom, "w") as handle:
    table.to_hdf5(handle, "MetaFlux finalize_otu_table")

# ── 5. otu_table.tsv, from the same rows and numbers ─────────────────────────────
# The layout this table has always had: a comment line, the '#OTU ID' header with the
# sample names and 'taxonomy', then one line per taxon, its slots joined by ';'.
# The first line, '# Constructed from biom file', is kept word for word as biom-format's
# own TSV export writes it, although this TSV is written here directly and not converted
# from the BIOM file, so that old and new tables compare line by line.
with open(out_tsv, "w") as handle:
    handle.write("# Constructed from biom file\n")
    handle.write("#OTU ID\t" + "\t".join(sample_names) + "\ttaxonomy\n")
    for taxid, counts in zip(kept_taxids, count_rows):
        count_texts = []
        for value in counts:
            count_texts.append(f"{value}")
        lineage_text = ";".join(lineage_for[taxid])
        line_fields = [taxid] + count_texts + [lineage_text]
        handle.write("\t".join(line_fields) + "\n")

logmsg(f"wrote {out_biom} and {out_tsv}")
log_fh.close()
