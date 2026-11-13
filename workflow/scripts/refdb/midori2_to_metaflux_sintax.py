#!/usr/bin/env python3
"""Reference-DB prep (workflow/scripts/refdb/): MIDORI2 COI SINTAX file -> MetaFlux SINTAX file.

One of the reference-preparation scripts. Each script here turns one upstream
database into the exact format a classifier wants; the rules in 10_refdb.smk just
fetch and call these. See workflow/scripts/refdb/README.md for the full list.

WHY THIS EXISTS
MIDORI2 (Leray et al. 2022) publishes its COI reference for SINTAX with headers like

  >KU746593.1.<1.>658;tax=k:Eukaryota_2759,p:Arthropoda_6656,c:Insecta_50557,o:Diptera_7147,f:Sarcophagidae_7381,g:Sarcophaga_7384,s:Sarcophaga munronis_1806722;

Five things in the file do not fit how MetaFlux classifies and reports:
  1. The DOMAIN (Eukaryota) sits under the letter k:. VSEARCH reads the letters as
     a fixed top-down order (d above k above p ...), and MetaFlux's rank model
     writes d__ for a domain and k__ only for a kingdom taken from the reference's
     own taxonomy (16S on SILVA 144, ITS on UNITE). MIDORI2's SINTAX file has no
     kingdom at all.
  2. Every name ends with its NCBI taxid (Arthropoda_6656). The taxid is useful
     (it is how the kingdom is found, below) but it is not part of the name: kept,
     it would end up in every taxonomy string and every filter token
     (keep: [k__Metazoa] would have to read k__Metazoa_33208).
  3. Where NCBI has no name at a rank, MIDORI2 fills the slot with a PLACEHOLDER
     made from the names below it: order_Lyonsiidae_120440, genus_Lepismatidae
     sp. ..., phylum_class_order_family_genus_Amoebozoa sp._1892891. A placeholder
     is not a taxon, but SINTAX would report it as one, and in the NUC_SP variant
     MetaFlux uses, half of the records added over NUC carry a placeholder genus
     (287,117 pseudo-genera against 58,627 real genera in GB273).
  4. The NUC_SP variant includes records whose species is an INTERIM name
     ("Sarcophaga sp. BIOUG12345-A01", "cf.", "aff.", BOLD codes). Reported as a
     species, such a name would look like an identification it is not.
  5. Some MIDORI2 records are bacterial COI homologues filed under a eukaryote's
     name (mostly Rickettsiales symbionts of insects; see the MetaFlux manual page
     "Bacterial COI homologues (reference set)"). Left in, they would pull the
     bacterial ASVs into animal taxa.

WHAT IT DOES (record by record, in file order; sequences are copied unchanged)
  - writes the domain under d: (always Eukaryota in MIDORI2);
  - adds the KINGDOM under k:, taken from MIDORI2's own NCBI taxonomy snapshot
    (nodes.dmp and names.dmp of MIDORI2_<release>_KRAKEN2_taxonomy.zip): from the
    record's species taxid, the first ancestor of rank "kingdom" (Metazoa, Fungi,
    Viridiplantae). Protists have no kingdom in NCBI, so their records get no k:
    field at all (SINTAX and the reference check both accept a missing rank);
  - removes the _<taxid> suffix from every name;
  - drops every placeholder rank (a value that starts with phylum_, class_,
    order_, family_ or genus_) and every interim species name (a species that is
    not a formal binomial, see FORMAL_SPECIES): the field is left out of the
    header, not written empty, because the reference check refuses empty names;
  - writes species as Genus_epithet (space -> underscore), as UNITE and PR2 do;
  - leaves out every record whose SEQUENCE is marked "removed" in the MIDORI2 list
    that ships with the set of bacterial COI homologues
    (bacterial_cox1_<set>_midori2_<release>_records.tsv). The list names each
    record together with the SHA-1 of its sequence (column sequence_sha1: upper
    case, line breaks removed), and the decision is made by sequence: MIDORI2 files
    the same sequence under several names, and a bacterial sequence must go under
    all of them, or the copies left behind would still pull bacterial ASVs into
    animal taxa. The build of the set decided per sequence too, so the list marks
    every copy it knows; a copy it does not name is left out all the same and
    counted in the report.

WHAT IT REFUSES (the script stops with a message rather than guess)
  - a header without ";tax=", or whose first field is not k:Eukaryota_<taxid>;
  - a species field without a numeric taxid suffix;
  - a taxid that is not in nodes.dmp;
  - a removal list without the column sequence_sha1 (a list of the bacterial set
    before gtdb-r232.v2, which removed records by name);
  - a removal list made from this very MIDORI2 file (its column
    midori2_file_sha256 equals the SHA-256 of the input file) whose removed
    sequences are not all in the file: the list and the file do not match.
  A list made from another MIDORI2 file does not stop the run: its removed
  sequences that are missing from this file are counted in the log and the report.

INPUT  : input.fasta   = MIDORI2_UNIQ_NUC_SP_<release>_CO1_SINTAX.fasta.gz (rule
                         fetch_midori2_coi extracts it from MIDORI2_<release>_SINTAX_sp.zip)
         input.nodes   = nodes.dmp, input.names = names.dmp (the taxonomy zip)
         input.removal = bacterial_cox1_<set>_midori2_<release>_records.tsv (the
                         MIDORI2 list of the set of bacterial COI homologues)
OUTPUT : output.fasta  = the MetaFlux SINTAX file (symbol midori2_coi_sintax in
                         COI.yaml), read by rule check_taxonomy_reference and then by
                         assign_taxonomy (VSEARCH --sintax) and by the COI screens
                         (identity E of the off-target guard)
         output.report = counts of everything this script changed (one line per
                         count, tab-separated), for a human and for the docs;
                         nothing downstream parses it
LOG    : log[0]        = a plain-text log for a human to check after the run: the
                         files read and written and the main counts, or the message
                         the script stopped with; nothing downstream reads it
"""
import gzip
import hashlib
import re
import sys
from collections import Counter
from pathlib import Path

sm = snakemake  # noqa: F821 (injected by Snakemake)

in_path = Path(sm.input.fasta)
nodes_path = Path(sm.input.nodes)
names_path = Path(sm.input.names)
removal_path = Path(sm.input.removal)
out_path = Path(sm.output.fasta)
report_path = Path(sm.output.report)
log_path = Path(sm.log[0])

# MIDORI2's SINTAX letters, in the order MIDORI2 writes them. k holds the domain.
MIDORI2_LETTERS = ["k", "p", "c", "o", "f", "g", "s"]
# Placeholders start with the name of the rank they stand in for.
PLACEHOLDER_PREFIXES = ("phylum_", "class_", "order_", "family_", "genus_")
# A formal species name: "Genus epithet", optionally a subspecies, or an epithet
# whose hyphen MIDORI2 turned into "_" (Nymphalis c_album). Anything else (sp.,
# cf., aff., nr., BOLD codes, hybrids with "x", lineage labels) is an interim name.
# The same rule as the development test of the COI marker (leave-out test).
FORMAL_SPECIES = re.compile(r"^[A-Z][a-z]+ [a-z]+(_[a-z]+)?( [a-z]+)?$")
# The tax part is split on the rank letters, not on commas: a few MIDORI2 species
# names contain a comma themselves.
FIELD = re.compile(r"(?:^|,)([kpcofgs]):(.*?)(?=,[kpcofgs]:|$)")


def stop(message):
    """Remove the half-written output, write the message to the log and stop the
    rule with it."""
    out_path.unlink(missing_ok=True)
    full_message = f"[midori2_to_metaflux_sintax] {message}"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(full_message + "\n")
    sys.exit(full_message)


def split_taxid(value):
    """'Arthropoda_6656' -> ('Arthropoda', '6656'); the taxid must be numeric."""
    name, _, taxid = value.rpartition("_")
    if not taxid.isdigit() or not name:
        return value, ""
    return name, taxid


# --- 1. NCBI tree of the MIDORI2 release: parent and rank of every taxid ----------
parent_of = {}
rank_of = {}
with nodes_path.open() as handle:
    for line in handle:
        fields = line.split("\t|\t")
        parent_of[fields[0]] = fields[1]
        rank_of[fields[0]] = fields[2]
kingdom_taxids = {taxid for taxid, rank in rank_of.items() if rank == "kingdom"}
kingdom_names = {}
with names_path.open() as handle:
    for line in handle:
        fields = line.rstrip("\t|\n").split("\t|\t")
        if fields[0] in kingdom_taxids and fields[3] == "scientific name":
            kingdom_names[fields[0]] = fields[1]

kingdom_cache = {}


def kingdom_of(taxid):
    """Name of the first ancestor of rank 'kingdom' (or the taxon itself), or ''."""
    if taxid in kingdom_cache:
        return kingdom_cache[taxid]
    path = []
    current = taxid
    found = ""
    while current in parent_of:
        if current in kingdom_cache:
            found = kingdom_cache[current]
            break
        path.append(current)
        if rank_of[current] == "kingdom":
            found = kingdom_names[current]
            break
        if parent_of[current] == current:      # the root (taxid 1) is its own parent
            break
        current = parent_of[current]
    else:
        stop(f"taxid {taxid} (or an ancestor) is not in {nodes_path}")
    for node in path:
        kingdom_cache[node] = found
    return found


# --- 2. Sequences to leave out: bacterial COI homologues under eukaryote names ----
# The list has one row per MIDORI2 record it examined, with the decision ("removed"
# or "listed") and the SHA-1 of the record's sequence. Every sequence of a removed
# row is left out, under whatever name MIDORI2 files it.
def sequence_sha1(sequence_lines):
    """SHA-1 of a record's sequence, the way the build of the bacterial set computes it:
    the sequence lines joined, line breaks removed, in upper case."""
    joined = "".join(line.rstrip("\r\n") for line in sequence_lines).upper()
    return hashlib.sha1(joined.encode("ascii")).hexdigest()


def file_sha256(path):
    """SHA-256 of a file's bytes, read in blocks of 1 MiB."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


removed_sequences = set()        # SHA-1 of every removed sequence
removed_records = set()          # the record IDs of the removed rows
list_made_from = set()           # midori2_file_sha256 of the rows (one value)
with removal_path.open() as handle:
    columns = handle.readline().rstrip("\n").split("\t")
    for needed in ("record_id", "decision", "sequence_sha1"):
        if needed not in columns:
            stop(f"{removal_path} has no column {needed!r}: it is a MIDORI2 list of the "
                 f"set of bacterial COI homologues before gtdb-r232.v2, which removed records by name; "
                 f"this MetaFlux removes them by sequence and needs the list of gtdb-r232.v2 "
                 f"or later")
    id_column = columns.index("record_id")
    decision_column = columns.index("decision")
    sha1_column = columns.index("sequence_sha1")
    made_from_column = columns.index("midori2_file_sha256") if "midori2_file_sha256" in columns else None
    for line in handle:
        values = line.rstrip("\n").split("\t")
        if made_from_column is not None and values[made_from_column]:
            list_made_from.add(values[made_from_column])
        if values[decision_column] == "removed":
            removed_sequences.add(values[sha1_column].lower())
            removed_records.add(values[id_column])

# --- 3. Rewrite the file ---------------------------------------------------------
# A record is its header line and the sequence lines after it. The sequence must be
# read before the record can be judged (its SHA-1), so each record is held until
# the next header, or the end of the file, and then written or left out.
counts = Counter()
seen_removed = set()
out_path.parent.mkdir(parents=True, exist_ok=True)


def write_record(fo, line_no, header, sequence_lines):
    """Leave out a record whose sequence the list removes; otherwise write it with
    its header converted and its sequence lines unchanged."""
    if ";tax=" not in header:
        stop(f"{in_path}, line {line_no}: no ';tax=' part: {header!r}")
    record_id, tax = header[1:].split(";tax=", 1)
    counts["records read"] += 1
    sha1 = sequence_sha1(sequence_lines)
    if sha1 in removed_sequences:
        seen_removed.add(sha1)
        counts["records left out: bacterial COI homologue under a eukaryote name"] += 1
        if record_id not in removed_records:
            counts["records left out: same sequence as a removed record, not named in the list"] += 1
        return
    fields = dict((letter, value) for letter, value in FIELD.findall(tax.rstrip(";")))
    domain, _ = split_taxid(fields.get("k", ""))
    if domain != "Eukaryota":
        stop(f"{in_path}, line {line_no}: the first field is {fields.get('k')!r}, "
             f"not k:Eukaryota_<taxid> as in MIDORI2: {header!r}")
    species_name, species_taxid = split_taxid(fields.get("s", ""))
    if not species_taxid:
        stop(f"{in_path}, line {line_no}: the species field has no taxid suffix: {header!r}")
    kingdom = kingdom_of(species_taxid)

    new_fields = ["d:Eukaryota"]
    if kingdom:
        new_fields.append(f"k:{kingdom}")
        counts[f"kingdom {kingdom}"] += 1
    else:
        counts["kingdom none (protists)"] += 1
    for letter in ("p", "c", "o", "f", "g"):
        name, _ = split_taxid(fields.get(letter, ""))
        if not name or name.startswith(PLACEHOLDER_PREFIXES):
            counts[f"placeholder dropped at {letter}:"] += 1
            continue
        new_fields.append(f"{letter}:{name}")
    if FORMAL_SPECIES.match(species_name):
        new_fields.append("s:" + species_name.replace(" ", "_"))
    else:
        counts["interim species name dropped"] += 1
    fo.write(f">{record_id};tax={','.join(new_fields)};\n")
    fo.writelines(sequence_lines)
    counts["records written"] += 1


with gzip.open(in_path, "rt") as fi, gzip.open(out_path, "wt", compresslevel=6) as fo:
    header = None            # the header of the record being read
    header_line_no = 0
    sequence_lines = []
    for line_no, line in enumerate(fi, 1):
        if line.startswith(">"):
            if header is not None:
                write_record(fo, header_line_no, header, sequence_lines)
            header = line.rstrip("\n")
            header_line_no = line_no
            sequence_lines = []
        else:
            sequence_lines.append(line)
    if header is not None:
        write_record(fo, header_line_no, header, sequence_lines)

# --- Removed sequences the file does not hold --------------------------------------
# A list made from this very file must find every removed sequence in it; anything
# else means the two do not belong together. A list made from another MIDORI2 file
# (another release) can name sequences this file lacks: they are counted, and the
# run goes on with the removals that apply.
missing = removed_sequences - seen_removed
input_sha256 = file_sha256(in_path)
list_matches_file = input_sha256 in list_made_from
counts["removed sequences in the list"] = len(removed_sequences)
counts["removed sequences of the list not in this MIDORI2 file"] = len(missing)
if missing and list_matches_file:
    stop(f"{len(missing)} removed sequences of {removal_path} are not in {in_path} "
         f"(for example SHA-1 {sorted(missing)[0]}), although the list was made from this "
         f"file (SHA-256 {input_sha256}): the list and the file do not match")

report_path.parent.mkdir(parents=True, exist_ok=True)
with report_path.open("w") as handle:
    handle.write("count\tvalue\n")
    for key in sorted(counts):
        handle.write(f"{key}\t{counts[key]}\n")

# --- 4. The log: the files and the main counts, for a human checking the reference
# after the run (every count is in the report above); nothing downstream reads it.
removed_count = counts["records left out: bacterial COI homologue under a eukaryote name"]
log_path.parent.mkdir(parents=True, exist_ok=True)
with log_path.open("w") as log:
    log.write(
        f"[midori2_to_metaflux_sintax] {counts['records read']} records read from {in_path}\n"
        f"[midori2_to_metaflux_sintax] {removed_count} left out as bacterial COI homologues "
        f"(list: {removal_path})\n"
        f"[midori2_to_metaflux_sintax] the list was made from "
        f"{'this MIDORI2 file' if list_matches_file else 'another MIDORI2 file'}; "
        f"{len(missing)} of its {len(removed_sequences)} removed sequences are not in it\n"
        f"[midori2_to_metaflux_sintax] {counts['records written']} records written to {out_path}\n"
        f"[midori2_to_metaflux_sintax] every count: {report_path}\n"
    )
