#!/usr/bin/env python3
"""A tiny Kraken2 "database" for the tests and the CI dry runs: its taxonomy file
taxo.k2d, and Kraken2- or Bracken-style reports written against the same tree.

WHAT IT IS FOR
  When a shotgun run starts, MetaFlux reads one file of the Kraken2 database:
  taxo.k2d, Kraken2's own copy of the taxonomy. From it MetaFlux learns which taxa
  the database ranks as a domain, a kingdom or a virus realm (_kraken2_rank_sets in
  workflow/rules/shared/00_common.smk), and the OTU-table step later fills the
  lineage slots of every row from that knowledge (45a_finalize_otu_table.py).
  A real database is tens of gigabytes, so the tests use this stand-in instead:
  about forty taxa, written in exactly the byte layout kraken2-build uses, so the
  very same reader code runs on it as on a real database. Kraken2 itself never runs
  on it (there is no k-mer index in it), which is fine for the tests of the Python
  parts and for dry runs, where no job runs at all.

WHERE IT IS USED
  - tests/test_kraken2_taxonomy_file.py reads taxo.k2d back with MetaFlux's reader.
  - tests/test_otu_table_builder.py feeds Bracken-style reports of this tree to
    45a_finalize_otu_table.py and checks the table it writes.
  - tests/test_extract_taxon_lookup.py runs the taxon-name lookup of the read
    extraction on Kraken2-style reports of this tree.
  - tests/test_token_check.py checks filter tokens against this tree's domains,
    kingdoms and realms.
  - .test/stubs/kraken2_db/ is the committed copy that the CI dry-run configs point
    references.kraken_db at. Regenerate it with
        python tests/make_test_kraken2_db.py .test/stubs/kraken2_db
    A test checks that the committed copy is what this script writes, so the two
    cannot drift apart.

THE TREE
  A cut-down NCBI taxonomy as NCBI has it since 2025: the domains Bacteria, Archaea
  and Eukaryota ranked "domain", with kingdoms below them (Pseudomonadati, Fungi); one
  protist lineage that NCBI places in no kingdom (Plasmodium falciparum); Viruses as an
  "acellular root" with one realm (Duplodnaviria) and its viral kingdom
  (Heunggongvirae); and one virus NCBI places in no realm (Pandoravirus salinus, under
  the unranked "unclassified viruses"). The genus Plasmodium has a subgenus of the
  same name, as in NCBI: the taxon-name lookup of the read extraction has a rule for
  exactly this case.
  The taxids of the well-known taxa are NCBI's (Bacteria 2, Escherichia coli 562,
  Fungi 4751, Viruses 10239 ...). The two Peduovirus nodes carry made-up numbers.
  The tree is a fixture, not a reference: nothing may be looked up in it.
  old_ncbi_nodes() gives the same tree as NCBI had it before 2025: the three domains
  and Viruses ranked "superkingdom", and no kingdom above the bacterial phyla. Kraken2
  databases built from that taxonomy are still in use, and MetaFlux must read them.

THE FILE LAYOUT (Kraken2 source, taxonomy.cc, Taxonomy::WriteToDisk)
  Every number is an unsigned 8-byte integer, least significant byte first:
    - bytes 0-7: the text K2TAXDAT;
    - three numbers: the number of node records, the size of the name block, the
      size of the rank block;
    - one 56-byte record per node (seven numbers): parent record, first child
      record, number of children, where the name starts in the name block, where
      the rank name starts in the rank block, the taxid, and one unused field;
    - the name block, then the rank block: texts, each ended by a zero byte.
  Record 0 is an empty placeholder (taxid 0); the real nodes follow in breadth-first
  order from the root, as kraken2-build stores them.

Standard library only. Run as a script to write the stub database folder:
    python tests/make_test_kraken2_db.py <folder>
"""
import struct
import sys
from pathlib import Path

# ── The tree: one row per taxon, (taxid, parent taxid, rank, name) ────────────────
# The root is its own parent, as in NCBI's nodes.dmp. Children are listed after
# their parent; the order of sister taxa here is the order they take in a report.
NODES = [
    (1,       1,       "no rank",        "root"),
    (131567,  1,       "no rank",        "cellular organisms"),
    (2,       131567,  "domain",         "Bacteria"),
    (3379134, 2,       "kingdom",        "Pseudomonadati"),
    (1224,    3379134, "phylum",         "Pseudomonadota"),
    (1236,    1224,    "class",          "Gammaproteobacteria"),
    (91347,   1236,    "order",          "Enterobacterales"),
    (543,     91347,   "family",         "Enterobacteriaceae"),
    (561,     543,     "genus",          "Escherichia"),
    (562,     561,     "species",        "Escherichia coli"),
    (2157,    131567,  "domain",         "Archaea"),
    (2759,    131567,  "domain",         "Eukaryota"),
    (33154,   2759,    "clade",          "Opisthokonta"),
    (4751,    33154,   "kingdom",        "Fungi"),
    (4890,    4751,    "phylum",         "Ascomycota"),
    (4891,    4890,    "class",          "Saccharomycetes"),
    (4892,    4891,    "order",          "Saccharomycetales"),
    (4893,    4892,    "family",         "Saccharomycetaceae"),
    (4930,    4893,    "genus",          "Saccharomyces"),
    (4932,    4930,    "species",        "Saccharomyces cerevisiae"),
    (2698737, 2759,    "clade",          "Sar"),
    (33630,   2698737, "clade",          "Alveolata"),
    (5794,    33630,   "phylum",         "Apicomplexa"),
    (422676,  5794,    "class",          "Aconoidasida"),
    (5819,    422676,  "order",          "Haemosporida"),
    (1639119, 5819,    "family",         "Plasmodiidae"),
    (5820,    1639119, "genus",          "Plasmodium"),
    (418103,  5820,    "subgenus",       "Plasmodium"),
    (5833,    418103,  "species",        "Plasmodium falciparum"),
    (10239,   1,       "acellular root", "Viruses"),
    (2731341, 10239,   "realm",          "Duplodnaviria"),
    (2731360, 2731341, "kingdom",        "Heunggongvirae"),
    (2731618, 2731360, "phylum",         "Uroviricota"),
    (2731619, 2731618, "class",          "Caudoviricetes"),
    (900001,  2731619, "genus",          "Peduovirus"),            # made-up taxid
    (900002,  900001,  "species",        "Peduovirus P37"),        # made-up taxid
    (12429,   10239,   "no rank",        "unclassified viruses"),
    (2060084, 12429,   "genus",          "Pandoravirus"),
    (1349410, 2060084, "species",        "Pandoravirus salinus"),
]

# The ranks Kraken2 gives a one-letter report code. Every other rank (clade,
# subgenus, realm, acellular root, no rank) gets the code of the nearest ranked
# ancestor plus its distance from it: R1 for "cellular organisms", D1 for a clade
# under a domain, G1 for a subgenus, R2 for a realm (Viruses is R1).
RANK_CODE = {
    "domain": "D", "superkingdom": "D", "kingdom": "K", "phylum": "P", "class": "C",
    "order": "O", "family": "F", "genus": "G", "species": "S",
}


def old_ncbi_nodes():
    """The same tree as NCBI had it before 2025: Bacteria, Archaea, Eukaryota and
    Viruses ranked "superkingdom" (Kraken2 codes all four D), and no kingdom between
    Bacteria and its phyla. The viral kingdom (Heunggongvirae) and Fungi existed then
    already, so they stay."""
    old_nodes = []
    for taxid, parent, rank, name in NODES:
        if taxid == 3379134:
            continue                          # no prokaryote kingdoms before 2025
        if parent == 3379134:
            parent = 2                        # the phylum hangs from Bacteria directly
        if rank in ("domain", "acellular root"):
            rank = "superkingdom"
        old_nodes.append((taxid, parent, rank, name))
    return old_nodes


def rank_sets_from_nodes(nodes):
    """What MetaFlux's reader must return for this tree: four {taxid: name} dicts
    (taxids as text) of its domains, kingdoms and virus realms, and of Viruses when
    the tree ranks it superkingdom. Worked out from the rank column alone, so the
    tests can compare the reader's answer with an independent one."""
    rank_sets = {"domain": {}, "kingdom": {}, "realm": {}, "virus_root": {}}
    for taxid, parent, rank, name in nodes:
        if rank in ("domain", "superkingdom"):
            if taxid == 10239:
                rank_sets["virus_root"][str(taxid)] = name
            else:
                rank_sets["domain"][str(taxid)] = name
        elif rank == "kingdom":
            rank_sets["kingdom"][str(taxid)] = name
        elif rank == "realm":
            rank_sets["realm"][str(taxid)] = name
    return rank_sets


def children_by_parent(nodes):
    """{parent taxid: [child taxids]} in the order of the node table; the root is
    not its own child."""
    children = {}
    for taxid, parent, rank, name in nodes:
        if taxid == parent:
            continue
        children.setdefault(parent, []).append(taxid)
    return children


def taxo_k2d_bytes(nodes=NODES):
    """The bytes of a taxo.k2d file for this tree, in kraken2-build's layout (see the
    module docstring)."""
    rank_of = {}
    name_of = {}
    parent_of = {}
    for taxid, parent, rank, name in nodes:
        rank_of[taxid] = rank
        name_of[taxid] = name
        parent_of[taxid] = parent
    children = children_by_parent(nodes)

    # Breadth-first order from the root, as kraken2-build stores the nodes. Record 0
    # is the empty placeholder, so the first real node gets record number 1.
    order = [1]
    position = 0
    while position < len(order):
        order.extend(children.get(order[position], []))
        position += 1
    record_of = {}
    for index, taxid in enumerate(order):
        record_of[taxid] = index + 1

    # The two text blocks: every name once, in record order; every rank name once.
    name_block = b""
    name_offset = {}
    for taxid in order:
        name_offset[taxid] = len(name_block)
        name_block += name_of[taxid].encode() + b"\0"
    rank_block = b""
    rank_offset = {}
    for rank in sorted(set(rank_of.values())):
        rank_offset[rank] = len(rank_block)
        rank_block += rank.encode() + b"\0"

    # One 56-byte record per node: seven 8-byte numbers ("<7Q": little-endian).
    records = [struct.pack("<7Q", 0, 0, 0, 0, 0, 0, 0)]    # record 0, the placeholder
    for taxid in order:
        if taxid == 1:
            parent_record = 0
        else:
            parent_record = record_of[parent_of[taxid]]
        child_taxids = children.get(taxid, [])
        if child_taxids:
            first_child_record = record_of[child_taxids[0]]
        else:
            first_child_record = 0
        records.append(struct.pack("<7Q", parent_record, first_child_record,
                                   len(child_taxids), name_offset[taxid],
                                   rank_offset[rank_of[taxid]], taxid, 0))
    header = b"K2TAXDAT" + struct.pack("<QQQ", len(records), len(name_block), len(rank_block))
    return header + b"".join(records) + name_block + rank_block


def report_lines(nodes, direct_reads, minimizer_columns=False, unclassified=0):
    """The lines of a report over this tree, as Kraken2 (and Bracken, which keeps the
    layout) writes them: one row per taxon that has reads in its clade, top to bottom
    in tree order, the name indented by two spaces per level below the root.

    direct_reads       {taxid: reads placed on that taxon itself}; a clade's count is
                       the sum over the taxon and everything below it.
    minimizer_columns  False: the standard six columns (percentage, clade reads,
                       direct reads, rank code, taxid, name), the layout of a Bracken
                       report. True: the eight columns Kraken2 writes with
                       --report-minimizer-data, as MetaFlux runs it: two minimizer
                       counts (written as 0 here) between the direct reads and the
                       rank code, and an "unclassified" row (code U, taxid 0) first.
    unclassified       the read count of that first row (eight columns only).
    """
    rank_of = {}
    name_of = {}
    for taxid, parent, rank, name in nodes:
        rank_of[taxid] = rank
        name_of[taxid] = name
    children = children_by_parent(nodes)

    def clade_reads(taxid):
        total = direct_reads.get(taxid, 0)
        for child in children.get(taxid, []):
            total += clade_reads(child)
        return total

    total_reads = clade_reads(1) + unclassified
    lines = []

    def one_row(code, taxid, name, depth, clade, direct):
        percent = 100.0 * clade / total_reads if total_reads else 0.0
        fields = [f"{percent:.2f}", str(clade), str(direct)]
        if minimizer_columns:
            fields += ["0", "0"]
        fields += [code, str(taxid), "  " * depth + name]
        lines.append("\t".join(fields))

    def write_subtree(taxid, depth, code_letter, code_distance):
        # A ranked node gets its own letter; an unranked one the letter of the
        # ranked node above it plus its distance from that node.
        if rank_of[taxid] in RANK_CODE:
            code_letter = RANK_CODE[rank_of[taxid]]
            code_distance = 0
        else:
            code_distance += 1
        if code_distance == 0:
            code = code_letter
        else:
            code = f"{code_letter}{code_distance}"
        clade = clade_reads(taxid)
        if clade == 0:
            return                            # Kraken2 lists only taxa that got reads
        one_row(code, taxid, name_of[taxid], depth, clade, direct_reads.get(taxid, 0))
        for child in children.get(taxid, []):
            write_subtree(child, depth + 1, code_letter, code_distance)

    if minimizer_columns:
        one_row("U", 0, "unclassified", 0, unclassified, unclassified)
    # The root is written R at depth 0: Kraken2 starts at R with distance -1 and the
    # unranked root adds one.
    write_subtree(1, 0, "R", -1)
    return lines


def write_database_stub(folder, nodes=NODES):
    """Write the stub database folder: taxo.k2d for this tree, and empty hash.k2d
    and opts.k2d files, so that the folder has the three files a Kraken2 database
    has (MetaFlux's start-up message names all three). Only taxo.k2d is ever read."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "taxo.k2d").write_bytes(taxo_k2d_bytes(nodes))
    (folder / "hash.k2d").write_bytes(b"")
    (folder / "opts.k2d").write_bytes(b"")
    return folder


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python tests/make_test_kraken2_db.py <folder>")
    written = write_database_stub(sys.argv[1])
    print(f"wrote {written}/taxo.k2d ({len(NODES)} taxa), hash.k2d and opts.k2d (empty)")
