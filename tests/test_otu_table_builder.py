"""The shotgun OTU-table builder, workflow/scripts/shotgun/45a_finalize_otu_table.py.

What is tested: rule finalize_otu_table reads one Bracken report per sample and
writes otu_table.tsv and otu_table.biom: one row per taxon at the Bracken rank, one
column per sample, and a lineage of eight slots d__ k__ p__ c__ o__ f__ g__ s__
taken from the reports' own rank codes. The tests run the whole script on reports
written over the tiny taxonomy of tests/make_test_kraken2_db.py and check:
  - the lineages: kingdom filled where the database has one (Pseudomonadati, Fungi),
    empty where it has none (a protist), the realm in d__ for a virus, empty d__ and
    k__ for a virus NCBI places in no realm, whole species names;
  - the counts: Bracken's clade reads, 0 for a sample without the taxon, rows in the
    order of first appearance; the BIOM file holds the same rows;
  - the stops: a report whose codes disagree with the database (an old Kraken2 wrote
    the domains as R1), a keep/discard filter that removes every taxon.
The script needs biom-format and numpy (workflow/envs/biom.yaml); without them the
tests here are skipped, not failed.
"""
import pytest

biom = pytest.importorskip("biom")

from helpers import SCRIPTS, fake_snakemake, run_script, write_text
from make_test_kraken2_db import NODES, old_ncbi_nodes, rank_sets_from_nodes, report_lines

SCRIPT = SCRIPTS / "shotgun" / "45a_finalize_otu_table.py"

# Reads placed directly on the species of the two test samples. Sample B has no
# yeast and no viruses, so those rows must read 0 in its column.
READS_A = {562: 400, 4932: 100, 5833: 50, 900002: 300, 1349410: 200}
READS_B = {562: 40, 5833: 5}

LINEAGE_ECOLI = ("d__Bacteria;k__Pseudomonadati;p__Pseudomonadota;c__Gammaproteobacteria;"
                 "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli")
LINEAGE_YEAST = ("d__Eukaryota;k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;"
                 "f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces cerevisiae")
LINEAGE_PLASMODIUM = ("d__Eukaryota;k__;p__Apicomplexa;c__Aconoidasida;o__Haemosporida;"
                      "f__Plasmodiidae;g__Plasmodium;s__Plasmodium falciparum")
LINEAGE_PHAGE = ("d__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;c__Caudoviricetes;"
                 "o__;f__;g__Peduovirus;s__Peduovirus P37")
LINEAGE_PANDORAVIRUS = "d__;k__;p__;c__;o__;f__;g__Pandoravirus;s__Pandoravirus salinus"


def build_table(tmp_path, nodes=NODES, reads_by_sample=None, filter_enabled=False,
                keep=(), discard=(), edit_report=None):
    """Write one Bracken report per sample over `nodes`, run the script on them with
    the rank sets of the same tree, and return the lines of otu_table.tsv. `edit_report`
    may change the text of sample A's report before the script reads it."""
    if reads_by_sample is None:
        reads_by_sample = {"A": READS_A, "B": READS_B}
    report_paths = []
    for sample, reads in reads_by_sample.items():
        text = "\n".join(report_lines(nodes, reads)) + "\n"
        if edit_report is not None and sample == "A":
            text = edit_report(text)
        report_paths.append(write_text(tmp_path / f"{sample}_report.txt", text))
    rank_sets = rank_sets_from_nodes(nodes)
    sm = fake_snakemake(
        input={"reports": report_paths},
        output={"biom": str(tmp_path / "otu_table.biom"), "tsv": str(tmp_path / "otu_table.tsv")},
        params={"samples": list(reads_by_sample), "tax_lev": "S",
                "domains": rank_sets["domain"], "kingdoms": rank_sets["kingdom"],
                "realms": rank_sets["realm"], "virus_root": rank_sets["virus_root"],
                "filter_enabled": filter_enabled, "filter_keep": list(keep),
                "filter_discard": list(discard)},
        log=[str(tmp_path / "finalize_otu_table.log")],
    )
    run_script(SCRIPT, sm)
    return (tmp_path / "otu_table.tsv").read_text().splitlines()


def test_lineages_counts_and_row_order(tmp_path):
    lines = build_table(tmp_path)
    assert lines[0] == "# Constructed from biom file"
    assert lines[1] == "#OTU ID\tA\tB\ttaxonomy"
    # Rows in the order the taxa first appear, reading sample A's report top to
    # bottom: bacteria, then the eukaryotes, then the viruses. Counts are written as
    # floating-point numbers, as the tables always were.
    assert lines[2:] == [
        f"562\t400.0\t40.0\t{LINEAGE_ECOLI}",
        f"4932\t100.0\t0.0\t{LINEAGE_YEAST}",
        f"5833\t50.0\t5.0\t{LINEAGE_PLASMODIUM}",
        f"900002\t300.0\t0.0\t{LINEAGE_PHAGE}",
        f"1349410\t200.0\t0.0\t{LINEAGE_PANDORAVIRUS}",
    ]


def test_biom_file_holds_the_same_rows_and_lineages(tmp_path):
    build_table(tmp_path)
    table = biom.load_table(str(tmp_path / "otu_table.biom"))
    assert list(table.ids(axis="observation")) == ["562", "4932", "5833", "900002", "1349410"]
    assert list(table.ids(axis="sample")) == ["A", "B"]
    first_row = table.metadata(axis="observation")[0]
    assert list(first_row["taxonomy"]) == LINEAGE_ECOLI.split(";")
    assert table.get_value_by_ids("562", "B") == 40.0
    assert table.get_value_by_ids("4932", "B") == 0.0


def test_pre_2025_database_still_gives_the_virus_its_realm(tmp_path):
    # On a database built from NCBI's older taxonomy Kraken2 codes Viruses D (it was
    # a superkingdom). That D must be skipped: the phage gets its realm in d__, never
    # d__Viruses, and the run must not stop on the code.
    lines = build_table(tmp_path, nodes=old_ncbi_nodes())
    rows = {line.split("\t")[0]: line.split("\t")[-1] for line in lines[2:]}
    assert rows["900002"] == LINEAGE_PHAGE
    assert rows["1349410"] == LINEAGE_PANDORAVIRUS
    # Without a kingdom rank above the bacterial phyla, k__ stays empty.
    assert rows["562"] == LINEAGE_ECOLI.replace("k__Pseudomonadati", "k__")


def test_keep_filter_keeps_only_matching_rows(tmp_path):
    lines = build_table(tmp_path, filter_enabled=True, keep=["d__Bacteria", "k__Fungi"])
    assert [line.split("\t")[0] for line in lines[2:]] == ["562", "4932"]
    lines = build_table(tmp_path, filter_enabled=True, discard=["g__Plasmodium"])
    assert "5833" not in [line.split("\t")[0] for line in lines[2:]]
    assert len(lines[2:]) == 4


def test_filter_is_ignored_when_switched_off(tmp_path):
    lines = build_table(tmp_path, filter_enabled=False, keep=["d__Bacteria"])
    assert len(lines[2:]) == 5


def test_a_filter_that_removes_every_taxon_stops_the_job(tmp_path):
    # k__Bacteria matches no slot (the domain is d__Bacteria): the table would be
    # empty, so the job stops and lists the d__ values the table does hold.
    with pytest.raises(SystemExit) as stop:
        build_table(tmp_path, filter_enabled=True, keep=["k__Bacteria"])
    message = str(stop.value)
    assert "removed every taxon of the table" in message
    assert "keep: [k__Bacteria]" in message
    assert "d__Bacteria (1 taxon)" in message
    assert "an empty d__ (1 taxon)" in message
    assert not (tmp_path / "otu_table.tsv").exists()


def test_a_report_of_an_old_kraken2_stops_on_the_domain_code(tmp_path):
    # Kraken2 before 2.17 wrote NCBI's 2025 domains as R1 instead of D. Such a
    # report would leave d__ empty everywhere, so the job stops and says why.
    def recode_bacteria(text):
        return text.replace("\tD\t2\t", "\tR1\t2\t")
    with pytest.raises(SystemExit) as stop:
        build_table(tmp_path, edit_report=recode_bacteria)
    message = str(stop.value)
    assert "Bacteria (taxid 2) is a domain in the Kraken2 database" in message
    assert "codes it R1, not D" in message
    assert "Kraken2 releases before 2.17" in message


def test_a_report_coding_a_non_domain_as_d_stops(tmp_path):
    # A report made with another database may code other taxa D.
    def recode_fungi(text):
        return text.replace("\tK\t4751\t", "\tD\t4751\t")
    with pytest.raises(SystemExit) as stop:
        build_table(tmp_path, edit_report=recode_fungi)
    assert "Fungi (taxid 4751) is coded D, but it is not a domain" in str(stop.value)


def test_a_table_with_no_taxon_at_the_rank_is_written_empty(tmp_path):
    # Reads only above the species rank: no row, but a table all the same.
    lines = build_table(tmp_path, reads_by_sample={"A": {561: 10}, "B": {561: 3}})
    assert lines == ["# Constructed from biom file", "#OTU ID\tA\tB\ttaxonomy"]
    log = (tmp_path / "finalize_otu_table.log").read_text()
    assert "WARNING: no sample has a taxon at rank S" in log
