"""The reader of Kraken2's taxonomy file, taxo.k2d (00_common.smk, _kraken2_rank_sets).

What is tested: when a shotgun run starts, MetaFlux reads taxo.k2d from the Kraken2
database to learn which taxa the database ranks as a domain, a kingdom or a virus
realm, and whether it still ranks Viruses as a superkingdom (NCBI did until 2025).
The OTU-table step fills the d__ and k__ slots from this, and checks the report codes
against it. The tests read the tiny taxonomy written by tests/make_test_kraken2_db.py
(the same layout kraken2-build writes) and compare the reader's answer with the
tree the file was written from. A damaged or foreign file must stop the run with a
message that says what is wrong.
"""
import pytest

from helpers import COMMON_SMK, TEST_DIR, load_functions
from make_test_kraken2_db import (NODES, old_ncbi_nodes, rank_sets_from_nodes,
                                  taxo_k2d_bytes)

reader = load_functions(COMMON_SMK, ["_kraken2_rank_sets", "_text_at"])


def test_modern_ncbi_tree_gives_its_domains_kingdoms_and_realms(tmp_path):
    taxo = tmp_path / "taxo.k2d"
    taxo.write_bytes(taxo_k2d_bytes(NODES))
    rank_sets = reader._kraken2_rank_sets(taxo)
    assert rank_sets == rank_sets_from_nodes(NODES)
    # Spelled out, so the expectation is visible here too:
    assert rank_sets["domain"] == {"2": "Bacteria", "2157": "Archaea", "2759": "Eukaryota"}
    assert rank_sets["kingdom"] == {"3379134": "Pseudomonadati", "4751": "Fungi",
                                    "2731360": "Heunggongvirae"}
    assert rank_sets["realm"] == {"2731341": "Duplodnaviria"}
    assert rank_sets["virus_root"] == {}        # Viruses is an "acellular root" now


def test_pre_2025_tree_puts_viruses_aside_and_counts_superkingdoms_as_domains(tmp_path):
    taxo = tmp_path / "taxo.k2d"
    taxo.write_bytes(taxo_k2d_bytes(old_ncbi_nodes()))
    rank_sets = reader._kraken2_rank_sets(taxo)
    assert rank_sets == rank_sets_from_nodes(old_ncbi_nodes())
    assert rank_sets["domain"] == {"2": "Bacteria", "2157": "Archaea", "2759": "Eukaryota"}
    assert rank_sets["virus_root"] == {"10239": "Viruses"}
    assert "3379134" not in rank_sets["kingdom"]  # no prokaryote kingdoms back then
    assert rank_sets["realm"] == {"2731341": "Duplodnaviria"}


def test_missing_file_stops_with_the_database_hint(tmp_path):
    with pytest.raises(SystemExit) as stop:
        reader._kraken2_rank_sets(tmp_path / "taxo.k2d")
    assert "not found" in str(stop.value)
    assert "hash.k2d, opts.k2d and taxo.k2d" in str(stop.value)


def test_a_file_without_the_kraken2_tag_is_refused(tmp_path):
    taxo = tmp_path / "taxo.k2d"
    taxo.write_bytes(b"NOTK2TAX" + taxo_k2d_bytes(NODES)[8:])
    with pytest.raises(SystemExit) as stop:
        reader._kraken2_rank_sets(taxo)
    assert "is not a Kraken2 taxonomy file" in str(stop.value)


def test_a_truncated_file_is_refused_with_the_expected_size(tmp_path):
    # A copy that stopped short (or a future layout) has a size its header does not
    # explain; its numbers cannot be trusted, so the run stops.
    complete = taxo_k2d_bytes(NODES)
    taxo = tmp_path / "taxo.k2d"
    taxo.write_bytes(complete[:-10])
    with pytest.raises(SystemExit) as stop:
        reader._kraken2_rank_sets(taxo)
    assert "does not have the layout MetaFlux reads" in str(stop.value)
    assert f"{len(complete)} expected" in str(stop.value)


def test_text_at_reads_up_to_the_zero_byte():
    assert reader._text_at(b"abc\0def\0", 0) == "abc"
    assert reader._text_at(b"abc\0def\0", 4) == "def"


def test_committed_stub_database_is_what_the_writer_writes():
    # .test/stubs/kraken2_db/taxo.k2d is what the CI dry runs read. It must stay in
    # step with the tree in make_test_kraken2_db.py; regenerate it with
    #   python tests/make_test_kraken2_db.py .test/stubs/kraken2_db
    committed = TEST_DIR / "stubs" / "kraken2_db" / "taxo.k2d"
    assert committed.is_file(), f"{committed} is missing"
    assert committed.read_bytes() == taxo_k2d_bytes(NODES)
