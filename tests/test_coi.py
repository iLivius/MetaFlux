"""The COI marker: the MIDORI2 conversion and the measurements of the COI screens.

What is tested
  - workflow/scripts/refdb/midori2_to_metaflux_sintax.py, on a tiny MIDORI2-shaped
    file with a tiny NCBI taxonomy: the domain moves from k: to d:, the kingdom is
    added from the taxonomy (and left out for a protist), the taxid suffixes go,
    placeholder ranks and interim species names are left out of the header, the
    records of the removal list are dropped, and the report counts all of it. The
    script stops on a removal list made for another MIDORI2 release.
  - the per-ASV measurements of workflow/scripts/amplicon/80d_coi_screens.py: the
    open reading frames (TAA/TAG only), the frame and length conditions of the
    frame check (an offset of the ASV's start allowed, a frameshift inside not),
    the A/C/G/T identity (N columns left out, internal gaps counted) and the
    species name read from a SINTAX label.
  - the shipped COI pack: SINTAX letters in VSEARCH's order, rank lists of equal
    length, reference_values naming ranks of the model.
  - the start-up check of amplicon.coi.length_tolerance in 00_common.smk
    (_coi_length_tolerance): whole codons pass as the same int, anything else stops.
"""
import gzip

import pytest

from helpers import COMMON_SMK, SCRIPTS, fake_snakemake, load_functions, marker_pack, run_script, write_text

CONVERTER = SCRIPTS / "refdb" / "midori2_to_metaflux_sintax.py"
SCREENS = SCRIPTS / "amplicon" / "80d_coi_screens.py"

# A tiny NCBI tree in nodes.dmp format: root 1 > Eukaryota 2759 > Metazoa 33208 >
# Arthropoda 6656 > Insecta 50557 > Diptera 7147 > Muscidae 7367 > Musca 7369 >
# Musca domestica 7370; and a protist line without a kingdom: 2759 > Amoebozoa 554915
# > its species 1892891.
NODES = "".join(
    f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
    for taxid, parent, rank in [
        ("1", "1", "no rank"), ("2759", "1", "domain"), ("33208", "2759", "kingdom"),
        ("6656", "33208", "phylum"), ("50557", "6656", "class"), ("7147", "50557", "order"),
        ("7367", "7147", "family"), ("7369", "7367", "genus"), ("7370", "7369", "species"),
        ("554915", "2759", "phylum"), ("1892891", "554915", "species"),
    ])
NAMES = "33208\t|\tMetazoa\t|\t\t|\tscientific name\t|\n2759\t|\tEukaryota\t|\t\t|\tscientific name\t|\n"
MIDORI2 = (
    ">A1.1.<1.>658;tax=k:Eukaryota_2759,p:Arthropoda_6656,c:Insecta_50557,o:Diptera_7147,"
    "f:Muscidae_7367,g:Musca_7369,s:Musca domestica_7370;\nACGTACGT\n"
    ">A2.1.<1.>658;tax=k:Eukaryota_2759,p:Arthropoda_6656,c:Insecta_50557,o:Diptera_7147,"
    "f:Muscidae_7367,g:genus_Muscidae sp. BIOUG123_7370,s:Muscidae sp. BIOUG123_7370;\nACGTACGA\n"
    ">P1.1.<1.>890;tax=k:Eukaryota_2759,p:phylum_class_order_family_genus_Amoebozoa sp._1892891,"
    "c:class_order_family_genus_Amoebozoa sp._1892891,o:order_family_genus_Amoebozoa sp._1892891,"
    "f:family_genus_Amoebozoa sp._1892891,g:genus_Amoebozoa sp._1892891,s:Amoebozoa sp._1892891;\nTTTT\n"
    ">B1.1.<1.>658;tax=k:Eukaryota_2759,p:Arthropoda_6656,c:Insecta_50557,o:Diptera_7147,"
    "f:Muscidae_7367,g:Musca_7369,s:Musca domestica_7370;\nGGGG\n"
)


def convert(tmp_path, removal_text):
    fasta = tmp_path / "midori2.fasta.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(MIDORI2)
    out = tmp_path / "out" / "converted.fasta.gz"
    report = tmp_path / "out" / "report.tsv"
    sm = fake_snakemake(
        input={"fasta": str(fasta), "nodes": write_text(tmp_path / "nodes.dmp", NODES),
               "names": write_text(tmp_path / "names.dmp", NAMES),
               "removal": write_text(tmp_path / "removal.tsv", removal_text)},
        output={"fasta": str(out), "report": str(report)})
    run_script(CONVERTER, sm)
    with gzip.open(out, "rt") as handle:
        headers = [line.rstrip("\n") for line in handle if line.startswith(">")]
    counts = dict(line.rstrip("\n").split("\t") for line in report.read_text().splitlines()[1:])
    return headers, counts


def test_conversion(tmp_path):
    headers, counts = convert(tmp_path, "record_id\tdecision\nB1.1.<1.>658\tremoved\nA1.1.<1.>658\tlisted\n")
    assert headers == [
        ">A1.1.<1.>658;tax=d:Eukaryota,k:Metazoa,p:Arthropoda,c:Insecta,o:Diptera,f:Muscidae,g:Musca,s:Musca_domestica;",
        ">A2.1.<1.>658;tax=d:Eukaryota,k:Metazoa,p:Arthropoda,c:Insecta,o:Diptera,f:Muscidae;",
        ">P1.1.<1.>890;tax=d:Eukaryota;",
    ]
    assert counts["records read"] == "4"
    assert counts["records written"] == "3"
    assert counts["records left out: bacterial COX1 under a eukaryote name"] == "1"
    assert counts["kingdom Metazoa"] == "2"
    assert counts["kingdom none (protists)"] == "1"
    assert counts["placeholder dropped at g:"] == "2"
    assert counts["interim species name dropped"] == "2"


def test_removal_list_of_another_release_stops(tmp_path):
    with pytest.raises(SystemExit) as stop:
        convert(tmp_path, "record_id\tdecision\nZZ999.1.<1.>658\tremoved\n")
    assert "belongs to another MIDORI2 release" in str(stop.value)


def test_screen_measurements():
    screens = load_functions(SCREENS, ["STOPS", "ACGT", "open_frames", "acgt_identity", "species_of"])
    # frame 0 reads TAA at the start; frames 1 and 2 have no TAA/TAG; TGA is no stop here
    assert screens.open_frames("TAAGGTGACCC") == [1, 2]
    # an N column is left out of both counts, an internal gap counts as a difference
    identity, columns = screens.acgt_identity("ACGTN-A", "ACGTAGA")
    assert columns == 6
    assert identity == pytest.approx(100 * 5 / 6)
    assert screens.species_of("X1;tax=d:Eukaryota,k:Metazoa,g:Musca,s:Musca_domestica;") == "Musca_domestica"
    assert screens.species_of("X2;tax=d:Eukaryota;") == ""


def test_frame_and_length():
    screens = load_functions(SCREENS, ["frame_and_length"])
    # run frame 1, reference length 313 (Leray-XT)
    assert screens.frame_and_length(313, [1], 1, 313, 9) == (True, True)
    # one base short at the start: open in frame 0, the frame before the run's
    assert screens.frame_and_length(312, [0], 1, 313, 9) == (True, True)
    # one base short at the end: still open in the run's frame
    assert screens.frame_and_length(312, [1], 1, 313, 9) == (True, True)
    # one base long at the start: open in frame 2
    assert screens.frame_and_length(314, [2], 1, 313, 9) == (True, True)
    # a frameshift inside: stops in both frames
    assert screens.frame_and_length(312, [2], 1, 313, 9) == (False, True)
    assert screens.frame_and_length(312, [], 1, 313, 9) == (False, True)
    # whole codons too far away
    assert screens.frame_and_length(301, [1], 1, 313, 9) == (True, False)


def test_coi_pack_rank_model():
    pack = marker_pack("COI")
    letters = pack["rank_letters"]
    order = ["d", "k", "p", "c", "o", "f", "g", "s", "t"]
    positions = [order.index(letter) for letter in letters]
    assert positions == sorted(positions)
    assert len(pack["tax_levels"]) == len(letters) == len(pack["rank_prefixes"])
    assert set(pack["reference_values"]) <= set(pack["tax_levels"])
    assert pack["taxonomy_refdb"] is None and pack["taxonomy_sintax_db"] == "midori2_coi_sintax"


# ── The length_tolerance check in 00_common.smk ─────────────────────────────────────

tolerance = load_functions(COMMON_SMK, ["_coi_length_tolerance"])


def test_length_tolerance_accepts_whole_codons_as_the_same_int():
    # A valid value must come out as the int the rule always received, or Snakemake
    # would see changed params and run coi_screens again.
    for value, expected in [(0, 0), (9, 9), (12, 12), (9.0, 9), ("9", 9), (" 6 ", 6)]:
        result = tolerance._coi_length_tolerance(value)
        assert result == expected
        assert type(result) is int


@pytest.mark.parametrize("value", [-3, 10, 9.5, "9.5", "nine", "", True, False, None])
def test_length_tolerance_stops_on_anything_else(value):
    with pytest.raises(SystemExit) as stop:
        tolerance._coi_length_tolerance(value)
    message = str(stop.value)
    assert "amplicon.coi.length_tolerance must be a whole number of bases" in message
    assert message.endswith(f"got {value!r}")
