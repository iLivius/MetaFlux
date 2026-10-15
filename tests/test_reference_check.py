"""The taxonomy reference check, workflow/scripts/refdb/check_reference_ranks.py.

What is tested: before an amplicon run classifies anything, rule
check_taxonomy_reference reads every header of the reference file and checks that
the file has exactly the ranks the marker pack declares, in the right order, with
names that can stand in those ranks (the pack's reference_values). A classifier
never checks what a rank means: a file with one rank more or less, or with two
ranks swapped, puts every name in the wrong column without an error. The tests run
the script on tiny reference files written on the spot, with the rank models and
value rules of the shipped marker packs, and check that
  - the right files pass and a report is written;
  - the wrong ones stop with a message that names the problem: the file of another
    SILVA release, SILVA's own download with its "--other" kingdoms, a lineage
    shorter than the model without a closing ';', PR2's own letter order, a prefix on
    a rank that has none, a domain in UNITE's kingdom rank;
  - a gzipped file is recognised by its first bytes, whatever its name.
"""
import gzip

import pytest

from helpers import SCRIPTS, fake_snakemake, marker_pack, run_script, write_text

SCRIPT = SCRIPTS / "refdb" / "check_reference_ranks.py"


def run_check(tmp_path, reference_text, ref_format, model, prefix_style,
              species_from_addspecies=False, gzipped=False, label="the test pack"):
    """Write `reference_text` as the reference file (plain or gzipped), run the
    script on it with the given rank model, and return the path of the report it
    writes when every check passes. `model` is the pack dict (or release block) with
    tax_levels, rank_letters, rank_prefixes and reference_values."""
    reference = tmp_path / "reference.fasta"
    if gzipped:
        with gzip.open(reference, "wt") as handle:
            handle.write(reference_text)
    else:
        reference.write_text(reference_text)
    report = tmp_path / "stats" / "taxonomy_reference_check.json"
    sm = fake_snakemake(
        input={"reference": str(reference)},
        output={"report": str(report)},
        params={"ref_format": ref_format, "tax_levels": model["tax_levels"],
                "rank_letters": model["rank_letters"], "rank_prefixes": model["rank_prefixes"],
                "prefix_style": prefix_style,
                "species_from_addspecies": species_from_addspecies,
                "reference_values": model["reference_values"], "label": label},
        log=[str(tmp_path / "logs" / "check.log")],
    )
    run_script(SCRIPT, sm)
    return report


def fasta(*lineages):
    """A FASTA text with one record per lineage (the sequence does not matter)."""
    records = []
    for lineage in lineages:
        records.append(f">{lineage}\nACGT\n")
    return "".join(records)


SILVA_144 = marker_pack("16S")["releases"]["144"]
SILVA_138 = marker_pack("16S")["releases"]["138.2"]

# SILVA-style DADA2 headers: bare names, a ';' after the last one, lineages that
# stop above genus are normal.
SILVA_144_HEADERS = fasta(
    "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
    "Enterobacteriaceae;Escherichia;",
    "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;",
    "Archaea;Methanobacteriati;Euryarchaeota;Methanobacteria;",
    "Bacteria;Incertae Sedis;Chloroflexota;",
)
SILVA_138_HEADERS = fasta(
    "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;",
    "Archaea;Euryarchaeota;Methanobacteria;",
)


# ── DADA2 training sets (method rdp) ─────────────────────────────────────────────

def test_silva_144_training_set_passes_and_writes_the_report(tmp_path):
    # 16S fills the Species rank from addSpecies, so the file stops at genus: seven
    # ranks for an eight-rank model.
    report = run_check(tmp_path, SILVA_144_HEADERS, "dada2", SILVA_144, "bare",
                       species_from_addspecies=True)
    assert report.is_file()
    text = report.read_text()
    assert '"result": "passed"' in text
    assert '"records": 4' in text
    assert '"expected_depth": 7' in text


def test_gzipped_file_is_recognised_by_its_bytes_not_its_name(tmp_path):
    report = run_check(tmp_path, SILVA_144_HEADERS, "dada2", SILVA_144, "bare",
                       species_from_addspecies=True, gzipped=True)
    assert '"result": "passed"' in report.read_text()


def test_other_release_file_stops_on_its_depth(tmp_path):
    # The SILVA 138.2 training set (six ranks, no kingdom) under the 144 model:
    # every name would land one column too high.
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, SILVA_138_HEADERS, "dada2", SILVA_144, "bare",
                  species_from_addspecies=True)
    message = str(stop.value)
    assert "does not hold the ranks declared by the test pack" in message
    assert "the deepest lineage has 6 ranks, but the rank model expects 7" in message
    # The phylum sits in the Kingdom column and fails the kingdom rule too.
    assert "Kingdom holds a name its rule does not allow" in message
    assert not (tmp_path / "stats" / "taxonomy_reference_check.json").exists()


def test_silvas_own_download_stops_on_the_other_suffix(tmp_path):
    # SILVA's download writes "Pseudomonadati--other" in the kingdom rank; the
    # training set MetaFlux builds removes the suffix. The raw file pinned in its
    # place must stop.
    headers = SILVA_144_HEADERS.replace("Bacteria;Pseudomonadati;", "Bacteria;Pseudomonadati--other;")
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, headers, "dada2", SILVA_144, "bare", species_from_addspecies=True)
    assert "'Pseudomonadati--other'" in str(stop.value)
    assert "ends_with [ati]" in str(stop.value)


def test_short_lineage_without_closing_semicolon_stops(tmp_path):
    # assignTaxonomy pads a short lineage rank by rank only after a closing ';'.
    # Without it the padding is glued onto the last name and every lineage after
    # this one is shifted, so the check stops.
    headers = SILVA_144_HEADERS.replace("Archaea;Methanobacteriati;Euryarchaeota;Methanobacteria;",
                                        "Archaea;Methanobacteriati;Euryarchaeota;Methanobacteria")
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, headers, "dada2", SILVA_144, "bare", species_from_addspecies=True)
    assert "without a ';' at its end" in str(stop.value)


def test_silva_138_training_set_passes_under_its_own_model(tmp_path):
    report = run_check(tmp_path, SILVA_138_HEADERS, "dada2", SILVA_138, "bare",
                       species_from_addspecies=True)
    assert '"result": "passed"' in report.read_text()


def test_gyrb_gene_tag_must_carry_no_prefix(tmp_path):
    # The DD7RZ8 lineage starts with a bare gene tag (gyrB or other), then
    # prefixed GTDB ranks down to a species binomial.
    gyrb = marker_pack("gyrB")
    good = fasta("gyrB;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;"
                 "f__Enterobacteriaceae;g__Escherichia;s__Escherichia_coli",
                 "other;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;"
                 "s__Bacillus_subtilis")
    report = run_check(tmp_path, good, "dada2", gyrb, "embedded")
    assert '"result": "passed"' in report.read_text()
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, good.replace(">gyrB;", ">g__gyrB;"), "dada2", gyrb, "embedded")
    assert "carries a rank prefix, but this rank has none in the pack" in str(stop.value)
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, good.replace(">other;", ">parE;"), "dada2", gyrb, "embedded")
    assert "Gene holds a name its rule does not allow" in str(stop.value)


def test_unite_layout_reads_the_lineage_from_the_fifth_field(tmp_path):
    # assignTaxonomy treats a file as the UNITE general release when each of its
    # first ten headers holds "FU|refs" or "FU|reps"; the lineage is then the 5th
    # '|'-separated field, with UNITE's own k__/p__ prefixes.
    its = marker_pack("ITS")
    records = []
    for number in range(12):
        records.append(f"Fungus_{number}|AB{number:06d}|SH{number:07d}.09FU|refs|k__Fungi;"
                       f"p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;"
                       f"g__Fusarium;s__Fusarium_sp{number}")
    report = run_check(tmp_path, fasta(*records), "dada2", its, "embedded")
    text = report.read_text()
    assert '"result": "passed"' in text
    assert '"unite_layout": true' in text
    # A domain in UNITE's kingdom rank would mean UNITE added a rank above it.
    bad = fasta(*records).replace("k__Fungi", "k__Eukaryota", 1)
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, bad, "dada2", its, "embedded")
    assert "Kingdom holds a name its rule does not allow" in str(stop.value)
    assert "'Eukaryota'" in str(stop.value)


# ── SINTAX files (method sintax) ─────────────────────────────────────────────────

def sintax_model(pack):
    """The sintax rank model of a pack that declares one of its own (18S)."""
    return {"tax_levels": pack["sintax_tax_levels"], "rank_letters": pack["sintax_rank_letters"],
            "rank_prefixes": pack["sintax_rank_prefixes"],
            "reference_values": pack["reference_values"]}


PR2_DOMAIN_FIRST = fasta(
    "AB001;tax=d:Eukaryota,k:TSAR,p:Alveolata-Dinoflagellata,c:Dinophyceae,o:Gonyaulacales,"
    "f:Gonyaulacaceae,g:Alexandrium,s:Alexandrium_minutum;",
    "AB002;tax=d:Eukaryota,k:Obazoa,p:Opisthokonta-Fungi,c:Ascomycota,o:Saccharomycetales,"
    "f:Saccharomycetaceae,g:Saccharomyces,s:Saccharomyces_cerevisiae;",
    "AB003;tax=d:Eukaryota:plas,k:Archaeplastida,p:Chlorophyta-Chlorophyta_X,c:Mamiellophyceae,"
    "o:Mamiellales,f:Bathycoccaceae,g:Ostreococcus,s:Ostreococcus_tauri;",
)
# PR2's own UTAX file: the domain under k: and the supergroup under d:, which VSEARCH
# reads as a fixed top-down order (so the supergroup becomes the top rank).
PR2_RAW = PR2_DOMAIN_FIRST.replace("d:Eukaryota:plas,k:Archaeplastida", "k:Eukaryota:plas,d:Archaeplastida") \
                          .replace("d:Eukaryota,k:TSAR", "k:Eukaryota,d:TSAR") \
                          .replace("d:Eukaryota,k:Obazoa", "k:Eukaryota,d:Obazoa")


def test_pr2_domain_first_file_passes_under_the_18s_sintax_model(tmp_path):
    model = sintax_model(marker_pack("18S"))
    report = run_check(tmp_path, PR2_DOMAIN_FIRST, "sintax", model, "bare")
    text = report.read_text()
    assert '"result": "passed"' in text
    assert '"records": 3' in text


def test_pr2s_own_letter_order_stops(tmp_path):
    model = sintax_model(marker_pack("18S"))
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, PR2_RAW, "sintax", model, "bare")
    message = str(stop.value)
    assert "letters out of the rank model's top-down order" in message
    # With the letters swapped, the supergroup letter holds the domain name.
    assert "Supergroup (letter k:) holds a name its rule does not allow" in message
    assert "Domain (letter d:) holds a name its rule does not allow" in message


def test_silva_138_sintax_file_under_the_144_model_lacks_the_kingdom_letter(tmp_path):
    # Six letters where the model has seven (k: missing): the SINTAX parser would
    # drop nothing and say nothing, so the check must notice the letter never occurs.
    headers = fasta("AB001;tax=d:Bacteria,p:Pseudomonadota,c:Gammaproteobacteria,"
                    "o:Enterobacterales,f:Enterobacteriaceae,g:Escherichia;")
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, headers, "sintax", SILVA_144, "bare", species_from_addspecies=True)
    assert "the letter k: (Kingdom) never occurs" in str(stop.value)
    # Under its own model the same file passes (the species letter is allowed to be
    # absent, because addSpecies fills that rank).
    report = run_check(tmp_path, headers, "sintax", SILVA_138, "bare", species_from_addspecies=True)
    assert '"result": "passed"' in report.read_text()


def test_a_letter_outside_the_model_stops(tmp_path):
    headers = fasta("AB001;tax=d:Bacteria,k:Pseudomonadati,p:Pseudomonadota,c:Gammaproteobacteria,"
                    "o:Enterobacterales,f:Enterobacteriaceae,g:Escherichia;")
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, headers, "sintax", SILVA_138, "bare", species_from_addspecies=True)
    assert "a letter that is not one of this rank model's letters (d,p,c,o,f,g,s)" in str(stop.value)


# ── Rules the script cannot test stop it before it reads anything ────────────────

def test_a_rule_for_a_rank_outside_the_model_stops_the_script(tmp_path):
    model = dict(SILVA_138)
    model["reference_values"] = {"Domain": {"one_of": ["Bacteria"]}, "Kingdom": {"ends_with": ["ati"]}}
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, SILVA_138_HEADERS, "dada2", model, "bare", species_from_addspecies=True)
    assert "a rule for the rank 'Kingdom', which is not in the rank model" in str(stop.value)
    assert "Fix the marker pack" in str(stop.value)


def test_a_rule_for_the_addspecies_rank_stops_the_script(tmp_path):
    model = dict(SILVA_138)
    model["reference_values"] = {"Species": {"not_one_of": ["sp."]}}
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, SILVA_138_HEADERS, "dada2", model, "bare", species_from_addspecies=True)
    assert "which addSpecies fills from a file of its own" in str(stop.value)


def test_an_empty_file_stops(tmp_path):
    with pytest.raises(SystemExit) as stop:
        run_check(tmp_path, "", "dada2", SILVA_138, "bare", species_from_addspecies=True)
    assert "the file holds no FASTA header" in str(stop.value)
