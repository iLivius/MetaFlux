"""The stop texts of rule amplicon_probe when no in-silico product is left to measure
(no_lengths_message in workflow/scripts/amplicon/40a_amplicon_probe.py).

What is tested, and why it matters:
  - The PCR-mode probe (16S, 18S, rpoB, COI) stops when no product is left for its length
    statistics, and each of the three reasons gets its own text:
      no product at all (no reference sequence had both primers);
      pass 3 left out every product (each holds a second forward-primer site);
      the products left all come from lineages the marker pack leaves out of the
      statistics (probe_exclude_taxa: for 16S Chloroplast, Eukaryota, Mitochondria),
      some perhaps after pass 3.
    The third case used to get the first text, which sent the user to check primers and
    orientation although the probe had found products.
  - The first two texts stay word for word as before.

The function is loaded on its own from the script (load_functions in helpers.py): the
probe itself runs cutadapt and seqtk, which the test environments do not have.
"""
from helpers import SCRIPTS, load_functions

probe = load_functions(SCRIPTS / "amplicon" / "40a_amplicon_probe.py", ["no_lengths_message"])

TAXA_16S = ["Chloroplast", "Eukaryota", "Mitochondria"]


def test_no_product_at_all_keeps_the_primer_text():
    text = probe.no_lengths_message(0, 0, 0, TAXA_16S)
    assert text == ("[amplicon_probe] ERROR: no reference sequences had both primers. "
                    "Check primer sequences against the reference DB orientation and "
                    "primer error tolerance.")


def test_every_product_with_a_second_site_keeps_its_text():
    text = probe.no_lengths_message(12, 12, 0, TAXA_16S)
    assert text == ("[amplicon_probe] ERROR: each of the 12 in-silico products holds a "
                    "second full-length forward-primer site, so none is left to measure "
                    "the amplicon length from. Set amplicon.expected_length to a number, "
                    "and length_filter.mode: manual with a range, for these primers.")


def test_products_all_from_left_out_lineages_say_so():
    text = probe.no_lengths_message(40, 0, 40, TAXA_16S)
    assert text.startswith("[amplicon_probe] ERROR: the primers gave 40 in-silico products, "
                           "but none is left to measure the amplicon length from: 40 come "
                           "from lineages the marker pack leaves out of the length statistics "
                           "(probe_exclude_taxa: Chloroplast, Eukaryota, Mitochondria).")
    assert "no reference sequences had both primers" not in text
    assert "second full-length" not in text


def test_left_out_lineages_after_pass_3_name_both_reasons():
    text = probe.no_lengths_message(40, 15, 25, TAXA_16S)
    assert "the primers gave 40 in-silico products" in text
    assert "25 come from lineages the marker pack leaves out" in text
    assert "and 15 hold a second full-length forward-primer site" in text
