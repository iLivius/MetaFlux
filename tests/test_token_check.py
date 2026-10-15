"""The start-up check of the keep/discard filter tokens (00_common.smk).

What is tested: a keep or discard token (such as d__Bacteria or k__Fungi) must be
able to match one whole slot of the taxonomy strings that the run will write.
MetaFlux checks that before any job starts, so that a token written for another
marker, method or MetaFlux version (k__Bacteria, the way the domain was written up
to v2.4.0) stops the run at once instead of silently removing every ASV or taxon at
the very end. The checks are pure functions in workflow/rules/shared/00_common.smk:
  _check_filter_tokens          the rules both modes share, against one table's slots
  _check_shotgun_filter_tokens  the same, plus what the Kraken2 database knows
                                (its domains, kingdoms and virus realms)
  _value_passes                 one name against one reference_values rule of a pack
  _as_taxon_list                how the config's keep/discard value becomes a list
The slot models come from the shipped marker packs, so a pack change shows up here.
"""
import pytest

from helpers import COMMON_SMK, load_functions, marker_pack
from make_test_kraken2_db import NODES, old_ncbi_nodes, rank_sets_from_nodes

# The functions under test, plus the constants and helpers they use.
checks = load_functions(COMMON_SMK, [
    "_as_taxon_list", "_value_passes", "_describe_value_rule", "_check_filter_tokens",
    "_check_shotgun_filter_tokens", "_DOMAIN_NAMES", "_TOKEN_PATTERN",
    "_KRAKEN2_SLOT_NAMES", "_KRAKEN2_SLOT_PREFIXES",
])


def silva_144_slots():
    """The rank names, prefixes and value rules of 16S on SILVA 144, from the pack:
    Domain, Kingdom, Phylum ... Species, with a rule for Domain and Kingdom."""
    release = marker_pack("16S")["releases"]["144"]
    return release["tax_levels"], release["rank_prefixes"], release["reference_values"]


# ── _as_taxon_list: the config value becomes a list ──────────────────────────────

def test_as_taxon_list_wraps_a_bare_token_and_empties_null():
    # YAML lets `keep: k__Fungi` through as plain text; it must become a one-item
    # list, never be split into letters. null and [] both mean "no token".
    assert checks._as_taxon_list("k__Fungi") == ["k__Fungi"]
    assert checks._as_taxon_list(None) == []
    assert checks._as_taxon_list([]) == []
    assert checks._as_taxon_list(["d__Bacteria", "d__Archaea"]) == ["d__Bacteria", "d__Archaea"]


# ── _value_passes: a name against a reference_values rule ────────────────────────

def test_value_passes_follows_the_three_kinds_of_rule():
    kingdom_rule = {"ends_with": ["ati"], "one_of": ["Incertae Sedis"]}
    assert checks._value_passes("Bacillati", kingdom_rule)            # allowed ending
    assert checks._value_passes("Incertae Sedis", kingdom_rule)       # allowed name
    assert not checks._value_passes("Bacillota", kingdom_rule)        # a phylum
    assert not checks._value_passes("Pseudomonadati--other", kingdom_rule)
    # not_one_of forbids names; without one_of/ends_with anything else passes.
    forbid_domains = {"not_one_of": ["Bacteria", "Archaea", "Eukaryota"]}
    assert checks._value_passes("Fungi", forbid_domains)
    assert not checks._value_passes("Eukaryota", forbid_domains)
    # A rule written out for a message lists its parts in a fixed order.
    assert checks._describe_value_rule(kingdom_rule) == "one_of [Incertae Sedis]; ends_with [ati]"


# ── _check_filter_tokens on 16S, SILVA 144 ───────────────────────────────────────

def test_good_16s_tokens_give_no_problem():
    names, prefixes, rules = silva_144_slots()
    tokens = ["d__Bacteria", "d__Archaea", "k__Pseudomonadati", "k__Incertae Sedis",
              "o__Chloroplast", "f__Mitochondria", "g__Ralstonia"]
    assert checks._check_filter_tokens(tokens, names, prefixes, rules, "keep") == []


def test_domain_under_k_names_the_token_to_write():
    # The v2.4.0 way of writing the domain. The message must name the replacement.
    names, prefixes, rules = silva_144_slots()
    problems = checks._check_filter_tokens(["k__Bacteria"], names, prefixes, rules,
                                           "amplicon.taxonomy.filter.keep")
    assert len(problems) == 1
    assert "puts a domain under k__" in problems[0]
    assert "write d__Bacteria" in problems[0]
    assert problems[0].startswith("amplicon.taxonomy.filter.keep:")


def test_a_phylum_in_the_kingdom_slot_fails_the_pack_rule():
    names, prefixes, rules = silva_144_slots()
    problems = checks._check_filter_tokens(["k__Bacillota"], names, prefixes, rules, "keep")
    assert len(problems) == 1
    assert "is not a Kingdom of this reference" in problems[0]
    assert "ends_with [ati]" in problems[0]


def test_other_suffix_is_read_as_the_plain_name():
    # SILVA 144 writes X--other for a clade next to X; the filter reads it as X.
    names, prefixes, rules = silva_144_slots()
    assert checks._check_filter_tokens(["k__Pseudomonadati--other"], names, prefixes,
                                       rules, "keep") == []


def test_unknown_prefix_missing_name_and_no_prefix_are_problems():
    names, prefixes, rules = silva_144_slots()
    problems = checks._check_filter_tokens(["x__Bacteria", "g__", "Bacteria"],
                                           names, prefixes, rules, "discard")
    assert len(problems) == 3
    assert "uses x__, which this table does not write" in problems[0]
    assert "has no name after its prefix" in problems[1]
    assert "has no rank prefix" in problems[2]
    assert "Bacteria is a domain: write d__Bacteria" in problems[2]


def test_a_value_yaml_did_not_read_as_text_is_a_problem():
    # `keep: [yes]` arrives as the value True; `keep: [""]` as an empty string.
    names, prefixes, rules = silva_144_slots()
    problems = checks._check_filter_tokens([True, ""], names, prefixes, rules, "keep")
    assert len(problems) == 2
    assert "is not text" in problems[0]
    assert "an empty token" in problems[1]


def test_prefix_note_explains_a_species_token_under_sintax():
    # Under method sintax the 16S species slot is left out (SILVA's SINTAX file
    # stops at genus); the run passes a note for the s__ prefix, quoted in the problem.
    names, prefixes, rules = silva_144_slots()
    without_species = names[:-1], prefixes[:-1]
    note = {"s__": "16S under method: sintax writes no species."}
    problems = checks._check_filter_tokens(["s__Escherichia_coli"], *without_species,
                                           rules, "keep", note)
    assert len(problems) == 1
    assert "uses s__, which this table does not write" in problems[0]
    assert "writes no species" in problems[0]


# ── Other markers: ITS (kingdom first), gyrB (a slot without prefix), 18S ─────────

def test_its_has_no_domain_slot_so_a_domain_token_is_explained():
    pack = marker_pack("ITS")
    problems = checks._check_filter_tokens(["d__Bacteria"], pack["tax_levels"],
                                           pack["rank_prefixes"], pack["reference_values"],
                                           "keep")
    assert len(problems) == 1
    assert "uses d__, which this table does not write" in problems[0]
    assert "its slots: Kingdom (k__), Phylum (p__)" in problems[0]
    # UNITE's kingdom slot holds Fungi; a domain name there is refused, and since
    # this table has no d__ slot the hint says so instead of naming a d__ token.
    assert checks._check_filter_tokens(["k__Fungi"], pack["tax_levels"],
                                       pack["rank_prefixes"], pack["reference_values"],
                                       "keep") == []
    problems = checks._check_filter_tokens(["k__Eukaryota"], pack["tax_levels"],
                                           pack["rank_prefixes"], pack["reference_values"],
                                           "keep")
    assert len(problems) == 1
    assert "puts a domain under k__" in problems[0]
    assert "this table has no domain rank: its first slot is Kingdom (k__)" in problems[0]


def test_gyrb_gene_tag_is_the_slot_without_a_prefix():
    # The DD7RZ8 reference starts every lineage with a bare gene tag, gyrB or
    # other (a co-amplified paralog). A bare token stands for that slot.
    pack = marker_pack("gyrB")
    model = (pack["tax_levels"], pack["rank_prefixes"], pack["reference_values"])
    assert checks._check_filter_tokens(["other"], *model, "discard") == []
    problems = checks._check_filter_tokens(["parE", "d__Bacteria"], *model, "discard")
    assert len(problems) == 2
    assert "parE is not a Gene of this reference" in problems[0]
    assert "one_of [gyrB, other]" in problems[0]
    assert "uses d__, which this table does not write" in problems[1]


def test_18s_sintax_refuses_a_domain_in_the_supergroup_slot():
    # PR2's own SINTAX file had the domain under the supergroup letter; a token
    # written from such a table (sg__Eukaryota) can never match the corrected one.
    pack = marker_pack("18S")
    model = (pack["sintax_tax_levels"], pack["sintax_rank_prefixes"], pack["reference_values"])
    assert checks._check_filter_tokens(["d__Eukaryota", "sg__TSAR"], *model, "keep") == []
    problems = checks._check_filter_tokens(["sg__Eukaryota"], *model, "keep")
    assert len(problems) == 1
    assert "Eukaryota is not a Supergroup of this reference" in problems[0]
    assert "write d__Eukaryota" in problems[0]


# ── _check_shotgun_filter_tokens: the eight Kraken2 slots and the database ───────

def test_shotgun_tokens_that_fit_the_database_pass_quietly():
    rank_sets = rank_sets_from_nodes(NODES)
    problems, warnings = checks._check_shotgun_filter_tokens(
        ["d__Bacteria", "d__Archaea", "k__Fungi", "d__Duplodnaviria"], ["g__Homo"],
        rank_sets, True)
    assert problems == []
    assert warnings == []


def test_shotgun_domain_under_k_is_a_problem():
    problems, warnings = checks._check_shotgun_filter_tokens(
        ["k__Bacteria"], [], rank_sets_from_nodes(NODES), True)
    assert len(problems) == 1
    assert "write d__Bacteria" in problems[0]
    assert problems[0].startswith("shotgun.taxonomy_filter.keep:")


def test_viruses_is_never_a_slot_value_and_the_realms_are_listed():
    # A virus carries its realm in d__; Viruses itself is written nowhere.
    problems, warnings = checks._check_shotgun_filter_tokens(
        [], ["d__Viruses"], rank_sets_from_nodes(NODES), False)
    assert len(problems) == 1
    assert "no slot of otu_table.tsv holds Viruses" in problems[0]
    assert "d__Duplodnaviria" in problems[0]
    # The same on a database built before 2025, where Viruses was a superkingdom.
    problems, warnings = checks._check_shotgun_filter_tokens(
        [], ["d__Viruses"], rank_sets_from_nodes(old_ncbi_nodes()), False)
    assert len(problems) == 1
    assert "no slot of otu_table.tsv holds Viruses" in problems[0]


def test_a_realm_under_k_names_the_d_token():
    problems, warnings = checks._check_shotgun_filter_tokens(
        ["k__Duplodnaviria"], [], rank_sets_from_nodes(NODES), True)
    assert len(problems) == 1
    assert "is a virus realm" in problems[0]
    assert "write d__Duplodnaviria" in problems[0]


def test_names_this_database_lacks_only_warn_unless_keep_cannot_match_at_all():
    rank_sets = rank_sets_from_nodes(NODES)
    # A realm name of another NCBI release: well formed, so only a warning in discard.
    problems, warnings = checks._check_shotgun_filter_tokens(
        [], ["d__Monodnaviria"], rank_sets, True)
    assert problems == []
    assert len(warnings) == 1
    assert "neither a domain nor a virus realm of this Kraken2 database" in warnings[0]
    # The same name as the only keep token, with the filter on: the table would be
    # empty, which is a problem.
    problems, warnings = checks._check_shotgun_filter_tokens(
        ["d__Monodnaviria"], [], rank_sets, True)
    assert len(warnings) == 1
    assert len(problems) == 1
    assert "otu_table.tsv would be empty" in problems[0]
    # With the filter switched off the same keep list only warns.
    problems, warnings = checks._check_shotgun_filter_tokens(
        ["d__Monodnaviria"], [], rank_sets, False)
    assert problems == []
    assert len(warnings) == 1


def test_kingdom_tokens_warn_on_a_database_without_kingdoms():
    # GTDB has no kingdom rank, so k__ is empty in every row of its tables.
    no_kingdoms = rank_sets_from_nodes(NODES)
    no_kingdoms["kingdom"] = {}
    problems, warnings = checks._check_shotgun_filter_tokens(
        [], ["k__Fungi"], no_kingdoms, True)
    assert problems == []
    assert len(warnings) == 1
    assert "has no kingdom rank" in warnings[0]
    # On a database with kingdoms, a name that is not one of them warns too.
    problems, warnings = checks._check_shotgun_filter_tokens(
        [], ["k__Metazoa"], rank_sets_from_nodes(NODES), True)
    assert problems == []
    assert "is not a kingdom of this Kraken2 database" in warnings[0]
