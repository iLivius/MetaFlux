"""The primer note of shotgun runs (_primer_fields_set_by_user in 00_common.smk).

What is tested, and why it matters:
  - A shotgun run prints "mode=shotgun but amplicon.primers fields are set — they will
    be ignored" only when amplicon.primers.fwd or .rev names a primer file of the
    user's own: the note is for a config that may have been meant as an amplicon run.
  - The template config/config.yaml is loaded under every run config, and its two
    primer fields hold the placeholders /path/to/fwd_primer.fasta and
    /path/to/rev_primer.fasta. Up to v2.6.1 those counted as set, so every shotgun run
    whose config did not mention primers printed the note. They no longer count, nor
    does an empty or null field, a missing primers block or a missing amplicon block.

The function is read from 00_common.smk and run on its own (load_functions in
helpers.py), on the amplicon block as the merged config holds it.
"""
from helpers import COMMON_SMK, load_functions

checks = load_functions(COMMON_SMK, ["_primer_fields_set_by_user"])

TEMPLATE_PRIMERS = {"fwd": "/path/to/fwd_primer.fasta", "rev": "/path/to/rev_primer.fasta",
                    "orientation": "fixed"}


def test_the_template_placeholders_get_no_note():
    """A run config that never mentions primers: the template's placeholders only."""
    assert checks._primer_fields_set_by_user({"primers": TEMPLATE_PRIMERS}) == []


def test_empty_null_or_missing_fields_get_no_note():
    assert checks._primer_fields_set_by_user({"primers": {"fwd": None, "rev": ""}}) == []
    assert checks._primer_fields_set_by_user({"primers": None}) == []
    assert checks._primer_fields_set_by_user({"type": "16S"}) == []
    assert checks._primer_fields_set_by_user(None) == []


def test_a_primer_file_of_the_users_own_gets_the_note():
    """One real path is enough; the placeholder left in the other field does not count."""
    one_real = dict(TEMPLATE_PRIMERS, fwd="/data/primers/515F.fasta")
    assert checks._primer_fields_set_by_user({"primers": one_real}) == ["fwd"]
    both_real = {"fwd": "primers/515F.fasta", "rev": "primers/806R.fasta"}
    assert checks._primer_fields_set_by_user({"primers": both_real}) == ["fwd", "rev"]
