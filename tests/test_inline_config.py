"""Settings typed inline with --config are read as in a config file.

What is tested
  - _as_in_a_config_file in 00_common.smk, on what Snakemake hands over for a block
    typed after --config, such as 'amplicon={phylogeny: {enabled: false}}': every
    value of such a block is text. true/false become booleans, whole numbers ints,
    decimals floats and null None, in blocks and lists at any depth; any other text (a
    marker name, a path, a primer, and 1e8, which YAML reads as text in a config file
    too) is kept as typed; a value that is not text (a top-level number Snakemake has
    already converted) is kept.
  - The result is what the same block gives when written in a config file. The block
    is read here with YAML's base loader, the loader Snakemake uses for --config, so
    the test starts from the very text that caused the bug ("false" switched a
    setting on) without needing Snakemake.

What is not tested here
  - that the result reaches the rules: the dry-run matrix (.test/run_dry_runs.sh)
    checks that end to end, with the phylogeny switch typed inline.
"""
import yaml

from helpers import COMMON_SMK, load_functions

inline = load_functions(COMMON_SMK, ["_as_in_a_config_file"])


def test_on_off_settings_become_booleans():
    block = {"phylogeny": {"enabled": "false", "support": "true"}, "coi": {"frame_check": "False"}}
    assert inline._as_in_a_config_file(block) == {
        "phylogeny": {"enabled": False, "support": True}, "coi": {"frame_check": False}}


def test_numbers_become_numbers_in_blocks_and_lists():
    block = {"dada2": {"max_ee": ["3", "2.5"], "learn_errors": {"max_consist": "999"}}}
    assert inline._as_in_a_config_file(block) == {
        "dada2": {"max_ee": [3, 2.5], "learn_errors": {"max_consist": 999}}}


def test_other_text_is_kept_as_typed():
    block = {"type": "ITS", "fastq_dir": "/data/run 1/fastq", "forward": "GTGYCAGCMGCCGCGGTAA",
             "nbases": "1e8", "note": "a: b"}
    assert inline._as_in_a_config_file(block) == block


def test_null_and_values_that_are_not_text():
    assert inline._as_in_a_config_file({"x": "null", "y": "~", "z": ""}) == {"x": None, "y": None, "z": None}
    assert inline._as_in_a_config_file({"mode": "shotgun", "threads": 8}) == {"mode": "shotgun", "threads": 8}


def test_same_result_as_the_block_written_in_a_config_file():
    text = "{phylogeny: {enabled: false}, dada2: {max_ee: [3, 3], learn_errors: {nbases: 1e8}}, type: COI}"
    from_a_config_file = yaml.safe_load(text)
    as_snakemake_hands_it_over = yaml.load(text, Loader=yaml.BaseLoader)
    # The text behind the bug: bool("false") is True.
    assert as_snakemake_hands_it_over["phylogeny"]["enabled"] == "false"
    assert as_snakemake_hands_it_over["dada2"]["max_ee"] == ["3", "3"]
    assert inline._as_in_a_config_file(as_snakemake_hands_it_over) == from_a_config_file
