"""A key written twice in the same block of a run config stops the run at start-up
(_check_no_repeated_keys in workflow/rules/shared/00_common.smk).

What is tested, and why it matters:
  - YAML keeps only the last of two keys written in one block, without a word, and
    Snakemake reads config files that way (PyYAML's FullLoader). Several messages of
    MetaFlux ask for one line under amplicon: -> dada2: -> learn_errors:. A user who
    pastes it as a new amplicon: block at the end of the run config would lose every
    other amplicon setting of the file (the marker, the primers, the taxonomy), which
    would fall back to the template. The check stops the run instead, naming the file,
    the key and the two lines: for two amplicon: blocks, for a key written twice inside
    a nested block, and in a JSON config too.
  - It does not stop on what is not a repeat: the same key in two different blocks
    (enabled: under phylogeny and under taxonomy.filter, for example), a key written
    next to `<<: *name` (YAML's way of copying a block in, where a key written next to
    it replaces the copied one on purpose), the template config/config.yaml and every
    config of the dry-run matrix. A file the check cannot read is left alone: Snakemake
    has read it already (a JSON file indented with tabs, which JSON allows and YAML does
    not).

The dry-run matrix checks the same end to end: the config
.test/config/expected_stop/16S_two_amplicon_blocks.yaml must stop at start-up with
this message (the test below also runs the check on that file).
"""
import pytest
import yaml

from helpers import COMMON_SMK, REPO_ROOT, TEST_DIR, load_functions, write_text

repeated = load_functions(COMMON_SMK, [
    "_RepeatedKey", "_LoaderRefusingRepeatedKeys", "_check_no_repeated_keys"])

# A run config with a second amplicon: block pasted at its end (lines 2 and 9).
TWO_AMPLICON_BLOCKS = """\
mode: amplicon
amplicon:
  type: 16S
  primers:
    fwd: primers/fwd.fasta
    rev: primers/rev.fasta
output:
  out_dir: /data/out
amplicon:
  dada2:
    learn_errors:
      error_model: loess
"""

# The dry-run config that must stop with this check's message.
EXPECTED_STOP_CONFIG = TEST_DIR / "config" / "expected_stop" / "16S_two_amplicon_blocks.yaml"


def stop_text(paths):
    """Run the check where it must stop; return the stop text."""
    with pytest.raises(SystemExit) as stop:
        repeated._check_no_repeated_keys(paths)
    return str(stop.value)


def test_snakemakes_loader_keeps_only_the_last_amplicon_block():
    # The silent loss the check exists for: without it, the marker and the primers of
    # the first block are gone, and the run would use the template's.
    config = yaml.load(TWO_AMPLICON_BLOCKS, Loader=yaml.FullLoader)
    assert config["amplicon"] == {"dada2": {"learn_errors": {"error_model": "loess"}}}


def test_two_amplicon_blocks_stop_the_run_with_the_file_the_key_and_the_lines(tmp_path):
    path = write_text(tmp_path / "run_config.yaml", TWO_AMPLICON_BLOCKS)
    text = stop_text([path])
    assert text.startswith(f"[MetaFlux] {path}: the key 'amplicon' is written twice in the "
                           "same block (lines 2 and 9).")
    assert "YAML keeps only the last of the two, without a word" in text
    assert "Merge the two into one block, each setting once" in text


def test_a_key_written_twice_inside_a_nested_block_stops_the_run(tmp_path):
    path = write_text(tmp_path / "run_config.yaml",
                      "amplicon:\n  dada2:\n    max_ee: [2, 5]\n    trunc_q: 2\n    max_ee: [3, 6]\n")
    assert "the key 'max_ee' is written twice in the same block (lines 3 and 5)" in stop_text([path])


def test_a_json_config_with_a_repeated_key_stops_the_run(tmp_path):
    # Snakemake also takes JSON config files; JSON keeps the last key as YAML does.
    path = write_text(tmp_path / "run_config.json",
                      '{"amplicon": {"type": "16S"}, "mode": "amplicon", "amplicon": {"x": 1}}')
    assert "the key 'amplicon' is written twice" in stop_text([path])


def test_every_file_is_checked_not_only_the_first(tmp_path):
    good = write_text(tmp_path / "good.yaml", "mode: amplicon\n")
    bad = write_text(tmp_path / "bad.yaml", TWO_AMPLICON_BLOCKS)
    assert stop_text([good, bad]).startswith(f"[MetaFlux] {bad}: ")


def test_the_same_key_in_two_different_blocks_is_no_repeat(tmp_path):
    path = write_text(tmp_path / "run_config.yaml",
                      "amplicon:\n  phylogeny:\n    enabled: true\n  taxonomy:\n    filter:\n"
                      "      enabled: false\n")
    assert repeated._check_no_repeated_keys([path]) is None


def test_a_key_next_to_a_merged_block_is_no_repeat(tmp_path):
    # `<<: *defaults` copies the block in; `backend: fasttree` replaces the copied value.
    text = ("defaults: &defaults\n  enabled: true\n  backend: iqtree\n"
            "amplicon:\n  phylogeny:\n    <<: *defaults\n    backend: fasttree\n")
    path = write_text(tmp_path / "run_config.yaml", text)
    assert repeated._check_no_repeated_keys([path]) is None
    loaded = yaml.load(text, Loader=yaml.FullLoader)
    assert loaded["amplicon"]["phylogeny"] == {"enabled": True, "backend": "fasttree"}


def test_a_file_yaml_cannot_read_is_left_alone(tmp_path):
    # A JSON config indented with tabs: Snakemake reads it as JSON, PyYAML cannot.
    path = write_text(tmp_path / "run_config.json",
                      '{\n\t"mode": "amplicon",\n\t"amplicon": {\n\t\t"type": "16S"\n\t}\n}\n')
    assert repeated._check_no_repeated_keys([path]) is None


def test_the_template_and_the_dry_run_configs_have_no_repeated_key():
    # The template is not checked at start-up (it is MetaFlux's own file), so it is
    # checked here, with every config of the dry-run matrix but the one that must stop.
    config_files = [REPO_ROOT / "config" / "config.yaml"]
    config_files += sorted((TEST_DIR / "config").glob("*.yaml"))
    config_files += [path for path in sorted((TEST_DIR / "config" / "expected_stop").glob("*.yaml"))
                     if path != EXPECTED_STOP_CONFIG]
    assert len(config_files) > 20
    assert repeated._check_no_repeated_keys(config_files) is None


def test_the_dry_run_config_with_two_amplicon_blocks_stops_with_its_expected_text():
    # The same message the dry-run matrix looks for (the "# EXPECT:" line of the file).
    expected = next(line for line in EXPECTED_STOP_CONFIG.read_text().splitlines()
                    if line.startswith("# EXPECT:"))[len("# EXPECT:"):].strip()
    assert expected in stop_text([EXPECTED_STOP_CONFIG])
