"""Raw reads: the path every rule reads, the start-up check of the R2 files, and the
sample names of the raw read count.

raw_fastq() (00_common.smk) builds a sample's raw-read path from the naming and the
extension that sample discovery found at start-up, so that .fq and .fq.gz files are
found as well as .fastq and .fastq.gz. _samples_without_r2() (00_common.smk) is the
start-up check that every sample also has its second-read file, with that same
extension. Both are loaded from 00_common.smk with load_functions (see helpers.py).
raw_fastq reads three values that discovery sets at start-up (FASTQ_DIR, FASTQ_NAMING
and EXTN); the tests set them by hand.

Rule count_reads_raw (amplicon 30_preprocess.smk) takes each sample's name from its
raw R1 file name; the last test runs that line of the rule with bash, as the rule does.
"""
import subprocess
from pathlib import Path

from helpers import COMMON_SMK, WORKFLOW, load_functions

paths = load_functions(COMMON_SMK, ["_samples_without_r2", "raw_fastq"])

# The four extensions the start-up check accepts (ALLOWED_EXTENSIONS in 00_common.smk).
EXTENSIONS = ["fastq.gz", "fq.gz", "fastq", "fq"]


def set_discovery(fastq_dir, naming, extn):
    """Give raw_fastq the three values that sample discovery sets at start-up."""
    paths.raw_fastq.__globals__.update(FASTQ_DIR=Path(fastq_dir), FASTQ_NAMING=naming, EXTN=extn)


def touch(folder, names):
    """Create empty files with these names in folder (only their existence matters)."""
    for name in names:
        (folder / name).write_text("")


def test_path_carries_the_extension_found_at_start_up(tmp_path):
    for extn in EXTENSIONS:
        set_discovery(tmp_path, "_R1", extn)
        assert paths.raw_fastq("soil01", 1) == str(tmp_path / f"soil01_R1.{extn}")
        assert paths.raw_fastq("soil01", 2) == str(tmp_path / f"soil01_R2.{extn}")


def test_sra_style_naming_gives_1_and_2(tmp_path):
    set_discovery(tmp_path, "_1", "fq.gz")
    assert paths.raw_fastq("SRR11605259", 1) == str(tmp_path / "SRR11605259_1.fq.gz")
    assert paths.raw_fastq("SRR11605259", 2) == str(tmp_path / "SRR11605259_2.fq.gz")


def test_complete_pairs_pass_the_r2_check(tmp_path):
    touch(tmp_path, ["a_R1.fq.gz", "a_R2.fq.gz", "b_R1.fq.gz", "b_R2.fq.gz"])
    assert paths._samples_without_r2(tmp_path, ["a", "b"], "_R1", "fq.gz") == []


def test_missing_r2_and_r2_with_another_extension_are_reported(tmp_path):
    # a: no R2 at all; b: an uncompressed R2 next to a gzipped R1; c: complete.
    touch(tmp_path, ["a_R1.fastq.gz", "b_R1.fastq.gz", "b_R2.fastq", "c_R1.fastq.gz", "c_R2.fastq.gz"])
    assert paths._samples_without_r2(tmp_path, ["a", "b", "c"], "_R1", "fastq.gz") == ["a", "b"]


def test_r2_check_follows_sra_style_naming(tmp_path):
    touch(tmp_path, ["SRR1_1.fq", "SRR1_2.fq", "SRR2_1.fq"])
    assert paths._samples_without_r2(tmp_path, ["SRR1", "SRR2"], "_1", "fq") == ["SRR2"]


def test_raw_read_count_names_the_sample_for_every_extension():
    # The line of rule count_reads_raw that turns a raw R1 file name into the sample
    # name. The sample name holds an underscore of its own: only the final "_R1" and
    # the extension may go.
    # The rule's text runs from "rule count_reads_raw:" to the next rule (the two count
    # rules after it read MetaFlux's own .fastq.gz files and have a line of their own).
    rule_file = WORKFLOW / "rules" / "amplicon" / "30_preprocess.smk"
    after_start = rule_file.read_text().split("rule count_reads_raw:", 1)[1]
    rule_text = after_start.split("\nrule ", 1)[0]
    name_lines = [line.strip() for line in rule_text.splitlines()
                  if line.strip().startswith('sample=$(basename "$fq"')]
    assert len(name_lines) == 1, "expected one sample-name line in rule count_reads_raw"
    for extn in EXTENSIONS:
        script = f'fq=/data/run/soil_01_R1.{extn}\n{name_lines[0]}\necho "$sample"\n'
        result = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True)
        assert result.stdout.strip() == "soil_01"
