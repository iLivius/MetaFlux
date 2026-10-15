"""The read-tracking tables: stats/read_tracking.txt in both modes.

What is tested: three scripts turn the per-sample counts of every pipeline step
into one table, one row per sample, one column per step, so that a sample that
loses most of its reads somewhere can be seen at a glance:
  amplicon  workflow/scripts/amplicon/80b_aggregate_read_counts.py, reading the
            count files of the preprocessing steps, DADA2's read.counts table, the
            sequence tables and the final asv_table.txt;
  shotgun   workflow/scripts/shotgun/90a_aggregate_read_counts.py (Kraken2 runs),
            reading BBDuk's PhiX log, BBMap's host-removal count, fastp's JSON and
            the Kraken2 report;
            workflow/scripts/shotgun/38b_sylph_read_tracking.py (sylph only), which
            keeps a word-for-word copy of 90a's parsers for the cleaning columns.
The tests feed each parser a tiny file in the real format and compare the numbers;
run the two shotgun scripts end to end on the same files, so their shared columns
must agree; and check that 38b's copies of the parsers are still identical to 90a's,
which the scripts ask for in a comment and nothing enforced until now.
"""
import json
from pathlib import Path

import pytest

from helpers import SCRIPTS, fake_snakemake, load_functions, run_script, source_segments, write_text

AMPLICON_SCRIPT = SCRIPTS / "amplicon" / "80b_aggregate_read_counts.py"
KRAKEN2_SCRIPT  = SCRIPTS / "shotgun" / "90a_aggregate_read_counts.py"
SYLPH_SCRIPT    = SCRIPTS / "shotgun" / "38b_sylph_read_tracking.py"

amplicon = load_functions(AMPLICON_SCRIPT, [
    "parse_stage_counts", "_unq", "parse_dada2_counts", "parse_seqtab_rowsums",
    "parse_asv_table_colsums",
])
shotgun = load_functions(KRAKEN2_SCRIPT, [
    "parse_bbduk_log", "parse_host_surviving", "parse_fastp_counts", "parse_kraken2_report",
    "_index",
])

# ── Fixture texts, in the formats the real tools write ───────────────────────────

# A count file of the amplicon preprocessing: one line per FASTQ, then a total.
COUNT_FILE = ("/run/01.reads/S1_R1.fastq.gz : 1000\n"
              "/run/01.reads/S2_R1.fastq.gz : 800\n"
              "Total read count : 1800\n")

# DADA2's read.counts, as R's write.table writes it (quoted names, empty first cell).
DADA2_COUNTS = ('""\t"stripped"\t"filtered"\t"denoised"\t"merged"\t"non_chimeric"\n'
                '"S1"\t990\t900\t890\t850\t840\n'
                '"S2"\t790\t700\t690\t650\t640\n')

# A sequence table: samples in rows, ASVs in columns.
SEQTAB = ('""\t"ASV_1"\t"ASV_2"\t"ASV_3"\n'
          '"S1"\t500\t300\t\n'
          '"S2"\t400\t\t200\n')

# The final ASV table: ASVs in rows, samples and the taxonomy string in columns.
ASV_TABLE = ('""\t"S1"\t"S2"\t"taxonomy"\n'
             '"ASV_1"\t500\t400\t"d__Bacteria;p__Pseudomonadota"\n'
             '"ASV_3"\t\t200\t"d__Bacteria;p__Bacillota"\n')

# The end of a BBDuk log: both lines count reads, so pairs are half of each.
BBDUK_LOG = ("Input is being processed as paired\n"
             "Input:                  600 reads           90000 bases.\n"
             "Contaminants:           2 reads (0.33%)     300 bases (0.33%)\n"
             "Result:                 598 reads (99.67%)  89700 bases (99.67%)\n")

DEHOST_STATS = "surviving_pairs\t290\n"

FASTP_JSON = json.dumps({
    "summary": {"before_filtering": {"total_reads": 598},
                "after_filtering": {"read1_mean_length": 147, "read2_mean_length": 146}},
    "filtering_result": {"passed_filter_reads": 590},
})

# A Kraken2 report written with --report-minimizer-data: eight columns, the rank
# code in the sixth. Unclassified first, then the root row with the classified total.
KRAKEN2_REPORT = (" 16.95\t50\t50\t0\t0\tU\t0\tunclassified\n"
                  " 83.05\t245\t0\t1200\t900\tR\t1\troot\n"
                  " 83.05\t245\t0\t1200\t900\tR1\t131567\t  cellular organisms\n"
                  " 83.05\t245\t245\t1200\t900\tD\t2\t    Bacteria\n")


# ── The amplicon parsers ─────────────────────────────────────────────────────────

def written(path, text):
    """Write a fixture file and return its path as a Path object, which is what the
    scripts hand their parsers."""
    return Path(write_text(path, text))


def test_amplicon_parsers_read_each_format(tmp_path):
    counts = amplicon.parse_stage_counts(written(tmp_path / "counts.txt", COUNT_FILE))
    assert counts == {"/run/01.reads/S1_R1.fastq.gz": 1000, "/run/01.reads/S2_R1.fastq.gz": 800}
    dada2 = amplicon.parse_dada2_counts(written(tmp_path / "read.counts", DADA2_COUNTS))
    assert dada2["S1"] == {"stripped": 990, "filtered": 900, "denoised": 890,
                           "merged": 850, "non_chimeric": 840}
    assert amplicon.parse_seqtab_rowsums(written(tmp_path / "seqtab.txt", SEQTAB)) == {
        "S1": 800, "S2": 600}
    assert amplicon.parse_asv_table_colsums(written(tmp_path / "asv_table.txt", ASV_TABLE)) == {
        "S1": 500, "S2": 600}
    assert amplicon._unq('"S1"') == "S1"


def test_amplicon_table_end_to_end_with_a_missing_sample(tmp_path):
    # The count files key their lines by file path, which here never equals a sample
    # name, so the raw/nophix/stripped cells read NA: the table never invents a 0 for
    # a sample a source file does not mention. The DADA2 and table columns are filled.
    sm = fake_snakemake(
        input={"raw_counts": write_text(tmp_path / "01.raw_reads.counts", COUNT_FILE),
               "nophix_counts": write_text(tmp_path / "02.nophix_reads.counts", COUNT_FILE),
               "stripped_counts": write_text(tmp_path / "03.stripped_reads.counts", COUNT_FILE),
               "dada2_counts": write_text(tmp_path / "read.counts", DADA2_COUNTS),
               "seqtab_lenfilt": write_text(tmp_path / "seqtab_lenfilt.txt", SEQTAB),
               "asv_table": write_text(tmp_path / "asv_table.txt", ASV_TABLE),
               "seqtab_extracted": []},
        output={"counts": str(tmp_path / "read_tracking.txt")},
        params={"samples": ["S1", "S2", "S3"], "extraction_enabled": False, "remove_phix": True},
        log=[str(tmp_path / "log.txt")],
    )
    run_script(AMPLICON_SCRIPT, sm)
    lines = (tmp_path / "read_tracking.txt").read_text().splitlines()
    assert lines[0] == ("sample\traw\tnophix\tstripped\tfiltered\tdenoised\tmerged\t"
                        "non_chimeric\tpost_length_filter\tpost_taxonomy_filter")
    assert lines[1] == "S1\tNA\tNA\tNA\t900\t890\t850\t840\t800\t500"
    assert lines[2] == "S2\tNA\tNA\tNA\t700\t690\t650\t640\t600\t600"
    assert lines[3] == "S3\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"


# ── The shotgun parsers ──────────────────────────────────────────────────────────

def test_shotgun_parsers_count_pairs(tmp_path):
    assert shotgun.parse_bbduk_log(written(tmp_path / "S1.log", BBDUK_LOG)) == (300, 299)
    assert shotgun.parse_host_surviving(written(tmp_path / "S1_dehost_stats.txt", DEHOST_STATS)) == 290
    assert shotgun.parse_fastp_counts(written(tmp_path / "S1_fastp.json", FASTP_JSON)) == (299, 295)
    assert shotgun.parse_kraken2_report(written(tmp_path / "S1_report.txt", KRAKEN2_REPORT)) == (245, 50)


def test_bbduk_log_without_its_summary_raises(tmp_path):
    # A log cut short (BBDuk did not finish) must not turn into a made-up number.
    with pytest.raises(ValueError):
        shotgun.parse_bbduk_log(written(tmp_path / "S1.log", "Input is being processed as paired\n"))


def test_index_strips_the_suffix_from_the_end_only(tmp_path):
    # A sample called Plot.log2 has the PhiX log Plot.log2.log; only the final
    # ".log" may go, or the sample could not be found.
    index = shotgun._index(["/x/Plot.log2.log", "/x/S1.log"], ".log")
    assert sorted(index) == ["Plot.log2", "S1"]
    assert str(index["Plot.log2"]).endswith("Plot.log2.log")


def shotgun_inputs(tmp_path):
    """The per-sample files of one sample, S1, as both shotgun scripts read them."""
    return {
        "fastp_jsons": [write_text(tmp_path / "S1_fastp.json", FASTP_JSON)],
        "phix_logs": [write_text(tmp_path / "S1.log", BBDUK_LOG)],
        "host_stats": [write_text(tmp_path / "S1_dehost_stats.txt", DEHOST_STATS)],
        "kraken_reports": [write_text(tmp_path / "S1_report.txt", KRAKEN2_REPORT)],
    }


def test_kraken2_and_sylph_tables_agree_on_the_cleaning_columns(tmp_path):
    inputs = shotgun_inputs(tmp_path)
    params = {"samples": ["S1"], "remove_phix": True, "remove_host": True}
    run_script(KRAKEN2_SCRIPT, fake_snakemake(
        input=inputs, output={"counts": str(tmp_path / "kraken2_table.txt")},
        params=params, log=[str(tmp_path / "k.log")]))
    run_script(SYLPH_SCRIPT, fake_snakemake(
        input=inputs, output={"counts": str(tmp_path / "sylph_table.txt")},
        params=params, log=[str(tmp_path / "s.log")]))
    kraken2_lines = (tmp_path / "kraken2_table.txt").read_text().splitlines()
    sylph_lines = (tmp_path / "sylph_table.txt").read_text().splitlines()
    assert kraken2_lines == ["sample\traw\tnophix\tno_host\ttrimmed\tclassified\tunclassified",
                             "S1\t300\t299\t290\t295\t245\t50"]
    assert sylph_lines == ["sample\traw\tnophix\tno_host\ttrimmed",
                           "S1\t300\t299\t290\t295"]


def test_without_phix_removal_raw_comes_from_fastp(tmp_path):
    # No BBDuk step: fastp is the first tool to see the reads, so its "before"
    # count is the raw count, and the nophix and no_host columns do not exist.
    inputs = shotgun_inputs(tmp_path)
    inputs["phix_logs"] = []
    inputs["host_stats"] = []
    run_script(KRAKEN2_SCRIPT, fake_snakemake(
        input=inputs, output={"counts": str(tmp_path / "table.txt")},
        params={"samples": ["S1"], "remove_phix": False, "remove_host": False},
        log=[str(tmp_path / "k.log")]))
    assert (tmp_path / "table.txt").read_text().splitlines() == [
        "sample\traw\ttrimmed\tclassified\tunclassified", "S1\t299\t295\t245\t50"]


def test_sylph_script_keeps_word_for_word_copies_of_the_parsers():
    # 38b cannot import 90a (90a builds its table as soon as it is loaded), so it
    # carries copies of four parsers. Both files ask, in a comment, that the copies
    # stay identical; this test is that check.
    shared = ["parse_bbduk_log", "parse_host_surviving", "parse_fastp_counts", "_index"]
    in_kraken2 = source_segments(KRAKEN2_SCRIPT, shared)
    in_sylph = source_segments(SYLPH_SCRIPT, shared)
    for name in shared:
        assert in_kraken2[name] == in_sylph[name], f"{name} differs between 90a and 38b"
