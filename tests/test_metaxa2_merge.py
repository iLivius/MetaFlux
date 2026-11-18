"""The 16S extraction step: ASVs that come out of Metaxa2 identical are merged.

Rule target_extract (amplicon 70_extract.smk) runs workflow/scripts/amplicon/
70a_metaxa2_extract.py. Metaxa2 seldom trims a 16S ASV, but when it does, two ASVs that
differed only in the bases it cut away come out with the same sequence. The script then
merges them into the lowest-numbered ASV (the most abundant: DADA2 numbers ASVs by
decreasing abundance) and adds up their counts, as the ITS step does after ITSx.

The tests run the script whole with run_script (see helpers.py). Metaxa2 itself is not
needed: a stand-in program named `metaxa2`, put first on PATH, copies into the folder
it runs in the output files each test prepared (the extraction FASTA, its results
table and the six per-origin FASTA files, all named seqs.*, as Metaxa2 names them for
the prefix "seqs"). Metaxa2 writes its records sorted by ID as text, so ASV_10 comes
before ASV_2; the inputs below follow that order.
"""
import os
import stat

from helpers import SCRIPTS, fake_snakemake, run_script, write_text

SCRIPT = SCRIPTS / "amplicon" / "70a_metaxa2_extract.py"

ORIGINS = ["mitochondria", "chloroplast", "bacteria", "archaea", "eukaryota", "uncertain"]

# DADA2's count table: one row per sample, one column per ASV, the header quoted by R.
# ASV_3 is the one Metaxa2 will not extract.
SEQTAB = ('""\t"ASV_1"\t"ASV_2"\t"ASV_3"\t"ASV_10"\n'
          '"S1"\t100\t40\t7\t5\n'
          '"S2"\t80\t0\t3\t9\n')
SEQS_IN = ">ASV_1\nAAAA\n>ASV_2\nCCCCGG\n>ASV_3\nTTTT\n>ASV_10\nCCCCTT\n"


def header(asv_id, length, letter="B"):
    """A header line as Metaxa2 writes it in its extraction FASTA."""
    return (f">{asv_id}|{letter} Predicted Bacterial  rRNA ({length} bp) "
            "From domain V4l to V4r on main strand")


def run_extract(tmp_path, monkeypatch, extracted, origin_of):
    """Run the script on SEQTAB, with the stand-in Metaxa2 returning `extracted`
    (a list of (ASV ID, extracted sequence), in Metaxa2's order) and the origin
    calls `origin_of` ({ASV ID: origin}). Returns the output folder."""
    fake_out = tmp_path / "fake_metaxa2_output"
    fasta = "".join(f"{header(aid, len(seq))}\n{seq}\n" for aid, seq in extracted)
    write_text(fake_out / "seqs.extraction.fasta", fasta)
    write_text(fake_out / "seqs.extraction.results", "results table\n")
    for origin in ORIGINS:
        records = "".join(f">{aid} SSU rRNA\n{seq}\n" for aid, seq in extracted
                          if origin_of.get(aid) == origin)
        write_text(fake_out / f"seqs.{origin}.fasta", records)

    # The stand-in `metaxa2`: copy the prepared files into the folder it runs in.
    bin_dir = tmp_path / "bin"
    program = bin_dir / "metaxa2"
    write_text(program, f'#!/bin/sh\ncp "{fake_out}"/seqs.* .\n')
    program.chmod(program.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")

    out = tmp_path / "05.dada2"
    sm = fake_snakemake(
        input={"seqs": write_text(tmp_path / "in" / "seqs.fasta", SEQS_IN),
               "seqtab_names": write_text(tmp_path / "in" / "seqtab_head_names.txt", SEQTAB)},
        output={"seqs": str(out / "seqs_extracted.fasta"),
                "seqtab_names": str(out / "seqtab_extracted_head_names.txt"),
                "results": str(out / "metaxa2_extraction.results.txt"),
                "origin": str(out / "metaxa2_origin.tsv"),
                "collapse_map": str(out / "metaxa2_collapse_map.tsv")},
        params={"prefix": str(out / "metaxa2" / "seqs")},
        log=[str(tmp_path / "logs" / "target_extract.log")],
    )
    sm.threads = 1
    run_script(SCRIPT, sm)
    return out


def test_nothing_is_merged_when_every_extracted_sequence_differs(tmp_path, monkeypatch):
    # The output is that of the step before the merge existed: Metaxa2's records in
    # its own order with plain IDs, and the count table without the ASV it dropped.
    out = run_extract(tmp_path, monkeypatch,
                      extracted=[("ASV_1", "AAAA"), ("ASV_10", "CCCCTT"), ("ASV_2", "CCCCGG")],
                      origin_of={"ASV_1": "bacteria", "ASV_10": "bacteria", "ASV_2": "chloroplast"})
    assert (out / "seqs_extracted.fasta").read_text() == (
        ">ASV_1\nAAAA\n>ASV_10\nCCCCTT\n>ASV_2\nCCCCGG\n")
    assert (out / "seqtab_extracted_head_names.txt").read_text() == (
        '""\t"ASV_1"\t"ASV_2"\t"ASV_10"\n'
        '"S1"\t100\t40\t5\n'
        '"S2"\t80\t0\t9\n')
    assert (out / "metaxa2_origin.tsv").read_text() == (
        "asv_id\torigin\nASV_1\tbacteria\nASV_10\tbacteria\nASV_2\tchloroplast\n")
    assert (out / "metaxa2_collapse_map.tsv").read_text().splitlines() == [
        "representative\tmerged_asvs\tn_merged\textracted_length",
        "ASV_1\tASV_1\t1\t4",
        "ASV_2\tASV_2\t1\t6",
        "ASV_10\tASV_10\t1\t6",
    ]


def test_asvs_trimmed_to_one_sequence_are_merged_into_the_most_abundant(tmp_path, monkeypatch):
    # Metaxa2 cut the last two bases of ASV_2 and ASV_10, the only ones in which they
    # differed. ASV_10 comes first in Metaxa2's order, but ASV_2 is the more abundant
    # (lower number), so ASV_2 stands for both, in its own place.
    out = run_extract(tmp_path, monkeypatch,
                      extracted=[("ASV_1", "AAAA"), ("ASV_10", "CCCC"), ("ASV_2", "CCCC")],
                      origin_of={"ASV_1": "bacteria", "ASV_10": "bacteria", "ASV_2": "bacteria"})
    assert (out / "seqs_extracted.fasta").read_text() == ">ASV_1\nAAAA\n>ASV_2\nCCCC\n"
    # ASV_2's column holds ASV_2 + ASV_10 (40 + 5, 0 + 9); every sample keeps its reads.
    assert (out / "seqtab_extracted_head_names.txt").read_text() == (
        '""\t"ASV_1"\t"ASV_2"\n'
        '"S1"\t100\t45\n'
        '"S2"\t80\t9\n')
    assert (out / "metaxa2_origin.tsv").read_text() == (
        "asv_id\torigin\nASV_1\tbacteria\nASV_2\tbacteria\n")
    assert (out / "metaxa2_collapse_map.tsv").read_text().splitlines() == [
        "representative\tmerged_asvs\tn_merged\textracted_length",
        "ASV_1\tASV_1\t1\t4",
        "ASV_2\tASV_2,ASV_10\t2\t4",
    ]
    log = (tmp_path / "logs" / "target_extract.log").read_text()
    assert "3 extracted ASV(s) -> 2 distinct sequence(s); 1 group(s)" in log
    assert "ASV_2 <- ASV_2, ASV_10 (4 bp after extraction)" in log
    assert "kept 2 / 4 ASV columns (1 not extracted by Metaxa2, 1 merged into another ASV)" in log


def test_a_merged_group_with_two_origin_calls_keeps_the_representatives(tmp_path, monkeypatch):
    # Not expected (the ASVs of a group share one sequence), but if Metaxa2 called them
    # differently, the representative's call is kept and the log says so.
    out = run_extract(tmp_path, monkeypatch,
                      extracted=[("ASV_1", "AAAA"), ("ASV_10", "CCCC"), ("ASV_2", "CCCC")],
                      origin_of={"ASV_1": "bacteria", "ASV_10": "chloroplast", "ASV_2": "bacteria"})
    assert "ASV_2\tbacteria" in (out / "metaxa2_origin.tsv").read_text().splitlines()
    log = (tmp_path / "logs" / "target_extract.log").read_text()
    assert ("WARNING: the ASVs merged into ASV_2 have different origin calls "
            "(ASV_2 bacteria, ASV_10 chloroplast); ASV_2 keeps its own, bacteria") in log
