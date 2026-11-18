"""The NCBI run info of each sample, written next to samples.tsv in BioProject mode.

Rule sra_runinfo (shotgun 30_preprocess.smk) runs workflow/scripts/shotgun/
30a_sra_runinfo.py on the run info that MetaFlux caches in fastq_dir at start-up
(<BioProject>_runinfo.csv, comma-separated as NCBI writes it). The tests run the
script whole with run_script (see helpers.py) on a small run-info file written on the
spot, shaped like NCBI's: a header line repeated further down, a blank line, a
single-end run that MetaFlux does not process, and quoted values with a comma and a
tab inside.
"""
import pytest

from helpers import SCRIPTS, fake_snakemake, run_script, write_text

SCRIPT = SCRIPTS / "shotgun" / "30a_sra_runinfo.py"

RUNINFO = (
    "Run,LibraryName,LibraryLayout,BioSample,SampleName,CenterName\n"
    'SRR3,lib_c,PAIRED,SAMN3,"soil, plot 3",AAFC\n'
    "SRR1,lib_a,PAIRED,SAMN1,soil plot 1,AAFC\n"
    "\n"
    "Run,LibraryName,LibraryLayout,BioSample,SampleName,CenterName\n"
    'SRR2,lib_b,PAIRED,SAMN2,"soil\tplot 2",AAFC\n'
    "SRR9,lib_z,SINGLE,SAMN9,single-end run,AAFC\n"
)


def run_runinfo(tmp_path, samples, runinfo=RUNINFO):
    """Run the script on `runinfo` for `samples`; return the output lines."""
    sm = fake_snakemake(
        input={"runinfo": write_text(tmp_path / "PRJNA1_runinfo.csv", runinfo)},
        output={"table": str(tmp_path / "out" / "sra_runinfo.tsv")},
        params={"samples": samples},
        log=[str(tmp_path / "out" / "logs" / "sra_runinfo.log")],
    )
    run_script(SCRIPT, sm)
    return (tmp_path / "out" / "sra_runinfo.tsv").read_text().splitlines()


def test_one_row_per_sample_in_the_order_given_with_ncbi_columns(tmp_path):
    lines = run_runinfo(tmp_path, ["SRR1", "SRR2", "SRR3"])
    assert lines[0] == "sample\tRun\tLibraryName\tLibraryLayout\tBioSample\tSampleName\tCenterName"
    assert [line.split("\t")[0] for line in lines[1:]] == ["SRR1", "SRR2", "SRR3"]
    assert lines[1] == "SRR1\tSRR1\tlib_a\tPAIRED\tSAMN1\tsoil plot 1\tAAFC"


def test_runs_not_processed_and_repeated_headers_are_left_out(tmp_path):
    lines = run_runinfo(tmp_path, ["SRR1", "SRR2", "SRR3"])
    assert len(lines) == 4                       # header + three samples
    assert not any("SRR9" in line for line in lines)


def test_a_comma_is_kept_and_a_tab_becomes_a_space(tmp_path):
    lines = run_runinfo(tmp_path, ["SRR1", "SRR2", "SRR3"])
    fields = {line.split("\t")[0]: line.split("\t") for line in lines[1:]}
    assert fields["SRR3"][5] == "soil, plot 3"
    assert fields["SRR2"][5] == "soil plot 2"
    assert all(len(line.split("\t")) == 7 for line in lines)


def test_a_sample_without_a_row_stops_the_step(tmp_path):
    with pytest.raises(SystemExit) as stop:
        run_runinfo(tmp_path, ["SRR1", "SRR4"])
    assert "1 sample(s) have no row" in str(stop.value)
    assert "SRR4" in str(stop.value)
    assert not (tmp_path / "out" / "sra_runinfo.tsv").exists()
