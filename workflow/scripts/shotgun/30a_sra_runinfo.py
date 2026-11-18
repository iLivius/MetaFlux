#!/usr/bin/env python3
"""Write the NCBI run info of this run's samples next to samples.tsv (BioProject mode).

What it is for: in BioProject mode MetaFlux asks NCBI, once, for the run info of every
SRA run in the BioProject and keeps it in fastq_dir as <BioProject>_runinfo.csv
(_fetch_bioproject_runs in 00_common.smk), where it serves to pick the paired runs.
This script copies the rows of the samples this run processes into the output folder,
so that each sample's run accession, BioSample, sample name, library and instrument
travel with the results: the start of a sample metadata table for the analysis that
follows. NCBI's run info holds no BioSample attributes (habitat, host, place, date),
so neither does this table.

Input (from rule sra_runinfo, shotgun 30_preprocess.smk):
  runinfo  the cached run info: comma-separated text with a header line, one row per
           run, columns named by NCBI (Run, ReleaseDate, ..., BioSample, SampleName,
           ...). NCBI can repeat the header line further down; such lines are skipped,
           like every row of a run this MetaFlux run does not process (single-end runs,
           runs left out by input.accession_list).
  samples  the run accessions MetaFlux processes (SAMPLES), in samples.tsv order.

Output: sra_runinfo.tsv, tab-separated, one row per sample in the order given: first
the column 'sample' (the run accession, the sample name of every MetaFlux table), then
every column of the run info, named and filled as NCBI writes them. A tab or a line
break inside a value (not seen in NCBI's run info, but not ruled out) becomes a space,
so every row stays one line with the same number of columns. Requested by rule all;
nothing in the pipeline reads it back.

Stops (exit message, nothing written) when a sample has no row in the run info. That
does not happen in a normal run, since the samples are read from this very file; it
would mean the cached file was replaced after the run started.
"""
import csv
import sys
from pathlib import Path

sm = snakemake  # noqa: F821

runinfo_path = Path(sm.input.runinfo)
table_path   = Path(sm.output.table)
samples      = list(sm.params.samples)
log_path     = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log = log_path.open("w")


def clean(value):
    """One value as it goes into the table: a tab or a line break becomes a space."""
    if value is None:
        return ""
    return value.replace("\t", " ").replace("\r", " ").replace("\n", " ")


# ── Read the run info: one dict per run, keyed by the run accession ──────────────
# The first row of each run accession is kept; a repeated header line (its Run cell
# reads "Run") is skipped. A second row for the same run is counted for the log.
rows_by_run = {}
n_repeated_runs = 0
with runinfo_path.open(newline="") as handle:
    reader = csv.DictReader(handle)
    columns = list(reader.fieldnames or [])
    for row in reader:
        run = row.get("Run") or ""
        if not run or run == "Run":
            continue
        if run in rows_by_run:
            n_repeated_runs += 1
            continue
        rows_by_run[run] = row

log.write(f"[sra_runinfo] {runinfo_path}: {len(rows_by_run)} runs, {len(columns)} columns\n")
if n_repeated_runs:
    log.write(f"[sra_runinfo] {n_repeated_runs} repeated row(s) of a run already read were skipped\n")

# ── Every sample must have its row ─────────────────────────────────────────────
missing = [sample for sample in samples if sample not in rows_by_run]
if missing:
    message = (f"[sra_runinfo] {len(missing)} sample(s) have no row in {runinfo_path}: "
               f"{', '.join(missing[:5])}. The samples of a BioProject run are read from "
               f"this file at start-up, so it was probably replaced after the run started; "
               f"start the run again.")
    log.write(message + "\n")
    log.close()
    sys.exit(message)

# ── Write the table ──────────────────────────────────────────────────────────────
table_path.parent.mkdir(parents=True, exist_ok=True)
with table_path.open("w") as out:
    out.write("\t".join(["sample"] + [clean(column) for column in columns]) + "\n")
    for sample in samples:
        row = rows_by_run[sample]
        values = [clean(row.get(column)) for column in columns]
        out.write("\t".join([sample] + values) + "\n")

log.write(f"[sra_runinfo] wrote {table_path}: {len(samples)} sample(s)\n")
log.close()
