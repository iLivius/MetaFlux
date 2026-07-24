#!/usr/bin/env python3
"""Run Metaxa2 on ASV sequences, collect all extracted SSU sequences,
and subset seqtab_head_names.txt to the retained ASV IDs.

Metaxa2 is run with -t all (all domains); the pipeline keeps every extracted
sequence regardless of domain — taxonomic gatekeeping happens post-assignment
via the contaminant filter in assign_taxonomy.R.
"""
import subprocess
import sys
from pathlib import Path


def subset_seqtab_by_ids(in_path: Path, kept_ids: set[str], out_path: Path, log) -> None:
    with in_path.open() as fh:
        raw_header = fh.readline().rstrip("\n")
        data_lines = fh.readlines()
    header = raw_header.split("\t")
    # R write.table wraps every column name in double quotes; strip them for matching.
    header_clean = [h.strip('"') for h in header]
    keep_idx = [0] + [i for i, h in enumerate(header_clean[1:], 1) if h in kept_ids]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("\t".join(header[i] for i in keep_idx) + "\n")
        for line in data_lines:
            cols = line.rstrip("\n").split("\t")
            fh.write("\t".join(cols[i] for i in keep_idx) + "\n")
    log.write(f"[metaxa2_extract] seqtab: kept {len(keep_idx)-1} / {len(header)-1} ASV columns\n")


sm = snakemake  # noqa: F821

input_fasta      = str(sm.input.seqs)
seqtab_names_in  = Path(sm.input.seqtab_names)
output_fasta     = Path(sm.output.seqs)
seqtab_names_out = Path(sm.output.seqtab_names)
results_out      = Path(sm.output.results)
prefix           = str(sm.params.prefix)
threads          = int(sm.threads)

log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log = log_path.open("w")

Path(prefix).parent.mkdir(parents=True, exist_ok=True)

cmd = [
    "metaxa2",
    "-i", input_fasta,
    "-o", prefix,
    "-t", "all",
    "--cpu", str(threads),
    "--fasta", "T",
    "--table", "T",
    "--plus", "T",
    "--silent", "T",
]
log.write(f"[metaxa2_extract] Running: {' '.join(cmd)}\n")
log.flush()
res = subprocess.run(cmd, stdout=log, stderr=log)
if res.returncode != 0:
    log.write(f"[metaxa2_extract] ERROR: metaxa2 exited with code {res.returncode}\n")
    log.close()
    sys.exit(res.returncode)

# Collect the combined extraction FASTA (all domains)
extraction_fasta = Path(f"{prefix}.extraction.fasta")
if not extraction_fasta.exists():
    log.write(f"[metaxa2_extract] ERROR: {extraction_fasta} not found\n")
    log.close()
    sys.exit(1)

# Write cleaned FASTA: keep only first token of each header as ASV ID
kept_ids: set[str] = set()
output_fasta.parent.mkdir(parents=True, exist_ok=True)
with extraction_fasta.open() as fi, output_fasta.open("w") as fo:
    current_seq: list[str] = []
    current_id: str | None = None
    for line in fi:
        if line.startswith(">"):
            if current_id is not None:
                fo.write(f">{current_id}\n{''.join(current_seq)}\n")
            current_id = line[1:].split()[0].split("|")[0]
            kept_ids.add(current_id)
            current_seq = []
        else:
            current_seq.append(line.rstrip("\n"))
    if current_id is not None:
        fo.write(f">{current_id}\n{''.join(current_seq)}\n")

log.write(f"[metaxa2_extract] {len(kept_ids)} ASV(s) extracted\n")

# Copy Metaxa2's results summary (a QC report, NOT consumed by any downstream
# rule). The real extraction output — the FASTA above — is already validated, so a
# missing summary must not fail an otherwise-complete run; but warn loudly and leave
# a self-explaining placeholder rather than a silent empty file.
# (Same missing-summary policy as 70b_itsx_extract.py — tolerant but loud.)
results_src = Path(f"{prefix}.extraction.results")
results_out.parent.mkdir(parents=True, exist_ok=True)
if results_src.exists():
    results_out.write_text(results_src.read_text())
else:
    log.write(
        f"[metaxa2_extract] WARNING: expected summary {results_src} not found; "
        "the extraction FASTA was produced, so continuing. Writing a placeholder.\n"
    )
    results_out.write_text(
        f"# Metaxa2 summary ({results_src.name}) was not produced by this run.\n"
        "# Extraction itself succeeded (see the extracted FASTA); this QC report is absent.\n"
    )

# Subset seqtab
subset_seqtab_by_ids(seqtab_names_in, kept_ids, seqtab_names_out, log)

# Keep all Metaxa2 raw outputs (seqs.summary.txt, seqs.hmmer.table, seqs.graph,
# seqs.extraction.fasta, seqs.extraction.fasta.1, seqs.extraction.results) under
# the prefix dir for QC / downstream inspection. Final dir layout:
#   <prefix_dir>/seqs.*  ← raw Metaxa2 working files (persisted)
#   5.dada2/seqs_extracted.fasta             ← promoted, cleaned headers
#   5.dada2/metaxa2_extraction.results.txt   ← promoted .extraction.results
log.write(f"[metaxa2_extract] Keeping raw outputs under: {Path(prefix).parent}\n")
log.write(f"[metaxa2_extract] DONE.\n")
log.close()
