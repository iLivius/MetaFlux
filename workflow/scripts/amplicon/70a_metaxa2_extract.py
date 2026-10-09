#!/usr/bin/env python3
"""Metaxa2: confirm and trim ASVs down to their actual SSU rRNA region.

This is the 16S "target_extract" step (see 70_extract.smk), run only when
amplicon.extraction.enabled is true. Metaxa2 uses HMM profiles — precomputed
statistical models of what a real SSU rRNA sequence looks like, built one per
domain (bacteria, archaea, eukaryote, plus mitochondrial and chloroplast SSU,
which are bacterial in origin) — to find the true rRNA boundaries inside each
input sequence and say which domain it best matches. We run it here because
DADA2's ASVs (amplicon sequence variants — DADA2's exact-sequence output unit)
are not guaranteed to all be genuine 16S: some can be non-specific amplification
or other artifacts that survived denoising. Metaxa2 re-confirms the sequence is
real rRNA and, where the ASV runs slightly past the true rRNA boundary (e.g.
into a flanking region), trims it back to just that region.

Input: seqs.fasta and seqtab_head_names.txt, straight from dada_seqtab.R — the
full, pre-length-filter ASV set (see rule target_extract in 70_extract.smk).

Run with -t all (all domains); every sequence Metaxa2 extracts is kept here,
whatever its origin — nothing is dropped for being non-bacterial. What is
removed is decided after taxonomy assignment, by the taxonomy filter
(80a_assign_taxonomy.R for the rdp path, 80c_parse_sintax.py for the sintax
path). For domains that filter goes by the classified lineage (the keep
tokens). For organelles it also uses Metaxa2's own origin call, which this
script writes to metaxa2_origin.tsv (amplicon.taxonomy.filter.organelles,
default metaxa2): on the 16S test set, 8 plant-mitochondrial ASVs whose
nearest SILVA reference is a mitochondrion at 97-99 % identity were
classified as bacteria (mostly Enterobacterales) and passed every token, while
Metaxa2 called all 8 mitochondrial.

Output: seqs_extracted.fasta (cleaned, trimmed sequences) and
seqtab_extracted_head_names.txt (the matching read-count table, subset to the
ASVs Metaxa2 could extract) are read next by dada_length_filter
(60d_dada_length_filter.py), which applies the ASV length window.
metaxa2_origin.tsv (one row per extracted ASV: asv_id, origin) is read by
assign_taxonomy (see above). Metaxa2's own raw output files are kept under the
working prefix directory for QC.
"""
import subprocess
import sys
from pathlib import Path


def subset_seqtab_by_ids(in_path: Path, kept_ids: set[str], out_path: Path, log) -> None:
    """Drop seqtab columns (ASVs) that Metaxa2 did not extract.

    in_path is the full seqtab_head_names.txt from dada_seqtab.R (every ASV,
    every sample's read count). kept_ids are the ASV IDs Metaxa2 did extract
    (built below, while parsing its output FASTA). Any ASV not in kept_ids —
    Metaxa2 could not confirm it as SSU rRNA in any domain — has its whole
    column removed, so its read counts drop out of every downstream table.
    """
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
origin_out       = Path(sm.output.origin)        # metaxa2_origin.tsv, see the end of this script
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
    "-t", "all",          # search all domain HMM profiles (bacteria/archaea/eukaryote)
    "--cpu", str(threads),
    "--fasta", "T",       # write the extracted-sequence FASTA (*.extraction.fasta)
    "--table", "T",       # write the per-sequence domain/taxonomy table
    "--plus", "T",        # also detect chloroplast/mitochondrial SSU (organelle-derived,
                          # bacterial in origin) rather than treating them as unclassifiable
    "--silent", "T",      # suppress Metaxa2's own console output (everything goes to our log)
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

# Write cleaned FASTA: Metaxa2 writes each header as the ASV ID immediately
# followed by "|<rRNA subunit tag>", then a space and a free-text description,
# e.g. ">ASV_12|SSU Bacterial 16S rRNA (253 bp) From domain 1 to 253 on +".
# Taking the first whitespace token gives "ASV_12|SSU"; splitting that on "|"
# isolates the plain ASV ID — this is what has to match the seqtab column
# headers for the subset below.
kept_ids: set[str] = set()
kept_order: list[str] = []      # the same IDs in file order, for metaxa2_origin.tsv
output_fasta.parent.mkdir(parents=True, exist_ok=True)
with extraction_fasta.open() as fi, output_fasta.open("w") as fo:
    current_seq: list[str] = []
    current_id: str | None = None
    for line in fi:
        if line.startswith(">"):
            if current_id is not None:
                fo.write(f">{current_id}\n{''.join(current_seq)}\n")
            current_id = line[1:].split()[0].split("|")[0]
            if current_id not in kept_ids:
                kept_order.append(current_id)
            kept_ids.add(current_id)
            current_seq = []
        else:
            current_seq.append(line.rstrip("\n"))
    if current_id is not None:
        fo.write(f">{current_id}\n{''.join(current_seq)}\n")

log.write(f"[metaxa2_extract] {len(kept_ids)} ASV(s) extracted\n")

# Copy Metaxa2's extraction results table: one row per input sequence, with
# extracted length, rRNA subunit type, completeness, match e-value/score, and
# extraction coordinates (plus a chimera flag where relevant). It is a QC
# reference only — NOT consumed by any downstream rule. The real extraction
# output — the FASTA above — is already validated, so a missing results file
# must not fail an otherwise-complete run; but warn loudly and leave a
# self-explaining placeholder rather than a silent empty file.
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

# ── Metaxa2's origin call per ASV → metaxa2_origin.tsv ──────────────────────
# After extracting, Metaxa2 aligns every extracted region against its own
# reference database (a BLAST step) and writes each sequence into one per-origin
# FASTA: <prefix>.mitochondria.fasta, .chloroplast.fasta, .bacteria.fasta,
# .archaea.fasta, .eukaryota.fasta, .uncertain.fasta. That second, alignment-
# based call is the one used here. The single letter in the extraction headers
# above ("ASV_1|M ...") is only Metaxa2's preliminary HMM guess, which calls far
# too many sequences mitochondrial (see "Metaxa2 reports origin twice" in
# docs/amplicon/markers/16S.md), so it is not used.
#
# The per-origin headers look like ">ASV_101 Mitochondrial SSU rRNA"; the first
# word is the ASV ID. The two organelle files are read first: if Metaxa2 ever
# listed an ASV in two files, the first call is kept (with a warning in the log),
# so an ASV listed as an organelle anywhere counts as an organelle.
#
# Output: one row per extracted ASV, in the order of Metaxa2's extraction FASTA
# (a lookup table: nothing depends on the row order), columns asv_id and origin
# (mitochondria, chloroplast, bacteria, archaea, eukaryota, uncertain; NA when
# Metaxa2 extracted the sequence but listed it in none of the files). Read by
# rule assign_taxonomy: the organelle filter (amplicon.taxonomy.filter.organelles)
# removes the mitochondria and chloroplast ASVs, and every ASV's call goes into the
# metaxa2_origin column of 6.taxonomy/asv_table_unfiltered.txt.
ORIGINS = ["mitochondria", "chloroplast", "bacteria", "archaea", "eukaryota", "uncertain"]
origin_of: dict[str, str] = {}
for origin in ORIGINS:
    origin_fasta = Path(f"{prefix}.{origin}.fasta")
    # Metaxa2 writes all six files on every -t all run, empty ones included. A
    # missing file means its output is not what this script expects, and an origin
    # table written without it could let organelles through unnoticed, so stop.
    if not origin_fasta.exists():
        log.write(
            f"[metaxa2_extract] ERROR: {origin_fasta} not found. Metaxa2 writes one FASTA "
            "per origin; without it the organelle filter cannot be built.\n"
        )
        log.close()
        sys.exit(1)
    with origin_fasta.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            asv_id = line[1:].split()[0].split("|")[0]
            if asv_id in origin_of:
                log.write(
                    f"[metaxa2_extract] WARNING: {asv_id} is listed as {origin_of[asv_id]} "
                    f"and as {origin}; keeping {origin_of[asv_id]}\n"
                )
                continue
            origin_of[asv_id] = origin

# An ID in an origin file that is not in the extraction FASTA should not happen
# (the files are built from the same extraction); say so if it ever does. Such an
# ID is not written, since that ASV was dropped above and never reaches taxonomy.
not_extracted = sorted(set(origin_of) - kept_ids)
if not_extracted:
    log.write(
        f"[metaxa2_extract] WARNING: {len(not_extracted)} ASV(s) have an origin call but "
        f"are not in {extraction_fasta.name} (first few: {not_extracted[:5]})\n"
    )

origin_out.parent.mkdir(parents=True, exist_ok=True)
n_per_origin = {origin: 0 for origin in ORIGINS}
n_no_call = 0
with origin_out.open("w") as fh:
    fh.write("asv_id\torigin\n")
    for asv_id in kept_order:
        origin = origin_of.get(asv_id, "NA")
        if origin == "NA":
            n_no_call += 1
        else:
            n_per_origin[origin] += 1
        fh.write(f"{asv_id}\t{origin}\n")
summary = ", ".join(f"{origin} {n}" for origin, n in n_per_origin.items())
log.write(f"[metaxa2_extract] origin calls of the {len(kept_order)} extracted ASV(s): "
          f"{summary}, no call {n_no_call} -> {origin_out.name}\n")

# Keep all Metaxa2 raw outputs (seqs.summary.txt, seqs.hmmer.table, seqs.graph,
# seqs.extraction.fasta, seqs.extraction.fasta.1, seqs.extraction.results) under
# the prefix dir for QC / downstream inspection. Final dir layout:
#   <prefix_dir>/seqs.*  ← raw Metaxa2 working files (persisted)
#   5.dada2/seqs_extracted.fasta             ← promoted, cleaned headers
#   5.dada2/metaxa2_extraction.results.txt   ← promoted .extraction.results
#   5.dada2/metaxa2_origin.tsv               ← origin call per ASV, from seqs.<origin>.fasta
log.write(f"[metaxa2_extract] Keeping raw outputs under: {Path(prefix).parent}\n")
log.write(f"[metaxa2_extract] DONE.\n")
log.close()
