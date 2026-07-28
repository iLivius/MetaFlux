#!/usr/bin/env python3
"""ITSx: trim ASVs down to the ITS sub-region actually used for fungal ID.

This is the ITS "target_extract" step (see 70_extract.smk), run only when
amplicon.extraction.enabled is true. ITS (Internal Transcribed Spacer) is the
standard fungal barcode marker: it sits between two conserved ribosomal genes
in the arrangement SSU-ITS1-5.8S-ITS2-LSU, and unlike 16S it varies a lot in
length between species, which is exactly why the ASV length filter downstream
(60d_dada_length_filter.py) needs a wide window rather than one fixed cut.
ITSx scans each input sequence against a library of HMM profiles — precomputed
statistical models of what each conserved flanking gene (SSU, 5.8S, LSU) looks
like, one profile per reference lineage — to find where those genes end and
the variable ITS1/ITS2 region begins, then extracts just that region.

Input: seqs.fasta and seqtab_head_names.txt, straight from dada_seqtab.R — the
full, pre-length-filter set of ASVs (amplicon sequence variants — DADA2's
exact-sequence output unit; see rule target_extract in 70_extract.smk).

Run with --only_full F so that partial detections are kept — mandatory for
targeted ITS1/ITS2 amplicon data, where a single short read only ever contains
ONE flanking gene (e.g. 5.8S but not LSU for ITS2), never both, so requiring a
"full" detection (both flanks present) would discard every real read. -t all
scans all eukaryote HMM profiles, not just fungi-specific ones — so a non-
fungal eukaryote that happened to co-amplify (e.g. plant or protist ITS) is
still recognised and extracted here rather than silently dropped as "not
found"; excluding it is left to the taxonomy keep/discard filter downstream
(amplicon.taxonomy.filter, typically keep: [k__Fungi] for this marker).

ITSx-output dedup: distinct DADA2 ASVs can yield identical extracted ITS
sequences when they differ only in flanking 5.8S/LSU bases (e.g. multi-primer
or multiplex amplicon designs). Such ASVs are collapsed into a single
representative (the lowest-numbered ASV ID per group, which DADA2 orders by
abundance) and their per-sample counts are summed in the seqtab.

Output: seqs_extracted.fasta (deduplicated extracted ITS1/ITS2 sequences),
seqtab_extracted_head_names.txt (matching, summed read counts), and
itsx_collapse_map.tsv (an audit trail of which original ASVs each output ASV
represents) are read next by dada_length_filter (60d_dada_length_filter.py).
ITSx's own raw output files are kept under the working prefix directory for QC.
"""
import shutil
import subprocess
import sys
from pathlib import Path


# ── Sort ASV IDs by their DADA2 abundance order (lower number = more abundant) ─
def asv_sort_key(asv_id: str) -> int:
    try:
        return int(asv_id.split("_")[-1])
    except ValueError:
        return 10**9  # ASVs without a numeric suffix sort last


# ── Sum seqtab columns that share an extracted sequence ─────────
def collapse_seqtab_by_groups(in_path: Path, collapse_map: dict, reps_in_order: list,
                              out_path: Path, log) -> None:
    """Write a seqtab where ASV columns sharing a sequence are summed.

    collapse_map: {original_asv_id → representative_asv_id}
    reps_in_order: representative ASV IDs in DADA2 abundance order
    ASV columns not present in collapse_map (e.g. not extracted by ITSx) are dropped.
    """
    with in_path.open() as fh:
        raw_header = fh.readline().rstrip("\n")
        data_lines = fh.readlines()
    header = raw_header.split("\t")
    header_clean = [h.strip('"') for h in header]

    sample_col = header[0]
    asv_cols = header_clean[1:]

    rep_to_source_idx: dict[str, list[int]] = {}
    for i, asv in enumerate(asv_cols, start=1):
        rep = collapse_map.get(asv)
        if rep is None:
            continue
        rep_to_source_idx.setdefault(rep, []).append(i)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        out_header = [sample_col] + [f'"{rep}"' for rep in reps_in_order]
        fh.write("\t".join(out_header) + "\n")
        for line in data_lines:
            cols = line.rstrip("\n").split("\t")
            out_row = [cols[0]]
            for rep in reps_in_order:
                total = sum(int(cols[i]) for i in rep_to_source_idx[rep])
                out_row.append(str(total))
            fh.write("\t".join(out_row) + "\n")

    n_merged_groups = sum(1 for rep in reps_in_order if len(rep_to_source_idx[rep]) > 1)
    n_dropped = len(asv_cols) - len(collapse_map)
    log.write(
        f"[itsx_extract] seqtab: {len(reps_in_order)} representative ASV column(s); "
        f"{n_merged_groups} group(s) summed >1 ASV; "
        f"{n_dropped} ASV column(s) dropped (not extracted by ITSx)\n"
    )


# ── Inputs, parameters, and log file ────────────────────────────
sm = snakemake  # noqa: F821

input_fasta       = str(sm.input.seqs)
seqtab_names_in   = Path(sm.input.seqtab_names)
output_fasta      = Path(sm.output.seqs)
seqtab_names_out  = Path(sm.output.seqtab_names)
results_out       = Path(sm.output.results)
collapse_map_out  = Path(sm.output.collapse_map)
its_region        = str(sm.params.its_region)   # "ITS1" or "ITS2"
prefix            = str(sm.params.prefix)
threads           = int(sm.threads)

log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log = log_path.open("w")

Path(prefix).parent.mkdir(parents=True, exist_ok=True)

# ── Run ITSx ────────────────────────────────────────────────────
cmd = [
    "ITSx",
    "-i", input_fasta,
    "-o", prefix,
    "--only_full", "F",
    "-t", "all",
    "--cpu", str(threads),
]
log.write(f"[itsx_extract] Running: {' '.join(cmd)}\n")
log.write(f"[itsx_extract] Extracting region: {its_region}\n")
log.flush()
res = subprocess.run(cmd, stdout=log, stderr=log)
if res.returncode != 0:
    log.write(f"[itsx_extract] ERROR: ITSx exited with code {res.returncode}\n")
    log.close()
    sys.exit(res.returncode)

region_fasta = Path(f"{prefix}.{its_region}.fasta")
if not region_fasta.exists() or region_fasta.stat().st_size == 0:
    log.write(
        f"[itsx_extract] ERROR: {region_fasta} missing or empty. "
        f"Check that your amplicon targets {its_region} and primers are correct.\n"
    )
    log.close()
    sys.exit(1)

# ── Parse ITSx region FASTA → {asv_id: extracted_sequence} ──────
# ITSx headers: ">ASV_1|F|ITS2 Extracted ..." → take first pipe-delimited token
extracted: dict[str, str] = {}
with region_fasta.open() as fi:
    current_id: str | None = None
    current_seq: list[str] = []
    for line in fi:
        if line.startswith(">"):
            if current_id is not None:
                extracted[current_id] = "".join(current_seq)
            raw_id = line[1:].strip()
            current_id = raw_id.split("|")[0].split()[0]
            current_seq = []
        else:
            current_seq.append(line.rstrip("\n"))
    if current_id is not None:
        extracted[current_id] = "".join(current_seq)

log.write(f"[itsx_extract] {len(extracted)} ASV(s) extracted for {its_region}\n")

# ── Group ASVs by identical extracted sequence ──────────────────
# Multi-primer / multiplex amplicon designs often produce distinct DADA2 ASVs
# whose only differences are in the flanking 5.8S/LSU regions. After ITSx trims
# those flanks away, the ASVs share an identical ITS sub-region and should be
# treated as the same marker. We pick the lowest-numbered ASV per group as the
# representative (DADA2 numbers ASVs in decreasing abundance order).
seq_to_ids: dict[str, list[str]] = {}
for asv_id, seq in extracted.items():
    seq_to_ids.setdefault(seq, []).append(asv_id)

seq_to_rep: dict[str, str] = {seq: min(ids, key=asv_sort_key)
                              for seq, ids in seq_to_ids.items()}
collapse_map: dict[str, str] = {asv_id: rep
                                for seq, ids in seq_to_ids.items()
                                for rep in [seq_to_rep[seq]]
                                for asv_id in ids}
reps_in_order: list[str] = sorted(seq_to_rep.values(), key=asv_sort_key)
rep_to_seq: dict[str, str] = {rep: seq for seq, rep in seq_to_rep.items()}

# ── Report grouping ─────────────────────────────────────────────
n_groups   = len(seq_to_ids)
n_collapsed = sum(1 for ids in seq_to_ids.values() if len(ids) > 1)
log.write(
    f"[itsx_extract] Deduplication: {len(extracted)} extracted ASVs → "
    f"{n_groups} unique sequence(s); {n_collapsed} group(s) collapse >1 ASV\n"
)
for ids in seq_to_ids.values():
    if len(ids) > 1:
        ids_sorted = sorted(ids, key=asv_sort_key)
        rep = ids_sorted[0]
        log.write(f"[itsx_extract]   {rep} <- {', '.join(ids_sorted)}\n")

# ── Write deduplicated FASTA (one entry per representative) ─────
output_fasta.parent.mkdir(parents=True, exist_ok=True)
with output_fasta.open("w") as fo:
    for rep in reps_in_order:
        fo.write(f">{rep}\n{rep_to_seq[rep]}\n")

# ── Write collapse map (audit trail: representative → merged ASVs) ─
# One row per representative ASV (singletons included) so downstream tools can
# always recover which originals each output ASV came from.
collapse_map_out.parent.mkdir(parents=True, exist_ok=True)
with collapse_map_out.open("w") as fo:
    fo.write("representative\tmerged_asvs\tn_merged\textracted_length\n")
    for rep in reps_in_order:
        seq = rep_to_seq[rep]
        members = sorted(seq_to_ids[seq], key=asv_sort_key)
        fo.write(f"{rep}\t{','.join(members)}\t{len(members)}\t{len(seq)}\n")
log.write(f"[itsx_extract] Collapse map written → {collapse_map_out}\n")

# ── Write deduplicated / summed seqtab ──────────────────────────
collapse_seqtab_by_groups(seqtab_names_in, collapse_map, reps_in_order,
                          seqtab_names_out, log)

# ── Copy ITSx's native summary file to declared output ──────────
# QC report, NOT consumed by any downstream rule. The real extraction output — the
# region FASTA — is already validated above, so a missing summary must not fail an
# otherwise-complete run; warn loudly and leave a self-explaining placeholder.
# (Same missing-summary policy as 70a_metaxa2_extract.py — tolerant but loud.)
summary_src = Path(f"{prefix}.summary.txt")
results_out.parent.mkdir(parents=True, exist_ok=True)
if summary_src.exists():
    shutil.copy2(summary_src, results_out)
    log.write(f"[itsx_extract] Copied {summary_src} → {results_out}\n")
else:
    log.write(
        f"[itsx_extract] WARNING: expected summary {summary_src} not found; "
        "the region FASTA was produced, so continuing. Writing a placeholder.\n"
    )
    results_out.write_text(
        f"# ITSx summary ({summary_src.name}) was not produced by this run.\n"
        "# Extraction itself succeeded (see the extracted FASTA); this QC report is absent.\n"
    )

log.write("[itsx_extract] DONE.\n")
log.close()
