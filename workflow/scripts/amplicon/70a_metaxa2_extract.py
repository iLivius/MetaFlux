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

Merging ASVs that come out identical: Metaxa2 usually returns a 16S ASV whole,
but now and then it trims one (on a V4 run of 144 field samples it cut 31 bases
off the 3' end of five ASVs). Two ASVs that differed only in the trimmed part
then come out with the same sequence. They are merged into one, as after ITSx
(70b_itsx_extract.py): the lowest-numbered ASV of the group represents it (DADA2
numbers ASVs by decreasing abundance, so it is the most abundant one) and the
group's per-sample counts are added up. Without the merge, the same sequence
would sit in two rows of every later table, as two ASVs of one organism.

Output: seqs_extracted.fasta (cleaned, trimmed sequences, one per ASV left after
merging) and seqtab_extracted_head_names.txt (the matching read-count table,
subset to the ASVs Metaxa2 could extract, with merged ASVs added up) are read
next by dada_length_filter (60h_dada_length_filter.py), which applies the ASV
length window. metaxa2_origin.tsv (one row per ASV left after merging: asv_id,
origin) is read by assign_taxonomy (see above). metaxa2_collapse_map.tsv (one row
per ASV left after merging, with the ASV IDs merged into it) is a record for the
user; no rule reads it. Metaxa2's own raw output files are kept under the working
prefix directory for QC.
"""
import subprocess
import sys
from pathlib import Path


def asv_number(asv_id: str) -> int:
    """The number of an ASV ID: 12 for ASV_12.

    DADA2 numbers its ASVs by decreasing abundance, so the lowest number of a group
    is its most abundant ASV, which becomes the group's representative. An ID
    without a number sorts last. The same rule as asv_sort_key in 70b_itsx_extract.py.
    """
    try:
        return int(asv_id.split("_")[-1])
    except ValueError:
        return 10**9


def merge_seqtab_columns(in_path: Path, representative_of: dict[str, str],
                         out_path: Path, log) -> None:
    """Write the read-count table of the extracted ASVs, merged ASVs added up.

    in_path is the full seqtab_head_names.txt from dada_seqtab.R: one row per
    sample, one column per ASV (every ASV DADA2 made), the header quoted by R.
    representative_of maps every ASV ID that Metaxa2 extracted to the ID of the
    ASV that stands for its sequence: itself, or the ASV it was merged into
    (built below from Metaxa2's output FASTA).

    Column by column, in the order of the input table:
      - an ASV that Metaxa2 did not extract is left out: Metaxa2 could not
        confirm it as SSU rRNA in any domain, so its reads drop out of every
        later table;
      - a representative keeps its column, in its place, with its own counts
        plus those of the ASVs merged into it;
      - an ASV merged into another one is left out: its counts are now in the
        representative's column.
    When nothing is merged, the result is the input without the columns of the
    ASVs Metaxa2 did not extract, written exactly as before the merge existed.
    """
    with in_path.open() as fh:
        raw_header = fh.readline().rstrip("\n")
        data_lines = fh.readlines()
    header = raw_header.split("\t")
    # R write.table wraps every column name in double quotes; strip them for matching.
    header_clean = [h.strip('"') for h in header]

    # The column positions of the members of each group, listed under the
    # representative's ID (a group of one is the representative alone).
    member_columns: dict[str, list[int]] = {}
    for i in range(1, len(header_clean)):
        rep = representative_of.get(header_clean[i])
        if rep is not None:
            member_columns.setdefault(rep, []).append(i)
    # The representatives' own columns, in the order of the input table.
    rep_columns = [i for i in range(1, len(header_clean))
                   if representative_of.get(header_clean[i]) == header_clean[i]]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("\t".join([header[0]] + [header[i] for i in rep_columns]) + "\n")
        for line in data_lines:
            cols = line.rstrip("\n").split("\t")
            out_row = [cols[0]]                    # the sample name
            for i in rep_columns:
                members = member_columns[header_clean[i]]
                if len(members) == 1:
                    out_row.append(cols[i])        # a group of one: the count as written
                else:
                    out_row.append(str(sum(int(cols[j]) for j in members)))
            fh.write("\t".join(out_row) + "\n")

    n_columns = len(header_clean) - 1
    n_extracted = sum(len(members) for members in member_columns.values())
    log.write(f"[metaxa2_extract] seqtab: kept {len(rep_columns)} / {n_columns} ASV columns "
              f"({n_columns - n_extracted} not extracted by Metaxa2, "
              f"{n_extracted - len(rep_columns)} merged into another ASV)\n")


sm = snakemake  # noqa: F821

input_fasta      = str(sm.input.seqs)
seqtab_names_in  = Path(sm.input.seqtab_names)
output_fasta     = Path(sm.output.seqs)
seqtab_names_out = Path(sm.output.seqtab_names)
results_out      = Path(sm.output.results)
origin_out       = Path(sm.output.origin)        # metaxa2_origin.tsv, see the end of this script
collapse_map_out = Path(sm.output.collapse_map)  # metaxa2_collapse_map.tsv, see the merge below
prefix           = str(sm.params.prefix)
threads          = int(sm.threads)

log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log = log_path.open("w")

# Metaxa2 runs inside the folder of its own outputs (prefix_dir), with only the last
# part of the prefix ("seqs") as -o. Why: Metaxa2 names its temporary folder
# "metaxa_temp_directory__" followed by the -o value with every "/" turned into "-",
# and creates it as ONE folder name in the folder it is started from. With the full
# prefix as -o, that name grew with out_dir, and once it passed 255 bytes, the most
# Linux file systems allow for one name (an out_dir longer than 210 bytes), the folder
# could not be made: Metaxa2 then extracted nothing, still exited 0, and the run
# failed later in dada_length_filter with a bare KeyError. Started in its output
# folder with -o "seqs", the name is always "metaxa_temp_directory__seqs", whatever
# out_dir is, and Metaxa2 deletes the folder when it ends. The output files are the
# same as before (<prefix>.extraction.fasta and the others, read below through the
# full prefix). The input path is made absolute so it still points to the file from
# inside prefix_dir.
prefix_dir  = Path(prefix).parent      # <out_dir>/05.dada2/metaxa2
prefix_name = Path(prefix).name        # "seqs"
prefix_dir.mkdir(parents=True, exist_ok=True)

cmd = [
    "metaxa2",
    "-i", str(Path(input_fasta).absolute()),
    "-o", prefix_name,
    "-t", "all",          # search all domain HMM profiles (bacteria/archaea/eukaryote)
    "--cpu", str(threads),
    "--fasta", "T",       # write the extracted-sequence FASTA (*.extraction.fasta)
    "--table", "T",       # write the per-sequence domain/taxonomy table
    "--plus", "T",        # also detect chloroplast/mitochondrial SSU (organelle-derived,
                          # bacterial in origin) rather than treating them as unclassifiable
    "--silent", "T",      # suppress Metaxa2's own console output (everything goes to our log)
]
log.write(f"[metaxa2_extract] Running in {prefix_dir}: {' '.join(cmd)}\n")
log.flush()
res = subprocess.run(cmd, stdout=log, stderr=log, cwd=prefix_dir)
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

# Read Metaxa2's extraction FASTA: one record per extracted ASV. Metaxa2 writes each
# header as the ASV ID immediately followed by "|<rRNA subunit tag>", then a space
# and a free-text description, e.g.
# ">ASV_12|B Predicted Bacterial  rRNA (253 bp) From domain V4l to V4r on main strand".
# Taking the first whitespace token gives "ASV_12|B"; splitting that on "|"
# isolates the plain ASV ID — this is what has to match the seqtab column
# headers below. `extracted` keeps the records in Metaxa2's order (a Python dict
# keeps the order in which keys were added), and that order is kept in the FASTA
# written below. An ID listed twice should not happen; if it does, the first record
# is kept and the log says so.
extracted: dict[str, str] = {}      # ASV ID -> extracted sequence
repeated_ids: list[str] = []


def add_record(asv_id: str | None, seq_lines: list[str]) -> None:
    """Store one FASTA record in `extracted`, unless its ID is there already."""
    if asv_id is None:
        return
    if asv_id in extracted:
        repeated_ids.append(asv_id)
        return
    extracted[asv_id] = "".join(seq_lines)


with extraction_fasta.open() as fi:
    current_id: str | None = None
    current_seq: list[str] = []
    for line in fi:
        if line.startswith(">"):
            add_record(current_id, current_seq)
            current_id = line[1:].split()[0].split("|")[0]
            current_seq = []
        else:
            current_seq.append(line.rstrip("\n"))
    add_record(current_id, current_seq)

log.write(f"[metaxa2_extract] {len(extracted)} ASV(s) extracted\n")
if repeated_ids:
    log.write(f"[metaxa2_extract] WARNING: {len(repeated_ids)} ASV ID(s) listed more than once "
              f"in {extraction_fasta.name}; the first record of each was kept "
              f"(first few: {repeated_ids[:5]})\n")

# ── Merge ASVs whose extracted sequences are identical ──────────────────────
# Group the extracted ASVs by sequence. Each group is represented by its
# lowest-numbered ASV, the most abundant one (asv_number, above), as after ITSx.
# Most groups hold one ASV: Metaxa2 seldom trims a 16S ASV, and two ASVs only
# come out identical when it trimmed away the only bases in which they differed.
ids_of_sequence: dict[str, list[str]] = {}
for asv_id, seq in extracted.items():
    ids_of_sequence.setdefault(seq, []).append(asv_id)

representative_of: dict[str, str] = {}      # every extracted ASV -> its group's representative
for ids in ids_of_sequence.values():
    rep = min(ids, key=asv_number)
    for asv_id in ids:
        representative_of[asv_id] = rep

# Every group as a list of IDs sorted by number, so its representative comes first;
# the groups in the order of their representatives' numbers. merged_groups are the
# groups of two or more ASVs.
groups = [sorted(ids, key=asv_number) for ids in ids_of_sequence.values()]
groups.sort(key=lambda group: asv_number(group[0]))
merged_groups = [group for group in groups if len(group) > 1]
log.write(f"[metaxa2_extract] merging: {len(extracted)} extracted ASV(s) -> "
          f"{len(ids_of_sequence)} distinct sequence(s); {len(merged_groups)} group(s) of "
          f"two or more ASVs merged\n")
for group in merged_groups:
    log.write(f"[metaxa2_extract]   {group[0]} <- {', '.join(group)} "
              f"({len(extracted[group[0]])} bp after extraction)\n")

# The cleaned FASTA: one record per representative, in Metaxa2's order, with the
# plain ASV ID as header. When nothing is merged, it holds every extracted record,
# exactly as before the merge existed.
output_fasta.parent.mkdir(parents=True, exist_ok=True)
with output_fasta.open("w") as fo:
    for asv_id, seq in extracted.items():
        if representative_of[asv_id] == asv_id:
            fo.write(f">{asv_id}\n{seq}\n")

# The collapse map, in the layout of ITSx's itsx_collapse_map.tsv: one row per
# representative, groups of one included, by ASV number, with the IDs merged into
# it (itself first), how many, and the length of the extracted sequence. A record
# for the user (which ASVs became one after extraction); no rule reads it.
collapse_map_out.parent.mkdir(parents=True, exist_ok=True)
with collapse_map_out.open("w") as fo:
    fo.write("representative\tmerged_asvs\tn_merged\textracted_length\n")
    for group in groups:
        fo.write(f"{group[0]}\t{','.join(group)}\t{len(group)}\t{len(extracted[group[0]])}\n")
log.write(f"[metaxa2_extract] collapse map written -> {collapse_map_out.name}\n")

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

# The read-count table: the extracted ASVs only, merged ASVs added up.
merge_seqtab_columns(seqtab_names_in, representative_of, seqtab_names_out, log)

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
# Output: one row per ASV left after merging (the representatives), in the order of
# Metaxa2's extraction FASTA (a lookup table: nothing depends on the row order),
# columns asv_id and origin (mitochondria, chloroplast, bacteria, archaea,
# eukaryota, uncertain; NA when Metaxa2 extracted the sequence but listed it in none
# of the files). A merged group takes the call of its representative. The ASVs of a
# group share one sequence, so Metaxa2 should give them all the same call; if it
# does not, the log says so. Read by rule assign_taxonomy: the organelle filter
# (amplicon.taxonomy.filter.organelles) removes the mitochondria and chloroplast
# ASVs, and every ASV's call goes into the metaxa2_origin column of
# 06.taxonomy/asv_table_unfiltered.txt.
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
not_extracted = sorted(set(origin_of) - set(extracted))
if not_extracted:
    log.write(
        f"[metaxa2_extract] WARNING: {len(not_extracted)} ASV(s) have an origin call but "
        f"are not in {extraction_fasta.name} (first few: {not_extracted[:5]})\n"
    )

# A merged group whose ASVs got different calls: keep the representative's, say so.
for group in merged_groups:
    calls = [origin_of.get(asv_id, "NA") for asv_id in group]
    if len(set(calls)) > 1:
        listed = ", ".join(f"{asv_id} {call}" for asv_id, call in zip(group, calls))
        log.write(f"[metaxa2_extract] WARNING: the ASVs merged into {group[0]} have different "
                  f"origin calls ({listed}); {group[0]} keeps its own, {calls[0]}\n")

origin_out.parent.mkdir(parents=True, exist_ok=True)
n_per_origin = {origin: 0 for origin in ORIGINS}
n_no_call = 0
n_rows = 0
with origin_out.open("w") as fh:
    fh.write("asv_id\torigin\n")
    for asv_id in extracted:
        if representative_of[asv_id] != asv_id:
            continue                    # merged into another ASV: no row of its own
        origin = origin_of.get(asv_id, "NA")
        if origin == "NA":
            n_no_call += 1
        else:
            n_per_origin[origin] += 1
        fh.write(f"{asv_id}\t{origin}\n")
        n_rows += 1
summary = ", ".join(f"{origin} {n}" for origin, n in n_per_origin.items())
log.write(f"[metaxa2_extract] origin calls of the {n_rows} ASV(s) left after merging: "
          f"{summary}, no call {n_no_call} -> {origin_out.name}\n")

# Keep all Metaxa2 raw outputs (seqs.summary.txt, seqs.hmmer.table, seqs.graph,
# seqs.extraction.fasta, seqs.extraction.fasta.1, seqs.extraction.results) under
# the prefix dir for QC / downstream inspection. Final dir layout:
#   <prefix_dir>/seqs.*  ← raw Metaxa2 working files (persisted)
#   05.dada2/seqs_extracted.fasta            ← promoted, cleaned headers, merged ASVs once
#   05.dada2/metaxa2_extraction.results.txt  ← promoted .extraction.results
#   05.dada2/metaxa2_origin.tsv              ← origin call per ASV, from seqs.<origin>.fasta
#   05.dada2/metaxa2_collapse_map.tsv        ← which ASVs were merged into which
log.write(f"[metaxa2_extract] Keeping raw outputs under: {Path(prefix).parent}\n")
log.write(f"[metaxa2_extract] DONE.\n")
log.close()
