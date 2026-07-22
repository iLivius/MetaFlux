#!/usr/bin/env python3
"""Parse vsearch --sintax tabbedout output into the same three table formats
produced by 80a_assign_taxonomy.R (RDP path).

VSEARCH --tabbedout columns:
  1  query ID
  2  full taxonomy with per-rank confidences: d:Bacteria(0.98),p:Proteobacteria(0.95),...
  3  strand (+/-)
  4  accepted taxonomy (ranks meeting --sintax_cutoff): d:Bacteria,p:Proteobacteria,...
     Empty when unclassified.

Output files (identical structure to the RDP path):
  asv_table.txt       — rows=ASV_IDs,   cols=[samples..., taxonomy]
  asv_table_seqs.txt  — rows=sequences, cols=[samples..., taxonomy]
  taxon_seq_table.txt — rows=ASV_IDs,   cols=[Kingdom..Species, sequence]

All outputs use R write.table quoting conventions (character values double-quoted,
integer counts unquoted) so downstream scripts are agnostic to the classifier.
"""
import gzip
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sm = snakemake  # noqa: F821

seqs_path     = Path(sm.input.seqs)
seqtab_path   = Path(sm.input.seqtab_names)
refdb_path    = Path(sm.input.refdb)
out_asv       = Path(sm.output.asv_table)
out_seqs      = Path(sm.output.asv_table_seqs)
out_taxon     = Path(sm.output.taxon_seq_table)

amp_type        = sm.params.amp_type
sintax_cutoff   = float(sm.params.sintax_cutoff)
filter_enabled  = bool(sm.params.filter_enabled)
threads         = int(sm.threads)
seed            = int(sm.params.seed)

# YAML null arrives as Python None or the string "null" depending on Snakemake version.
def _none_or_str(v) -> str | None:
    return None if (v is None or str(v).strip().lower() in ("null", "none", "")) else str(v)

include_pattern = _none_or_str(sm.params.include_pattern)
exclude_pattern = _none_or_str(sm.params.exclude_pattern)

if include_pattern is None:
    include_pattern = "k__Archaea|k__Bacteria" if amp_type == "16S" else "k__Fungi"
if exclude_pattern is None:
    exclude_pattern = "chloroplast|mitochondria" if amp_type == "16S" else None

log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log = log_path.open("w")


def logmsg(msg: str) -> None:
    log.write(f"[assign_taxonomy/sintax] {msg}\n")
    log.flush()


# ── Rank tables ────────────────────────────────────────────────────────────
RANK_LETTERS = ("d", "p", "c", "o", "f", "g", "s")
RANK_NAMES   = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
RANK_PREFIX  = ("k__", "p__", "c__", "o__", "f__", "g__", "s__")


def _strip_existing_prefix(val: str) -> str:
    """Remove any leading rank prefix that UNITE sometimes embeds in values."""
    for pfx in RANK_PREFIX:
        if val.startswith(pfx):
            return val[len(pfx):]
    return val


def parse_sintax_column(field: str) -> dict[str, str]:
    """Parse VSEARCH tabbedout column 4 (accepted taxonomy at cutoff).
    Format: 'd:Bacteria,p:Proteobacteria,...'  — no confidence values here.
    Returns {rank_letter: value} for each present rank."""
    result: dict[str, str] = {}
    if not field or field in ("*", ""):
        return result
    for token in field.split(","):
        token = token.strip()
        if ":" not in token:
            continue
        letter, value = token.split(":", 1)
        letter = letter.strip()
        value  = _strip_existing_prefix(value.strip().strip('"'))
        if value:
            result[letter] = value
    return result


def build_tax_string(ranks: dict[str, str]) -> str:
    """Build 'k__X;p__Y;...' taxonomy string from {letter: value} dict."""
    parts = []
    for letter, prefix in zip(RANK_LETTERS, RANK_PREFIX):
        if letter in ranks:
            parts.append(prefix + ranks[letter])
    return ";".join(parts)


# ── Read FASTA ─────────────────────────────────────────────────────────────
def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current_id: str | None = None
    buf: list[str] = []
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(buf)
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(buf)
    return seqs


fasta = read_fasta(seqs_path)
asv_ids = list(fasta.keys())
logmsg(f"FASTA: {len(fasta)} ASVs loaded from {seqs_path}")

# ── Read seqtab_head_names ─────────────────────────────────────────────────
with seqtab_path.open() as fh:
    header_line = fh.readline().rstrip("\n")
    data_lines  = fh.readlines()

seqtab_cols    = header_line.split("\t")
seqtab_asv_ids = [c.strip('"') for c in seqtab_cols[1:]]
samples: list[str] = []
counts: dict[str, list[str]] = {asv: [] for asv in seqtab_asv_ids}

for line in data_lines:
    cols = line.rstrip("\n").split("\t")
    samples.append(cols[0].strip('"'))
    for i, asv in enumerate(seqtab_asv_ids, 1):
        counts[asv].append(cols[i] if i < len(cols) else "0")

shared_ids = [aid for aid in asv_ids if aid in set(seqtab_asv_ids)]
logmsg(f"Seqtab: {len(shared_ids)} ASVs shared with FASTA")

# ── Run VSEARCH SINTAX ─────────────────────────────────────────────────────
out_asv.parent.mkdir(parents=True, exist_ok=True)
tabbedout = out_asv.parent / "_sintax_raw.tsv"

# The UNITE SINTAX DB is .gz; VSEARCH reads gzipped FASTA natively.
# --randseed: vsearch's sintax bootstrap confidence defaults to seed 0, i.e.
# "use a random data source". Pinning it removes that as a source of drift,
# but does NOT make multithreaded --sintax fully reproducible: vsearch races
# threads on one RNG stream, so per-rank confidence near sintax_cutoff can
# still shift a little between runs at threads > 1 (verified empirically).
# Full byte-reproducibility requires threads == 1 — see README Troubleshooting.
cmd = [
    "vsearch", "--sintax", str(seqs_path),
    "--db",            str(refdb_path),
    "--sintax_cutoff", str(sintax_cutoff),
    "--tabbedout",     str(tabbedout),
    "--strand",        "both",
    "--threads",       str(threads),
    "--randseed",      str(seed),
]
logmsg(f"Running: {' '.join(cmd)}")
res = subprocess.run(cmd, stdout=log, stderr=log)
if res.returncode != 0:
    logmsg(f"ERROR: vsearch exited with code {res.returncode}")
    log.close()
    sys.exit(res.returncode)

# ── Parse tabbedout ────────────────────────────────────────────────────────
sintax_result: dict[str, dict[str, str]] = {}
with tabbedout.open() as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if not parts:
            continue
        query_id      = parts[0]
        accepted_col  = parts[3] if len(parts) >= 4 else (parts[1] if len(parts) >= 2 else "")
        sintax_result[query_id] = parse_sintax_column(accepted_col)

logmsg(f"SINTAX: {len(sintax_result)} query results parsed")

# ── Build per-ASV taxonomy structures ─────────────────────────────────────
rank_rows:   dict[str, dict[str, str]] = {}
tax_strings: dict[str, str]            = {}

for asv_id in shared_ids:
    ranks = sintax_result.get(asv_id, {})
    tax_strings[asv_id] = build_tax_string(ranks)
    rank_rows[asv_id] = {
        col: ranks.get(letter, "")
        for letter, col in zip(RANK_LETTERS, RANK_NAMES)
    }

# ── Contaminant filter ─────────────────────────────────────────────────────
filtered_ids = list(shared_ids)
n_before = len(filtered_ids)

if filter_enabled:
    if include_pattern:
        filtered_ids = [aid for aid in filtered_ids
                        if re.search(include_pattern, tax_strings[aid], re.IGNORECASE)]
        logmsg(f"Include filter ({include_pattern}): {len(filtered_ids)} / {n_before} ASVs retained")
    if exclude_pattern:
        n_pre_excl   = len(filtered_ids)
        filtered_ids = [aid for aid in filtered_ids
                        if not re.search(exclude_pattern, tax_strings[aid], re.IGNORECASE)]
        logmsg(f"Exclude filter ({exclude_pattern}): {n_pre_excl - len(filtered_ids)} ASV(s) removed")

logmsg(f"Final ASV count: {len(filtered_ids)}")


# ── Write helpers — R write.table quoting conventions ─────────────────────
def q(s: str) -> str:
    """Wrap string in double quotes (R character column quoting)."""
    return f'"{s}"'


def write_count_table(path: Path, row_ids: list[str], row_names: list[str]) -> None:
    """Write asv_table.txt or asv_table_seqs.txt.
    Counts are unquoted integers; taxonomy is a quoted string."""
    with path.open("w") as fh:
        fh.write("\t".join([""] + [q(s) for s in samples] + [q("taxonomy")]) + "\n")
        for asv_id, row_name in zip(row_ids, row_names):
            vals = counts.get(asv_id, ["0"] * len(samples))
            fh.write("\t".join([q(row_name)] + vals + [q(tax_strings[asv_id])]) + "\n")


def write_taxon_table(path: Path, row_ids: list[str]) -> None:
    """Write taxon_seq_table.txt — all rank + sequence values are quoted strings."""
    with path.open("w") as fh:
        fh.write("\t".join([""] + [q(c) for c in RANK_NAMES] + [q("sequence")]) + "\n")
        for asv_id in row_ids:
            row  = rank_rows[asv_id]
            vals = [q(row.get(col, "")) for col in RANK_NAMES]
            fh.write("\t".join([q(asv_id)] + vals + [q(fasta[asv_id])]) + "\n")


# ── Output (a): asv_table.txt ──────────────────────────────────────────────
write_count_table(out_asv, filtered_ids, filtered_ids)

# ── Output (b): asv_table_seqs.txt ────────────────────────────────────────
seq_row_names = [fasta[aid] for aid in filtered_ids]
write_count_table(out_seqs, filtered_ids, seq_row_names)

# ── Output (c): taxon_seq_table.txt ───────────────────────────────────────
write_taxon_table(out_taxon, filtered_ids)

logmsg(f"DONE. Outputs written to {out_asv.parent}")
log.close()
