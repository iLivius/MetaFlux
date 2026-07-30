#!/usr/bin/env python3
"""Assign a taxonomic lineage to each ASV using VSEARCH's SINTAX classifier.

This is one of the two implementations of the assign_taxonomy rule (see rule
assign_taxonomy in 80_taxonomy.smk), selected when config amplicon.taxonomy.
method is "sintax" (the default) rather than "rdp". SINTAX is a k-mer-based
classifier: instead of aligning each ASV to reference sequences, it compares
short fixed-length substrings (k-mers) of the ASV against those of every
reference sequence and picks the best-matching lineage, which is much cheaper
than alignment and scales well to large reference databases. This is a
different algorithm from the RDP naive Bayesian classifier the other path
(80a_assign_taxonomy.R) uses, but this script is written to produce byte-
identical table STRUCTURE, so every script downstream of assign_taxonomy is
agnostic to which classifier actually ran.

Input: seqs (the length-filtered ASV FASTA, seqs_lenfilt.fasta) and
seqtab_names (its matching read-count table, seqtab_lenfilt_head_names.txt) —
both produced by dada_length_filter (60d_dada_length_filter.py), i.e. this is
the FINAL ASV set, nothing upstream of this step still changes which ASVs
exist. refdb is the marker's SINTAX-formatted reference database (e.g. SILVA
for 16S, UNITE for ITS — see taxonomy_sintax_db in workflow/markers/*.yaml).

Mechanism: runs `vsearch --sintax`, which reports, for each ASV, a lineage at
every rank together with a bootstrap confidence per rank — bootstrap here means
vsearch repeatedly resamples the ASV's k-mers and re-classifies, and the
confidence is the fraction of those resamples that agreed on that rank's call.
Only ranks whose confidence clears --sintax_cutoff (config
amplicon.taxonomy.sintax_cutoff) are kept; deeper ranks are left blank rather
than reported unreliably.

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
These three files are final pipeline deliverables (requested by `rule all`);
the only thing that reads one of them back in is aggregate_read_counts.py,
which sums asv_table.txt's sample columns to report the post-taxonomy-filter
read count per sample.

All outputs use R write.table quoting conventions (character values double-quoted,
integer counts unquoted) so downstream scripts are agnostic to the classifier.
"""
import gzip
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

# Rank model from the marker profile (see 00_common.smk MARKERS), replacing the
# previously hardcoded 7-rank tuples. The SINTAX path always normalises values
# to bare letters then re-adds the prefixes, so it needs no prefix_style.
RANK_LETTERS = tuple(sm.params.rank_letters)
RANK_NAMES   = tuple(sm.params.tax_levels)
RANK_PREFIX  = tuple(sm.params.rank_prefixes)

# Contaminant keep/discard: rank-token lists from config, no per-marker default.
filter_keep    = list(sm.params.filter_keep) or []
filter_discard = list(sm.params.filter_discard) or []

log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log = log_path.open("w")


def logmsg(msg: str) -> None:
    log.write(f"[assign_taxonomy/sintax] {msg}\n")
    log.flush()


def _strip_existing_prefix(val: str) -> str:
    """Remove any leading rank prefix that UNITE sometimes embeds in values.

    Some reference databases (UNITE for ITS) store the rank prefix as part of
    the taxon name itself, e.g. the genus field literally reads "g__Fusarium"
    rather than "Fusarium". If we didn't strip that off first, build_tax_string
    (below) would add its own prefix on top and produce "g__g__Fusarium". This
    runs unconditionally on every value regardless of marker, so a reference
    that does NOT embed prefixes is simply unaffected (nothing to strip).
    """
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


logmsg(f"amp_type={amp_type}")
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

# Subset to ASVs present in the FASTA (length filter / extraction may have removed some)
shared_ids = [aid for aid in asv_ids if aid in set(seqtab_asv_ids)]
logmsg(f"Seqtab: {len(shared_ids)} ASVs shared with FASTA")

# ── Run VSEARCH SINTAX ─────────────────────────────────────────────────────
out_asv.parent.mkdir(parents=True, exist_ok=True)
tabbedout = out_asv.parent / "_sintax_raw.tsv"

# The SINTAX DB (SILVA for 16S, PR2/UTAX for 18S, UNITE for ITS) is .gz; VSEARCH
# reads gzipped FASTA natively.
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
# Every line --tabbedout writes has exactly 4 columns, even for a query with no
# hit at all (columns 2-4 are then just empty rather than absent) — verified
# directly against vsearch, including a forced no-hit case. A line with a
# different shape means something upstream changed (a vsearch version with a
# different --tabbedout format, a corrupted file, ...), not a case to guess
# through silently, so it fails here with the line itself in the message.
sintax_result: dict[str, dict[str, str]] = {}
with tabbedout.open() as fh:
    for lineno, line in enumerate(fh, start=1):
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 4:
            sys.exit(
                f"[assign_taxonomy] {tabbedout}:{lineno}: expected 4 tab-separated "
                f"columns from vsearch --tabbedout, got {len(parts)}: {line!r}"
            )
        query_id, accepted_col = parts[0], parts[3]
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

# ── Contaminant filter (rank-aware keep/discard lists; config-only) ─────────
# Each token (e.g. k__Bacteria, o__Chloroplast) matches whole ';'-delimited
# segments of the taxonomy string — so g__Bacillus hits only a genus named
# exactly that, never a substring elsewhere. keep runs first (defines the
# universe), then discard prunes. Empty list = that direction is a no-op.
def _seg_match(tax_string: str, tokens: list[str]) -> bool:
    segs = {s.strip() for s in tax_string.split(";")}
    return any(tok in segs for tok in tokens)


filtered_ids = list(shared_ids)
n_before = len(filtered_ids)

if filter_enabled:
    if filter_keep:
        filtered_ids = [aid for aid in filtered_ids
                        if _seg_match(tax_strings[aid], filter_keep)]
        logmsg(f"Keep filter ({','.join(filter_keep)}): {len(filtered_ids)} / {n_before} ASVs retained")
    if filter_discard:
        n_pre_excl   = len(filtered_ids)
        filtered_ids = [aid for aid in filtered_ids
                        if not _seg_match(tax_strings[aid], filter_discard)]
        logmsg(f"Discard filter ({','.join(filter_discard)}): {n_pre_excl - len(filtered_ids)} ASV(s) removed")

logmsg(f"Final ASV count: {len(filtered_ids)}")

# Fail loud rather than write a header-only table: a filter that removes 100% of
# ASVs is almost always a marker/config mismatch (e.g. 16S keep tokens on an ITS run).
if filter_enabled and n_before > 0 and len(filtered_ids) == 0:
    msg = (f"the keep/discard filter removed ALL {n_before} ASVs (amp_type={amp_type}). "
           "The keep/discard lists likely do not match this marker: for ITS use keep: [k__Fungi]; "
           "for 16S use keep: [k__Bacteria, k__Archaea]. Fix amplicon.taxonomy.filter, "
           "or set amplicon.taxonomy.filter.enabled: false.")
    logmsg("ERROR: " + msg)
    log.close()
    sys.exit(f"[assign_taxonomy/sintax] {msg}")


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
