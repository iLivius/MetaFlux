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

For 16S, metaxa2_origin is also read: Metaxa2's own origin call per ASV
(5.dada2/metaxa2_origin.tsv, from target_extract — present only when Metaxa2
ran). With organelles on (params.organelle_filter, from
amplicon.taxonomy.filter.organelles) the ASVs it calls mitochondria or
chloroplast are removed together with keep/discard (see "Organelle filter").

Output files (identical structure to the RDP path):
  asv_table.txt            — rows=ASV_IDs,   cols=[samples..., taxonomy]
  asv_table_seqs.txt       — rows=sequences, cols=[samples..., taxonomy]
  taxon_seq_table.txt      — rows=ASV_IDs,   cols=[Kingdom..Species, sequence]
  asv_table_unfiltered.txt — rows=ASV_IDs, every ASV that reached the filter,
                             cols=[samples..., taxonomy, metaxa2_origin, filter_status]
The first three hold only the ASVs the filter kept. All four are final pipeline
deliverables (requested by `rule all`); aggregate_read_counts.py sums
asv_table.txt's sample columns to report the post-taxonomy-filter read count
per sample, and the phylogeny (85a_phylo_input.py) reads asv_table.txt and
taxon_seq_table.txt. Nothing reads asv_table_unfiltered.txt back in: it is for
the user's own analyses, e.g. comparing how many organelle reads two primer
sets pick up.

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
out_unfilt    = Path(sm.output.asv_table_unfiltered)
# Metaxa2's per-ASV origin table: a path when Metaxa2 ran in this run (16S with
# extraction on), absent (None) for every other marker — see rule assign_taxonomy.
origin_input  = sm.input.get("metaxa2_origin")
origin_path   = Path(origin_input) if origin_input else None

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
# Best-hit check of the discard filter (see "Discard filter, best-hit check" below):
# minimum support at the token's rank, as a percentage (0-100) like min_boot, or
# None when the check is switched off. SINTAX reports confidences as 0-1, so the
# value is divided by 100 when it is compared.
discard_min_support = sm.params.discard_min_support
# Organelle filter (see "Organelle filter" below): True when the ASVs Metaxa2 calls
# mitochondria or chloroplast are to be removed. Decided in 00_common.smk: it is
# only ever True for 16S with extraction on, filter.enabled true and
# filter.organelles: metaxa2.
organelle_filter = bool(sm.params.organelle_filter)

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


def parse_sintax_full(field: str) -> dict[str, tuple[str, float]]:
    """Parse VSEARCH tabbedout column 2 (the best hit's full lineage).

    Format: 'd:Bacteria(1.00),p:Pseudomonadota(0.97),...,g:Pseudomonas(0.41)' —
    every rank of the top-scoring reference, each with its bootstrap confidence,
    including ranks below sintax_cutoff. Returns {rank_letter: (value, confidence)}.
    The confidence is the last "(...)" of each token, because some SILVA names
    contain brackets themselves (e.g. "Marinimicrobia (SAR406 clade)").
    """
    result: dict[str, tuple[str, float]] = {}
    if not field or field in ("*", ""):
        return result
    for token in field.split(","):
        token = token.strip()
        if ":" not in token or not token.endswith(")") or "(" not in token:
            continue
        letter, rest = token.split(":", 1)
        value, conf = rest.rsplit("(", 1)
        try:
            confidence = float(conf[:-1])
        except ValueError:
            continue
        value = _strip_existing_prefix(value.strip().strip('"'))
        if value:
            result[letter.strip()] = (value, confidence)
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
# --randseed: without it vsearch uses seed 0, which means "take the seed from a
# random source", and the bootstrap confidence changes from run to run. With a
# fixed seed, VSEARCH 2.32 and later (pinned in envs/vsearch.yaml) give every ASV
# its own random stream, worked out from the seed and the ASV's place in the input
# FASTA, so the result is the same at any --threads value. Older releases shared
# one stream between all threads and were only reproducible on one thread. A seed
# of 0 in amplicon.seed would therefore switch reproducibility off for this step.
# See docs/troubleshooting.md for the reproducibility of the whole amplicon path.
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
# vsearch writes each ASV's line as soon as a thread has finished it, so above one
# thread the LINE ORDER of _sintax_raw.tsv changes from run to run (the lines
# themselves do not). Everything below is looked up by ASV id and written in the
# FASTA order (shared_ids), so the output tables do not depend on that order.
sintax_result: dict[str, dict[str, str]] = {}
sintax_best:   dict[str, dict[str, tuple[str, float]]] = {}
with tabbedout.open() as fh:
    for lineno, line in enumerate(fh, start=1):
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 4:
            sys.exit(
                f"[assign_taxonomy] {tabbedout}:{lineno}: expected 4 tab-separated "
                f"columns from vsearch --tabbedout, got {len(parts)}: {line!r}"
            )
        query_id, full_col, accepted_col = parts[0], parts[1], parts[3]
        sintax_result[query_id] = parse_sintax_column(accepted_col)
        sintax_best[query_id]   = parse_sintax_full(full_col)

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
#
# SILVA 144 "--other" names (amplicon.taxonomy.silva_other: verbatim) count as the
# plain name: the suffix is removed from both the segments and the tokens before
# comparing, so g__Bradyrhizobium also removes "g__Bradyrhizobium--other". A filter
# token is meant to name a taxon, and SILVA considers the "--other" clade part of
# that taxon. With silva_other: strip, or SILVA 138.2, no name carries the suffix
# and this changes nothing. Same rule as seg_match in 80a_assign_taxonomy.R.
def _drop_other(name: str) -> str:
    return name[: -len("--other")] if name.endswith("--other") else name


def _seg_match(tax_string: str, tokens: list[str]) -> bool:
    segs = {_drop_other(s.strip()) for s in tax_string.split(";")}
    return any(_drop_other(tok) in segs for tok in tokens)


# ── Metaxa2 origin per ASV (16S only; empty otherwise) ──
# metaxa2_origin.tsv has one row per ASV Metaxa2 extracted: asv_id, origin
# (mitochondria, chloroplast, bacteria, archaea, eukaryota, uncertain, or NA for
# no call). Used by the organelle filter below and written, for every ASV, to the
# metaxa2_origin column of asv_table_unfiltered.txt. An ASV with no entry (or for
# a marker where Metaxa2 does not run) has no origin, written as NA.
metaxa2_origin: dict[str, str] = {}
if origin_path is not None:
    with origin_path.open() as fh:
        fh.readline()                                   # header: asv_id, origin
        for line in fh:
            asv_id, origin = line.rstrip("\n").split("\t")
            if origin != "NA":
                metaxa2_origin[asv_id] = origin
    logmsg(f"Metaxa2 origin calls read for {len(metaxa2_origin)} ASVs from {origin_path}")
# The rule only switches the filter on when Metaxa2 ran, so the table is always
# there when it is needed; this guards against the two drifting apart.
if organelle_filter and origin_path is None:
    logmsg("ERROR: the organelle filter is on but no Metaxa2 origin table was given")
    log.close()
    sys.exit("[assign_taxonomy/sintax] the organelle filter is on but no Metaxa2 origin "
             "table was given (input metaxa2_origin of rule assign_taxonomy)")


def reads_of(asv_ids: list[str]) -> int:
    """Total reads, over all samples, of the given ASVs (from the seqtab counts)."""
    return sum(int(c) for aid in asv_ids for c in counts[aid])


filtered_ids = list(shared_ids)
n_before = len(filtered_ids)

# What the filter did to each ASV, for the filter_status column of
# asv_table_unfiltered.txt. The steps below run in a fixed order — keep, discard,
# best-hit check, organelle filter — and each one only looks at the ASVs the
# earlier steps kept. So a removed ASV is labelled with the FIRST step that
# removed it: removed_keep, removed_discard, removed_best_hit or
# removed_organelle. removed_organelle therefore means "only Metaxa2 caught it".
filter_status: dict[str, str] = {aid: "kept" for aid in shared_ids}

if filter_enabled:
    if filter_keep:
        for aid in filtered_ids:
            if not _seg_match(tax_strings[aid], filter_keep):
                filter_status[aid] = "removed_keep"
        filtered_ids = [aid for aid in filtered_ids if filter_status[aid] == "kept"]
        logmsg(f"Keep filter ({','.join(filter_keep)}): {len(filtered_ids)} / {n_before} ASVs retained")
    if filter_discard:
        n_pre_excl   = len(filtered_ids)
        for aid in filtered_ids:
            if _seg_match(tax_strings[aid], filter_discard):
                filter_status[aid] = "removed_discard"
        filtered_ids = [aid for aid in filtered_ids if filter_status[aid] == "kept"]
        logmsg(f"Discard filter ({','.join(filter_discard)}): {n_pre_excl - len(filtered_ids)} ASV(s) removed")

    # ── Discard filter, best-hit check ──
    # The accepted taxonomy stops at the last rank whose confidence reaches
    # sintax_cutoff. An ASV whose best hit is Mitochondria at family level, but
    # with a confidence below the cutoff there, is reported only down to class
    # and would slip past f__Mitochondria. So each discard token is also compared
    # with the best hit's name at the token's own rank (column 2 of the vsearch
    # output), and the ASV is removed when the names match and the confidence at
    # that rank is at least discard_min_support / 100. Same rule as the rdp path
    # (80a_assign_taxonomy.R); the keep list is not affected.
    if filter_discard and discard_min_support is not None:
        prefix_to_letter = dict(zip(RANK_PREFIX, RANK_LETTERS))
        weak_hits: dict[str, str] = {}
        for aid in filtered_ids:
            best = sintax_best.get(aid, {})
            for tok in filter_discard:
                # longest matching prefix, so "sg__" is not read as "s__"
                matches = [p for p in RANK_PREFIX if p and tok.startswith(p)]
                if not matches:
                    continue
                pfx = max(matches, key=len)
                letter = prefix_to_letter[pfx]
                if letter not in best:
                    continue
                name, confidence = best[letter]
                if (_drop_other(name) == _drop_other(tok[len(pfx):])
                        and confidence * 100 >= float(discard_min_support)):
                    weak_hits[aid] = f"{tok} at {confidence:.2f} confidence"
                    break
        for aid, why in weak_hits.items():
            logmsg(f"  {aid}: best hit {why} (reported as {tax_strings[aid]})")
            filter_status[aid] = "removed_best_hit"
        filtered_ids = [aid for aid in filtered_ids if aid not in weak_hits]
        logmsg(f"Discard filter, best-hit check (support >= {discard_min_support}%): "
               f"{len(weak_hits)} more ASV(s) removed")

    # ── Organelle filter (Metaxa2's own origin call; 16S) ──
    # Removes every ASV that Metaxa2 placed in its mitochondria or chloroplast FASTA
    # (amplicon.taxonomy.filter.organelles: metaxa2, the default). It runs last, so it
    # only removes what keep, discard and the best-hit check let through. That is
    # exactly the case it is for: plant mitochondria whose classifier best hit is a
    # bacterium, so no token can see them. On the 16S test set, 8 ASVs with a
    # mitochondrion as nearest SILVA reference (97-99 % identity) were classified as
    # bacteria with SINTAX on SILVA 138.2 (best hit Enterobacterales / Pantoea) and
    # passed every token, while Metaxa2 called all 8 mitochondrial.
    # Metaxa2's other calls (eukaryota, uncertain, the domains) remove nothing; the
    # keep tokens decide domains. The log reports how many organelle ASVs and reads
    # Metaxa2 found, how many of them the earlier steps had already removed, and lists
    # each ASV that only Metaxa2 caught.
    if organelle_filter:
        organelle_ids = [aid for aid in shared_ids
                         if metaxa2_origin.get(aid) in ("mitochondria", "chloroplast")]
        n_mito = sum(1 for aid in organelle_ids if metaxa2_origin[aid] == "mitochondria")
        n_chloro = len(organelle_ids) - n_mito
        # Organelle ASVs the earlier steps had already removed, per step.
        already = {step: [aid for aid in organelle_ids if filter_status[aid] == step]
                   for step in ("removed_keep", "removed_discard", "removed_best_hit")}
        already_ids = [aid for ids in already.values() for aid in ids]
        # Organelle ASVs still in the table: only Metaxa2 caught these.
        only_metaxa2 = [aid for aid in organelle_ids if filter_status[aid] == "kept"]
        for aid in only_metaxa2:
            filter_status[aid] = "removed_organelle"
        filtered_ids = [aid for aid in filtered_ids if filter_status[aid] == "kept"]

        logmsg(f"Organelle filter (Metaxa2 origin): {len(organelle_ids)} ASV(s), "
               f"{reads_of(organelle_ids)} reads, called mitochondria ({n_mito}) or "
               f"chloroplast ({n_chloro})")
        logmsg(f"  already removed above by the tokens or the best-hit check: "
               f"{len(already_ids)} ASV(s), {reads_of(already_ids)} reads "
               f"(keep {len(already['removed_keep'])}, discard {len(already['removed_discard'])}, "
               f"best-hit check {len(already['removed_best_hit'])})")
        logmsg(f"  removed only by the organelle filter: {len(only_metaxa2)} ASV(s), "
               f"{reads_of(only_metaxa2)} reads")
        for aid in only_metaxa2:
            reported = tax_strings[aid] or "unclassified"
            logmsg(f"    {aid}: {metaxa2_origin[aid]} by Metaxa2, {reads_of([aid])} reads "
                   f"(reported as {reported})")

logmsg(f"Final ASV count: {len(filtered_ids)}")

# Fail loud rather than write a header-only table: a filter that removes 100% of
# ASVs is almost always a marker/config mismatch (e.g. 16S keep tokens on an ITS run).
if filter_enabled and n_before > 0 and len(filtered_ids) == 0:
    msg = (f"the keep/discard filter removed ALL {n_before} ASVs (amp_type={amp_type}). "
           "The keep/discard lists likely do not match this marker: for ITS use keep: [k__Fungi]; "
           "for 16S use keep: [k__Bacteria, k__Archaea]. Fix amplicon.taxonomy.filter, "
           "or set amplicon.taxonomy.filter.enabled: false.")
    # If the organelle filter took part, say how much, because then the lists may
    # well be right and the sample simply holds little but organelle 16S.
    n_organelle_only = sum(1 for s in filter_status.values() if s == "removed_organelle")
    if n_organelle_only > 0:
        msg += (f" {n_organelle_only} of them were removed only because Metaxa2 called them "
                "mitochondria or chloroplast; amplicon.taxonomy.filter.organelles: off keeps those.")
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


def write_unfiltered_table(path: Path, row_ids: list[str]) -> None:
    """Write asv_table_unfiltered.txt: the layout of asv_table.txt for EVERY ASV
    that reached the filter (organelles and other removed ASVs included), plus two
    columns: metaxa2_origin (Metaxa2's call, or NA) and filter_status (kept, or
    the step that removed the ASV). NA is written unquoted, as R's write.table
    writes a missing value, so read.table and pandas both read it as missing."""
    with path.open("w") as fh:
        fh.write("\t".join([""] + [q(s) for s in samples]
                           + [q("taxonomy"), q("metaxa2_origin"), q("filter_status")]) + "\n")
        for asv_id in row_ids:
            vals = counts.get(asv_id, ["0"] * len(samples))
            origin = metaxa2_origin.get(asv_id)
            origin_cell = q(origin) if origin else "NA"
            fh.write("\t".join([q(asv_id)] + vals
                               + [q(tax_strings[asv_id]), origin_cell, q(filter_status[asv_id])]) + "\n")


# ── Output (d): asv_table_unfiltered.txt ──────────────────────────────────
# Every ASV that entered the filter, in the same order as asv_table.txt, so the
# kept rows of this table are asv_table.txt itself plus the two extra columns.
write_unfiltered_table(out_unfilt, shared_ids)

logmsg(f"DONE. Outputs written to {out_asv.parent}")
log.close()
