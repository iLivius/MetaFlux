#!/usr/bin/env python3
"""Put SILVA 144's Kingdom rank back next to each ASV's final taxonomy (16S only).

WHY THIS STEP EXISTS
SILVA 144 added a Kingdom rank between Domain and Phylum (Bacillati, Pseudomonadati,
Methanobacteriati, ...). MetaFlux classifies on six ranks instead — the Kingdom field
is removed from the training set beforehand (rule convert_silva_train,
workflow/scripts/refdb/silva_drop_kingdom.py) — so that every taxonomy string,
filter token (k__Bacteria, o__Chloroplast, f__Mitochondria) and downstream table keeps
the shape it had with SILVA 138.2. This script restores the kingdom afterwards.

That works because in SILVA 144 the kingdom follows from the phylum: every phylum
sits under one kingdom only, so knowing an ASV's (domain, phylum) is enough. The
lookup table comes from the same reference file the ASVs were classified against,
written by convert_silva_train next to the 6-rank training set.

The existing tables (asv_table.txt, asv_table_seqs.txt, taxon_seq_table.txt) are left
exactly as they are. The kingdom goes into a separate SIDECAR file, so nothing that
already reads those tables can be affected by it.

INPUT (rule restore_kingdom, 80_taxonomy.smk):
  taxon_seq_table = 6.taxonomy/taxon_seq_table.txt from assign_taxonomy (either the
                    rdp or the sintax implementation — both write the same layout):
                    one row per ASV left after the keep/discard filter, one column
                    per rank (Kingdom..Species, where MetaFlux's "Kingdom" column
                    holds the DOMAIN) plus the sequence. Written with R write.table
                    conventions: values in double quotes; a missing rank is "" (sintax)
                    or an unquoted NA (rdp).
  kingdom_map     = the (domain, phylum) -> kingdom TSV from convert_silva_train.
  params.tax_levels = the rank column names used in taxon_seq_table.
OUTPUT:
  6.taxonomy/taxonomy_lineage_full.tsv — tab-separated, one row per ASV, same ASVs
  and order as taxon_seq_table.txt, with these columns:
    asv_id        ASV_N, as in every other table
    Domain        the value MetaFlux's tables call "Kingdom" (Bacteria / Archaea)
    Kingdom       SILVA 144's kingdom, restored from the phylum (may be empty, see below)
    Phylum .. Genus, Species   copied unchanged from taxon_seq_table.txt
    lineage       the SILVA-style path Domain;Kingdom;Phylum;...;Genus; (SILVA's own
                  tax_slv layout: bare names, ';' after each), down to the deepest
                  rank the ASV was assigned. A rank left empty in the middle keeps an
                  empty slot (";;"), so a name's position always tells its rank.
    kingdom_note  why the kingdom is empty, or a remark about it; empty when the
                  kingdom was restored without any caveat
  Nothing in the pipeline reads this file back; it is a deliverable for the user's
  own analysis, joined to the abundance table by asv_id.

WHEN THE KINGDOM IS LEFT EMPTY (each case says so in kingdom_note)
  - the ASV has no domain, or only a domain: the kingdom is read off the phylum, so
    without a phylum there is nothing to look it up from (guessing from the domain
    alone would be wrong — Bacteria spans several kingdoms);
  - the phylum is SILVA's placeholder "Incertae Sedis", which occurs under more than
    one kingdom (marked "ambiguous" in the table); picking one would invent a
    placement SILVA itself does not make;
  - the domain is Eukaryota: SILVA 144's DADA2 file carries eukaryotes at domain level
    only, so it holds no eukaryote kingdom to restore.
A (domain, phylum) pair that is missing from the table altogether is NOT left empty:
it means the table and the reference the ASVs were classified against do not match
(e.g. a reference path pinned by hand), and the script stops with an error.
"""
import sys
from pathlib import Path

sm = snakemake  # noqa: F821 (injected by Snakemake)

taxon_table_path = Path(sm.input.taxon_seq_table)
kingdom_map_path = Path(sm.input.kingdom_map)
out_path         = Path(sm.output.lineage_full)
log_path         = Path(sm.log[0])
tax_levels       = list(sm.params.tax_levels)   # e.g. Kingdom, Phylum, ..., Species

out_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)
log = log_path.open("w")


def logmsg(msg: str) -> None:
    log.write(f"[restore_kingdom] {msg}\n")
    log.flush()


def cell_value(raw: str) -> str:
    """One field of an R write.table file -> its plain value ("" when missing).

    R writes a missing value as a bare NA (no quotes) and every real character value
    inside double quotes, so a bare NA is the only "missing" marker; the sintax
    implementation writes missing ranks as "" instead.
    """
    if raw == "NA":
        return ""
    if len(raw) >= 2 and raw.startswith('"') and raw.endswith('"'):
        return raw[1:-1]
    return raw


# ── Read the (domain, phylum) -> kingdom table ──
kingdom_of: dict[tuple[str, str], dict[str, str]] = {}
with kingdom_map_path.open() as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        row = dict(zip(header, line.rstrip("\n").split("\t")))
        kingdom_of[(row["domain"], row["phylum"])] = row
logmsg(f"kingdom table: {len(kingdom_of)} (domain, phylum) pairs from {kingdom_map_path}")

# ── Read the final per-rank taxonomy table ──
with taxon_table_path.open() as fh:
    header = [cell_value(c) for c in fh.readline().rstrip("\n").split("\t")]
    rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]

missing_cols = [c for c in tax_levels if c not in header]
if missing_cols:
    sys.exit(f"[restore_kingdom] {taxon_table_path} has no column(s) {missing_cols}; "
             f"header is {header}")
col_index = {name: header.index(name) for name in tax_levels}

domain_col  = tax_levels[0]      # MetaFlux calls the domain slot "Kingdom"
phylum_col  = tax_levels[1]
species_col = tax_levels[-1]
# SILVA ranks between domain and species, in order (Phylum..Genus).
lineage_cols = tax_levels[1:-1]

out_header = (["asv_id", "Domain", "Kingdom"] + lineage_cols + [species_col]
              + ["lineage", "kingdom_note"])

n_restored = 0
n_empty    = 0
with out_path.open("w") as fo:
    fo.write("\t".join(out_header) + "\n")

    for fields in rows:
        asv_id = cell_value(fields[0])
        values = {name: cell_value(fields[idx]) for name, idx in col_index.items()}
        domain = values[domain_col]
        phylum = values[phylum_col]

        # ── Decide the kingdom (see WHEN THE KINGDOM IS LEFT EMPTY above) ──
        kingdom = ""
        note    = ""
        if not domain:
            note = "not classified to domain"
        elif domain == "Eukaryota":
            note = "SILVA 144's DADA2 file gives eukaryotes the domain only, no kingdom"
        elif not phylum:
            note = ("classified to domain only: the kingdom is read off the phylum, "
                    "so none can be given")
        else:
            entry = kingdom_of.get((domain, phylum))
            if entry is None:
                msg = (f"{asv_id}: ({domain}, {phylum}) is not in the kingdom table "
                       f"{kingdom_map_path}. The table and the reference this run was "
                       f"classified against do not match.")
                logmsg("ERROR: " + msg)
                log.close()
                sys.exit(f"[restore_kingdom] {msg}")
            kingdom = entry["kingdom"]   # empty for an ambiguous phylum
            note    = entry["note"]

        if kingdom:
            n_restored += 1
        else:
            n_empty += 1

        # ── SILVA-style lineage: Domain;Kingdom;Phylum;...;Genus; ──
        # Stop at the deepest assigned rank; keep empty slots before it so every
        # name stays at its rank's position.
        path = [domain, kingdom] + [values[c] for c in lineage_cols]
        while path and path[-1] == "":
            path.pop()
        lineage = "".join(name + ";" for name in path)

        out_row = ([asv_id, domain, kingdom]
                   + [values[c] for c in lineage_cols]
                   + [values[species_col], lineage, note])
        fo.write("\t".join(out_row) + "\n")

logmsg(f"{len(rows)} ASVs: kingdom restored for {n_restored}, left empty for {n_empty} "
       f"(reason in the kingdom_note column)")
logmsg(f"DONE. Wrote {out_path}")
log.close()
