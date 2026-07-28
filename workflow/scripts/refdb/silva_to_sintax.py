#!/usr/bin/env python3
"""Reference-DB prep (workflow/scripts/refdb/): SILVA DADA2 trainset -> VSEARCH SINTAX.

One of the reference-preparation scripts. Each script here turns one upstream
database into the exact format a classifier wants; the rules in 10_refdb.smk just
fetch and call these. See workflow/scripts/refdb/README.md for the full list.

SCOPE: this converter is SILVA/16S-specific — it assumes a BARE 6-rank lineage
(kingdom..genus, no species) with no rank prefixes and no accession field. It is
NOT a general DADA2->SINTAX tool: a marker with embedded prefixes (k__/p__),
a species rank, or a gene tag (gyrB) would need a rank-model-aware converter. When
a second marker needs SINTAX, generalise this using the pack's rank metadata
(tax_levels / rank_letters / rank_prefixes / prefix_style) rather than copying it.

The SILVA DADA2 trainset (toGenus) uses the full lineage as the FASTA header
with no separate accession field:

  >Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Vibrionaceae;Vibrio;

VSEARCH SINTAX requires each record to have a unique ID and a ;tax= tag:

  >1;tax=d:Bacteria,p:Pseudomonadota,c:Gammaproteobacteria,o:Enterobacterales,f:Vibrionaceae,g:Vibrio

A sequential integer is used as the unique ID — database IDs are only used
internally by VSEARCH and do not appear in the tabbedout classification output.
Trailing empty ranks (from the trailing semicolon) and genuinely missing
higher-level ranks are silently skipped; VSEARCH handles partial lineages fine.
"""
import gzip
from pathlib import Path

sm = snakemake  # noqa: F821

# Paths come from rule convert_silva_sintax in 10_refdb.smk:
#   input.fasta  = the SILVA DADA2 trainset (fetch_silva_train's downloaded file)
#   output.fasta = the SINTAX-formatted FASTA this script writes (cached as
#                  silva_sintax in 16S.yaml; consumed next by VSEARCH's sintax
#                  command when taxonomy.method: sintax is chosen for 16S)
#   log[0]       = a plain-text log for a human to check after the run; nothing
#                  downstream reads it
in_path  = Path(sm.input.fasta)
out_path = Path(sm.output.fasta)
log_path = Path(sm.log[0])

# SINTAX rank-letter codes, one per lineage position in order: domain, phylum,
# class, order, family, genus. Matches 16S.yaml's rank_letters entry.
RANK_LETTERS = ("d", "p", "c", "o", "f", "g")

# Running totals for the log summary at the end, plus seq_id: the sequential
# integer used as each record's new SINTAX ID (see the WHY note above).
n_converted = 0
n_sequences = 0
seq_id      = 0

out_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

with (gzip.open(in_path,  "rt") as fi,
      gzip.open(out_path, "wt") as fo,
      log_path.open("w")        as log):

    for line in fi:
        line = line.rstrip("\n")

        if not line.startswith(">"):
            # Sequence line (the DNA letters) — passed through unchanged, just
            # re-adding the newline that was stripped above.
            fo.write(line + "\n")
            n_sequences += 1
            continue

        # Header is the full lineage: >Kingdom;Phylum;Class;Order;Family;Genus;
        # Strip leading ">" and split on ";", dropping empty tokens.
        ranks = [r.strip() for r in line[1:].split(";") if r.strip()]

        # Pair each rank letter with the lineage value at the same position
        # (RANK_LETTERS[0]=d with ranks[0]=Kingdom, and so on). If this lineage
        # has fewer than 6 ranks — a shallow or partial SILVA classification —
        # stop at whatever is actually present instead of inventing empty ranks
        # for the rest (see the module docstring's note on missing ranks).
        pairs = [
            f"{RANK_LETTERS[i]}:{ranks[i]}"
            for i in range(min(len(ranks), len(RANK_LETTERS)))
        ]

        seq_id += 1
        fo.write(f">{seq_id};tax={','.join(pairs)}\n")
        n_converted += 1

    # Summary for a human checking the refdb rebuild; nothing downstream parses
    # this log.
    log.write(
        f"[silva_to_sintax] {n_converted} headers converted, "
        f"{n_sequences} sequence lines written\n"
        f"[silva_to_sintax] Output: {out_path}\n"
    )
