#!/usr/bin/env python3
"""Convert the SILVA DADA2 RDP trainset to VSEARCH SINTAX format.

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

in_path  = Path(sm.input.fasta)
out_path = Path(sm.output.fasta)
log_path = Path(sm.log[0])

RANK_LETTERS = ("d", "p", "c", "o", "f", "g")

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
            fo.write(line + "\n")
            n_sequences += 1
            continue

        # Header is the full lineage: >Kingdom;Phylum;Class;Order;Family;Genus;
        # Strip leading ">" and split on ";", dropping empty tokens.
        ranks = [r.strip() for r in line[1:].split(";") if r.strip()]

        pairs = [
            f"{RANK_LETTERS[i]}:{ranks[i]}"
            for i in range(min(len(ranks), len(RANK_LETTERS)))
        ]

        seq_id += 1
        fo.write(f">{seq_id};tax={','.join(pairs)}\n")
        n_converted += 1

    log.write(
        f"[silva_to_sintax] {n_converted} headers converted, "
        f"{n_sequences} sequence lines written\n"
        f"[silva_to_sintax] Output: {out_path}\n"
    )
