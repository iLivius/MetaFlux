#!/usr/bin/env python3
"""Reference-DB prep (workflow/scripts/refdb/): SILVA 144 DADA2 trainset -> the 7-rank
training set MetaFlux classifies with, kingdom names without SILVA's "--other".

One of the reference-preparation scripts. Each script here turns one upstream
database into the exact format a classifier wants; the rules in 10_refdb.smk just
fetch and call these. See workflow/scripts/refdb/README.md for the full list.

WHY THIS STEP EXISTS
SILVA 144 introduced a seven-rank prokaryote taxonomy, adding Kingdom between Domain
and Phylum (following the 2024 change to the prokaryote code that made kingdom a
formal rank). The DADA2 training set published by SILVA therefore reads

  >Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;

where SILVA 138.2 read

  >Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;

MetaFlux classifies on exactly these seven ranks (16S.yaml, release "144": Domain,
Kingdom, Phylum ... Genus, plus the Species slot that addSpecies fills), so the
domain lands in d__ and the kingdom in k__, each with its own confidence. The file
is still not used exactly as published, because of SILVA's "--other" suffix: SILVA
writes "X--other" for a clade that does not contain the sequence of the
nomenclatural type of X, but whose sequences SILVA still counts as part of X. So
this script, before any classifier is trained:
  1. writes every kingdom WITHOUT the suffix, always, whatever
     amplicon.taxonomy.silva_other says. At kingdom level the suffix is a curation
     note about SILVA's tree, not a separate taxon: it splits one kingdom into two
     labels (Bacteroidota, Acidobacteriota, Planctomycetota and Verrucomicrobiota
     all sit in "Pseudomonadati--other", Pseudomonadota in "Pseudomonadati";
     Methanobacteriati is split the same way), and with two labels fewer ASVs reach
     the confidence cutoff at kingdom level, while no other rank changes. The
     keep/discard filter treats "X--other" as X anyway;
  2. removes the suffix from the other ranks (phylum to genus) only when
     amplicon.taxonomy.silva_other is "strip"; with "verbatim" they keep it.

INPUT (rule convert_silva_train):
  input.fasta        = silva_nr99_v144_toGenus_trainset.fa.gz exactly as published
                       by SILVA (fetched and md5-checked by fetch_silva_train_raw)
  params.strip_other = True when amplicon.taxonomy.silva_other is "strip"
OUTPUT:
  output.fasta       = the 7-rank training set (gzip). Used by assignTaxonomy (rdp
                       path), the amplicon probe (in-silico PCR, which only reads
                       the sequences) and convert_silva_sintax (sintax path).
  log[0]             = plain-text summary for a human; nothing reads it back.

Sequence lines are never altered — only the FASTA headers change, so the ASVs are
compared against exactly the sequences SILVA published.

STRICT SHAPE CHECK
The script refuses (exits with a message) rather than guesses when a header does
not have the shape of a SILVA 144 toGenus header — so it cannot silently mangle a
SILVA 138.2 file, a toSpecies file (8 ranks) or a corrupted download:
  - every header ends with ";" and has no empty rank in the middle;
  - the first rank is Bacteria, Archaea or Eukaryota;
  - Eukaryota records carry the domain only (the DADA2 build keeps 500 of them, at
    domain level, as an outgroup);
  - Bacteria/Archaea records have between 1 and 7 ranks, and the deepest record in
    the file has exactly 7;
  - the second rank of a prokaryote record looks like a kingdom: kingdom names end
    in "-ati" under the prokaryote code (Bacillati, Pseudomonadati, ...), possibly
    followed by "--other", or it is "Incertae Sedis". A 138.2 file, whose second
    rank is a phylum (Pseudomonadota, Bacillota, ...), fails here.
"""
import gzip
import sys
from collections import Counter
from pathlib import Path

sm = snakemake  # noqa: F821 (injected by Snakemake)

in_path     = Path(sm.input.fasta)
out_fasta   = Path(sm.output.fasta)
log_path    = Path(sm.log[0])
strip_other = bool(sm.params.strip_other)

OTHER_SUFFIX     = "--other"
DOMAINS          = ("Bacteria", "Archaea", "Eukaryota")
PROKARYOTES      = ("Bacteria", "Archaea")
EXPECTED_MAX_DEPTH = 7     # Domain;Kingdom;Phylum;Class;Order;Family;Genus

out_fasta.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)


def strip_other_suffix(name: str) -> str:
    """Remove SILVA's "--other" suffix from a rank name.

    SILVA 144 also merges two intermixed taxa into one label joined by "+", each
    part with its own suffix (e.g. "Enterobacter--other+Klebsiella--other"), so the
    suffix is removed from every "+"-separated part: "Enterobacter+Klebsiella".
    """
    parts = []
    for part in name.split("+"):
        if part.endswith(OTHER_SUFFIX):
            part = part[: -len(OTHER_SUFFIX)]
        parts.append(part)
    return "+".join(parts)


def looks_like_kingdom(name: str) -> bool:
    """True if a rank name has the form of a SILVA 144 kingdom (see module docstring)."""
    base = name[: -len(OTHER_SUFFIX)] if name.endswith(OTHER_SUFFIX) else name
    return base == "Incertae Sedis" or base.endswith("ati")


def fail(message: str) -> None:
    """Stop the rule with a message in both the log and Snakemake's console output."""
    log.write(f"[silva_prepare_trainset] ERROR: {message}\n")
    log.close()
    sys.exit(f"[silva_prepare_trainset] {message}")


# Running totals for the summary in the log.
n_records          = 0
n_sequence_lines   = 0
n_kingdom_renamed  = 0     # kingdom names written without "--other" (always done)
n_other_renamed    = 0     # phylum-to-genus names that lost it (silva_other: strip only)
depth_counts       = Counter()
domain_counts      = Counter()
kingdom_counts     = Counter()   # kingdom names as written to the output

# compresslevel 6 instead of Python's default 9: the file is a local cache that is
# rebuilt on demand, and level 9 makes this step several times slower for a few
# percent smaller output.
with (gzip.open(in_path,   "rt") as fi,
      gzip.open(out_fasta, "wt", compresslevel=6) as fo,
      log_path.open("w")         as log):

    for line_number, line in enumerate(fi, start=1):
        line = line.rstrip("\n")

        if not line.startswith(">"):
            # Sequence line — passed through unchanged.
            fo.write(line + "\n")
            n_sequence_lines += 1
            continue

        n_records += 1
        body = line[1:]

        # ── Shape check for this header ──
        if not body.endswith(";"):
            fail(f"{in_path}:{line_number}: header does not end with ';' — this is not a "
                 f"SILVA DADA2 training-set header: {line!r}")
        ranks = body[:-1].split(";")
        if any(r == "" for r in ranks):
            fail(f"{in_path}:{line_number}: header has an empty rank: {line!r}")
        domain = ranks[0]
        if domain not in DOMAINS:
            fail(f"{in_path}:{line_number}: first rank {domain!r} is not one of {DOMAINS}: {line!r}")
        if domain == "Eukaryota" and len(ranks) != 1:
            fail(f"{in_path}:{line_number}: a Eukaryota record with more than the domain "
                 f"rank — not the shape of the SILVA 144 DADA2 file: {line!r}")
        if domain in PROKARYOTES and len(ranks) > EXPECTED_MAX_DEPTH:
            fail(f"{in_path}:{line_number}: {len(ranks)} ranks, more than the "
                 f"{EXPECTED_MAX_DEPTH} of a toGenus file (a toSpecies file?): {line!r}")
        depth_counts[len(ranks)] += 1
        domain_counts[domain] += 1

        # ── Kingdom (second rank of a prokaryote record): always without "--other" ──
        # Domain-only records (incl. the Eukaryota outgroup) have no kingdom field
        # and are written as they are.
        if domain in PROKARYOTES and len(ranks) >= 2:
            kingdom = ranks[1]
            if not looks_like_kingdom(kingdom):
                fail(f"{in_path}:{line_number}: second rank {kingdom!r} is not a kingdom "
                     f"name — this does not look like a SILVA 144 file (a SILVA 138.2 "
                     f"file, which has no kingdom, has the phylum there): {line!r}")
            plain_kingdom = strip_other_suffix(kingdom)
            if plain_kingdom != kingdom:
                n_kingdom_renamed += 1
            ranks[1] = plain_kingdom
            kingdom_counts[plain_kingdom] += 1

        # ── Phylum to genus: "--other" removed only with silva_other: strip ──
        # The domain (position 0) never carries the suffix (checked above), and the
        # kingdom (position 1) has just been handled, so the loop starts at the phylum.
        if strip_other:
            for position in range(2, len(ranks)):
                stripped = strip_other_suffix(ranks[position])
                if stripped != ranks[position]:
                    n_other_renamed += 1
                    ranks[position] = stripped

        fo.write(">" + ";".join(ranks) + ";\n")

    # ── File-level shape checks ──
    if n_records == 0:
        fail(f"{in_path} contains no FASTA headers.")
    max_depth = max(depth_counts)
    if max_depth != EXPECTED_MAX_DEPTH:
        fail(f"the deepest lineage in {in_path} has {max_depth} ranks; a SILVA 144 toGenus "
             f"training set has {EXPECTED_MAX_DEPTH} (Domain..Genus).")

    # ── Summary for a human checking the refdb rebuild ──
    log.write(
        f"[silva_prepare_trainset] input:  {in_path}\n"
        f"[silva_prepare_trainset] {n_records} headers, {n_sequence_lines} sequence lines "
        f"(sequences unchanged)\n"
        f"[silva_prepare_trainset] records per domain: {dict(domain_counts)}\n"
        f"[silva_prepare_trainset] ranks per record: {dict(sorted(depth_counts.items()))}\n"
        f"[silva_prepare_trainset] records per kingdom (as written): "
        f"{dict(kingdom_counts.most_common())}\n"
        f"[silva_prepare_trainset] kingdom names written without '--other': {n_kingdom_renamed}\n"
        f"[silva_prepare_trainset] silva_other={'strip' if strip_other else 'verbatim'}; "
        f"phylum-to-genus names with '--other' removed: {n_other_renamed}\n"
        f"[silva_prepare_trainset] output: {out_fasta}\n"
    )
