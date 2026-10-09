#!/usr/bin/env python3
"""Reference-DB prep (workflow/scripts/refdb/): SILVA 144 DADA2 trainset -> 6-rank trainset
plus a phylum -> kingdom table.

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

MetaFlux's 16S rank model has six slots before species (domain written as k__,
then p__ ... g__). Fed the 7-rank file, every name lands one slot too deep, the
o__Chloroplast / f__Mitochondria filter tokens stop matching, and nothing errors.
So this script removes the Kingdom field before any classifier is trained, which
keeps every downstream table and filter exactly as it was with 138.2, and records
which kingdom each phylum belongs to so the kingdom can be put back afterwards
(rule restore_kingdom, 80_taxonomy.smk, writes it to a separate sidecar table).

The kingdom can be recovered from the phylum because, in SILVA 144, a phylum sits
under one kingdom only. The one exception is the placeholder phylum "Incertae
Sedis" (sequences SILVA could not place in a named phylum), which occurs under more
than one kingdom; the table marks it "ambiguous" and leaves its kingdom empty
rather than picking one.

INPUT (rule convert_silva_train):
  input.fasta        = silva_nr99_v144_toGenus_trainset.fa.gz exactly as published
                       by SILVA (fetched and md5-checked by fetch_silva_train_raw)
  params.strip_other = True when amplicon.taxonomy.silva_other is "strip"
OUTPUT:
  output.fasta       = the 6-rank training set (gzip). Used by assignTaxonomy (rdp
                       path), the amplicon probe (in-silico PCR, which only reads
                       the sequences) and convert_silva_sintax (sintax path).
  output.kingdom_map = TSV, one row per (domain, phylum): the kingdom, how many
                       reference sequences support it, and whether it is unique.
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
from collections import Counter, defaultdict
from pathlib import Path

sm = snakemake  # noqa: F821 (injected by Snakemake)

in_path     = Path(sm.input.fasta)
out_fasta   = Path(sm.output.fasta)
out_map     = Path(sm.output.kingdom_map)
log_path    = Path(sm.log[0])
strip_other = bool(sm.params.strip_other)

OTHER_SUFFIX     = "--other"
DOMAINS          = ("Bacteria", "Archaea", "Eukaryota")
PROKARYOTES      = ("Bacteria", "Archaea")
EXPECTED_MAX_DEPTH = 7     # Domain;Kingdom;Phylum;Class;Order;Family;Genus

out_fasta.parent.mkdir(parents=True, exist_ok=True)
out_map.parent.mkdir(parents=True, exist_ok=True)
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
    log.write(f"[silva_drop_kingdom] ERROR: {message}\n")
    log.close()
    sys.exit(f"[silva_drop_kingdom] {message}")


# Running totals for the log, and the raw material of the kingdom table:
# kingdoms_of[(domain, phylum)] counts how many reference sequences put that
# phylum under each kingdom.
n_records          = 0
n_sequence_lines   = 0
n_names_stripped   = 0
n_kingdom_only     = 0     # records that stop at kingdom (no phylum)
depth_before       = Counter()
domain_counts      = Counter()
kingdom_counts     = Counter()
kingdoms_of        = defaultdict(Counter)

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
        depth_before[len(ranks)] += 1
        domain_counts[domain] += 1

        # ── Split off the kingdom (prokaryotes with at least 2 ranks) ──
        if domain in PROKARYOTES and len(ranks) >= 2:
            kingdom = ranks[1]
            if not looks_like_kingdom(kingdom):
                fail(f"{in_path}:{line_number}: second rank {kingdom!r} is not a kingdom "
                     f"name — this does not look like a SILVA 144 file (a 6-rank SILVA "
                     f"138.2 file has the phylum there): {line!r}")
            kept = [domain] + ranks[2:]
        else:
            # Domain-only records (incl. the Eukaryota outgroup) have no kingdom field.
            kingdom = None
            kept = ranks

        # ── Optional "--other" removal (amplicon.taxonomy.silva_other: strip) ──
        if strip_other:
            stripped = [strip_other_suffix(r) for r in kept]
            n_names_stripped += sum(1 for a, b in zip(kept, stripped) if a != b)
            kept = stripped
            if kingdom is not None:
                kingdom = strip_other_suffix(kingdom)

        # ── Record the phylum -> kingdom pair for the lookup table ──
        if kingdom is not None:
            kingdom_counts[kingdom] += 1
            if len(kept) >= 2:
                phylum = kept[1]
                kingdoms_of[(domain, phylum)][kingdom] += 1
            else:
                # e.g. "Bacteria;Pseudomonadati;" — becomes domain-only after the
                # drop, so there is no phylum to attach the kingdom to.
                n_kingdom_only += 1

        fo.write(">" + ";".join(kept) + ";\n")

    # ── File-level shape checks ──
    if n_records == 0:
        fail(f"{in_path} contains no FASTA headers.")
    max_depth = max(depth_before)
    if max_depth != EXPECTED_MAX_DEPTH:
        fail(f"the deepest lineage in {in_path} has {max_depth} ranks; a SILVA 144 toGenus "
             f"training set has {EXPECTED_MAX_DEPTH} (Domain..Genus).")

    # ── Write the phylum -> kingdom table ──
    # One row per (domain, phylum). "unique": every reference with this phylum sits
    # under the same kingdom, which is the kingdom restored later. "ambiguous": the
    # phylum occurs under more than one kingdom (in SILVA 144 only the placeholder
    # phylum "Incertae Sedis"), so no kingdom is restored for it and the candidates
    # are listed in the note column instead.
    n_unique = n_ambiguous = 0
    with out_map.open("w") as fm:
        fm.write("domain\tphylum\tkingdom\tn_refs\tstatus\tnote\n")
        for (domain, phylum) in sorted(kingdoms_of):
            counts = kingdoms_of[(domain, phylum)]
            n_refs = sum(counts.values())
            if len(counts) == 1:
                kingdom = next(iter(counts))
                status  = "unique"
                note    = ""
                if kingdom == "Incertae Sedis":
                    note = "SILVA 144 places this phylum in no named kingdom"
                n_unique += 1
            else:
                kingdom = ""
                status  = "ambiguous"
                candidates = ", ".join(f"{k} ({n})" for k, n in counts.most_common())
                note = f"placed under more than one kingdom in SILVA 144: {candidates}"
                n_ambiguous += 1
            fm.write(f"{domain}\t{phylum}\t{kingdom}\t{n_refs}\t{status}\t{note}\n")

    # ── Summary for a human checking the refdb rebuild ──
    log.write(
        f"[silva_drop_kingdom] input:  {in_path}\n"
        f"[silva_drop_kingdom] {n_records} headers, {n_sequence_lines} sequence lines "
        f"(sequences unchanged)\n"
        f"[silva_drop_kingdom] records per domain: {dict(domain_counts)}\n"
        f"[silva_drop_kingdom] ranks per record before the drop: {dict(sorted(depth_before.items()))}\n"
        f"[silva_drop_kingdom] records per kingdom: {dict(kingdom_counts.most_common())}\n"
        f"[silva_drop_kingdom] {n_kingdom_only} records stop at kingdom and became domain-only\n"
        f"[silva_drop_kingdom] silva_other={'strip' if strip_other else 'verbatim'}; "
        f"rank names with '--other' removed: {n_names_stripped}\n"
        f"[silva_drop_kingdom] kingdom table: {len(kingdoms_of)} (domain, phylum) pairs, "
        f"{n_unique} unique, {n_ambiguous} ambiguous\n"
        f"[silva_drop_kingdom] output: {out_fasta}\n"
        f"[silva_drop_kingdom] output: {out_map}\n"
    )
