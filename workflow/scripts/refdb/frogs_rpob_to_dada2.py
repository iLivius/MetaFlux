#!/usr/bin/env python3
"""Reference-DB prep (workflow/scripts/refdb/): FROGS rpoB release -> DADA2 trainset.

One of the reference-preparation scripts. Each script here turns one upstream
database into the exact format a classifier wants; the rules in 10_refdb.smk just
fetch/extract and call these. See workflow/scripts/refdb/README.md for the full list.

WHY THIS EXISTS
The FROGS RefSeq rpoB release ships the sequences for FROGS' own affiliation tool,
so every FASTA header is an accession followed by a FROGS lineage in which every
node carries a "[id: N]" database tag, under a constant "Root":

  >WP_095092576.1 Root;k__Bacteria [id: 1];p__Bacillota [id: 2];...;s__Staphylococcus_simiae [id: 7]

DADA2's assignTaxonomy instead wants the header to BE the lineage — semicolon
separated, no accession, no per-node tags — and MetaFlux writes the first rank,
the domain, as d__ (edit 3 below):

  >d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus_simiae

WHAT IT DOES (three edits per header, sequence lines untouched)
  1. drop the leading ">ACCESSION Root;"  -> the lineage now starts at k__Bacteria
  2. strip every " [id: N]" database tag
  3. write the first rank as a domain: k__Bacteria -> d__Bacteria. FROGS puts the
     domain under the k__ prefix. This release carries NCBI's taxonomy as it was
     when FROGS built it (July 2024), before NCBI added kingdoms for prokaryotes
     (Pseudomonadati, Bacillati, ...), so every lineage has a domain and no
     kingdom. In MetaFlux's amplicon tables d__ is the domain (16S, 18S) and k__
     the kingdom rank (16S on SILVA 144, ITS), so the rpoB pack
     (workflow/markers/rpoB.yaml) names the first rank Domain, with the prefix
     d__, and this edit makes the file say the same. Only k__Bacteria and
     k__Archaea are rewritten; a header whose first rank is anything else stops
     the script (see HEADER-SHAPE VALIDATION below).

The result is a 7-rank NCBI RefSeq lineage with prefixes already embedded and a
binomial species — consumed by the rpoB marker pack as taxonomy_refdb (prefix_style
"embedded"). Input and output are both gzipped; sequence lines pass through verbatim
so ASV-facing sequences are never altered.

HEADER-SHAPE VALIDATION
The three regexes are anchored to the exact shape of the 2024-07-07 FROGS release
(see workflow/markers/rpoB.yaml's references.rpob_frogs.url, which pins that dated
filename — not a "latest" pointer). If a future re-download ever points at a newer
FROGS release with a different header shape, a silent no-op here would ship
corrupted taxonomy strings (still-present accession/Root/id-tag text mistaken for
a rank value, or a kingdom written as a domain) with no error at all. So each edit
is checked: if a header doesn't start with ">ACCESSION Root;", still contains
"[id:" once the tags are stripped, or then does not start with k__Bacteria or
k__Archaea (for example a first rank that is a kingdom, such as k__Pseudomonadati,
or one already written d__), the script stops and names the offending line rather
than passing it through.

INPUT  : the raw FROGS FASTA, extracted from the .tar.gz by rule fetch_rpob_archive.
OUTPUT : the DADA2 trainset FASTA, used by rule assign_taxonomy for rpoB.
"""
import gzip
import re
import sys
from pathlib import Path

sm = snakemake  # noqa: F821 (injected by Snakemake)

# Paths come from rule convert_rpob_to_dada2 in 10_refdb.smk:
#   input.fasta  = the raw FROGS FASTA (fetch_rpob_archive's output, still gzipped,
#                  still FROGS-shaped)
#   output.fasta = the DADA2 trainset FASTA this script writes (cached as rpob_frogs
#                  in rpoB.yaml; consumed next by rule assign_taxonomy for rpoB)
#   log[0]       = a plain-text log for a human to check after the run; nothing
#                  downstream reads it
in_path  = Path(sm.input.fasta)
out_path = Path(sm.output.fasta)
log_path = Path(sm.log[0])

# Edit 1: ">ACCESSION Root;"  ->  ">"   (accession = first run of non-space chars)
DROP_ACCESSION_ROOT = re.compile(r"^>[^ ]+ Root;")
# Edit 2: remove each " [id: N]" tag wherever it appears in the header
DROP_ID_TAG = re.compile(r" \[id: \d+\]")
# Edit 3: the first rank, FROGS' k__ followed by a domain name, is written as a
# domain, d__ (see "WHAT IT DOES" above). The name is captured so it is kept as it
# is, and the ";" after it makes the pattern match the whole first rank only.
KINGDOM_LABELLED_DOMAIN = re.compile(r"^>k__(Bacteria|Archaea);")


def _open_maybe_gzip(path: Path, mode: str):
    """Open .gz transparently; plain text otherwise (so the script is testable
    on an un-gzipped FASTA outside the pipeline)."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


# Running totals, used only for the summary line written to the log at the
# end — DADA2 doesn't read the log, so these exist purely so a human can
# sanity-check the conversion afterwards (right number of sequences, etc.)
n_headers = 0
n_seq_lines = 0

out_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

with (_open_maybe_gzip(in_path, "rt") as fi,
      gzip.open(out_path, "wt") as fo,
      log_path.open("w") as log):

    # Walk the FASTA one line at a time. A header line (">...") is rewritten
    # below; every other line is sequence data and is copied through as-is.
    for line_no, line in enumerate(fi, 1):
        if line.startswith(">"):
            # Apply the three edits, in order, to the header (newline preserved), and
            # confirm each one actually fired — a header that doesn't match means
            # the upstream FROGS format has changed since this script was written.
            new_line, n_root_sub = DROP_ACCESSION_ROOT.subn(">", line)
            if n_root_sub == 0:
                log.write(
                    f"[frogs_rpob_to_dada2] ERROR: header on line {line_no} does not "
                    f"match the expected '>ACCESSION Root;...' shape: {line.strip()!r}\n"
                )
                sys.exit(
                    f"[frogs_rpob_to_dada2] {in_path}, line {line_no}: header does not "
                    f"match the expected FROGS shape '>ACCESSION Root;...' ({line.strip()!r}). "
                    "The upstream FROGS rpoB release may have changed its header format — "
                    "compare a few headers from the new download against this script's "
                    "docstring, then update DROP_ACCESSION_ROOT/DROP_ID_TAG to match."
                )
            new_line, _ = DROP_ID_TAG.subn("", new_line)
            if "[id:" in new_line:
                log.write(
                    f"[frogs_rpob_to_dada2] ERROR: header on line {line_no} still has an "
                    f"unremoved '[id: N]' tag after conversion: {new_line.strip()!r}\n"
                )
                sys.exit(
                    f"[frogs_rpob_to_dada2] {in_path}, line {line_no}: an '[id: N]' tag "
                    f"survived the conversion ({new_line.strip()!r}). The tag format may "
                    "have changed — update DROP_ID_TAG to match the new shape."
                )
            # Edit 3: k__Bacteria -> d__Bacteria (and k__Archaea -> d__Archaea). After
            # edits 1 and 2 the header starts with its first rank. If that rank is not
            # k__ followed by a domain name, the FROGS release has changed its rank
            # model (a first rank that is a kingdom, or one already written d__):
            # stop rather than guess what the first rank holds.
            new_line, n_domain_sub = KINGDOM_LABELLED_DOMAIN.subn(r">d__\1;", new_line)
            if n_domain_sub == 0:
                log.write(
                    f"[frogs_rpob_to_dada2] ERROR: header on line {line_no} does not "
                    f"start with k__Bacteria or k__Archaea after conversion: "
                    f"{new_line.strip()!r}\n"
                )
                sys.exit(
                    f"[frogs_rpob_to_dada2] {in_path}, line {line_no}: the first rank is "
                    f"not k__Bacteria or k__Archaea ({new_line.strip()!r}). FROGS writes "
                    "the domain under k__ and this script rewrites it to d__; any other "
                    "first rank means the FROGS release has changed its rank model. "
                    "Compare a few headers from the new download with this script's "
                    "docstring and with the rank model in workflow/markers/rpoB.yaml, "
                    "then update KINGDOM_LABELLED_DOMAIN and the pack together."
                )
            line = new_line
            n_headers += 1
        else:
            n_seq_lines += 1
        fo.write(line)

    # Final tally for whoever checks the log after a refdb rebuild — not parsed
    # by Snakemake or any downstream rule.
    log.write(
        f"[frogs_rpob_to_dada2] {n_headers} headers rewritten, "
        f"{n_seq_lines} sequence line(s) passed through\n"
        f"[frogs_rpob_to_dada2] Output: {out_path}\n"
    )
