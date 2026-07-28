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
separated, no accession, no per-node tags:

  >k__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus_simiae

WHAT IT DOES (two edits per header, sequence lines untouched)
  1. drop the leading ">ACCESSION Root;"  -> the lineage now starts at k__Bacteria
  2. strip every " [id: N]" database tag

The result is a 7-rank NCBI RefSeq lineage with prefixes already embedded and a
binomial species — consumed by the rpoB marker pack as taxonomy_refdb (prefix_style
"embedded"). Input and output are both gzipped; sequence lines pass through verbatim
so ASV-facing sequences are never altered.

HEADER-SHAPE VALIDATION
Both regexes are anchored to the exact shape of the 2024-07-07 FROGS release (see
workflow/markers/rpoB.yaml's references.rpob_frogs.url, which pins that dated
filename — not a "latest" pointer). If a future re-download ever points at a newer
FROGS release with a different header shape, a silent no-op here would ship
corrupted taxonomy strings (still-present accession/Root/id-tag text mistaken for
a rank value) with no error at all. So each edit is checked: if a header doesn't
start with ">ACCESSION Root;", or still contains "[id:" after both substitutions,
the script stops and names the offending line rather than passing it through.

INPUT  : the raw FROGS FASTA, extracted from the .tar.gz by rule fetch_rpob_archive.
OUTPUT : the DADA2 trainset FASTA, used by rule assign_taxonomy for rpoB.
"""
import gzip
import re
import sys
from pathlib import Path

sm = snakemake  # noqa: F821 (injected by Snakemake)

in_path  = Path(sm.input.fasta)
out_path = Path(sm.output.fasta)
log_path = Path(sm.log[0])

# Edit 1: ">ACCESSION Root;"  ->  ">"   (accession = first run of non-space chars)
DROP_ACCESSION_ROOT = re.compile(r"^>[^ ]+ Root;")
# Edit 2: remove each " [id: N]" tag wherever it appears in the header
DROP_ID_TAG = re.compile(r" \[id: \d+\]")


def _open_maybe_gzip(path: Path, mode: str):
    """Open .gz transparently; plain text otherwise (so the script is testable
    on an un-gzipped FASTA outside the pipeline)."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


n_headers = 0
n_seq_lines = 0

out_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

with (_open_maybe_gzip(in_path, "rt") as fi,
      gzip.open(out_path, "wt") as fo,
      log_path.open("w") as log):

    for line_no, line in enumerate(fi, 1):
        if line.startswith(">"):
            # Apply both edits, in order, to the header (newline preserved), and
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
            line = new_line
            n_headers += 1
        else:
            n_seq_lines += 1
        fo.write(line)

    log.write(
        f"[frogs_rpob_to_dada2] {n_headers} headers rewritten, "
        f"{n_seq_lines} sequence line(s) passed through\n"
        f"[frogs_rpob_to_dada2] Output: {out_path}\n"
    )
