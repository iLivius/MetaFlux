#!/usr/bin/env python3
"""Reference-DB prep (workflow/scripts/refdb/): PR2 UTAX file -> SINTAX file, domain first.

One of the reference-preparation scripts. Each script here turns one upstream
database into the exact format a classifier wants; the rules in 10_refdb.smk just
fetch and call these. See workflow/scripts/refdb/README.md for the full list.

WHY THIS EXISTS
PR2 (v5.1.1) writes the headers of its UTAX/SINTAX file with the DOMAIN under the
letter k: and the SUPERGROUP under the letter d:

  >AB353770.1.1740_U;tax=k:Eukaryota,d:TSAR,p:Alveolata-Dinoflagellata,c:Dinophyceae,...

VSEARCH's --sintax does not read these letters as rank names. It reads them as a
fixed top-down order, d above k above p, c, o, f, g, s, and it counts the
confidence of a rank only over the bootstrap rounds that agreed on every rank above
it. With PR2's file as published, the supergroup is therefore VSEARCH's top rank
and the domain hangs below it: the domain can never be more certain than the
supergroup. Whenever the supergroup stays below sintax_cutoff, the domain is
dropped too, even when nearly every bootstrap round agreed on it, and the ASV ends
up with no taxonomy at all.

WHAT IT DOES (headers only; sequence lines are copied unchanged)
It swaps the first two values, so that the letters follow the ranks top-down:

  >AB353770.1.1740_U;tax=d:Eukaryota,k:TSAR,p:Alveolata-Dinoflagellata,c:Dinophyceae,...

VSEARCH then classifies the domain first. The supergroup now sits under the letter
k:, which for VSEARCH is only a position, the second rank. MetaFlux still writes it
as sg__ in the Supergroup column: 18S.yaml pairs the SINTAX letter k with
Supergroup and sg__ (sintax_rank_letters / sintax_rank_prefixes). Only the raw
VSEARCH output, 06.taxonomy/_sintax_raw.tsv, shows the supergroup as k:<supergroup>.
Everything after the first two values is copied character for character.

WHAT IT REFUSES (the script stops with a message rather than guess)
  - a header without a ";tax=" part;
  - a header whose letters are not exactly k,d,p,c,o,f,g,s in this order, the
    shape of PR2 5.1.1. A PR2 release that adds or drops a rank, or that already
    lists the domain under d:, stops here, and the 18S marker pack and this script
    have to be checked against the new file;
  - a k: value that is not one of PR2's domain values (Eukaryota, Bacteria,
    Archaea, or the organelle tags Eukaryota:plas, :mito, :nucl, :apic, :chro).
Together these make a second swap impossible: a file this script already wrote
starts with d:, not k:, so feeding it back in stops at the letter check.

INPUT  : input.fasta  = PR2's UTAX file (symbol pr2_utax in 18S.yaml): the download,
         gzipped as PR2 publishes it (rule fetch_pr2_utax), or the file a run config
         pins with references.pr2.utax, which may also be plain text. Before this
         step existed VSEARCH read that file directly, and VSEARCH reads either form,
         so this script does too (open_utax below).
OUTPUT : output.fasta = the domain-first SINTAX file (symbol pr2_sintax in 18S.yaml),
         read next by rule assign_taxonomy (VSEARCH --sintax, via
         80c_parse_sintax.py) when an 18S run uses taxonomy.method: sintax. The rdp
         path never reads it: it uses PR2's separate DADA2 file.
"""
import gzip
import sys
from collections import Counter
from pathlib import Path

sm = snakemake  # noqa: F821 (injected by Snakemake)

# Paths come from rule convert_pr2_sintax in 10_refdb.smk:
#   input.fasta  = PR2's UTAX FASTA: gzipped as downloaded, or plain when the run
#                  config pins an uncompressed copy (references.pr2.utax)
#   output.fasta = the domain-first SINTAX FASTA this script writes (gzipped)
#   log[0]       = a plain-text log for a human to check after the run; nothing
#                  downstream reads it
in_path  = Path(sm.input.fasta)
out_path = Path(sm.output.fasta)
log_path = Path(sm.log[0])

# The letters of every PR2 5.1.1 UTAX header, in the order PR2 writes them:
# k = Domain, d = Supergroup, p = Division-Subdivision, then class to species.
EXPECTED_LETTERS = ["k", "d", "p", "c", "o", "f", "g", "s"]

# The values PR2 puts in its domain field. Besides the three domains it tags
# organelle and endosymbiont-derived 18S/16S genes with a suffix: plastid (plas),
# mitochondrion (mito), nucleomorph (nucl), apicoplast (apic) and the Paulinella
# chromatophore (chro).
DOMAIN_VALUES = {
    "Eukaryota", "Bacteria", "Archaea",
    "Eukaryota:plas", "Eukaryota:mito", "Eukaryota:nucl",
    "Eukaryota:apic", "Eukaryota:chro",
}


def open_utax(path):
    """Open PR2's UTAX file for reading as text, whether it is gzipped or not.
    The choice is made by the file's first two bytes (every gzip file starts with
    1f 8b), not by its name, as VSEARCH and check_reference_ranks.py do: a pinned
    file may be plain even when its name ends in .gz, and the other way round."""
    with open(path, "rb") as raw:
        first_bytes = raw.read(2)
    if first_bytes == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path)


def stop(message, log, fo):
    """Write the problem to the log, remove the half-written output, and stop the
    rule with the same message, so no truncated reference is ever left behind."""
    log.write(f"[pr2_utax_domain_first] ERROR: {message}\n")
    log.close()
    fo.close()
    out_path.unlink(missing_ok=True)
    sys.exit(f"[pr2_utax_domain_first] {message}")


# Running totals for the summary at the end of the log, so a human can check the
# rebuild (every header rewritten, no sequence lost, the domain counts as expected).
n_headers = 0
n_seq_lines = 0
records_per_domain = Counter()

out_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

# compresslevel=6 is gzip's usual level. Python's default, 9, makes this file only
# marginally smaller and takes several times longer to write.
with (open_utax(in_path) as fi,
      gzip.open(out_path, "wt", compresslevel=6) as fo,
      log_path.open("w") as log):

    # Walk the FASTA one line at a time. A header line (">...") is rewritten below;
    # every other line is sequence data and is copied through as it is.
    for line_no, line in enumerate(fi, 1):
        if not line.startswith(">"):
            fo.write(line)
            n_seq_lines += 1
            continue

        header = line.rstrip("\n")
        if ";tax=" not in header:
            stop(f"{in_path}, line {line_no}: no ';tax=' part in the header: {header!r}",
                 log, fo)

        # ">AB353770.1.1740_U" and "k:Eukaryota,d:TSAR,p:...,s:..."
        record_id, tax = header.split(";tax=", 1)
        tokens = tax.split(",")

        # Each token is "<letter>:<value>". Split at the FIRST colon only, because
        # PR2's values can contain one themselves (Eukaryota:plas, TSAR:apic).
        letters = []
        values = []
        for token in tokens:
            if ":" not in token:
                stop(f"{in_path}, line {line_no}: the part {token!r} has no "
                     f"'<letter>:<value>' shape: {header!r}", log, fo)
            letter, value = token.split(":", 1)
            letters.append(letter)
            values.append(value)

        if letters != EXPECTED_LETTERS:
            stop(f"{in_path}, line {line_no}: the letters are {','.join(letters)}, not "
                 f"{','.join(EXPECTED_LETTERS)} as in PR2 5.1.1. PR2 may have changed "
                 f"its UTAX format, or this file was converted already; check the new "
                 f"file against workflow/markers/18S.yaml and this script: {header!r}",
                 log, fo)

        domain = values[0]
        supergroup = values[1]
        if domain not in DOMAIN_VALUES:
            stop(f"{in_path}, line {line_no}: k: holds {domain!r}, which is not a PR2 "
                 f"domain value ({', '.join(sorted(DOMAIN_VALUES))}): {header!r}",
                 log, fo)
        records_per_domain[domain] += 1

        # The swap: domain under d:, supergroup under k:. The tokens from the third
        # one on (p: to s:) are reused exactly as they were read.
        new_tokens = [f"d:{domain}", f"k:{supergroup}"] + tokens[2:]
        fo.write(f"{record_id};tax={','.join(new_tokens)}\n")
        n_headers += 1

    # Final tally for whoever checks the log after a refdb rebuild; not parsed by
    # Snakemake or any downstream rule.
    log.write(
        f"[pr2_utax_domain_first] {n_headers} headers rewritten (domain moved to d:, "
        f"supergroup to k:), {n_seq_lines} sequence lines passed through\n"
    )
    log.write("[pr2_utax_domain_first] records per domain value:\n")
    for domain, count in records_per_domain.most_common():
        log.write(f"[pr2_utax_domain_first]   {domain}: {count}\n")
    log.write(f"[pr2_utax_domain_first] Output: {out_path}\n")
