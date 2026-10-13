#!/usr/bin/env python3
"""Reference check (workflow/scripts/refdb/): does the taxonomy reference hold the
ranks its marker pack declares? Stop the run before classification when it does not.

Unlike the other scripts in this folder, this one changes no file. It reads the
reference the classifier is about to use and writes a small report.

WHY THIS STEP EXISTS
A classifier never checks what a rank MEANS. DADA2's assignTaxonomy names the
columns of its result from the marker pack's tax_levels, one after the other, and
VSEARCH reads the SINTAX letters (d:, k:, p:, ...) as a fixed top-down order. If a
reference has one rank more or less than the pack expects, or two ranks swapped,
every name lands in the wrong column and nothing errors. This has happened twice:
SILVA 144 added a kingdom rank between domain and phylum, and PR2's UTAX file puts
its domain under the letter k:, below the supergroup.

Counting ranks or letters is not enough. PR2's own UTAX file uses the same eight
letters as the corrected file MetaFlux builds, only in another order, and a SILVA
138.2 SINTAX file (no kingdom) read under the SILVA 144 model only lacks the letter
k:, which the SINTAX parser (80c_parse_sintax.py) skips without a word. So this
script also checks the NAMES in the ranks where such a shift shows: the domain rank
must hold domain names, SILVA 144's kingdom rank kingdom names, and so on. Each
marker pack declares these rules in reference_values.

A stop here means that the run config pins a wrong file, or that the reference
itself changed upstream (a new release with another rank). The fix is then the
right file, or a new release block in the marker pack.

WHAT IS CHECKED, on every header of the file
  DADA2 training set (method rdp):
    - the lineage is read as assignTaxonomy reads it: the whole header, or, for the
      UNITE general release, the 5th '|' field;
    - no rank inside a lineage is empty;
    - no lineage is deeper than the rank model (minus the Species rank when
      addSpecies fills it, as for 16S), and the deepest lineage is exactly that
      deep. Shallower lineages are normal when the header ends with ';', as
      SILVA's do: assignTaxonomy then pads them with empty ranks. Without that
      ';' it glues its padding onto the last name and its table of ranks no
      longer lines up, so such a lineage stops the run;
    - prefix_style embedded: every name starts with the pack's prefix for its rank,
      and a rank whose prefix is "" (gyrB's gene tag) carries no prefix at all;
      prefix_style bare: no name carries a rank prefix.
  SINTAX file (method sintax):
    - every letter is one of the pack's letters, and the letters of a header follow
      the pack's order top-down (a rank may be missing: UNITE leaves out ranks it
      has no name for);
    - every letter the pack declares occurs somewhere in the file, except the
      species letter when the species comes from addSpecies (16S: SILVA's SINTAX
      file stops at genus).
  Both:
    - reference_values: every name of a rank the pack gives a rule for (taken
      without its prefix) must pass that rule. A rule has up to three kinds of
      entry:
        one_of     - the names allowed in this rank;
        ends_with  - endings that are allowed (SILVA 144 kingdoms end in "ati");
        not_one_of - names that must never appear in this rank.
      A name fails when it is in not_one_of. Otherwise, when one_of or ends_with is
      given, it must be in one_of or end with one of the endings. When the run
      starts, 00_common.smk has already checked that the active pack's rules name
      ranks of this rank model and use only these three kinds, each a list of at
      least one name.

INPUT (rule check_taxonomy_reference, 10_refdb.smk):
  input.reference = the exact file the classifier reads (rdp: the DADA2 training
                    set; sintax: the SINTAX file), a path pinned in the run config
                    included. It comes from a download rule or a conversion rule in
                    10_refdb.smk, or is used as it is when it already exists.
  params          = ref_format ("dada2" or "sintax"); the rank model of the active
                    method (tax_levels, rank_letters, rank_prefixes); prefix_style;
                    species_from_addspecies; reference_values; label (how messages
                    name the check, e.g. "the 16S marker pack, release 144, method
                    sintax").
OUTPUT:
  output.report   = stats/taxonomy_reference_check.json, a record of what was
                    checked: records, ranks per record, letters seen (sintax) and
                    the names seen in each checked rank. Written only when every
                    check passes. Rule assign_taxonomy (80_taxonomy.smk) takes it as
                    an input, so classification cannot start on a reference that
                    failed; nothing reads its content.
  log[0]          = the same summary, then every problem found, or "OK".

Standard library only.
"""
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

sm = snakemake  # noqa: F821 (injected by Snakemake)

# ── Inputs and settings from the rule ─────────────────────────────────────────
reference_path = Path(sm.input.reference)
report_path    = Path(sm.output.report)
log_path       = Path(sm.log[0])
ref_format     = sm.params.ref_format           # "dada2" (rdp) or "sintax"
tax_levels     = list(sm.params.tax_levels)     # rank names, top-down
rank_letters   = list(sm.params.rank_letters)   # SINTAX letters, one per rank
rank_prefixes  = list(sm.params.rank_prefixes)  # d__, k__, ... one per rank
prefix_style   = sm.params.prefix_style         # "bare" or "embedded"
species_from_addspecies = bool(sm.params.species_from_addspecies)
value_rules    = dict(sm.params.reference_values or {})   # rank name -> rule
label          = sm.params.label

# A rank prefix as it appears at the start of a name: d__, k__, p__, sg__, sbd__ ...
PREFIX_PATTERN = re.compile(r"^[a-z]{1,3}__")

# How much goes into the stop message, so that it stays readable on a terminal.
# The log lists every problem.
MAX_PROBLEMS_IN_MESSAGE = 10
MAX_EXAMPLE_HEADERS     = 3     # headers quoted per problem
MAX_NAMES_PER_PROBLEM   = 8     # offending names quoted per problem
MAX_NAMES_IN_REPORT     = 60    # names listed per checked rank in the JSON report

report_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)


# ── Problems found, in the order first seen ───────────────────────────────────
# One entry per kind of problem (its description): how many headers show it, which
# names or letters caused it, and the first few headers as examples. A problem of
# the whole file (the deepest lineage, a letter that never occurs) has no header;
# its count stays 0.
problems = {}
n_records = 0      # headers read so far; also the number of the current header


def note_problem(description, header=None, found=None):
    """Remember one header (or the whole file, header=None) that shows a problem,
    and the name or letter that caused it, if there is one. A header that shows the
    same problem twice (two empty ranks, say) is counted and quoted once."""
    if description not in problems:
        problems[description] = {"headers": 0, "found": Counter(), "examples": [],
                                 "last_header": 0}
    entry = problems[description]
    if found is not None:
        entry["found"][found] += 1
    if header is None:
        return
    if entry["last_header"] == n_records:
        return
    entry["last_header"] = n_records
    entry["headers"] += 1
    if len(entry["examples"]) < MAX_EXAMPLE_HEADERS:
        entry["examples"].append(header)


def value_passes(name, rule):
    """True when a name passes its rank's reference_values rule (see the module
    docstring). Example, SILVA 144's Kingdom rule {ends_with: [ati], one_of:
    [Incertae Sedis]}: "Bacillati" and "Incertae Sedis" pass; "Bacillota" (a phylum)
    and "Pseudomonadati--other" fail."""
    forbidden_names = rule.get("not_one_of") or []
    if name in forbidden_names:
        return False
    allowed_names   = rule.get("one_of") or []
    allowed_endings = rule.get("ends_with") or []
    if not allowed_names and not allowed_endings:
        return True
    if name in allowed_names:
        return True
    for ending in allowed_endings:
        if name.endswith(ending):
            return True
    return False


def describe_rule(rule):
    """A rule written out for messages, e.g. "one_of [Incertae Sedis]; ends_with [ati]"."""
    parts = []
    for kind in ("one_of", "ends_with", "not_one_of"):
        entries = rule.get(kind) or []
        if entries:
            parts.append(f"{kind} [{', '.join(entries)}]")
    return "; ".join(parts)


def read_headers(path):
    """Yield every FASTA header of the file, without its '>' and without the spaces
    and line ending around it (assignTaxonomy trims them too). A byte that is not
    valid text is replaced, never a crash.
    The file may be gzipped. DADA2 and VSEARCH read a gzipped file whatever its name,
    so here too this is decided by the file's first two bytes (every gzip file starts
    with 1f 8b), not by a ".gz" at the end of its name."""
    with open(path, "rb") as raw:
        first_bytes = raw.read(2)
    if first_bytes == b"\x1f\x8b":
        handle = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    else:
        handle = open(path, encoding="utf-8", errors="replace")
    with handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip()


def is_unite_layout(path):
    """True when assignTaxonomy treats the file as the UNITE general release: each
    of its first ten headers contains "FU|refs" or "FU|reps" (dada2's own test, also
    used by 80a_assign_taxonomy.R). It then reads the lineage from the 5th '|' field
    of every header instead of the whole header."""
    first_headers = []
    for header in read_headers(path):
        first_headers.append(header)
        if len(first_headers) == 10:
            break
    if len(first_headers) < 10:
        return False
    for header in first_headers:
        if not re.search(r"FU\|re[pf]s", header):
            return False
    return True


def strip_pack_prefix(name):
    """A SINTAX name without a leading rank prefix of this pack, the way the SINTAX
    parser (80c_parse_sintax.py) reads it: a file that writes g:g__Fusarium ends up
    as g__Fusarium in the tables, so the rule is tested on Fusarium."""
    for prefix in rank_prefixes:
        if prefix and name.startswith(prefix):
            return name[len(prefix):]
    return name


# ── Rules by rank name -> by position ─────────────────────────────────────────
# The pack names its rules by rank (Domain, Kingdom ...); the checks below walk a
# lineage by position. 00_common.smk has already stopped the run if a rule names a
# rank this model does not have; the test is repeated here so the script never
# silently skips a rule when it is run with other settings.
for rank_name in value_rules:
    if rank_name not in tax_levels:
        sys.exit(
            f"[check_taxonomy_reference] reference_values has a rule for the rank "
            f"{rank_name!r}, which is not in the rank model of {label} "
            f"({', '.join(tax_levels)}). Fix the marker pack."
        )

# ── Running totals for the report (n_records is counted above) ────────────────
depth_counts = Counter()       # ranks per record -> number of records
letters_seen = Counter()       # SINTAX letter -> number of times it occurs
names_seen   = {}              # checked rank -> Counter of its names
for rank_name in value_rules:
    names_seen[rank_name] = Counter()


def check_name(rank_name, name, header, letter=None):
    """Count a name of a checked rank and test it against that rank's rule."""
    names_seen[rank_name][name] += 1
    rule = value_rules[rank_name]
    if value_passes(name, rule):
        return
    if letter is None:
        where = rank_name
    else:
        where = f"{rank_name} (letter {letter}:)"
    note_problem(f"{where} holds a name its rule does not allow ({describe_rule(rule)})",
                 header, name)


# ── DADA2 training set (method rdp) ───────────────────────────────────────────
if ref_format == "dada2":
    # For 16S the Species rank is filled later by addSpecies from its own file, so
    # the training set stops one rank short of the rank model; every other marker's
    # training set covers all of its ranks.
    expected_depth = len(tax_levels)
    species_note = ""
    if species_from_addspecies:
        expected_depth = len(tax_levels) - 1
        species_note = " (the Species rank comes from addSpecies)"
    unite_layout = is_unite_layout(reference_path)

    for header in read_headers(reference_path):
        n_records += 1

        # The lineage, as assignTaxonomy reads it.
        if unite_layout:
            fields = header.split("|")
            if len(fields) < 5:
                note_problem("the header has no 5th '|' field, where assignTaxonomy "
                             "reads a UNITE lineage", header)
                continue
            lineage = fields[4]
        else:
            lineage = header
        # assignTaxonomy ignores one ';' at the end of a lineage (R's strsplit drops
        # a last empty field), so a SILVA header "Bacteria;...;Pseudomonas;" has as
        # many ranks as names. Whether the ';' was there is remembered, because it
        # decides how assignTaxonomy pads a short lineage (see below).
        ends_with_semicolon = lineage.endswith(";")
        if ends_with_semicolon:
            lineage = lineage[:-1]
        names = lineage.split(";")
        depth_counts[len(names)] += 1

        if len(names) > expected_depth:
            note_problem(f"a lineage deeper than the {expected_depth} ranks the rank "
                         f"model expects in this file{species_note}", header)
            continue

        # A lineage shallower than the rank model is normal (SILVA has lineages that
        # end above genus), but only when it ends with ';'. assignTaxonomy brings
        # every lineage to the depth of the deepest one by appending
        # "_DADA2_UNSPECIFIED;" once per missing rank. After a closing ';' each piece
        # is a rank of its own, which DADA2 reports as NA. Without the ';' the first
        # piece is glued onto the last name (p__Bacillota_DADA2_UNSPECIFIED), and the
        # lineage stays one rank short. assignTaxonomy then cuts the lineages, read one
        # after the other, into rows of equal length, so the ranks of every lineage
        # that follows this one in the file are shifted, and the bootstrap support of
        # those genera is counted on the wrong ranks, without an error.
        if len(names) < expected_depth and not ends_with_semicolon:
            note_problem(f"a lineage shallower than the {expected_depth} ranks the rank "
                         f"model expects{species_note}, without a ';' at its end: "
                         f"assignTaxonomy would glue its padding (_DADA2_UNSPECIFIED) onto "
                         f"the last name, which shifts the ranks of every lineage that "
                         f"follows it in the file", header)

        # Each name in turn: is it empty, does it carry the right prefix, and, for
        # the ranks the pack gives a rule for, is it a name that rank may hold?
        for position, name in enumerate(names):
            rank_name = tax_levels[position]
            prefix = rank_prefixes[position]
            if name == "":
                note_problem("an empty rank inside the lineage", header)
                continue
            if prefix_style == "embedded":
                if prefix == "":
                    # A rank without a prefix in this pack (gyrB's gene tag): its
                    # names must carry none either.
                    if PREFIX_PATTERN.match(name):
                        note_problem(f"{rank_name} (rank {position + 1}) carries a rank "
                                     f"prefix, but this rank has none in the pack",
                                     header, name)
                        continue
                    bare_name = name
                else:
                    if not name.startswith(prefix):
                        note_problem(f"{rank_name} (rank {position + 1}) does not start "
                                     f"with {prefix}", header, name)
                        continue
                    bare_name = name[len(prefix):]
            else:
                # Bare style: MetaFlux adds the prefixes itself, so a prefix already
                # in the file would be written twice (d__d__Bacteria).
                if PREFIX_PATTERN.match(name):
                    note_problem(f"{rank_name} (rank {position + 1}) carries a rank "
                                 f"prefix, but the pack declares this reference bare "
                                 f"(prefix_style: bare)", header, name)
                    continue
                bare_name = name
            if rank_name in value_rules:
                check_name(rank_name, bare_name, header)

    # assignTaxonomy makes one column per rank of the deepest lineage in the file.
    if n_records > 0:
        deepest = max(depth_counts)
        if deepest != expected_depth:
            note_problem(f"the deepest lineage has {deepest} ranks, but the rank model "
                         f"expects {expected_depth} in this file{species_note}")

# ── SINTAX file (method sintax) ───────────────────────────────────────────────
elif ref_format == "sintax":
    position_of_letter = {}
    for position, letter in enumerate(rank_letters):
        position_of_letter[letter] = position
    letters_text = ",".join(rank_letters)

    for header in read_headers(reference_path):
        n_records += 1
        if ";tax=" not in header:
            note_problem("the header has no ';tax=' field", header)
            continue
        # VSEARCH reads the tax= field up to the next ';' of the header.
        tax = header.split(";tax=", 1)[1]
        tax = tax.split(";", 1)[0]
        fields = tax.split(",")
        depth_counts[len(fields)] += 1

        # Each "letter:name" field in turn. last_position is the rank of the field
        # before, so a letter that is not below it shows a header out of order.
        last_position = -1
        for field in fields:
            if ":" not in field:
                note_problem("a taxonomy field that is not 'letter:name'", header, field)
                continue
            letter, name = field.split(":", 1)
            if name == "":
                note_problem("a taxonomy field with a letter but no name", header, field)
                continue
            letters_seen[letter] += 1
            if letter not in position_of_letter:
                note_problem(f"a letter that is not one of this rank model's letters "
                             f"({letters_text}); the SINTAX parser would drop that rank "
                             f"without a word", header, f"{letter}:")
                continue
            position = position_of_letter[letter]
            if position <= last_position:
                note_problem(f"letters out of the rank model's top-down order "
                             f"({letters_text})", header)
            last_position = position
            rank_name = tax_levels[position]
            if rank_name in value_rules:
                check_name(rank_name, strip_pack_prefix(name), header, letter)

    # Every declared letter must occur somewhere. The species letter is the last one
    # in every pack; under addSpecies (16S) the SINTAX file has no species.
    last_position = len(rank_letters) - 1
    for position, letter in enumerate(rank_letters):
        if position == last_position and species_from_addspecies:
            continue
        if letters_seen[letter] == 0:
            note_problem(f"the letter {letter}: ({tax_levels[position]}) never occurs, "
                         f"although the rank model declares it")

else:
    sys.exit(f"[check_taxonomy_reference] unknown ref_format {ref_format!r} "
             f"(expected 'dada2' or 'sintax')")

if n_records == 0:
    note_problem("the file holds no FASTA header")


# ── Summary of what was checked (the report, and the top of the log) ──────────
ranks_per_record = {}
for depth in sorted(depth_counts):
    ranks_per_record[str(depth)] = depth_counts[depth]

names_report = {}
for rank_name, counts in names_seen.items():
    most_common = {}
    for name, count in counts.most_common(MAX_NAMES_IN_REPORT):
        most_common[name] = count
    names_report[rank_name] = {"distinct_names": len(counts), "most_common": most_common}

summary = {
    "reference": str(reference_path),
    "checked_against": label,
    "format": ref_format,
    "records": n_records,
    "ranks_per_record": ranks_per_record,
    "rank_model": tax_levels,
    "rank_prefixes": rank_prefixes,
    "prefix_style": prefix_style,
    "reference_values": value_rules,
    "names_seen": names_report,
}
if ref_format == "dada2":
    summary["expected_depth"] = expected_depth
    summary["unite_layout"] = unite_layout
else:
    summary["rank_letters"] = rank_letters
    summary["letters_seen"] = dict(letters_seen)


def describe_problem(description, entry):
    """The lines that describe one problem: how often, which names, example headers."""
    if entry["headers"] > 0:
        line = f"  - {description}: {entry['headers']:,} header(s)"
    else:
        line = f"  - {description}"
    found = entry["found"]
    if found:
        quoted = []
        for name, count in found.most_common(MAX_NAMES_PER_PROBLEM):
            quoted.append(f"{name!r} ({count:,})")
        line += "; found: " + ", ".join(quoted)
        if len(found) > MAX_NAMES_PER_PROBLEM:
            line += f" and {len(found) - MAX_NAMES_PER_PROBLEM:,} more"
    lines = [line]
    for example in entry["examples"]:
        lines.append(f"      e.g. >{example[:200]}")
    return lines


with log_path.open("w") as log:
    log.write(json.dumps(summary, indent=1) + "\n")

    if problems:
        heading = (f"[check_taxonomy_reference] {reference_path} does not hold the ranks "
                   f"declared by {label} ({', '.join(tax_levels)}):")
        advice = ("  The classifier would put names on the wrong ranks, so the run stops "
                  "before any classification. Either the run config points at a wrong "
                  "file (references:, and for 16S the SILVA release chosen with "
                  "references.silva.release), or the reference changed upstream and its "
                  "marker pack needs a new release block (see the reference_values notes "
                  "in workflow/markers/16S.yaml).")

        # The log gets every problem.
        log.write(heading + "\n")
        for description, entry in problems.items():
            for line in describe_problem(description, entry):
                log.write(line + "\n")
        log.write(advice + "\n")

        # The stop message gets the first MAX_PROBLEMS_IN_MESSAGE of them.
        message_lines = [heading]
        shown = 0
        for description, entry in problems.items():
            if shown == MAX_PROBLEMS_IN_MESSAGE:
                break
            message_lines.extend(describe_problem(description, entry))
            shown += 1
        not_shown = len(problems) - shown
        if not_shown > 0:
            message_lines.append(f"  ... and {not_shown} more kind(s) of problem, "
                                 f"listed in {log_path}")
        message_lines.append(advice)
        sys.exit("\n".join(message_lines))

    # Every check passed: write the report that lets assign_taxonomy start.
    summary["result"] = "passed"
    report_path.write_text(json.dumps(summary, indent=1) + "\n")
    log.write(f"[check_taxonomy_reference] OK: {n_records:,} records of "
              f"{reference_path} match {label}\n")
