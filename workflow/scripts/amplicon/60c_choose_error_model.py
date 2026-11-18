#!/usr/bin/env python3
"""Choose the error function DADA2 learns each error model with, from the quality
values the filtered reads carry.

What this is for: DADA2 learns one error model per read direction for the whole run
(rule learn_errors, next in the chain): for each quality score, how often each base is
read as each other base. In every round of that learning an "error function" turns the
misreads counted at each quality score into error rates. Which function fits depends
on the kind of quality values the reads carry:
  - continuous: a MiSeq writes nearly every value from 2 to 38 (filtered MiSeq reads
    carry 28-32 distinct values above Q2, measured on 56 files of six test sets).
    DADA2's default loess fit draws a smooth curve through the rates; it is the fit
    MetaFlux used for every run up to v2.6.1.
  - binned: a NovaSeq, a NextSeq 1000/2000, a MiSeq i100 and other recent instruments
    write only a few values ("bins", for example 9, 20 and 34). DADA2's author
    recommends a binned fit for them (makeBinnedQualErrfun), which joins the rates
    counted at those values with straight lines.
The values depend on the instrument, its control software and even its kit, so they
are never taken from a table: dada_filter counted them in the first 100,000 filtered
reads of every sample (the reads learnErrors learns from), and this script pools those
counts and decides. The function names are those of the registry in
60_error_functions.R; SUITED_CLASS and ON_OTHER_CLASS below hold one entry per name.

Inputs (see rule choose_error_model in 60_dada2.smk):
  quality_values : stats/dada2/{sample}.quality_values.json of every sample (written by
                   dada_filter, 60b_dada_filter.R): per read direction, the reads
                   scanned, the distinct values and the bases at each value.
  params         : error_model_requested, the setting amplicon.dada2.learn_errors.
                   error_model as checked at start-up (auto, loess, binned, or the
                   experimental custom or scam); trunc_q, only to say in a message that a
                   raised trunc_q removed the lower quality values.

Decision, for R1 and R2 separately (decide_direction below):
  1. Voters: samples with at least 1,000 reads scanned. A blank with a few dozen reads
     can lack values by chance, so smaller samples are listed but do not vote and
     never stop anything. If no sample has 1,000 reads, every sample votes (warning).
  2. The quality values of the run are the union of the voters' values, Q2 and below
     included (they are the anchors the binned fit is handed); the class is read from
     the number of distinct values above Q2 (Q0-Q2 are the no-call values of every
     Illumina table, and trunc_q: 2 removes them): fewer than 3 = no usable quality
     information; 3-10 = binned, when no two of them are neighbours (29 and 30, say);
     3-10 that include neighbours, and 11-19 = undetermined (seen on no instrument: no
     instrument table has neighbouring values, while MiSeq reads under a trunc_q of 28
     to 30 keep only the 8 to 10 values above it, all neighbours); 20 or more =
     continuous. Neighbours that only two binned tables make together (a NovaSeq X
     writes 40 on some versions of its software and 41 on others) leave the run binned,
     with two tables (C6 below).
  3. Stops, in this order (the first that applies ends the decision of the direction):
     one voter without quality information among voters that have it, or every voter
     without it while their values together reach three (C3, any setting); no quality
     information at all (C2, any setting); a function set in the config on the kind
     of values where that stops, binned on continuous values (C7), checked before the
     warnings that follow so that a run about to stop prints no "going on ... as set"
     warning; voters of different classes (C5 under auto, a warning when a function is
     set); binned voters with two different sets of values (C6 under auto, a warning
     when set); an undetermined class (C4 under auto, a warning when set). So loess set
     in the config turns C4-C6 into warnings and goes on; binned does the same for C4
     and C6, but on binned and continuous voters together their values make the run
     continuous, and binned stops with C7.
  4. The function: auto takes loess for continuous and binned for binned values. A
     function set in the config is used; loess set on binned values warns and goes on
     (W1; it is the fit of v2.6.1 and earlier).
  5. A non-voter carrying values the voters do not gets a note (W4): its extra values
     are neither used for the choice nor as anchors.
  error_model: custom (EXPERIMENTAL): the user's own error function, named in
  amplicon.dada2.learn_errors.custom_error_function, for both directions, decided by
  the setting. The class of the values is still read and written down (the record and
  the report show it), and the stops that describe a problem of the data itself still
  apply, because no error function can fix them: no usable quality information (C2),
  one sample without it (C3). What only serves to choose between loess and binned does
  not apply: mixed classes, two binned sets and an undetermined class give the
  warnings any function set in the config gets (the warnings that replace C5, C6 and
  C4), and the function-on-the-other-class rule (C7, W1) and the note on non-voters
  (W4, about the choice and the anchors of the binned fit) are skipped.
Both directions are decided before anything stops, so every stop text is printed at
once. Warnings and notes go to the log and the terminal as they come; a stop text goes
to the log and then, once, to the terminal as the text of the SystemExit (Python
prints it itself, so it is not printed before as well), and nothing is written, so
learn_errors cannot start.

Outputs:
  stats/dada2/error_model_choice.json : per read direction the class, the function,
      how it was decided, the quality values (the anchors of a binned fit), the bases
      at each value summed over the voters, and who voted; the warnings; the values of
      every sample. Read by learn_errors (which function to learn with, on which
      values); error_model_report takes it as an input, so that the report is drawn
      again when the choice changes. The same inputs give the same bytes, so a rerun
      that changes nothing makes no later step run again.
  stats/dada2/quality_values.tsv : one row per sample and read direction, for people to
      read; nothing in the workflow reads it.
"""
import json
import sys
from collections import Counter
from pathlib import Path

# ── The thresholds of the decision, each with the measurement behind it ──────────
# Q0-Q2: the no-call values of every Illumina quality table; the default trunc_q: 2
# cuts each read at its first base at or below Q2, so they rarely reach the filtered
# reads, and when they do (trunc_q 0 or 1) they say nothing about the instrument.
NO_CALL_MAX = 2
# Fewer distinct values above Q2 than this: no usable quality information. No Illumina
# table has fewer than three values above Q2; one or two is the mark of a file whose
# values were rewritten (ENA copies and SRA Lite files carry Q30, and Q3 for rejected
# bases, at every base).
MIN_QUALITY_VALUES = 3
# 3 to 10 values above Q2, no two of them neighbours: binned. Instrument tables carry 3
# to 7 values above Q2 (NovaSeq 6000: 11, 25, 37; NextSeq 500/550: 14, 21, 27, 32, 36;
# HiSeq 2500 binned: up to 7), and binned files on disk carried 3-5. No table has two
# neighbouring values (the NovaSeq X table ends in 40, or in 41 on some versions of its
# software, never in both). 3 to 10 values that include neighbours are continuous values
# with the lower ones cut away (MiSeq reads under trunc_q 28-30 keep only the 8-10
# values above it), and count as undetermined (see classify).
MAX_BINNED_VALUES = 10
# 20 or more values above Q2: continuous. Filtered MiSeq reads carry 28-32 (measured on
# 56 files of the six MiSeq test sets; the fewest, 28, in gyrB R1). 11-19 is seen on no
# instrument, and MetaFlux does not guess there.
MIN_CONTINUOUS_VALUES = 20
# Smaller samples are listed but do not vote, and do not decide anything: a blank with
# a few dozen reads can lack values by chance. MiSeq samples with about 5,000 filtered
# reads still showed all 32 values, and every value of a file appears within its first
# 1,000 reads.
MIN_READS_TO_VOTE = 1000
# The function auto takes for each class.
AUTO_FUNCTION = {"continuous": "loess", "binned": "binned"}
# One entry per function of the registry (60_error_functions.R): the kind of quality
# values it is made for, and what happens when the config sets it on the other kind.
# loess on binned values: a warning (it is what MetaFlux did up to v2.6.1, and a user
# may want to keep it); binned on continuous values: a stop (it would join some 30
# counted rates with straight lines instead of fitting a curve).
SUITED_CLASS = {"loess": "continuous", "binned": "binned", "scam": "both"}
ON_OTHER_CLASS = {"loess": "warn", "binned": "stop", "scam": "none"}
# The setting's word for the user's own error function (EXPERIMENTAL; error_model:
# custom). Not a registry entry, so it has no entry in the two tables above: it is made
# for no particular kind of quality values, and learn_errors adds it to the registry
# when it runs, from the file named in amplicon.dada2.learn_errors.custom_error_function.
CUSTOM_FUNCTION = "custom"

# A stop text lists at most this many sample rows; the log has them all.
MAX_ROWS_IN_STOP = 20
DIRECTIONS = ("R1", "R2")


# ── Reading the counts ────────────────────────────────────────────────────────────

def read_quality_values(paths):
    """Read the quality values of every sample from its quality_values.json.

    Parameters
    ----------
    paths : list of str, the stats/dada2/{sample}.quality_values.json files, as
        dada_filter writes them: {"sample": name, "scan_limit": 100000, "R1":
        {"reads_scanned": n, "values": [...], "bases_per_value": {"9": n, ...}}, "R2":
        {...}}.

    Returns
    -------
    A list sorted by sample name, one dict per sample: {"sample": name, "R1":
    {"reads_scanned": int, "values": [int, ...] ascending, "bases_per_value": {int:
    int}}, "R2": {...}}.
    """
    samples = []
    for path in paths:
        data = json.loads(Path(path).read_text())
        sample = {"sample": data["sample"]}
        for direction in DIRECTIONS:
            part = data[direction]
            bases = {}
            for value, count in part["bases_per_value"].items():
                bases[int(value)] = int(count)
            sample[direction] = {
                "reads_scanned": int(part["reads_scanned"]),
                "values": sorted(int(value) for value in part["values"]),
                "bases_per_value": bases,
            }
        samples.append(sample)
    samples.sort(key=lambda sample: sample["sample"])
    return samples


def above_q2(values):
    """The values above the no-call values Q0-Q2, ascending."""
    return sorted(value for value in set(values) if value > NO_CALL_MAX)


def first_neighbours(values):
    """The first two neighbouring values of a set of quality values, such as (29, 30),
    or None when no two values are neighbours (one apart).

    No instrument table has neighbouring values (see MAX_BINNED_VALUES), while
    continuous values have them everywhere, also when a raised trunc_q has left only the
    few values above it."""
    ordered = sorted(set(values))
    for lower, upper in zip(ordered, ordered[1:]):
        if upper - lower == 1:
            return (lower, upper)
    return None


def classify(values_above_q2):
    """Say what kind of quality values a set of values above Q2 is.

    Parameters
    ----------
    values_above_q2 : iterable of int, the distinct quality values above Q2.

    Returns
    -------
    str: "none" (fewer than 3 values: no usable quality information), "binned" (3-10,
    no two of them neighbours), "undetermined" (3-10 that include neighbours, or 11-19)
    or "continuous" (20 or more).
    """
    values = sorted(set(values_above_q2))
    count = len(values)
    if count < MIN_QUALITY_VALUES:
        return "none"
    if count <= MAX_BINNED_VALUES:
        # Neighbours, such as 29 and 30, are no binned table but continuous values with
        # the lower ones cut away: MiSeq reads under trunc_q 28-30 keep 8-10 such values,
        # and the binned fit would join their rates with straight lines without a word.
        if first_neighbours(values) is not None:
            return "undetermined"
        return "binned"
    if count < MIN_CONTINUOUS_VALUES:
        return "undetermined"
    return "continuous"


# ── Texts ─────────────────────────────────────────────────────────────────────────

def values_text(values):
    """Quality values for a message: "9, 20, 34", or "32 values, 7-38" when there are
    more than ten; "none" for no value."""
    values = sorted(set(values))
    if not values:
        return "none"
    if len(values) > MAX_BINNED_VALUES:
        return f"{len(values)} values, {values[0]}-{values[-1]}"
    return ", ".join(str(value) for value in values)


def names_text(names):
    """Sample names for a message, at most MAX_ROWS_IN_STOP of them, then how many more."""
    shown = ", ".join(names[:MAX_ROWS_IN_STOP])
    if len(names) > MAX_ROWS_IN_STOP:
        shown += f" and {len(names) - MAX_ROWS_IN_STOP} more"
    return shown


def falco_path(out, sample, direction):
    """The Falco report of a sample's filtered reads (rule falco, stage "filtered")."""
    return f"{out}/stats/falco/{sample}_{direction}_filtered/fastqc_report.html"


def sample_table(direction, samples, voters, limit=None, log_path=None):
    """The per-sample table of one read direction, as lines of text.

    Parameters
    ----------
    direction : "R1" or "R2".
    samples : from read_quality_values.
    voters : set of the names of the voting samples of this direction.
    limit : at most this many sample rows, then a line saying where the rest is
        (None: every row).
    log_path : this step's log, which holds every row; named in that last line. It is
        the absolute path Snakemake hands the job (the output folder is absolute), as
        every path of the stop texts: the user's shell would read a relative path
        against the MetaFlux folder, where MetaFlux is started, not against the output
        folder. Needed only with a limit.

    Returns
    -------
    list of str.
    """
    width = max([len("sample")] + [len(sample["sample"]) for sample in samples])
    lines = [f"  {direction}  {'sample'.ljust(width)}  reads scanned  votes  "
             "quality values above Q2"]
    for count, sample in enumerate(samples):
        if limit is not None and count == limit:
            lines.append(f"  ... and {len(samples) - limit} more in {log_path}")
            break
        part = sample[direction]
        votes = "yes" if sample["sample"] in voters else "no"
        lines.append(f"  {direction}  {sample['sample'].ljust(width)}  "
                     f"{str(part['reads_scanned']).rjust(13)}  {votes.ljust(5)}  "
                     f"{values_text(above_q2(part['values']))}")
    return lines


def trunc_q_sentence_c2(trunc_q):
    """The C2 sentence for a trunc_q above 2, or "" (see stop_c2)."""
    if trunc_q <= NO_CALL_MAX:
        return ""
    return (f" amplicon.dada2.trunc_q is {trunc_q}: filterAndTrim cut each read at its first "
            f"base at or below Q{trunc_q}, so the lower values are gone from the filtered "
            "reads; set it back to 2, the default, first.")


def stop_c2(direction, declared, trunc_q, out, table_lines):
    """C2: no usable quality information in the run (any setting)."""
    return (
        f"[choose_error_model] STOP: the filtered {direction} reads carry fewer than three "
        "distinct quality values above Q2 (values seen in the first 100,000 reads of each "
        f"voting sample: {values_text(declared)}). DADA2 cannot learn an error model without "
        "quality information, whatever error_model says."
        + trunc_q_sentence_c2(trunc_q)
        + " Reads downloaded from ENA, and SRA Lite files, often carry one made-up value (Q30, "
        "and Q3 for rejected bases) instead of the sequencer's own; fetch such reads from NCBI "
        "with sra-tools instead (see \"Are the quality values real?\" in "
        "docs/about/test-datasets.md). Each sample's filtered reads: "
        f"{falco_path(out, '<sample>', direction)}.\n" + "\n".join(table_lines)
    )


def stop_c3(direction, none_voters, others_values, out):
    """C3: one or more voters without quality information among voters that have it.

    others_values are the values above Q2 of the other voters, those with quality
    information. When there are none, every voter lacks quality information, yet their
    values together reach three (two files rewritten with different made-up values, for
    example 3, 30 and 3, 40), so the run is not classed "none" and C2 does not apply. The
    usual text would then compare the samples with "the other voting samples (none)",
    which do not exist, so that case has a text of its own."""
    if not others_values:
        described = [f"{name} ({part['reads_scanned']:,} reads scanned: "
                     f"{values_text(above_q2(part['values']))})" for name, part in none_voters]
        return (
            "[choose_error_model] STOP: no voting sample carries three or more distinct "
            f"quality values above Q2 in its {direction} reads: {names_text(described)}. "
            "DADA2 cannot learn an error model without quality information, whatever "
            "error_model says. Reads downloaded from ENA, and SRA Lite files, often carry "
            "made-up values (Q30, and Q3 for rejected bases) instead of the sequencer's own; "
            "fetch such reads from NCBI with sra-tools instead (see \"Are the quality values "
            "real?\" in docs/about/test-datasets.md). Their filtered reads: "
            f"{falco_path(out, '<sample>', direction)}."
        )
    if len(none_voters) == 1:
        name, part = none_voters[0]
        return (
            f"[choose_error_model] STOP: sample {name} ({part['reads_scanned']:,} {direction} "
            "reads scanned) carries fewer than three distinct quality values above Q2 "
            f"({values_text(above_q2(part['values']))}), unlike the other voting samples "
            f"({values_text(others_values)}). One sample whose quality values were rewritten "
            "(for example downloaded from ENA) would make the one error model of the run wrong "
            "for every sample. Fetch its reads again from NCBI with sra-tools, or leave the "
            f"sample out. Its filtered reads: {falco_path(out, name, direction)}."
        )
    described = [f"{name} ({part['reads_scanned']:,} reads scanned: "
                 f"{values_text(above_q2(part['values']))})" for name, part in none_voters]
    return (
        f"[choose_error_model] STOP: {len(none_voters)} samples carry fewer than three "
        f"distinct quality values above Q2 in their {direction} reads, unlike the other voting "
        f"samples ({values_text(others_values)}): {names_text(described)}. One sample whose "
        "quality values were rewritten (for example downloaded from ENA) would make the one "
        "error model of the run wrong for every sample. Fetch their reads again from NCBI with "
        "sra-tools, or leave these samples out. Their filtered reads: "
        f"{falco_path(out, '<sample>', direction)}."
    )


def stop_c4(direction, declared_above, trunc_q, out, table_lines):
    """C4: an undetermined kind of quality values under auto, in one of two forms.

    11-19 distinct values above Q2: seen on no instrument. 3-10 values that include
    neighbours (29 and 30, say): no instrument table has neighbouring values, but
    continuous values with the lower ones cut away do, as MiSeq reads under a raised
    trunc_q keep only the values above it; this form names trunc_q when it is above 2,
    and the likely causes otherwise. Both forms end with the two lines that set the
    function and with the per-sample table (table_lines)."""
    listed = ", ".join(str(value) for value in declared_above)
    falco = falco_path(out, "<sample>", direction)
    how_to_go_on = (
        " If you know what the data are, add one of these lines under amplicon: -> dada2: -> "
        "learn_errors: in your run config, inside the amplicon: block it already has (do not "
        "add a second amplicon: block; YAML keeps only the last one), and start MetaFlux again:\n"
        "    error_model: loess   # continuous quality values (MiSeq-like)\n"
        "    error_model: binned  # binned quality values\n" + "\n".join(table_lines)
    )
    neighbours = first_neighbours(declared_above)
    if len(declared_above) <= MAX_BINNED_VALUES and neighbours is not None:
        lower, upper = neighbours
        if trunc_q > NO_CALL_MAX:
            cause = (f" amplicon.dada2.trunc_q is {trunc_q}: filterAndTrim cut each read at "
                     f"its first base at or below Q{trunc_q}, so only the values above "
                     f"Q{trunc_q} remain. Set it back to 2, the default, and the choice sees "
                     "every value the reads carry.")
        else:
            cause = (" Reads trimmed at a quality threshold before they reached MetaFlux (by "
                     "the sequencing facility or another tool), or quality values rewritten "
                     "by a facility, are the likely causes.")
        return (
            f"[choose_error_model] STOP: the filtered {direction} reads carry "
            f"{len(declared_above)} distinct quality values above Q2 ({listed}), and some of "
            f"them are neighbours ({lower} and {upper}). MetaFlux reads 3 to 10 values as "
            "binned quality scores only when there are gaps between them, as in every "
            "instrument table (NovaSeq 6000: 11, 25, 37); neighbouring values are what "
            "continuous ones (MiSeq) look like once the lower values are gone, and MetaFlux "
            "does not guess."
            + cause
            + f" The per-sample values are below, and each sample's filtered reads are in {falco}."
            + how_to_go_on
        )
    trunc_q_sentence = ""
    if trunc_q > NO_CALL_MAX:
        trunc_q_sentence = (f" amplicon.dada2.trunc_q is {trunc_q}, so only values above "
                            f"Q{trunc_q} remain, which can make continuous values look few.")
    return (
        f"[choose_error_model] STOP: the filtered {direction} reads carry {len(declared_above)} "
        f"distinct quality values above Q2 ({listed}). MetaFlux reads 3 to 10 values as binned "
        "quality scores (NovaSeq, NextSeq, MiSeq i100 and other recent instruments) and 20 or "
        "more as continuous ones (MiSeq), has seen nothing in between, and does not guess. A "
        "mix of instruments, or quality values rewritten by a sequencing facility, are the "
        "likely causes; the per-sample values are below, and each sample's filtered reads are "
        f"in {falco}."
        + trunc_q_sentence
        + how_to_go_on
    )


def class_clauses(voters, voter_class):
    """The C5 description of each class among the voters: "3 samples carry binned
    values (9, 20, 34)", "1 carries continuous ones (32 values, 7-38)", ...; the verb
    agrees with the number of samples ("1 sample carries", "2 samples carry")."""
    clauses = []
    for kind in ("binned", "continuous", "undetermined"):
        members = [part for name, part in voters if voter_class[name] == kind]
        if not members:
            continue
        values = set()
        for part in members:
            values.update(above_q2(part["values"]))
        values = sorted(values)
        count = len(members)
        verb = "carries" if count == 1 else "carry"
        # The first clause names the samples ("1 sample", "3 samples"), the others give
        # only the number.
        if clauses:
            who = f"{count}"
        elif count == 1:
            who = "1 sample"
        else:
            who = f"{count} samples"
        if kind == "binned":
            clauses.append(f"{who} {verb} binned values ({values_text(values)})")
        elif kind == "continuous":
            clauses.append(f"{who} {verb} continuous ones ({len(values)} values, "
                           f"{values[0]}-{values[-1]})")
        else:
            clauses.append(f"{who} {verb} {len(values)} values, neither binned nor continuous "
                           f"({values_text(values)})")
    if len(clauses) == 1:
        return clauses[0]
    return ", ".join(clauses[:-1]) + " and " + clauses[-1]


def stop_c5(direction, voters, voter_class, out, table_lines):
    """C5: voters of different classes under auto.

    The text offers error_model: loess only. binned cannot pass the mix C5 describes in
    the usual case: a sample with continuous values carries 20 or more values on its own,
    so the voters' values together make the run continuous, and binned set in the config
    then stops with C7. Like C4, it says where the line goes (inside the amplicon: block
    the run config has: YAML keeps only the last of two amplicon: blocks), and like C2-C4
    where the Falco reports of the filtered reads are."""
    return (
        "[choose_error_model] STOP: the voting samples do not share one kind of quality values "
        f"in the {direction} reads: {class_clauses(voters, voter_class)}; they were probably "
        "sequenced on different instruments. DADA2 learns one error model per read direction "
        "for all samples of a run, so they cannot be denoised as one run. Run the groups in "
        "separate output folders, or, knowing that one model then serves both kinds of reads, "
        "add error_model: loess under amplicon: -> dada2: -> learn_errors: in your run config, "
        "inside the amplicon: block it already has (do not add a second amplicon: block; YAML "
        "keeps only the last one). Each sample's filtered reads: "
        f"{falco_path(out, '<sample>', direction)}.\n" + "\n".join(table_lines)
    )


def stop_c6(direction, common, extra_voters, out):
    """C6: binned voters with two different sets of values under auto. The verb agrees
    with the samples named ("B1 also carries 9", "S4, S5 also carry 9")."""
    extra_values = set()
    for name, extra in extra_voters:
        extra_values.update(extra)
    verb = "carries" if len(extra_voters) == 1 else "carry"
    return (
        "[choose_error_model] STOP: the voting samples carry two different sets of binned "
        f"quality values in the {direction} reads: most carry {values_text(common)}, but "
        f"{names_text([name for name, extra in extra_voters])} also {verb} "
        f"{values_text(extra_values)}. Two instruments, or two versions of the sequencer's "
        "software, made these files, and one error model for both would be fitted to neither. "
        "Run the groups in separate output folders. error_model: loess lets the run go on with "
        "one curve for all values. error_model: binned is not advised here: learnErrors learns "
        "from a few whole samples, and the values of a group it does not draw would be fitted "
        "on averaged qualities of the other group. Each sample's filtered reads: "
        f"{falco_path(out, '<sample>', direction)}."
    )


def stop_c7(direction, function, run_class, declared_above, only_together=False,
            table_lines=None):
    """C7: a function set on the other kind of quality values, where that stops.

    only_together : True when the voters' values reach the continuous class only
        together, because not every voter carries continuous values on its own (binned
        and continuous samples in one run, for example). The text then does not say "as a
        MiSeq writes them", which would describe none of the samples, names the mix
        instead, gives the advice of C5 (separate output folders, or loess) and ends with
        the per-sample table (table_lines)."""
    if function == "binned" and run_class == "continuous":
        count = len(declared_above)
        if only_together:
            return (
                "[choose_error_model] STOP: amplicon.dada2.learn_errors.error_model is binned, "
                f"but together the voting samples' filtered {direction} reads carry {count} "
                f"distinct quality values above Q2 ({declared_above[0]}-{declared_above[-1]}), "
                "which MetaFlux reads as continuous, although not every sample carries that "
                "many on its own: the samples were probably sequenced on different instruments "
                "(the per-sample values are below). The binned fit would join the error rates "
                f"of {count} values with straight lines instead of fitting a curve. Run the "
                "groups in separate output folders, or set error_model to loess, knowing that "
                "one model then serves all the samples.\n" + "\n".join(table_lines or [])
            )
        return (
            "[choose_error_model] STOP: amplicon.dada2.learn_errors.error_model is binned, but "
            f"the filtered {direction} reads carry continuous quality values ({count} distinct "
            f"values above Q2, {declared_above[0]}-{declared_above[-1]}), as a MiSeq writes "
            f"them. The binned fit would join the error rates of {count} values with straight "
            "lines instead of fitting a curve. Set error_model to auto or loess, or delete the "
            "line from your run config."
        )
    # A function added to the registry later, set on the kind of values it is not for.
    return (
        f"[choose_error_model] STOP: amplicon.dada2.learn_errors.error_model is {function}, "
        f"which is made for {SUITED_CLASS[function]} quality values, but the filtered "
        f"{direction} reads carry {run_class} ones ({values_text(declared_above)}). Set "
        "error_model to auto, or delete the line from your run config."
    )


def warning_w1(direction, function, run_class, declared_above):
    """W1: a function set on the other kind of quality values, where that warns."""
    if function == "loess" and run_class == "binned":
        return (
            f"[choose_error_model] WARNING: {direction}: error_model is loess, but the filtered "
            f"reads carry binned quality values ({values_text(declared_above)}). DADA2's author "
            "recommends the binned fit (makeBinnedQualErrfun) for such data; loess is the fit "
            "MetaFlux used up to v2.6.1. Going on with loess, as set."
        )
    return (
        f"[choose_error_model] WARNING: {direction}: error_model is {function}, which is made "
        f"for {SUITED_CLASS[function]} quality values, but the filtered reads carry {run_class} "
        f"ones ({values_text(declared_above)}). Going on with {function}, as set."
    )


def warning_w2(direction, requested, declared_above, trunc_q):
    """W2: an undetermined kind of quality values with a function set in the config
    (custom included); the run goes on with it. One form per form of C4: 11-19 values,
    or 3-10 values that include neighbours, which names trunc_q when it is above 2."""
    neighbours = first_neighbours(declared_above)
    if len(declared_above) <= MAX_BINNED_VALUES and neighbours is not None:
        lower, upper = neighbours
        trunc_q_sentence = ""
        if trunc_q > NO_CALL_MAX:
            trunc_q_sentence = (f" amplicon.dada2.trunc_q is {trunc_q}, so only the values "
                                f"above Q{trunc_q} remain.")
        return (
            f"[choose_error_model] WARNING: {direction}: {len(declared_above)} distinct quality "
            f"values above Q2 ({values_text(declared_above)}), and some are neighbours "
            f"({lower} and {upper}), which binned quality values never are: neither binned nor "
            f"continuous as MetaFlux knows them.{trunc_q_sentence} Going on with {requested}, "
            "as set."
        )
    return (
        f"[choose_error_model] WARNING: {direction}: {len(declared_above)} distinct quality "
        "values above Q2, neither binned nor continuous as MetaFlux knows them. Going on with "
        f"{requested}, as set."
    )


# ── The decision ──────────────────────────────────────────────────────────────────

def voting_samples(direction, samples):
    """Which samples vote in one read direction.

    Parameters
    ----------
    direction : "R1" or "R2".
    samples : from read_quality_values.

    Returns
    -------
    (voters, everyone_votes): the set of names of the samples with at least
    MIN_READS_TO_VOTE reads scanned, or of every sample when none has that many (then
    everyone_votes is True).
    """
    voters = {sample["sample"] for sample in samples
              if sample[direction]["reads_scanned"] >= MIN_READS_TO_VOTE}
    if voters:
        return voters, False
    return {sample["sample"] for sample in samples}, True


def most_common_set(set_counts):
    """The set of values carried by most voters; a tie goes to the larger set, then to the
    set whose sorted list comes first (9, 20, 34 before 11, 25, 37).

    Parameters
    ----------
    set_counts : Counter of tuples of values above Q2 -> number of voters.

    Returns
    -------
    tuple of int.
    """
    return min(set_counts, key=lambda values: (-set_counts[values], -len(values), values))


def decide_direction(direction, samples, requested, trunc_q, out, log_path):
    """Decide the error function of one read direction (steps 1-5 of the module docstring).

    Parameters
    ----------
    direction : "R1" or "R2".
    samples : from read_quality_values.
    requested : the setting error_model: "auto", a registry name, or "custom" (the
        user's own function, experimental).
    trunc_q : amplicon.dada2.trunc_q, for the stop texts.
    out : the run's output folder, for the paths in the stop texts.
    log_path : this step's log, named where a stop text's table of samples is cut short.

    Returns
    -------
    (decision, warnings, stops). decision is the dict that goes into
    error_model_choice.json under this direction, or None when a stop applies; warnings
    and stops are lists of message texts.
    """
    warnings = []
    stops = []
    # A function set in the config, the user's own (custom) included: the stops that
    # exist only because auto has to choose become warnings below.
    forced = requested != "auto"
    custom = requested == CUSTOM_FUNCTION

    # Step 1: who votes.
    voter_names, everyone_votes = voting_samples(direction, samples)
    if everyone_votes:
        warnings.append(f"[choose_error_model] WARNING: no sample has 1,000 or more filtered "
                        f"{direction} reads; every sample votes.")
    voters = [(sample["sample"], sample[direction]) for sample in samples
              if sample["sample"] in voter_names]
    non_voters = [(sample["sample"], sample[direction]) for sample in samples
                  if sample["sample"] not in voter_names]
    table_lines = sample_table(direction, samples, voter_names, limit=MAX_ROWS_IN_STOP,
                               log_path=log_path)

    # Step 2: the values of the run (the voters' union) and its class.
    declared = set()
    for name, part in voters:
        declared.update(part["values"])
    declared = sorted(declared)
    declared_above = above_q2(declared)
    run_class = classify(declared_above)
    voter_class = {name: classify(above_q2(part["values"])) for name, part in voters}
    # Two binned tables can hold neighbouring values together although neither does on
    # its own: a NovaSeq X writes 40 on some versions of its software and 41 on others,
    # so samples of two such runs carry 40 and 41 between them. Such neighbours are no
    # sign of continuous values cut by trunc_q, which every sample would show on its own,
    # so the run stays binned, and C6 below names the two tables (with a function set in
    # the config, its warning).
    undetermined_by_neighbours = (run_class == "undetermined"
                                  and len(declared_above) <= MAX_BINNED_VALUES)
    every_voter_binned = all(kind == "binned" for kind in voter_class.values())
    if undetermined_by_neighbours and every_voter_binned:
        run_class = "binned"

    # Step 3, the stops. One voter without quality information among voters that have it.
    # This stop and the next (no quality information at all) describe the data, not the
    # choice, so they apply under every setting, custom included: an error model is a
    # rate per quality score, and no function can learn one from values that were made
    # up on the way.
    none_voters = [(name, part) for name, part in voters if voter_class[name] == "none"]
    if none_voters and run_class != "none":
        others = set()
        for name, part in voters:
            if voter_class[name] != "none":
                others.update(above_q2(part["values"]))
        stops.append(stop_c3(direction, none_voters, sorted(others), out))
        return None, warnings, stops

    # No quality information at all.
    if run_class == "none":
        stops.append(stop_c2(direction, declared, trunc_q, out, table_lines))
        return None, warnings, stops

    # A function set in the config on the kind of values where that stops: binned on
    # continuous values (C7; ON_OTHER_CLASS says "stop"). Checked here, before the
    # warnings below, because those say "going on with one ... model ..., as set": a run
    # that is about to stop must not print them first. The case that needs this order is
    # binned set on binned and continuous voters together: their values make the run
    # continuous, so it stops here, and C7 then describes the mix rather than "as a MiSeq
    # writes them" (only_together: not every voter carries continuous values on its own).
    # The user's own function (custom) suits no particular kind and never stops here.
    if forced and not custom:
        on_other_kind = (run_class in ("binned", "continuous")
                         and SUITED_CLASS[requested] != "both"
                         and run_class != SUITED_CLASS[requested])
        if on_other_kind and ON_OTHER_CLASS[requested] == "stop":
            only_together = any(voter_class[name] != run_class for name, part in voters)
            stops.append(stop_c7(direction, requested, run_class, declared_above,
                                 only_together, table_lines))
            return None, warnings, stops

    # Voters of different classes: two kinds of instrument in one run. Under auto this
    # stops, because auto cannot pick one function for both kinds; a function set in
    # the config (custom included) goes on, with a warning, since one model then serves
    # reads of two instruments.
    classes_present = sorted(set(voter_class.values()))
    if len(classes_present) > 1:
        if not forced:
            stops.append(stop_c5(direction, voters, voter_class, out, table_lines))
            return None, warnings, stops
        warnings.append(f"[choose_error_model] WARNING: {direction}: the voting samples carry "
                        "different kinds of quality values; going on with one "
                        f"{requested} model for all of them, as set.")
    elif classes_present == ["binned"]:
        # Binned voters with two different sets of values: two quality tables. A voter
        # that carries fewer values than the most common set is fine (a value can be
        # rare); a voter with a value outside it is not.
        set_counts = Counter(tuple(above_q2(part["values"])) for name, part in voters)
        common = most_common_set(set_counts)
        extra_voters = []
        for name, part in voters:
            extra = sorted(set(above_q2(part["values"])) - set(common))
            if extra:
                extra_voters.append((name, extra))
        if extra_voters:
            if not forced:
                stops.append(stop_c6(direction, common, extra_voters, out))
                return None, warnings, stops
            # The sets of the voters with extra values, next to the common one (a voter
            # whose set is a subset of the common one is not a second table).
            other_sets = set()
            for name, part in voters:
                if set(above_q2(part["values"])) - set(common):
                    other_sets.add(tuple(above_q2(part["values"])))
            sets_text = "; ".join([values_text(common)]
                                  + [values_text(values) for values in sorted(other_sets)])
            if requested == "binned":
                warnings.append(
                    f"[choose_error_model] WARNING: {direction}: the voting samples carry two "
                    f"sets of binned values ({sets_text}); going on with one binned model on all "
                    f"of them ({values_text(declared_above)}), as set. learnErrors learns from a "
                    "few whole samples, and the values of a group it does not draw are fitted "
                    "on averaged qualities of the other group.")
            else:
                warnings.append(
                    f"[choose_error_model] WARNING: {direction}: the voting samples carry two "
                    f"sets of binned values ({sets_text}); going on with one {requested} model "
                    "for all of them, as set.")

    # An undetermined kind of values: 11-19 values, or 3-10 that include neighbours. The
    # stop only says that auto cannot tell loess from binned; with a function set in the
    # config (custom included) the run goes on, with a warning that the values are
    # unusual.
    if run_class == "undetermined":
        if not forced:
            stops.append(stop_c4(direction, declared_above, trunc_q, out, table_lines))
            return None, warnings, stops
        warnings.append(warning_w2(direction, requested, declared_above, trunc_q))

    # Step 4: the function.
    if not forced:
        function = AUTO_FUNCTION[run_class]
    elif custom:
        # The user's own function: made for no particular kind of quality values, so
        # there is nothing to compare the class with (no C7 stop, no W1 warning).
        function = CUSTOM_FUNCTION
    else:
        function = requested
        # On the other kind of values: a function that stops there has stopped above,
        # before the warnings (C7); one that warns (loess on binned values) warns here.
        if (SUITED_CLASS[function] != "both"
                and run_class in ("binned", "continuous")
                and run_class != SUITED_CLASS[function]):
            if ON_OTHER_CLASS[function] == "warn":
                warnings.append(warning_w1(direction, function, run_class, declared_above))

    # Step 5: non-voters with values the voters do not carry. The note says that their
    # values are used neither for the choice nor as anchors of the binned fit; under
    # custom nothing is chosen and nothing is anchored (the user's function gets DADA2's
    # counts of whatever samples learnErrors draws), so it is left out.
    if custom:
        non_voters_to_note = []
    else:
        non_voters_to_note = non_voters
    for name, part in non_voters_to_note:
        extra = sorted(set(above_q2(part["values"])) - set(declared_above))
        if extra:
            warnings.append(
                f"[choose_error_model] NOTE: {name} ({part['reads_scanned']:,} {direction} reads "
                "scanned) has fewer than 1,000 filtered reads and does not vote; it carries "
                f"quality values the voting samples do not ({values_text(extra)}). They are not "
                "used to choose the error function, nor as anchors of the binned fit. If "
                "learnErrors draws this sample under the binned fit, a value outside the "
                "anchors stops the learning with the reason.")

    # The bases at each value of the run, summed over the voters: how much each value
    # weighs in the reads the model is learned from.
    bases_per_value = {}
    for value in declared:
        total = 0
        for name, part in voters:
            total += part["bases_per_value"].get(value, 0)
        bases_per_value[str(value)] = total

    decision = {
        "quality_class": run_class,
        "error_function": function,
        "decided_by": "setting" if forced else "auto",
        "quality_values": declared,
        "values_above_q2": declared_above,
        "n_values_above_q2": len(declared_above),
        "bases_per_value": bases_per_value,
        "samples_voting": len(voters),
        "samples_not_voting": [name for name, part in non_voters],
    }
    return decision, warnings, stops


def decision_line(direction, decision):
    """The log line that says what was decided, and why, such as
    "[choose_error_model] R2: binned quality values (9, 20, 34 above Q2) in 78 voting
    sample(s) -> binned fit, chosen by auto". The values handed to learning are added
    when they include values at or below Q2 (they then differ from those shown). The
    user's own function is named "your own error function (custom, experimental)"."""
    how = "chosen by auto" if decision["decided_by"] == "auto" else "as set in error_model"
    if decision["error_function"] == CUSTOM_FUNCTION:
        function_text = "your own error function (custom, experimental)"
    else:
        function_text = f"{decision['error_function']} fit"
    line = (f"[choose_error_model] {direction}: {decision['quality_class']} quality values "
            f"({values_text(decision['values_above_q2'])} above Q2) in "
            f"{decision['samples_voting']} voting sample(s) -> {function_text}, {how}")
    if decision["quality_values"] != decision["values_above_q2"]:
        line += f"; values handed to learning: {values_text(decision['quality_values'])}"
    return line


# ── Writing the outputs ───────────────────────────────────────────────────────────

def write_choice(path, requested, trunc_q, decisions, warnings, samples, voters):
    """Write stats/dada2/error_model_choice.json.

    Parameters
    ----------
    path : the output file.
    requested : the setting error_model.
    trunc_q : amplicon.dada2.trunc_q.
    decisions : {direction: decision} from decide_direction.
    warnings : every warning and note, in the order they were given.
    samples : from read_quality_values.
    voters : {direction: set of voting sample names}.

    Returns
    -------
    None. The file is indented JSON with a final line end; keys keep this order and
    samples are sorted, so the same inputs give the same bytes.
    """
    sample_entries = []
    for sample in samples:
        entry = {"sample": sample["sample"]}
        for direction in DIRECTIONS:
            part = sample[direction]
            entry[direction] = {"reads_scanned": part["reads_scanned"],
                                "votes": sample["sample"] in voters[direction],
                                "quality_values": part["values"]}
        sample_entries.append(entry)
    choice = {
        "requested": requested,
        "trunc_q": trunc_q,
        "thresholds": {"min_values": MIN_QUALITY_VALUES,
                       "max_binned_values": MAX_BINNED_VALUES,
                       "min_continuous_values": MIN_CONTINUOUS_VALUES,
                       "min_reads_to_vote": MIN_READS_TO_VOTE},
        "R1": decisions["R1"],
        "R2": decisions["R2"],
        "warnings": warnings,
        "samples": sample_entries,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(choice, indent=2) + "\n")


def write_quality_tsv(path, samples, voters):
    """Write stats/dada2/quality_values.tsv: one row per sample and read direction.

    Columns: sample, read_direction, reads_scanned, votes (yes/no), n_values_above_q2,
    values (comma-separated, Q2 and below included), quality_class (of this sample on
    its own), bases_at_or_below_q2. For people to read; nothing in the workflow reads it.
    """
    header = ["sample", "read_direction", "reads_scanned", "votes", "n_values_above_q2",
              "values", "quality_class", "bases_at_or_below_q2"]
    lines = ["\t".join(header)]
    for sample in samples:
        for direction in DIRECTIONS:
            part = sample[direction]
            values_above = above_q2(part["values"])
            low_bases = 0
            for value, count in part["bases_per_value"].items():
                if value <= NO_CALL_MAX:
                    low_bases += count
            row = [sample["sample"], direction, str(part["reads_scanned"]),
                   "yes" if sample["sample"] in voters[direction] else "no",
                   str(len(values_above)), ",".join(str(value) for value in part["values"]),
                   classify(values_above), str(low_bases)]
            lines.append("\t".join(row))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


# ── Main ──────────────────────────────────────────────────────────────────────────

def main():
    sm = snakemake  # noqa: F821 (injected by Snakemake)

    paths = sm.input.quality_values
    if isinstance(paths, str):
        paths = [paths]
    choice_path = Path(sm.output.choice)
    tsv_path = Path(sm.output.quality_tsv)
    # The output folder: the choice file is <out>/stats/dada2/error_model_choice.json.
    out = choice_path.parents[2]
    requested = sm.params.error_model_requested
    trunc_q = int(sm.params.trunc_q)
    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_lines = []

    def log(text):
        """A note, a warning or a decision: into the log and onto the terminal."""
        log_lines.append(text)
        print(text, file=sys.stderr, flush=True)

    def stop(text):
        """A stop: into the log, then once onto the terminal as the SystemExit text."""
        log_lines.append(text)
        log_path.write_text("\n".join(log_lines) + "\n")
        sys.exit(text)

    samples = read_quality_values(paths)

    log_lines.append(f"[choose_error_model] error_model: {requested}; {len(samples)} sample(s); "
                     f"a sample votes with at least {MIN_READS_TO_VOTE} reads scanned")
    voters = {}
    for direction in DIRECTIONS:
        voters[direction], everyone_votes = voting_samples(direction, samples)
        log_lines.extend(sample_table(direction, samples, voters[direction]))

    # Both directions are decided before anything stops, so all stop texts come at once.
    decisions = {}
    all_warnings = []
    all_stops = []
    for direction in DIRECTIONS:
        decision, warnings, stops = decide_direction(direction, samples, requested, trunc_q, out,
                                                     log_path)
        for warning in warnings:
            log(warning)
        all_warnings.extend(warnings)
        all_stops.extend(stops)
        decisions[direction] = decision
    if all_stops:
        stop("\n\n".join(all_stops))

    classes = (decisions["R1"]["quality_class"], decisions["R2"]["quality_class"])
    if classes[0] != classes[1] and set(classes) == {"binned", "continuous"}:
        warning = (f"[choose_error_model] WARNING: R1 carries {classes[0]} and R2 {classes[1]} "
                   "quality values, which is unusual for one sequencing run; check the input "
                   "files.")
        log(warning)
        all_warnings.append(warning)

    for direction in DIRECTIONS:
        log(decision_line(direction, decisions[direction]))

    write_choice(choice_path, requested, trunc_q, decisions, all_warnings, samples, voters)
    write_quality_tsv(tsv_path, samples, voters)
    log_lines.append("[choose_error_model] DONE. error_model_choice.json and quality_values.tsv written.")
    log_path.write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
