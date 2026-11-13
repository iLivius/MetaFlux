#!/usr/bin/env python3
"""Stop before denoising when a DADA2 error model did not converge and was not accepted.

DADA2 learns one error model per read direction for the whole run (rule learn_errors,
60d_learn_errors.R): for each quality score, how often each base is read as each
other base. It repeats two steps until a round gives back a table it has seen before,
which MetaFlux calls "converged". A model can also end without converging: at the cap
on rounds (amplicon.dada2.learn_errors.max_consist), or by MetaFlux's flat-rounds
rule (max_flat_rounds), which ends a loop that has stopped improving. Such a model is
often fine, but the numbers alone cannot tell it from one fitted to faulty reads: a
library with reads of both strands in one file looked exactly like a slow, healthy
one. So the run does not go on with it silently. This step stops the run before
denoising, says what to check, and prints the line that accepts this one model by its
code, such as R2:0a55be94 (read direction, a colon, the first 8 characters of an md5
sum of the model's numbers).

Inputs (see rule check_error_model in 60_dada2.smk):
  record : stats/dada2/error_model.json from learn_errors: per read direction the
           rounds, how the learning ended (stopped_by: dada2_repeat, flat_rounds or
           cap), the change in every round and the model code.
  report : stats/dada2/error_model_report.json from error_model_report (60e), taken so
           that the report is on disk before this step can stop the run (the stop text
           sends the user to its figures). Only its list of warnings is read: when the
           report had trouble (a figure it could not draw, a refit that does not give
           back the model), a note says so on the terminal.
  config : amplicon.dada2.learn_errors.accept_unconverged, the accepted codes, read
           from the config inside the job rather than from params, so that adding or
           removing a code makes no step run again (00_common.smk checks and
           normalises the list at start-up; accepted_codes() below normalises it
           again, in case it ever arrives as written in the run config).
  params : orientation, only to say whether stats/orientation/ exists.

Decision, per read direction:
  converged       stopped_by is dada2_repeat;
  accepted        not converged, and the model code is in accept_unconverged;
  not accepted    otherwise.
With a direction not accepted, the stop text (stop_text below) goes to the log and
then, once, to the terminal through sys.exit (Python prints the text of SystemExit
itself, so it is not printed before as well), and no verdict is written, so
dada_seqtab cannot run. Otherwise the verdict is written: one line, the same for the
same models, such as
    R1:0a55be94 converged; R2:7ebe3bed accepted (not converged: flat-rounds rule, round 42)
with a warning on the terminal for each accepted model, a note for listed codes
that name no model of this run, and a note when the report has warnings.

Output: stats/dada2/error_model_check.txt, the verdict, which dada_seqtab takes as an
input and copies into its log. In a finished run it records which models the ASV
tables were made with, and whether one of them was accepted without converging.
"""
import json
import sys
from pathlib import Path

# MetaFlux's default cap on rounds (config/config.yaml); a cap below it gets its own
# first paragraph in the stop text, because raising it is then the likely fix.
DEFAULT_MAX_CONSIST = 200


def rounds_text(count):
    """Return "1 round" or "200 rounds", for the messages."""
    return "1 round" if count == 1 else f"{count} rounds"


def how_it_ended(part, max_consist):
    """Say in a few words how a model that did not converge ended.

    Parameters
    ----------
    part : the record's entry of one read direction.
    max_consist : the cap on rounds.

    Returns
    -------
    str, such as "flat-rounds rule, round 42" or "cap of 200 rounds"; used in the
    verdict and in the warning for an accepted model.
    """
    if part["stopped_by"] == "flat_rounds":
        return f"flat-rounds rule, round {part['flat_stop_round']}"
    return f"cap of {rounds_text(max_consist)}"


def number(value):
    """A change value with three significant digits, as the learning log writes it."""
    return f"{value:.3g}"


def accepted_codes(learn_cfg):
    """The accepted model codes, as a set, from the learn_errors block of the config.

    00_common.smk checks amplicon.dada2.learn_errors.accept_unconverged at start-up and
    writes it back into the config as a sorted list of normalised codes, so this step
    normally receives exactly that list. The normalisation is repeated here so that the
    check still accepts the right models if that write-back were ever lost: a single
    code written without brackets would otherwise arrive as text, which set() would
    split into its characters, and a code typed in other capitals ("r2:0A55BE94") would
    match no model.

    Parameters
    ----------
    learn_cfg : dict, the config's amplicon -> dada2 -> learn_errors block (or {}).

    Returns
    -------
    set of str, each code with "R1"/"R2" in capitals and the md5 characters in lower
    case, as the stop message prints them. Items that are not text are left out (the
    start-up check stops the run on them before any job runs).
    """
    listed = learn_cfg.get("accept_unconverged") or []
    if isinstance(listed, str):
        listed = [listed]                    # one code written without the brackets
    codes = set()
    for item in listed:
        if not isinstance(item, str):
            continue
        text = item.strip()
        codes.add(text[:2].upper() + text[2:].lower())
    return codes


def decide(record, accepted):
    """Sort the two read directions into converged, accepted and not accepted.

    Parameters
    ----------
    record : the learning record (dict).
    accepted : the set of accepted codes.

    Returns
    -------
    dict {direction: "converged" | "accepted" | "unaccepted"}, for R1 and R2.
    """
    status = {}
    for direction in ("R1", "R2"):
        part = record[direction]
        if part["converged"]:
            status[direction] = "converged"
        elif part["model_code"] in accepted:
            status[direction] = "accepted"
        else:
            status[direction] = "unaccepted"
    return status


def verdict_line(record, status):
    """The one-line verdict written for dada_seqtab.

    Parameters
    ----------
    record : the learning record.
    status : from decide(); no direction may be "unaccepted".

    Returns
    -------
    str, without the line end.
    """
    parts = []
    for direction in ("R1", "R2"):
        code = record[direction]["model_code"]
        if status[direction] == "converged":
            parts.append(f"{code} converged")
        else:
            how = how_it_ended(record[direction], record["max_consist"])
            parts.append(f"{code} accepted (not converged: {how})")
    return "; ".join(parts)


def direction_paragraph(direction, part, record):
    """The paragraph of the stop text that says how one direction's learning ended.

    Parameters
    ----------
    direction : "R1" or "R2".
    part : the record's entry of that direction (not converged).
    record : the whole record (for the cap and N).

    Returns
    -------
    str: one paragraph, wrapped by hand at about 95 characters.
    """
    changes = part["change_per_round"]
    if part["stopped_by"] == "flat_rounds":
        stop_round = part["flat_stop_round"]
        n_flat = record["max_flat_rounds"]
        second = f" (round 2: {number(changes[1])})" if len(changes) > 1 else ""
        # The window of the rule: "the last round" for N = 1, "the last 20 rounds" otherwise.
        window = "the last round" if n_flat == 1 else f"the last {n_flat} rounds"
        return (
            f"{direction}: MetaFlux ended the learning in round {stop_round} because the model "
            f"had stopped improving\n(amplicon.dada2.learn_errors.max_flat_rounds: {n_flat}): "
            f"it changed by {number(changes[stop_round - 1])} in that round, less\nthan 1/1000 "
            f"of round 1's {number(changes[0])}{second}, and the lowest change so far fell by "
            f"less\nthan 10 % over {window}. The model kept is the table that "
            f"went into round {stop_round}."
        )
    cap = record["max_consist"]
    allowed = f"all {cap} rounds it is allowed" if cap > 1 else "the one round it is allowed"
    if len(changes) >= 3:
        movement = (f"by {number(changes[-1])} in its last round, against {number(changes[0])} "
                    f"and {number(changes[1])} in\nrounds 1 and 2.")
    elif len(changes) == 2:
        movement = (f"by {number(changes[1])} in its last round, round 2, against "
                    f"{number(changes[0])} in round 1.")
    else:
        movement = f"by {number(changes[0])} in its only round."
    return (
        f"{direction}: learnErrors used {allowed} (amplicon.dada2.learn_errors.max_consist)\n"
        f"and the model was still changing: {movement}"
    )


def stop_text(record, status, accepted, out, log_path, orientation):
    """The full stop text: what happened, what to check, and how to go on.

    Parameters
    ----------
    record : the learning record.
    status : from decide(), with at least one "unaccepted" direction.
    accepted : the accepted codes of this run (a set), echoed at the end.
    out : the run's output folder (absolute), so every path printed is absolute: the
        user's shell reads a relative path against the MetaFlux folder, where
        MetaFlux is started, not against the output folder.
    log_path : this step's log, which holds the same text.
    orientation : "fixed" or "mixed".

    Returns
    -------
    str: the paragraphs, separated by blank lines; the first and the last line each
    say where the text is kept, so that either one, read alone, leads to the log.
    """
    unaccepted = [d for d in ("R1", "R2") if status[d] == "unaccepted"]
    unaccepted_codes = [record[d]["model_code"] for d in unaccepted]
    # The codes this run needs to go on: every model that did not converge, accepted
    # already or not, so that the printed line can replace the one in the run config.
    needed_codes = sorted(record[d]["model_code"] for d in ("R1", "R2")
                          if not record[d]["converged"])
    stopped_by = {record[d]["stopped_by"] for d in unaccepted}
    cap = record["max_consist"]

    model_words = "error model" if len(unaccepted) == 1 else "error models"
    this_model = "this model" if len(unaccepted) == 1 else "these models"
    # The first line names the model codes, so that it identifies the models on its own
    # (in a scrolled terminal, or quoted in a question).
    paragraphs = [
        f"[check_error_model] STOP: the DADA2 {model_words} of the {' and '.join(unaccepted)} "
        f"reads did not converge ({', '.join(unaccepted_codes)}),\nso MetaFlux stops before "
        f"denoising. What to check, and how to go on with {this_model}, follows;\nthe same "
        f"text is in {log_path}."
    ]

    if "cap" in stopped_by and cap < DEFAULT_MAX_CONSIST:
        paragraphs.append(
            f"Your limit of {rounds_text(cap)} (amplicon.dada2.learn_errors.max_consist) is below "
            f"MetaFlux's default\nof {DEFAULT_MAX_CONSIST}, and some libraries need 20 to 100 "
            "rounds. Raise it, or delete the line from your run\nconfig, and start MetaFlux "
            "again; a model that converges within the new limit needs no code."
        )

    for direction in unaccepted:
        paragraphs.append(direction_paragraph(direction, record[direction], record))

    # Each path ends its line and the explanation follows, indented, so that the text
    # reads the same however long the path of the output folder is.
    first_direction = unaccepted[0]
    read_checks = [
        "These numbers cannot tell a usable model from one fitted to faulty reads: a library "
        "with reads\nof both strands in one file looked just like this. Check the reads first:",
        f"  - {out}/logs/pick_trunclen.log\n"
        "    a \"very short primary cuts\" warning can mean that the quality rises along part "
        "of the\n    reads, as it does when some of them are reversed;",
        f"  - {out}/stats/falco/<sample>_{first_direction}_stripped/fastqc_report.html\n"
        "    the per-base quality should fall along the read; a rise from low values at the "
        "start is\n    the sign of reversed reads, and flat lines at 2 to 4 values mean "
        "binned quality scores;",
    ]
    if orientation == "mixed":
        read_checks.append(f"  - {out}/stats/orientation/\n"
                           "    how many pairs each pass kept;")
    read_checks.append(
        f"  - {out}/03.stripped/\n"
        "    the start of the reads: \"DADA2's error model never converges\" in "
        "docs/troubleshooting.md\n    gives a one-line check for primer bases left on the "
        "reads;"
    )
    read_checks.append(
        f"  - {out}/stats/dada2/error_model_quality_values.png\n"
        "    few quality values (binned) need error_model: binned or auto, not loess."
    )
    read_checks.append(
        "If the reads are clean, slow learning is common in very diverse samples (soil, "
        "sediment), and\nDADA2's author advises using such a model."
    )
    if "flat_rounds" in stopped_by:
        read_checks.append(
            "On the libraries replayed to set this rule, a model stopped by it changed at "
            "most 0.04 % of the\nreads compared with the model the loop ended with (loess "
            "fit), and at most 0.02 % (binned fit)."
        )
    paragraphs.append("\n".join(read_checks))

    model_checks = ["Then look at the model:"]
    for direction in unaccepted:
        model_checks.append(
            f"  - {out}/stats/dada2/error_model_fit_{direction}.png\n"
            "    the solid line should follow the filled points and fall as quality rises;"
        )
    model_checks.append(
        f"  - {out}/stats/dada2/error_model_convergence.png\n"
        "    the change per round, the stop and the cap;"
    )
    model_checks.append(
        f"  - {out}/stats/dada2/error_model_fit.tsv and {out}/stats/dada2/error_model.json\n"
        "    the numbers."
    )
    # The documentation page on the error model explains every figure and table, with
    # the values of the models MetaFlux was tested on; the troubleshooting entry walks
    # through the checks of this text.
    model_checks.append(
        "How to read them: docs/amplicon/error-model.md, and \"DADA2's error model never "
        "converges\" in\ndocs/troubleshooting.md."
    )
    paragraphs.append("\n".join(model_checks))

    codes_list = ", ".join(f'"{code}"' for code in needed_codes)
    paragraphs.append(
        f"To go on with {this_model}, add this line to the run config of this output folder "
        "(the file you\ngive with --configfile, not config/config.yaml):\n"
        f"    accept_unconverged: [{codes_list}]\n"
        "It belongs under amplicon: -> dada2: -> learn_errors:. Put it inside the amplicon: "
        "block your\nrun config already has (replace an accept_unconverged line that is "
        "already there); do not paste a\nsecond amplicon: block, because YAML keeps only "
        "the last of two blocks with the same name, and your\nother amplicon settings would "
        "be lost. In a run config that has none of these blocks yet:\n"
        "    amplicon:\n"
        "      dada2:\n"
        "        learn_errors:\n"
        f"          accept_unconverged: [{codes_list}]\n"
        "Then start MetaFlux again with the same command. The model is not learned again: "
        "dada_seqtab and\nthe steps after it run. A code names one model: keep the line, and "
        "a later rerun of this folder\ngoes on; other reads, settings or software, and "
        "possibly another type of computer, give another\nmodel with another code, and the "
        "run stops here again.\n"
        f"(accept_unconverged in this run: {json.dumps(sorted(accepted))}.)"
    )
    paragraphs.append(f"[check_error_model] This text is also in {log_path}.")
    return "\n\n".join(paragraphs)


def report_warnings_note(report_path, out):
    """The note that repeats the warnings of the error-model report, or None.

    The report (rule error_model_report, 60e) never fails because of a figure or a
    table: it draws a placeholder instead and lists the problem in the "warnings" of
    error_model_report.json. This note makes sure such a problem is seen before the
    ASV tables are made with the model.

    Parameters
    ----------
    report_path : stats/dada2/error_model_report.json.
    out : the run's output folder (absolute), for the paths in the note.

    Returns
    -------
    str or None: the note, naming how many warnings there are and the first one; a
    note that the file could not be read; or None when the report has no warnings.
    """
    try:
        warnings = json.loads(Path(report_path).read_text()).get("warnings") or []
    except (OSError, ValueError, AttributeError) as error:
        return (
            f"[check_error_model] NOTE: the warnings of the error-model report could not be "
            f"read from {report_path} ({error}); its log is {out}/logs/error_model_report.log."
        )
    if not warnings:
        return None
    first = str(warnings[0]).rstrip(".")
    return (
        f"[check_error_model] NOTE: the error-model report has {len(warnings)} warning(s), "
        f"the first: {first}. All of them are in {out}/stats/dada2/error_model_report.json "
        f"and {out}/logs/error_model_report.log."
    )


def main():
    sm = snakemake  # noqa: F821 (injected by Snakemake)

    record_path = Path(sm.input.record)
    record = json.loads(record_path.read_text())
    # The output folder: the record is <out>/stats/dada2/error_model.json.
    out = record_path.parents[2]
    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    verdict_path = Path(sm.output.verdict)
    orientation = sm.params.orientation

    # Read from the config, not from params (see the module docstring); every level
    # may be missing in a stand-in config, which then accepts nothing.
    learn_cfg = sm.config.get("amplicon", {}).get("dada2", {}).get("learn_errors", {})
    accepted = accepted_codes(learn_cfg)

    status = decide(record, accepted)

    # A model that did not converge and is not accepted: the stop text goes into the
    # log, then onto the terminal once, as the text of SystemExit. No verdict is
    # written, so dada_seqtab cannot run.
    if "unaccepted" in status.values():
        text = stop_text(record, status, accepted, out, log_path, orientation)
        log_path.write_text(text + "\n")
        sys.exit(text)

    # Every model converged or was accepted: say on the terminal (and in the log) what
    # the user should know, then write the verdict.
    log_lines = []
    for direction in ("R1", "R2"):
        if status[direction] == "accepted":
            part = record[direction]
            log_lines.append(
                f"[check_error_model] WARNING: going on with {part['model_code']}, which did "
                f"not converge ({how_it_ended(part, record['max_consist'])}), because "
                "amplicon.dada2.learn_errors.accept_unconverged lists its code."
            )
    codes_of_run = {record[d]["model_code"] for d in ("R1", "R2")}
    stale = sorted(accepted - codes_of_run)
    if stale:
        log_lines.append(
            f"[check_error_model] NOTE: accept_unconverged also lists {', '.join(stale)}, "
            "which name no model of this run."
        )
    # The report's own warnings (a figure it could not draw, for example).
    report_note = report_warnings_note(sm.input.report, out)
    if report_note:
        log_lines.append(report_note)
    for message in log_lines:
        print(message, file=sys.stderr)

    line = verdict_line(record, status)
    verdict_path.parent.mkdir(parents=True, exist_ok=True)
    verdict_path.write_text(line + "\n")
    log_lines.append(f"[check_error_model] {line}")
    log_path.write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
