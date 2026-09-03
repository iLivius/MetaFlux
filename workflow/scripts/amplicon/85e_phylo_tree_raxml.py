#!/usr/bin/env python3
"""Infer the ASV tree with RAxML-NG (alternative backend).

Third step of the optional 16S phylogeny (see 85_phylogeny.smk), for users who choose
phylogeny.backend: raxml-ng. Full maximum likelihood, run in two steps.

WHY TWO STEPS — THIS IS NOT AN OPTIMISATION
-------------------------------------------
RAxML-NG requires a minimum number of alignment patterns per thread (its own
rule of thumb is 200-1000) and it ENFORCES that: exceed roughly twice its
recommendation and it prints "Too few patterns per thread! RAxML-NG will terminate now
to avoid wasting resources" and stops. A 16S amplicon alignment has only ~430 distinct
patterns, which works out at about 2 recommended threads and a hard failure somewhere
above 4. Handing it a normal Snakemake allocation of 8 or 16 threads does not run
slowly — it crashes the rule outright.

So the first step, --parse, is run precisely to ask RAxML-NG how many threads IT thinks
this alignment supports. It also validates the alignment, reports an expected memory
requirement (logged, useful for tuning resources.mem_mb), and writes a binary .rba copy
that the second step loads quickly. The second step, --search, then runs the actual
inference with min(what Snakemake allocated, what RAxML-NG recommended).

--workers is likewise an explicit integer rather than `auto`. Workers parallelise
across independent starting trees rather than across alignment sites, which is the
scaling axis that actually works for a short alignment — but `auto` reads the machine,
and this workflow's rule is that the same input must not behave differently on
different hardware.

--seed is mandatory here in a way it is not for other tools: RAxML-NG's default seed is
the current wall-clock time, so an unseeded run is non-reproducible by design. It comes
from amplicon.seed, the same seed as the rest of the amplicon path.

--extra seq-dup-remove is never exposed. It drops duplicate sequences from the tree,
which would leave the tree describing a different feature set than the abundance table
— the one outcome this module must never produce. The defaults keep every tip.

THE RERUN TRAP
--------------
Like IQ-TREE, RAxML-NG refuses to overwrite the outputs of a previous run at the same
prefix, and it resumes automatically from a leftover .raxml.ckp checkpoint. Both
prefixes (the --parse one and the --search one) are therefore cleared before starting.
Snakemake handles the DECLARED outputs; this clean covers everything else RAxML-NG
scatters around its prefix.

Inputs (see rule phylo_tree, raxml-ng branch, in 85_phylogeny.smk)
  aln : 7.phylogeny/asv_16s.aln.fasta, from phylo_align.

Outputs
  tree     : <backend dir>/asv_16s.raxml.bestTree — the best-scoring ML tree, read
             next by phylo_export.
  run_json : <backend dir>/tree_run.json — resolved model, the recommended/allocated/
             effective thread counts, memory estimate, version, seed and both command
             lines. Folded into phylogeny.params.json by phylo_export.
"""
import json
import re
import shlex
import subprocess
import sys
from pathlib import Path


def clear_stale_prefix_files(prefix: Path, log) -> None:
    """Remove leftovers from a previous run at this prefix.

    RAxML-NG will not overwrite existing output files, and it silently resumes from a
    leftover .raxml.ckp — so an interrupted earlier run would otherwise either abort
    this one or contaminate it with stale state.
    """
    removed = []
    for stale in sorted(prefix.parent.glob(prefix.name + ".*")):
        if stale.is_file():
            stale.unlink()
            removed.append(stale.name)
    if removed:
        log(f"[phylo_tree] Cleared {len(removed)} file(s) from a previous run: "
            f"{', '.join(removed)}")


def raxml_version(log) -> str:
    """Record which RAxML-NG build produced this tree."""
    try:
        completed = subprocess.run(
            ["raxml-ng", "--version"],
            capture_output=True, text=True, check=False,
        )
        for line in (completed.stdout or completed.stderr).splitlines():
            if "RAxML-NG" in line:
                return line.strip()
        return "unknown"
    except Exception as exc:                                    # noqa: BLE001
        log(f"[phylo_tree] WARNING: could not read RAxML-NG's version ({exc})")
        return "unknown"


def parse_thread_recommendation(parse_output: str, log) -> int | None:
    """Read the thread count RAxML-NG recommends out of its --parse output.

    The line looks like:  "* Recommended number of threads / MPI processes: 2"
    Returning None means the line was not found; the caller then falls back to a
    conservative fixed value rather than guessing high, because guessing high is the
    failure mode that crashes the search.
    """
    match = re.search(
        r"Recommended number of threads[^:]*:\s*(\d+)", parse_output, re.IGNORECASE
    )
    if match:
        return int(match.group(1))
    log("[phylo_tree] WARNING: could not find a thread recommendation in the --parse "
        "output.")
    return None


def parse_memory_estimate(parse_output: str) -> str | None:
    """Read RAxML-NG's own memory estimate from --parse, for the run record."""
    match = re.search(r"Estimated memory requirements[^:]*:\s*(.+)", parse_output)
    return match.group(1).strip() if match else None


def main() -> int:
    sm = snakemake  # noqa: F821

    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = log_path.open("w")

    def log(message: str) -> None:
        print(message, file=log_fh, flush=True)
        print(message, file=sys.stderr, flush=True)

    aln_in       = Path(sm.input.aln)
    tree_out     = Path(sm.output.tree)
    run_json_out = Path(sm.output.run_json)

    prefix          = Path(sm.params.prefix)
    parse_prefix    = Path(sm.params.parse_prefix)
    model           = str(sm.params.model)
    model_requested = str(sm.params.model_requested)
    seed            = int(sm.params.seed)
    extra_args      = str(sm.params.extra_args or "").strip()
    allocated_threads = int(sm.threads)

    prefix.parent.mkdir(parents=True, exist_ok=True)
    clear_stale_prefix_files(parse_prefix, log)
    clear_stale_prefix_files(prefix, log)

    # ── Step 1: --parse ───────────────────────────────────────────────────────
    parse_command = [
        "raxml-ng", "--parse",
        "--msa", str(aln_in),
        "--model", model,
        "--prefix", str(parse_prefix),
    ]
    log(f"[phylo_tree] Backend raxml-ng, model {model} (requested: {model_requested}), "
        f"seed {seed}")
    log(f"[phylo_tree] Step 1/2 — asking RAxML-NG how many threads this alignment supports:")
    log(f"[phylo_tree] Running: {' '.join(parse_command)}")

    parse_run = subprocess.run(parse_command, capture_output=True, text=True)
    parse_output = (parse_run.stdout or "") + (parse_run.stderr or "")

    log("[phylo_tree] --- raxml-ng --parse output ---")
    log(parse_output.rstrip())
    log("[phylo_tree] -------------------------------")

    if parse_run.returncode != 0:
        log(f"[phylo_tree] ERROR: raxml-ng --parse exited with status {parse_run.returncode}")
        log_fh.close()
        return parse_run.returncode

    recommended = parse_thread_recommendation(parse_output, log)
    memory_estimate = parse_memory_estimate(parse_output)
    if memory_estimate:
        log(f"[phylo_tree] RAxML-NG's own memory estimate: {memory_estimate}")

    # Take the smaller of what Snakemake allocated and what RAxML-NG recommends. The
    # recommendation is the binding constraint: exceeding roughly twice it is a hard
    # error, not a slowdown. When the recommendation could not be read, fall back to a
    # deliberately conservative 2 rather than trusting the allocation.
    if recommended is None:
        effective_threads = min(allocated_threads, 2)
        log(f"[phylo_tree] Falling back to {effective_threads} thread(s) — RAxML-NG "
            "hard-errors when given too many threads for a short alignment, so a low "
            "guess is the safe one.")
    else:
        effective_threads = max(1, min(allocated_threads, recommended))
        log(f"[phylo_tree] Threads: {effective_threads} "
            f"(Snakemake allocated {allocated_threads}, RAxML-NG recommends {recommended})")

    # Workers parallelise across independent starting trees rather than across alignment
    # sites, which is the scaling axis that actually helps a short alignment. Explicit
    # integer, never `auto`, so the run does not depend on the machine.
    #
    # The count MUST divide the thread count exactly: RAxML-NG refuses to start
    # otherwise ("The specified number of threads (5) is not a multiple of the number of
    # parallel tree searches (2)"). Halving is therefore not enough on its own — 5 // 2
    # is 2, which does not divide 5. Take the largest divisor that is at most half the
    # threads, so each worker gets at least two threads, and fall back to a single
    # worker for thread counts with no such divisor (1, 2, 3, and every prime).
    workers = 1
    for candidate in range(effective_threads // 2, 1, -1):
        if effective_threads % candidate == 0:
            workers = candidate
            break
    log(f"[phylo_tree] Workers: {workers} "
        f"({effective_threads} threads / {workers} worker(s) = "
        f"{effective_threads // workers} thread(s) each; RAxML-NG requires an exact divisor)")

    parse_binary = Path(str(parse_prefix) + ".raxml.rba")
    if not parse_binary.exists():
        log(f"[phylo_tree] ERROR: --parse did not produce {parse_binary.name}")
        log_fh.close()
        return 1

    # ── Step 2: --search ──────────────────────────────────────────────────────
    search_command = [
        "raxml-ng", "--search",
        "--msa", str(parse_binary),
        "--model", model,
        "--seed", str(seed),
        "--threads", str(effective_threads),
        "--workers", str(workers),
        "--prefix", str(prefix),
    ]
    if extra_args:
        search_command += shlex.split(extra_args)
        log(f"[phylo_tree] Appending extra_args.raxml_ng: {extra_args}")

    log(f"[phylo_tree] Step 2/2 — inference:")
    log(f"[phylo_tree] Running: {' '.join(search_command)}")

    search_run = subprocess.run(search_command, capture_output=True, text=True)
    search_output = (search_run.stdout or "") + (search_run.stderr or "")

    log("[phylo_tree] --- raxml-ng --search output ---")
    log(search_output.rstrip())
    log("[phylo_tree] --------------------------------")

    if search_run.returncode != 0:
        log(f"[phylo_tree] ERROR: raxml-ng --search exited with status "
            f"{search_run.returncode}")
        log_fh.close()
        return search_run.returncode

    if not tree_out.exists():
        log(f"[phylo_tree] ERROR: raxml-ng reported success but {tree_out.name} is missing.")
        log_fh.close()
        return 1

    # RAxML-NG writes .raxml.bestModel after EVERY search, not only when model selection
    # ran. Its contents are the fitted parameters of whatever model was used — e.g.
    # "GTR{0.97/1.08/...}+FU{0.24/...}+G4m{99.8}" for a plain GTR+G run — not the name of
    # a chosen model. Reading it unconditionally into model_used would therefore replace
    # a pinnable model name ("GTR+G") with a parameter dump, and claim a selection step
    # that never happened.
    #
    # So: the fitted parameters are recorded under their own key, always, because they
    # are genuine provenance. model_used only changes when MOOSE was actually asked for
    # (model: DNA), which is the one case where the model itself was not known in advance.
    best_model_file = Path(str(prefix) + ".raxml.bestModel")
    model_parameters = None
    if best_model_file.exists():
        lines = best_model_file.read_text(errors="replace").strip().splitlines()
        if lines:
            model_parameters = lines[0].split(",")[0].strip()

    moose_requested = model.strip().upper() == "DNA"
    if moose_requested and model_parameters:
        model_used = model_parameters
        log(f"[phylo_tree] MOOSE model selection chose: {model_used}")
    else:
        model_used = model
        if model_parameters:
            log(f"[phylo_tree] Fitted parameters for {model}: {model_parameters}")

    record = {
        "tool": "raxml-ng",
        "version": raxml_version(log),
        "backend": "raxml-ng",
        "model_requested": model_requested,
        "model_passed": model,
        "model_used": model_used,
        # The optimised rate/frequency parameters RAxML-NG fitted. Recorded separately
        # from model_used so the model stays a name a user can pin.
        "model_parameters": model_parameters,
        "model_selection_ran": moose_requested,
        "threads_allocated": allocated_threads,
        "threads_recommended_by_parse": recommended,
        "threads_effective": effective_threads,
        "workers": workers,
        "memory_estimate": memory_estimate,
        "seed": seed,
        "extra_args": extra_args,
        "command_parse": " ".join(parse_command),
        "command_search": " ".join(search_command),
    }
    run_json_out.parent.mkdir(parents=True, exist_ok=True)
    run_json_out.write_text(json.dumps(record, indent=2) + "\n")

    log(f"[phylo_tree] DONE. Tree written to {tree_out.name}")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
