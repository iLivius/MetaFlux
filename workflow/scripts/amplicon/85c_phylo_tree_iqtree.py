#!/usr/bin/env python3
"""Infer the ASV tree with IQ-TREE 3 (the default backend).

Third step of the optional 16S phylogeny (see 85_phylogeny.smk). Reads the MAFFT
alignment and runs full maximum-likelihood inference. IQ-TREE is the default of the
three backends because it is a genuine ML implementation with real model selection
available, at a runtime that is perfectly tolerable for the few hundred to few
thousand ASVs a typical amplicon run produces.

THE MODEL, AND WHY IT IS PINNED
-------------------------------
The shipped default is GTR+F+G4, not IQ-TREE's own default of -m MFP. MFP runs
ModelFinder on every execution, crossing 22 base substitution models with frequency
variants and several kinds of rate heterogeneity. That is a lot of work whose cost
depends on the data, so runtime becomes unpredictable. On the 16S test data the search
picked TPM3u+R5 on real ASVs and GTR+F+I+R5 on reference fragments — GTR-family
matrices whose choice barely moves the distances diversity metrics use. IQ-TREE's own
documentation recommends exactly this pattern: select once, then pin the result for
subsequent analyses.

Model selection is still available: setting phylogeny.model to MFP runs ModelFinder,
and the model it chooses is written into the .iqtree report and recorded in this
script's run JSON. That is the distinction this module draws between kinds of
automatic behaviour — one that writes down what it picked is fine, one that silently
depends on the machine is not.

THREADS
-------
Modest on purpose. IQ-TREE parallelises the likelihood calculation across alignment
COLUMNS, so its parallel efficiency improves with alignment LENGTH. A 16S amplicon
alignment is ~250-430 columns, which is short; IQ-TREE's own documentation warns that
too many cores can actually slow a short alignment down. -T AUTO is never used: it
picks a thread count by timing the machine it is on, under whatever load it is under
at that moment, so the same input would run differently on different hardware — or on
the same hardware on a busy afternoon.

THE RERUN TRAP
--------------
IQ-TREE refuses to re-run an analysis that already completed, and — more subtly — it
restores the random seed from its .ckp.gz checkpoint file. So a leftover checkpoint
from an earlier run would either abort this rule outright or silently override a
changed seed. This script therefore clears everything sharing the output prefix before
starting. Note the division of labour: Snakemake already removes a rule's DECLARED
outputs before re-running it, so the manual clean here is what handles the files
IQ-TREE writes that are NOT declared (.log, .bionj, .mldist, .uniqueseq.phy, and the
extra files produced under -B or MFP).

Inputs (see rule phylo_tree, iqtree branch, in 85_phylogeny.smk)
  aln : 7.phylogeny/asv_16s.aln.fasta, from phylo_align.

Outputs
  tree     : <backend dir>/asv_16s.treefile — the ML tree, read next by phylo_export,
             which validates it and writes the canonical unrooted Newick.
  report   : <backend dir>/asv_16s.iqtree — IQ-TREE's own report, including the model
             it used or selected. Inspection only; nothing reads it back.
  ckp      : <backend dir>/asv_16s.ckp.gz — IQ-TREE's checkpoint. Declared so
             Snakemake cleans it between runs.
  run_json : <backend dir>/tree_run.json — resolved model, threads, version, seed and
             command line, folded into phylogeny.params.json by phylo_export.
"""
import json
import shlex
import subprocess
import sys
from pathlib import Path


def clear_stale_prefix_files(prefix: Path, log) -> None:
    """Remove every file IQ-TREE may have left from a previous run of this prefix.

    Necessary because IQ-TREE treats a completed analysis as something not to repeat,
    and reads its seed back from .ckp.gz — meaning a stale checkpoint would quietly
    win over the seed this run was given. Removing the files is preferred over passing
    --redo because it also guarantees no half-written artifact from a crashed run is
    mistaken for output of this one.
    """
    removed = []
    for stale in sorted(prefix.parent.glob(prefix.name + ".*")):
        if stale.is_file():
            stale.unlink()
            removed.append(stale.name)
    if removed:
        log(f"[phylo_tree] Cleared {len(removed)} file(s) from a previous run: "
            f"{', '.join(removed)}")


def iqtree_version(log) -> str:
    """Record which IQ-TREE build produced this tree."""
    try:
        completed = subprocess.run(
            ["iqtree3", "--version"],
            capture_output=True, text=True, check=False,
        )
        lines = (completed.stdout or completed.stderr).strip().splitlines()
        return lines[0] if lines else "unknown"
    except Exception as exc:                                    # noqa: BLE001
        log(f"[phylo_tree] WARNING: could not read IQ-TREE's version ({exc})")
        return "unknown"


def selected_model_from_report(report: Path) -> str | None:
    """Pull the model IQ-TREE actually used out of its .iqtree report.

    Only interesting when the user asked for model selection (MFP): the whole point of
    permitting that kind of automatic behaviour is that the choice is recoverable
    afterwards, so it is copied into the run JSON rather than left buried in a report
    nobody reads.
    """
    if not report.exists():
        return None
    for line in report.read_text(errors="replace").splitlines():
        stripped = line.strip()
        # IQ-TREE writes e.g. "Model of substitution: GTR+F+G4"
        if stripped.startswith("Model of substitution:"):
            return stripped.split(":", 1)[1].strip()
    return None


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
    report_out   = Path(sm.output.report)
    run_json_out = Path(sm.output.run_json)

    prefix          = Path(sm.params.prefix)
    model           = str(sm.params.model)
    model_requested = str(sm.params.model_requested)
    support         = bool(sm.params.support)
    seed            = int(sm.params.seed)
    extra_args      = str(sm.params.extra_args or "").strip()
    threads         = int(sm.threads)

    prefix.parent.mkdir(parents=True, exist_ok=True)
    clear_stale_prefix_files(prefix, log)

    command = [
        "iqtree3",
        "-s", str(aln_in),
        "--seqtype", "DNA",
        "-m", model,
        # Both the thread count and its ceiling are explicit integers: IQ-TREE
        # otherwise reserves the right to scale up on its own.
        "-T", str(threads),
        "--threads-max", str(threads),
        "--seed", str(seed),
        "--prefix", str(prefix),
        # Safe mode uses a slower but numerically more robust likelihood kernel. It
        # engages by itself past 2000 sequences; passing it explicitly means the
        # behaviour does not change when a run happens to cross that count.
        "--safe",
    ]

    if support:
        # Ultrafast bootstrap plus the SH-aLRT test. Off by default because neither
        # UniFrac nor Faith's PD reads support values, and computing them dominates
        # the runtime. IQ-TREE's documented thresholds for a single-gene tree like
        # this one are UFBoot >= 95 and SH-aLRT >= 80.
        command += ["-B", "1000", "--alrt", "1000"]
        log("[phylo_tree] Support enabled: adding -B 1000 --alrt 1000 "
            "(interpretation thresholds for a single-gene tree: UFBoot >= 95, SH-aLRT >= 80)")

    if extra_args:
        command += shlex.split(extra_args)
        log(f"[phylo_tree] Appending extra_args.iqtree: {extra_args}")

    log(f"[phylo_tree] Backend iqtree, model {model} "
        f"(requested: {model_requested}), {threads} thread(s), seed {seed}")
    log(f"[phylo_tree] Running: {' '.join(command)}")

    completed = subprocess.run(command, capture_output=True, text=True)

    # IQ-TREE's console output goes into this rule's log, so a failure is diagnosable
    # even though the tool's own .log file is left undeclared (Snakemake deletes a
    # failed job's declared outputs, which would take that file with it).
    if completed.stdout:
        log("[phylo_tree] --- IQ-TREE output ---")
        log(completed.stdout.rstrip())
        log("[phylo_tree] ----------------------")
    if completed.stderr:
        log("[phylo_tree] --- IQ-TREE stderr ---")
        log(completed.stderr.rstrip())
        log("[phylo_tree] ----------------------")

    if completed.returncode != 0:
        log(f"[phylo_tree] ERROR: IQ-TREE exited with status {completed.returncode}")
        log_fh.close()
        return completed.returncode

    if not tree_out.exists():
        log(f"[phylo_tree] ERROR: IQ-TREE reported success but {tree_out.name} is missing.")
        log_fh.close()
        return 1

    selected_model = selected_model_from_report(report_out)
    if selected_model and selected_model != model:
        log(f"[phylo_tree] Model selection chose: {selected_model}")

    record = {
        "tool": "iqtree3",
        "version": iqtree_version(log),
        "backend": "iqtree",
        "model_requested": model_requested,
        "model_passed": model,
        "model_used": selected_model or model,
        "support": support,
        "threads": threads,
        "seed": seed,
        "extra_args": extra_args,
        "command": " ".join(command),
    }
    run_json_out.parent.mkdir(parents=True, exist_ok=True)
    run_json_out.write_text(json.dumps(record, indent=2) + "\n")

    log(f"[phylo_tree] DONE. Tree written to {tree_out.name}")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
