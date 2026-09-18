#!/usr/bin/env python3
"""Infer the ASV tree with FastTree 2 (alternative backend).

Third step of the optional 16S phylogeny (see 85_phylogeny.smk), for users who choose
phylogeny.backend: fasttree. FastTree computes an approximate maximum-likelihood tree
— heuristic topology search plus ML branch lengths — and is dramatically faster than
IQ-TREE or RAxML-NG: about 93 seconds and under 200 MB for 10,000 ASVs on one core.
It is the backend QIIME 2 and nf-core/ampliseq use by default.

WHY -gtr IS PASSED EXPLICITLY
-----------------------------
FastTree's default nucleotide model is Jukes-Cantor, which assumes every substitution
is equally likely and all four bases equally frequent — not defensible for 16S. It has
no model selection of any kind; the only choices are JC and GTR. Worth knowing: QIIME 2
invokes FastTree with no model flags at all, so the default 16S tree from QIIME 2 and
ampliseq is a Jukes-Cantor tree, and there is no supported way to change it there.
MetaFlux passes -gtr, and honours phylogeny.model: jc if you deliberately want JC.

SINGLE-THREADED, ALWAYS
-----------------------
There is a parallel build, FastTreeMP, and this script never calls it. Its own
documentation states both that "FastTreeMP will not give exactly the same results as
FastTree" and that "FastTreeMP is not deterministic" — and its speedup does not extend
to the maximum-likelihood phase beyond about three cores. Given the serial binary is
already fast, there is nothing to buy and reproducibility to lose. OMP_NUM_THREADS is
pinned to 1 as well, so even a FastTree built with OpenMP support cannot quietly
spread out.

-gamma AND WHAT IT MEANS DOWNSTREAM
-----------------------------------
-gamma rescales the tree's branch lengths by a single global factor — on the 16S test
set total tree length went from 12.15 without it to 15.30 with it, +26% — and costs
about 5% more runtime. The consequence is worth
stating plainly, because it lands differently on the two metrics people use this tree
for: a global scale factor CANCELS in UniFrac, which is a ratio of branch lengths, but
does NOT cancel in Faith's PD, which is an absolute sum.

-noml is never used: it is documented to return slightly negative branch lengths,
which would break both metrics.

Inputs (see rule phylo_tree, fasttree branch, in 85_phylogeny.smk)
  aln : 7.phylogeny/asv_16s.aln.fasta, from phylo_align.

Outputs
  tree     : <backend dir>/asv_16s.nwk — read next by phylo_export, which validates it
             and writes the canonical unrooted Newick.
  run_json : <backend dir>/tree_run.json — resolved model flags, version, threads,
             folded into phylogeny.params.json by phylo_export.
"""
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path


def assert_double_precision(log) -> str | None:
    """Refuse to run against a single-precision FastTree build.

    This matters more than it sounds, and it matters specifically for what this module
    exports. A single-precision build cannot represent very short branches: it floors
    them at 5e-4 where a double-precision build reaches 5e-9. On a measured 2,000-ASV
    test that put 44% of branches on the floor and added about 13% to the total tree
    length — and total tree length IS Faith's PD, by definition. The error grows with
    how many closely related sequences a sample carries, so it is largest exactly where
    amplicon data is richest.

    Bioconda has built FastTree with double precision since 2017 and 2.2.0 makes it
    the default, so the shipped environment is safe; this guards the case where the
    binary comes from somewhere else.

    FastTree prints its banner to STDERR, so stderr must be captured — checking stdout
    alone would find an empty string and reject every build.

    Returns the banner line, which doubles as FastTree's version statement (it has no
    --version flag), or None if the check failed. The caller records it in the run JSON:
    the rule's log is truncated on every re-run, so pointing at the log would leave the
    provenance file with no version at all.
    """
    try:
        completed = subprocess.run(
            ["FastTree", "-help"],
            capture_output=True, text=True, check=False,
        )
    except FileNotFoundError:
        log("[phylo_tree] ERROR: the FastTree binary was not found in this environment.")
        return None

    banner = (completed.stderr or "") + (completed.stdout or "")
    first_line = banner.strip().splitlines()[0] if banner.strip() else ""
    log(f"[phylo_tree] FastTree reports: {first_line}")

    if "double precision" not in first_line.lower():
        log("[phylo_tree] ERROR: this FastTree build is not double precision. Short "
            "branch lengths would be floored at 5e-4 instead of 5e-9, inflating total "
            "tree length (and therefore Faith's PD) by roughly 13% on a typical ASV set. "
            "Use the bioconda build (workflow/envs/fasttree.yaml), which is compiled "
            "with -DUSE_DOUBLE.")
        return None
    return first_line


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

    model           = str(sm.params.model).lower()
    model_requested = str(sm.params.model_requested)
    support         = bool(sm.params.support)
    seed            = int(sm.params.seed)
    extra_args      = str(sm.params.extra_args or "").strip()

    version_banner = assert_double_precision(log)
    if version_banner is None:
        log_fh.close()
        return 1

    command = ["FastTree", "-nt"]

    # GTR unless the user explicitly asked for Jukes-Cantor. Honouring `jc` matters:
    # silently substituting GTR for a requested JC would hand back a tree that is not
    # what was asked for, with nothing to show it.
    if model == "gtr":
        command.append("-gtr")
    else:
        log("[phylo_tree] model: jc — running FastTree's Jukes-Cantor default "
            "(-gtr omitted, as requested)")

    # Applied under both models: it is a branch-length rescaling, independent of the
    # substitution model. See the module docstring for what a global rescaling does to
    # UniFrac versus Faith's PD.
    command.append("-gamma")

    if not support:
        # SH-like local support values are FastTree's only use of randomness, and they
        # cost about 19% of the runtime. Nothing downstream in this module's intended
        # use reads them, so they are off by default — which also makes the run
        # completely deterministic.
        command.append("-nosupport")
    else:
        # Support values come from resampling alignment sites, which is the one place
        # FastTree uses its random number generator. Seed it from amplicon.seed rather
        # than leaving FastTree on its own built-in default, so this run is reproducible
        # for the same reason every other backend is. (-seed is documented under
        # `FastTree -expert`.)
        command += ["-seed", str(seed)]
        log("[phylo_tree] Support enabled: keeping FastTree's SH-like local support "
            f"values (NOT bootstrap; 0-1 scale), seeded with {seed}")

    if extra_args:
        command += shlex.split(extra_args)
        log(f"[phylo_tree] Appending extra_args.fasttree: {extra_args}")

    command.append(str(aln_in))

    # The serial binary is used, so this is belt-and-braces — but a FastTree built with
    # OpenMP would otherwise read this variable from the environment and could spread
    # across cores, which is the documented non-deterministic case.
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"

    log(f"[phylo_tree] Backend fasttree, model {model} (requested: {model_requested}), "
        "single-threaded, OMP_NUM_THREADS=1")
    log(f"[phylo_tree] Running: {' '.join(command)} > {tree_out.name}")

    tree_out.parent.mkdir(parents=True, exist_ok=True)
    # FastTree writes the Newick to stdout and its progress to stderr.
    with tree_out.open("w") as tree_fh:
        completed = subprocess.run(
            command, stdout=tree_fh, stderr=subprocess.PIPE, text=True, env=env
        )

    if completed.stderr:
        log("[phylo_tree] --- FastTree output ---")
        log(completed.stderr.rstrip())
        log("[phylo_tree] -----------------------")

    if completed.returncode != 0:
        log(f"[phylo_tree] ERROR: FastTree exited with status {completed.returncode}")
        log_fh.close()
        return completed.returncode

    if tree_out.stat().st_size == 0:
        log(f"[phylo_tree] ERROR: FastTree reported success but {tree_out.name} is empty.")
        log_fh.close()
        return 1

    record = {
        "tool": "FastTree",
        # FastTree has no --version flag; the -help banner IS its version statement, and
        # the double-precision assert above already had to read it.
        # The banner ends in a colon ("FastTree 2.2.0 Double precision:"); drop it.
        "version": version_banner.rstrip(":").strip(),
        "backend": "fasttree",
        "model_requested": model_requested,
        "model_used": "GTR+gamma" if model == "gtr" else "JC+gamma",
        "support": support,
        "threads": 1,
        "seed": seed,
        # With -nosupport (the default) FastTree uses no randomness at all, so no seed is
        # passed and none is needed — the run is deterministic regardless. The seed is
        # only actually given to FastTree when support values are computed.
        "seed_used": bool(support),
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
