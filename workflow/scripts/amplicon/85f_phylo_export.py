#!/usr/bin/env python3
"""Validate the inferred tree and write MetaFlux's canonical phylogeny outputs.

Fourth step of the optional 16S phylogeny (see 85_phylogeny.smk). Whichever backend
ran, this script produces the same small set of files, so downstream use in R does not
depend on which tool built the tree.

THE GATES COME FIRST
--------------------
Nothing canonical is written until the tree has passed every check below. A tree that
is subtly wrong is worse than no tree at all, because it looks perfectly usable and
the error only shows up as a diversity number nobody can explain:

  1. Every ASV that went in is a tip, exactly once, and there are no extra tips. This
     is the property the whole module rests on — the tree and the abundance table must
     describe the same feature set, or every metric computed from the pair is
     meaningless.
  2. The Newick parses.
  3. Every branch length is present, finite, and not negative. Negative branch lengths
     break both Faith's PD and UniFrac; missing ones make them undefined.

THE TREE IS EXPORTED UNROOTED, AND ONLY UNROOTED
------------------------------------------------
This is a decision, not an oversight. Faith's PD and UniFrac both depend on where the
root sits. In practice users prune the tree before computing them — decontam against
negative-control samples, abundance filtering of rare ASVs — and a root placed on the
FULL tree (a midpoint in particular) is no longer valid once tips are removed. Rooting
therefore belongs in the user's R session, after pruning and immediately before the
metric (phytools::midpoint.root() or phangorn::midpoint()). No rooted copy is offered:
a file that is only correct until the user prunes — which they will — is a trap.

Inputs (see rule phylo_export in 85_phylogeny.smk)
  aln        : the alignment, from phylo_align — its headers are checked against the
               tree's tips.
  tree       : the backend's own tree file, from phylo_tree.
  input_json : from phylo_input — carries the authoritative ASV ID list, which was
               already verified against 6.taxonomy/asv_table.txt when it was written.
  aln_json,
  tree_json  : the per-step run records, merged into phylogeny.params.json here.

Outputs
  unrooted    : 7.phylogeny/asv_16s.unrooted.nwk — THE deliverable. Nothing in
                MetaFlux reads it; it is for your R session.
  params_json : 7.phylogeny/phylogeny.params.json — the full provenance record
                (resolved strategy and model, thread counts, tool versions, seed,
                input checksums, extra_args). Read by phylo_qc for its report header.
"""
import json
import math
import sys
from collections import Counter
from pathlib import Path

import dendropy
from dendropy.dataio.newickreader import NewickReader


def read_fasta_ids(path: Path) -> list[str]:
    """Return the sequence IDs of a FASTA, in file order."""
    ids: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                label = line[1:].strip()
                if label:
                    ids.append(label.split()[0])
    return ids


def load_tree(path: Path, log) -> dendropy.Tree:
    """Read a Newick tree, or stop with a message naming the file that failed.

    preserve_underscores is essential: without it DendroPy follows the old NEXUS
    convention of turning underscores into spaces, which would silently rewrite every
    ASV_1 label as "ASV 1" and break the match against the abundance table.
    """
    try:
        return dendropy.Tree.get(
            path=str(path),
            schema="newick",
            preserve_underscores=True,
        )
    except NewickReader.NewickReaderDuplicateTaxonError as exc:
        # DendroPy refuses a tree with two tips of the same name before check_tips
        # ever sees it, so report it as what it is — a tip-accounting failure — rather
        # than as a generic parse error.
        log(f"[phylo_export] ERROR: the tree in {path} repeats a tip label; no tree with "
            f"duplicate tips can be matched to the abundance table. {exc}")
        raise SystemExit(1) from exc
    except Exception as exc:                                    # noqa: BLE001
        log(f"[phylo_export] ERROR: could not parse {path} as Newick: {exc}")
        raise SystemExit(1) from exc


def check_tips(tree: dendropy.Tree, expected_ids: list[str], source: str, log) -> bool:
    """Gate 1 — the tree's tips must be exactly the expected ASV set, each once."""
    tip_labels = [leaf.taxon.label for leaf in tree.leaf_node_iter() if leaf.taxon is not None]
    unlabelled = sum(1 for leaf in tree.leaf_node_iter() if leaf.taxon is None)

    if unlabelled:
        log(f"[phylo_export] ERROR: the tree has {unlabelled} unlabelled tip(s).")
        return False

    expected = set(expected_ids)
    found    = set(tip_labels)

    if len(tip_labels) != len(found):
        # Normally unreachable: DendroPy already rejects duplicate tip labels while
        # reading the file (load_tree reports that case). Kept as the second line of
        # defence in case a future reader lets them through.
        duplicates = sorted(
            label for label, count in Counter(tip_labels).items() if count > 1
        )[:10]
        log(f"[phylo_export] ERROR: the tree repeats tip label(s) (up to 10): {duplicates}")
        return False

    if found != expected:
        missing = sorted(expected - found)[:10]
        extra   = sorted(found - expected)[:10]
        log(f"[phylo_export] ERROR: the tree's tips do not match the {source} ASV set.")
        log(f"[phylo_export]   expected {len(expected)} tips, found {len(found)}")
        if missing:
            log(f"[phylo_export]   missing from the tree (up to 10): {missing}")
        if extra:
            log(f"[phylo_export]   present in the tree but not in the ASV set (up to 10): {extra}")
        return False

    log(f"[phylo_export] Tip accounting OK: all {len(expected)} ASVs appear exactly once.")
    return True


def check_branch_lengths(tree: dendropy.Tree, log) -> bool:
    """Gate 3 — every branch length present, finite, and non-negative.

    The root's own edge legitimately has no length (there is nothing above it), so it
    is skipped. Everything else must carry a usable number: a missing or negative
    branch length makes Faith's PD and UniFrac either undefined or wrong.
    """
    missing = 0
    negative: list[str] = []
    non_finite: list[str] = []

    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        length = node.edge_length
        label = node.taxon.label if node.taxon is not None else "<internal node>"

        if length is None:
            missing += 1
        elif not math.isfinite(length):
            non_finite.append(label)
        elif length < 0:
            negative.append(label)

    if missing:
        log(f"[phylo_export] ERROR: {missing} branch(es) have no length. Faith's PD and "
            "UniFrac both need branch lengths on every edge.")
    if non_finite:
        log(f"[phylo_export] ERROR: non-finite branch length(s) at (up to 10): {non_finite[:10]}")
    if negative:
        log(f"[phylo_export] ERROR: negative branch length(s) at (up to 10): {negative[:10]}. "
            "These break both Faith's PD and UniFrac.")

    if missing or non_finite or negative:
        return False

    log("[phylo_export] Branch lengths OK: all present, finite and non-negative.")
    return True


def write_newick(tree: dendropy.Tree, path: Path) -> None:
    """Write a Newick file for downstream use in R.

    suppress_rooting keeps DendroPy's [&R]/[&U] rooting comment out of the file: ape's
    read.tree copes with it, but plenty of other readers do not, and the rooting state
    is documented by the filename and the docs rather than by a comment.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(
        path=str(path),
        schema="newick",
        suppress_rooting=True,
        unquoted_underscores=True,
    )


def main() -> int:
    sm = snakemake  # noqa: F821

    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = log_path.open("w")

    def log(message: str) -> None:
        print(message, file=log_fh, flush=True)
        print(message, file=sys.stderr, flush=True)

    aln_in        = Path(sm.input.aln)
    tree_in       = Path(sm.input.tree)
    input_json_in = Path(sm.input.input_json)
    aln_json_in   = Path(sm.input.aln_json)
    tree_json_in  = Path(sm.input.tree_json)

    unrooted_out    = Path(sm.output.unrooted)
    params_json_out = Path(sm.output.params_json)

    backend    = str(sm.params.backend)
    phylo_dir  = Path(sm.params.phylo_dir)
    extra_args = dict(sm.params.extra_args)

    input_record = json.loads(input_json_in.read_text())
    expected_ids = input_record["asv_ids"]

    log(f"[phylo_export] Validating the {backend} tree against {len(expected_ids)} eligible ASVs")

    tree = load_tree(tree_in, log)

    # Gate 1, run twice against two independent records of what the ASV set is: the
    # alignment that was actually fed to the backend, and the ID list recorded by
    # phylo_input (already checked against 6.taxonomy/asv_table.txt when written).
    # Agreement across both is what lets the docs promise that the tree's tips and the
    # abundance table are the same feature set.
    aligned_ids = read_fasta_ids(aln_in)
    ok = check_tips(tree, expected_ids, "6.taxonomy", log)
    if ok and sorted(aligned_ids) != sorted(expected_ids):
        log("[phylo_export] ERROR: the alignment and the recorded ASV set disagree — "
            "the alignment step and the taxonomy tables are out of step.")
        ok = False

    # Gate 3 (gate 2, "does it parse", was passed by load_tree above).
    if ok:
        ok = check_branch_lengths(tree, log)

    if not ok:
        log("[phylo_export] Refusing to write the canonical tree: it failed at least one "
            "correctness check above.")
        log_fh.close()
        return 1

    # ── Canonical output ──────────────────────────────────────────────────────
    # Unrooted, and only unrooted — no rooted copy is offered. Users prune the tree in
    # R first (decontam against negative controls, abundance filtering of rare ASVs),
    # and any root placed on the full ASV set stops being valid the moment a tip is
    # removed. A file that is only correct until the user does the thing they will
    # certainly do would be a trap, not a convenience.
    write_newick(tree, unrooted_out)
    log(f"[phylo_export] Wrote the tree: {unrooted_out.name} (UNROOTED — the only form "
        "MetaFlux exports; root it in R after you have pruned to your final ASV set)")

    # Switching backend leaves the previous one's directory in
    # place. That is deliberate — it is useful for comparison — but say so, so nobody
    # mistakes it for output of this run.
    for candidate in sorted(phylo_dir.glob("*")):
        if candidate.is_dir() and candidate.name in {"iqtree", "fasttree", "raxml-ng"} \
                and candidate.name != backend:
            log(f"[phylo_export] NOTE: {candidate.name}/ is left over from an earlier run "
                f"with a different backend. The canonical outputs above come from "
                f"{backend}.")

    # ── Provenance ────────────────────────────────────────────────────────────
    # Everything resolved at run time, in one file: the point is that any choice this
    # module made for you can be read back and pinned explicitly later.
    params_record = {
        "module": "phylogeny",
        "marker": "16S",
        "backend": backend,
        "input": input_record,
        "alignment": json.loads(aln_json_in.read_text()),
        "tree": json.loads(tree_json_in.read_text()),
        "extra_args": extra_args,
        "outputs": {
            "unrooted": unrooted_out.name,
        },
        "note": (
            "MetaFlux computes no diversity statistics. Compute Faith's PD, UniFrac or "
            "anything else from these files in your own downstream tools."
        ),
    }
    params_json_out.parent.mkdir(parents=True, exist_ok=True)
    params_json_out.write_text(json.dumps(params_record, indent=2) + "\n")

    log(f"[phylo_export] DONE. Provenance recorded in {params_json_out.name}")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
