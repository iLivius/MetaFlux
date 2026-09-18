#!/usr/bin/env python3
"""Write the phylogeny QC report — or record why no tree was built.

Last step of the optional 16S phylogeny (see 85_phylogeny.smk), and the only one whose
output `rule all` asks for directly. Two jobs:

  1. When a tree was built, summarise it in a form that answers "should I trust this?"
     without opening the tree in R.
  2. When it was skipped because there were too few ASVs, say so, in a file, rather
     than leaving a silently missing tree for the user to puzzle over.

THE LONG-BRANCH REPORT IS THE POINT
-----------------------------------
Of everything in here, the long-branch report is the number worth looking at. The known
failure mode of a de novo tree built from short amplicon reads is that a few sequences
only distantly related to the rest — an off-target amplification, a contaminant the
taxonomy filter did not catch — end up on very long branches. Faith's PD is an absolute
sum of branch length, so a long branch adds to it directly; unweighted UniFrac is even
more sensitive, and a couple of such tips can manufacture a separation between sample
groups that looks exactly like real biology. This is documented in the literature: in
one published case three low-abundance archaeal ASVs on a long branch produced a fake
UniFrac cluster separation that disappeared once the tree was built properly.

So: any tip sitting on an unusually long branch is flagged here, and the longest few
are always listed whether or not they cross the threshold. Two measures are reported —
the pendant edge (the branch leading to that tip alone) as the primary flag, and
root-to-tip alongside it. See long_branch_report below for why the obvious measure is
not the reliable one.

Flagged, never removed. Dropping tips would desynchronise the tree from the abundance
table, which is the one thing this module must never do. What to do about a flagged tip
is a judgement call — check what it was classified as, consider whether the contaminant
filter should have caught it, and rerun your analysis without it to see whether your
conclusion depends on it.

Root-to-tip distance needs a root, and the exported tree is deliberately unrooted, so a
midpoint-rooted COPY is made in memory purely for this measurement. That copy is never
written anywhere; it does not affect the exported tree.

Inputs (see rule phylo_qc in 85_phylogeny.smk)
  input_json  : from phylo_input — always present. Its n_eligible is what decided
                whether the rest of these inputs exist at all.
  aln,
  unrooted,
  params_json : present only when a tree was actually built (4 or more eligible ASVs).

Output
  qc_json : stats/phylogeny/phylogeny_qc.json — for a human to read after a run.
            Nothing else in MetaFlux reads it back.
"""
import json
import math
import statistics
import sys
from pathlib import Path

import dendropy


def alignment_stats(path: Path) -> dict:
    """Summarise the alignment: how many sequences, how wide, how gappy, how variable.

    Gap fraction and the count of variable columns are the quick sanity check on
    whether the alignment is doing its job. A same-region 16S ASV set should align
    tightly, so a high gap fraction usually means something got in that does not belong
    — which is the same story the long-branch report tells from the other end.
    """
    labels: list[str] = []
    sequences: list[str] = []
    current: list[str] = []

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
                labels.append(line[1:].strip().split()[0])
            else:
                current.append(line.strip())
    if current:
        sequences.append("".join(current))

    if not sequences:
        return {"n_sequences": 0}

    width = len(sequences[0])
    total_cells = width * len(sequences)
    gap_count = sum(seq.count("-") for seq in sequences)

    variable_columns = 0
    constant_columns = 0
    for column_index in range(width):
        column = {seq[column_index].upper() for seq in sequences}
        bases = column - {"-", "N", "?"}
        if len(bases) > 1:
            variable_columns += 1
        elif len(bases) == 1 and "-" not in column:
            constant_columns += 1

    return {
        "n_sequences": len(sequences),
        "alignment_length": width,
        "gap_fraction": round(gap_count / total_cells, 4) if total_cells else 0.0,
        "variable_columns": variable_columns,
        "constant_columns": constant_columns,
    }


def long_branch_report(tree_path: Path, multiple: float, pendant_fence_iqr: float, log) -> dict:
    """Flag tips on unusually long branches, by two measures.

    TWO MEASURES, AND WHY THE SECOND ONE EXISTS
    -------------------------------------------
    The obvious measure is root-to-tip distance, and it is reported here. But it has a
    weakness that shows up on exactly the case this report is meant to catch, and the
    weakness is structural rather than bad luck.

    The exported tree is unrooted, so measuring anything "from the root" requires
    picking one, and midpoint is the natural choice. Midpoint rooting places the root
    halfway along the LONGEST path in the tree — and when one tip is on a very long
    branch, that longest path runs straight through it. So the offending branch gets
    split across the root, half on each side, and every tip ends up at a broadly
    similar distance from it. The measure flattens out precisely where the signal is.

    Measured on this module's own synthetic test set, with one deliberately divergent
    sequence planted among 30 ordinary ones:

        root-to-tip    1.05x the median   (would NOT be flagged at 1.5x)
        pendant edge  22.71x the median   (unmistakable)

    So the PENDANT EDGE — the length of the single branch leading to that tip, and
    nothing else — is the primary flag here. It needs no root at all, which makes it
    the right measure for an unrooted tree, and it is the direct expression of what a
    contaminant looks like: a sequence far from everything else, joined to the rest of
    the tree by one long branch.

    Root-to-tip is still reported alongside it, for continuity with how this artifact
    is described in the literature and because a clade of several related off-target
    sequences shows up better there — each member's own pendant edge can be short even
    though the whole group hangs off one long stem.

    THE TWO THRESHOLD RULES DIFFER, ON PURPOSE
    ------------------------------------------
    Root-to-tip distances are all measured from the same root, so on a tree with a
    clear centre they cluster and a small multiple of the median separates signal from
    noise — 1.5x, which is roughly the ratio in the published example (an artificial
    branch at 1.43 against a mean tip-to-root of 0.94). How tightly they cluster
    depends on the backend and on which optimum the search found: on the same 211
    ASVs the rule flagged 9 tips on one IQ-TREE tree, 17 on another, 18 on FastTree's
    and 58 on RAxML-NG's, while the pendant fence flagged 6-7 every time. That is why
    root-to-tip is the secondary measure: read the pendant flags first.

    Pendant edges CANNOT use a multiple of the median, and this was learned on real
    data, not reasoned out. A real ASV set contains many near-identical sequences —
    on the 211-ASV 16S test set a quarter of all pendant edges are ~1e-6 — so the
    median is close to zero and any multiple of it lands deep inside normal
    variation: "5x the median" flagged 81 of 211 tips there, which is a useless
    report. The rule used instead is the standard boxplot outlier fence,

        threshold = Q3 + k * IQR        (k = pendant_fence_iqr, default 3)

    which does not care how small the lower half of the distribution is. On the same
    real data it flags 7 tips, and they are the right ones: singleton divergent
    lineages (an Acidobacteriota, a Bdellovibrionota) in a mostly-Pseudomonadota
    community — genuinely long branches, worth a look, not artifacts of the
    threshold. On the synthetic test set the fence still isolates the planted
    contaminant (0.33 against a fence of ~0.05) exactly as the old rule did.

    Small trees are a blind spot of the fence: with six or fewer tips the largest
    value takes part in the Q3 interpolation itself, so Q3 + 3*IQR always lies above
    it and nothing can be flagged (checked numerically: one planted value of 10 among
    values of 1e-6 is not flagged at n <= 6 and is flagged from n = 7). Below 7 tips
    the log says so instead of reporting "no outliers", and the longest_tips list and
    the root-to-tip rule carry the diagnosis.

    Both rules only decide what gets called out in the log; the five longest tips by
    each measure are always listed in the report either way, so nothing is hidden by
    a threshold being too strict.
    """
    tree = dendropy.Tree.get(
        path=str(tree_path), schema="newick", preserve_underscores=True
    )

    # Pendant edges first: read straight off the unrooted tree, no rooting involved.
    pendant: dict[str, float] = {
        leaf.taxon.label: (leaf.edge_length or 0.0)
        for leaf in tree.leaf_node_iter()
        if leaf.taxon is not None
    }

    if not pendant:
        return {"n_tips": 0}

    # Root-to-tip needs a root, so make a midpoint-rooted copy in memory. The exported
    # tree is untouched by this — nothing here is written back.
    rooted = dendropy.Tree.get(
        path=str(tree_path), schema="newick", preserve_underscores=True
    )
    rooted.reroot_at_midpoint(update_bipartitions=False)

    root_to_tip: dict[str, float] = {}
    for leaf in rooted.leaf_node_iter():
        if leaf.taxon is None:
            continue
        distance = 0.0
        node = leaf
        while node.parent_node is not None:
            distance += node.edge_length or 0.0
            node = node.parent_node
        root_to_tip[leaf.taxon.label] = distance

    def summarise(values: dict[str, float], key_name: str, threshold: float,
                  rule: str) -> dict:
        median = statistics.median(values.values())
        flagged = sorted(
            ((label, value) for label, value in values.items() if value > threshold),
            key=lambda pair: pair[1], reverse=True,
        )
        longest = sorted(values.items(), key=lambda pair: pair[1], reverse=True)[:5]
        return {
            "median": round(median, 6),
            "flag_rule": rule,
            "flag_threshold": round(threshold, 6) if math.isfinite(threshold) else None,
            "n_flagged": len(flagged),
            # The list is capped at 50 entries so the report stays readable on a tree
            # where a rule misfires; n_flagged is always the true count and the flag
            # below says when the list is shorter than it.
            "flagged_tips_truncated": len(flagged) > 50,
            "flagged_tips": [
                {"asv_id": label, key_name: round(value, 6)}
                for label, value in flagged[:50]
            ],
            "longest_tips": [
                {"asv_id": label, key_name: round(value, 6)}
                for label, value in longest
            ],
        }

    # Pendant fence: Q3 + k*IQR (see the docstring for why not a multiple of the
    # median). With every pendant edge identical (all-duplicate input) the IQR is 0
    # and the fence collapses to Q3 — then nothing exceeds it, which is right.
    pendant_values = sorted(pendant.values())
    q1 = statistics.quantiles(pendant_values, n=4)[0] if len(pendant_values) > 1 else pendant_values[0]
    q3 = statistics.quantiles(pendant_values, n=4)[2] if len(pendant_values) > 1 else pendant_values[0]
    pendant_fence = q3 + pendant_fence_iqr * (q3 - q1)

    r2t_median = statistics.median(root_to_tip.values())
    r2t_threshold = r2t_median * multiple if r2t_median > 0 else float("inf")

    pendant_summary = summarise(pendant, "pendant_edge", pendant_fence,
                                f"Q3 + {pendant_fence_iqr}*IQR")
    r2t_summary     = summarise(root_to_tip, "root_to_tip", r2t_threshold,
                                f"{multiple}x median")

    if pendant_summary["n_flagged"]:
        worst = pendant_summary["flagged_tips"][0]
        log(f"[phylo_qc] LONG BRANCHES: {pendant_summary['n_flagged']} tip(s) sit on a "
            f"pendant edge above the outlier fence "
            f"(Q3 + {pendant_fence_iqr}*IQR = {pendant_fence:.4f}). Worst: "
            f"{worst['asv_id']} at {worst['pendant_edge']}.")
        log("[phylo_qc] These are the tips that contribute most branch length per tip "
            "to Faith's PD and weigh most in unweighted UniFrac. Check what they were "
            "classified as: in an environmental sample they are often genuine divergent "
            "lineages with a single representative; an off-target or contaminant the "
            "taxonomy filter missed looks the same. Reported only; nothing has been "
            "removed.")
    elif len(pendant) < 7:
        log(f"[phylo_qc] Only {len(pendant)} tips: the pendant-edge fence cannot flag "
            "anything at this size (see long_branch_report). Read the longest_tips list "
            "in the report instead.")
    else:
        log(f"[phylo_qc] No pendant edge above the outlier fence "
            f"(Q3 + {pendant_fence_iqr}*IQR = {pendant_fence:.4f}) — no long-branch "
            "outliers.")

    if r2t_summary["n_flagged"]:
        log(f"[phylo_qc] ({r2t_summary['n_flagged']} tip(s) also exceed {multiple}x the "
            "median root-to-tip distance.)")

    return {
        "n_tips": len(pendant),
        # The primary measure: rooting-independent, and the one that actually catches a
        # single divergent sequence. See this function's docstring for the numbers.
        "pendant_edge": pendant_summary,
        # Reported alongside; better at spotting a whole off-target CLADE hanging off
        # one long stem, weaker at spotting a single long branch (midpoint rooting
        # splits that branch across the root).
        "root_to_tip_midpoint_rooted": r2t_summary,
        "n_flagged": pendant_summary["n_flagged"],
        "note": (
            "Pendant edge = the length of the single branch leading to that tip, read "
            "off the unrooted tree. Root-to-tip is measured on a midpoint-rooted copy "
            "made in memory for this report only; the exported tree stays unrooted. "
            "Flagged tips are reported, never removed — removing them here would "
            "desynchronise the tree from the abundance table. With fewer than 7 tips "
            "the pendant-edge fence cannot flag anything; read longest_tips instead."
        ),
    }


def main() -> int:
    sm = snakemake  # noqa: F821

    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = log_path.open("w")

    def log(message: str) -> None:
        print(message, file=log_fh, flush=True)
        print(message, file=sys.stderr, flush=True)

    input_json_in = Path(sm.input.input_json)
    qc_json_out   = Path(sm.output.qc_json)

    min_asvs = int(sm.params.min_asvs)
    backend  = str(sm.params.backend)
    multiple = float(sm.params.long_branch_multiple)
    pendant_fence_iqr = float(sm.params.long_branch_pendant_fence_iqr)

    input_record = json.loads(input_json_in.read_text())
    n_eligible = input_record["n_eligible"]

    # The skip case. phylo_qc's input function did not ask for a tree, so none was
    # built and none of the other inputs exist. Record why, so the absence of a tree is
    # explained in a file rather than being something the user has to work out.
    if n_eligible < min_asvs:
        log(f"[phylo_qc] No tree was built: {n_eligible} eligible ASV(s), minimum is "
            f"{min_asvs}.")
        report = {
            "skipped": True,
            "reason": (
                f"Only {n_eligible} eligible ASV(s) after taxonomy and contaminant "
                f"filtering; at least {min_asvs} are needed. With fewer than 4 tips there "
                "is only one possible unrooted topology, so there is nothing to infer."
            ),
            "n_eligible": n_eligible,
            "min_asvs": min_asvs,
            "backend_configured": backend,
        }
        qc_json_out.parent.mkdir(parents=True, exist_ok=True)
        qc_json_out.write_text(json.dumps(report, indent=2) + "\n")
        log(f"[phylo_qc] DONE. Skip recorded in {qc_json_out.name}")
        log_fh.close()
        return 0

    # The normal case.
    aln_in         = Path(sm.input.aln)
    unrooted_in    = Path(sm.input.unrooted)
    params_json_in = Path(sm.input.params_json)

    params_record = json.loads(params_json_in.read_text())

    log(f"[phylo_qc] Reporting on the {backend} tree of {n_eligible} ASVs")

    aln_summary = alignment_stats(aln_in)
    log(f"[phylo_qc] Alignment: {aln_summary.get('n_sequences')} sequences x "
        f"{aln_summary.get('alignment_length')} columns, "
        f"gap fraction {aln_summary.get('gap_fraction')}, "
        f"{aln_summary.get('variable_columns')} variable columns")

    branches = long_branch_report(unrooted_in, multiple, pendant_fence_iqr, log)

    # Tip accounting is re-stated here rather than merely assumed: phylo_export already
    # refuses to write a tree whose tips do not match, so this is the report's record of
    # a check that has passed, in the one file a user actually opens after a run.
    tip_accounting = {
        "n_eligible_asvs": n_eligible,
        "n_tips": branches.get("n_tips", 0),
        "all_asvs_present_exactly_once": branches.get("n_tips", 0) == n_eligible,
        "note": (
            "Enforced by phylo_export before the tree was written: it refuses to produce "
            "asv_16s.unrooted.nwk unless every ASV in 6.taxonomy/asv_table.txt appears "
            "exactly once as a tip."
        ),
    }

    report = {
        "skipped": False,
        "marker": "16S",
        "backend": backend,
        "aligner_strategy": params_record.get("alignment", {}).get("strategy_resolved"),
        "model_used": params_record.get("tree", {}).get("model_used"),
        "tip_accounting": tip_accounting,
        "alignment": aln_summary,
        "long_branches": branches,
        "reproducibility_note": (
            "FastTree here is single-threaded with no random component: the same "
            "alignment file gives a byte-identical tree. The alignment's record order "
            "still depends on ASV numbering, which can differ between runs for equally "
            "abundant ASVs; see docs/amplicon/phylogeny.md, Reproducibility."
            if backend == "fasttree" else
            "Exactly reproducible only for the same input file, the same seed and one "
            "thread. A different thread count, or a different order of the input "
            "sequences (which happens between runs when equally abundant ASVs are "
            "numbered differently), sends the search to a different, about equally good "
            "tree: on the 16S test set RF 0.37 and patristic r 0.93 between two such "
            "trees, per-sample PD r 0.999. Compare trees by patristic distance or PD, "
            "never by diff; see docs/amplicon/phylogeny.md, Reproducibility."
        ),
        "downstream_note": (
            "MetaFlux computes no diversity statistics. Prune this tree to your filtered "
            "ASV set in R (ape::drop.tip or phyloseq::prune_taxa), then root it, then "
            "compute Faith's PD or UniFrac. Root after pruning, not before."
        ),
    }
    qc_json_out.parent.mkdir(parents=True, exist_ok=True)
    qc_json_out.write_text(json.dumps(report, indent=2) + "\n")

    log(f"[phylo_qc] DONE. Report written to {qc_json_out.name}")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
