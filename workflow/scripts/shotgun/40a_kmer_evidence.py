#!/usr/bin/env python3
"""Score every Kraken2 species call by how much *distinct sequence* supports it,
and optionally remove the ones that fail.

WHY THIS EXISTS
---------------
Kraken2's read count tells you how many reads it chose to label with a taxon. It
does not tell you how much sequence unique to that organism actually turned up in
the sample, and those are different questions. A species that is really present
contributes many DIFFERENT species-specific sequences. A spurious call is typically
a small set of sequences — conserved, repetitive, or shared with a relative that
really is abundant — matched over and over. In a read count the two look identical.

Running Kraken2 with --report-minimizer-data (rule kraken2, 35_kraken2.smk) gives
the second number. A minimizer is the representative k-mer Kraken2 picks from each
window, so counting distinct minimizers is a sampled way of counting distinct
k-mers — see the conversion note under "THRESHOLD UNITS" below, the units matter.

WHAT THE MINIMIZER COLUMNS ACTUALLY COUNT — read before interpreting them
------------------------------------------------------------------------
The obvious reading is wrong, and it changes what the numbers mean. The columns are
NOT "the minimizers carried by the reads assigned to this taxon". Kraken2 attributes
each minimizer OCCURRENCE it looks up to whichever taxon that minimizer maps to in
the database (its lowest common ancestor), pooled over every read in the sample,
regardless of where each read finally got assigned. The report then rolls those
counts up the tree cumulatively, the same way it does the read counts.

The give-away in any real report is taxa with a couple of assigned reads carrying
hundreds of thousands of minimizer counts. In the MetaFlux test data one Streptomyces
row has 2 clade reads and 345,108 total minimizers — 737x more than the ~468 k-mer
positions two read pairs physically contain. Under the "this taxon's own reads"
reading that is impossible; under LCA attribution it is unremarkable.

So column 5 for a species means: how many DISTINCT database minimizers specific to
that species (or to something below it) were seen anywhere in this sample. That is a
statement about how much species-specific sequence is present — which is exactly the
quantity KrakenUniq calls unique k-mers, so the published thresholds below do
transfer once the units are converted. It is NOT a per-read genome-coverage measure,
and describing it as "how much of the genome the reads covered" would be wrong.

Column 4 divided by column 5 (call it duplicity) is then the average number of times
each of those species-specific minimizers was seen. Low means the evidence is spread
across much distinct sequence; high means a small set of positions was hit again and
again. Useful to eyeball, and unlike the raw counts it barely tracks abundance
(measured Pearson r = +0.20 against log reads, versus +0.71 for the distinct count) —
but it is a WORSE filter, so it is reported and not thresholded (see below).

WHAT THIS SCRIPT DOES
---------------------
Input:  one sample's Kraken2 report, 8 columns because of --report-minimizer-data
        (rule kraken2, 35_kraken2.smk).

Output, all three written on every run whether or not the gate is switched on:
  1. {sample}_species_evidence.tsv — one row per taxon at the target rank, with its
     read count, its minimizer counts, the derived ratios, and the verdict the
     current thresholds give it. This is the table to look at BEFORE trusting a
     threshold on your own data.
  2. {sample}_report_gated.txt — the Kraken2 report with the failing taxa removed.
     Written unconditionally so you can diff it against the original and see exactly
     what the gate would do. Whether Bracken actually reads this file or the original
     is decided by shotgun.kmer_evidence.enabled (see rule bracken, 45_abundance.smk).
  3. {sample}_evidence.png — reads against distinct minimizers, log-log, one point per
     taxon, coloured by verdict, with the threshold drawn on. An outlier becomes
     something you see rather than something you compute.

HOW A TAXON IS REMOVED — AND WHY IT IS NOT JUST A DELETED LINE
--------------------------------------------------------------
Two things have to be right or the pruned report quietly corrupts Bracken's arithmetic.

(a) Whole subtree, never a bare row. Bracken (est_abundance.py) rebuilds the taxonomy
    tree from the INDENTATION of the report, attaching each row to the last row seen
    one indent level above it. Delete a species row but leave its strain (S1) children
    behind and those children get re-parented onto whichever species happens to
    precede them at the same depth — silently inflating an innocent neighbour. So when
    a taxon fails, every deeper-indented row beneath it goes with it.

(b) The reads are handed back to the parent rather than deleted. Column 2 of a Kraken2
    report is the clade total (this taxon plus everything under it) and column 3 is
    the reads sitting on this node alone. Bracken redistributes a parent's column-3
    reads down onto the target rank. So when we remove a species we add its whole
    clade total to its parent's column 3, and Bracken redistributes those reads among
    the siblings that DID pass. The parent's column 2 is untouched (the reads never
    left the clade), so the report's own arithmetic still balances: for every node,
    column 2 == column 3 + sum of the children's column 2. Verified on the test data —
    zero violations, zero indentation gaps, and no column other than a parent's
    column 3 ever changes.

    Two honest caveats, both measured rather than assumed:

    * The reads do not all end up with close relatives. The rationale for crediting
      the parent is that a spurious call usually borrows its evidence from a
      neighbour, so the neighbour is the likely true source. That holds for most but
      not all of the traffic: on the four test samples 70% of the reads gained by
      surviving taxa can be matched to a loss inside the same genus, and the other
      30% crosses genus or family boundaries.

    * Some reads ARE lost, despite the credit. If every species under a genus fails
      the check, Bracken has nothing left at the target rank to redistribute that
      genus's reads to, and drops them. On the test data that cost 14,606 of
      9,299,837 reads (-0.16%), with 27-38 genera per sample emptied completely.
      Deleting the rows without crediting the parent would lose far more; crediting
      is the better of the two, not a guarantee of conservation.

    Discarding the reads outright was deliberately not implemented: it would break the
    column-2 balance at every ancestor up to the root, and throw away real sequence
    over what is a naming problem at one rank.

THRESHOLD UNITS — THE EASY MISTAKE
----------------------------------
Every published threshold (KrakenUniq 2018, aMeta 2023, and Oskolkov 2026 in
Front. Microbiol. 17:1603339) is quoted in UNIQUE K-MERS. Kraken2's column 5 counts DISTINCT MINIMIZERS, which are
a sample of the k-mers, not all of them. Each k-mer of length k spans w = k - l + 1
candidate minimizers of length l, and the standard expected density of selected
minimizers is 2 / (w + 1). For the usual nucleotide settings (k = 35, l = 31) that is
w = 5 and a density of 1/3 — roughly one distinct minimizer per three k-mers. So a
published "1,000 unique k-mers" is about 333 here. Feeding the published number
straight into the minimizer column makes the filter three times harsher than its
authors intended. Check k and l on the first line of your database's inspect.txt.

CONSUMED BY
-----------
The gated report feeds rule bracken (45_abundance.smk) when the gate is on, and from
there the correction flows through kraken_biom into the final otu_table.biom /
otu_table.tsv — so the pipeline's normal output table IS the corrected table, not a
side artefact. The TSV and PNG are terminal outputs: nothing downstream reads them.
"""
import matplotlib

matplotlib.use("Agg")  # render to file only — no display on a compute node
import matplotlib.pyplot as plt

sm = snakemake  # noqa: F821

# Everything comes from rule kmer_evidence in 40_kmer_evidence.smk, which in turn
# reads the shotgun.kmer_evidence block resolved in 00_common.smk.
report_path    = str(sm.input.report)
out_evidence   = str(sm.output.evidence)
out_gated      = str(sm.output.gated)
out_plot       = str(sm.output.plot)

sample         = str(sm.params.sample)
target_rank    = str(sm.params.rank)              # shotgun.bracken.tax_lev — the rank Bracken will report
min_distinct   = int(sm.params.min_distinct)      # distinct-minimizer floor
min_reads      = int(sm.params.min_reads)         # clade-read floor, applied together with the above
protect_frac   = float(sm.params.protect_frac)    # never prune a taxon holding more than this share of classified reads
gate_enabled   = bool(sm.params.enabled)          # labels the log and the figure title; the gated
                                                  # report itself is written either way. What this
                                                  # flag really decides lives in rule bracken.

log = open(sm.log[0], "w")


def logmsg(msg: str) -> None:
    """One line to the log, flushed immediately. Most of the logging happens in the
    summary block at the end of the script; the flush matters for the few warnings
    raised during parsing, which would otherwise sit in the buffer and be lost if a
    later step died."""
    log.write(f"[kmer_evidence] {msg}\n")
    log.flush()


# ───────────────────────── Read the report ─────────────────────────
# A Kraken2 report with --report-minimizer-data has eight tab-separated columns:
#   0 percentage of classified reads in this clade
#   1 reads in this clade   (cumulative: this taxon + everything below it)
#   2 reads assigned directly to this taxon and no deeper
#   3 total minimizers in the read data associated with this clade
#   4 distinct minimizers in the read data associated with this clade (HyperLogLog estimate)
#   5 rank code (U, R, D, K, P, C, O, F, G, S, and the ...1/...2 sub-ranks)
#   6 NCBI taxid
#   7 scientific name, indented two spaces per level of the taxonomy tree
#
# Depth comes from that indentation and nothing else — it is what Bracken uses to
# rebuild the tree, so it is what we use to find parents and subtrees.
COL_CLADE_READS  = 1
COL_DIRECT_READS = 2
COL_TOTAL_MIN    = 3
COL_DISTINCT_MIN = 4
COL_RANK         = 5
COL_TAXID        = 6
COL_NAME         = 7
N_COLS           = 8

rows = []  # one dict per report line, in file order — the order is the tree traversal

with open(report_path) as fh:
    for lineno, line in enumerate(fh, start=1):
        line = line.rstrip("\n")
        if not line.strip():
            continue
        fields = line.split("\t")
        # Fail loud rather than guess. Six columns means kraken2 ran without
        # --report-minimizer-data, in which case there is no evidence to score and
        # silently continuing would produce a meaningless table.
        if len(fields) != N_COLS:
            log.close()
            raise SystemExit(
                f"[kmer_evidence] {report_path}:{lineno}: expected {N_COLS} tab-separated "
                f"columns (kraken2 --report-minimizer-data), got {len(fields)}. "
                f"Line: {line!r}"
            )
        raw_name = fields[COL_NAME]
        indent = len(raw_name) - len(raw_name.lstrip(" "))
        rows.append({
            "fields":       fields,                       # kept verbatim so the gated report can be rewritten in place
            "depth":        indent // 2,                  # two spaces per taxonomy level
            "clade_reads":  int(fields[COL_CLADE_READS]),
            "direct_reads": int(fields[COL_DIRECT_READS]),
            "total_min":    int(fields[COL_TOTAL_MIN]),
            "distinct_min": int(fields[COL_DISTINCT_MIN]),
            "rank":         fields[COL_RANK],
            "taxid":        fields[COL_TAXID],
            "name":         raw_name.strip(),
        })

if not rows:
    log.close()
    raise SystemExit(f"[kmer_evidence] {report_path} is empty — nothing to score.")

# ───────────────────── Link every row to its parent ─────────────────────
# The report is a depth-first walk of the taxonomy, so the parent of a row at depth d
# is simply the most recent row seen at depth d-1. A running stack indexed by depth
# gives that in one pass. Rows at depth 0 ("unclassified" and "root") have no parent.
parent_of = {}          # row index -> parent row index
stack = {}              # depth -> index of the most recent row at that depth

for i, row in enumerate(rows):
    d = row["depth"]
    if d > 0 and (d - 1) in stack:
        parent_of[i] = stack[d - 1]
    stack[d] = i
    # Anything deeper than this row belongs to a branch we have now left.
    for deeper in [k for k in stack if k > d]:
        del stack[deeper]

# Total classified reads — the clade total on the root row. This is the denominator
# for the "how big a share of the sample is this taxon" safety rail below.
total_classified = 0
for row in rows:
    if row["taxid"] == "1":
        total_classified = row["clade_reads"]
        break
if total_classified == 0:
    logmsg(f"WARNING: no root row (taxid 1) or zero classified reads in {report_path}; "
           f"the abundance safety rail cannot fire.")

# ─────────────────── Score every taxon at the target rank ───────────────────
# Candidates are rows whose rank code is EXACTLY the Bracken target rank — 'S', not
# 'S1'. Sub-rank rows (strains) are not scored in their own right; they are carried
# along with whatever species they sit under.
#
# A candidate passes when it clears BOTH floors (this is how aMeta pairs its two
# filters), or when it is too abundant to remove on a heuristic. Order matters only
# for the reason string; the outcome is the same either way.
candidates = [i for i, row in enumerate(rows) if row["rank"] == target_rank]

verdicts = {}   # row index -> (decision, reason)

for i in candidates:
    row = rows[i]
    share = (row["clade_reads"] / total_classified) if total_classified else 0.0

    if row["distinct_min"] >= min_distinct and row["clade_reads"] >= min_reads:
        verdicts[i] = ("keep", "passed both floors")
    elif share > protect_frac:
        # Safety rail against a pathological case: a taxon to which Kraken2 itself
        # directly assigned a large share of the sample is never removed on thin
        # evidence alone. If something that dominant really is spurious, that is a
        # database problem to look at by hand, not something to delete silently.
        #
        # Read `share` carefully: it is this taxon's Kraken2 CLADE reads over ALL
        # classified reads. It is NOT the taxon's share of the final Bracken table,
        # which is typically far larger because Bracken pushes higher-rank reads down
        # onto species. So the rail does not, and cannot, protect a taxon that is big
        # only after redistribution — and on the four test samples it never fired once.
        # That is the intended behaviour rather than a gap: the clearest false positive
        # in the test data (Pseudomonas sp. JS425 — 373 Kraken2 reads on 189 distinct
        # minimizers, which Bracken inflated to 64,961 reads, 3% of the sample) is
        # exactly the kind of taxon a Bracken-based rail would have rescued.
        verdicts[i] = ("keep", f"protected: {share * 100:.2f}% of classified reads")
    elif row["distinct_min"] < min_distinct:
        verdicts[i] = ("prune", f"distinct minimizers {row['distinct_min']} < {min_distinct}")
    else:
        verdicts[i] = ("prune", f"clade reads {row['clade_reads']} < {min_reads}")

pruned = {i for i in candidates if verdicts[i][0] == "prune"}

# ───────────────── Work out which rows leave with each pruned taxon ─────────────────
# Everything deeper-indented that follows a pruned row is part of its subtree and goes
# with it — see point (a) in the module docstring for what happens otherwise.
drop = set()
for i in sorted(pruned):
    drop.add(i)
    base_depth = rows[i]["depth"]
    j = i + 1
    while j < len(rows) and rows[j]["depth"] > base_depth:
        drop.add(j)
        j += 1

# Hand each pruned taxon's clade reads back to its parent's own-node count, so Bracken
# redistributes them among the siblings that survived — see point (b) in the docstring.
# Only the top of each pruned subtree contributes; its children's reads are already
# inside its clade total.
credit = {}     # parent row index -> reads gained
orphaned = 0    # reads from pruned taxa with no parent in the report (should not happen)

for i in sorted(pruned):
    # Climb to the nearest ancestor that is still in the report. A pruned taxon's
    # parent is whatever row sits one level above it, which in practice is usually a
    # 'G1' node ("unclassified <Genus>") rather than the genus itself — on the test
    # data G1 rows took 61% of all credited reads. The climb also means that pointing
    # the gate at some other rank cannot silently credit a row that is itself on its
    # way out, which would lose the reads rather than move them.
    p = parent_of.get(i)
    while p is not None and p in drop:
        p = parent_of.get(p)
    if p is None:
        orphaned += rows[i]["clade_reads"]
        continue
    credit[p] = credit.get(p, 0) + rows[i]["clade_reads"]

if orphaned:
    logmsg(f"WARNING: {orphaned} reads came from pruned taxa with no parent row and "
           f"could not be handed back — they will be missing from the Bracken output.")

# How much of that credit can Bracken actually place? It can only push a node's reads
# down onto taxa still present at the target rank beneath it. If pruning emptied a
# parent completely, its credited reads — AND the reads it already held itself — have
# nowhere to go and Bracken drops them. That is the single biggest reason the final
# table loses reads, so it is worth measuring here rather than discovering it later by
# diffing two Bracken runs. This is a report-only estimate: Bracken additionally needs
# a k-mer distribution entry for the node and applies its own -t threshold, so the real
# figure is slightly worse than this.
live_parents = 0
dead_credit = 0
dead_own = 0
for p, gained in credit.items():
    base_depth = rows[p]["depth"]
    j = p + 1
    has_survivor = False
    while j < len(rows) and rows[j]["depth"] > base_depth:
        if j not in drop and rows[j]["rank"] == target_rank:
            has_survivor = True
            break
        j += 1
    if has_survivor:
        live_parents += 1
    else:
        dead_credit += gained
        dead_own += rows[p]["direct_reads"]

# ───────────────────────── Write the gated report ─────────────────────────
# Column 1 (the percentage) and column 2 (the clade total) are left alone on every
# surviving row: removing a species does not change how many reads its genus's clade
# contains, only how they are labelled inside it. Only the parents' column 3 moves.
with open(out_gated, "w") as fh:
    for i, row in enumerate(rows):
        if i in drop:
            continue
        fields = list(row["fields"])
        if i in credit:
            fields[COL_DIRECT_READS] = str(row["direct_reads"] + credit[i])
        fh.write("\t".join(fields) + "\n")

# ───────────────────────── Write the evidence table ─────────────────────────
# One row per candidate taxon, sorted by clade reads so the organisms that matter most
# to the sample are at the top. Two derived numbers are included because they are what
# people actually eyeball:
#   duplicity            = total / distinct minimizers. The average number of times each
#                          of this taxon's own database minimizers was seen in the
#                          sample. Low means the evidence is spread over much distinct
#                          sequence; high means a small set of positions was hit again
#                          and again. On the MetaFlux test data the middle half of taxa
#                          sit between about 15 and 75, with the clearest false
#                          positives in the hundreds — but the absolute scale moves with
#                          sequencing depth and database, so compare taxa WITHIN a
#                          sample rather than against a remembered number.
#   minimizers_per_read  = distinct minimizers / clade reads. Included because it is the
#                          obvious thing to compute, and because it is a trap: it
#                          anti-correlates strongly with abundance (measured Pearson
#                          r = -0.74 against log reads), so thresholding it removes the
#                          abundant taxa and keeps the rare ones. Diagnostic only.
def safe_ratio(numerator: int, denominator: int) -> str:
    return f"{numerator / denominator:.3f}" if denominator else "NA"


with open(out_evidence, "w") as fh:
    fh.write("\t".join([
        "taxid", "name", "rank", "clade_reads", "direct_reads",
        "total_minimizers", "distinct_minimizers", "duplicity",
        "minimizers_per_read", "pct_of_classified", "decision", "reason",
    ]) + "\n")
    for i in sorted(candidates, key=lambda k: rows[k]["clade_reads"], reverse=True):
        row = rows[i]
        decision, reason = verdicts[i]
        share = (row["clade_reads"] / total_classified * 100) if total_classified else 0.0
        fh.write("\t".join([
            row["taxid"], row["name"], row["rank"],
            str(row["clade_reads"]), str(row["direct_reads"]),
            str(row["total_min"]), str(row["distinct_min"]),
            safe_ratio(row["total_min"], row["distinct_min"]),
            safe_ratio(row["distinct_min"], row["clade_reads"]),
            f"{share:.4f}", decision, reason,
        ]) + "\n")

# ───────────────────────── Plot ─────────────────────────
# Reads against distinct minimizers, both on log scales because both span several
# orders of magnitude across a metagenome. Real organisms sit on a broad diagonal
# band; taxa held up by a handful of shared k-mers sit well below it. The horizontal
# line is the threshold, so what the gate removes is the region beneath it.
kept_x  = [rows[i]["clade_reads"]  for i in candidates if verdicts[i][0] == "keep"]
kept_y  = [rows[i]["distinct_min"] for i in candidates if verdicts[i][0] == "keep"]
cut_x   = [rows[i]["clade_reads"]  for i in candidates if verdicts[i][0] == "prune"]
cut_y   = [rows[i]["distinct_min"] for i in candidates if verdicts[i][0] == "prune"]

fig, ax = plt.subplots(figsize=(7, 5.5))
# Points are drawn semi-transparent: in a deep metagenome hundreds of low-abundance
# taxa overlap near the origin, and solid markers would hide how dense that corner is.
ax.scatter(kept_x, kept_y, s=18, alpha=0.6, color="#2a7f62", label=f"keep (n={len(kept_x)})")
ax.scatter(cut_x,  cut_y,  s=18, alpha=0.6, color="#c1442e", label=f"prune (n={len(cut_x)})")
ax.axhline(min_distinct, color="#444444", linestyle="--", linewidth=1,
           label=f"threshold = {min_distinct}")
if min_reads > 0:
    ax.axvline(min_reads, color="#444444", linestyle=":", linewidth=1,
               label=f"read floor = {min_reads}")
# Log axes need at least one positive value on each axis. A report with no rows at the
# target rank (a very shallow database, or tax_lev pointed at a rank the report does not
# use) leaves both scatters empty, and matplotlib then raises inside tight_layout rather
# than drawing an empty frame. Draw an explicit "nothing to show" panel in that case —
# the run should not die over a figure.
if kept_x or cut_x:
    ax.set_xscale("log")
    ax.set_yscale("log")
else:
    ax.text(0.5, 0.5, f"no taxa at rank {target_rank} in this report",
            ha="center", va="center", transform=ax.transAxes, fontsize=11, color="#666666")
ax.set_xlabel("Reads in clade")
ax.set_ylabel("Distinct minimizers in clade")
ax.set_title(f"{sample} — k-mer evidence per taxon (rank {target_rank})\n"
             f"gate {'ON' if gate_enabled else 'OFF (preview only)'}")
ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.5)
fig.tight_layout()
fig.savefig(out_plot, dpi=150)
plt.close(fig)

# ───────────────────────── Summarise to the log ─────────────────────────
n_keep       = len(candidates) - len(pruned)
reads_pruned = sum(rows[i]["clade_reads"] for i in pruned)
# Reads Kraken2 actually placed AT the gated rank. This, not the total classified
# count, is the population the gate can touch at all — in a metagenome most reads
# stop at genus or above and are never candidates. Quoting the pruned reads against
# the total classified count makes the gate look ~30x more harmless than it is.
reads_at_rank = sum(rows[i]["clade_reads"] for i in candidates)
share_rank    = (reads_pruned / reads_at_rank * 100) if reads_at_rank else 0.0
share_all     = (reads_pruned / total_classified * 100) if total_classified else 0.0
protected     = [i for i in candidates if verdicts[i][1].startswith("protected")]

logmsg(f"sample={sample} rank={target_rank} "
       f"min_distinct_minimizers={min_distinct} min_reads={min_reads} "
       f"never_prune_above_fraction={protect_frac}")
logmsg(f"taxa at rank {target_rank}: {len(candidates)} -> kept {n_keep}, pruned {len(pruned)}")
logmsg(f"rows removed from the report (taxa plus their sub-rank children): {len(drop)}")
logmsg(f"reads under pruned taxa: {reads_pruned} — {share_rank:.2f}% of the "
       f"{reads_at_rank} reads Kraken2 placed at rank {target_rank}, "
       f"{share_all:.2f}% of all {total_classified} classified reads. "
       f"Handed back to {len(credit)} parent taxa.")
logmsg(f"of those {len(credit)} parents, {live_parents} still have a surviving taxon at "
       f"rank {target_rank} beneath them and can actually be redistributed; "
       f"{len(credit) - live_parents} were emptied by the pruning. Bracken cannot place "
       f"reads under an emptied parent, so an estimated {dead_credit} credited reads plus "
       f"{dead_own} reads the parents already held will be dropped from the final table.")
# The counts above are Kraken2's. Bracken has not run yet, and it will have been
# assigning far more reads to some of these species than Kraken2 did, by pushing
# genus-level reads down onto them. So the change to the final table is bigger than
# the percentages above — on the MetaFlux test data, 6-19x bigger. Say so here rather
# than let a reassuring "0.3%" be the last word anyone reads.
logmsg("NOTE: those are Kraken2's own counts. Bracken redistributes higher-rank reads "
       "down onto species, so a pruned taxon may have been holding many more reads in "
       "the final table than it holds here. Compare 03.abundance with the gate off to "
       "see the real effect on your data.")
if protected:
    logmsg(f"{len(protected)} taxa kept by the abundance safety rail despite thin evidence: "
           + ", ".join(rows[i]["name"] for i in protected))
if gate_enabled:
    logmsg("gate is ON — Bracken will read the gated report.")
else:
    logmsg("gate is OFF — Bracken will read the original Kraken2 report. The gated "
           "report above is written anyway so you can diff it and see what enabling "
           "the gate would change.")

log.close()
