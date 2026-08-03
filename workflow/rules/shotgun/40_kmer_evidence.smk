# Species-level k-mer evidence — scores how much distinct sequence supports each
# taxon Kraken2 called, and optionally removes the ones that fail before Bracken
# ever sees them.
#
# WHERE THIS SITS
#   kraken2 (35_kraken2.smk)  ->  kmer_evidence (here)  ->  bracken (45_abundance.smk)
#
# It is a step BETWEEN the two, not an alternative to either. Switching the gate on
# does not skip Bracken; it hands Bracken a cleaned report instead of the raw one.
# Because the correction happens upstream, it flows through bracken -> kraken_biom ->
# finalize_otu_table, so the pipeline's ordinary otu_table.biom / otu_table.tsv become
# the corrected tables. There is no separate "filtered" table to reconcile afterwards.
#
# WHY NOT FILTER THE FINISHED TABLE INSTEAD
#   Deleting a cell from the final table throws the reads away. Removing the taxon
#   beforehand hands its reads back to its parent, so Bracken redistributes them among
#   the relatives that do have evidence — which is where they most likely came from,
#   since the reason we distrusted the call is that its k-mers are shared with a
#   neighbour. Same intent, very different arithmetic.
#
# WHAT ALWAYS HAPPENS vs WHAT THE FLAG CONTROLS
#   This rule runs on every shotgun run. It always writes the evidence table, the plot
#   AND the gated report — the last one so you can diff it against the Kraken2 report
#   and see exactly what the gate would remove from your own samples before trusting
#   it. shotgun.kmer_evidence.enabled decides one thing only: whether rule bracken
#   reads the gated report or the original. With the gate off the pipeline behaves
#   exactly as it did before this rule existed.
#
# Outputs and who reads them:
#   02b.evidence/{sample}_species_evidence.tsv -> nothing; a final artefact for you
#   02b.evidence/{sample}_report_gated.txt     -> rule bracken, but only when the
#                                                 gate is on (see 45_abundance.smk)
#   stats/kmer_evidence/{sample}_evidence.png  -> nothing; a final artefact for you
#
# All the parsing, scoring and subtree-pruning logic lives in
# scripts/shotgun/40a_kmer_evidence.py — see that script's docstring for the two
# details that have to be right (prune whole subtrees; hand the reads to the parent)
# and for the k-mer/minimizer unit conversion behind the default threshold.
#
# priority: 3, the same band as kraken2 — this is cheap bookkeeping on a file that
# already exists, and running it promptly keeps bracken unblocked.

rule kmer_evidence:
    input:
        report = OUT / "02.classification" / "{sample}_report.txt",
    output:
        evidence = OUT / "02b.evidence" / "{sample}_species_evidence.tsv",
        gated    = OUT / "02b.evidence" / "{sample}_report_gated.txt",
        plot     = OUT / "stats" / "kmer_evidence" / "{sample}_evidence.png",
    params:
        sample       = "{sample}",
        # Gate at whatever rank Bracken is going to report, so the rows being scored
        # are exactly the rows that become table rows. All resolved in 00_common.smk.
        rank         = TAX_LEV,
        min_distinct = KMER_MIN_DISTINCT,
        min_reads    = KMER_MIN_READS,
        protect_frac = KMER_PROTECT_FRAC,
        enabled      = KMER_GATE_ENABLED,
    threads: lambda wc: threads_for("kmer_evidence")
    resources:
        mem_mb = lambda wc: mem_mb_for("kmer_evidence"),
    conda:
        "../../envs/python_utils.yaml"
    log:
        LOGS / "kmer_evidence" / "{sample}.log"
    priority: 3
    script:
        "../../scripts/shotgun/40a_kmer_evidence.py"
