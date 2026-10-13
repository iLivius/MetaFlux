# Bracken abundance re-estimation + cross-sample OTU table assembly. This is
# the last stretch of the shotgun pipeline: it takes Kraken2's per-read calls
# (35_kraken2.smk) and turns them into the one table of taxon x sample
# abundances (otu_table.tsv) a user actually opens for analysis.
#
# bracken            : per-sample re-estimation of read abundance at a target
#                rank. Kraken2 classifies many reads at internal nodes ABOVE
#                that rank (e.g. only to Genus, because a read wasn't
#                distinctive enough to call to Species) — Bracken uses a
#                Bayesian model, built from each Kraken2 database's own
#                pre-computed k-mer distributions, to redistribute those
#                reads down onto the target rank, giving a more complete
#                (and more accurate) count at that rank than Kraken2's raw
#                report alone. Bracken DBs ship distributions for several
#                read lengths (50, 75, 100, 150, 200, 250, 300); we pick the
#                one closest to the *post-trim* mean read length parsed from
#                the fastp JSON. Mis-using a 150 distrib on 100 bp reads biases
#                the species-level redistribution.
# finalize_otu_table : build the final OTU table straight from the per-sample
#                Bracken reports: one row per taxon at the rank Bracken
#                re-estimated at (shotgun.bracken.tax_lev), one column per
#                sample, and a lineage of eight fixed slots d__ k__ p__ c__ o__
#                f__ g__ s__ taken from the reports' own rank codes (a virus
#                carries its realm in d__). Applies the optional rank-aware taxon
#                keep/discard filter, then writes otu_table.biom (BIOM HDF5, the
#                standard Biological Observation Matrix container most microbiome
#                tools read/write) and otu_table.tsv from the same numbers.

# Takes this sample's Kraken2 report (per-rank read counts, 35_kraken2.smk)
# and its fastp JSON (post-trim read-length summary, 30_preprocess.smk) as
# input; writes a re-estimated Bracken report plus a plain-text per-taxon
# table. finalize_otu_table (below) builds one table from the reports of all
# samples.
#
# WHICH REPORT THIS READS. One of two files, decided once when the workflow is
# parsed by shotgun.kmer_evidence.enabled (KMER_GATE_ENABLED, 00_common.smk):
#   gate off (default) — Kraken2's own report, exactly as before the k-mer
#                        evidence rule existed. Nothing changes.
#   gate on            — the gated report from rule kmer_evidence
#                        (40_kmer_evidence.smk): the same report with taxa that
#                        failed the k-mer evidence check removed and their reads
#                        handed back to their parent, for Bracken to redistribute
#                        among the relatives that passed. Mostly it can: where a
#                        parent still has a surviving taxon at the target rank the
#                        reads move, but where pruning emptied the parent completely
#                        Bracken has no destination and drops them. The per-sample
#                        kmer_evidence log estimates how many reads that is before
#                        Bracken runs.
# Bracken itself runs either way — the gate filters Bracken's input, it does not
# replace or bypass Bracken. The corrected counts therefore propagate through
# finalize_otu_table into the final otu_table, which is the whole point of putting
# the gate here rather than on the finished table.
#
# priority: 2, continuing the pipeline-order countdown started in
# 30_preprocess.smk — see that file's header for why.
rule bracken:
    input:
        report = (OUT / "02b.evidence" / "{sample}_report_gated.txt") if KMER_GATE_ENABLED
                 else (OUT / "02.classification" / "{sample}_report.txt"),
        json   = OUT / "01.preprocessing" / "{sample}_fastp.json",
    output:
        report = OUT / "03.abundance" / "{sample}_report.txt",
        out    = OUT / "03.abundance" / "{sample}_output.txt",
    params:
        db      = str(KRAKEN_DB),
        tax_lev = TAX_LEV,
        thresh  = BRACKEN_THRESH,     # a number, or the literal word "auto"
        alpha   = BRACKEN_ALPHA,      # used only when thresh == "auto"
        min_t   = BRACKEN_MIN_T,      # floor for the auto value
    threads: lambda wc: threads_for("bracken")
    resources:
        mem_mb = lambda wc: mem_mb_for("bracken"),
    conda:
        "../../envs/bracken.yaml"
    log:
        LOGS / "bracken" / "{sample}.log"
    priority: 2
    shell:
        r"""
        # mean read length (R1+R2)/2 after filtering, from fastp json
        mean=$(jq -r '(.summary.after_filtering.read1_mean_length + .summary.after_filtering.read2_mean_length) / 2 | floor' {input.json})
        # Closest available Bracken read-length database (the kmer_distrib files are
        # built for these read lengths) by squared distance to the mean read length.
        # Computed with jq, which the bracken env provides, rather than with the
        # host's awk (or bc), so the result never depends on what the machine happens
        # to have installed. min_by keeps the FIRST of two equally distant lengths
        # (mean 125 -> 100), exactly as the earlier awk loop did.
        closest=$(jq -n --argjson mean "$mean" \
            '[50, 75, 100, 150, 200, 250, 300] | min_by((. - $mean) * (. - $mean))')
        echo "{wildcards.sample}: mean=$mean kmer=$closest" > {log}

        # Minimum-read threshold. With a fixed number in the config it is used as given.
        # With "auto" it is derived here, at run time, because it depends on how many
        # read pairs THIS sample actually got classified — a number that does not exist
        # until Kraken2 has finished, i.e. long after the workflow was parsed.
        #
        # The classified count is the clade total on the report's root row (taxid 1).
        # Columns are counted from the RIGHT so this works on both report layouts: the
        # standard 6-column one and the 8-column one produced by
        # --report-minimizer-data. From the end, the fields are always
        # ... rank_code, taxid, name, so the taxid is the second-to-last field, and
        # the clade total is the second field in both layouts.
        thresh="{params.thresh}"
        if [ "$thresh" = "auto" ]; then
            # jq reads the report as raw text lines (-R), splits each on tabs and
            # prints column 2 of the first line whose taxid column is exactly "1".
            classified=$(jq -nrR \
                'first(inputs | split("\t") | select(.[-2] == "1") | .[1]) // empty' \
                {input.report})
            if [ -z "$classified" ] || [ "$classified" -le 0 ] 2>/dev/null; then
                echo "[bracken] WARNING: could not read a classified-read count from" \
                     "{input.report}; falling back to the floor {params.min_t}" >> {log}
                thresh={params.min_t}
            else
                # threshold = max(min_t, alpha * classified), rounded to a whole read
                # (x.5 rounds up). jq, from the bracken env, does the arithmetic in
                # the same double precision the earlier host-awk call used.
                thresh=$(jq -n --argjson c "$classified" --argjson a {params.alpha} \
                    --argjson m {params.min_t} '[($a * $c + 0.5 | floor), $m] | max')
                echo "[bracken] auto threshold: {params.alpha} x $classified classified" \
                     "read pairs -> -t $thresh" >> {log}
            fi
        else
            echo "[bracken] fixed threshold from config: -t $thresh" >> {log}
        fi

        bracken -d {params.db} -i {input.report} -r "$closest" \
            -l {params.tax_lev} -t "$thresh" \
            -o {output.out} -w {output.report} >> {log} 2>&1
        """


# Builds the pipeline's final shotgun table, otu_table.tsv and otu_table.biom, straight
# from every sample's Bracken report (rule bracken above). This is the shotgun
# equivalent of the amplicon asv_table.txt / taxon_seq_table.txt: the last file a user
# opens. The docstring of 45a_finalize_otu_table.py explains every rule; in short:
#   rows     one per taxon at the Bracken rank (shotgun.bracken.tax_lev), one column
#            per sample, counts = Bracken's re-estimated reads;
#   lineage  eight fixed slots, d__ k__ p__ c__ o__ f__ g__ s__, taken from the rank
#            codes of the reports themselves (D K P C O F G S), so every Kraken2 table
#            has the same eight slots. On a GTDB database, which has no kingdom rank,
#            k__ is empty. A virus has no domain: its d__ slot holds its realm, the top
#            rank of virus taxonomy, or stays empty when the virus has none;
#   species  the whole name the database gives it;
#   filter   the optional keep/discard filter (shotgun.taxonomy_filter);
#   checks   the job stops when the D or K codes of a report disagree with the
#            database (reports of a Kraken2 older than 2.17, or of another database),
#            when one taxon has two lineages, and when the filter removes every taxon
#            of a table that had some (an empty table would look like a finished run).
#
# Up to MetaFlux v2.4.0 a separate rule ran kraken-biom here, and this rule cleaned its
# output. kraken-biom dropped the kingdom rows, wrote the domain as k__, cut the genus
# off species names, and, depending on which report it read first, could fill a slot
# with a name from a neighbouring branch of the tree. Reading the reports directly
# fixes all four, with the same counts: Bracken's estimates at the Bracken rank, as
# kraken-biom took them.
#
# Inputs: the Bracken report of every sample (03.abundance/{sample}_report.txt), in the
# order of SAMPLES. Kraken2 and Bracken are untouched by this rule: changing it re-runs
# only this step.
# Params: the sample names (the column order), the Bracken rank, what the Kraken2
# database ranks as a domain, a kingdom or a virus realm (read from its taxo.k2d when
# the workflow starts, KRAKEN2_RANK_SETS in 00_common.smk), and the
# shotgun.taxonomy_filter settings.
# Outputs: otu_table.biom and otu_table.tsv, final targets in _shotgun_targets()
# (00_common.smk); nothing else in the pipeline reads them.
# Runs in the biom env (envs/biom.yaml), which it shares with sylph_tables
# (38_sylph.smk): python, numpy and biom-format, nothing else.
rule finalize_otu_table:
    input:
        reports = expand(OUT / "03.abundance" / "{s}_report.txt", s=SAMPLES),
    output:
        biom = OUT / "03.abundance" / "otu_table.biom",
        tsv  = OUT / "03.abundance" / "otu_table.tsv",
    params:
        samples        = SAMPLES,      # the column order, the same as input.reports
        tax_lev        = TAX_LEV,      # the Bracken rank: the rank of every row
        # {taxid: name} of what the Kraken2 database ranks as a domain, a kingdom and a
        # virus realm, and of Viruses itself when the database ranks it superkingdom
        # (built from NCBI taxonomy before 2025); empty dicts are possible.
        domains        = KRAKEN2_RANK_SETS["domain"],
        kingdoms       = KRAKEN2_RANK_SETS["kingdom"],
        realms         = KRAKEN2_RANK_SETS["realm"],
        virus_root     = KRAKEN2_RANK_SETS["virus_root"],
        filter_enabled = OTU_FILTER_ENABLED,
        filter_keep    = OTU_FILTER_KEEP,
        filter_discard = OTU_FILTER_DISCARD,
    threads: lambda wc: threads_for("finalize_otu_table")
    resources:
        mem_mb = lambda wc: mem_mb_for("finalize_otu_table"),
    conda:
        "../../envs/biom.yaml"
    log:
        LOGS / "finalize_otu_table.log"
    priority: 1
    script:
        "../../scripts/shotgun/45a_finalize_otu_table.py"
