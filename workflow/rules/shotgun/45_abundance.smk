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
# kraken_biom        : assemble per-sample Bracken reports into one raw OTU table
#                (BIOM HDF5, the standard Biological Observation Matrix
#                container most microbiome tools read/write), one row per taxon
#                at the rank Bracken re-estimated at (shotgun.bracken.tax_lev).
# finalize_otu_table : normalise the taxonomy strings (unify to ';'-no-space and
#                a Genus+species binomial, matching the amplicon tables), strip
#                the "_report" suffix kraken-biom puts on sample names, and apply
#                the optional rank-aware taxon keep/discard filter. Writes the
#                final otu_table.biom + otu_table.tsv, kept consistent because
#                both are rendered from the same filtered BIOM table.

# Takes this sample's Kraken2 report (per-rank read counts, 35_kraken2.smk)
# and its fastp JSON (post-trim read-length summary, 30_preprocess.smk) as
# input; writes a re-estimated Bracken report plus a plain-text per-taxon
# table. kraken_biom (below) assembles the reports across all samples into
# one table.
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
# kraken_biom into the final otu_table, which is the whole point of putting the
# gate here rather than on the finished table.
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


# Assembles every sample's Bracken report (input.reports, from rule bracken
# above) into ONE BIOM table — one row per taxon, one column per sample. This
# raw table still carries kraken-biom's own sample-naming and taxonomy-string
# conventions; finalize_otu_table (below) is what cleans those up into the
# form the rest of MetaFlux (and its users) expect. Marked temp() because
# only the finalized table is meant to be read further.
#
# WHICH ROWS END UP IN THE TABLE. kraken-biom keeps only taxa whose rank lies
# between --max and --min. At the --min rank it takes the clade count (report
# column 2: reads on the taxon plus everything below it); above --min it takes
# only the reads placed directly on the taxon (column 3). Both flags are set to
# the Bracken rank, so the table holds exactly the taxa Bracken re-estimated,
# each with Bracken's estimate: the same numbers as {sample}_output.txt.
#
# kraken-biom's defaults (--max O --min S) assume a species report. At S they
# are harmless: Bracken writes 0 direct reads on every row above species, so the
# extra ranks add nothing (otu_table.tsv is byte-identical either way). At any
# other rank they failed without an error (MetaFlux v2.4.0 and earlier):
#   G, F, O - the Bracken rank falls in the "direct reads" branch, and Bracken
#             3.1 writes 0 there for every taxon with child rows in the Kraken2
#             report (upstream Bracken issue #315): empty or near-empty table;
#   C, P, D - above --max O, so nothing was kept. Raising only --min would stop
#             kraken-biom with "Max and Min ranks are out of order", hence both.
rule kraken_biom:
    input:
        reports = expand(OUT / "03.abundance" / "{s}_report.txt", s=SAMPLES),
    output:
        # Raw table straight from kraken-biom: sample ids still carry the
        # "_report" suffix and taxonomy is kraken-biom's "; "-joined form.
        # finalize_otu_table cleans both. Marked temp — only the final table ships.
        raw_biom = temp(OUT / "03.abundance" / "otu_table.raw.biom"),
    params:
        # The same rank rule bracken was given (shotgun.bracken.tax_lev, one
        # capital letter, checked when the workflow is parsed in 00_common.smk).
        tax_lev = TAX_LEV,
    threads: lambda wc: threads_for("kraken_biom")
    resources:
        mem_mb = lambda wc: mem_mb_for("kraken_biom"),
    conda:
        "../../envs/kraken-biom.yaml"
    log:
        LOGS / "kraken_biom.log"
    priority: 1
    shell:
        # Kept on one line: in this plain (not r"...") string Python would remove a
        # backslash line break before the shell ever saw it.
        """
        kraken-biom {input.reports} --max {params.tax_lev} --min {params.tax_lev} --fmt hdf5 -o {output.raw_biom} > {log} 2>&1
        """


# Cleans up kraken_biom's raw table and writes the pipeline's final OTU
# table — both otu_table.biom and otu_table.tsv, rendered from the same
# filtered in-memory table so the two never disagree. This is the shotgun
# equivalent of the amplicon pipeline's asv_table.txt / taxon_seq_table.txt:
# the last file a user actually opens. Sample-id cleanup, taxonomy-string
# normalization (kraken-biom's bare species epithet -> a Genus+species
# binomial where that's unambiguous), and the optional taxon keep/discard
# filter (shotgun.taxonomy_filter in the config) all happen in
# 45a_finalize_otu_table.py — see that script's own docstring for the exact
# rules, including a documented limitation around viral strain names.
rule finalize_otu_table:
    input:
        raw_biom = OUT / "03.abundance" / "otu_table.raw.biom",
    output:
        biom = OUT / "03.abundance" / "otu_table.biom",
        tsv  = OUT / "03.abundance" / "otu_table.tsv",
    params:
        filter_enabled = OTU_FILTER_ENABLED,
        filter_keep    = OTU_FILTER_KEEP,
        filter_discard = OTU_FILTER_DISCARD,
    threads: lambda wc: threads_for("finalize_otu_table")
    resources:
        mem_mb = lambda wc: mem_mb_for("finalize_otu_table"),
    conda:
        "../../envs/kraken-biom.yaml"
    log:
        LOGS / "finalize_otu_table.log"
    priority: 1
    script:
        "../../scripts/shotgun/45a_finalize_otu_table.py"
