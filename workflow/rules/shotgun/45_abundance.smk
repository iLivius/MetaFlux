# Bracken abundance re-estimation + cross-sample OTU table assembly.
#
# bracken            : per-sample re-estimation of read abundance at a target
#                rank. Bracken DBs ship pre-computed k-mer distributions for
#                several read lengths (50, 75, 100, 150, 200, 250, 300); we pick
#                the one closest to the *post-trim* mean read length parsed from
#                the fastp JSON. Mis-using a 150 distrib on 100 bp reads biases
#                the species-level redistribution.
# kraken_biom        : assemble per-sample Bracken reports into one raw OTU table
#                (BIOM HDF5).
# finalize_otu_table : normalise the taxonomy strings (unify to ';'-no-space and
#                a Genus+species binomial, matching the amplicon tables), strip
#                the "_report" suffix kraken-biom puts on sample names, and apply
#                the optional rank-aware taxon keep/discard filter. Writes the
#                final otu_table.biom + otu_table.tsv, kept consistent because
#                both are rendered from the same filtered BIOM table.

rule bracken:
    input:
        report = OUT / "02.classification" / "{sample}_report.txt",
        json   = OUT / "01.preprocessing" / "{sample}_fastp.json",
    output:
        report = OUT / "03.abundance" / "{sample}_report.txt",
        out    = OUT / "03.abundance" / "{sample}_output.txt",
    params:
        db      = str(KRAKEN_DB),
        tax_lev = TAX_LEV,
        thresh  = BRACKEN_THRESH,
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
        # closest available bracken k-mer length by squared distance
        kmers=(50 75 100 150 200 250 300)
        closest=$(for k in "${{kmers[@]}}"; do echo "$k $(echo "($mean - $k)^2" | bc)"; done | sort -k2n | head -n1 | cut -d' ' -f1)
        echo "{wildcards.sample}: mean=$mean kmer=$closest" > {log}
        bracken -d {params.db} -i {input.report} -r "$closest" \
            -l {params.tax_lev} -t {params.thresh} \
            -o {output.out} -w {output.report} >> {log} 2>&1
        """


rule kraken_biom:
    input:
        reports = expand(OUT / "03.abundance" / "{s}_report.txt", s=SAMPLES),
    output:
        # Raw table straight from kraken-biom: sample ids still carry the
        # "_report" suffix and taxonomy is kraken-biom's "; "-joined form.
        # finalize_otu_table cleans both. Marked temp — only the final table ships.
        raw_biom = temp(OUT / "03.abundance" / "otu_table.raw.biom"),
    threads: lambda wc: threads_for("kraken_biom")
    resources:
        mem_mb = lambda wc: mem_mb_for("kraken_biom"),
    conda:
        "../../envs/kraken-biom.yaml"
    log:
        LOGS / "kraken_biom.log"
    priority: 1
    shell:
        """
        kraken-biom {input.reports} --fmt hdf5 -o {output.raw_biom} > {log} 2>&1
        """


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
