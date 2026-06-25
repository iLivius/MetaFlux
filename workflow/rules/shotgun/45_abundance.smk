# Bracken abundance re-estimation + cross-sample OTU table assembly.
#
# bracken      : per-sample re-estimation of read abundance at a target taxonomic
#                rank. Bracken DBs ship pre-computed k-mer distributions for
#                several read lengths (50, 75, 100, 150, 200, 250, 300); we pick
#                the one closest to the *post-trim* mean read length parsed from
#                the fastp JSON. Mis-using a 150 distrib on 100 bp reads biases
#                the species-level redistribution.
# kraken_biom  : assemble per-sample Bracken reports into one OTU/abundance table
#                (BIOM HDF5 + TSV). The awk pass strips the "_report" suffix that
#                kraken-biom otherwise puts on column names.

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
        biom = OUT / "03.abundance" / "otu_table.biom",
        tsv  = OUT / "03.abundance" / "otu_table.tsv",
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
        kraken-biom {input.reports} --fmt hdf5 -o {output.biom} > {log} 2>&1
        tmp=$(mktemp)
        biom convert -i {output.biom} -o $tmp --to-tsv --header-key taxonomy >> {log} 2>&1
        awk 'NR==2 {{gsub("_report", "", $0)}} {{print}}' $tmp > {output.tsv} 2>> {log}
        rm -f $tmp
        """
