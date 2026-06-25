# Per-taxon read extraction via KrakenTools.
#
# wildcards.taxon is the slug (e.g. "bacteria"); we map it back to the real
# scientific name via params.target so we can grep the Bracken report for the
# NCBI taxid at rank D. --include-children means we get everything under the
# requested taxon, not just reads classified exactly at that rank.
#
# Only activated when shotgun.extract_taxa is non-empty (the corresponding
# targets are only added to rule all in that case — see _shotgun_targets()
# in 00_common.smk).

rule extract_taxon_reads:
    input:
        report = OUT / "03.abundance" / "{sample}_report.txt",
        kraken = OUT / "02.classification" / "{sample}_output.txt",
        r1     = OUT / "02.classification" / "{sample}_R1.fastq.gz",
        r2     = OUT / "02.classification" / "{sample}_R2.fastq.gz",
    output:
        r1 = temp(OUT / "04.extracted_reads" / "{sample}_{taxon}_R1.fastq"),
        r2 = temp(OUT / "04.extracted_reads" / "{sample}_{taxon}_R2.fastq"),
    params:
        target = lambda wc: {t.lower().replace(" ", "_"): t for t in EXTRACT_TAXA}[wc.taxon],
    threads: lambda wc: threads_for("extract_taxon_reads")
    resources:
        mem_mb = lambda wc: mem_mb_for("extract_taxon_reads"),
    conda:
        "../../envs/krakentools.yaml"
    log:
        LOGS / "extract_taxon_reads" / "{sample}_{taxon}.log"
    shell:
        r"""
        target="{params.target}"
        # Find the NCBI taxid for the requested taxon name in the bracken report.
        taxid=$(awk -v t="$target" -F '\t' '$6 ~ t {{print $5; exit}}' {input.report})
        echo "{wildcards.sample}: taxon='$target' taxid=$taxid" > {log}
        if [ -n "$taxid" ]; then
            extract_kraken_reads.py -k {input.kraken} \
                -1 {input.r1} -2 {input.r2} \
                -t "$taxid" --include-children --fastq-output \
                -o {output.r1} -o2 {output.r2} \
                -r {input.report} >> {log} 2>&1
        else
            # Taxon absent from this sample — emit empty files to satisfy the DAG.
            echo "no taxid found for '$target'" >> {log}
            : > {output.r1}; : > {output.r2}
        fi
        """


# extract_kraken_reads.py doesn't write .gz natively, so we pigz afterwards.
# Kept in the kraken2 env since pigz is already there.
rule compress_extracted:
    input:
        r1 = OUT / "04.extracted_reads" / "{sample}_{taxon}_R1.fastq",
        r2 = OUT / "04.extracted_reads" / "{sample}_{taxon}_R2.fastq",
    output:
        r1 = OUT / "04.extracted_reads" / "{sample}_{taxon}_R1.fastq.gz",
        r2 = OUT / "04.extracted_reads" / "{sample}_{taxon}_R2.fastq.gz",
    threads: DL_THREADS
    conda:
        "../../envs/kraken2.yaml"
    log:
        LOGS / "compress_extracted" / "{sample}_{taxon}.log"
    shell:
        """
        pigz -p {threads} -c {input.r1} > {output.r1} 2> {log}
        pigz -p {threads} -c {input.r2} > {output.r2} 2>> {log}
        """
