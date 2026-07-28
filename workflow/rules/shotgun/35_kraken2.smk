# Kraken2 taxonomic classification — the core classification step of the
# shotgun pipeline. Takes the fastp-trimmed reads (trim_adapters's output,
# 30_preprocess.smk) and, for every read, either assigns it to a node in the
# NCBI taxonomy tree or leaves it unclassified, by matching the read's k-mers
# (short fixed-length substrings) against the Kraken2 database (KRAKEN_DB, a
# directory of pre-built k-mer -> taxon lookup tables the user supplies via
# config.references.kraken_db — MetaFlux does not build this database itself).
#
# --confidence (KRAKEN_CONF) is a stringency threshold: raising it demands
# that a bigger share of a read's k-mers agree with the assigned taxon before
# Kraken2 accepts the call, trading raw classified-read count for accuracy
# (see the shotgun.kraken.confidence comment in config.yaml for the benchmark
# literature behind the chosen default). --minimum-hit-groups (KRAKEN_HIT) is
# a related stringency knob, counting distinct groups of consecutive matching
# k-mers rather than individual k-mers. --report-minimizer-data adds two
# extra columns to the report; those extra columns are why
# 90a_aggregate_read_counts.py (80_stats.smk) reads the rank code from column
# 6 rather than the standard report's column 4 — don't drop this flag without
# checking that script too.
#
# Outputs and who reads them downstream:
#   {sample}_report.txt      -> bracken (45_abundance.smk) re-estimates
#                                abundance from it; aggregate_read_counts
#                                and MultiQC (80_stats.smk / 90_report.smk)
#                                also read it, for classified/unclassified
#                                read counts.
#   {sample}_output.txt      -> one line per input read (classified AND
#                                unclassified), recording which taxid, if
#                                any, it was assigned to. extract_taxon_reads
#                                (55_extract_taxa.smk) reads this, alongside
#                                the R1/R2 fastqs below, to pull out every
#                                read under a requested taxon.
#   {sample}_R1/_R2.fastq.gz -> ONLY the reads Kraken2 classified (built from
#                                --classified-out, not --unclassified-out) —
#                                the source sequences extract_taxon_reads
#                                extracts from.
#
# The `--classified-out` template uses # as a placeholder that kraken2 expands
# to _1 / _2 (the leading underscore is included by kraken2 itself), so this
# rule produces {sample}_1.fastq and _2 files — we rename to _R1/_R2 below.
#
# With --memory-mapping the DB is mmap'd and shared across concurrent kraken2
# processes via the OS file cache, so parallel kraken2 jobs can run safely;
# without it each process privately loads the full DB (~hash.k2d size each).
# The per-process mem_mb is computed at parse time in 00_common.smk based on
# shotgun.kraken.memory_mapping and the actual hash.k2d size.
#
# priority: 3, continuing the pipeline-order countdown started in
# 30_preprocess.smk — see that file's header for why.

rule kraken2:
    input:
        r1 = OUT / "01.preprocessing" / "{sample}_trim_R1.fastq.gz",
        r2 = OUT / "01.preprocessing" / "{sample}_trim_R2.fastq.gz",
    output:
        c_r1   = OUT / "02.classification" / "{sample}_R1.fastq.gz",
        c_r2   = OUT / "02.classification" / "{sample}_R2.fastq.gz",
        report = OUT / "02.classification" / "{sample}_report.txt",
        out    = OUT / "02.classification" / "{sample}_output.txt",
    params:
        db             = str(KRAKEN_DB),
        confidence     = KRAKEN_CONF,
        hit_groups     = KRAKEN_HIT,
        # kraken2 v2.1.3 doesn't auto-gzip --classified-out even when the
        # filename ends in .gz — it just writes plain fastq with a misleading
        # name. So we ask for plain .fastq and pigz them ourselves into the
        # final _R1/_R2.fastq.gz outputs.
        classified     = str(OUT / "02.classification" / "{sample}#.fastq"),
        c_r1_native    = str(OUT / "02.classification" / "{sample}_1.fastq"),
        c_r2_native    = str(OUT / "02.classification" / "{sample}_2.fastq"),
        memory_mapping = "--memory-mapping" if KRAKEN_MMAP else "",
    threads: KRAKEN2_THREADS
    resources:
        mem_mb = KRAKEN2_MEM_MB,            # computed in 00_common.smk
    conda:
        "../../envs/kraken2.yaml"
    log:
        LOGS / "kraken2" / "{sample}.log"
    priority: 3
    shell:
        """
        mkdir -p $(dirname {output.report})
        kraken2 --db {params.db} --confidence {params.confidence} \
            {params.memory_mapping} \
            --minimum-hit-groups {params.hit_groups} --use-names \
            --threads {threads} --paired --gzip-compressed \
            --classified-out {params.classified} \
            --report {output.report} --report-minimizer-data \
            --output {output.out} {input.r1} {input.r2} > {log} 2>&1
        pigz -p {threads} -c {params.c_r1_native} > {output.c_r1}
        pigz -p {threads} -c {params.c_r2_native} > {output.c_r2}
        rm -f {params.c_r1_native} {params.c_r2_native}
        """
