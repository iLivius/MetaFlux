# Kraken2 taxonomic classification.
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
        # kraken2 v2.17 doesn't auto-gzip --classified-out even when the
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
