# aggregate_read_counts : joins per-stage read counts (BBDuk, fastp, Kraken2)
# into one wide TSV — one row per sample, one column per pipeline stage.
# Columns: raw → [nophix] → [no_host] → trimmed → classified → unclassified.
# Optional columns are present only when the corresponding decontamination step
# is active (REMOVE_PHIX / HOST_GENOMES flags resolved in 00_common.smk).

rule aggregate_read_counts:
    input:
        fastp_jsons    = expand(OUT / "01.preprocessing" / "{s}_fastp.json",    s=SAMPLES),
        kraken_reports = expand(OUT / "02.classification" / "{s}_report.txt",   s=SAMPLES),
        **( {"phix_stats": expand(OUT / "01.preprocessing" / "{s}_dephix_stats.txt", s=SAMPLES)}
            if REMOVE_PHIX else {} ),
        **( {"host_stats": expand(OUT / "01.preprocessing" / "{s}_dehost_stats.txt", s=SAMPLES)}
            if HOST_GENOMES else {} ),
    output:
        counts = OUT / "stats" / "read_tracking.txt",
    params:
        samples     = SAMPLES,
        remove_phix = REMOVE_PHIX,
        remove_host = bool(HOST_GENOMES),
    log:
        LOGS / "aggregate_read_counts.log",
    conda:
        "../../envs/python_utils.yaml"
    script:
        "../../scripts/shotgun/90a_aggregate_read_counts.py"
