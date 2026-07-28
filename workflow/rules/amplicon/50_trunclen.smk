# Auto-pick truncLen by parsing falco's per-base quality table.
# 16S/18S/gyrB/rpoB: first bin where median-across-samples Q1 < q_threshold, then
# enforce the overlap constraint (truncR1 + truncR2 >= amplicon_length + min_overlap)
# via resolve_policy. ITS: the same first-drop cut is computed and recorded, but the
# overlap constraint is skipped — truncLen is overridden to c(0,0) in dada_filter.
# Manual mode bypasses Q analysis, using config-provided truncLen values.

rule pick_trunclen:
    input:
        falco_r1 = expand(
            OUT / "stats" / "falco" / "{sample}_R1_stripped" / "fastqc_data.txt",
            sample=SAMPLES,
        ),
        falco_r2 = expand(
            OUT / "stats" / "falco" / "{sample}_R2_stripped" / "fastqc_data.txt",
            sample=SAMPLES,
        ),
        # probe_json is conditionally included: only when expected_length="auto",
        # so pick_trunclen depends on amplicon_probe, which in turn depends on that
        # marker's probe reference — e.g. fetch_silva_train (16S) or fetch_uchime (ITS).
        # Otherwise [] = no dependency.
        probe_json = [PROBE_JSON] if WANTS_PROBE else [],
    output:
        json = OUT / "stats" / "trunclen.json",
    params:
        amplicon_type   = AMPLICON_TYPE,
        expected_length = config["amplicon"]["expected_length"],
        min_overlap     = config["amplicon"]["min_overlap"],
        mode            = config["amplicon"]["trunc_len"]["mode"],
        q_threshold     = config["amplicon"]["trunc_len"]["q_threshold"],
        q_floor         = config["amplicon"]["trunc_len"]["q_floor"],
        resolve_policy  = config["amplicon"]["trunc_len"]["resolve_policy"],
        manual_r1       = config["amplicon"]["trunc_len"]["manual_r1"],
        manual_r2       = config["amplicon"]["trunc_len"]["manual_r2"],
        probe_stat      = PROBE_LENGTH_STAT,
        min_read_coverage_pct = config["amplicon"]["trunc_len"]["min_read_coverage_pct"],
    log:
        LOGS / "pick_trunclen.log",
    threads: lambda wc: threads_for("pick_trunclen")
    resources:
        mem_mb = lambda wc: mem_mb_for("pick_trunclen"),
    script:
        "../../scripts/amplicon/50a_pick_trunclen.py"
