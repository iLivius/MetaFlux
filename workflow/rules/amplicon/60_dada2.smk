# DADA2 rules: quality plots → per-sample filter → all-sample denoising chain
# → length filter.
#
# dada_quality_plots  : plotQualityProfile on stripped reads (aggregate diagnostic).
# dada_filter         : per-sample filterAndTrim using truncLen from trunclen.json.
# dada_seqtab         : learnErrors → dada → mergePairs → makeSeqtab → removeBimeras.
# dada_length_filter  : ASV length window applied AFTER target extraction (see 70_extract.smk).


rule dada_length_filter:
    input:
        # seqs_pre: always the pre-extraction DADA2 FASTA, used to record the
        # pre-extraction length distribution for diagnostic comparison.
        seqs_pre      = OUT / "5.dada2" / "seqs.fasta",
        # seqs / seqtab_names: the actual filter source — extracted outputs when
        # extraction is enabled, else the raw DADA2 outputs.
        seqs          = (OUT / "5.dada2" / "seqs_extracted.fasta") if EXTRACTION_ENABLED
                        else (OUT / "5.dada2" / "seqs.fasta"),
        seqtab_names  = (OUT / "5.dada2" / "seqtab_extracted_head_names.txt") if EXTRACTION_ENABLED
                        else (OUT / "5.dada2" / "seqtab_head_names.txt"),
        trunclen_json = OUT / "stats" / "trunclen.json",
        # probe_json included only when expected_length=auto
        probe_json    = [PROBE_JSON] if WANTS_PROBE else [],
    output:
        seqs         = OUT / "5.dada2" / "seqs_lenfilt.fasta",
        seqtab_names = OUT / "5.dada2" / "seqtab_lenfilt_head_names.txt",
        seqtab_seqs  = OUT / "5.dada2" / "seqtab_lenfilt_head_seqs.txt",
        stats_json   = OUT / "stats" / "dada2" / "asv_length_stats.json",
        hist_png     = OUT / "stats" / "dada2" / "asv_length_hist.png",
    params:
        mode               = config["amplicon"]["length_filter"]["mode"],
        window_margin      = config["amplicon"]["length_filter"]["window_margin"],
        manual_range       = config["amplicon"]["length_filter"]["range"],
        amp_type           = AMPLICON_TYPE,
        extraction_enabled = EXTRACTION_ENABLED,
    log:
        LOGS / "dada_length_filter.log",
    conda:
        "../../envs/python_utils.yaml"
    threads: lambda wc: threads_for("dada_length_filter")
    resources:
        mem_mb = lambda wc: mem_mb_for("dada_length_filter"),
    script:
        "../../scripts/amplicon/60d_dada_length_filter.py"


rule dada_quality_plots:
    input:
        fnFs = expand(OUT / "3.stripped" / "{sample}_R1.fastq.gz", sample=SAMPLES),
        fnRs = expand(OUT / "3.stripped" / "{sample}_R2.fastq.gz", sample=SAMPLES),
    output:
        plot_r1     = OUT / "stats" / "dada2" / "stripped_read_R1_qual_plot.png",
        plot_r2     = OUT / "stats" / "dada2" / "stripped_read_R2_qual_plot.png",
        plot_r1_pdf = OUT / "stats" / "dada2" / "stripped_read_R1_qual_plot.pdf",
        plot_r2_pdf = OUT / "stats" / "dada2" / "stripped_read_R2_qual_plot.pdf",
    log:
        LOGS / "dada_quality_plots.log",
    conda:
        "../../envs/dada2.yaml"
    threads: lambda wc: threads_for("dada_quality_plots")
    resources:
        mem_mb = lambda wc: mem_mb_for("dada_quality_plots"),
    script:
        "../../scripts/amplicon/60a_dada_quality_plots.R"


rule dada_filter:
    input:
        r1            = OUT / "3.stripped" / "{sample}_R1.fastq.gz",
        r2            = OUT / "3.stripped" / "{sample}_R2.fastq.gz",
        trunclen_json = OUT / "stats" / "trunclen.json",
        # probe_json present only in auto mode; used to derive min_len from q1
        probe_json    = [PROBE_JSON] if WANTS_PROBE else [],
    output:
        r1    = OUT / "4.filtered" / "{sample}_R1_filt.fastq.gz",
        r2    = OUT / "4.filtered" / "{sample}_R2_filt.fastq.gz",
        stats = OUT / "stats" / "dada2" / "{sample}.filter_stats.json",
    params:
        max_ee       = config["amplicon"]["dada2"]["max_ee"],
        trunc_q      = config["amplicon"]["dada2"]["trunc_q"],
        min_len      = config["amplicon"]["dada2"]["filter"]["min_len"],
        max_len      = config["amplicon"]["dada2"]["filter"]["max_len"],
        amp_type     = AMPLICON_TYPE,
        min_len_stat = config["amplicon"]["dada2"]["filter"]["min_len_stat"],
    log:
        LOGS / "dada_filter" / "{sample}.log",
    conda:
        "../../envs/dada2.yaml"
    threads: lambda wc: threads_for("dada_filter")
    resources:
        mem_mb = lambda wc: mem_mb_for("dada_filter"),
    script:
        "../../scripts/amplicon/60b_dada_filter.R"


rule dada_seqtab:
    input:
        r1           = expand(OUT / "4.filtered" / "{sample}_R1_filt.fastq.gz", sample=SAMPLES),
        r2           = expand(OUT / "4.filtered" / "{sample}_R2_filt.fastq.gz", sample=SAMPLES),
        filter_stats = expand(OUT / "stats" / "dada2" / "{sample}.filter_stats.json", sample=SAMPLES),
    output:
        seqtab_seqs     = OUT / "5.dada2" / "seqtab_head_seqs.txt",
        seqtab_names    = OUT / "5.dada2" / "seqtab_head_names.txt",
        seqs_fasta      = OUT / "5.dada2" / "seqs.fasta",
        read_counts     = OUT / "5.dada2" / "read.counts",
        err_plot_r1     = OUT / "stats" / "dada2" / "filtered_read_R1_error_plot.png",
        err_plot_r2     = OUT / "stats" / "dada2" / "filtered_read_R2_error_plot.png",
        err_plot_r1_pdf = OUT / "stats" / "dada2" / "filtered_read_R1_error_plot.pdf",
        err_plot_r2_pdf = OUT / "stats" / "dada2" / "filtered_read_R2_error_plot.pdf",
    params:
        pool                     = config["amplicon"]["dada2"]["pool"],
        learn_errors_nbases      = config["amplicon"]["dada2"]["learn_errors"]["nbases"],
        learn_errors_max_consist = config["amplicon"]["dada2"]["learn_errors"]["max_consist"],
        omega_a                  = config["amplicon"]["dada2"]["dada"]["omega_a"],
        merge_min_overlap        = config["amplicon"]["dada2"]["merge"]["min_overlap"],
        merge_max_mismatch       = config["amplicon"]["dada2"]["merge"]["max_mismatch"],
        merge_just_concatenate   = config["amplicon"]["dada2"]["merge"]["just_concatenate"],
        merge_trim_overhang      = config["amplicon"]["dada2"]["merge"]["trim_overhang"],
        chimera_method           = config["amplicon"]["dada2"]["chimera"]["method"],
        chimera_min_fold_parent  = config["amplicon"]["dada2"]["chimera"]["min_fold_parent"],
        chimera_allow_one_off    = config["amplicon"]["dada2"]["chimera"]["allow_one_off"],
        seed                     = config["amplicon"]["seed"],
    log:
        LOGS / "dada_seqtab.log",
    conda:
        "../../envs/dada2.yaml"
    threads: lambda wc: threads_for("dada_seqtab")
    resources:
        mem_mb = lambda wc: mem_mb_for("dada_seqtab"),
    script:
        "../../scripts/amplicon/60c_dada_seqtab.R"
