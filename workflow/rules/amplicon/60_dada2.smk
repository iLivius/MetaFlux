# DADA2 rules: this file turns the primer-trimmed FASTQs in 3.stripped/ (written
# by trim_primers, 30_preprocess.smk) into a table of ASVs — DADA2's term for
# "amplicon sequence variant", an inferred true biological sequence resolved to
# single-nucleotide precision, as opposed to a fuzzy-radius OTU cluster. Order:
#
# dada_quality_plots  : plotQualityProfile on stripped reads (aggregate diagnostic
#                       only — does not feed any decision; pick_trunclen decides
#                       truncLen separately, from falco's tables, in 50_trunclen.smk).
# dada_filter         : per-sample filterAndTrim using truncLen from trunclen.json.
# dada_seqtab         : learnErrors → dada → mergePairs → makeSeqtab → removeBimeras,
#                       run once across ALL samples together so ASVs are called
#                       consistently for the whole run, not sample by sample.
# dada_length_filter  : ASV length window applied AFTER target extraction (see 70_extract.smk).


# Applies a final length window to the ASVs coming out of DADA2 (and, if
# extraction is enabled, out of target_extract). This is a plausibility filter,
# not a taxonomic one: even after denoising and chimera removal, off-target
# amplification, primer dimers, or extraction artifacts can leave a few ASVs at
# lengths that don't make sense for the marker being run. The keep window is
# centred on the amplicon_probe length distribution (or a manual range from
# config — see 60d_dada_length_filter.py for the exact priority order between
# manual range, probe-based auto, and its fallback).
# Reads BOTH seqs_pre (dada_seqtab's raw seqs.fasta, kept only so the stats JSON
# can show the pre- vs post-extraction length distributions side by side) and
# the actual filter source (seqs_extracted.fasta when extraction ran, otherwise
# seqs.fasta again). Writes the length-filtered FASTA/seqtab that assign_taxonomy
# classifies (via _seqs_for_taxonomy / _seqtab_for_taxonomy in 00_common.smk) and
# that aggregate_read_counts folds into the run's read-tracking table.
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


# Plots the aggregate per-base quality profile (both directions, every sample
# together) of the primer-trimmed reads in 3.stripped/. This is a picture for a
# human to look at — it does NOT feed pick_trunclen's truncation decision, which
# is computed separately and automatically from falco's per-sample quality
# tables (50_trunclen.smk). Nothing downstream reads these plots either; they
# are a QC artifact bundled with the run's other outputs, not an input to
# anything else in the pipeline.
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


# Runs DADA2's filterAndTrim on one sample's primer-trimmed reads (3.stripped/).
# truncLen comes from trunclen.json — pick_trunclen's decision — and cuts every
# read to that length, discarding any read shorter than it. ITS is the one
# exception: 60b_dada_filter.R overrides truncLen to c(0, 0) for it, because a
# fixed cut would discard genuinely short fungal amplicons (see 50_trunclen.smk
# for the full reasoning). maxEE/truncQ/min_len/max_len come straight from
# config for every marker, except min_len, which for ITS only can instead be
# derived from the probe JSON when dada2.filter.min_len_stat is set (see the
# script for the trade-off that makes). Writes the filtered FASTQ pair that
# dada_seqtab reads below, plus a small per-sample JSON of reads in/out that
# dada_seqtab folds into the run's read-tracking counts.
rule dada_filter:
    input:
        r1            = OUT / "3.stripped" / "{sample}_R1.fastq.gz",
        r2            = OUT / "3.stripped" / "{sample}_R2.fastq.gz",
        trunclen_json = OUT / "stats" / "trunclen.json",
        # probe_json present only in auto mode; used by dada_filter to derive min_len
        # from q1, but for ITS only — every other marker keeps the config min_len.
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


# The core DADA2 step, run ONCE across every sample's filtered reads together —
# unlike dada_filter above, which runs per sample — because that is what lets
# DADA2 build one shared error model and call ASVs consistently across the
# whole run:
#   learnErrors  → fits a per-direction error model from a subsample of the
#                  filtered reads (learn_errors.nbases in config).
#   dada         → the actual denoising step: infers which true biological
#                  sequence (ASV) each read most likely came from, correcting
#                  for the sequencing errors the model above describes. `pool`
#                  controls whether samples borrow strength from each other
#                  when doing this (false/pseudo/true — see config for the
#                  trade-off).
#   mergePairs   → stitches each sample's denoised R1 and R2 into one
#                  full-length ASV wherever they overlap by at least
#                  merge.min_overlap bp.
#   makeSequenceTable → assembles every sample's merged ASVs into one
#                  sample x ASV count matrix.
#   removeBimeraDenovo → flags and drops chimeric ASVs: sequences that look
#                  like the front half of one real ASV spliced to the back half
#                  of another — a PCR artifact, not a real organism.
# Writes seqs.fasta (the ASV sequences), two seqtab flavours (ASV-ID-keyed and
# sequence-keyed), a per-sample read-tracking table, and the error-model
# diagnostic plots. Consumed next by target_extract (if extraction is enabled,
# 70_extract.smk) or directly by dada_length_filter above.
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
