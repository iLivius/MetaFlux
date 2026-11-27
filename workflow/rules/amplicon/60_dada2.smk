# DADA2 rules: this file turns the primer-trimmed FASTQs in 03.stripped/ (written
# by trim_primers, 30_preprocess.smk) into a table of ASVs — DADA2's term for
# "amplicon sequence variant", an inferred true biological sequence resolved to
# single-nucleotide precision, as opposed to a fuzzy-radius OTU cluster. Order:
#
# dada_quality_plots  : plotQualityProfile on stripped reads (aggregate diagnostic
#                       only — does not feed any decision; pick_trunclen decides
#                       truncLen separately, from falco's tables, in 50_trunclen.smk).
# dada_filter         : per-sample filterAndTrim using truncLen from trunclen.json;
#                       also counts the quality values the filtered reads carry.
# choose_error_model  : pools those counts and picks, per read direction, the error
#                       function DADA2 learns with: loess for continuous quality
#                       values (MiSeq), binned for binned ones (NovaSeq, NextSeq),
#                       or the one set in the config. A local rule.
# learn_errors        : learnErrors for R1, then R2, over ALL samples' filtered reads,
#                       with the chosen functions (or, under the experimental
#                       error_model: custom, the user's own): the two error models
#                       (05.dada2/*.rds) and a record of how their learning ended
#                       (stats/dada2/error_model.json).
# error_model_report  : the report on the two models: how closely each follows the
#                       error rates it was learned from (every error function
#                       refitted on the same counts), how the learning ended, the
#                       quality values behind the choice; tables, figures and the
#                       MultiQC section, on disk even when the check below stops the
#                       run.
# check_error_model   : stops the run before denoising when a model did not
#                       converge and its code is not in accept_unconverged;
#                       otherwise writes a one-line verdict. A local rule.
# dada_seqtab         : dada → mergePairs → makeSeqtab → removeBimeras with the saved
#                       models, run once across ALL samples together so ASVs are called
#                       consistently for the whole run, not sample by sample.
# dada_length_filter  : ASV length window applied AFTER target extraction (see 70_extract.smk).
#
# Why no input in this file is wrapped in ancient(). Three inputs here are files a step
# only waits for: the Falco reports of the filtered reads (choose_error_model), the
# error-model report (check_error_model) and the verdict of the check (dada_seqtab).
# ancient() would say "the file must exist, but a newer copy is no reason to run this
# step again", which would spare the denoising when only the report is drawn again. In
# Snakemake 9 that is not safe. Once a job has to run, Snakemake stops checking the
# jobs after it for reasons of their own (a changed script, setting or environment)
# and runs them only because they read that job's output (update_needrun in
# Snakemake's dag.py; checked with Snakemake 9.23.1). An ancient() input cuts that
# chain, so a step further down with a reason of its own, such as a new taxonomy
# script, would be left out of that run without a word and run only on the next one.
# So these are normal inputs, and the price is said where they are declared: a step
# made again makes every step after it run again.


# Applies a final length window to the ASVs coming out of DADA2 (and, if
# extraction is enabled, out of target_extract). This is a plausibility filter,
# not a taxonomic one: even after denoising and chimera removal, off-target
# amplification, primer dimers, or extraction artifacts can leave a few ASVs at
# lengths that don't make sense for the marker being run. The keep window is
# centred on the amplicon_probe length distribution (or a manual range from
# config — see 60h_dada_length_filter.py for the exact priority order between
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
        seqs_pre      = OUT / "05.dada2" / "seqs.fasta",
        # seqs / seqtab_names: the actual filter source — extracted outputs when
        # extraction is enabled, else the raw DADA2 outputs.
        seqs          = (OUT / "05.dada2" / "seqs_extracted.fasta") if EXTRACTION_ENABLED
                        else (OUT / "05.dada2" / "seqs.fasta"),
        seqtab_names  = (OUT / "05.dada2" / "seqtab_extracted_head_names.txt") if EXTRACTION_ENABLED
                        else (OUT / "05.dada2" / "seqtab_head_names.txt"),
        trunclen_json = OUT / "stats" / "trunclen.json",
        # probe_json included only when expected_length=auto
        probe_json    = [PROBE_JSON] if WANTS_PROBE else [],
    output:
        seqs         = OUT / "05.dada2" / "seqs_lenfilt.fasta",
        seqtab_names = OUT / "05.dada2" / "seqtab_lenfilt_head_names.txt",
        seqtab_seqs  = OUT / "05.dada2" / "seqtab_lenfilt_head_seqs.txt",
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
        "../../scripts/amplicon/60h_dada_length_filter.py"


# Plots the aggregate per-base quality profile (both directions, every sample
# together) of the primer-trimmed reads in 03.stripped/. This is a picture for a
# human to look at — it does NOT feed pick_trunclen's truncation decision, which
# is computed separately and automatically from falco's per-sample quality
# tables (50_trunclen.smk). Nothing downstream reads these plots either; they
# are a QC artifact bundled with the run's other outputs, not an input to
# anything else in the pipeline.
rule dada_quality_plots:
    input:
        fnFs = expand(OUT / "03.stripped" / "{sample}_R1.fastq.gz", sample=SAMPLES),
        fnRs = expand(OUT / "03.stripped" / "{sample}_R2.fastq.gz", sample=SAMPLES),
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


# Runs DADA2's filterAndTrim on one sample's primer-trimmed reads (03.stripped/).
# truncLen comes from trunclen.json — pick_trunclen's decision — and cuts every
# read to that length, discarding any read shorter than it. ITS is the one
# exception: 60b_dada_filter.R overrides truncLen to c(0, 0) for it, because a
# fixed cut would discard genuinely short fungal amplicons (see 50_trunclen.smk
# for the full reasoning). maxEE/truncQ/min_len/max_len come straight from
# config for every marker, except min_len, which for ITS only can instead be
# derived from the probe JSON when dada2.filter.min_len_stat is set (see the
# script for the trade-off that makes). Writes the filtered FASTQ pair, which
# learn_errors reads first (the error models are learned from it) and dada_seqtab
# then denoises; rule falco (20_qc.smk, stage "filtered") also profiles each file for
# the MultiQC report. Plus two small per-sample JSON files: filter_stats.json, the
# reads in/out, which dada_seqtab reads for the run's read-tracking table (unchanged
# since v2.6.1), and quality_values.json, the quality values found in the first
# 100,000 filtered reads of each file, which choose_error_model reads. The quality
# values have a file of their own so that an output folder made before it existed
# runs this step again by itself, under any --rerun-triggers: choose_error_model,
# which every run asks for, needs the file, and Snakemake makes a missing input of a
# step that has to run. (A field inside filter_stats.json would leave that to the
# script's date, which Snakemake checks only when --rerun-triggers includes code.)
# When no read pair passes, filterAndTrim writes no files and the script stops with a
# message that names the filtering settings.
rule dada_filter:
    input:
        r1            = OUT / "03.stripped" / "{sample}_R1.fastq.gz",
        r2            = OUT / "03.stripped" / "{sample}_R2.fastq.gz",
        trunclen_json = OUT / "stats" / "trunclen.json",
        # probe_json present only in auto mode; used by dada_filter to derive min_len
        # from q1, but for ITS only — every other marker keeps the config min_len.
        probe_json    = [PROBE_JSON] if WANTS_PROBE else [],
    output:
        r1             = OUT / "04.filtered" / "{sample}_R1_filt.fastq.gz",
        r2             = OUT / "04.filtered" / "{sample}_R2_filt.fastq.gz",
        stats          = OUT / "stats" / "dada2" / "{sample}.filter_stats.json",
        quality_values = OUT / "stats" / "dada2" / "{sample}.quality_values.json",
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


# choose_error_model and check_error_model run on the machine Snakemake itself runs
# on, never as cluster jobs, for two reasons: each is a few lines of standard-library
# Python that should not wait in a cluster queue behind heavy jobs, and each can stop
# the run with a text meant for the user. A local job prints that text on the terminal
# Snakemake runs in, also under a cluster profile, where the output of a cluster job
# goes to the scheduler's log.
localrules: choose_error_model, check_error_model


# Picks the error function DADA2 learns each error model with (rule learn_errors
# below), from the quality values dada_filter counted in every sample's filtered reads
# (stats/dada2/{sample}.quality_values.json). Why it matters: in every round of
# learning, an error function turns the misreads counted at each quality score into
# error rates. DADA2's default, the loess fit, draws a smooth curve and is made for the
# continuous quality values a MiSeq writes (some 30 distinct values); NovaSeq, NextSeq
# and other recent instruments write only 3 to 6 "binned" values, for which DADA2's
# author recommends the binned fit (makeBinnedQualErrfun).
# With amplicon.dada2.learn_errors.error_model: auto (the default) the script decides
# per read direction from the number of distinct values above Q2 in the samples with
# at least 1,000 filtered reads (3-10 values with gaps between them: binned; 20 or
# more: continuous); loess, binned or scam set in the config is used as set, with a warning
# or a stop when it does not suit the reads. It stops the run, before any learning, when the
# reads carry no usable quality information (fewer than 3 values, as in files
# re-encoded by ENA), when one sample does, when samples from different instruments are
# mixed, or when 11-19 values, or 3-10 that include neighbours (as MiSeq reads keep
# under a trunc_q of 28 to 30), leave the kind open; each stop says how to go on. The
# experimental error_model: custom picks the user's own function for both directions;
# the stops about the data itself still apply. See 60c_choose_error_model.py for the
# full rules.
# Writes stats/dada2/error_model_choice.json (the function and the quality values per
# direction: read by learn_errors, and an input of error_model_report, so that the
# report is drawn again when the choice changes) and
# stats/dada2/quality_values.tsv (one row per sample and direction, for people; nothing
# reads it, and the run does not ask for it, so that deleting it costs no new learning:
# see _amplicon_targets() in 00_common.smk). trunc_q is a param only for the messages
# (a raised trunc_q removes the lower quality values); changing it reruns dada_filter,
# and with it this step, anyway. The step also waits for the Falco reports of every
# sample's filtered reads (input falco_filtered), without reading them: its stop
# messages send the user there.
# As with every input, a Falco report made again makes this step run again, and with
# it the learning and every DADA2 step after it (see the comment at the input).
rule choose_error_model:
    input:
        # The quality values of every sample's filtered reads (dada_filter). In an output
        # folder made before these files existed they are missing, so dada_filter runs
        # again there for every sample, whatever --rerun-triggers says.
        quality_values = expand(OUT / "stats" / "dada2" / "{sample}.quality_values.json", sample=SAMPLES),
        # The Falco report of each sample's filtered R1 and R2 reads (rule falco,
        # stage "filtered", 20_qc.smk). The stop messages of this step, and the
        # out-of-range stop of learn_errors after it, point to these reports, so they
        # must exist when this step runs. Without this input, Snakemake can run this
        # quick step before some of those Falco jobs have started, and a stop here then
        # ends the run without them (unless Snakemake runs with --keep-going): the
        # order in which Snakemake starts jobs that are ready at the same time is not
        # fixed.
        # What it costs: this is a normal input (not ancient(), see the note at the top
        # of this file), so deleting a sample's filtered-stage Falco folder
        # (stats/falco/<sample>_R1_filtered/, a target of the run) makes Falco run again
        # for that sample, then this step, the learning of the error models and every
        # DADA2 step after it. Deleting one file inside the folder changes nothing at
        # first: Snakemake looks for the missing outputs of a step that is not a target
        # only when a step after it has to run anyway (update_needrun in Snakemake's
        # dag.py; checked with Snakemake 9.23.1). The file is made again when a step that
        # reads it runs for another reason (MultiQC reads fastqc_data.txt of every stage),
        # and then the learning and every DADA2 step after it run too. The stripped stage
        # behaves the same way through pick_trunclen (50_trunclen.smk). A new Falco
        # version costs nothing more than before: Falco then runs again at every stage,
        # and its stripped-stage reports already lead to pick_trunclen, dada_filter and
        # the learning.
        falco_filtered = expand(
            OUT / "stats" / "falco" / "{sample}_R{r}_filtered" / "fastqc_report.html",
            sample=SAMPLES, r=[1, 2]),
    output:
        choice      = OUT / "stats" / "dada2" / "error_model_choice.json",
        quality_tsv = OUT / "stats" / "dada2" / "quality_values.tsv",
    params:
        # Checked at start-up (00_common.smk): auto, a name of the registry, or custom
        # (the user's own function, experimental).
        error_model_requested = ERROR_MODEL_SETTING,
        trunc_q               = config["amplicon"]["dada2"]["trunc_q"],
    log:
        LOGS / "choose_error_model.log",
    conda:
        "../../envs/python_utils.yaml"
    script:
        "../../scripts/amplicon/60c_choose_error_model.py"


# Learns DADA2's error model of this run, once for the R1 reads and once for the R2
# reads, from every sample's filtered pair in 04.filtered/ (written by dada_filter
# above). The model says, for each quality score, how often each base is read as
# another one; dada_seqtab needs it to tell a real sequence variant from a
# sequencing error of a more abundant one. learnErrors reads whole samples in random
# order until learn_errors.nbases bases are in (the order is drawn with R's random
# numbers, from amplicon.seed), then repeats two steps on those reads: denoise them
# with the current error rates, and estimate the rates again from the result. The
# loop ends in one of three ways, which the script records:
#   - a round gives back a table that already went into a round: DADA2's own stop,
#     which MetaFlux calls "converged" (the table that went into the same round, or
#     into one a few rounds earlier: a cycle);
#   - the flat-rounds rule (learn_errors.max_flat_rounds): the model has stopped
#     improving, and MetaFlux ends the loop early;
#   - the cap of learn_errors.max_consist rounds.
# Each direction is learned with the error function choose_error_model picked
# (error_model_choice.json), looked up in the registry 60_error_functions.R; the binned
# fit is anchored on the quality values of the voting samples' filtered reads. Under
# error_model: custom (EXPERIMENTAL) both directions are learned with the user's own
# function, read from the R file named in amplicon.dada2.learn_errors.
# custom_error_function, and every table it returns is checked before DADA2 uses it.
# Writes the two models as DADA2 returned them (05.dada2/error_model_R1.rds and
# _R2.rds, read by error_model_report and dada_seqtab) and the record of how their
# learning ended (stats/dada2/error_model.json: the function and why, rounds, the
# change per round, the stop, the model code; read by error_model_report and
# check_error_model).
# The step always finishes unless DADA2 or the error function itself fails, or the
# script's own check of the rounds (stop text L7 in 60d_learn_errors.R) finds that its
# record of the rounds does not match DADA2's tables; that check stops the job after
# DADA2 has finished, without writing the models or the record, and the reason is in the
# log. Otherwise the models and the record are on disk whatever the check after it
# decides.
# Conda environment: workflow/envs/dada2.yaml, like every other DADA2 step, except under
# error_model: custom. The models are then learned with the user's own function, which
# may load the R packages mgcv or scam (a GAM fit: a smooth curve of error rate against
# quality score), so the step runs in dada2_custom.yaml, the same environment with those
# two packages added. They are kept out of dada2.yaml because conda-forge has no r-scam
# for Apple Silicon Macs or Linux on ARM: in dada2.yaml it would stop every amplicon run
# there (the comments in dada2_custom.yaml say more). The choice is made when Snakemake
# reads the workflow, from the setting checked at start-up (ERROR_MODEL_SETTING,
# 00_common.smk), so a run that learns with loess or binned never builds
# dada2_custom.yaml; a run that learns with scam or custom does. Snakemake counts a change of environment as new software, but
# switching error_model to or from custom makes this step run again anyway: the choice
# file it reads changes.
# Threads and memory: dada_seqtab's values unless the run config names
# learn_errors under resources (see threads_for_or in 00_common.smk).
rule learn_errors:
    input:
        r1 = expand(OUT / "04.filtered" / "{sample}_R1_filt.fastq.gz", sample=SAMPLES),
        r2 = expand(OUT / "04.filtered" / "{sample}_R2_filt.fastq.gz", sample=SAMPLES),
        # The function and the quality values per direction. It also carries the
        # setting error_model, so a changed setting changes this file and learning
        # runs again; no error_model param is needed.
        choice = OUT / "stats" / "dada2" / "error_model_choice.json",
        # The registry of error functions and DADA2's vendored binned fit, which the
        # script sources. As inputs, a change to either makes learning run again (a
        # script is otherwise judged by its own date only, and these are not this
        # rule's script). See the header of 60_error_functions.R for what that costs.
        error_functions = workflow.basedir + "/scripts/amplicon/60_error_functions.R",
        vendored_binned = workflow.basedir + "/scripts/amplicon/vendor/dada2_makeBinnedQualErrfun.R",
        # The R file of the user's own error function (experimental), only under
        # error_model: custom (checked at start-up: CUSTOM_ERROR_FUNCTION in
        # 00_common.smk). As an input, an edit of the file makes the models be learned
        # again, with everything after them. Under any other setting it is no input
        # here: a function named only for comparison is read by the report alone, so
        # editing it never makes the models be learned again.
        custom_function_file = [CUSTOM_ERROR_FUNCTION["file"]] if ERROR_MODEL_SETTING == "custom" else [],
    output:
        model_r1 = OUT / "05.dada2" / "error_model_R1.rds",
        model_r2 = OUT / "05.dada2" / "error_model_R2.rds",
        record   = OUT / "stats" / "dada2" / "error_model.json",
    params:
        nbases          = config["amplicon"]["dada2"]["learn_errors"]["nbases"],
        # Both checked at start-up (00_common.smk): a whole number >= 1 and >= 0.
        max_consist     = LEARN_ERRORS_MAX_CONSIST,
        max_flat_rounds = LEARN_ERRORS_MAX_FLAT_ROUNDS,
        seed            = config["amplicon"]["seed"],
        # Under error_model: custom, the function's name in that file and the path as
        # written in the config (for the record); "" under any other setting, so that
        # a function named only for comparison never makes the learning run again.
        custom_function_name               = CUSTOM_ERROR_FUNCTION["name"] if ERROR_MODEL_SETTING == "custom" else "",
        custom_function_file_as_configured = CUSTOM_ERROR_FUNCTION["file_as_configured"] if ERROR_MODEL_SETTING == "custom" else "",
    log:
        LOGS / "learn_errors.log",
    conda:
        # dada2_custom.yaml (dada2.yaml plus mgcv and scam) only when the models are
        # learned with the user's own function; see "Conda environment" above the rule.
        "../../envs/dada2_custom.yaml" if ERROR_MODEL_SETTING in ("custom", "scam") else "../../envs/dada2.yaml"
    threads: lambda wc: threads_for_or("learn_errors", threads_for("dada_seqtab"))
    resources:
        mem_mb = lambda wc: mem_mb_for_or("learn_errors", mem_mb_for("dada_seqtab")),
    script:
        "../../scripts/amplicon/60d_learn_errors.R"


# The report on the two error models learn_errors wrote, made in every amplicon run
# (60e_error_model_report.R). For each read direction it rebuilds the transition counts
# the model was fitted to (from the record), refits every error function of the
# registry on them (the smoothing step only, not a new learning: the refit of the
# function used must give back the model exactly), and measures how closely each
# fitted table follows the counted rates, separately at the quality values the reads
# carry and at the averaged qualities of identical reads. It writes:
#   - tables: error_model_fit.tsv (the measures per direction and function),
#     error_model_rates.tsv (counted, nominal and fitted rates), error_model_rounds.tsv
#     (the rounds of learning) and error_model_report.json (the same, plus how the
#     learning ended and the report's own warnings);
#   - figures (PNG and PDF): the fit of each model, the convergence of both, the
#     quality values of the filtered reads, and DADA2's own plotErrors picture of each
#     model (the names it has had since v1.0.0);
#   - the "DADA2 error model" section of the MultiQC report: two tables and copies of
#     the four figures in stats/dada2/multiqc/, which rule multiqc takes as inputs.
# It is a step of its own, ahead of check_error_model, because the report is the
# user's evidence for accepting a model that did not converge: the check takes the
# report JSON as an input, so the report is on disk even when the check stops the
# run. The price of that link: a report drawn again (after a change of this rule's
# script, for example) makes the check, dada_seqtab and every step after it run again
# too, though the models stay the same (see the note on ancient() at the top of this
# file). A figure or table that cannot be made becomes a placeholder and a warning,
# never a failed job (a failed report would block the denoising); the job fails only
# when the record or a model cannot be read. Its files are listed once, in
# ERROR_MODEL_REPORT (00_common.smk), and asked for as targets.
# A function of the user's own (EXPERIMENTAL), named in amplicon.dada2.learn_errors.
# custom_error_function, is refitted and drawn with the others whenever it is named,
# labelled "custom (<name>)": as the model used under error_model: custom, as a
# comparison under any other setting. Its file is then an input of this rule, so an edit
# of it draws the report again; through the link described above, the check,
# dada_seqtab and every step after it then run again too, with the same models.
# Conda environment: dada2.yaml, or dada2_custom.yaml (dada2.yaml plus the R packages
# mgcv and scam; see learn_errors above) whenever a function of the user's own is named,
# because the report then calls that function too: as the model used under
# error_model: custom, as a comparison under any other error_model. Chosen when
# Snakemake reads the workflow, from CUSTOM_ERROR_FUNCTION (00_common.smk). With a
# function named only for comparison, the models are learned in dada2.yaml and refitted
# here in dada2_custom.yaml. On Linux (x86-64) that changes nothing: both environments
# are installed from pin files that hold the same builds of R, dada2 and the maths
# libraries, so the refit of the function used gives back the model bit for bit
# (refit_reproduces_model; checked for loess on the 16S test set). On machines without
# pin files, conda solves the two environments separately and may pick other builds;
# the report may then warn that the refit does not give back the model, while the
# models themselves stay as learned. Naming or removing a function makes this step run
# again anyway, because the function's file is then an input of it or no longer one.
# Threads and memory: 1 thread and the memory of dada_quality_plots, the other DADA2
# step that only draws figures (4000 MB in the template), unless the run config names
# error_model_report under resources (see threads_for_or and mem_mb_for_or in
# 00_common.smk). The report draws figures and refits a few curves in one thread,
# whatever it is given (10 to 15 seconds on the runs checked, and about 1 GB of memory
# at its peak). dada_seqtab's 8 threads and 16 GB would hold eight cores of a
# workstation for that one-thread job and, on a cluster, make it wait for an 8-core,
# 16-GB slot before the check after it could print a stop.
rule error_model_report:
    input:
        record   = OUT / "stats" / "dada2" / "error_model.json",
        model_r1 = OUT / "05.dada2" / "error_model_R1.rds",
        model_r2 = OUT / "05.dada2" / "error_model_R2.rds",
        # The choice of error function (its bases per quality value are drawn) and the
        # registry with DADA2's vendored binned fit (every function in it is refitted):
        # inputs so that the report is made again whenever they change.
        choice          = OUT / "stats" / "dada2" / "error_model_choice.json",
        error_functions = workflow.basedir + "/scripts/amplicon/60_error_functions.R",
        vendored_binned = workflow.basedir + "/scripts/amplicon/vendor/dada2_makeBinnedQualErrfun.R",
        # The R file of the user's own error function (experimental), whenever one is
        # named, whatever error_model says (CUSTOM_ERROR_FUNCTION in 00_common.smk).
        custom_function_file = [CUSTOM_ERROR_FUNCTION["file"]] if CUSTOM_ERROR_FUNCTION else [],
    output:
        **ERROR_MODEL_REPORT,
    params:
        # The function's name in that file, and the path as written in the config (for
        # the report's notes); "" when no function of the user's own is named.
        custom_function_name               = CUSTOM_ERROR_FUNCTION["name"] if CUSTOM_ERROR_FUNCTION else "",
        custom_function_file_as_configured = CUSTOM_ERROR_FUNCTION["file_as_configured"] if CUSTOM_ERROR_FUNCTION else "",
    log:
        LOGS / "error_model_report.log",
    conda:
        # dada2_custom.yaml whenever a function of the user's own is named, whatever
        # error_model says; see "Conda environment" above the rule.
        "../../envs/dada2_custom.yaml" if (CUSTOM_ERROR_FUNCTION or ERROR_MODEL_SETTING == "scam") else "../../envs/dada2.yaml"
    threads: lambda wc: threads_for_or("error_model_report", 1)
    resources:
        mem_mb = lambda wc: mem_mb_for_or("error_model_report", mem_mb_for("dada_quality_plots")),
    script:
        "../../scripts/amplicon/60e_error_model_report.R"


# The check between learning and denoising. Reads the record of learn_errors and,
# for each read direction, decides: converged (DADA2's own stop), accepted (the
# model did not converge, but its code, such as R2:0a55be94, is listed in
# amplicon.dada2.learn_errors.accept_unconverged), or neither. With a model of the
# third kind it stops the run before denoising, and prints, on the terminal and in
# its log, what to check in the reads, where the model's figures are, and the exact
# line that accepts this model. Otherwise it writes the verdict, one line such as
#   R1:0a55be94 converged; R2:7ebe3bed accepted (not converged: flat-rounds rule, round 42)
# which dada_seqtab takes as an input and copies into its log, and repeats the
# warnings of the report, if it has any, as a note.
# Why a step of its own, before dada_seqtab: Snakemake deletes the outputs of a job
# when the job starts, so a stop inside dada_seqtab would first remove the finished
# 05.dada2 tables of an earlier run; a stop here deletes nothing.
# accept_unconverged is read inside the job, from the config (00_common.smk says why),
# so adding or removing a code makes no step run again.
rule check_error_model:
    input:
        record = OUT / "stats" / "dada2" / "error_model.json",
        # The report must exist before the check can stop the run, because it is
        # what the stop text sends the user to; the check also reads its list of
        # warnings, for a note. One file stands for the whole report: Snakemake fails a
        # job that does not write all of its outputs. A normal input (not ancient(),
        # see the note at the top of this file): a report made again makes the check
        # run again, and with it dada_seqtab and every step after it.
        report = ERROR_MODEL_REPORT["report_json"],
    output:
        verdict = OUT / "stats" / "dada2" / "error_model_check.txt",
    params:
        # Only to decide whether the stop text points to stats/orientation/, which
        # exists under orientation: mixed alone.
        orientation = ORIENTATION,
    log:
        LOGS / "check_error_model.log",
    conda:
        "../../envs/python_utils.yaml"
    script:
        "../../scripts/amplicon/60f_check_error_model.py"


# The core DADA2 step, run ONCE across every sample's filtered reads together —
# unlike dada_filter above, which runs per sample — because that is what lets
# DADA2 call ASVs consistently across the whole run. It reads the two error models
# that learn_errors saved, once check_error_model has written its verdict, and runs:
#   dada         → the actual denoising step: infers which true biological
#                  sequence (ASV) each read most likely came from, correcting
#                  for the sequencing errors the error model describes. `pool`
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
# sequence-keyed) and a per-sample read-tracking table. Consumed next by
# target_extract (if extraction is enabled, 70_extract.smk) or directly by
# dada_length_filter above.
rule dada_seqtab:
    input:
        r1           = expand(OUT / "04.filtered" / "{sample}_R1_filt.fastq.gz", sample=SAMPLES),
        r2           = expand(OUT / "04.filtered" / "{sample}_R2_filt.fastq.gz", sample=SAMPLES),
        filter_stats = expand(OUT / "stats" / "dada2" / "{sample}.filter_stats.json", sample=SAMPLES),
        model_r1     = OUT / "05.dada2" / "error_model_R1.rds",
        model_r2     = OUT / "05.dada2" / "error_model_R2.rds",
        # The verdict of check_error_model: it must exist (the check let the run go
        # on); dada_seqtab only copies it into its log. A verdict that does not exist
        # yet, as after a stop followed by an acceptance, makes the check run and then
        # this step. A normal input (not ancient(), see the note at the top of this
        # file): a check that runs again, even on the same models, makes this step and
        # every step after it run again.
        verdict      = OUT / "stats" / "dada2" / "error_model_check.txt",
    output:
        seqtab_seqs     = OUT / "05.dada2" / "seqtab_head_seqs.txt",
        seqtab_names    = OUT / "05.dada2" / "seqtab_head_names.txt",
        seqs_fasta      = OUT / "05.dada2" / "seqs.fasta",
        read_counts     = OUT / "05.dada2" / "read.counts",
    params:
        pool                     = config["amplicon"]["dada2"]["pool"],
        omega_a                  = config["amplicon"]["dada2"]["dada"]["omega_a"],
        merge_min_overlap        = config["amplicon"]["dada2"]["merge"]["min_overlap"],
        merge_max_mismatch       = config["amplicon"]["dada2"]["merge"]["max_mismatch"],
        merge_just_concatenate   = config["amplicon"]["dada2"]["merge"]["just_concatenate"],
        merge_trim_overhang      = config["amplicon"]["dada2"]["merge"]["trim_overhang"],
        chimera_method           = config["amplicon"]["dada2"]["chimera"]["method"],
        chimera_min_fold_parent  = config["amplicon"]["dada2"]["chimera"]["min_fold_parent"],
        chimera_allow_one_off    = config["amplicon"]["dada2"]["chimera"]["allow_one_off"],
    log:
        LOGS / "dada_seqtab.log",
    conda:
        "../../envs/dada2.yaml"
    threads: lambda wc: threads_for("dada_seqtab")
    resources:
        mem_mb = lambda wc: mem_mb_for("dada_seqtab"),
    script:
        "../../scripts/amplicon/60g_dada_seqtab.R"
