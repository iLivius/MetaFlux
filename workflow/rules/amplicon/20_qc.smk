# Per-sample, per-direction, per-stage read QC via falco.
#
# Falco profiles the reads of every sample, both read directions, at four points of
# the amplicon path, the "stages" (FALCO_STAGES in 00_common.smk lists the stages a
# run has; _falco_input below finds the file of each one):
#   raw       the input FASTQ files, as they came in;
#   nophix    after PhiX removal (02.no_phix/), only when
#             amplicon.decontamination.remove_phix is on;
#   stripped  after primer trimming (03.stripped/);
#   filtered  after DADA2's filterAndTrim (04.filtered/, rule dada_filter in
#             60_dada2.smk): the reads the error models are learned from
#             (learn_errors) and then denoised (dada_seqtab).
# Running the same QC at each stage shows where reads were lost and how the quality
# along the reads changed, instead of only the last snapshot.
#
# Who reads the reports: pick_trunclen (50_trunclen.smk) reads the per-base quality
# tables of the stripped stage to choose where to truncate the reads. Nothing reads
# the reports of the other stages except MultiQC, which shows every stage of every
# sample in its section "Falco" (Falco writes FastQC's format, which MultiQC reads
# with its FastQC module; checked with MultiQC 1.35). choose_error_model
# (60_dada2.smk) waits for the reports of the filtered stage without reading them,
# because its stop messages send the user to them. They are normal inputs of that
# step, so deleting a sample's filtered-stage Falco folder (a target of the run) makes
# the error models be learned again, as deleting its stripped-stage folder already does
# through pick_trunclen. Deleting one file inside such a folder changes nothing at
# first; the file is made again when a step that reads it runs for another reason
# (MultiQC reads fastqc_data.txt of every stage), and then the error models are learned
# again too (see the comment at the input falco_filtered in 60_dada2.smk).
#
# What the filtered stage shows: the quality of the reads DADA2 learns its error
# models from. In the MultiQC report its rows give the reads that passed
# filterAndTrim and the mean quality at each position of the read. A mean hides
# binned quality scores (the few values recent instruments write, for example 9, 20
# and 34), because the mean of a few values can fall anywhere between them. Each
# sample's own Falco page (stats/falco/<sample>_R1_filtered/fastqc_report.html) draws
# the quality along the read as box plots, built from the values the reads carry, so
# binned reads show as flat lines on their common values (where Falco pools
# neighbouring positions into one box, it averages them, which can fall between two
# values). On the NextSeq 2000 reads tested, the lines sat at 34 along nearly the
# whole read and at 20 only near the end of R2; 9, carried by fewer than one base in
# ten, never showed.
# Every value, with its bases, is in the error-model report
# (stats/dada2/error_model_quality_values.png, also in the MultiQC section "DADA2
# error model").
#
# What keeps MultiQC from merging the stages of one sample is NOT the output-dir
# name — MultiQC names each amplicon sample from the "Filename" field inside
# fastqc_data.txt (see the shell comment below), and falco is fed a stage-renamed
# symlink so that field carries the stage. The output-dir leaf name only makes the
# Snakemake outputs unique.

def _falco_input(wildcards):
    """The FASTQ file rule falco profiles for one sample, read direction and stage.

    Parameters
    ----------
    wildcards : Snakemake wildcards
        `sample`, `r` (1 or 2) and `stage` (raw, nophix, stripped or filtered).

    Returns
    -------
    str
        raw: the sample's input file (raw_fastq() in 00_common.smk); nophix: its
        PhiX-filtered file in 02.no_phix/; stripped: its primer-trimmed file in
        03.stripped/; filtered: its quality-filtered file in 04.filtered/, the name
        rule dada_filter writes ({sample}_R{r}_filt.fastq.gz), so that the filtered
        stage runs after dada_filter and profiles exactly the reads DADA2 goes on with.
    """
    if wildcards.stage == "raw":
        return raw_fastq(wildcards.sample, int(wildcards.r))
    if wildcards.stage == "nophix":
        return str(OUT / "02.no_phix" / f"{wildcards.sample}_R{wildcards.r}.fastq.gz")
    if wildcards.stage == "stripped":
        return str(OUT / "03.stripped" / f"{wildcards.sample}_R{wildcards.r}.fastq.gz")
    if wildcards.stage == "filtered":
        return str(OUT / "04.filtered" / f"{wildcards.sample}_R{wildcards.r}_filt.fastq.gz")
    raise ValueError(f"unknown falco stage: {wildcards.stage}")


# Runs falco (a fast, drop-in re-implementation of FastQC — same reports, much
# quicker on large FASTQs) on one sample, one read direction, one pipeline stage.
# Input comes from _falco_input above: the raw FASTQ, the PhiX-filtered FASTQ
# (02.no_phix/, only when PhiX removal is on), the primer-trimmed FASTQ
# (03.stripped/) or the quality-filtered FASTQ (04.filtered/, written by
# dada_filter), depending on which {stage} wildcard Snakemake is building. A
# filtered file is never empty: when filterAndTrim keeps no read of a sample,
# dada_filter writes no file and stops the run with the reason.
# Output is the standard FastQC-style trio (fastqc_data.txt, fastqc_report.html,
# summary.txt) in a per-sample/direction/stage folder. The "stripped" stage's
# fastqc_data.txt is not just a report: pick_trunclen (50_trunclen.smk) reads its
# per-base quality table to decide where to truncate reads before DADA2. Every
# stage's output also feeds the MultiQC aggregate report at the end of the run;
# nothing else reads the other stages (choose_error_model only waits for the
# "filtered" reports, see the header; a filtered report made again still makes the
# error models be learned again).
rule falco:
    input:
        fastq = _falco_input,
    output:
        outdir  = directory(OUT / "stats" / "falco" / "{sample}_R{r}_{stage}"),
        data    = OUT / "stats" / "falco" / "{sample}_R{r}_{stage}" / "fastqc_data.txt",
        report  = OUT / "stats" / "falco" / "{sample}_R{r}_{stage}" / "fastqc_report.html",
        summary = OUT / "stats" / "falco" / "{sample}_R{r}_{stage}" / "summary.txt",
    log:
        LOGS / "falco" / "{sample}_R{r}_{stage}.log",
    conda:
        "../../envs/falco.yaml"
    threads: lambda wc: threads_for("falco")
    wildcard_constraints:
        stage = "raw|nophix|stripped|filtered",
        r = "1|2",
    shell:
        # Falco writes the input filename into the "Filename" field of
        # fastqc_data.txt, which MultiQC then uses as the sample name. To
        # disambiguate raw/nophix/stripped (all share the same basename
        # CDRa1_R1.fastq.gz), feed falco a stage-renamed symlink. Result:
        # Filename = "CDRa1_R1_raw.fastq.gz" → MultiQC sample = "CDRa1_R1_raw".
        # The filtered stage's file has a name of its own (CDRa1_R1_filt.fastq.gz),
        # but it goes through the same renamed link, so that its sample name has the
        # form of the others ("CDRa1_R1_filtered"), which multiqc_config.yaml renames
        # to "QC | quality-filtered reads | CDRa1_R1".
        # These comments are not part of the command Snakemake records (only the
        # quoted text below is; checked with Snakemake 9.23.1), so editing them
        # makes no Falco job run again.
        """
        mkdir -p {output.outdir}
        tmpdir=$(mktemp -d -p $(dirname {output.outdir}) .falco_input.XXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        link="$tmpdir/{wildcards.sample}_R{wildcards.r}_{wildcards.stage}.fastq.gz"
        ln -sf "$(realpath {input.fastq})" "$link"
        falco "$link" -o {output.outdir} > {log} 2>&1
        """
