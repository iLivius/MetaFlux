# Per-sample, per-direction, per-stage read QC via falco.
# What keeps MultiQC from deduplicating raw/nophix/stripped is NOT the output-dir
# name — MultiQC names each amplicon sample from the "Filename" field inside
# fastqc_data.txt (see the shell comment below), and falco is fed a stage-renamed
# symlink so that field carries the stage. The output-dir leaf name only makes the
# Snakemake outputs unique.

def _falco_input(wildcards):
    """Resolve falco input based on the stage wildcard."""
    if wildcards.stage == "raw":
        return raw_fastq(wildcards.sample, int(wildcards.r))
    if wildcards.stage == "nophix":
        return str(OUT / "2.no_phix" / f"{wildcards.sample}_R{wildcards.r}.fastq.gz")
    if wildcards.stage == "stripped":
        return str(OUT / "3.stripped" / f"{wildcards.sample}_R{wildcards.r}.fastq.gz")
    raise ValueError(f"unknown falco stage: {wildcards.stage}")


# Runs falco (a fast, drop-in re-implementation of FastQC — same reports, much
# quicker on large FASTQs) on one sample, one read direction, one pipeline stage.
# Input comes from _falco_input above: the raw FASTQ, the PhiX-filtered FASTQ
# (2.no_phix/, only when PhiX removal is on), or the primer-trimmed FASTQ
# (3.stripped/) depending on which {stage} wildcard Snakemake is building.
# Output is the standard FastQC-style trio (fastqc_data.txt, fastqc_report.html,
# summary.txt) in a per-sample/direction/stage folder. The "stripped" stage's
# fastqc_data.txt is not just a report: pick_trunclen (50_trunclen.smk) reads its
# per-base quality table to decide where to truncate reads before DADA2. Every
# stage's output also feeds the MultiQC aggregate report at the end of the run.
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
        stage = "raw|nophix|stripped",
        r = "1|2",
    shell:
        # Falco writes the input filename into the "Filename" field of
        # fastqc_data.txt, which MultiQC then uses as the sample name. To
        # disambiguate raw/nophix/stripped (all share the same basename
        # CDRa1_R1.fastq.gz), feed falco a stage-renamed symlink. Result:
        # Filename = "CDRa1_R1_raw.fastq.gz" → MultiQC sample = "CDRa1_R1_raw".
        """
        mkdir -p {output.outdir}
        tmpdir=$(mktemp -d -p $(dirname {output.outdir}) .falco_input.XXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        link="$tmpdir/{wildcards.sample}_R{wildcards.r}_{wildcards.stage}.fastq.gz"
        ln -sf "$(realpath {input.fastq})" "$link"
        falco "$link" -o {output.outdir} > {log} 2>&1
        """
