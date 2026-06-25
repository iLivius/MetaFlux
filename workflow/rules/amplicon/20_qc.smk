# Per-sample, per-direction, per-stage read QC via falco.
# Output dirs encode the stage in the leaf name so MultiQC reports each stage
# as a distinct sample (raw/nophix/stripped won't deduplicate).

def _falco_input(wildcards):
    """Resolve falco input based on the stage wildcard."""
    if wildcards.stage == "raw":
        return raw_fastq(wildcards.sample, int(wildcards.r))
    if wildcards.stage == "nophix":
        return str(OUT / "2.no_phix" / f"{wildcards.sample}_R{wildcards.r}.fastq.gz")
    if wildcards.stage == "stripped":
        return str(OUT / "3.stripped" / f"{wildcards.sample}_R{wildcards.r}.fastq.gz")
    raise ValueError(f"unknown falco stage: {wildcards.stage}")


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
