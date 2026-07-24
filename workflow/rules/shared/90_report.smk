# Aggregate QC across upstream rule outputs via MultiQC.
# Inputs and scan dirs are mode-dependent: amplicon path adds falco, cutadapt,
# and bowtie2 (PhiX) logs; shotgun path adds fastp, kraken, bracken, BBDuk PhiX
# stats, and (when host_genomes is set) BBMap host-removal stats.

def _multiqc_inputs():
    """Resolve MultiQC inputs based on MODE — only the files actually produced
    by the active workflow are listed here."""
    if MODE == "amplicon":
        out = {
            "falco": expand(
                OUT / "stats" / "falco" / "{sample}_R{r}_{stage}" / "fastqc_data.txt",
                stage=FALCO_STAGES, sample=SAMPLES, r=[1, 2],
            ),
            "cutadapt_json": expand(
                OUT / "stats" / "cutadapt" / "{sample}.{pass_}.cutadapt.json",
                sample=SAMPLES, pass_=["passA_5prime", "passA_3prime"],
            ),
        }
        # bowtie2 PhiX logs exist only when PhiX removal runs.
        if REMOVE_PHIX:
            out["rm_phix_logs"] = expand(LOGS / "rm_phix" / "{sample}_phixFilter.log", sample=SAMPLES)
    else:  # shotgun
        out = {
            "fastp":   expand(OUT / "01.preprocessing" / "{s}_fastp.json", s=SAMPLES),
            "kraken":  expand(OUT / "02.classification" / "{s}_report.txt", s=SAMPLES),
            "bracken": expand(OUT / "03.abundance" / "{s}_report.txt", s=SAMPLES),
        }
        if REMOVE_PHIX:
            out["dephix_stats"] = expand(OUT / "01.preprocessing" / "{s}_dephix_stats.txt", s=SAMPLES)
        if HOST_GENOMES:
            out["dehost_stats"] = expand(OUT / "01.preprocessing" / "{s}_dehost_stats.txt", s=SAMPLES)
    return out


def _multiqc_scan_dirs():
    """Directories MultiQC walks. Kept distinct from _multiqc_inputs because
    MultiQC accepts either explicit files or directory trees on the CLI; we use
    directories for cleaner sample naming and let inputs[] enforce the DAG."""
    if MODE == "amplicon":
        dirs = [
            OUT / "stats" / "falco",
            OUT / "stats" / "cutadapt",
        ]
        if REMOVE_PHIX:
            dirs.append(LOGS / "rm_phix")
        return dirs
    return [
        OUT / "01.preprocessing",
        OUT / "02.classification",
        OUT / "03.abundance",
    ]


rule multiqc:
    input:
        **_multiqc_inputs(),
        mqc_config = workflow.basedir + "/multiqc_config.yaml",
    output:
        report = OUT / "multiqc" / "multiqc_report.html",
        data   = directory(OUT / "multiqc" / "multiqc_data"),
    params:
        outdir     = OUT / "multiqc",
        scan_dirs  = _multiqc_scan_dirs(),
        mqc_config = workflow.basedir + "/multiqc_config.yaml",
        # -d --dirs-depth 1 prefixes every sample name with its immediate parent
        # directory (e.g. '01.preprocessing | SRR123_fastp'), making the regex
        # patterns in multiqc_config.yaml unambiguous and insensitive to whether
        # a module derives its name from the filename or from file content.
        # Amplicon mode does not need this because its modules use filename-based
        # naming consistently and the existing regexes already work without -d.
        dirs_flag  = "-d --dirs-depth 1" if MODE == "shotgun" else "",
    log:
        LOGS / "multiqc.log",
    conda:
        "../../envs/multiqc.yaml"
    threads: lambda wc: threads_for("multiqc")
    resources:
        mem_mb = lambda wc: mem_mb_for("multiqc"),
    shell:
        # Falco output is FastQC-format, so MultiQC reports it under the
        # "FastQC" module — expected behavior. The amplicon side uses regex
        # renames in multiqc_config.yaml to give each stage a readable label.
        """
        mkdir -p {params.outdir}
        multiqc {params.scan_dirs} \
                {params.dirs_flag} \
                --config {params.mqc_config} \
                --outdir {params.outdir} \
                --force \
                > {log} 2>&1
        """
