# Aggregate QC across upstream rule outputs via MultiQC — the very last rule
# in the pipeline, for EITHER mode, and shared by both.
#
# This file is included last in the Snakefile, after whichever mode's rule
# files (rules/amplicon/*.smk or rules/shotgun/*.smk) were pulled in — see the
# comment in the Snakefile itself. That ordering matters because the two
# functions below (_multiqc_inputs / _multiqc_scan_dirs) branch on MODE and
# list the exact stats/log files the ACTIVE mode's rules produce (falco,
# cutadapt, bowtie2 for amplicon; fastp, kraken2, bracken, BBDuk/BBMap for
# shotgun) — those rules need to already exist in the DAG for Snakemake to
# know how to build their outputs, which is only guaranteed once this file
# comes after them. Separately, at RUN time, the MultiQC tool itself also
# walks the directories named in scan_dirs looking for anything else it
# recognizes: inputs[] is what tells Snakemake's dependency graph which files
# must exist first; scan_dirs is what tells MultiQC where to look once they do.
#
# Because this file works for either mode, nothing below should assume an
# amplicon-only or shotgun-only file exists outside of the two MODE branches
# in _multiqc_inputs/_multiqc_scan_dirs — that's the only place the two modes
# are told apart.
#
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


# Runs the MultiQC tool once over every stats/log directory listed in
# scan_dirs, and writes one aggregate HTML report (plus its accompanying
# multiqc_data/ directory of parsed tables) — the single-file summary a
# human actually looks at after a run finishes. This is the last rule in the
# DAG for either mode: rule all (all_targets() in 00_common.smk) lists
# multiqc_report.html as a target for both amplicon and shotgun runs, so
# nothing downstream of this rule exists.
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
