# Shotgun preprocessing: optional SRA fetch (BioProject mode), samples manifest,
# PhiX removal (BBDuk), host removal (BBDuk against concatenated host reference),
# and adapter/quality trimming (fastp).
#
# Stage chain (each step optional): raw → (phix) → (host) → fastp. The
# fastp_upstream() / host_upstream() helpers in 00_common.smk pick the right
# upstream stage so a consumer rule doesn't have to know what ran before it.

import shutil
import subprocess
from urllib.parse import urlparse


# ──────────────────────── SRA fetch (BioProject mode) ──────────
# Only triggered when input.bioproject is set. In local-FASTQ mode, the
# dependency on these output paths is satisfied by user-provided files.
# NOTE: prefetch puts the .sra under {tmp}/{accession}/{accession}.sra,
# hence the nested path in the fasterq-dump call below.
rule download_sra:
    output:
        r1 = FASTQ_DIR / "{sample}_R1.fastq.gz",
        r2 = FASTQ_DIR / "{sample}_R2.fastq.gz",
    params:
        tmp = lambda wc: str(FASTQ_DIR / f".tmp_{wc.sample}"),
    threads: lambda wc: threads_for("download_sra")
    conda:
        "../../envs/sra-tools.yaml"
    log:
        LOGS / "download_sra" / "{sample}.log"
    priority: 10
    shell:
        r"""
        mkdir -p {params.tmp}
        prefetch --output-directory {params.tmp} {wildcards.sample} > {log} 2>&1
        fasterq-dump --split-files --skip-technical --threads {threads} \
            --temp {params.tmp} --outdir {params.tmp} \
            {params.tmp}/{wildcards.sample}/{wildcards.sample}.sra >> {log} 2>&1
        pigz -p {threads} -c {params.tmp}/{wildcards.sample}_1.fastq > {output.r1} 2>> {log}
        pigz -p {threads} -c {params.tmp}/{wildcards.sample}_2.fastq > {output.r2} 2>> {log}
        rm -rf {params.tmp}
        """


# ──────────────────────── Samples manifest ─────────────────────
# Small TSV listing every sample's source + FASTQ paths. Useful for
# downstream R/Python scripts that don't want to re-parse the config or
# re-glob fastq_dir.
rule samples_manifest:
    output:
        manifest = OUT / "samples.tsv",
    params:
        samples    = SAMPLES,
        mode       = INPUT_MODE,
        bioproject = BIOPROJECT or "",
        fastq_dir  = str(FASTQ_DIR),
        naming     = FASTQ_NAMING,
        extn       = EXTN,
    run:
        m2 = "_R2" if params.naming == "_R1" else "_2"
        with open(output.manifest, "w") as f:
            f.write("sample\tsource\tbioproject\tr1\tr2\n")
            for s in params.samples:
                r1 = f"{params.fastq_dir}/{s}{params.naming}.{params.extn}"
                r2 = f"{params.fastq_dir}/{s}{m2}.{params.extn}"
                f.write(f"{s}\t{params.mode}\t{params.bioproject}\t{r1}\t{r2}\n")


# ──────────────────────── PhiX removal (BBDuk) ─────────────────
# BBDuk ships phix174_ill.ref.fa.gz inside the bbmap conda package, so
# `ref=phix` works without any download/index-build step.
# k=31 hdist=1 is the BBTools-recommended setting for PhiX (see their
# decontamination guide). stats= produces a small text file that MultiQC
# picks up automatically.
rule decontam_phix:
    input:
        r1 = lambda wc: raw_r1(wc.sample),
        r2 = lambda wc: raw_r2(wc.sample),
    output:
        r1    = temp(OUT / "01.preprocessing" / "{sample}_dephix_R1.fastq.gz"),
        r2    = temp(OUT / "01.preprocessing" / "{sample}_dephix_R2.fastq.gz"),
        stats = OUT / "01.preprocessing" / "{sample}_dephix_stats.txt",
    # PhiX is a 5 kb reference — BBDuk scales poorly past ~4–6 worker threads,
    # so we cap here and let Snakemake run multiple samples in parallel instead.
    threads: lambda wc: threads_for("decontam_phix")
    resources:
        mem_mb = lambda wc: mem_mb_for("decontam_phix"),
    conda:
        "../../envs/bbtools.yaml"
    log:
        LOGS / "decontam_phix" / "{sample}.log"
    priority: 6
    shell:
        """
        bbduk.sh ref=phix in1={input.r1} in2={input.r2} \
            out1={output.r1} out2={output.r2} stats={output.stats} \
            k=31 hdist=1 threads={threads} overwrite=t > {log} 2>&1
        """


# ──────────────────────── Host reference assembly ──────────────
# Resolves shotgun.decontamination.host_genomes into one concatenated
# host_combined.fna.gz. Download URLs, copy local paths, gzip plain FASTAs,
# then binary-concat (gzipped files concatenate to a valid gzip stream, so
# we avoid re-gzipping a 3 Gb human genome).
rule fetch_host_refs:
    output:
        combined = temp(OUT / "01.preprocessing" / "refs" / "host_combined.fna.gz"),
    params:
        sources = HOST_GENOMES,
    log:
        LOGS / "fetch_host_refs.log"
    priority: 8
    run:
        outdir = os.path.dirname(output.combined)
        os.makedirs(outdir, exist_ok=True)
        downloaded = []
        with open(log[0], "w") as lh:
            for i, src in enumerate(params.sources):
                dest = os.path.join(outdir, f"host_{i:02d}.fna.gz")
                if urlparse(str(src)).scheme in ("http", "https", "ftp"):
                    lh.write(f"download {src}\n"); lh.flush()
                    urllib.request.urlretrieve(src, dest)
                else:
                    if not os.path.exists(src):
                        raise FileNotFoundError(f"host reference missing: {src}")
                    if str(src).endswith((".gz", ".bgz")):
                        shutil.copy(src, dest)
                    else:
                        subprocess.check_call(f"gzip -c {src} > {dest}", shell=True)
                downloaded.append(dest)
            # Concat (binary) — gzip streams are catable as-is.
            with open(output.combined, "wb") as out:
                for d in downloaded:
                    with open(d, "rb") as f:
                        shutil.copyfileobj(f, out)
            for d in downloaded:
                os.remove(d)


# ──────────────────────── Host index (BBMap) ──────────────────
# Build the BBMap index ONCE from the concatenated host reference, then reuse it
# for every sample (far cheaper than rebuilding per sample). build=1 writes the
# index under {path}/ref/. Sized for a human-masked reference (~24-28 Gb).
rule build_host_index:
    input:
        ref = OUT / "01.preprocessing" / "refs" / "host_combined.fna.gz",
    output:
        idx = directory(OUT / "01.preprocessing" / "refs" / "bbmap_index" / "ref"),
    params:
        path = lambda wc: str(OUT / "01.preprocessing" / "refs" / "bbmap_index"),
    threads: lambda wc: threads_for("build_host_index")
    resources:
        mem_mb = lambda wc: mem_mb_for("build_host_index"),
    conda:
        "../../envs/bbtools.yaml"
    log:
        LOGS / "build_host_index.log"
    priority: 7
    shell:
        r"""
        xmx=$(( {resources.mem_mb} * 85 / 100 ))
        bbmap.sh ref={input.ref} path={params.path} build=1 \
            threads={threads} -Xmx${{xmx}}m > {log} 2>&1
        """


# ──────────────────────── Host removal (BBMap) ─────────────────
# Alignment-based host removal — the BBTools-authors' recommended approach for
# divergent host reads (cf. removehuman.sh). Reads that do NOT map to the host
# index (outu) are kept; mapped reads are dropped. minid is the identity cutoff.
# The surviving-pair count is written to stats for read-tracking (BBMap, unlike
# BBDuk, emits no #Total/#Matched table, so we count the kept R1 records).
rule decontam_host:
    input:
        unpack(host_upstream),
        idx = rules.build_host_index.output.idx,
    output:
        r1    = temp(OUT / "01.preprocessing" / "{sample}_dehost_R1.fastq.gz"),
        r2    = temp(OUT / "01.preprocessing" / "{sample}_dehost_R2.fastq.gz"),
        stats = OUT / "01.preprocessing" / "{sample}_dehost_stats.txt",
    params:
        path   = lambda wc: str(OUT / "01.preprocessing" / "refs" / "bbmap_index"),
        min_id = HOST_MIN_ID,
    threads: lambda wc: threads_for("decontam_host")
    resources:
        mem_mb = lambda wc: mem_mb_for("decontam_host"),
    conda:
        "../../envs/bbtools.yaml"
    log:
        LOGS / "decontam_host" / "{sample}.log"
    priority: 5
    shell:
        r"""
        xmx=$(( {resources.mem_mb} * 85 / 100 ))
        bbmap.sh path={params.path} in1={input.r1} in2={input.r2} \
            outu1={output.r1} outu2={output.r2} \
            minid={params.min_id} maxindel=3 bwr=0.16 bw=12 \
            quickmatch fast minhits=2 qtrim=rl trimq=10 untrim \
            threads={threads} -Xmx${{xmx}}m > {log} 2>&1
        # Surviving (non-host) read PAIRS = R1 records in the kept output.
        kept=$(( $(zcat {output.r1} | wc -l) / 4 ))
        printf 'surviving_pairs\t%s\n' "$kept" > {output.stats}
        """


# ──────────────────────── Adapter + quality trim (fastp) ───────
# --cut_front / --cut_right is sliding-window trimming from both ends;
# --length_required drops short reads after trimming. fastp emits .fastq.gz
# natively when the output filename ends in .gz.
rule trim_adapters:
    input:
        unpack(fastp_upstream),
    output:
        r1   = temp(OUT / "01.preprocessing" / "{sample}_trim_R1.fastq.gz"),
        r2   = temp(OUT / "01.preprocessing" / "{sample}_trim_R2.fastq.gz"),
        html = OUT / "01.preprocessing" / "{sample}_fastp.html",
        json = OUT / "01.preprocessing" / "{sample}_fastp.json",
    params:
        min_len = MIN_READ_LEN,
    # fastp scales poorly past ~16 threads; no point asking for more.
    threads: lambda wc: threads_for("fastp")
    resources:
        mem_mb = lambda wc: mem_mb_for("fastp"),
    conda:
        "../../envs/fastp.yaml"
    log:
        LOGS / "fastp" / "{sample}.log"
    priority: 4
    shell:
        """
        fastp --detect_adapter_for_pe --length_required {params.min_len} \
            --cut_front --cut_right --thread {threads} --verbose \
            -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -j {output.json} -h {output.html} > {log} 2>&1
        """
