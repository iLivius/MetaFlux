# Shotgun preprocessing: turns raw paired-end reads into the trimmed FASTQs that
# 35_kraken2.smk classifies. Five things happen here, in order:
#   1. download_sra     — only in BioProject mode: fetch FASTQs from SRA instead
#                          of expecting them already in fastq_dir.
#   2. samples_manifest  — write a small TSV recording where each sample's reads
#                          came from.
#   3. decontam_phix     — remove PhiX spike-in reads (BBDuk). Optional.
#   4. fetch_host_refs / build_host_index / decontam_host — remove host (or any
#                          other unwanted) DNA by mapping against a reference
#                          with BBMap. Optional.
#   5. trim_adapters     — adapter and quality trimming (fastp). Always runs,
#                          and its output is what 35_kraken2.smk actually
#                          classifies — every earlier step here exists only to
#                          clean reads up before that point, not to change what
#                          gets classified.
#
# Stage chain (each of steps 3-4 is switched on/off by shotgun.decontamination
# in the config): raw → (phix) → (host) → fastp. The fastp_upstream() /
# host_upstream() helper functions (defined in 00_common.smk) work out which of
# these stages actually ran for this config and hand back the right upstream
# files, so a consumer rule never has to know whether the step before it ran or
# was skipped.
#
# priority: the values in this file count down from 10 (download_sra) to 4
# (trim_adapters), continuing on through 35_kraken2.smk and 45_abundance.smk
# (3, 2, 1) — roughly in pipeline order. When several samples have jobs ready
# to run at once and there are not enough cores for all of them, Snakemake
# uses priority to break the tie, preferring the higher number — so an
# earlier-stage job tends to run before a later-stage one, and the batch of
# samples moves through the pipeline together instead of a few samples racing
# ahead while others are still waiting on their first step.

import shutil
import subprocess
from urllib.parse import urlparse


# ──────────────────────── SRA fetch (BioProject mode) ──────────
# Downloads one SRA run's reads and writes them as this pipeline's standard
# {sample}_R1/_R2.fastq.gz pair. wildcards.sample here is an SRA run accession
# (e.g. SRR12345678) — the run accessions for a BioProject are resolved once,
# at parse time, by _fetch_bioproject_runs() in 00_common.smk and become
# SAMPLES for the rest of the run.
#
# Only triggered when input.bioproject is set (BioProject mode). In
# local-FASTQ mode this rule is never invoked — samples are instead discovered
# by globbing fastq_dir, and the dependency on these same output paths is
# satisfied by the user's own files already sitting there. Either way,
# downstream rules read raw reads through raw_r1()/raw_r2() (00_common.smk),
# so they never need to know which of the two happened.
#
# prefetch downloads the raw .sra archive; fasterq-dump extracts it to paired
# plain FASTQ; pigz then compresses to the .gz this pipeline expects
# everywhere else.
# NOTE: prefetch puts the .sra under {tmp}/{accession}/{accession}.sra (an
# SRA-tools convention, not something chosen here), hence the nested path in
# the fasterq-dump call below.
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
# Writes one TSV row per sample: sample name, where it came from (local FASTQ
# vs BioProject download), the BioProject accession (blank for local runs),
# and the R1/R2 paths. Everything it needs (SAMPLES, INPUT_MODE, BIOPROJECT,
# FASTQ_DIR, FASTQ_NAMING, EXTN) is a global already resolved once in
# 00_common.smk, so this rule has no real file input — it just renders those
# globals to disk. Nothing else in the pipeline reads samples.tsv back in; it
# exists purely as a human/external-script-friendly record of what this run
# actually processed, without having to re-parse the config or re-glob
# fastq_dir.
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
# First decontamination stage — only in the DAG when
# shotgun.decontamination.remove_phix is true (REMOVE_PHIX, 00_common.smk).
# Takes the sample's raw FASTQs (raw_r1/raw_r2) and drops any read that
# matches PhiX, a synthetic virus genome Illumina runs routinely spike into
# a lane as a sequencing control — real, but not part of the sample, so it
# has no business being counted as an organism in the final results.
#
# BBDuk ships phix174_ill.ref.fa.gz inside the bbmap conda package, so
# `ref=phix` works without any download/index-build step.
# k=31 hdist=1 is the BBTools-recommended setting for PhiX (see their
# decontamination guide): match on 31-mers, allowing up to 1 mismatch
# (hdist = hamming distance) per k-mer. stats= produces a small text file
# that MultiQC picks up automatically, and that aggregate_read_counts
# (80_stats.smk) also parses for the raw/nophix read-tracking columns.
#
# Output feeds fastp_upstream()/host_upstream() (00_common.smk), which route
# decontam_host or trim_adapters to these _dephix_R1/R2 files whenever this
# stage ran, instead of the raw FASTQs.
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
# Turns shotgun.decontamination.host_genomes (a list of local paths and/or
# URLs from the config, resolved to HOST_GENOMES in 00_common.smk) into one
# concatenated host_combined.fna.gz that build_host_index (below) then
# indexes. Only reached when host_genomes is non-empty: fastp_upstream()/
# host_upstream() (00_common.smk) never point at the host-removal stage
# otherwise, so on a run with no host_genomes configured, this rule and the
# two below it simply never get pulled into the DAG.
#
# Each source is downloaded (URL) or copied/gzipped (local path) into its own
# temp file, then binary-concatenated in that order — gzipped files
# concatenate to a valid gzip stream, so this avoids re-gzipping a
# multi-gigabyte human genome just to merge it with the other references.
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
# Build the BBMap index ONCE from fetch_host_refs's concatenated host
# reference, then reuse it for every sample (far cheaper than rebuilding per
# sample) — decontam_host (below) reads this directory
# (rules.build_host_index.output.idx) as one of its inputs for every sample.
# build=1 writes the index under {path}/ref/. mem_mb is sized in the config
# (resources.mem_mb.build_host_index) for a human-masked reference
# (~24-28 Gb); a smaller/different host reference would need less.
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
# divergent host reads (cf. removehuman.sh), i.e. reads that don't match the
# host reference exactly but still come from it. Reads that do NOT map to the
# host index (outu) are kept; mapped reads are dropped. minid (params.min_id,
# from shotgun.decontamination.host_min_id in the config) is the identity
# cutoff — how similar a read has to be to the host reference before BBMap
# calls it host DNA and drops it.
# quickmatch/fast are BBMap's own speed presets for this kind of large-scale
# decontamination mapping; minhits/qtrim/trimq/untrim are further tuning
# flags from the same BBMap guide (minimum seed hits before attempting an
# alignment; quality-trim both read ends at Q10 before mapping, then restore
# the trimmed bases in the output so read length isn't silently changed).
# Flag order here is not incidental: `fast` is a macro that adjusts several
# alignment parameters at once, including bandwidth-related ones, so it could
# in principle clash with the explicit bwr=0.16/bw=12 set just before it.
# BBTools' own removehuman.sh and removecatdogmousehuman.sh (shipped in the
# same bbtools conda package) use this identical flag sequence — bwr/bw before
# quickmatch/fast — for exactly this purpose, so this ordering is the
# tool authors' own tested pattern, not a coincidence worth "cleaning up".
# (quickmatch itself is a plain boolean toggle, not a macro, and does not
# touch bwr/bw at all.)
# The surviving-pair count is written to stats for read-tracking (BBMap, unlike
# BBDuk, emits no #Total/#Matched table, so we count the kept R1 records).
# Feeds fastp_upstream() (00_common.smk), which routes trim_adapters to these
# _dehost_R1/R2 files whenever host removal ran.
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
# Last preprocessing stage — always runs, regardless of which earlier stages
# were skipped. fastp_upstream() (00_common.smk) already worked out which of
# raw / PhiX-removed / host-removed reads is the last stage that actually ran
# for this config, so this rule doesn't repeat that REMOVE_PHIX/HOST_GENOMES
# logic itself. Its _trim_R1/_trim_R2 output is what kraken2 (35_kraken2.smk)
# actually classifies.
#
# --cut_front / --cut_right is sliding-window quality trimming from both ends
# of each read; --length_required (min_len, from
# shotgun.fastp.min_read_length in the config) drops any read too short to be
# useful after trimming. fastp emits .fastq.gz natively when the output
# filename ends in .gz — no separate pigz step needed here, unlike the
# BBDuk/BBMap/Kraken2 outputs elsewhere in this file.
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
