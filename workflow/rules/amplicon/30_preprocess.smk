# Amplicon preprocessing: PhiX removal (bowtie2), primer revcomp, primer
# trimming (Cutadapt), and per-stage read counts.

# ─────────────────── 1.reads/ symlinks (visibility) ─────────────
# Symlinks input FASTQs into 1.reads/. Raw files stay in fastq_dir; downstream
# rules (rm_phix, falco_raw, count_reads_raw) read from there directly.
# These symlinks are listed in rule all for visibility only.
rule link_reads:
    input:
        fastq = lambda wc: raw_fastq(wc.sample, int(wc.r)),
    output:
        link = OUT / "1.reads" / "{sample}_R{r}.fastq.gz",
    wildcard_constraints:
        r = "1|2",
    shell:
        r"""
        mkdir -p $(dirname {output.link})
        ln -sfn {input.fastq} {output.link}
        """


# ──────────────────────── Read counts ──────────────────────────
# Counts reads (R1) at each pipeline stage. Format = "{sample} : {n}";
# trailing line "Total read count : N" is skipped by parsers.

rule count_reads_raw:
    input:
        r1 = expand("{fq}", fq=[raw_fastq(s, 1) for s in SAMPLES]),
    output:
        counts = OUT / "stats" / "1.raw_reads.counts",
    log:
        LOGS / "count_reads" / "raw.log",
    shell:
        r"""
        mkdir -p $(dirname {output.counts})
        : > {output.counts}
        total=0
        for fq in {input.r1}; do
          if [[ "$fq" == *.gz ]]; then
            c=$(zcat "$fq" | wc -l)
          else
            c=$(wc -l < "$fq")
          fi
          n=$((c / 4))
          total=$((total + n))
          sample=$(basename "$fq" | sed 's/_R1\.fastq\.gz$//' | sed 's/_R1\.fastq$//')
          echo "$sample : $n" >> {output.counts}
        done
        echo "Total read count : $total" >> {output.counts}
        """


rule count_reads_nophix:
    input:
        r1 = expand(OUT / "2.no_phix" / "{sample}_R1.fastq.gz", sample=SAMPLES),
    output:
        counts = OUT / "stats" / "2.nophix_reads.counts",
    log:
        LOGS / "count_reads" / "nophix.log",
    shell:
        r"""
        mkdir -p $(dirname {output.counts})
        : > {output.counts}
        total=0
        for fq in {input.r1}; do
          c=$(zcat "$fq" | wc -l)
          n=$((c / 4))
          total=$((total + n))
          sample=$(basename "$fq" | sed 's/_R1\.fastq\.gz$//')
          echo "$sample : $n" >> {output.counts}
        done
        echo "Total read count : $total" >> {output.counts}
        """


rule count_reads_stripped:
    input:
        r1 = expand(OUT / "3.stripped" / "{sample}_R1.fastq.gz", sample=SAMPLES),
    output:
        counts = OUT / "stats" / "3.stripped_reads.counts",
    log:
        LOGS / "count_reads" / "stripped.log",
    shell:
        r"""
        mkdir -p $(dirname {output.counts})
        : > {output.counts}
        total=0
        for fq in {input.r1}; do
          c=$(zcat "$fq" | wc -l)
          n=$((c / 4))
          total=$((total + n))
          sample=$(basename "$fq" | sed 's/_R1\.fastq\.gz$//')
          echo "$sample : $n" >> {output.counts}
        done
        echo "Total read count : $total" >> {output.counts}
        """


# ──────────────────────── PhiX removal ─────────────────────────
rule rm_phix:
    input:
        r1  = lambda wc: raw_fastq(wc.sample, 1),
        r2  = lambda wc: raw_fastq(wc.sample, 2),
        idx = rules.build_phix_index.output.idx,
    output:
        r1 = OUT / "2.no_phix" / "{sample}_R1.fastq.gz",
        r2 = OUT / "2.no_phix" / "{sample}_R2.fastq.gz",
    params:
        prefix = str(REFDB_PHIX_PREFIX),
    log:
        # MultiQC's bowtie2 module extracts the sample name from the log
        # basename — suffix '_phixFilter' is stripped back to '{sample}' by
        # the regex rule in multiqc_config.yaml.
        LOGS / "rm_phix" / "{sample}_phixFilter.log",
    benchmark:
        BENCH / "rm_phix" / "{sample}.tsv"
    conda:
        "../../envs/bowtie2.yaml"
    threads: lambda wc: threads_for("bowtie2")
    resources:
        mem_mb = lambda wc: mem_mb_for("bowtie2"),
    shell:
        r"""
        mkdir -p $(dirname {output.r1})
        # bowtie2 --un-conc-gz writes "<basename>.<n>.gz" — use a temp basename then rename.
        tmp=$(mktemp -u -p $(dirname {output.r1}) {wildcards.sample}.unmapped.XXXX)
        bowtie2 -x {params.prefix} \
                -1 {input.r1} -2 {input.r2} \
                --threads {threads} \
                --un-conc-gz ${{tmp}}_R%.fastq.gz \
                -S /dev/null \
                > {log} 2>&1
        mv ${{tmp}}_R1.fastq.gz {output.r1}
        mv ${{tmp}}_R2.fastq.gz {output.r2}
        """


# ──────────────────────── Primer revcomps ──────────────────────
rule revcomp_primers:
    input:
        fwd = PRIMER_FWD,
        rev = PRIMER_REV,
    output:
        fwd_rc = OUT / "_aux" / "primers" / "fwd_primer_rc.fasta",
        rev_rc = OUT / "_aux" / "primers" / "rev_primer_rc.fasta",
    log:
        LOGS / "revcomp_primers.log",
    conda:
        "../../envs/seqtk.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.fwd_rc})
        seqtk seq -r {input.fwd} > {output.fwd_rc} 2> {log}
        seqtk seq -r {input.rev} > {output.rev_rc} 2>> {log}
        """


# ──────────────────────── Primer trimming ──────────────────────
# Cutadapt with cross-paired 3' adapters:
#   Pass A: -g fwd -G rev   +   -a rev_rc -A fwd_rc
#   Pass B (swap, only if orientation=mixed):
#           -g rev -G fwd   +   -a fwd_rc -A rev_rc
#   Then revcomp pass-B output and concat with pass-A.
rule trim_primers:
    input:
        # r1/r2 resolve to 2.no_phix when REMOVE_PHIX else the raw FASTQs.
        unpack(primer_trim_upstream),
        fwd    = PRIMER_FWD,
        rev    = PRIMER_REV,
        fwd_rc = rules.revcomp_primers.output.fwd_rc,
        rev_rc = rules.revcomp_primers.output.rev_rc,
    output:
        r1       = OUT / "3.stripped" / "{sample}_R1.fastq.gz",
        r2       = OUT / "3.stripped" / "{sample}_R2.fastq.gz",
        json_pa5 = OUT / "stats" / "cutadapt" / "{sample}.passA_5prime.cutadapt.json",
        json_pa3 = OUT / "stats" / "cutadapt" / "{sample}.passA_3prime.cutadapt.json",
    params:
        orientation = ORIENTATION,
        max_err     = config["amplicon"]["cutadapt"]["max_error_rate"],
        min_len     = config["amplicon"]["cutadapt"]["min_length"],
    log:
        LOGS / "trim_primers" / "{sample}.log",
    benchmark:
        BENCH / "trim_primers" / "{sample}.tsv"
    conda:
        "../../envs/cutadapt.yaml"
    threads: lambda wc: threads_for("cutadapt")
    shell:
        # MultiQC's cutadapt module reads the JSON's input.path1 field for
        # sample naming. To distinguish 5' and 3' passes (and to keep them
        # from clashing with bowtie2's bare {sample}), each cutadapt call
        # reads from a stage-named symlink, so the JSON records
        # '{sample}_5prime_R1' / '{sample}_3prime_R1' as inputs. The regex
        # in multiqc_config.yaml strips those suffixes back to '{sample}'
        # behind a "Primer trim 5'|3'" prefix.
        r"""
        mkdir -p $(dirname {output.r1}) $(dirname {output.json_pa5})
        tmpdir=$(mktemp -d -p $(dirname {output.r1}) {wildcards.sample}.cutadapt.XXXX)
        trap 'rm -rf "$tmpdir"' EXIT

        # Stage-named symlinks for cutadapt to read from (drives MultiQC sample naming).
        ln -sf "$(realpath {input.r1})" $tmpdir/{wildcards.sample}_5prime_R1.fastq.gz
        ln -sf "$(realpath {input.r2})" $tmpdir/{wildcards.sample}_5prime_R2.fastq.gz

        # ---- Pass A 5': assume R1=fwd, R2=rev ----
        cutadapt -j {threads} \
                 -g file:{input.fwd} -G file:{input.rev} \
                 -e {params.max_err} --pair-filter=any \
                 --minimum-length {params.min_len} --discard-untrimmed \
                 --json {output.json_pa5} \
                 -o $tmpdir/{wildcards.sample}_after5p_R1.fq.gz \
                 -p $tmpdir/{wildcards.sample}_after5p_R2.fq.gz \
                 $tmpdir/{wildcards.sample}_5prime_R1.fastq.gz \
                 $tmpdir/{wildcards.sample}_5prime_R2.fastq.gz \
                 > {log} 2>&1

        # Rename 5' output → 3' input so the 3' pass's JSON records '_3prime' as sample.
        mv $tmpdir/{wildcards.sample}_after5p_R1.fq.gz $tmpdir/{wildcards.sample}_3prime_R1.fq.gz
        mv $tmpdir/{wildcards.sample}_after5p_R2.fq.gz $tmpdir/{wildcards.sample}_3prime_R2.fq.gz

        # ---- Pass A 3': cross-paired 3' adapters ----
        cutadapt -j {threads} \
                 -a file:{input.rev_rc} -A file:{input.fwd_rc} \
                 -e {params.max_err} --pair-filter=any \
                 --minimum-length {params.min_len} \
                 --json {output.json_pa3} \
                 -o $tmpdir/{wildcards.sample}_after3p_R1.fq.gz \
                 -p $tmpdir/{wildcards.sample}_after3p_R2.fq.gz \
                 $tmpdir/{wildcards.sample}_3prime_R1.fq.gz \
                 $tmpdir/{wildcards.sample}_3prime_R2.fq.gz \
                 >> {log} 2>&1

        if [ "{params.orientation}" = "mixed" ]; then
            echo "[trim_primers] orientation=mixed → running swap pass" >> {log}

            # ---- Pass B (swap) 5': assume R1=rev, R2=fwd ----
            # (Swap-pass cutadapt logs go through the same JSON-less log file;
            # not exposed to MultiQC. Tmp paths use sample name to avoid clashes.)
            cutadapt -j {threads} \
                     -g file:{input.rev} -G file:{input.fwd} \
                     -e {params.max_err} --pair-filter=any \
                     --minimum-length {params.min_len} --discard-untrimmed \
                     -o $tmpdir/{wildcards.sample}_pb5_R1.fq.gz \
                     -p $tmpdir/{wildcards.sample}_pb5_R2.fq.gz \
                     {input.r1} {input.r2} \
                     >> {log} 2>&1

            # ---- Pass B (swap) 3': cross-paired 3' adapters ----
            cutadapt -j {threads} \
                     -a file:{input.fwd_rc} -A file:{input.rev_rc} \
                     -e {params.max_err} --pair-filter=any \
                     --minimum-length {params.min_len} \
                     -o $tmpdir/{wildcards.sample}_pb3_R1.fq.gz \
                     -p $tmpdir/{wildcards.sample}_pb3_R2.fq.gz \
                     $tmpdir/{wildcards.sample}_pb5_R1.fq.gz \
                     $tmpdir/{wildcards.sample}_pb5_R2.fq.gz \
                     >> {log} 2>&1

            # Reverse-complement pass-B output so it matches pass-A orientation.
            zcat $tmpdir/{wildcards.sample}_pb3_R1.fq.gz | seqtk seq -r - | gzip \
                > $tmpdir/{wildcards.sample}_pb3_R1_rc.fq.gz
            zcat $tmpdir/{wildcards.sample}_pb3_R2.fq.gz | seqtk seq -r - | gzip \
                > $tmpdir/{wildcards.sample}_pb3_R2_rc.fq.gz

            # Concat pass-A + reverse-complemented pass-B (same orientation).
            cat $tmpdir/{wildcards.sample}_after3p_R1.fq.gz \
                $tmpdir/{wildcards.sample}_pb3_R1_rc.fq.gz > {output.r1}
            cat $tmpdir/{wildcards.sample}_after3p_R2.fq.gz \
                $tmpdir/{wildcards.sample}_pb3_R2_rc.fq.gz > {output.r2}
        else
            mv $tmpdir/{wildcards.sample}_after3p_R1.fq.gz {output.r1}
            mv $tmpdir/{wildcards.sample}_after3p_R2.fq.gz {output.r2}
        fi
        """
