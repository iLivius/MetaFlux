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
# These three counts files (raw / nophix / stripped) are read back later by
# aggregate_read_counts (80_taxonomy.smk), which joins them with the DADA2 and
# taxonomy counts into one read-tracking table for the whole run.
#
# count_reads_raw handles both .gz and plain .fastq input because the run's raw
# FASTQs (in fastq_dir) can be either; count_reads_nophix and count_reads_stripped
# only ever zcat, because their inputs are rm_phix's and trim_primers' outputs,
# which this pipeline always writes gzipped.

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
# Aligns each sample's raw read pair against the PhiX bowtie2 index (built by
# build_phix_index, 10_refdb.smk) and keeps only the pairs that DON'T align —
# bowtie2's --un-conc-gz writes the pairs that fail to align concordantly (i.e.
# as a proper pair) to a gzipped FASTQ pair, which is exactly the "not PhiX"
# read set this rule wants. Runs only when amplicon.decontamination.remove_phix
# is true; when it's false, trim_primers reads straight from the raw FASTQs
# instead (see primer_trim_upstream in 00_common.smk). Output feeds trim_primers
# and the "nophix" falco QC stage.
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
# Reverse-complements the fwd/rev primer FASTAs once, up front, with seqtk.
# trim_primers (below) needs these to catch primer read-through at the 3' end:
# if a read runs past the amplicon and into the OTHER primer site, that primer
# shows up there in its reverse-complement orientation, not its original one.
# Computing the revcomps here (once) rather than inside trim_primers (once per
# sample) avoids repeating the same seqtk call for every sample in the run.
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
# Removes the fwd/rev primers from each read pair with Cutadapt, using
# cross-paired 3' adapters:
#   Pass A: -g fwd -G rev   +   -a rev_rc -A fwd_rc
#   Pass B (swap, only if orientation=mixed):
#           -g rev -G fwd   +   -a fwd_rc -A rev_rc
#   Then revcomp pass-B output and concat with pass-A.
#
# Why pass B exists: some library preps come back with reads in mixed
# orientation — part of a sample's R1 file is actually a reverse read, and part
# of R2 is actually a forward read. Trimming with the forward primer alone
# (pass A only) simply loses those flipped pairs as "primer not found". Setting
# amplicon.primers.orientation: mixed adds pass B, which re-tries with the
# primers swapped, then reverse-complements its output so it lines back up with
# pass A's orientation before the two are concatenated into one output file.
# With orientation: fixed (the default), only pass A runs.
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
        # Why the symlinks below, instead of just handing cutadapt the real files:
        #
        # Each sample goes through cutadapt twice here — once to cut the forward
        # primer off the front of the reads (5' pass), once to cut the reverse
        # primer's leftover off the back (3' pass). MultiQC has to show both in the
        # report as two separate rows, but it does not take the sample name from
        # the output filename we choose — it digs it out of cutadapt's own JSON
        # stats file, from the path of the INPUT file cutadapt was given. So
        # whatever name we hand cutadapt as input is the name that ends up in the
        # report, whether we intended that or not.
        #
        # If both passes were simply given "{sample}", both JSON files would report
        # the same name, and the report would show only one row where two are
        # needed (and it could even be confused with bowtie2's own PhiX-removal
        # stats, which already use the bare "{sample}" name).
        #
        # The fix: before each cutadapt call, we make a symlink — a renamed pointer
        # to the same file, nothing is copied — that spells out which pass this is,
        # e.g. "{sample}_5prime_R1.fastq.gz". Cutadapt reads through the symlink
        # and writes THAT name into its JSON. multiqc_config.yaml then has a rule
        # that recognises the "_5prime"/"_3prime" tag and rewrites it back to a
        # clean "Primer trim 5'|3' | {sample}" label in the final report. Only the
        # report's bookkeeping passes through this detour — the real output files
        # keep their normal "{sample}_R1.fastq.gz" names.
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
