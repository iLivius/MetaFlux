# Amplicon preprocessing: PhiX removal (bowtie2), primer revcomp, primer
# trimming (Cutadapt), and per-stage read counts.

# ─────────────────── 01.reads/ symlinks (visibility) ───────────
# Symlinks input FASTQs into 01.reads/. Raw files stay in fastq_dir; downstream
# rules (rm_phix, falco_raw, count_reads_raw) read from there directly.
# These symlinks are listed in rule all for visibility only.
rule link_reads:
    input:
        fastq = lambda wc: raw_fastq(wc.sample, int(wc.r)),
    output:
        link = OUT / "01.reads" / "{sample}_R{r}.fastq.gz",
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
# Each file name starts with the two-digit number of the folder whose reads it
# counts, so stats/ sorts in pipeline order like the output folders do:
# 01 = the raw reads (linked in 01.reads/), 02 = after PhiX removal (02.no_phix/),
# 03 = after primer trimming (03.stripped/).
#
# count_reads_raw handles both .gz and plain .fastq input because the run's raw
# FASTQs (in fastq_dir) can be either; count_reads_nophix and count_reads_stripped
# only ever zcat, because their inputs are rm_phix's and trim_primers' outputs,
# which this pipeline always writes gzipped.

rule count_reads_raw:
    input:
        r1 = expand("{fq}", fq=[raw_fastq(s, 1) for s in SAMPLES]),
    output:
        counts = OUT / "stats" / "01.raw_reads.counts",
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
        r1 = expand(OUT / "02.no_phix" / "{sample}_R1.fastq.gz", sample=SAMPLES),
    output:
        counts = OUT / "stats" / "02.nophix_reads.counts",
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
        r1 = expand(OUT / "03.stripped" / "{sample}_R1.fastq.gz", sample=SAMPLES),
    output:
        counts = OUT / "stats" / "03.stripped_reads.counts",
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
#
# WHY --reorder, AND WHY IT IS NOT OPTIONAL HERE
# On several threads bowtie2 writes each pair as soon as the thread handling it
# finishes, and reads take different amounts of work, so the surviving pairs come
# out in a different order every run. Nothing is lost or duplicated, but that order
# reaches DADA2, which numbers ASVs by decreasing abundance and keeps the input
# order for ties — so two ASVs with exactly the same total count can trade IDs
# between two runs of the same data. That was measured on the 16S test set (three
# runs, one swapped pair), and it propagates: a different ASV numbering means a
# different row order in the phylogeny alignment, and the tree search then settles
# on a different, equally good tree. --reorder makes bowtie2 emit the pairs in the
# order they were read, which is what a single-threaded run would produce. Measured on
# one test sample of ~50,000 pairs over six alternating repeats: no reliable difference
# in wall time (1.62 s either way) and about 4-6 MB more memory, because a thread that
# finishes early has to hold its output until the reads ahead of it are done.
rule rm_phix:
    input:
        r1  = lambda wc: raw_fastq(wc.sample, 1),
        r2  = lambda wc: raw_fastq(wc.sample, 2),
        idx = rules.build_phix_index.output.idx,
    output:
        r1 = OUT / "02.no_phix" / "{sample}_R1.fastq.gz",
        r2 = OUT / "02.no_phix" / "{sample}_R2.fastq.gz",
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
                --reorder \
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
#   Then swap the two mates of every pass-B pair (its R2 becomes R1) and
#   concatenate with pass A. Under orientation=mixed the two 5' calls (-g/-G)
#   also need at least 10 primer bases (--overlap 10; see below).
#
# Why pass B exists: some library preps come back with reads in mixed
# orientation — part of a sample's R1 file is actually a reverse read, and part
# of R2 is actually a forward read. Trimming with the forward primer alone
# (pass A only) simply loses those flipped pairs as "primer not found". Setting
# amplicon.primers.orientation: mixed adds pass B, which re-tries with the
# primers swapped. A pair that pass B keeps still has its mates the other way
# round, so its two mates are swapped back before the two passes are
# concatenated into one pair of output files: the recovered pair then has its
# forward read in R1, starting right after the forward primer, and its reverse
# read in R2, like the pairs of pass A, which is what DADA2's merging step
# expects. Nothing is reverse-complemented; the comment at the concatenation
# step below says why that would be wrong. With orientation: fixed (the
# default), only pass A runs.
#
# Why the 5' calls need 10 primer bases under mixed. Both passes read every
# input pair, not only the ones pass A dropped. By default Cutadapt counts a 5'
# primer as found as soon as a read starts with the primer's last three bases
# (minimum overlap 3, because a read may begin partway into its primer), or
# when the whole primer turns up anywhere in the read within the error rate.
# With some primer pairs the start of each primer resembles the end of the
# other: the COI primer BF3 ends in HCG, which a BR2 read matches whenever its
# degenerate third base is a G (TCG), and BR2 ends in YCA, which a BF3 read
# matches whenever its third base is an A (CCA). A pair could then also pass
# the call meant for the other orientation, with only those three bases cut,
# and was written twice: once right, and once with its mates the wrong way
# round and most of a primer still on. Every version up to v2.5.0 did this.
# Measured with Cutadapt's default, 13.3 % and 17.1 % of the pairs of two
# BF3/BR2 libraries of an insect mock were written twice, 0.14 % with the COI
# primers fwhF2/fwhR2n, 0.31 % in the rpoB test set, 0.04 % in the gyrB set,
# 2 pairs in the 18S set and none in the 16S and ITS sets. With 10 bases, 1
# pair is left in one BF3/BR2 library and 32 in the gyrB set, whose primers are
# highly degenerate: pairs in which the other orientation's primer still
# matches 10 or more of its bases; none elsewhere. The 3' calls keep
# Cutadapt's default: a read that runs through a short amplicon into the far
# primer can legitimately end a few bases into it. A primer shorter than 10
# bases needs all of its bases (Cutadapt caps the minimum at the primer's
# length).
#
# What the 10 bases cost: the pairs in which a primer is found only from 3 to 9
# of its bases at the start of a read, which orientation: fixed keeps. That is
# 0.006-0.15 % of the pairs in the test sets of the five shipped markers and at
# most 0.06 % in the four COI libraries measured. Hardly any of them are
# cleanly trimmed amplicon reads: the read cut short carries the rest of a
# primer that differs from the given one by more than the error rate allows, or
# it is a read of something else, such as a low-complexity read. The amplicon
# overview gives the measurements. With orientation: fixed the one pass keeps
# Cutadapt's default and runs exactly as in v2.5.0.
#
# What the swap means for DADA2 in a truly mixed library: the output R1 file
# holds sequencer read 1 of the pairs that came in the right way round and
# sequencer read 2 of the swapped ones (and R2 the other way round). DADA2
# learns one error model for all R1 reads and one for all R2 reads
# (dada_seqtab), and pick_trunclen reads one quality profile per direction, so
# both see a blend of the two sequencing reads instead of one read each. When
# every pair is swapped, the pass-B part of the output is what orientation:
# fixed gives on the two input files swapped, less the few pairs the 10-base
# rule drops, and pass A adds almost nothing.
rule trim_primers:
    input:
        # r1/r2 resolve to 02.no_phix when REMOVE_PHIX else the raw FASTQs.
        unpack(primer_trim_upstream),
        fwd    = PRIMER_FWD,
        rev    = PRIMER_REV,
        fwd_rc = rules.revcomp_primers.output.fwd_rc,
        rev_rc = rules.revcomp_primers.output.rev_rc,
    output:
        r1       = OUT / "03.stripped" / "{sample}_R1.fastq.gz",
        r2       = OUT / "03.stripped" / "{sample}_R2.fastq.gz",
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

        # Under orientation: mixed both 5' calls need at least 10 bases of their
        # primer (--overlap 10) instead of Cutadapt's default of 3, so that a pair
        # is not kept by both passes through a chance match of a few bases (see
        # the rule header). Under orientation: fixed the variable stays empty, so
        # pass A gets exactly the options it had in v2.5.0. It is used unquoted on
        # purpose: empty, it adds nothing to the command line (quoted, it would
        # hand Cutadapt an empty argument); set, it gives the two words --overlap
        # and 10.
        overlap_5p=""
        if [ "{params.orientation}" = "mixed" ]; then
            overlap_5p="--overlap 10"
        fi

        # ---- Pass A 5': assume R1=fwd, R2=rev ----
        cutadapt -j {threads} \
                 -g file:{input.fwd} -G file:{input.rev} \
                 -e {params.max_err} --pair-filter=any $overlap_5p \
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
                     -e {params.max_err} --pair-filter=any $overlap_5p \
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

            # Concatenate pass A and pass B, swapping the mates of every pass-B
            # pair. A pass-B pair came in the other way round: its R1 is the reverse
            # read (it started with the reverse primer) and its R2 the forward read.
            # So the output R1 gets pass-B R2 and the output R2 gets pass-B R1, and
            # the recovered pair has its forward read in R1, starting right after
            # the forward primer, as DADA2 expects. Both output files keep the same
            # pair order, and that order is all DADA2 uses to pair the mates
            # (dada_filter leaves filterAndTrim's matchIDs at its default, FALSE),
            # so an R1 record that keeps the name it had as read 2 is harmless. A
            # pair that both passes kept is written twice here; the 10-base rule of
            # the 5' calls makes that rare (see the rule header).
            #
            # MetaFlux up to v2.5.0 reverse-complemented each pass-B mate instead
            # (seqtk seq -r). That gives the right strand but the wrong pair
            # geometry: R1 then covers the 3' end of the amplicon and R2 the 5' end,
            # and mergePairs with trimOverhang keeps only the stretch where the two
            # overlap. (Their quality strings also run backwards, so truncation cuts
            # the best bases, but the geometry alone is enough: the result was the
            # same without truncation.) On the 418-bp BF3/BR2 COI amplicon of an
            # insect mock in which every R1 carries the reverse primer (SRA run
            # SRR8082166), none of the merged reads came within 9 bp of the full
            # 418 bp that way, against 99.4 % with the mates swapped (99.1 % at
            # exactly 418 bp).
            cat $tmpdir/{wildcards.sample}_after3p_R1.fq.gz \
                $tmpdir/{wildcards.sample}_pb3_R2.fq.gz > {output.r1}
            cat $tmpdir/{wildcards.sample}_after3p_R2.fq.gz \
                $tmpdir/{wildcards.sample}_pb3_R1.fq.gz > {output.r2}
        else
            mv $tmpdir/{wildcards.sample}_after3p_R1.fq.gz {output.r1}
            mv $tmpdir/{wildcards.sample}_after3p_R2.fq.gz {output.r2}
        fi
        """
