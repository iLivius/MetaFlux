# sylph genome profiling: the second shotgun profiler, used when shotgun.profiler is
# 'sylph' or 'both'. The Snakefile leaves this file out under 'kraken2' (see
# SKIP_RULE_FILES in 00_common.smk), so with the default none of these rules exist.
#
# HOW sylph DIFFERS FROM KRAKEN2
#   Kraken2 labels every read pair, one by one. sylph labels no read at all. It reduces
#   each sample to a "sketch" (a fixed 1-in-200 subset of the sample's 31-bp k-mers)
#   and compares it with sketches of whole reference genomes. From how many of a
#   genome's k-mers turn up, and how often, it estimates the genome's identity to what
#   is in the sample (only genomes at 95 % ANI or more, roughly the same species, are
#   reported) and its coverage. Coverage x genome size then gives the share of the
#   sample's sequence that genome explains. So sylph's numbers are ESTIMATES derived
#   from coverage, not counts of classified reads. With no per-read output, read
#   extraction (which picks reads by the taxon Kraken2 gave each one) has nothing to
#   work from, and the k-mer evidence gate filters the Kraken2 report, which only
#   Kraken2 writes: both are switched off under profiler: sylph (00_common.smk).
#
# THE RULES, in data-flow order:
#   sylph_sketch        per sample: trimmed reads
#                         -> 02.classification/sylph/sketches/{sample}.paired.sylsp
#   sylph_profile       ONE job: every sketch + every database
#                         -> 03.abundance/sylph_profile.tsv (one row per sample x genome)
#                         -> 02.classification/sylph/sylph_profile.per_sample.txt
#                            (sylph's per-sample lines)
#   sylph_taxprof       ONE job: sylph-tax adds each genome's lineage
#                         -> 02.classification/sylph/taxprof/{sample}.sylphmpa
#   sylph_tables        ONE job: the tables a user opens
#                         -> 03.abundance/sylph_table.tsv (+ .biom, .genome_copies.tsv),
#                            03.abundance/sylph_genomes.tsv, stats/sylph/ (MultiQC)
#   sylph_read_tracking profiler: sylph only -> stats/read_tracking.txt (under kraken2
#                         and both, aggregate_read_counts in 80_stats.smk writes it)
#
# THE READS: exactly the trimmed read pairs Kraken2 gets (trim_adapters in
# 30_preprocess.smk: after PhiX removal, host removal and fastp), so the two profilers
# can be compared on identical reads. Under profiler: both, Snakemake keeps those
# temporary files until both kraken2 and sylph_sketch have read them. Nothing in the
# Kraken2 chain reads the files made here, and nothing here reads Kraken2's: the two
# chains run side by side and their tables are never merged.
#
# THE SKETCHES are kept (not temporary): 3-7 MB per sample for 5 million read pairs,
# but about 0.5 GB for a very deep soil library, so they count on disk. Each sketch
# also records the full path of the trimmed R1 file it was made from. That changes no
# result, but it means two sketches of the same reads made in two different out_dirs
# are not byte-identical: compare profiles and tables between runs, not sketches.
#
# priority: 3 for the sketch (the same as kraken2, so the temporary trimmed reads are
# freed early), then 2 and 1, continuing the countdown started in 30_preprocess.smk.


# ──────────────────────── Sketch one sample ────────────────────
# Turns one sample's trimmed read pairs into a sylph sketch.
#   -c 200  keep about 1 k-mer in 200 (SYLPH_SKETCH_C, fixed in 00_common.smk). A
#           sample sketch can only be compared with databases built at the same or a
#           coarser c: c 200 and c 1000 databases work, a c 100 one does not (sylph
#           then skips every sample; see sylph_profile). 200 is also the density of
#           the prebuilt c 200 databases, which find low-coverage genomes best.
#   -S      store the MetaFlux sample name in the sketch, so sylph's profile says
#           'SRR29240428' rather than a file path, and sylph-tax names its file after
#           it. It also names the sketch file: {sample}.paired.sylsp.
# ONE thread by default (resources.threads.sylph_sketch), and that is a
# reproducibility setting, not a speed one: with 2 or more threads sylph removes
# duplicate k-mers in parallel and the sketch differs slightly from run to run
# (measured on the Zymo mock at 8 threads: 4,395,232 vs 4,395,237 k-mer counts, and
# abundances moved in the 4th decimal). At one thread the sketch is identical every
# time, and it is still fast (about 6 s for 4.7 million read pairs); Snakemake
# sketches several samples at once. That holds with PhiX or host removal on too:
# BBDuk and BBMap write the reads they keep in the order they read them (ordered=t in
# 30_preprocess.smk), so the reads reach fastp, and then this rule, in the same order
# on every run.
# Two failures are turned into a stop with a plain message:
#   - sylph cannot sketch (exit code not 0). The usual cause is a sample with no read
#     pairs left after fastp: sylph cannot sketch an empty file.
#   - sylph wrote a WARN or ERROR line (in sylph 1.0.0 the sketch step warns only
#     about invalid read records and odd file names), so nothing is skipped silently.
#     sylph writes the level as a word between spaces ('... WARN  [sylph::sketch] ...'),
#     which is what the grep below looks for, here and in sylph_profile.
# Output: the sketch, read by sylph_profile next.
rule sylph_sketch:
    input:
        r1 = OUT / "01.preprocessing" / "{sample}_trim_R1.fastq.gz",
        r2 = OUT / "01.preprocessing" / "{sample}_trim_R2.fastq.gz",
    output:
        sketch = OUT / "02.classification" / "sylph" / "sketches" / "{sample}.paired.sylsp",
    params:
        # The output's folder, written out because sylph -d takes a folder and names the
        # sketch file itself (from -S). snakemake --lint flags it as a hard-coded output
        # prefix; that is harmless here.
        outdir = str(OUT / "02.classification" / "sylph" / "sketches"),
        c      = SYLPH_SKETCH_C,
    threads: lambda wc: threads_for_or("sylph_sketch", 1)
    resources:
        mem_mb = lambda wc: mem_mb_for("sylph_sketch"),
    conda:
        "../../envs/sylph.yaml"
    log:
        LOGS / "sylph_sketch" / "{sample}.log"
    priority: 3
    shell:
        r"""
        sylph --version > {log}
        if ! sylph sketch -1 {input.r1} -2 {input.r2} -S {wildcards.sample} -c {params.c} \
                -t {threads} -d {params.outdir} >> {log} 2>&1; then
            echo "[MetaFlux] sylph_sketch: sylph could not sketch the trimmed reads of" \
                 "{wildcards.sample} (details in {log}). The usual cause is a sample with no read" \
                 "pairs left after fastp (see 01.preprocessing/{wildcards.sample}_fastp.html): sylph" \
                 "cannot sketch an empty file. Remove the sample from fastq_dir, or run it with" \
                 "shotgun.profiler: kraken2." >&2
            exit 1
        fi
        if grep -E " (WARN|ERROR) " {log} > /dev/null; then
            echo "[MetaFlux] sylph_sketch: sylph warned while sketching {wildcards.sample}" \
                 "(for example, it skipped an invalid read record); see {log}." >&2
            exit 1
        fi
        """


# ──────────────────────── Profile all samples ──────────────────
# Compares every sample's sketch with every database in references.sylph_databases,
# in ONE sylph run. One run matters for two reasons: the databases are loaded into
# memory once for the whole run (GTDB r232 c200 alone is 26 GB) rather than once per
# sample; and genomes from different databases compete for the same k-mers, so an
# organism present in two databases is reported once, by its closer genome, instead of
# being counted twice. That holds within a sample. When two databases hold the SAME
# genome (GlobDB holds every GTDB genome, under its accession without the version
# suffix), the two copies tie, and which one sylph reports can change from sample to
# sample: the genome then ends up on two rows of the tables. So 00_common.smk refuses
# GTDB and GlobDB in the same run, before any job starts.
#
# -u (estimate unknown) is always on. Without it, sylph stretches each sample's
# Sequence_abundance to add up to 100 %, as if the database explained every read.
# With it, Sequence_abundance is a share of ALL the sample's reads, and what is left
# over is the unknown part that sylph_tables writes as the 'unassigned' row. One
# exception: when the coverage estimates of a sample add up to all of its sequence or
# more, sylph caps the explained share at 100 % and shares that among the detected
# genomes as it does without -u, so every genome of that sample is scaled down, by an
# amount sylph does not report. -u does not change which genomes are reported or
# their Taxonomic_abundance. It also renames column 6 from Eff_cov to True_cov, which
# sylph_tables uses to check -u was on.
# No thresholds are passed, so sylph 1.0.0's defaults apply. The thread count changes
# no number (tested at 4 and 16 threads). One thing is not fixed from run to run: the
# order in which sylph lists genomes with exactly the same abundance (the Zymo mock,
# whose strains are mixed in equal amounts, has several such ties). That touches only
# the raw sylph_profile.tsv and sylph-tax's files; the MetaFlux tables sort their own
# rows, so they come out byte-identical. To compare two raw profiles, sort them first.
#
# The databases are inputs, not params, so a replaced database file re-runs this step
# (and the steps after it), but not the sketches. sylph holds them in memory: the
# memory request (SYLPH_PROFILE_MEM_MB) is worked out at parse time from their size on
# disk, in 00_common.smk.
# The benchmark file (benchmarks/sylph_profile.tsv) under-reports this job's memory, so
# it must not be used to size resources.mem_mb.sylph_profile. Snakemake checks memory
# every 0.5 s for its first 30 checks (about the first 15-20 s), then only every 30 s,
# while sylph's memory climbs as it reads the databases and peaks in the last seconds
# of the job. The recorded max_rss therefore usually misses the peak (measured: max_rss
# 18,845 MB recorded, against a real peak of 25,671 MB in the same units, 26.9 GB).
# Size it from the database files instead: about 1.1 x their size, plus headroom.
#
# sylph exits with code 0 in three situations that must not pass silently, so each
# one is turned into a stop here:
#   - A database built at a lower c than the sketches (c 100, e.g. some viral
#     databases): sylph writes an ERROR line, skips EVERY sample and leaves a profile
#     with a header only.
#   - The same genomes loaded from two databases: sylph writes a WARN line, and which
#     copy of each genome it reports then changes from run to run.
#   - Fewer per-sample result lines than samples: some sample was not profiled.
# A sample in which no genome passes sylph's thresholds is NOT an error: with -u sylph
# still writes its result line ("0.00% of reads detected"), and sylph_tables gives it a
# column of zeros with all its read pairs unassigned.
# Outputs: the profile (read by sylph_taxprof and sylph_tables) and a copy of sylph's
# per-sample log lines (read by sylph_tables: genomes per sample, % of reads detected,
# and which samples got sylph's fixed 99.5 % read accuracy because their genomes were
# covered too thinly for its own estimate).
rule sylph_profile:
    input:
        sketches  = expand(OUT / "02.classification" / "sylph" / "sketches" / "{s}.paired.sylsp", s=SAMPLES),
        databases = SYLPH_DB_PATHS,
    output:
        profile    = OUT / "03.abundance" / "sylph_profile.tsv",
        per_sample = OUT / "02.classification" / "sylph" / "sylph_profile.per_sample.txt",
    params:
        n_samples = len(SAMPLES),
        c         = SYLPH_SKETCH_C,
    threads: lambda wc: threads_for("sylph_profile")
    resources:
        mem_mb = SYLPH_PROFILE_MEM_MB,            # computed in 00_common.smk
    conda:
        "../../envs/sylph.yaml"
    log:
        LOGS / "sylph_profile.log"
    benchmark:
        BENCH / "sylph_profile.tsv"
    priority: 2
    shell:
        r"""
        sylph --version > {log}
        if ! sylph profile {input.databases} {input.sketches} -u -t {threads} \
                -o {output.profile} >> {log} 2>&1; then
            echo "[MetaFlux] sylph_profile: sylph stopped (its message is at the end of {log})." \
                 "If it names a database file, check that the download is complete." >&2
            exit 1
        fi
        if grep -E " (WARN|ERROR) " {log} > /dev/null; then
            echo "[MetaFlux] sylph_profile: sylph wrote warnings or errors (below; full log {log})." >&2
            echo "  Usual causes: a database built at a lower c than MetaFlux's sketches (c {params.c})," >&2
            echo "  which makes sylph skip every sample and still exit 0; or the same genomes loaded" >&2
            echo "  from two databases, which makes its output change from run to run." >&2
            grep -E " (WARN|ERROR) " {log} >&2
            exit 1
        fi
        n_done=$(grep -c -F "% of reads detected in database by profile" {log} || true)
        if [ "$n_done" -ne {params.n_samples} ]; then
            echo "[MetaFlux] sylph_profile: sylph reported results for $n_done of {params.n_samples}" \
                 "samples (see {log}); stopping so that no sample goes missing from the tables" \
                 "without notice." >&2
            exit 1
        fi
        grep -F -e "genomes passing profiling threshold" \
                -e "% of reads detected in database by profile" \
                -e "as read accuracy estimate" {log} > {output.per_sample}
        """


# ──────────────────────── Add lineages (sylph-tax) ─────────────
# sylph reports genome FILES (gtdb_genomes_reps_r232/.../GCF_000025985.1_genomic.fna.gz).
# sylph-tax looks each one up in the lineage tables named by the 'taxonomy' entries of
# references.sylph_databases (GTDB_r232, FungiRefSeq-latest, ... or your own table)
# and writes one {sample}.sylphmpa per sample: a row for every rank from domain to
# species, then one row per genome (t__<accession>).
#   --no-config --taxonomy-dir DIR  read the built-in tables from DIR
#        (references.sylph_taxonomy_dir) and nothing in the user's home folder.
#        Without --no-config, sylph-tax reads, and creates if missing,
#        ~/.config/sylph-tax/config.json, so the run would depend on whatever an
#        earlier manual sylph-tax session left there. Both options must come BEFORE
#        the word taxprof (placed after it, sylph-tax stops with exit code 2).
#   -t   the taxonomies, one per database, in config order (each named once).
#   -o   a file-name prefix, not a folder: '<folder>/' gives <folder>/{sample}.sylphmpa.
#        sylph-tax does not create the folder, hence the mkdir.
# The output is a whole folder rather than one file per sample, because a sample in
# which sylph detected no genome gets no file at all. Snakemake empties the folder
# before this step runs, so files from an earlier run never mix with new ones.
# The lineage tables are inputs: a table whose content changes re-runs this step and
# sylph_tables, but not sylph_profile. (Snakemake compares small inputs, up to 1 MB,
# by content, so merely touching the small fungal table does not re-run anything.)
# A genome missing from the named tables gets a NO_TAXONOMY lineage here, with only a
# warning and exit code 0; sylph_tables stops on it with a message that says what to fix.
# Output: the folder of .sylphmpa files, read by sylph_tables.
rule sylph_taxprof:
    input:
        profile  = OUT / "03.abundance" / "sylph_profile.tsv",
        metadata = SYLPH_TAXONOMY_FILES,
    output:
        taxprof = directory(OUT / "02.classification" / "sylph" / "taxprof"),
    params:
        # '--taxonomy-dir <folder>', or nothing when every taxonomy is a file path.
        taxonomy_dir_option = SYLPH_TAXONOMY_DIR_OPTION,
        # sylph-tax names (GTDB_r232) or table paths, in config order, each once.
        taxonomies          = SYLPH_TAXONOMY_ARGS,
    threads: lambda wc: threads_for("sylph_taxprof")
    resources:
        mem_mb = lambda wc: mem_mb_for("sylph_taxprof"),
    conda:
        "../../envs/sylph.yaml"
    log:
        LOGS / "sylph_taxprof.log"
    priority: 1
    shell:
        r"""
        sylph-tax --version > {log}
        mkdir -p {output.taxprof}
        if ! sylph-tax --no-config {params.taxonomy_dir_option} taxprof {input.profile} \
                -t {params.taxonomies} -o {output.taxprof}/ >> {log} 2>&1; then
            echo "[MetaFlux] sylph_taxprof: sylph-tax stopped; its message is in {log}." >&2
            exit 1
        fi
        """


# ──────────────────────── Build the tables ─────────────────────
# Turns sylph's per-genome percentages into the tables a user opens. The arithmetic,
# the checks and the exact columns are explained in scripts/shotgun/38a_sylph_tables.py:
#   estimated read pairs = Sequence_abundance / 100 x read pairs after QC
#   'unassigned'         = read pairs after QC - all genome rows
# so every column of sylph_table.tsv adds up to the 'trimmed' count of
# stats/read_tracking.txt (one logged exception: a sample sylph explains to 100 %).
# Inputs: the sylph-tax folder (lineages and the two abundances of every genome), the
# raw profile (to check that -u was on, and each genome's FASTA header for the name
# notes), sylph's per-sample lines (genomes per sample, which samples got sylph's fixed
# read-accuracy value) and the fastp JSONs (read pairs after QC). The params also
# carry the lineage tables sylph_taxprof used; the script opens them only to name the
# table a wrong lineage came from, when it has to stop.
# Outputs: the five files listed in _shotgun_targets() (00_common.smk), plus the
# MultiQC copy of the summary, which rule multiqc (90_report.smk) shows as a table.
#
# Runs in the biom env (envs/biom.yaml), which it shares with finalize_otu_table, the
# rule that writes otu_table.biom: both need biom-format, and this script needs nothing
# else beyond the Python standard library and numpy. As with every env file, any edit
# to envs/biom.yaml, even a comment, counts as new software: in finished runs this step
# runs again, and with it MultiQC, which shows the summary this step writes; where
# Kraken2 runs too (profiler: both), finalize_otu_table also runs again.
rule sylph_tables:
    input:
        taxprof     = OUT / "02.classification" / "sylph" / "taxprof",
        profile     = OUT / "03.abundance" / "sylph_profile.tsv",
        per_sample  = OUT / "02.classification" / "sylph" / "sylph_profile.per_sample.txt",
        fastp_jsons = expand(OUT / "01.preprocessing" / "{s}_fastp.json", s=SAMPLES),
    output:
        table   = OUT / "03.abundance" / "sylph_table.tsv",
        biom    = OUT / "03.abundance" / "sylph_table.biom",
        copies  = OUT / "03.abundance" / "sylph_table.genome_copies.tsv",
        genomes = OUT / "03.abundance" / "sylph_genomes.tsv",
        summary = OUT / "stats" / "sylph" / "sylph_summary.tsv",
        mqc     = OUT / "stats" / "sylph" / "sylph_summary_mqc.tsv",
    params:
        samples   = SAMPLES,            # fixes the column order
        # 'file name (taxonomy)' of every database, written into line 1 of the tables,
        # so a table always says which databases its names came from.
        databases = SYLPH_DB_LABELS,
        # The taxonomies sylph_taxprof used and the lineage table behind each (same
        # order), so that a stop on an odd lineage can name the table it came from.
        taxonomy_names  = SYLPH_TAXONOMY_ARGS,
        taxonomy_tables = SYLPH_TAXONOMY_FILES,
    threads: lambda wc: threads_for("sylph_tables")
    resources:
        mem_mb = lambda wc: mem_mb_for("sylph_tables"),
    conda:
        "../../envs/biom.yaml"
    log:
        LOGS / "sylph_tables.log"
    priority: 1
    script:
        "../../scripts/shotgun/38a_sylph_tables.py"


# ──────────────────────── Read tracking (profiler: sylph only) ─
# stats/read_tracking.txt shows how many read pairs each sample kept at every
# cleaning step: raw -> [nophix] -> [no_host] -> trimmed. Under kraken2 and both,
# aggregate_read_counts (80_stats.smk) writes this file, with two Kraken2 columns
# added. Under profiler: sylph that rule file is not read, so this rule writes the
# same file without the Kraken2 columns. It exists only when Kraken2 does not run,
# so the two rules never compete for the same output.
# Inputs: the fastp JSONs, plus BBDuk's PhiX log and the BBMap host stats when
# those steps ran (REMOVE_PHIX, HOST_GENOMES). The parsing is a copy of the Kraken2
# rule's script (scripts/shotgun/38b_sylph_read_tracking.py says why).
# Output: stats/read_tracking.txt, a final file for the user; nothing reads it back.
if not RUN_KRAKEN2:

    rule sylph_read_tracking:
        input:
            fastp_jsons = expand(OUT / "01.preprocessing" / "{s}_fastp.json", s=SAMPLES),
            **( {"phix_logs": expand(LOGS / "decontam_phix" / "{s}.log", s=SAMPLES)}
                if REMOVE_PHIX else {} ),
            **( {"host_stats": expand(OUT / "01.preprocessing" / "{s}_dehost_stats.txt", s=SAMPLES)}
                if HOST_GENOMES else {} ),
        output:
            counts = OUT / "stats" / "read_tracking.txt",
        params:
            samples     = SAMPLES,
            remove_phix = REMOVE_PHIX,
            remove_host = bool(HOST_GENOMES),
        threads: lambda wc: threads_for("sylph_read_tracking")
        resources:
            mem_mb = lambda wc: mem_mb_for("sylph_read_tracking"),
        conda:
            "../../envs/python_utils.yaml"
        log:
            LOGS / "sylph_read_tracking.log"
        script:
            "../../scripts/shotgun/38b_sylph_read_tracking.py"
