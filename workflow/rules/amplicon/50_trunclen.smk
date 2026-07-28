# Decides how far to trim R1 and R2 before DADA2 sees them, and writes that decision
# to trunclen.json. Nothing downstream re-derives it: the overlap check further down
# this file, dada_filter (which applies the trim), and dada_length_filter all just
# read this rule's output.
#
# In auto mode, R1 and R2 are each cut using falco's per-base quality table for that
# sample set, in three steps:
#   1. Read-coverage cap — each sample is first capped at the read length that
#      min_read_coverage_pct of its reads still reach. Falco's quality table runs to
#      the LONGEST read in a file however few reads got that far, so without this cap
#      a clean run could pick a position almost no read actually reaches. See
#      50a_pick_trunclen.py for the full reasoning and a worked example.
#   2. Quality-based cut — within what step 1 left standing, cut just before the
#      first position where quality (the median Q1 across samples) drops below
#      q_threshold.
#   3. Overlap check (16S, 18S, gyrB, rpoB only) — the two cuts must still leave
#      enough overlap for DADA2 to merge R1 and R2:
#      truncR1 + truncR2 >= amplicon_length + min_overlap. If they don't,
#      resolve_policy decides whether to extend the cuts, relax q_threshold, or stop.
#
# ITS never reaches step 3: ITS amplicons vary too much in length for a single fixed
# cut to be sensible, so dada_filter overrides truncLen to c(0,0) regardless of what
# this rule computes. The step-2 cut is still written to trunclen.json for the
# record, it just has no effect on an ITS run.
#
# Manual mode (trunc_len.mode: manual) skips all of the above and takes truncLen
# straight from manual_r1 / manual_r2 in the config.

rule pick_trunclen:
    input:
        falco_r1 = expand(
            OUT / "stats" / "falco" / "{sample}_R1_stripped" / "fastqc_data.txt",
            sample=SAMPLES,
        ),
        falco_r2 = expand(
            OUT / "stats" / "falco" / "{sample}_R2_stripped" / "fastqc_data.txt",
            sample=SAMPLES,
        ),
        # probe_json is conditionally included: only when expected_length="auto",
        # so pick_trunclen depends on amplicon_probe, which in turn depends on that
        # marker's probe reference — e.g. fetch_silva_train (16S) or fetch_uchime (ITS).
        # Otherwise [] = no dependency.
        probe_json = [PROBE_JSON] if WANTS_PROBE else [],
    output:
        json = OUT / "stats" / "trunclen.json",
    params:
        amplicon_type   = AMPLICON_TYPE,
        expected_length = config["amplicon"]["expected_length"],
        min_overlap     = config["amplicon"]["min_overlap"],
        mode            = config["amplicon"]["trunc_len"]["mode"],
        q_threshold     = config["amplicon"]["trunc_len"]["q_threshold"],
        q_floor         = config["amplicon"]["trunc_len"]["q_floor"],
        resolve_policy  = config["amplicon"]["trunc_len"]["resolve_policy"],
        manual_r1       = config["amplicon"]["trunc_len"]["manual_r1"],
        manual_r2       = config["amplicon"]["trunc_len"]["manual_r2"],
        probe_stat      = PROBE_LENGTH_STAT,
        min_read_coverage_pct = config["amplicon"]["trunc_len"]["min_read_coverage_pct"],
    log:
        LOGS / "pick_trunclen.log",
    threads: lambda wc: threads_for("pick_trunclen")
    resources:
        mem_mb = lambda wc: mem_mb_for("pick_trunclen"),
    script:
        "../../scripts/amplicon/50a_pick_trunclen.py"
