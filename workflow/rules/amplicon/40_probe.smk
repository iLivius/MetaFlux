# Amplicon length probe — only when amplicon.expected_length == "auto".
#
# 16S (pcr mode):    two-pass cutadapt against SILVA — keep reference sequences
#                    that carry both the fwd primer at 5' and revcomp(rev) at 3'.
#                    Surviving trimmed sequences = the in-silico PCR products.
#
# ITS (direct mode): measure lengths directly from UNITE UCHIME pre-extracted
#                    ITS1 or ITS2 sequences. No primer binding sites needed;
#                    cutadapt is not used. This avoids the unreliable in-silico
#                    PCR approach for ITS, where primers (ITS1F, ITS4) bind in
#                    flanking rRNA regions absent from UNITE sequences.
#
# Both modes compute the same length distribution statistics and write a JSON;
# downstream pick_trunclen reads probe_length_stat from config to select which
# statistic to use as expected_length.
#
# Output is cached under refdb/cache/ and survives across pipeline reruns.
# Filename encodes amplicon type, reference tag, and a hash of both primer
# sequences, so the probe step is skipped whenever those three inputs are unchanged.

rule amplicon_probe:
    input:
        ref_fasta = PROBE_REF_FASTA,
        fwd       = PRIMER_FWD,
        rev       = PRIMER_REV,
    output:
        json      = PROBE_JSON,
        amplicons = PROBE_AMPLICONS_FA,
    params:
        probe_mode  = PROBE_MODE,
        max_err     = config["amplicon"]["cutadapt"]["max_error_rate"],
        pcr_min_len = 30,                       # broad floor for PCR mode; full distribution reported anyway
        amp_type    = AMPLICON_TYPE,
        primer_hash = PRIMER_HASH,
        ref_tag     = PROBE_REF_TAG,
    log:
        LOGS / "amplicon_probe.log",
    conda:
        "../../envs/cutadapt.yaml"
    threads: lambda wc: threads_for("cutadapt")
    resources:
        mem_mb = lambda wc: mem_mb_for("cutadapt"),
    script:
        "../../scripts/amplicon/40a_amplicon_probe.py"
