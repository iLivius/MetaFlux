# Taxonomy assignment and final output table generation.
#
# assign_taxonomy       : two implementations selected at parse time by
#                         amplicon.taxonomy.method (rdp | sintax):
#                           rdp    — DADA2 assignTaxonomy (RDP classifier) + addSpecies (16S only)
#                           sintax — vsearch --sintax (memory-efficient; ~2–3 GB regardless of ASV count)
#                         Both paths apply the same contaminant filter and write
#                         the same four output tables: the three filtered ones,
#                         and asv_table_unfiltered.txt with every ASV that reached
#                         the filter (organelles included), its Metaxa2 origin and
#                         what the filter did to it. The ranks of every table are
#                         the marker pack's rank model, i.e. the ranks of the
#                         marker's reference (16S: those of the chosen SILVA
#                         release, Domain, Kingdom, Phylum ... on SILVA 144 and
#                         Domain, Phylum ... on SILVA 138.2). Either path starts
#                         only once rule check_taxonomy_reference (10_refdb.smk)
#                         has found that the reference holds exactly those ranks.
# aggregate_read_counts : joins all per-stage read-count sources into one wide
#                         TSV (one row per sample, one column per stage).

_tax_cfg = config["amplicon"]["taxonomy"]

if TAXONOMY_METHOD == "rdp":
    # RDP path: DADA2's assignTaxonomy is a naive Bayesian classifier trained on
    # k-mers from the reference; min_boot is the minimum bootstrap confidence
    # (out of 100 resamples) a rank must clear to be reported, otherwise that
    # rank is left blank rather than guessed. addSpecies does a separate,
    # near-exact sequence match against a species-level reference and runs for
    # 16S only — species_db arrives as an empty list for every other marker, so
    # 80a_assign_taxonomy.R skips that step for them entirely.
    rule assign_taxonomy:
        input:
            seqs         = _seqs_for_taxonomy,
            seqtab_names = _seqtab_for_taxonomy,
            refdb        = TAXONOMY_REFDB,
            species_db   = TAXONOMY_SPECIES_DB,    # addSpecies for 16S only; empty list for all other markers (ITS, 18S, gyrB, rpoB)
            # The report of the reference check (check_taxonomy_reference,
            # 10_refdb.smk), written only when refdb holds the ranks the marker pack
            # declares. The script does not read it: it is listed so that
            # classification waits for the check and never starts on a reference
            # that failed it.
            reference_check = OUT / "stats" / "taxonomy_reference_check.json",
            # Metaxa2's origin call per ASV (target_extract), for the organelle filter
            # and the metaxa2_origin column. Only when Metaxa2 runs (16S with extraction
            # on); for every other marker the rule's inputs are unchanged.
            **({"metaxa2_origin": METAXA2_ORIGIN_TSV} if METAXA2_RUNS else {}),
        output:
            asv_table       = OUT / "06.taxonomy" / "asv_table.txt",
            asv_table_seqs  = OUT / "06.taxonomy" / "asv_table_seqs.txt",
            taxon_seq_table = OUT / "06.taxonomy" / "taxon_seq_table.txt",
            # Every ASV that reached the filter, with metaxa2_origin and filter_status.
            asv_table_unfiltered = OUT / "06.taxonomy" / "asv_table_unfiltered.txt",
        params:
            amp_type        = AMPLICON_TYPE,
            min_boot        = _tax_cfg["min_boot"],
            try_rc          = _tax_cfg["try_rc"],
            tax_levels      = TAX_LEVELS,
            rank_prefixes   = RANK_PREFIXES,
            prefix_style    = PREFIX_STYLE,
            filter_enabled  = FILTER_ENABLED,
            filter_keep     = FILTER_KEEP,
            filter_discard  = FILTER_DISCARD,
            discard_min_support = FILTER_DISCARD_MIN_SUPPORT,   # best-hit check; None = off
            organelle_filter = ORGANELLE_FILTER,   # remove Metaxa2's mitochondria/chloroplast calls
            seed            = config["amplicon"]["seed"],
        log:
            LOGS / "assign_taxonomy.log",
        conda:
            "../../envs/taxonomy.yaml"
        threads: lambda wc: threads_for("assign_taxonomy")
        resources:
            mem_mb = lambda wc: mem_mb_for("assign_taxonomy"),
        script:
            "../../scripts/amplicon/80a_assign_taxonomy.R"

else:  # sintax
    # SINTAX path (the shipped default): vsearch --sintax is a k-mer-based
    # classifier that reports a bootstrap-like confidence per rank; sintax_cutoff
    # (0-1) is the confidence threshold a rank must clear to be reported, playing
    # the same role as min_boot above but on a 0-1 scale instead of 0-100
    # (config note: 0.8 is roughly equivalent to min_boot 80). It is the default
    # here mainly for memory: vsearch holds ~2-3 GB regardless of ASV count,
    # where assignTaxonomy's RAM scales with reference size and ASV count.
    rule assign_taxonomy:
        input:
            seqs         = _seqs_for_taxonomy,
            seqtab_names = _seqtab_for_taxonomy,
            refdb        = TAXONOMY_SINTAX_DB,
            # The report of the reference check — same as the rdp rule above.
            reference_check = OUT / "stats" / "taxonomy_reference_check.json",
            # Metaxa2's origin call per ASV — same as the rdp rule above.
            **({"metaxa2_origin": METAXA2_ORIGIN_TSV} if METAXA2_RUNS else {}),
        output:
            asv_table       = OUT / "06.taxonomy" / "asv_table.txt",
            asv_table_seqs  = OUT / "06.taxonomy" / "asv_table_seqs.txt",
            taxon_seq_table = OUT / "06.taxonomy" / "taxon_seq_table.txt",
            # Every ASV that reached the filter, with metaxa2_origin and filter_status.
            asv_table_unfiltered = OUT / "06.taxonomy" / "asv_table_unfiltered.txt",
        params:
            amp_type        = AMPLICON_TYPE,
            sintax_cutoff   = _tax_cfg["sintax_cutoff"],
            # sintax rank model — differs from the rdp one for markers whose two
            # reference files disagree on depth (see MARKERS in 00_common.smk).
            tax_levels      = SINTAX_TAX_LEVELS,
            rank_letters    = SINTAX_RANK_LETTERS,
            rank_prefixes   = SINTAX_RANK_PREFIXES,
            filter_enabled  = FILTER_ENABLED,
            filter_keep     = FILTER_KEEP,
            filter_discard  = FILTER_DISCARD,
            discard_min_support = FILTER_DISCARD_MIN_SUPPORT,   # best-hit check; None = off
            organelle_filter = ORGANELLE_FILTER,   # remove Metaxa2's mitochondria/chloroplast calls
            seed            = config["amplicon"]["seed"],
        log:
            LOGS / "assign_taxonomy.log",
        conda:
            "../../envs/vsearch.yaml"
        threads: lambda wc: threads_for("assign_taxonomy")
        resources:
            mem_mb = lambda wc: mem_mb_for("assign_taxonomy"),
        script:
            "../../scripts/amplicon/80c_parse_sintax.py"


# Joins every per-stage read-count source produced earlier in the run (the
# preprocessing .counts files, DADA2's own read.counts, the length-filtered
# seqtab, and the final asv_table) into one wide table: one row per sample, one
# column per stage, in pipeline order (raw → ... → post_taxonomy_filter — see
# 80b_aggregate_read_counts.py for the exact column list, which varies with
# remove_phix / extraction_enabled). post_taxonomy_filter is summed from the
# filtered asv_table.txt, so it counts every removal of the taxonomy filter: keep,
# discard, the best-hit check and the organelle filter. This is the run's
# end-of-line diagnostic: read_tracking.txt is requested directly by rule all and
# meant for a human to scan for where a sample lost most of its reads; no other
# rule reads it back.
rule aggregate_read_counts:
    input:
        raw_counts       = OUT / "stats" / "1.raw_reads.counts",
        stripped_counts  = OUT / "stats" / "3.stripped_reads.counts",
        dada2_counts     = OUT / "05.dada2" / "read.counts",
        seqtab_lenfilt   = OUT / "05.dada2" / "seqtab_lenfilt_head_names.txt",
        # seqtab_extracted only present when extraction is enabled
        seqtab_extracted = [OUT / "05.dada2" / "seqtab_extracted_head_names.txt"] if EXTRACTION_ENABLED else [],
        asv_table        = OUT / "06.taxonomy" / "asv_table.txt",
        # nophix_counts only present when PhiX removal runs
        **({"nophix_counts": OUT / "stats" / "2.nophix_reads.counts"} if REMOVE_PHIX else {}),
    output:
        counts = OUT / "stats" / "read_tracking.txt",
    params:
        extraction_enabled = EXTRACTION_ENABLED,
        remove_phix        = REMOVE_PHIX,
        samples            = SAMPLES,
    log:
        LOGS / "aggregate_read_counts.log",
    # Standard library only; the env pins the interpreter instead of using the host's.
    conda:
        "../../envs/python_utils.yaml"
    script:
        "../../scripts/amplicon/80b_aggregate_read_counts.py"
