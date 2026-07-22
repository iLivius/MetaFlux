# Taxonomy assignment and final output table generation.
#
# assign_taxonomy       : two implementations selected at parse time by
#                         amplicon.taxonomy.method (rdp | sintax):
#                           rdp    — DADA2 assignTaxonomy (RDP classifier) + addSpecies (16S only)
#                           sintax — vsearch --sintax (memory-efficient; ~2–3 GB regardless of ASV count)
#                         Both paths apply the same contaminant filter and write
#                         the same three output tables.
# aggregate_read_counts : joins all per-stage read-count sources into one wide
#                         TSV (one row per sample, one column per stage).

_tax_cfg = config["amplicon"]["taxonomy"]

if TAXONOMY_METHOD == "rdp":
    rule assign_taxonomy:
        input:
            seqs         = _seqs_for_taxonomy,
            seqtab_names = _seqtab_for_taxonomy,
            refdb        = TAXONOMY_REFDB,
            species_db   = TAXONOMY_SPECIES_DB,    # addSpecies for 16S only; empty list for ITS
        output:
            asv_table       = OUT / "6.taxonomy" / "asv_table.txt",
            asv_table_seqs  = OUT / "6.taxonomy" / "asv_table_seqs.txt",
            taxon_seq_table = OUT / "6.taxonomy" / "taxon_seq_table.txt",
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
    rule assign_taxonomy:
        input:
            seqs         = _seqs_for_taxonomy,
            seqtab_names = _seqtab_for_taxonomy,
            refdb        = TAXONOMY_SINTAX_DB,
        output:
            asv_table       = OUT / "6.taxonomy" / "asv_table.txt",
            asv_table_seqs  = OUT / "6.taxonomy" / "asv_table_seqs.txt",
            taxon_seq_table = OUT / "6.taxonomy" / "taxon_seq_table.txt",
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


rule aggregate_read_counts:
    input:
        raw_counts       = OUT / "stats" / "1.raw_reads.counts",
        stripped_counts  = OUT / "stats" / "3.stripped_reads.counts",
        dada2_counts     = OUT / "5.dada2" / "read.counts",
        seqtab_lenfilt   = OUT / "5.dada2" / "seqtab_lenfilt_head_names.txt",
        # seqtab_extracted only present when extraction is enabled
        seqtab_extracted = [OUT / "5.dada2" / "seqtab_extracted_head_names.txt"] if EXTRACTION_ENABLED else [],
        asv_table        = OUT / "6.taxonomy" / "asv_table.txt",
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
    script:
        "../../scripts/amplicon/80b_aggregate_read_counts.py"
