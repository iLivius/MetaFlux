# Target extraction rules: Metaxa2 (16S) or ITSx (ITS).
# Activated only when amplicon.extraction.enabled is true.
# Both rules are named target_extract; only one is defined depending on the
# marker profile's extractor (resolved at parse time), so the DAG always sees a
# single rule. A marker whose extractor is "none" (18S, gyrB, rpoB) defines no
# rule here and must run with extraction disabled; 00_common.smk's guard forces
# EXTRACTION_ENABLED off (with a warning) if such a marker's config still has
# extraction.enabled: true, so a "none" marker never reaches this file expecting
# a rule it doesn't get. The two branches below stay exhaustive for the same
# reason: a "none" marker takes neither.
#
# Metaxa2 : -t all (all domains). All extracted sequences pass through here. Metaxa2's
#           own origin call per ASV (bacteria, archaea, eukaryota, chloroplast,
#           mitochondria, uncertain) is written to metaxa2_origin.tsv; the taxonomy
#           filter (80a/80c) uses its organelle calls to remove mitochondria and
#           chloroplasts (amplicon.taxonomy.filter.organelles), while the keep/discard
#           tokens on the classified lineage still decide everything else.
# ITSx    : -t all, --only_full F (mandatory for single-region ITS1/ITS2 amplicons).
#           The region file matching amplicon.its_region is kept.

if MARKER_EXTRACTOR == "metaxa2":
    # Runs Metaxa2 on the full ASV set from dada_seqtab (seqs.fasta +
    # seqtab_head_names.txt, both from 60_dada2.smk). A 16S amplicon sits
    # entirely inside the 16S gene — every standard primer binds within it — so
    # there is little or nothing for Metaxa2 to actually trim; here it mainly
    # acts as a filter, dropping ASVs it cannot recognise as SSU rRNA at all
    # (off-target amplification). That is why -t all (every domain Metaxa2
    # recognises, not just bacteria/archaea) keeps every extracted sequence
    # regardless of domain: nothing is pruned here. Writes the extracted FASTA
    # and a seqtab subset to the ASVs Metaxa2 could resolve, both consumed next
    # by dada_length_filter (60_dada2.smk). `results` is Metaxa2's own summary
    # report — useful for inspection, but not read by any other rule.
    #
    # `origin` is Metaxa2's own origin call for every extracted ASV, one row per
    # ASV, collected from its per-origin FASTAs (seqs.mitochondria.fasta, ...).
    # Rule assign_taxonomy reads it: with amplicon.taxonomy.filter.organelles:
    # metaxa2 (the default) the ASVs called mitochondria or chloroplast are
    # removed there, next to the keep/discard tokens, and every ASV's call is
    # written to the metaxa2_origin column of 06.taxonomy/asv_table_unfiltered.txt.
    # Removing them on the classified lineage alone is not enough: on the 16S test
    # set, plant mitochondria whose nearest SILVA reference is a mitochondrion at
    # 97-99 % identity were reported as bacteria (mostly Enterobacterales) and
    # escaped every token.
    rule target_extract:
        input:
            seqs         = OUT / "05.dada2" / "seqs.fasta",
            seqtab_names = OUT / "05.dada2" / "seqtab_head_names.txt",
        output:
            seqs         = OUT / "05.dada2" / "seqs_extracted.fasta",
            seqtab_names = OUT / "05.dada2" / "seqtab_extracted_head_names.txt",
            results      = OUT / "05.dada2" / "metaxa2_extraction.results.txt",
            origin       = METAXA2_ORIGIN_TSV,       # 05.dada2/metaxa2_origin.tsv
        params:
            # All raw Metaxa2 outputs (seqs.*) are persisted under this dir for
            # downstream QC / inspection; see 70a_metaxa2_extract.py.
            prefix = lambda wc: str(OUT / "05.dada2" / "metaxa2" / "seqs"),
        log:
            LOGS / "target_extract.log",
        conda:
            "../../envs/metaxa2.yaml"
        threads: lambda wc: threads_for("target_extract")
        resources:
            mem_mb = lambda wc: mem_mb_for("target_extract"),
        script:
            "../../scripts/amplicon/70a_metaxa2_extract.py"

elif MARKER_EXTRACTOR == "itsx":
    # Runs ITSx on the full ASV set from dada_seqtab to cut away the flanking
    # rRNA fragments (5.8S, SSU, LSU) an ITS ASV carries outside the marker
    # region proper, keeping only the sub-region named by amplicon.its_region
    # (ITS1 or ITS2). --only_full F keeps partial detections, which is
    # mandatory here: a targeted ITS1/ITS2 amplicon never contains the FULL ITS
    # region (flanked by both SSU and LSU) that ITSx's "full" detection mode
    # expects. -t all scans across all of ITSx's eukaryote lineage profiles.
    # Side effect of cutting the flanks away: distinct DADA2 ASVs that only
    # differed in those flanking bases can come out with an IDENTICAL extracted
    # sub-region. The script collapses each such group into one representative
    # ASV and sums their per-sample counts (collapse_map, below) — otherwise
    # the same organism would be double-counted under two ASV IDs. Writes the
    # deduplicated extracted FASTA and summed seqtab, both consumed next by
    # dada_length_filter (60_dada2.smk); `results` is ITSx's own summary (QC
    # only); `collapse_map` is the audit trail of which original ASV IDs were
    # merged into which representative — not read by any other rule, but
    # useful for checking why the ASV count dropped after extraction.
    rule target_extract:
        input:
            seqs         = OUT / "05.dada2" / "seqs.fasta",
            seqtab_names = OUT / "05.dada2" / "seqtab_head_names.txt",
        output:
            seqs         = OUT / "05.dada2" / "seqs_extracted.fasta",
            seqtab_names = OUT / "05.dada2" / "seqtab_extracted_head_names.txt",
            results      = OUT / "05.dada2" / "itsx_extraction.summary.txt",
            collapse_map = OUT / "05.dada2" / "itsx_collapse_map.tsv",
        params:
            its_region = ITS_REGION,
            # Raw ITSx outputs (seqs.ITS1/ITS2.fasta, seqs.positions.txt, seqs.summary.txt,
            # ...) are persisted under this dir for QC/inspection, same as Metaxa2's above —
            # the leading underscore is only a naming convention, nothing deletes this dir.
            prefix     = lambda wc: str(OUT / "05.dada2" / "_itsx_tmp" / "seqs"),
        log:
            LOGS / "target_extract.log",
        conda:
            "../../envs/itsx.yaml"
        threads: lambda wc: threads_for("target_extract")
        resources:
            mem_mb = lambda wc: mem_mb_for("target_extract"),
        script:
            "../../scripts/amplicon/70b_itsx_extract.py"
