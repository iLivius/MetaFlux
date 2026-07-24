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
# Metaxa2 : -t all (all domains). All extracted sequences pass through; domain-
#           level filtering happens post-taxonomy via the contaminant filter
#           in assign_taxonomy.R.
# ITSx    : -t all, --only_full F (mandatory for single-region ITS1/ITS2 amplicons).
#           The region file matching amplicon.its_region is kept.

if MARKER_EXTRACTOR == "metaxa2":
    rule target_extract:
        input:
            seqs         = OUT / "5.dada2" / "seqs.fasta",
            seqtab_names = OUT / "5.dada2" / "seqtab_head_names.txt",
        output:
            seqs         = OUT / "5.dada2" / "seqs_extracted.fasta",
            seqtab_names = OUT / "5.dada2" / "seqtab_extracted_head_names.txt",
            results      = OUT / "5.dada2" / "metaxa2_extraction.results.txt",
        params:
            # All raw Metaxa2 outputs (seqs.*) are persisted under this dir for
            # downstream QC / inspection; see 70a_metaxa2_extract.py.
            prefix = lambda wc: str(OUT / "5.dada2" / "metaxa2" / "seqs"),
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
    rule target_extract:
        input:
            seqs         = OUT / "5.dada2" / "seqs.fasta",
            seqtab_names = OUT / "5.dada2" / "seqtab_head_names.txt",
        output:
            seqs         = OUT / "5.dada2" / "seqs_extracted.fasta",
            seqtab_names = OUT / "5.dada2" / "seqtab_extracted_head_names.txt",
            results      = OUT / "5.dada2" / "itsx_extraction.summary.txt",
            collapse_map = OUT / "5.dada2" / "itsx_collapse_map.tsv",
        params:
            its_region = ITS_REGION,
            # Raw ITSx outputs (seqs.ITS1/ITS2.fasta, seqs.positions.txt, seqs.summary.txt,
            # ...) are persisted under this dir for QC/inspection, same as Metaxa2's above —
            # the leading underscore is only a naming convention, nothing deletes this dir.
            prefix     = lambda wc: str(OUT / "5.dada2" / "_itsx_tmp" / "seqs"),
        log:
            LOGS / "target_extract.log",
        conda:
            "../../envs/itsx.yaml"
        threads: lambda wc: threads_for("target_extract")
        resources:
            mem_mb = lambda wc: mem_mb_for("target_extract"),
        script:
            "../../scripts/amplicon/70b_itsx_extract.py"
