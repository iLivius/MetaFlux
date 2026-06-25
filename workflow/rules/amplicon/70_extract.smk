# Target extraction rules: Metaxa2 (16S) or ITSx (ITS).
# Activated only when amplicon.extraction.enabled is true.
# Both rules are named target_extract; only one is defined depending on
# AMPLICON_TYPE (resolved at parse time), so the DAG always sees a single rule.
#
# Metaxa2 : -t all (all domains). All extracted sequences pass through; domain-
#           level filtering happens post-taxonomy via the contaminant filter
#           in assign_taxonomy.R.
# ITSx    : -t all, --only_full F (mandatory for single-region ITS1/ITS2 amplicons).
#           The region file matching amplicon.its_region is kept.

if AMPLICON_TYPE == "16S":
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

else:  # ITS
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
