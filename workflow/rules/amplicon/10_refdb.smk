# Reference DB management (amplicon-only). Each rule runs only when its output is
# missing — DBs are cached under refdb/.
#
# REFERENCE PREPARATION MAP — every reference falls into one of three prep tiers:
#
#   TIER 0  plain download, ready for its tool (NO code)
#           silva_species, unite_sintax, pr2_dada2, pr2_utax, silva_euk,
#           gyrb_dada2, plus silva_train (SILVA 138.2) or silva_train_raw
#           (SILVA 144)  → handled by the GENERATED fetch loop at the bottom of
#           this file, straight from the marker pack's url, md5-checked when the
#           pack lists a checksum.
#
#   TIER 1  download + unpack an archive (bespoke rule, extraction only)
#           unite_fasta (.tgz)         → fetch_unite
#           unite_uchime_ITS1/2 (.zip) → fetch_uchime
#           rpob_frogs (.tar.gz)       → fetch_rpob_archive  (raw FASTA, still FROGS-shaped)
#
#   TIER 2  reshape the FASTA into the format a classifier wants (bespoke rule +
#           a named script in workflow/scripts/refdb/ — that folder is the ONE place
#           the database-specific transform logic lives)
#           silva_train + silva_kingdom_map
#                         ← convert_silva_train     (refdb/silva_drop_kingdom.py; SILVA 144 only)
#           silva_sintax  ← convert_silva_sintax   (refdb/silva_to_sintax.py)
#           rpob_frogs    ← convert_rpob_to_dada2   (refdb/frogs_rpob_to_dada2.py)
#           pr2_sintax    ← convert_pr2_sintax      (refdb/pr2_utax_domain_first.py)
#
# PhiX is separate: it is a decontamination reference, not taxonomy — fetched and
# bowtie2-indexed below.

# Downloads the PhiX genome — a sequencing-control genome Illumina runs routinely
# spike into the flow cell — from the NCBI URL in the config. This raw FASTA feeds
# build_phix_index below; rm_phix (30_preprocess.smk) then aligns every sample's
# reads against that index to strip out PhiX contamination before primer trimming.
rule fetch_phix:
    output:
        fasta = REFDB_PHIX_FASTA,
    params:
        url = config["references"]["phix"]["fetch_url"],
    log:
        LOGS / "refdb" / "fetch_phix.log",
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -O {output.fasta}.gz {params.url} > {log} 2>&1
        gunzip -f {output.fasta}.gz >> {log} 2>&1
        """


# Builds the bowtie2 index from the PhiX FASTA fetched above (a one-time cost —
# the index is reused by every sample and every future run against this reference).
# Consumed by rm_phix (30_preprocess.smk), which aligns each sample's reads
# against it and keeps only the reads that do NOT map, i.e. the non-PhiX reads.
rule build_phix_index:
    input:
        fasta = REFDB_PHIX_FASTA,
    output:
        idx = multiext(
            str(REFDB_PHIX_PREFIX),
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
            ".rev.1.bt2", ".rev.2.bt2",
        ),
    params:
        prefix = str(REFDB_PHIX_PREFIX),
    log:
        LOGS / "refdb" / "build_phix_index.log",
    conda:
        "../../envs/bowtie2.yaml"
    threads: lambda wc: threads_for("bowtie2")
    shell:
        """
        bowtie2-build {input.fasta} {params.prefix} > {log} 2>&1
        """


# ───────────────────── UNITE general-release fetch (ITS) ───────────────
# Taxonomy reference ONLY (taxonomy_refdb in ITS.yaml) — never the probe. The
# probe uses the separate, pre-extracted UCHIME ITS1/ITS2 subregion instead
# (fetch_uchime, below), since ITS's probe_mode is "direct": it needs sequences
# already cut to amplicon length, not the full ITS region this file provides.
rule fetch_unite:
    output:
        fasta = UNITE_FASTA,
    params:
        url = PACK_REF_URLS["unite_fasta"],
    log:
        LOGS / "refdb" / "fetch_unite.log",
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.fasta})
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        wget -O $tmp/unite.tgz {params.url} > {log} 2>&1
        tar -xzf $tmp/unite.tgz -C $tmp >> {log} 2>&1
        # Pick the "dynamic_s_all" FASTA from the archive (UNITE convention)
        src=$(find $tmp -type f -name 'sh_general_release_dynamic_s_all_*.fasta' | head -1)
        if [ -z "$src" ]; then
            echo "Could not find sh_general_release_dynamic_s_all_*.fasta in UNITE archive" >> {log}
            exit 1
        fi
        mv "$src" {output.fasta}
        """


# ───────────────────── SILVA 144: drop the Kingdom rank (16S) ──────────────
# SILVA 144 added a Kingdom rank (Bacillati, Pseudomonadati, ...) between Domain
# and Phylum, so its lineages have 7 ranks where MetaFlux's 16S rank model (and
# SILVA 138.2) has 6 before the species slot. Fed as-is, every name would land one
# rank too deep: genus in the species slot, Chloroplast at family, Mitochondria at
# genus — and the o__Chloroplast / f__Mitochondria filter tokens would quietly stop
# matching anything.
#
# This rule reads the file exactly as SILVA publishes it (silva_train_raw) and
# writes two files next to each other:
#   - the 6-rank training set (Domain;Phylum;...;Genus), sequences untouched —
#     the file assignTaxonomy, the amplicon probe and convert_silva_sintax use;
#   - a phylum -> kingdom table built from the same headers, which rule
#     restore_kingdom (80_taxonomy.smk) uses to put the Kingdom back next to the
#     final taxonomy.
# amplicon.taxonomy.silva_other: strip also removes SILVA's "--other" suffix here,
# before any classifier sees the names (the output file names carry the setting).
#
# Defined only when the active SILVA release ships a Kingdom (144). Under 138.2 the
# trainset is a plain download and this rule does not exist, so two rules can never
# claim the same output file.
if SILVA_DROPS_KINGDOM:

    rule convert_silva_train:
        input:
            fasta = SILVA_TRAIN_RAW,
        output:
            fasta       = SILVA_TRAIN,
            kingdom_map = SILVA_KINGDOM_MAP,
        params:
            strip_other = SILVA_OTHER == "strip",
        log:
            LOGS / "refdb" / "convert_silva_train.log",
        conda:
            "../../envs/python_utils.yaml"
        script:
            "../../scripts/refdb/silva_drop_kingdom.py"


# ───────────────────── SILVA SINTAX conversion (16S) ──────────────────────
# Reformats the SILVA RDP trainset (the 6-rank one MetaFlux classifies with) into
# VSEARCH SINTAX header format (>id;tax=d:...,p:...,g:...). No download — depends on
# fetch_silva_train (SILVA 138.2) or convert_silva_train (SILVA 144). Cached at the
# silva_sintax path the marker pack gives for the active release.
rule convert_silva_sintax:
    input:
        fasta = SILVA_TRAIN,
    output:
        fasta = SILVA_SINTAX,
    log:
        LOGS / "refdb" / "convert_silva_sintax.log",
    conda:
        "../../envs/python_utils.yaml"
    script:
        "../../scripts/refdb/silva_to_sintax.py"


# ───────────────────── UNITE UCHIME ITS1/ITS2 subregions (ITS) ─────────────
# The UNITE UCHIME reference dataset (.zip): ITS1 and ITS2 already cut to their
# amplicon-length subregions. MetaFlux uses these as the ITS PROBE substrate
# (ITS.yaml probe_ref: unite_uchime_{region}, direct mode) — lengths are read
# straight off them, no in-silico PCR. This is NOT the SINTAX taxonomy DB: that
# is unite_sintax, a plain .gz download handled by the generated fetch loop at
# the bottom of this file.
rule fetch_uchime:
    output:
        its1 = UNITE_UCHIME_ITS1_FA,
        its2 = UNITE_UCHIME_ITS2_FA,
    params:
        url = PACK_REF_URLS["unite_uchime_ITS1"],
    log:
        LOGS / "refdb" / "fetch_uchime.log",
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.its1})
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        wget -O $tmp/unite_uchime.zip {params.url} > {log} 2>&1
        unzip -q $tmp/unite_uchime.zip -d $tmp >> {log} 2>&1
        src1=$(find $tmp -type f -name '*_ITS1.fasta' | head -1)
        src2=$(find $tmp -type f -name '*_ITS2.fasta' | head -1)
        if [ -z "$src1" ] || [ -z "$src2" ]; then
            echo "Could not find ITS1/ITS2 FASTA files in UNITE UCHIME archive" >> {log}
            exit 1
        fi
        mv "$src1" {output.its1}
        mv "$src2" {output.its2}
        """


# ───────────────────── rpoB FROGS → DADA2 (fetch + convert) ────────────────
# The FROGS RefSeq rpoB release is a .tar.gz bundling the FASTA, a BLAST index and
# an RDP-classifier tree, with FROGS-shaped headers:
#   >WP_095092576.1 Root;k__Bacteria [id: 1];p__Bacillota [id: 2];...;s__X [id: 7]
# Prepared in two steps that mirror the SILVA fetch+convert pattern above:
#   fetch_rpob_archive   — download + unpack, keep only the raw FASTA (Tier 1)
#   convert_rpob_to_dada2 — rewrite headers to a DADA2 trainset (Tier 2), the
#                           reshape logic living in refdb/frogs_rpob_to_dada2.py

# Raw FROGS FASTA (still FROGS-shaped); consumed by the converter, then deleted.
RPOB_FROGS_RAW = RPOB_FROGS.parent / "rpob_frogs_raw.fasta.gz"

rule fetch_rpob_archive:
    output:
        fasta = temp(RPOB_FROGS_RAW),
    params:
        url = PACK_REF_URLS["rpob_frogs"],
    log:
        LOGS / "refdb" / "fetch_rpob_archive.log",
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.fasta})
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        wget -O $tmp/rpob.tar.gz "{params.url}" > {log} 2>&1
        tar -xzf $tmp/rpob.tar.gz -C $tmp >> {log} 2>&1
        # The single true FASTA in the archive (the BLAST index files end in
        # .fasta.n**, and there is a .fasta.properties — exclude both).
        src=$(find $tmp -type f -name '*.fasta' ! -name '*.properties' | head -1)
        if [ -z "$src" ]; then
            echo "Could not find the rpoB FASTA in the FROGS archive" >> {log}
            exit 1
        fi
        gzip -c "$src" > {output.fasta}
        """


rule convert_rpob_to_dada2:
    input:
        fasta = RPOB_FROGS_RAW,
    output:
        fasta = RPOB_FROGS,
    log:
        LOGS / "refdb" / "convert_rpob_to_dada2.log",
    conda:
        "../../envs/python_utils.yaml"
    script:
        "../../scripts/refdb/frogs_rpob_to_dada2.py"


# ── 18S references ────────────────────────────────────────────────────────
# SILVA-Euk is the probe substrate (in-silico primer recovery 96-97% vs 66-75%
# on PR2); PR2 is the taxonomy reference. PR2 ships two files with DIFFERENT
# rank depths — the DADA2 one (9 ranks) feeds the rdp path, the UTAX one
# (8 ranks) feeds sintax after the conversion below — so both are fetched
# independently.
#
# Three downloads plus one derived file. None of the three downloads gets a
# bespoke rule here: they are plain files (no `archive` key in 18S.yaml), so they
# fall through to the generated fetch loop at the bottom of this file (Tier 0),
# which builds fetch_silva_euk, fetch_pr2_dada2 and fetch_pr2_utax automatically.
# The derived file, pr2_sintax, is made from pr2_utax by rule convert_pr2_sintax
# right below (Tier 2).

# Rewrites PR2's UTAX file into the SINTAX file the 18S sintax path classifies
# with. PR2 puts the domain under the SINTAX letter k: and the supergroup under
# d:, but VSEARCH reads the letters as a fixed top-down order (d, k, p, c, ...)
# and counts a rank's confidence only over the bootstrap rounds that agreed on
# every rank above it. On the file as published the supergroup was therefore the
# top rank, and an ASV lost its domain whenever its supergroup was uncertain. The
# script swaps the two letters, headers only (sequences untouched), and stops on
# any file without PR2 5.1.1's letter order or domain values, so a file can never
# be swapped twice.
#   in : PR2_UTAX   (pr2_utax, downloaded by the generated rule fetch_pr2_utax)
#   out: PR2_SINTAX (pr2_sintax), read by rule assign_taxonomy (80_taxonomy.smk)
#        on the sintax path. Built only when an 18S run uses method: sintax; the
#        rdp path reads PR2's DADA2 file and never asks for it.
rule convert_pr2_sintax:
    input:
        fasta = PR2_UTAX,
    output:
        fasta = PR2_SINTAX,
    log:
        LOGS / "refdb" / "convert_pr2_sintax.log",
    conda:
        "../../envs/python_utils.yaml"
    script:
        "../../scripts/refdb/pr2_utax_domain_first.py"


# ── Generated reference fetches ───────────────────────────────────────────
# One rule per plain-download reference declared by the ACTIVE marker pack
# (workflow/markers/<type>.yaml → references.<symbol>.url), so adding such a
# reference — or a whole new marker — is a data change with no rule to write
# here. References needing archive extraction (UNITE .tgz/.zip, PhiX .gz) or a
# transform (silva_sintax, derived from silva_train) keep the bespoke rules above.
#
# When the pack lists an md5 for the file, the download is checked against it
# before anything uses it: a truncated download, a changed file behind the same URL
# or a proxy error page would otherwise be trained on without complaint. On a
# mismatch the file is deleted and the rule fails, saying which checksum it got.
# References without an md5 entry are downloaded exactly as before.
for _sym, _spec in PACK_FETCHABLE.items():

    rule:
        name:   f"fetch_{_sym}"
        output: _spec["path"]
        params:
            url = _spec["url"],
            md5 = _spec["md5"],
        log:
            LOGS / "refdb" / f"fetch_{_sym}.log",
        conda:
            "../../envs/bowtie2.yaml"
        shell:
            """
            mkdir -p $(dirname {output})
            wget -O {output} "{params.url}" > {log} 2>&1
            if [ -n "{params.md5}" ]; then
                got=$(md5sum {output} | cut -d' ' -f1)
                if [ "$got" != "{params.md5}" ]; then
                    echo "md5 mismatch for {output}: expected {params.md5}, got $got" >> {log}
                    rm -f {output}
                    exit 1
                fi
                echo "md5 OK: {params.md5}" >> {log}
            fi
            """
