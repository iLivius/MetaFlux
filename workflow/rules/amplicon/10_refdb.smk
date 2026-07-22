# Reference DB management (amplicon-only): PhiX fetch + bowtie2 index;
# SILVA training-set + species (16S); UNITE general release + UCHIME (ITS).
# Each rule only runs when its output is missing — DBs are cached under refdb/.

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


# ───────────────────── SILVA training-set fetch (16S) ──────────────────
# Triggered on demand: amplicon_probe (when expected_length=auto) and the
# DADA2 taxonomy step. Cached at the path declared in config.
rule fetch_silva_train:
    output:
        fasta = SILVA_TRAIN,
    params:
        url = config["references"]["silva"]["fetch_url_train"],
    log:
        LOGS / "refdb" / "fetch_silva_train.log",
    conda:
        "../../envs/bowtie2.yaml"                  # provides wget
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -O {output.fasta} {params.url} > {log} 2>&1
        """


rule fetch_silva_species:
    output:
        fasta = SILVA_SPECIES,
    params:
        url = config["references"]["silva"]["fetch_url_species"],
    log:
        LOGS / "refdb" / "fetch_silva_species.log",
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -O {output.fasta} {params.url} > {log} 2>&1
        """


# ───────────────────── UNITE fetch (ITS) ───────────────────────────────
# Two separate downloads:
#   fetch_unite  : full UNITE general release (used by assignTaxonomy)
#   fetch_uchime : UNITE UCHIME release with pre-extracted ITS1/ITS2 sequences
#                  (used for the auto length-probe distribution)
rule fetch_unite:
    output:
        fasta = UNITE_FASTA,
    params:
        url = config["references"]["unite"]["fetch_url"],
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


# ───────────────────── SILVA SINTAX conversion (16S) ──────────────────────
# Reformats the already-downloaded SILVA RDP trainset into VSEARCH SINTAX
# header format (>id;tax=d:...,p:...,g:...). No download — depends on
# fetch_silva_train. Cached at refdb/silva/silva_nr99_v138.2_sintax.fa.gz.
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
        "../../scripts/amplicon/10a_silva_to_sintax.py"


# ───────────────────── UNITE SINTAX download (ITS) ─────────────────────────
# Official UNITE+INSD VSEARCH/SINTAX release (February 2025, UNITE v10).
# Already in SINTAX header format — no conversion needed.
rule fetch_unite_sintax:
    output:
        fasta = UNITE_SINTAX,
    params:
        url = config["references"]["unite"]["fetch_url_sintax"],
    log:
        LOGS / "refdb" / "fetch_unite_sintax.log",
    conda:
        "../../envs/bowtie2.yaml"                  # provides wget
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -O {output.fasta} {params.url} > {log} 2>&1
        """


rule fetch_uchime:
    output:
        its1 = UNITE_UCHIME_ITS1_FA,
        its2 = UNITE_UCHIME_ITS2_FA,
    params:
        url = config["references"]["unite"]["fetch_url_uchime"],
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


# ── 18S references ────────────────────────────────────────────────────────
# SILVA-Euk is the probe substrate (in-silico primer recovery 96-97% vs 66-75%
# on PR2); PR2 is the taxonomy reference. PR2 ships two files with DIFFERENT
# rank depths — the DADA2 one (9 ranks) feeds the rdp path, the UTAX one
# (8 ranks) feeds sintax — so both are fetched independently.

rule fetch_silva_euk:
    output:
        fasta = SILVA_EUK,
    params:
        url = config["references"].get("silva_euk", {}).get(
            "fetch_url",
            "https://zenodo.org/records/1447330/files/silva_132.18s.99_rep_set.dada2.fa.gz?download=1",
        ),
    log:
        LOGS / "refdb" / "fetch_silva_euk.log",
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -O {output.fasta} "{params.url}" > {log} 2>&1
        """


rule fetch_pr2_dada2:
    output:
        fasta = PR2_DADA2,
    params:
        url = config["references"].get("pr2", {}).get(
            "fetch_url_dada2",
            "https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_dada2.fasta.gz",
        ),
    log:
        LOGS / "refdb" / "fetch_pr2_dada2.log",
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -O {output.fasta} "{params.url}" > {log} 2>&1
        """


rule fetch_pr2_utax:
    output:
        fasta = PR2_UTAX,
    params:
        url = config["references"].get("pr2", {}).get(
            "fetch_url_utax",
            "https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_UTAX.fasta.gz",
        ),
    log:
        LOGS / "refdb" / "fetch_pr2_utax.log",
    conda:
        "../../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -O {output.fasta} "{params.url}" > {log} 2>&1
        """
