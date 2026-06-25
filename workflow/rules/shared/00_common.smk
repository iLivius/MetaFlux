# Shared utilities: mode dispatch, path resolution, sample discovery,
# per-rule resource lookups, and mode-specific globals.
#
# This module is included first in the Snakefile and resolves everything
# that downstream rule modules consume at parse time.

import csv
import hashlib
import json
import os
import sys
import urllib.parse
import urllib.request
from pathlib import Path

from snakemake.io import glob_wildcards


# ────────────────────────── Mode dispatch ──────────────────────
MODE = config["mode"].lower()
if MODE not in ("amplicon", "shotgun"):
    sys.exit(f"config.mode must be 'amplicon' or 'shotgun' (got: {config['mode']!r})")


# ────────────────────────── Path resolution ────────────────────
# All filesystem paths are absolute Path objects, resolved at parse time.
# No `workdir:` directive — shell rules execute from the invocation dir.
FASTQ_DIR = Path(config["input"]["fastq_dir"]).resolve()
OUT       = Path(config["output"]["out_dir"]).resolve()
LOGS      = OUT / "logs"
BENCH     = OUT / "benchmarks"

# PhiX is shared (both modes); amplicon also uses the index basename.
REFDB_PHIX_FASTA  = Path(config["references"]["phix"]["fasta"]).resolve()
REFDB_PHIX_PREFIX = REFDB_PHIX_FASTA.with_suffix("")   # strip .fna → bowtie2 index basename

# Project-internal cache that survives across runs and out_dirs.
PROJECT_DIR     = Path(workflow.basedir).parent.resolve()
PROBE_CACHE_DIR = PROJECT_DIR / "refdb" / "cache"

# Safe in both modes: BioProject mode needs fastq_dir to be writable.
FASTQ_DIR.mkdir(parents=True, exist_ok=True)


# ────────────────────── Resource lookups ───────────────────────
# Single source of truth for per-rule threads / mem_mb is config["resources"].
# Rules reference these via lambda wc: threads_for("<rule>"), etc.
def threads_for(rule_name: str) -> int:
    cfg = config.get("resources", {})
    return int(cfg.get("threads", {}).get(rule_name, cfg.get("threads_default", 1)))


def mem_mb_for(rule_name: str) -> int:
    cfg = config.get("resources", {})
    return int(cfg.get("mem_mb", {}).get(rule_name, cfg.get("mem_mb_default", 1000)))


# ────────────────────── Amplicon-only globals ──────────────────
if MODE == "amplicon":
    amp_cfg = config["amplicon"]

    # ── References ──
    SILVA_TRAIN          = Path(config["references"]["silva"]["train"]).resolve()
    SILVA_SPECIES        = Path(config["references"]["silva"]["species"]).resolve()
    UNITE_FASTA          = Path(config["references"]["unite"]["fasta"]).resolve()
    UNITE_UCHIME_ITS1_FA = Path(config["references"]["unite"]["uchime_its1"]).resolve()
    UNITE_UCHIME_ITS2_FA = Path(config["references"]["unite"]["uchime_its2"]).resolve()

    # ── Primers (hard error if missing — mode=amplicon requires them) ──
    PRIMER_FWD = Path(amp_cfg["primers"]["fwd"]).resolve()
    PRIMER_REV = Path(amp_cfg["primers"]["rev"]).resolve()
    for _p in (PRIMER_FWD, PRIMER_REV):
        if not _p.exists():
            sys.exit(f"[MetaFlux] mode=amplicon requires primer FASTA but not found: {_p}")

    AMPLICON_TYPE = amp_cfg["type"].upper()
    if AMPLICON_TYPE not in ("16S", "ITS"):
        sys.exit(f"amplicon.type must be '16S' or 'ITS' (got: {amp_cfg['type']!r})")

    ITS_REGION = amp_cfg.get("its_region", "ITS2").upper()
    if AMPLICON_TYPE == "ITS" and ITS_REGION not in ("ITS1", "ITS2"):
        sys.exit(f"amplicon.its_region must be 'ITS1' or 'ITS2' (got: {ITS_REGION!r})")

    ORIENTATION = amp_cfg["primers"].get("orientation", "fixed").lower()
    if ORIENTATION not in ("fixed", "mixed"):
        sys.exit(f"amplicon.primers.orientation must be 'fixed' or 'mixed' (got: {ORIENTATION!r})")

    # ── PhiX decontamination toggle (symmetric with shotgun) ──
    # When False, bowtie2 PhiX removal is skipped and Cutadapt reads the raw
    # FASTQs directly. The 'nophix' falco stage / counts / read-tracking column
    # exist only when this is True.
    REMOVE_PHIX  = bool(amp_cfg.get("decontamination", {}).get("remove_phix", True))
    FALCO_STAGES = ["raw", "nophix", "stripped"] if REMOVE_PHIX else ["raw", "stripped"]

    # ── Amplicon-probe cache ──
    # Cache key = sha256(fwd_primer + "|" + rev_primer)[:12]. Changing primers
    # invalidates the cache automatically; same primers across runs hit the same key.
    def _primer_pair_hash(fwd: Path, rev: Path) -> str:
        h = hashlib.sha256()
        h.update(fwd.read_bytes())
        h.update(b"|")
        h.update(rev.read_bytes())
        return h.hexdigest()[:12]

    PRIMER_HASH = _primer_pair_hash(PRIMER_FWD, PRIMER_REV)

    if AMPLICON_TYPE == "16S":
        # 16S: in-silico PCR against SILVA using two-pass cutadapt
        PROBE_REF_FASTA = SILVA_TRAIN
        PROBE_REF_TAG   = "silva_v138.2_toGenus"
        PROBE_MODE      = "pcr"
    else:
        # ITS: direct length measurement from UNITE UCHIME pre-extracted subregion
        # sequences; no primer binding sites needed, no cutadapt involved.
        PROBE_REF_FASTA = UNITE_UCHIME_ITS1_FA if ITS_REGION == "ITS1" else UNITE_UCHIME_ITS2_FA
        PROBE_REF_TAG   = f"unite_uchime_{ITS_REGION}"
        PROBE_MODE      = "direct"

    PROBE_JSON         = PROBE_CACHE_DIR / f"probe_{AMPLICON_TYPE}_{PROBE_REF_TAG}_{PRIMER_HASH}.json"
    PROBE_AMPLICONS_FA = PROBE_CACHE_DIR / f"probe_{AMPLICON_TYPE}_{PROBE_REF_TAG}_{PRIMER_HASH}.amplicons.fa.gz"

    def _wants_probe(cfg) -> bool:
        v = cfg.get("expected_length")
        return isinstance(v, str) and v.strip().lower() == "auto"

    WANTS_PROBE = _wants_probe(amp_cfg)

    # ── Extraction / taxonomy globals ──
    EXTRACTION_ENABLED  = amp_cfg["extraction"].get("enabled", True)

    # RDP path (method: rdp)
    TAXONOMY_REFDB      = SILVA_TRAIN if AMPLICON_TYPE == "16S" else UNITE_FASTA
    TAXONOMY_SPECIES_DB = [str(SILVA_SPECIES)] if AMPLICON_TYPE == "16S" else []

    # SINTAX path (method: sintax)
    SILVA_SINTAX        = Path(config["references"]["silva"]["sintax"]).resolve()
    UNITE_SINTAX        = Path(config["references"]["unite"]["sintax"]).resolve()
    TAXONOMY_SINTAX_DB  = SILVA_SINTAX if AMPLICON_TYPE == "16S" else UNITE_SINTAX

    # Taxonomy method — validated at parse time
    TAXONOMY_METHOD = amp_cfg["taxonomy"].get("method", "rdp").lower()
    if TAXONOMY_METHOD not in ("rdp", "sintax"):
        sys.exit(
            f"amplicon.taxonomy.method must be 'rdp' or 'sintax' "
            f"(got: {amp_cfg['taxonomy']['method']!r})"
        )

    # Taxonomy always reads the length-filtered outputs (dada_length_filter
    # runs after target_extract — see amplicon/60_dada2.smk).
    def _seqs_for_taxonomy(wildcards=None):
        return OUT / "5.dada2" / "seqs_lenfilt.fasta"

    def _seqtab_for_taxonomy(wildcards=None):
        return OUT / "5.dada2" / "seqtab_lenfilt_head_names.txt"


# ────────────────────── Shotgun-only globals ───────────────────
else:  # MODE == "shotgun"
    sht_cfg = config["shotgun"]

    # ── Reference DB (user-supplied; not auto-fetched) ──
    KRAKEN_DB = Path(config["references"]["kraken_db"]).resolve()
    if not KRAKEN_DB.exists():
        sys.exit(f"[MetaFlux] mode=shotgun requires kraken_db but not found: {KRAKEN_DB}")

    # ── Primer warning (not an error — amplicon section can be filled out
    # and ignored when running shotgun). ──
    _amp_primers = config.get("amplicon", {}).get("primers", {})
    if _amp_primers.get("fwd") or _amp_primers.get("rev"):
        sys.stderr.write(
            "[MetaFlux] warning: mode=shotgun but amplicon.primers fields are set — they will be ignored\n"
        )

    # ── Decontamination toggles ──
    _decontam   = sht_cfg.get("decontamination", {}) or {}
    REMOVE_PHIX  = bool(_decontam.get("remove_phix", True))
    HOST_GENOMES = list(_decontam.get("host_genomes", []) or [])
    HOST_MIN_ID  = float(_decontam.get("host_min_id", 0.95))   # BBMap minid for host removal

    # ── fastp ──
    MIN_READ_LEN = int(sht_cfg["fastp"]["min_read_length"])

    # ── Kraken2 / Bracken ──
    KRAKEN_CONF      = float(sht_cfg["kraken"]["confidence"])
    KRAKEN_HIT       = int(sht_cfg["kraken"]["hit_groups"])
    KRAKEN_MMAP      = bool(sht_cfg["kraken"].get("memory_mapping", True))
    KRAKEN2_THREADS  = sht_cfg["kraken"].get("threads") or threads_for("kraken2")
    TAX_LEV          = sht_cfg["bracken"]["tax_lev"]
    BRACKEN_THRESH   = int(sht_cfg["bracken"]["threshold"])
    EXTRACT_TAXA     = list(sht_cfg.get("extract_taxa", []) or [])
    DL_THREADS       = int(sht_cfg.get("download_threads", 4))

    # Kraken2 mem_mb is derived from the mmap mode + actual DB size at parse time:
    #   mmap on   → ~20 GB private workspace per process (DB is shared via the OS
    #               file cache, so it doesn't count against this rule's budget).
    #   mmap off  → hash.k2d size + 10% headroom for buffers / read I/O.
    #               Self-adjusts to the actual DB (16 GB pluspf_16gb, ~100 GB pluspf,
    #               200+ GB core_nt, etc).
    def _kraken_hash_size_mb(db_path: Path):
        h = db_path / "hash.k2d"
        return (h.stat().st_size // (1024 * 1024)) if h.exists() else None

    _hash_mb = _kraken_hash_size_mb(KRAKEN_DB)
    if KRAKEN_MMAP:
        KRAKEN2_MEM_MB = 20000
    elif _hash_mb:
        KRAKEN2_MEM_MB = int(_hash_mb * 1.1)
    else:
        # only reached if the DB path doesn't exist yet — conservative ~100 GB default
        KRAKEN2_MEM_MB = 110000

    # ── BioProject SRA fetch ──
    _inp           = config.get("input", {}) or {}
    BIOPROJECT     = _inp.get("bioproject") or None
    ACCESSION_LIST = _inp.get("accession_list") or None


# ────────────────────── Sample discovery ───────────────────────
ALLOWED_EXTENSIONS = ("fastq", "fq", "fastq.gz", "fq.gz")
BAD_CHARS = set("*#@%^/! ?&:;|<>")             # underscore allowed — SRR ids don't use it


def _fetch_bioproject_runs(bioproject: str, cache_path: Path):
    """Resolve PAIRED SRA runs for a BioProject via NCBI e-utils.
    Caches the runinfo CSV so subsequent invocations are offline."""
    if not cache_path.exists():
        sys.stderr.write(f"[MetaFlux] resolving {bioproject} via NCBI...\n")
        esearch = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=sra&term={urllib.parse.quote(bioproject)}[BioProject]"
            "&retmax=10000&retmode=json"
        )
        with urllib.request.urlopen(esearch, timeout=60) as r:
            ids = json.load(r).get("esearchresult", {}).get("idlist", [])
        if not ids:
            sys.exit(f"[MetaFlux] no SRA records for {bioproject}")
        efetch = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=sra&id={','.join(ids)}&rettype=runinfo&retmode=text"
        )
        with urllib.request.urlopen(efetch, timeout=120) as r:
            cache_path.write_text(r.read().decode())

    # Only keep PAIRED runs — single-end metagenomes don't fit the rest of the pipeline.
    runs = []
    with open(cache_path) as f:
        for row in csv.DictReader(f):
            if row.get("Run") and row.get("LibraryLayout", "").upper() == "PAIRED":
                runs.append(row["Run"])
    return sorted(set(runs))


def _discover_local(fastq_dir: Path):
    """Glob paired FASTQs. Prefers {sample}_R1/_R2 (amplicon + our default),
    falls back to {sample}_1/_2 (raw fasterq-dump convention).
    Returns (samples, extensions, naming) where naming is '_R1' or '_1'."""
    sR, eR = glob_wildcards(str(fastq_dir / "{sample}_R1.{extn}"))
    if sR:
        return list(sR), list(eR), "_R1"
    s1, e1 = glob_wildcards(str(fastq_dir / "{sample}_1.{extn}"))
    return list(s1), list(e1), "_1"


# Resolve SAMPLES + extension + naming at parse time.
if MODE == "shotgun" and BIOPROJECT:
    _cache = FASTQ_DIR / f"{BIOPROJECT}_runinfo.csv"
    SAMPLES = _fetch_bioproject_runs(BIOPROJECT, _cache)
    if ACCESSION_LIST:
        SAMPLES = [s for s in SAMPLES if s in set(ACCESSION_LIST)]
    if not SAMPLES:
        sys.exit("[MetaFlux] no samples after filtering accession_list")
    EXTN         = "fastq.gz"                  # download_sra always produces .gz
    FASTQ_NAMING = "_R1"                       # download_sra renames SRA's _1/_2 to _R1/_R2
    INPUT_MODE   = "bioproject"
else:
    _samples, _exts, FASTQ_NAMING = _discover_local(FASTQ_DIR)
    if not _samples:
        sys.exit(f"[MetaFlux] no *{FASTQ_NAMING}.* FASTQs found in {FASTQ_DIR}")
    if len(set(_exts)) != 1 or _exts[0] not in ALLOWED_EXTENSIONS:
        sys.exit(f"[MetaFlux] mixed or unsupported FASTQ extensions: {set(_exts)}")
    SAMPLES    = sorted(_samples)
    EXTN       = _exts[0]
    INPUT_MODE = "local"

    if MODE == "amplicon" and FASTQ_NAMING != "_R1":
        sys.exit(
            f"[MetaFlux] mode=amplicon requires {{sample}}_R1/_R2 naming; "
            f"found {{sample}}_1/_2 in {FASTQ_DIR}"
        )

# Sample-name sanity check — snakemake gets unhappy when wildcards contain
# special characters; fail fast with a clear message.
for _s in SAMPLES:
    if any(c in _s for c in BAD_CHARS):
        sys.exit(f"[MetaFlux] sample '{_s}' has disallowed characters")

sys.stderr.write(
    f"[MetaFlux v1.0.0] mode={MODE}, input={INPUT_MODE}, "
    f"{len(SAMPLES)} sample(s) in {FASTQ_DIR}\n"
)


# ────────────────────── Path helpers (shared) ──────────────────
def raw_fastq(sample: str, read: int) -> str:
    """Resolve a sample's raw FASTQ path for read 1 or 2.
    Honors the run-level naming convention (_R1/_R2 or _1/_2) discovered above,
    and prefers .gz over plain .fastq."""
    if read == 1:
        mate = FASTQ_NAMING
    else:
        mate = "_R2" if FASTQ_NAMING == "_R1" else "_2"
    gz    = FASTQ_DIR / f"{sample}{mate}.fastq.gz"
    plain = FASTQ_DIR / f"{sample}{mate}.fastq"
    # Pre-download (BioProject) the .gz doesn't exist yet — return the .gz path
    # so download_sra's output declaration aligns with downstream consumers.
    if gz.exists():
        return str(gz)
    if plain.exists():
        return str(plain)
    return str(gz)


def raw_r1(sample: str) -> str:
    return raw_fastq(sample, 1)


def raw_r2(sample: str) -> str:
    return raw_fastq(sample, 2)


# ────────────────── Shotgun-only chain helpers ─────────────────
# Preprocessing chain: raw → (phix) → (host) → fastp. Each function picks the
# correct upstream stage so the consumer rule doesn't have to know what ran before.
if MODE == "shotgun":
    def fastp_upstream(wildcards):
        s = wildcards.sample
        if HOST_GENOMES:
            return {
                "r1": str(OUT / "01.preprocessing" / f"{s}_dehost_R1.fastq.gz"),
                "r2": str(OUT / "01.preprocessing" / f"{s}_dehost_R2.fastq.gz"),
            }
        if REMOVE_PHIX:
            return {
                "r1": str(OUT / "01.preprocessing" / f"{s}_dephix_R1.fastq.gz"),
                "r2": str(OUT / "01.preprocessing" / f"{s}_dephix_R2.fastq.gz"),
            }
        return {"r1": raw_r1(s), "r2": raw_r2(s)}

    def host_upstream(wildcards):
        s = wildcards.sample
        if REMOVE_PHIX:
            return {
                "r1": str(OUT / "01.preprocessing" / f"{s}_dephix_R1.fastq.gz"),
                "r2": str(OUT / "01.preprocessing" / f"{s}_dephix_R2.fastq.gz"),
            }
        return {"r1": raw_r1(s), "r2": raw_r2(s)}


# ────────────────── Amplicon-only chain helper ─────────────────
# Cutadapt reads from the PhiX-filtered stage when PhiX removal is on, otherwise
# straight from the raw FASTQs — so trim_primers doesn't hard-code the upstream.
if MODE == "amplicon":
    def primer_trim_upstream(wildcards):
        s = wildcards.sample
        if REMOVE_PHIX:
            return {
                "r1": str(OUT / "2.no_phix" / f"{s}_R1.fastq.gz"),
                "r2": str(OUT / "2.no_phix" / f"{s}_R2.fastq.gz"),
            }
        return {"r1": raw_fastq(s, 1), "r2": raw_fastq(s, 2)}


# ────────────────────── rule all targets ───────────────────────
def _amplicon_targets():
    """Default target list for mode=amplicon. Mirrors DADAism v3 rule-all."""
    targets = [
        # Raw FASTQ symlinks (1.reads/ — visibility only)
        *expand(OUT / "1.reads" / "{sample}_R{r}.fastq.gz", sample=SAMPLES, r=[1, 2]),
    ]
    # PhiX-removed FASTQs and their read count exist only when PhiX removal runs.
    if REMOVE_PHIX:
        targets += [
            *expand(OUT / "2.no_phix" / "{sample}_R{r}.fastq.gz", sample=SAMPLES, r=[1, 2]),
            OUT / "stats" / "2.nophix_reads.counts",
        ]
    targets += [
        # Preprocessing outputs
        *expand(OUT / "3.stripped" / "{sample}_R{r}.fastq.gz", sample=SAMPLES, r=[1, 2]),
        # Read count summaries
        OUT / "stats" / "1.raw_reads.counts",
        OUT / "stats" / "3.stripped_reads.counts",
        # Falco QC (per sample, per direction, per stage — 'nophix' only if REMOVE_PHIX)
        *expand(OUT / "stats" / "falco" / "{sample}_R{r}_{stage}",
                stage=FALCO_STAGES, sample=SAMPLES, r=[1, 2]),
        # TruncLen decision
        OUT / "stats" / "trunclen.json",
        # DADA2 quality diagnostic plots
        OUT / "stats" / "dada2" / "stripped_read_R1_qual_plot.png",
        OUT / "stats" / "dada2" / "stripped_read_R2_qual_plot.png",
        OUT / "stats" / "dada2" / "stripped_read_R1_qual_plot.pdf",
        OUT / "stats" / "dada2" / "stripped_read_R2_qual_plot.pdf",
        # DADA2 filterAndTrim
        *expand(OUT / "4.filtered" / "{sample}_R1_filt.fastq.gz", sample=SAMPLES),
        *expand(OUT / "4.filtered" / "{sample}_R2_filt.fastq.gz", sample=SAMPLES),
        # DADA2 core
        OUT / "5.dada2" / "seqs.fasta",
        OUT / "5.dada2" / "seqtab_head_seqs.txt",
        OUT / "5.dada2" / "seqtab_head_names.txt",
        OUT / "5.dada2" / "read.counts",
        OUT / "stats" / "dada2" / "filtered_read_R1_error_plot.png",
        OUT / "stats" / "dada2" / "filtered_read_R2_error_plot.png",
        OUT / "stats" / "dada2" / "filtered_read_R1_error_plot.pdf",
        OUT / "stats" / "dada2" / "filtered_read_R2_error_plot.pdf",
        # ASV length distribution stats + filter
        OUT / "stats" / "dada2" / "asv_length_stats.json",
        OUT / "stats" / "dada2" / "asv_length_hist.png",
        # Taxonomy
        OUT / "6.taxonomy" / "asv_table.txt",
        OUT / "6.taxonomy" / "asv_table_seqs.txt",
        OUT / "6.taxonomy" / "taxon_seq_table.txt",
        # Combined read-tracking table
        OUT / "stats" / "read_tracking.txt",
        # MultiQC aggregate
        OUT / "multiqc" / "multiqc_report.html",
    ]
    return targets


def _shotgun_targets():
    """Default target list for mode=shotgun. Mirrors MetaFlux v0.1.0 rule-all."""
    out = [
        *expand(OUT / "01.preprocessing" / "{s}_fastp.json",          s=SAMPLES),
        *expand(OUT / "02.classification" / "{s}_report.txt",         s=SAMPLES),
        *expand(OUT / "02.classification" / "{s}_output.txt",         s=SAMPLES),
        *expand(OUT / "02.classification" / "{s}_R1.fastq.gz",        s=SAMPLES),
        *expand(OUT / "02.classification" / "{s}_R2.fastq.gz",        s=SAMPLES),
        *expand(OUT / "03.abundance" / "{s}_report.txt",              s=SAMPLES),
        *expand(OUT / "03.abundance" / "{s}_output.txt",              s=SAMPLES),
        OUT / "03.abundance" / "otu_table.tsv",
        OUT / "stats" / "read_tracking.txt",
        OUT / "multiqc" / "multiqc_report.html",
        OUT / "samples.tsv",
    ]
    # Per-taxon extracted reads — only added when extract_taxa is non-empty.
    for taxon in EXTRACT_TAXA:
        # Slugify for clean wildcards; map back to the scientific name at runtime.
        slug = taxon.lower().replace(" ", "_")
        out.extend(expand(OUT / "04.extracted_reads" / f"{{s}}_{slug}_R1.fastq.gz", s=SAMPLES))
        out.extend(expand(OUT / "04.extracted_reads" / f"{{s}}_{slug}_R2.fastq.gz", s=SAMPLES))
    return out


def all_targets():
    return _amplicon_targets() if MODE == "amplicon" else _shotgun_targets()
