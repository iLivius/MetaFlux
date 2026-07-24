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
import yaml
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


# ────────────────────── Marker registry ────────────────────────
# One profile per amplicon marker. Adding a marker means adding an entry here
# (plus its reference-fetch rule in amplicon/10_refdb.smk) — every choice that
# used to be a scattered `if AMPLICON_TYPE == "16S" else ...` now lives in one
# place, so a new marker can no longer silently inherit the ITS branch.
#
#   probe_mode          "pcr"    — in-silico PCR of the primers against a
#                                  full-length reference (16S/SILVA).
#                       "direct" — the reference is already the trimmed
#                                  subregion; measure its lengths (ITS/UNITE).
#   probe_ref           reference symbol to probe (resolved via REF_PATHS below);
#                       "{region}" is filled from amplicon.its_region.
#   probe_ref_tag       label embedded in the probe-cache filename.
#   probe_stat_key      which amplicon.probe_length_stat.<key> to read (only the
#                       pcr path uses it; the direct path ignores the value).
#   taxonomy_refdb      reference symbol for the DADA2/RDP training set.
#   taxonomy_species_db reference symbol for DADA2 addSpecies, or None.
#   taxonomy_sintax_db  reference symbol for the VSEARCH SINTAX database.
#   extractor           "metaxa2" | "itsx" | "none" — target-region extractor.
#   tax_levels          ordered rank column names emitted by the classifier
#                       (assignTaxonomy columns / SINTAX ranks).
#   rank_letters        the SINTAX one-letter rank codes paired to tax_levels.
#   rank_prefixes       the output prefixes (k__/p__/…) paired to tax_levels.
#   sintax_tax_levels / sintax_rank_letters / sintax_rank_prefixes
#                       OPTIONAL per-method overrides. The rdp and sintax
#                       references usually share a rank model (16S, ITS), but not
#                       always: PR2 ships 9 ranks in its DADA2 file and only 8 in
#                       its UTAX file, where Division+Subdivision are merged into
#                       one 'p:' field. When absent these fall back to the fields
#                       above. Prefixes are keyed to rank MEANING, so the same
#                       keep/discard token works under either method.
#   prefix_style        "bare"     — classifier values carry no rank prefix, so
#                                    the builder adds rank_prefixes (SILVA/16S).
#                                    The species slot becomes a Genus+species
#                                    binomial.
#                       "embedded" — values already carry k__/p__ prefixes, so
#                                    the builder emits them verbatim (UNITE/ITS).
#                       (rdp path only; the SINTAX path always normalises to
#                       bare letters then re-adds rank_prefixes.)
#
# The contaminant keep/discard lists are NOT a pack field: they live only in
# config (amplicon.taxonomy.filter.keep / .discard), documented per marker in
# the config template + README — there is no hidden default to reconcile.
#
# These fields are DATA, in workflow/markers/<type>.yaml (see _load_marker_packs).
# Adding a marker means adding a pack file, not editing this module.
MARKERS_SHIPPED_DIR = Path(workflow.basedir) / "markers"
MARKERS_USER_DIR    = Path("config") / "markers"


def _load_marker_packs():
    """Load marker packs from YAML into their OWN namespace.

    Packs are DATA, never merged into `config`: a pack is read as MARKERS[type][...]
    and a run setting as config["amplicon"][...], so no key has two homes and there
    is no precedence rule to remember. A pack therefore cannot silently supply a
    user preference, and a run config cannot accidentally shadow a marker fact.

    workflow/markers/*.yaml ships with the pipeline; config/markers/*.yaml is the
    user's own, shadowing a shipped pack of the same name. Adding a marker is a
    data-only change — no code edit.
    """
    packs = {}
    for d in (MARKERS_SHIPPED_DIR, MARKERS_USER_DIR):
        if not d.is_dir():
            continue
        for f in sorted(d.glob("*.yaml")):
            try:
                with open(f) as fh:
                    pack = yaml.safe_load(fh) or {}
            except Exception as e:
                sys.exit(f"[MetaFlux] could not read marker pack {f}: {e}")
            if not isinstance(pack, dict):
                sys.exit(f"[MetaFlux] marker pack {f} must be a YAML mapping")
            packs[f.stem] = pack
    if not packs:
        sys.exit(
            f"[MetaFlux] no marker packs found in {MARKERS_SHIPPED_DIR} "
            f"(or {MARKERS_USER_DIR}). The installation looks incomplete."
        )
    return packs


MARKERS = _load_marker_packs()

def _as_taxon_list(v):
    """Normalize a keep/discard config value to a list of tokens.
    YAML lets `keep: [k__Fungi]` and `keep: k__Fungi` both parse; the bare scalar
    arrives as a str and MUST NOT reach the scripts unwrapped (Python's list("k__Fungi")
    would char-split it, silently emptying the table). null/empty → []."""
    if v is None:
        return []
    if isinstance(v, str):
        return [v]
    return list(v)


# ────────────────────── Amplicon-only globals ──────────────────
if MODE == "amplicon":
    amp_cfg = config["amplicon"]

    # ── References (resolved paths + symbol table) ──
    # A path comes from the run config when given (so existing configs are
    # unchanged), otherwise from the owning marker pack's references.<symbol>.file
    # relative to references.refdb_root. Packs therefore carry every marker's DB
    # locations and URLs, and the config need only name a reference to RELOCATE it.
    REFDB_ROOT = Path(config["references"].get("refdb_root", "refdb"))

    _PACK_REF_FILES = {}
    PACK_REF_URLS = {}
    for _pack in MARKERS.values():
        for _sym, _spec in (_pack.get("references") or {}).items():
            if not isinstance(_spec, dict):
                continue
            if _spec.get("file"):
                _PACK_REF_FILES.setdefault(_sym, _spec["file"])
            if _spec.get("url"):
                PACK_REF_URLS.setdefault(_sym, _spec["url"])

    def _ref_path(symbol, *cfg_keys):
        """Resolve a reference path: explicit config entry first, else pack default."""
        node = config["references"]
        for k in cfg_keys:
            node = node.get(k) if isinstance(node, dict) else None
        if node:
            return Path(node).resolve()
        if symbol in _PACK_REF_FILES:
            return (REFDB_ROOT / _PACK_REF_FILES[symbol]).resolve()
        sys.exit(
            f"[MetaFlux] reference '{symbol}' has no path. Either set references."
            f"{'.'.join(cfg_keys)} in the config, or give the marker pack a "
            f"references.{symbol}.file entry."
        )

    SILVA_TRAIN          = _ref_path("silva_train",       "silva", "train")
    SILVA_SPECIES        = _ref_path("silva_species",     "silva", "species")
    SILVA_SINTAX         = _ref_path("silva_sintax",      "silva", "sintax")
    UNITE_FASTA          = _ref_path("unite_fasta",       "unite", "fasta")
    UNITE_SINTAX         = _ref_path("unite_sintax",      "unite", "sintax")
    UNITE_UCHIME_ITS1_FA = _ref_path("unite_uchime_ITS1", "unite", "uchime_its1")
    UNITE_UCHIME_ITS2_FA = _ref_path("unite_uchime_ITS2", "unite", "uchime_its2")
    PR2_DADA2            = _ref_path("pr2_dada2",         "pr2",   "dada2")
    PR2_UTAX             = _ref_path("pr2_utax",          "pr2",   "utax")
    SILVA_EUK            = _ref_path("silva_euk",         "silva_euk", "fasta")
    GYRB_DADA2           = _ref_path("gyrb_dada2",        "gyrb",  "dada2")
    RPOB_FROGS           = _ref_path("rpob_frogs",        "rpob",  "frogs")

    # Marker profiles reference DBs by symbol; resolve symbol → Path here.
    REF_PATHS = {
        "silva_train":       SILVA_TRAIN,
        "silva_species":     SILVA_SPECIES,
        "silva_sintax":      SILVA_SINTAX,
        "unite_fasta":       UNITE_FASTA,
        "unite_sintax":      UNITE_SINTAX,
        "unite_uchime_ITS1": UNITE_UCHIME_ITS1_FA,
        "unite_uchime_ITS2": UNITE_UCHIME_ITS2_FA,
        "pr2_dada2":         PR2_DADA2,
        "pr2_utax":          PR2_UTAX,
        "silva_euk":         SILVA_EUK,
        "gyrb_dada2":        GYRB_DADA2,
        "rpob_frogs":        RPOB_FROGS,
    }

    # ── Primers (hard error if missing — mode=amplicon requires them) ──
    PRIMER_FWD = Path(amp_cfg["primers"]["fwd"]).resolve()
    PRIMER_REV = Path(amp_cfg["primers"]["rev"]).resolve()
    for _p in (PRIMER_FWD, PRIMER_REV):
        if not _p.exists():
            sys.exit(f"[MetaFlux] mode=amplicon requires primer FASTA but not found: {_p}")

    # Match amplicon.type to a marker pack case-INSENSITIVELY, then adopt the pack's
    # own name as canonical. 16S/ITS/18S are uppercase so this is a no-op for them;
    # mixed-case gene markers (gyrB, rpoB) keep their conventional casing rather than
    # being flattened to GYRB/RPOB in logs and cache filenames. A user may write
    # gyrb / GYRB / gyrB and all resolve to the shipped pack 'gyrB'.
    _marker_by_ci = {name.upper(): name for name in MARKERS}
    _req_type = amp_cfg["type"]
    if _req_type.upper() not in _marker_by_ci:
        sys.exit(f"amplicon.type must be one of {sorted(MARKERS)} (got: {_req_type!r})")
    AMPLICON_TYPE = _marker_by_ci[_req_type.upper()]

    ITS_REGION = amp_cfg.get("its_region", "ITS2").upper()
    if AMPLICON_TYPE == "ITS" and ITS_REGION not in ("ITS1", "ITS2"):
        sys.exit(f"amplicon.its_region must be 'ITS1' or 'ITS2' (got: {ITS_REGION!r})")

    # The active marker's profile — the single source for every dispatch below.
    PROFILE = MARKERS[AMPLICON_TYPE]

    # Symbols the ACTIVE pack can fetch with a plain download. 10_refdb.smk turns
    # each into a generated fetch rule. Entries carrying an `archive` (UNITE .tgz /
    # .zip) need extraction and keep their bespoke rule there instead, as does
    # silva_sintax, which is derived from silva_train rather than downloaded.
    PACK_FETCHABLE = {
        _sym: {"path": REF_PATHS[_sym], "url": _spec["url"]}
        for _sym, _spec in (PROFILE.get("references") or {}).items()
        if isinstance(_spec, dict)
        and _spec.get("url")
        and not _spec.get("archive")
        and _sym in REF_PATHS
    }

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

    # ── Amplicon-probe wiring (from the marker profile) ──
    #   pcr mode    (16S, 18S, rpoB) → in-silico PCR of the primers against a
    #               full-length reference (SILVA / SILVA-Euk / FROGS).
    #   direct mode (ITS, gyrB)      → lengths read straight from a reference that
    #               is already amplicon-length (UNITE UCHIME / DD7RZ8).
    PROBE_MODE        = PROFILE["probe_mode"]
    if PROBE_MODE not in ("pcr", "direct"):
        sys.exit(
            f"[MetaFlux] marker pack {AMPLICON_TYPE!r} has an invalid probe_mode "
            f"{PROBE_MODE!r} — must be 'pcr' or 'direct'."
        )
    PROBE_REF_FASTA   = REF_PATHS[PROFILE["probe_ref"].format(region=ITS_REGION)]
    PROBE_REF_TAG     = PROFILE["probe_ref_tag"].format(region=ITS_REGION)
    # Every marker except ITS can reach resolve_expected_length() in pick_trunclen.py
    # (ITS returns early and never resolves a stat) — so gate this validation on
    # amp_type, not probe_mode. A direct-mode marker like gyrB still needs a valid
    # probe_length_stat entry if it is ever run with expected_length: auto, even
    # though its own probe never runs an in-silico PCR.
    _stat_key = PROFILE["probe_stat_key"]
    _stat_cfg = config["amplicon"].get("probe_length_stat", {}) or {}
    if AMPLICON_TYPE != "ITS" and _stat_key not in _stat_cfg:
        sys.exit(
            f"[MetaFlux] amplicon.probe_length_stat is missing the '{_stat_key}' key that "
            f"marker {AMPLICON_TYPE} needs. Add it, e.g.\n"
            f"    probe_length_stat:\n"
            f"      {_stat_key}: p95"
        )
    PROBE_LENGTH_STAT = _stat_cfg.get(_stat_key)

    PROBE_JSON         = PROBE_CACHE_DIR / f"probe_{AMPLICON_TYPE}_{PROBE_REF_TAG}_{PRIMER_HASH}.json"
    PROBE_AMPLICONS_FA = PROBE_CACHE_DIR / f"probe_{AMPLICON_TYPE}_{PROBE_REF_TAG}_{PRIMER_HASH}.amplicons.fa.gz"

    def _wants_probe(cfg) -> bool:
        v = cfg.get("expected_length")
        return isinstance(v, str) and v.strip().lower() == "auto"

    WANTS_PROBE = _wants_probe(amp_cfg)

    # ITS's length_filter never resolves expected_length through pick_trunclen.py
    # (its ITS branch always writes expected_length: null to trunclen.json, since
    # truncLen is fixed to c(0,0) for ITS regardless of amplicon length — see
    # 50a_pick_trunclen.py). dada_length_filter's mode=auto has two ways to size its
    # window: the probe JSON (only produced when WANTS_PROBE, i.e. expected_length:
    # 'auto'), or a fallback that reads expected_length out of trunclen.json. For ITS
    # with a manual expected_length, WANTS_PROBE is False (no probe JSON) and the
    # fallback hits that permanent null, crashing on int(None) — AFTER the full DADA2
    # run. Catch the combination here instead, before anything expensive runs.
    if AMPLICON_TYPE == "ITS" and amp_cfg["length_filter"]["mode"] == "auto" and not WANTS_PROBE:
        sys.exit(
            "[MetaFlux] ITS with amplicon.length_filter.mode: auto needs "
            "amplicon.expected_length: 'auto' too, so a probe-based length window can be "
            "built — a manual (non-'auto') expected_length gives the ITS length filter "
            "nothing to size its window from. Either set expected_length: 'auto', or set "
            "length_filter.mode: manual with an explicit length_filter.range."
        )

    # ── Extraction / taxonomy globals (from the marker profile) ──
    EXTRACTION_ENABLED  = amp_cfg["extraction"].get("enabled", True)
    MARKER_EXTRACTOR    = PROFILE["extractor"]           # metaxa2 | itsx | none
    if MARKER_EXTRACTOR not in ("metaxa2", "itsx", "none"):
        sys.exit(
            f"[MetaFlux] marker pack {AMPLICON_TYPE!r} has an invalid extractor "
            f"{MARKER_EXTRACTOR!r} — must be 'metaxa2', 'itsx', or 'none'."
        )

    # A marker whose profile has no region extractor (18S, gyrB, rpoB) defines no
    # target_extract rule in 70_extract.smk. Leaving extraction.enabled true would
    # make the DAG ask for seqs_extracted.fasta, which nothing produces — a
    # confusing missing-input error. Force it off and say so instead.
    if MARKER_EXTRACTOR == "none" and EXTRACTION_ENABLED:
        sys.stderr.write(
            f"[MetaFlux] warning: marker {AMPLICON_TYPE} has no target-region extractor; "
            "ignoring amplicon.extraction.enabled: true (no extraction step will run)\n"
        )
        EXTRACTION_ENABLED = False

    # RDP path (method: rdp)
    TAXONOMY_REFDB      = REF_PATHS[PROFILE["taxonomy_refdb"]]
    TAXONOMY_SPECIES_DB = (
        [str(REF_PATHS[PROFILE["taxonomy_species_db"]])]
        if PROFILE["taxonomy_species_db"] else []
    )

    # SINTAX path (method: sintax). Some markers have no SINTAX build at all
    # (gyrB's DD7RZ8, rpoB's FROGS release ship a DADA2/rdp trainset only), so the
    # pack leaves taxonomy_sintax_db null and this stays None; the method guard
    # below refuses method: sintax for such a marker rather than failing obscurely.
    _sintax_sym         = PROFILE["taxonomy_sintax_db"]
    TAXONOMY_SINTAX_DB  = REF_PATHS[_sintax_sym] if _sintax_sym else None

    # ── Taxonomy rank model (from the marker profile) ──
    # Drives the taxonomy-string builder in both 80a (rdp) and 80c (sintax),
    # replacing the previously hardcoded 7-rank tuples. 16S/ITS/gyrB/rpoB use the
    # Linnaean 7; 18S/PR2 is the one 9-rank marker today.
    # rdp path rank model
    TAX_LEVELS    = PROFILE["tax_levels"]
    RANK_LETTERS  = PROFILE["rank_letters"]
    RANK_PREFIXES = PROFILE["rank_prefixes"]
    PREFIX_STYLE  = PROFILE["prefix_style"]

    # sintax path rank model — identical to the rdp one unless the marker declares
    # otherwise (18S does: PR2's UTAX file has 8 ranks to the DADA2 file's 9).
    SINTAX_TAX_LEVELS    = PROFILE.get("sintax_tax_levels",    TAX_LEVELS)
    SINTAX_RANK_LETTERS  = PROFILE.get("sintax_rank_letters",  RANK_LETTERS)
    SINTAX_RANK_PREFIXES = PROFILE.get("sintax_rank_prefixes", RANK_PREFIXES)

    # These three (and the sintax_* triad) are positionally matched, hand-maintained
    # lists in the marker pack — nothing else in the pipeline checks they line up.
    # A miscounted pack silently corrupts taxonomy strings instead of erroring: R's
    # out-of-range vector indexing returns NA (80a renders "NAvalue"), and Python's
    # zip() silently truncates to the shortest list (80c drops trailing ranks) —
    # neither crashes. Catch the malformed pack here instead, at parse time.
    if not (len(TAX_LEVELS) == len(RANK_LETTERS) == len(RANK_PREFIXES)):
        sys.exit(
            f"[MetaFlux] marker pack {AMPLICON_TYPE!r} has mismatched rank-model list "
            f"lengths: tax_levels={len(TAX_LEVELS)}, rank_letters={len(RANK_LETTERS)}, "
            f"rank_prefixes={len(RANK_PREFIXES)} — all three must be the same length."
        )
    if not (len(SINTAX_TAX_LEVELS) == len(SINTAX_RANK_LETTERS) == len(SINTAX_RANK_PREFIXES)):
        sys.exit(
            f"[MetaFlux] marker pack {AMPLICON_TYPE!r} has mismatched sintax rank-model "
            f"list lengths: sintax_tax_levels={len(SINTAX_TAX_LEVELS)}, "
            f"sintax_rank_letters={len(SINTAX_RANK_LETTERS)}, "
            f"sintax_rank_prefixes={len(SINTAX_RANK_PREFIXES)} — all three must be the same length."
        )

    # ── Contaminant keep/discard (config-only; no per-marker code default) ──
    # Lists of rank-prefixed tokens (e.g. k__Bacteria, o__Chloroplast) matched
    # against whole taxonomy-string segments. Empty/null list → that direction
    # is a no-op. The recommended per-marker values ship in the config template.
    _filter_cfg      = amp_cfg["taxonomy"]["filter"]
    # Fail loud on the pre-step-3 schema instead of silently ignoring it: the old
    # include_pattern/exclude_pattern regex keys used to carry per-type defaults;
    # a reused config that still has them (and no keep/discard) would otherwise run
    # with the filter silently disabled.
    if "include_pattern" in _filter_cfg or "exclude_pattern" in _filter_cfg:
        sys.exit(
            "[MetaFlux] amplicon.taxonomy.filter uses the retired include_pattern/"
            "exclude_pattern regex schema. Replace them with keep/discard token lists, e.g.\n"
            "    keep:    [k__Bacteria, k__Archaea]\n"
            "    discard: [o__Chloroplast, f__Mitochondria]   # 16S\n"
            "  or keep: [k__Fungi]   discard: []              # ITS"
        )
    FILTER_ENABLED   = bool(_filter_cfg.get("enabled", True))
    FILTER_KEEP      = _as_taxon_list(_filter_cfg.get("keep"))
    FILTER_DISCARD   = _as_taxon_list(_filter_cfg.get("discard"))

    # Taxonomy method — validated at parse time
    TAXONOMY_METHOD = amp_cfg["taxonomy"].get("method", "rdp").lower()
    if TAXONOMY_METHOD not in ("rdp", "sintax"):
        sys.exit(
            f"amplicon.taxonomy.method must be 'rdp' or 'sintax' "
            f"(got: {amp_cfg['taxonomy']['method']!r})"
        )
    # A marker with no SINTAX reference (gyrB, rpoB) can only be classified by rdp.
    # Catch method: sintax here with a clear message instead of a None-path crash
    # deeper in the sintax rule.
    if TAXONOMY_METHOD == "sintax" and TAXONOMY_SINTAX_DB is None:
        sys.exit(
            f"[MetaFlux] marker {AMPLICON_TYPE} has no SINTAX reference database, so "
            f"amplicon.taxonomy.method: sintax is not available for it. Use method: rdp."
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

    # ── OTU-table taxon filter (shotgun counterpart of amplicon keep/discard) ──
    # Same rank-aware token matching; applied to the final otu_table. Config-only,
    # no code default — default OFF (see config comment).
    _sht_filter        = sht_cfg.get("taxonomy_filter", {}) or {}
    OTU_FILTER_ENABLED = bool(_sht_filter.get("enabled", False))
    OTU_FILTER_KEEP    = _as_taxon_list(_sht_filter.get("keep"))
    OTU_FILTER_DISCARD = _as_taxon_list(_sht_filter.get("discard"))

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
        # reached when hash.k2d isn't present under an existing kraken_db dir yet
        # (e.g. the DB is still downloading/being built) — conservative ~100 GB default
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
    f"[MetaFlux v2.1.0] mode={MODE}, input={INPUT_MODE}, "
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
    """Default target list for mode=amplicon — the files `rule all` requests for an amplicon run."""
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
    """Default target list for mode=shotgun — the files `rule all` requests for a shotgun run."""
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
