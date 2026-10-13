# Shared utilities: mode dispatch, path resolution, sample discovery,
# per-rule resource lookups, and mode-specific globals.
#
# This module is included first in the Snakefile (before any amplicon or shotgun
# rule file), and it runs entirely at PARSE time — before any actual analysis rule
# executes — because everything in it either validates the config outright (calling
# sys.exit with a specific message the moment something is wrong) or builds a Python
# value that later rule files read as a plain variable. A rule elsewhere in the
# workflow that writes `params: threads_for("some_rule")` or reads `AMPLICON_TYPE` is
# relying on this file having already run and set those names at the module level.
#
# Concretely, this file is where these things get decided, once, for the whole run:
#   - which mode is active (MODE) and the absolute paths every rule reads/writes
#     (FASTQ_DIR, OUT, LOGS, BENCH, ...)
#   - per-rule CPU/RAM lookups (threads_for, mem_mb_for) that every `resources:` and
#     `threads:` block in every rule file calls instead of hard-coding a number
#   - for amplicon mode: which marker is running (AMPLICON_TYPE), that marker's full
#     profile of facts (PROFILE, loaded from workflow/markers/<type>.yaml — see the
#     "Marker registry" section below), and every resolved reference-database path
#   - for shotgun mode: which profiler turns the cleaned reads into taxon tables
#     (shotgun.profiler: kraken2, sylph or both), the Kraken2/Bracken settings, the
#     sylph databases and lineage tables, and the PhiX/host decontamination toggles
#   - for sylph: sample names it cannot use, and out_dirs that already hold another
#     profiler's results (both checked right after sample discovery)
#   - which rule files the Snakefile leaves out of this run (SKIP_RULE_FILES)
#   - which samples exist and how their FASTQ files are named (SAMPLES, EXTN,
#     FASTQ_NAMING), discovered once here rather than re-globbed by every rule
#   - the final target file list `rule all` builds (all_targets(), at the bottom)
#
# If you are reading a rule file and wondering where some ALL_CAPS name or helper
# function it uses comes from, it is almost certainly defined somewhere in this file.

import csv
import hashlib
import json
import os
import re
import sys
import urllib.parse
import urllib.request
import yaml
from pathlib import Path

from snakemake.io import glob_wildcards


# ────────────────────────── Mode dispatch ──────────────────────
# The single config key that decides which half of the workflow exists at all. Every
# `if MODE == "amplicon":` / `else:` block further down in this file — and the choice
# of which rules/<mode>/*.smk directory the Snakefile includes — reads this one
# variable. Checked and failed here, at the very top, so a typo in the config
# (`mode: amplikon`) stops the run immediately with a clear message instead of
# surfacing later as a confusing missing-variable error deep in some rule file.
MODE = config["mode"].lower()
if MODE not in ("amplicon", "shotgun"):
    sys.exit(f"config.mode must be 'amplicon' or 'shotgun' (got: {config['mode']!r})")

# Rule files the Snakefile leaves out of this run (see the include loop there). It
# stays empty in amplicon mode. In shotgun mode the shotgun block below fills it from
# shotgun.profiler, so the Kraken2 rule files are not even read when only sylph runs,
# and the sylph rule file is not read when only Kraken2 runs.
SKIP_RULE_FILES = []


# ────────────────────────── Path resolution ────────────────────
# Every directory a rule reads from or writes to, computed once here as an absolute
# Path so that no rule file has to know or care what directory `snakemake` happened to
# be launched from. There is deliberately no `workdir:` directive anywhere in this
# workflow — shell commands run from wherever you invoked snakemake, and it is
# `.resolve()` here, not a workdir change, that turns the config's (possibly relative)
# `fastq_dir` / `out_dir` into paths that work regardless.
FASTQ_DIR = Path(config["input"]["fastq_dir"]).resolve()   # where this run's raw reads live
OUT       = Path(config["output"]["out_dir"]).resolve()    # where every output of this run goes
LOGS      = OUT / "logs"
BENCH     = OUT / "benchmarks"

# PhiX is shared (both modes); amplicon also uses the index basename.
REFDB_PHIX_FASTA  = Path(config["references"]["phix"]["fasta"]).resolve()
REFDB_PHIX_PREFIX = REFDB_PHIX_FASTA.with_suffix("")   # strip .fna → bowtie2 index basename

# A cache directory INSIDE the repository (not under OUT), so its contents survive
# across different runs and different output directories. Currently holds only the
# amplicon-probe JSON/FASTA cache (see PROBE_JSON below): the in-silico PCR of a given
# marker's primers against its reference is the same regardless of which samples you
# are running, so it is computed once here and reused by every subsequent run that
# shares the same marker + primer pair, rather than being tied to one out_dir.
PROJECT_DIR     = Path(workflow.basedir).parent.resolve()
PROBE_CACHE_DIR = PROJECT_DIR / "refdb" / "cache"

# Created immediately, before any rule needs it, because BioProject mode (shotgun's
# SRA-fetch path) writes downloaded FASTQs into this directory — it has to exist
# before the first download rule can run, and a plain "does this directory exist" read
# earlier in this file would otherwise fail on a first-ever run.
FASTQ_DIR.mkdir(parents=True, exist_ok=True)


# ────────────────────── Resource lookups ───────────────────────
# Every rule in every rule file asks for its CPU and memory allotment through these
# two functions rather than writing a number directly in its own `threads:` /
# `resources:` block. That means a single YAML block — config["resources"] — is the
# only place per-rule compute limits are ever set; changing how many threads
# `dada_seqtab` gets means editing one line in the config, not hunting through rule
# files for a hard-coded `threads: 8`. A rule not named in the config falls back to
# `threads_default` / `mem_mb_default`, so adding a new rule never requires touching
# the config first.
def threads_for(rule_name: str) -> int:
    cfg = config.get("resources", {})
    return int(cfg.get("threads", {}).get(rule_name, cfg.get("threads_default", 1)))


# Same lookup, but with a fallback chosen by the rule instead of threads_default.
# For almost every rule threads_default is the right thing to land on; the tree
# search is the exception, because there the thread count changes the RESULT and not
# just the runtime, so its safe value is 1 rather than whatever the machine's default
# happens to be. A run config that does name phylo_tree still wins.
def threads_for_or(rule_name: str, fallback: int) -> int:
    cfg = config.get("resources", {})
    return int(cfg.get("threads", {}).get(rule_name, fallback))


def mem_mb_for(rule_name: str) -> int:
    cfg = config.get("resources", {})
    return int(cfg.get("mem_mb", {}).get(rule_name, cfg.get("mem_mb_default", 1000)))


# ────────────────────── Marker registry ────────────────────────
# "A marker" here means everything that differs between 16S, ITS, 18S, gyrB and rpoB
# but is a FACT about that marker rather than a choice about your particular run: which
# reference database to measure amplicon length against, which database to classify
# against, whether a region-extraction tool exists for it, and how many taxonomic
# ranks its reference carries. One profile (a "pack") holds all of that per marker.
# Adding a marker means adding an entry here (plus its reference-fetch rule in
# amplicon/10_refdb.smk) — every choice that used to be a scattered
# `if AMPLICON_TYPE == "16S" else ...` now lives in one place, so a new marker can no
# longer silently inherit the ITS branch by falling through an incomplete if-chain.
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
#                       (assignTaxonomy columns / SINTAX ranks): the ranks of the
#                       marker's own reference files (16S on SILVA 144: Domain,
#                       Kingdom, Phylum ... Genus, plus the Species slot that
#                       addSpecies fills). A pack that ships several reference
#                       releases (16S) declares these three fields per release
#                       instead (see _select_pack_releases below).
#   rank_letters        the SINTAX one-letter rank codes paired to tax_levels.
#   rank_prefixes       the output prefixes (d__/k__/p__/…) paired to tax_levels.
#   sintax_tax_levels / sintax_rank_letters / sintax_rank_prefixes
#                       OPTIONAL per-method overrides. The rdp and sintax
#                       references usually share a rank model (16S, ITS), but not
#                       always: PR2 ships 9 ranks in its DADA2 file and only 8 in
#                       its UTAX file, where Division+Subdivision are merged into
#                       one 'p:' field. When absent these fall back to the fields
#                       above. Prefixes are keyed to rank MEANING, so a keep/
#                       discard token works the same under either method at every
#                       rank the two models share (18S: d__, sg__ and c__ ... s__;
#                       see 18S.yaml for dv__ and sbd__).
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


# ── SILVA '--other' handling (16S, SILVA 144 onwards) ──
# SILVA 144 marks some clades with the suffix "--other": a clade that does NOT
# contain the nomenclatural type of a taxon, but whose sequences SILVA still
# considers part of that taxon (e.g. "Bacillaceae--other" next to "Bacillaceae").
# amplicon.taxonomy.silva_other chooses what MetaFlux does with it:
#   verbatim (default) — keep the names exactly as SILVA writes them. The keep/
#                        discard filter and the addSpecies genus check then treat
#                        "X--other" as "X" (see 80a/80c), so a token such as
#                        g__Bradyrhizobium still removes g__Bradyrhizobium--other.
#   strip              — remove the suffix from every name before the classifier is
#                        trained, so "X" and "X--other" become one taxon.
# The kingdom is the one exception: it is always written without the suffix,
# whichever value is set (rule convert_silva_train, 10_refdb.smk), because at
# kingdom level the suffix only splits one kingdom into two labels.
# The value only matters for a release that ships the suffix (144); 138.2 has none.
# It is read here, before the marker packs are finalised, because it is part of the
# derived reference FILE NAMES below (the two settings give different files).
_silva_other_cfg = (
    ((config.get("amplicon") or {}).get("taxonomy") or {}).get("silva_other", "verbatim")
)
SILVA_OTHER = str(_silva_other_cfg).strip().lower()
if SILVA_OTHER not in ("verbatim", "strip"):
    sys.exit(
        f"[MetaFlux] amplicon.taxonomy.silva_other must be 'verbatim' or 'strip' "
        f"(got: {_silva_other_cfg!r})"
    )
# The tag that replaces "{other}" in a pack's reference file names.
_SILVA_OTHER_TAG = {"verbatim": "other-verbatim", "strip": "other-stripped"}[SILVA_OTHER]


def _select_pack_releases(packs):
    """Resolve marker packs that ship more than one reference release.

    A pack may carry a `releases:` block (today only 16S, for SILVA 144 and 138.2).
    The run config picks one with references.<release_selector>.release (for 16S:
    references.silva.release); without it the pack's default_release is used. The
    chosen block's fields are copied over the pack's top-level fields — for SILVA:
    probe_ref_tag, the rank model (tax_levels, rank_letters, rank_prefixes, which
    differ between 144 and 138.2) and `references` — so every later line of this
    file sees one ordinary pack and does not need to know releases exist.

    "{other}" in a reference file name is filled with the silva_other tag above.

    A pack without `releases:` is returned unchanged.
    """
    for name, pack in packs.items():
        releases = pack.get("releases")
        if not releases:
            continue
        selector = pack.get("release_selector")
        default  = str(pack.get("default_release"))
        # YAML reads an unquoted 138.2 as a number; str() makes "138.2" and 138.2 the same.
        ref_cfg   = (config.get("references") or {}).get(selector) or {}
        requested = ref_cfg.get("release") if isinstance(ref_cfg, dict) else None
        chosen    = default if requested is None else str(requested).strip()
        available = {str(k): v for k, v in releases.items()}
        if chosen not in available:
            sys.exit(
                f"[MetaFlux] references.{selector}.release = {requested!r} is not a release "
                f"the {name} marker pack ships. Available: {sorted(available)} "
                f"(default {default!r})."
            )
        block = available[chosen] or {}
        # Fill "{other}" in every reference file name of the chosen release.
        refs = {}
        for sym, spec in (block.get("references") or {}).items():
            spec = dict(spec)
            if "file" in spec:
                spec["file"] = spec["file"].format(other=_SILVA_OTHER_TAG)
            refs[sym] = spec
        for key, value in block.items():
            if key != "references":
                pack[key] = value
        # Release entries are added to (and win over) any references the pack
        # declares at top level, so a pack can mix release-specific and shared files.
        merged_refs = dict(pack.get("references") or {})
        merged_refs.update(refs)
        pack["references"] = merged_refs
        # Kept on the pack so logs and rules can say which release is running.
        pack["release"] = chosen
        del pack["releases"]
    return packs


MARKERS = _select_pack_releases(MARKERS)

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
# Runs only when mode: amplicon. Everything from here down to the matching `else:`
# (shotgun globals) is normally never even evaluated on a shotgun run — this whole
# block exists to turn a config plus a marker pack into the roughly twenty ALL_CAPS
# variables the amplicon rule files (rules/amplicon/*.smk) use directly: which marker
# is active, every resolved reference-database path, the taxonomy rank model, the
# region-extractor choice, and so on. The block is ordered roughly the way a run
# actually needs these facts: references first, then which marker was requested and
# its profile, then the probe (amplicon-length measurement) setup, then extraction and
# taxonomy. Several `sys.exit(...)` calls are mixed in along the way — each one is a
# config or marker-pack problem that would otherwise only surface as a confusing
# crash much later, once DADA2 or a taxonomy rule is already running.
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
    PR2_SINTAX           = _ref_path("pr2_sintax",        "pr2",   "sintax")  # built, not downloaded
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
        "pr2_sintax":        PR2_SINTAX,
        "silva_euk":         SILVA_EUK,
        "gyrb_dada2":        GYRB_DADA2,
        "rpob_frogs":        RPOB_FROGS,
    }

    # ── SILVA training set (only for a release whose download needs preparing) ──
    # SILVA 144 lineages have seven ranks, Domain;Kingdom;Phylum;...;Genus, where
    # SILVA 138.2 had six (no kingdom). MetaFlux classifies on the ranks of the
    # chosen release, but SILVA 144's file is not used exactly as published: rule
    # convert_silva_train (10_refdb.smk) first writes every kingdom without SILVA's
    # "--other" suffix (and, with silva_other: strip, removes it from the other
    # ranks too). The pack signals this by declaring the raw download
    # (silva_train_raw); a release without it (138.2) is used as downloaded, has no
    # conversion step, and SILVA_TRAIN_RAW stays None.
    SILVA_CONVERTS_TRAIN = "silva_train_raw" in _PACK_REF_FILES
    SILVA_TRAIN_RAW      = None
    if SILVA_CONVERTS_TRAIN:
        SILVA_TRAIN_RAW = _ref_path("silva_train_raw", "silva", "train_raw")
        REF_PATHS["silva_train_raw"] = SILVA_TRAIN_RAW
    SILVA_RELEASE = MARKERS.get("16S", {}).get("release")

    # SINTAX letters for the ranks of the SILVA training set, top-down, read by rule
    # convert_silva_sintax (10_refdb.smk) to write the 16S SINTAX file. They are the
    # 16S pack's letters for the chosen release without the species letter "s": the
    # SILVA toGenus files stop at genus (the Species slot is filled by addSpecies on
    # the rdp path and stays empty under sintax).
    #   SILVA 144   -> d, k, p, c, o, f, g
    #   SILVA 138.2 -> d, p, c, o, f, g
    # The 16S pack is read here, not PROFILE (the active marker, set further down),
    # because the rule exists in every amplicon run; only a 16S sintax run uses it.
    _silva_pack = MARKERS.get("16S", {})
    if "sintax_rank_letters" in _silva_pack:
        _silva_letters = _silva_pack["sintax_rank_letters"]
    else:
        _silva_letters = _silva_pack.get("rank_letters") or []
    SILVA_SINTAX_LETTERS = []
    for _letter in _silva_letters:
        if _letter != "s":
            SILVA_SINTAX_LETTERS.append(_letter)

    # references.silva.kingdom_map pointed at the phylum -> kingdom table MetaFlux
    # once wrote for SILVA 144. MetaFlux now classifies the kingdom like every other
    # rank and writes no such table, so the setting is no longer read. A warning, not
    # a stop: the run is unaffected, but the user should know the line does nothing.
    _silva_ref_cfg = config["references"].get("silva")
    if isinstance(_silva_ref_cfg, dict) and "kingdom_map" in _silva_ref_cfg:
        sys.stderr.write(
            "[MetaFlux] warning: references.silva.kingdom_map is no longer used and is "
            "ignored. MetaFlux classifies the SILVA 144 kingdom like every other rank "
            "and writes no phylum -> kingdom table; remove the line from your run "
            "config.\n"
        )

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
    # "md5" is the checksum the pack lists for the file (empty string when it lists
    # none); the generated fetch rule checks the download against it.
    PACK_FETCHABLE = {
        _sym: {"path": REF_PATHS[_sym], "url": _spec["url"], "md5": str(_spec.get("md5") or "")}
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
    # Reference lineages the probe leaves out of its length statistics (16S:
    # Chloroplast, Mitochondria, Eukaryota — see 40a_amplicon_probe.py). The list
    # changes the result, so it is written into the cache file name as well: a
    # probe computed with a different list (or none, as before v2.5.0) is never
    # picked up by mistake. Markers without the field keep their old cache names.
    PROBE_EXCLUDE_TAXA = list(PROFILE.get("probe_exclude_taxa") or [])
    if PROBE_EXCLUDE_TAXA:
        PROBE_REF_TAG += "_excl-" + "-".join(sorted(t.lower() for t in PROBE_EXCLUDE_TAXA))
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
    # replacing the previously hardcoded 7-rank tuples. The ranks are those of the
    # marker's reference: 16S has eight slots on SILVA 144 (Domain, Kingdom, Phylum
    # ... Species) and seven on SILVA 138.2 (no kingdom), declared per release in
    # 16S.yaml; ITS, gyrB and rpoB have seven; 18S/PR2 has nine under rdp and eight
    # under sintax.
    #
    # Guard first: a pack with releases (16S) gets its rank model only from the
    # chosen release block, merged in by _select_pack_releases above. A pack that
    # declares tax_levels, rank_letters or rank_prefixes nowhere, neither at its top
    # level nor in that block (an incomplete pack of one's own in config/markers/,
    # for example), would otherwise stop on the first line below with a bare Python
    # KeyError. Stop here with a sentence that names the pack and the missing fields.
    _missing_rank_fields = []
    for _field in ("tax_levels", "rank_letters", "rank_prefixes"):
        if not PROFILE.get(_field):
            _missing_rank_fields.append(_field)
    if _missing_rank_fields:
        if PROFILE.get("release"):
            _where = f"at its top level or in its release block {PROFILE['release']!r}"
        else:
            _where = "at its top level"
        sys.exit(
            f"[MetaFlux] marker pack {AMPLICON_TYPE!r} has no complete rank model: "
            f"{', '.join(_missing_rank_fields)} not declared {_where}. A pack in "
            f"config/markers/ replaces the shipped one completely; compare it with "
            f"workflow/markers/{AMPLICON_TYPE}.yaml."
        )

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
    # Lists of rank-prefixed tokens (e.g. d__Bacteria, o__Chloroplast) matched
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
            "    keep:    [d__Bacteria, d__Archaea]\n"
            "    discard: [o__Chloroplast, f__Mitochondria]   # 16S\n"
            "  or keep: [k__Fungi]   discard: []              # ITS"
        )
    FILTER_ENABLED   = bool(_filter_cfg.get("enabled", True))
    FILTER_KEEP      = _as_taxon_list(_filter_cfg.get("keep"))
    FILTER_DISCARD   = _as_taxon_list(_filter_cfg.get("discard"))

    # Best-hit check of the discard filter: an ASV is also removed when its
    # classifier best hit carries a discard token at the token's rank with at least
    # this much support (percent, 0-100; rdp bootstrap, or SINTAX confidence x 100),
    # even if that rank is below min_boot / sintax_cutoff and so missing from the
    # reported taxonomy. It exists for organelle 16S (see 80a/80c). null or false
    # switches it off, which gives the filter of MetaFlux v2.4.0 and earlier.
    #
    # YAML reads a bare off / no / false as the boolean False, so those switch it off
    # too. A bare on / yes / true arrives as True, and Python counts True as the
    # number 1 (float(True) is 1.0), so without the check below the run would quietly
    # use a 1 % threshold. The value has to be a number, so True is refused.
    _dms = _filter_cfg.get("discard_min_support", 50)
    if _dms is True:
        sys.exit(
            "[MetaFlux] amplicon.taxonomy.filter.discard_min_support must be a number "
            "from 0 to 100, or null (also a bare false / off / no) to switch the "
            "best-hit check off (got: True; YAML reads on / yes / true as True, so "
            "write the number, e.g. 50)"
        )
    if _dms is None or _dms is False:
        FILTER_DISCARD_MIN_SUPPORT = None
    else:
        try:
            FILTER_DISCARD_MIN_SUPPORT = float(_dms)
        except (TypeError, ValueError):
            FILTER_DISCARD_MIN_SUPPORT = -1.0
        if not 0 <= FILTER_DISCARD_MIN_SUPPORT <= 100:
            sys.exit(
                f"[MetaFlux] amplicon.taxonomy.filter.discard_min_support must be a number "
                f"from 0 to 100, or null (also a bare false / off / no) to switch the "
                f"best-hit check off (got: {_dms!r})"
            )

    # ── Organelle filter (16S: Metaxa2's own origin call) ──
    # Metaxa2 (rule target_extract, 70_extract.smk) sorts every ASV it extracts into
    # one origin: bacteria, archaea, eukaryota, chloroplast, mitochondria or uncertain.
    # With organelles: metaxa2 (the default) the taxonomy filter also removes every ASV
    # Metaxa2 called mitochondria or chloroplast, after keep/discard and the best-hit
    # check. It catches plant mitochondria that no token can reach, because the
    # classifier's best hit for them is a bacterium: on the 16S test set, 8 ASVs whose
    # nearest SILVA reference is a mitochondrion (97-99 % identity) are classified as
    # bacteria, mostly Enterobacterales (best hit Pantoea / Erwiniaceae), with SINTAX on
    # SILVA 138.2 and with rdp on either release, and Metaxa2 calls all 8 mitochondrial.
    # Only the organelle call is used. Metaxa2's other calls (eukaryota, uncertain, the
    # domains) remove nothing; the keep tokens still decide domains.
    #
    # The filter can only act when Metaxa2 actually ran: the marker's extractor is
    # metaxa2 (only 16S) and extraction is enabled. config/config.yaml sets the key for
    # every run, because it is the base layer under every --configfile, so for the other
    # markers, or with extraction off, the setting is simply inactive. That gets one
    # info line below, never an error.
    #
    # YAML reads a bare `off` (also no / false) as the boolean False, so False means
    # off here, and so does null, as for discard_min_support. A missing key means the
    # default, metaxa2. `on` / `yes` / `true` arrive as True and are refused: the value
    # names where the organelle call comes from, and today that can only be metaxa2.
    _organelles_cfg = _filter_cfg.get("organelles", "metaxa2")
    if _organelles_cfg is None or _organelles_cfg is False:
        ORGANELLE_SOURCE = "off"
    else:
        ORGANELLE_SOURCE = str(_organelles_cfg).strip().lower()
    if ORGANELLE_SOURCE not in ("metaxa2", "off"):
        # A bare on / yes / true reaches us as True; say so, or "(got: True)" puzzles.
        _yaml_hint = ""
        if _organelles_cfg is True:
            _yaml_hint = " — YAML reads on / yes / true as True; write metaxa2"
        sys.exit(
            f"[MetaFlux] amplicon.taxonomy.filter.organelles must be 'metaxa2' (remove the "
            f"ASVs Metaxa2 calls mitochondria or chloroplast) or 'off' "
            f"(got: {_organelles_cfg!r}{_yaml_hint})"
        )

    # Metaxa2 runs in this run, and so writes its per-ASV origin table
    # (5.dada2/metaxa2_origin.tsv), only for a metaxa2 marker with extraction on.
    # Whenever it does, rule assign_taxonomy reads the table, also with organelles: off,
    # so that the metaxa2_origin column of 6.taxonomy/asv_table_unfiltered.txt is filled
    # in either case. For every other marker the rule's inputs stay as they were.
    METAXA2_RUNS = MARKER_EXTRACTOR == "metaxa2" and EXTRACTION_ENABLED
    METAXA2_ORIGIN_TSV = OUT / "5.dada2" / "metaxa2_origin.tsv"

    # Whether ASVs are actually removed on Metaxa2's organelle call. filter.enabled is
    # the master switch of the whole taxonomy filter, so false switches this off too.
    ORGANELLE_FILTER = ORGANELLE_SOURCE == "metaxa2" and METAXA2_RUNS and FILTER_ENABLED

    # One line on every amplicon run, saying what the organelle filter does and why.
    if ORGANELLE_FILTER:
        _organelle_note = ("on — ASVs that Metaxa2 calls mitochondria or chloroplast "
                           "are removed")
    elif ORGANELLE_SOURCE == "off":
        _organelle_note = "off (amplicon.taxonomy.filter.organelles: off)"
    elif MARKER_EXTRACTOR != "metaxa2":
        _organelle_note = (f"inactive — marker {AMPLICON_TYPE} does not run Metaxa2, "
                           "so there is no organelle call to use")
    elif not EXTRACTION_ENABLED:
        _organelle_note = ("inactive — amplicon.extraction.enabled is false, so Metaxa2 "
                           "does not run")
    else:
        _organelle_note = "inactive — amplicon.taxonomy.filter.enabled is false"
    sys.stderr.write(f"[MetaFlux] organelle filter: {_organelle_note}\n")

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

    # One line on every 16S run naming the SILVA release, the silva_other setting
    # and the training-set file, so the log shows which reference classified the run.
    if AMPLICON_TYPE == "16S":
        sys.stderr.write(
            f"[MetaFlux] 16S reference: SILVA {SILVA_RELEASE} "
            f"(silva_other: {SILVA_OTHER}; trainset {TAXONOMY_REFDB.name})\n"
        )

    # ── Optional 16S phylogeny (stage 85) ──────────────────────────────────
    # Opt-in de novo tree for the run's ASVs, so users can compute phylogenetic
    # diversity (Faith's PD, UniFrac) in their own R session. MetaFlux's job ends
    # when the Newick is written — it computes no diversity statistics itself.
    #
    # Everything here is resolved at parse time so that rules/amplicon/85_phylogeny.smk
    # can be defined (or not defined at all) without re-reading the config. Every key
    # is read with an explicit default, following the shotgun kmer_evidence block
    # above: a config written before this feature existed — including one passed with
    # --configfile from outside the repo, where config/config.yaml is NOT loaded as a
    # base layer (see the Snakefile) — must still parse untouched. The guards below
    # can therefore only ever fire on a value the user actually wrote, never on a
    # missing key.
    #
    # WHY 16S ONLY. A de novo tree is meaningful only for a marker that can be
    # globally aligned across divergent taxa. ITS is not reliably alignable above
    # about genus level; gyrB co-amplifies its paralogue parE, which enters a tree as
    # a long out-group branch and inflates Faith's PD; rpoB is mutationally saturated
    # at all three codon positions at phylum/domain scope. 18S is the one marker worth
    # revisiting, but branch lengths from short 18S amplicons are documented to be
    # underestimated and nobody has validated UniFrac/PD from them. See
    # docs/amplicon/phylogeny.md for the citations.
    _phylo_cfg = amp_cfg.get("phylogeny", {}) or {}

    PHYLO_ENABLED          = bool(_phylo_cfg.get("enabled", False))
    PHYLO_BACKEND          = str(_phylo_cfg.get("backend", "iqtree")).strip().lower()
    PHYLO_ALIGNER_STRATEGY = str(_phylo_cfg.get("aligner_strategy", "auto")).strip().lower()
    PHYLO_MODEL            = str(_phylo_cfg.get("model", "auto")).strip()
    PHYLO_SUPPORT          = bool(_phylo_cfg.get("support", False))
    # extra_args must be a mapping (tool -> flag string). A bare string here would
    # otherwise die inside dict() with a Python traceback instead of a sentence.
    _phylo_extra_raw = _phylo_cfg.get("extra_args", {}) or {}
    if not isinstance(_phylo_extra_raw, dict):
        sys.exit(
            "[MetaFlux] amplicon.phylogeny.extra_args must be a mapping of tool -> flag "
            "string, with keys mafft, iqtree, fasttree and/or raxml_ng, e.g.\n"
            "  extra_args:\n    iqtree: \"--fast\"\n"
            f"(got: {_phylo_extra_raw!r})"
        )
    PHYLO_EXTRA_ARGS       = dict(_phylo_extra_raw)

    # There is deliberately NO rooting option. The tree is exported unrooted, full stop:
    # users prune it in R first (decontam against negative controls, abundance
    # filtering), and any root placed here would stop being valid the moment a tip is
    # removed. A `root:` key in the config is from a pre-release draft of this module —
    # refuse it with a pointer rather than silently ignoring a setting the user
    # believes is doing something.
    if "root" in _phylo_cfg:
        sys.exit(
            "[MetaFlux] amplicon.phylogeny.root is not a setting. MetaFlux exports the "
            "tree UNROOTED only — root it in R after pruning to your final ASV set "
            "(see docs/amplicon/phylogeny.md). Remove the 'root' key."
        )

    # Marker guard, in the same "force off with a warning" style as the
    # EXTRACTION_ENABLED guard above: an unsupported marker is a recoverable config
    # combination, not a typo, so coerce it to the safe state and say so rather than
    # aborting a run whose taxonomy half is perfectly valid.
    if PHYLO_ENABLED and AMPLICON_TYPE != "16S":
        sys.stderr.write(
            f"[MetaFlux] warning: amplicon.phylogeny is 16S-only; marker {AMPLICON_TYPE} "
            "is not supported (see docs/amplicon/phylogeny.md for the per-marker "
            "reasons). Ignoring amplicon.phylogeny.enabled: true — no tree will be built.\n"
        )
        PHYLO_ENABLED = False

    # The remaining guards validate values the user typed. They run even when the
    # feature is disabled, so a typo is caught the first time it is written rather
    # than only on the run where it finally matters.
    if PHYLO_BACKEND not in ("iqtree", "fasttree", "raxml-ng"):
        sys.exit(
            f"[MetaFlux] amplicon.phylogeny.backend must be 'iqtree', 'fasttree' or 'raxml-ng' "
            f"(got: {PHYLO_BACKEND!r})"
        )
    if PHYLO_ALIGNER_STRATEGY not in ("auto", "linsi", "fftns2"):
        sys.exit(
            f"[MetaFlux] amplicon.phylogeny.aligner_strategy must be 'auto', 'linsi' or 'fftns2' "
            f"(got: {PHYLO_ALIGNER_STRATEGY!r})"
        )
    # Unknown extra_args keys are a silent no-op otherwise: writing `iqtre: "-x"`
    # would simply never reach any tool, and the user would spend an afternoon
    # wondering why their flag did nothing.
    _PHYLO_EXTRA_KEYS = ("mafft", "iqtree", "fasttree", "raxml_ng")
    _bad_extra = sorted(set(PHYLO_EXTRA_ARGS) - set(_PHYLO_EXTRA_KEYS))
    if _bad_extra:
        sys.exit(
            f"[MetaFlux] amplicon.phylogeny.extra_args has unknown key(s) {_bad_extra}. "
            f"Valid keys are: {', '.join(_PHYLO_EXTRA_KEYS)}"
        )

    # RAxML-NG support values: bootstrapping dominates this backend's runtime, and
    # support values are not inputs to UniFrac or Faith's PD, so the setting buys
    # nothing here. Force it off with a warning rather than silently ignoring it.
    if PHYLO_ENABLED and PHYLO_SUPPORT and PHYLO_BACKEND == "raxml-ng":
        sys.stderr.write(
            "[MetaFlux] warning: amplicon.phylogeny.support is not available for the "
            "raxml-ng backend (bootstrapping dominates its runtime and support values "
            "are not used by any downstream diversity metric). Proceeding without support.\n"
        )
        PHYLO_SUPPORT = False

    # Model resolution. `auto` is MetaFlux's own choice per backend, pinned rather
    # than delegated to the tool's search: IQ-TREE's own default (-m MFP) re-searches
    # 22 base models on every run, which makes runtime data-dependent. On the 16S test
    # data the selection picked TPM3u+R5 (real ASVs) and GTR+F+I+R5 (reference
    # fragments) — GTR-family matrices whose choice barely moves the distances the
    # diversity metrics use (see docs/amplicon/phylogeny.md). IQ-TREE itself documents
    # this select-once-then-pin pattern.
    #
    # A tool's own model search is still reachable — it is the class of automatic
    # behaviour that WRITES DOWN what it picked (IQ-TREE's choice lands in .iqtree,
    # RAxML-NG's in .raxml.bestModel), so provenance survives. Set model: MFP or
    # model: DNA to use it. What is never reachable is automatic behaviour that
    # depends on the machine rather than the data (-T AUTO, --threads auto), because
    # the same input must not behave differently on different hardware.
    _PHYLO_MODEL_AUTO = {
        "iqtree":   "GTR+F+G4",
        "fasttree": "gtr",       # FastTree has only JC and GTR; see the guard below
        "raxml-ng": "GTR+G",
    }
    if PHYLO_MODEL.lower() == "auto":
        PHYLO_MODEL_RESOLVED = _PHYLO_MODEL_AUTO[PHYLO_BACKEND]
    else:
        PHYLO_MODEL_RESOLVED = PHYLO_MODEL

    # FastTree is the one backend with no model selection at all: its nucleotide
    # options are Jukes-Cantor (the default, which is why -gtr is passed explicitly)
    # and GTR, and nothing else. An unrecognised model string would otherwise be
    # dropped on the floor and the user would get a GTR tree while believing they had
    # asked for something else — quietly wrong output, which this workflow treats as
    # worse than a loud failure.
    if PHYLO_BACKEND == "fasttree" and PHYLO_MODEL_RESOLVED.lower() not in ("jc", "gtr"):
        sys.exit(
            f"[MetaFlux] amplicon.phylogeny.model {PHYLO_MODEL!r} is not available for the "
            "fasttree backend, which implements only Jukes-Cantor and GTR (no model "
            "selection). Use 'auto' (= gtr), 'gtr', or 'jc' — or switch to "
            "backend: iqtree / raxml-ng, which do support model selection."
        )

    # Advisory only — never forces anything. The whole point of feeding the tree the
    # POST-taxonomy ASV set is that MetaFlux's contaminant filter has already removed
    # chloroplast, mitochondrial and wrong-domain sequences, which are exactly the
    # off-target ASVs that show up in a de novo tree as long branches and distort
    # UniFrac and Faith's PD. When the filter removes little or nothing, the tree
    # still gets built, but from a poorly filtered set — so say so before the run
    # rather than leaving the user to discover it in the long-branch QC report.
    #
    # Two situations get different messages, because they protect the tree to a
    # different degree:
    #   - nothing filters at all: filter.enabled is false, or the keep/discard lists
    #     are empty and the organelle filter is not acting either (organelles: off,
    #     or extraction off so Metaxa2 never ran);
    #   - the keep/discard lists are empty but the organelle filter still acts
    #     (the shipped 16S default): Metaxa2's mitochondria and chloroplast calls
    #     are removed, but eukaryotic, unclassified and other off-target ASVs are not.
    _filter_lists_empty = not (FILTER_KEEP or FILTER_DISCARD)
    if PHYLO_ENABLED and not FILTER_ENABLED:
        sys.stderr.write(
            "[MetaFlux] warning: phylogeny is enabled but the contaminant filter "
            "(amplicon.taxonomy.filter.enabled) is off. The tree will include any "
            "chloroplast, mitochondrial or off-target ASVs, which typically appear as "
            "long branches and inflate Faith's PD. Check "
            "stats/phylogeny/phylogeny_qc.json's long-branch report.\n"
        )
    elif PHYLO_ENABLED and _filter_lists_empty and not ORGANELLE_FILTER:
        sys.stderr.write(
            "[MetaFlux] warning: phylogeny is enabled but the contaminant filter "
            "(amplicon.taxonomy.filter) has empty keep/discard lists and the organelle "
            "filter is not acting, so no ASV is removed. The tree will include any "
            "chloroplast, mitochondrial or off-target ASVs, which typically appear as "
            "long branches and inflate Faith's PD. Check "
            "stats/phylogeny/phylogeny_qc.json's long-branch report.\n"
        )
    elif PHYLO_ENABLED and _filter_lists_empty:
        sys.stderr.write(
            "[MetaFlux] warning: phylogeny is enabled but amplicon.taxonomy.filter has "
            "empty keep/discard lists. The organelle filter still removes the ASVs "
            "Metaxa2 calls mitochondria or chloroplast, but eukaryotic, unclassified "
            "and other off-target ASVs will reach the tree, where they typically appear "
            "as long branches and inflate Faith's PD. Consider "
            "keep: [d__Bacteria, d__Archaea] and discard: [o__Chloroplast, "
            "f__Mitochondria]. Check stats/phylogeny/phylogeny_qc.json's long-branch "
            "report.\n"
        )

    # Taxonomy always reads the length-filtered outputs (dada_length_filter
    # runs after target_extract — see amplicon/60_dada2.smk).
    def _seqs_for_taxonomy(wildcards=None):
        return OUT / "5.dada2" / "seqs_lenfilt.fasta"

    def _seqtab_for_taxonomy(wildcards=None):
        return OUT / "5.dada2" / "seqtab_lenfilt_head_names.txt"


# ────────────────────── Shotgun-only globals ───────────────────
# The shotgun counterpart of the amplicon block above — runs only when mode: shotgun,
# and only one of the two blocks ever executes for a given run. Shotgun mode has no
# marker packs and no probe/truncation machinery to resolve, so this block is
# considerably shorter. It first reads which profiler turns the cleaned reads into
# taxon tables (shotgun.profiler) and, from that, which rule files the Snakefile reads.
# It then validates the user-supplied Kraken2 database (only when Kraken2 runs),
# resolves the decontamination and classification settings rules/shotgun/*.smk read
# directly, and works out how much memory to give the Kraken2 rule from the database's
# own size on disk (see _kraken_hash_size_mb below) rather than a fixed guess. Last, it
# checks the sylph databases and lineage tables (only when sylph runs) and works out
# sylph_profile's memory the same way, from the database files' size.
else:  # MODE == "shotgun"
    sht_cfg = config["shotgun"]

    # ── Which profiler turns the cleaned reads into taxon tables (shotgun.profiler) ──
    # The read cleaning (PhiX, host, fastp) is the same whichever profiler is chosen:
    #   kraken2  Kraken2 + Bracken, exactly as in earlier versions (the default)
    #   sylph    sylph + sylph-tax instead: genomes detected from k-mer sketches of the
    #            reads, abundance estimated from genome coverage
    #   both     both chains on the same trimmed reads, side by side; never merged
    # Read with a default, so a config written before this key existed still runs as
    # kraken2, with nothing changed in a run that has already finished.
    _profiler_raw = sht_cfg.get("profiler", "kraken2")
    PROFILER = str(_profiler_raw).strip().lower()

    # Stop on any other value here, before a single rule file is read, so that a typo
    # such as "sylp" is never quietly treated as one of the three choices.
    if PROFILER not in ("kraken2", "sylph", "both"):
        sys.exit(
            f"[MetaFlux] shotgun.profiler must be kraken2, sylph or both (got: {_profiler_raw!r}).\n"
            "  kraken2  Kraken2 + Bracken, as in earlier MetaFlux versions (the default)\n"
            "  sylph    sylph genome profiling against references.sylph_databases\n"
            "  both     both on the same reads, side by side; their tables are never merged"
        )
    RUN_KRAKEN2 = PROFILER in ("kraken2", "both")
    RUN_SYLPH   = PROFILER in ("sylph", "both")

    # Whole rule files belong to one profiler, and the Snakefile skips the files of a
    # profiler that does not run (see SKIP_RULE_FILES at the top of this file).
    # Every file in the Kraken2 list reads Kraken2 output; 80_stats.smk is there because
    # its read-tracking script reads the Kraken2 reports. A new rule file that needs
    # Kraken2 output belongs in this list too. 30_preprocess.smk (read cleaning) is in
    # neither list, because every profiler starts from the reads it cleans.
    KRAKEN2_RULE_FILES = [
        "35_kraken2.smk",
        "40_kmer_evidence.smk",
        "45_abundance.smk",
        "55_extract_taxa.smk",
        "80_stats.smk",
    ]
    SYLPH_RULE_FILES = [
        "38_sylph.smk",
    ]
    if not RUN_KRAKEN2:
        SKIP_RULE_FILES += KRAKEN2_RULE_FILES
    if not RUN_SYLPH:
        SKIP_RULE_FILES += SYLPH_RULE_FILES

    # ── Kraken2 database (user-supplied, never downloaded; checked only when Kraken2 runs) ──
    # A sylph-only run may keep the template's /path/to/kraken2_db placeholder, and that
    # must not stop it. When Kraken2 runs, the value and the message are exactly as
    # before: rules kraken2 and bracken pass str(KRAKEN_DB) as a parameter (a changed
    # value would re-run finished runs), and three docs pages quote the message
    # (quick-start.md, installation.md, shotgun/databases.md).
    _kraken_db_cfg = (config.get("references") or {}).get("kraken_db")
    if _kraken_db_cfg:
        KRAKEN_DB = Path(_kraken_db_cfg).resolve()
    else:
        KRAKEN_DB = None
    if RUN_KRAKEN2 and (KRAKEN_DB is None or not KRAKEN_DB.exists()):
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

    # The rank Bracken re-estimates at (rule bracken), and the rank rule kraken_biom
    # builds the OTU table at (both in 45_abundance.smk). Only these seven letters
    # work end to end: Bracken matches the value against the rank codes in the
    # Kraken2 report, and kraken-biom knows no others. Without this check a wrong
    # value failed late or not at all: a lowercase letter crashed Bracken, but only
    # after Kraken2 had finished, and "K" went through Bracken and then gave an
    # empty table. Stop here instead, before anything runs.
    if TAX_LEV not in ("D", "P", "C", "O", "F", "G", "S"):
        sys.exit(f"[MetaFlux] shotgun.bracken.tax_lev must be one of D, P, C, O, F, G, S "
                 f"(a single capital letter), got {TAX_LEV!r}")

    # Bracken's minimum-read threshold. Either a fixed integer, or the string "auto",
    # in which case the rule works it out per sample from that sample's own classified
    # read count (see the config comment, and the shell block in 45_abundance.smk).
    # It cannot be resolved here: the number of classified reads is only known after
    # Kraken2 has run, which is long after the workflow is parsed.
    _bracken_thresh = sht_cfg["bracken"]["threshold"]
    if isinstance(_bracken_thresh, str) and _bracken_thresh.strip().lower() == "auto":
        BRACKEN_THRESH = "auto"
    else:
        BRACKEN_THRESH = int(_bracken_thresh)
        if BRACKEN_THRESH < 0:
            sys.exit("[MetaFlux] shotgun.bracken.threshold must be >= 0, or the word 'auto'")
    BRACKEN_ALPHA  = float(sht_cfg["bracken"].get("threshold_alpha", 5.0e-05))
    BRACKEN_MIN_T  = int(sht_cfg["bracken"].get("threshold_min", 10))
    if not 0 < BRACKEN_ALPHA < 1:
        sys.exit(f"[MetaFlux] shotgun.bracken.threshold_alpha must be a small fraction "
                 f"(e.g. 5.0e-05), got {BRACKEN_ALPHA}")
    EXTRACT_TAXA     = list(sht_cfg.get("extract_taxa", []) or [])

    # Read extraction (55_extract_taxa.smk) picks reads by the taxon Kraken2 gave each
    # one, and sylph labels no single read. Under profiler: sylph there is nothing to
    # pick reads by, so warn and switch extraction off. Under profiler: both Kraken2
    # runs, so extract_taxa is kept exactly as configured.
    if PROFILER == "sylph" and EXTRACT_TAXA:
        sys.stderr.write(
            f"[MetaFlux] warning: shotgun.extract_taxa needs Kraken2's per-read output, which "
            f"profiler: sylph does not produce. Ignoring extract_taxa {EXTRACT_TAXA}: no reads "
            f"will be extracted (profiler: both keeps it).\n"
        )
        EXTRACT_TAXA = []

    # ── Species-level k-mer evidence (rule kmer_evidence, 40_kmer_evidence.smk) ──
    # The rule itself always runs and always writes its evidence table, plot and
    # gated report. KMER_GATE_ENABLED decides one thing only: whether rule bracken
    # takes the gated report or Kraken2's original (see 45_abundance.smk). Defaults
    # are duplicated here so a config written before this block existed still parses.
    _kmer_cfg          = sht_cfg.get("kmer_evidence", {}) or {}
    KMER_GATE_ENABLED  = bool(_kmer_cfg.get("enabled", False))
    KMER_MIN_DISTINCT  = int(_kmer_cfg.get("min_distinct_minimizers", 333))
    KMER_MIN_READS     = int(_kmer_cfg.get("min_reads", 0))
    KMER_PROTECT_FRAC  = float(_kmer_cfg.get("never_prune_above_fraction", 0.01))

    if KMER_MIN_DISTINCT < 0 or KMER_MIN_READS < 0:
        sys.exit("[MetaFlux] shotgun.kmer_evidence: min_distinct_minimizers and "
                 "min_reads must be >= 0")
    if not 0 < KMER_PROTECT_FRAC <= 1:
        sys.exit("[MetaFlux] shotgun.kmer_evidence.never_prune_above_fraction must be "
                 f"a fraction in (0, 1], got {KMER_PROTECT_FRAC}")

    # The k-mer evidence gate filters Kraken2's report before Bracken, and profiler:
    # sylph runs neither, so the setting would do nothing. Warn and switch it off. This
    # sits before the rank warning below, so that warning never fires for a gate that
    # does not run. Under profiler: both the gate works as configured.
    if PROFILER == "sylph" and KMER_GATE_ENABLED:
        sys.stderr.write(
            "[MetaFlux] warning: shotgun.kmer_evidence filters the Kraken2 report before "
            "Bracken, and profiler: sylph runs neither. Ignoring kmer_evidence.enabled: true.\n"
        )
        KMER_GATE_ENABLED = False

    # The default threshold is converted from species-level literature (see the
    # config comment). Pointed at any coarser rank it is far too permissive to mean
    # anything — a genus pools every species under it — so say so rather than let it
    # look like a considered choice.
    if KMER_GATE_ENABLED and TAX_LEV != "S":
        sys.stderr.write(
            f"[MetaFlux] warning: shotgun.kmer_evidence is enabled at rank {TAX_LEV}, but "
            f"min_distinct_minimizers={KMER_MIN_DISTINCT} was derived for SPECIES. Check the "
            f"evidence tables before trusting the pruned output at this rank.\n"
        )

    # ── OTU-table taxon filter (shotgun counterpart of amplicon keep/discard) ──
    # Same rank-aware token matching; applied to the final otu_table. Config-only,
    # no code default — default OFF (see config comment).
    _sht_filter        = sht_cfg.get("taxonomy_filter", {}) or {}
    OTU_FILTER_ENABLED = bool(_sht_filter.get("enabled", False))
    OTU_FILTER_KEEP    = _as_taxon_list(_sht_filter.get("keep"))
    OTU_FILTER_DISCARD = _as_taxon_list(_sht_filter.get("discard"))

    # The taxon filter acts on otu_table.tsv (Kraken2 + Bracken) only. The sylph tables
    # are written unfiltered, so that each column still adds up to the read pairs after
    # QC. Say so, so nobody assumes the filter also cleaned the sylph tables.
    if RUN_SYLPH and OTU_FILTER_ENABLED:
        sys.stderr.write(
            "[MetaFlux] warning: shotgun.taxonomy_filter applies to otu_table.tsv (Kraken2 + "
            "Bracken) only; the sylph tables in 03.abundance/ are written unfiltered, so each "
            "column still adds up to the read pairs after QC.\n"
        )

    # Kraken2 mem_mb is derived from the mmap mode + actual DB size at parse time:
    #   mmap on   → ~20 GB private workspace per process (DB is shared via the OS
    #               file cache, so it doesn't count against this rule's budget).
    #   mmap off  → hash.k2d size + 10% headroom for buffers / read I/O.
    #               Self-adjusts to the actual DB (16 GB pluspf_16gb, ~100 GB pluspf,
    #               200+ GB core_nt, etc).
    def _kraken_hash_size_mb(db_path: Path):
        h = db_path / "hash.k2d"
        return (h.stat().st_size // (1024 * 1024)) if h.exists() else None

    # No kraken_db at all is possible only when Kraken2 does not run (checked above);
    # then there is no database size to read, and the memory figure is never used.
    if KRAKEN_DB is not None:
        _hash_mb = _kraken_hash_size_mb(KRAKEN_DB)
    else:
        _hash_mb = None
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

    # ── sylph databases and lineage tables (profiler sylph | both; user-supplied) ──
    # Everything the sylph rules (38_sylph.smk) need is resolved and checked here,
    # before any job starts. That matters more for sylph than for Kraken2, because
    # sylph and sylph-tax accept most mistakes with exit code 0 and carry on:
    #   - a database file with an unexpected suffix is skipped with only a warning;
    #   - a database sketched more densely (c 100) than MetaFlux sketches the reads
    #     makes sylph skip every sample;
    #   - a lineage table from another release silently gives genomes wrong names.
    # What can only be seen at run time is checked in the rules themselves.
    # The values built here (left empty under profiler: kraken2):
    #   SYLPH_DB_PATHS        the database files, in config order (sylph_profile inputs)
    #   SYLPH_TAXONOMY_ARGS   what sylph-tax -t receives: built-in names or table paths,
    #                         each once (sylph_taxprof)
    #   SYLPH_TAXONOMY_FILES  the lineage tables behind them (sylph_taxprof inputs)
    #   SYLPH_DB_LABELS       'file (taxonomy)' per database: start-up line, table line 1
    #   SYLPH_TAXONOMY_DIR    references.sylph_taxonomy_dir as a full path, or None
    #   SYLPH_PROFILE_MEM_MB  memory for the one sylph_profile job

    # How densely sylph_sketch samples the reads: 1 k-mer in 200. Fixed by design. A
    # sketch can be compared with databases built at the same or a coarser c, so the
    # prebuilt c 200 and c 1000 databases both work, and c 200 finds low-coverage
    # genomes best.
    SYLPH_SKETCH_C = 200

    # sylph-tax 1.9.2's built-in taxonomy names, and the lineage table each one reads
    # from references.sylph_taxonomy_dir. Copied from sylph-tax 1.9.2
    # (sylph_tax/metadata_files.py); update together with the pin in
    # workflow/envs/sylph.yaml. Used here only to check, before anything runs, that
    # the table a run needs is really in the folder.
    SYLPH_TAX_FILES = {
        "FungiRefSeq-latest":     "fungi_refseq_2025-10-11_metadata.tsv.gz",
        "FungiRefSeq-2024-07-25": "fungi_refseq_2024-07-25_metadata.tsv.gz",
        "GTDB_r214":              "gtdb_r214_metadata.tsv.gz",
        "GTDB_r220":              "gtdb_r220_metadata.tsv.gz",
        "GTDB_r226":              "gtdb_r226_metadata.tsv.gz",
        "GTDB_r232":              "gtdb_r232_metadata.tsv.gz",
        "GlobDB_r226":            "globdb_r226_sylph_tax.tsv.gz",
        "GlobDB_r232":            "globdb_r232_taxonomy.tsv.gz",
        "OceanDNA":               "ocean_dna_metadata.tsv.gz",
        "SoilSMAG":               "soil_smag_metadata.tsv.gz",
        "TaraEukaryoticSMAG":     "tara_SMAGs_metadata.tsv.gz",
    }
    # sylph-tax also knows three viral taxonomies. Their ranks (realm, viral kingdom,
    # vFAM, ...) do not fit the seven d__..s__ slots (domain ... species) of the sylph
    # tables, so they are refused instead of producing a table phyloseq would read wrongly.
    SYLPH_TAX_VIRAL = {"IMGVR_4.1", "UHGV_default", "UHGV_ictv"}
    # A database paired with the lineage table of another release gets wrong genome
    # names without any warning (sylph-tax looks each genome up by accession, and most
    # accessions are in both releases), so each pair is checked against what the
    # database file name says. The prebuilt GTDB and GlobDB file names carry their
    # release (gtdb-r232-c200-dbv1.syldb, globdb_r232_sylph_c200.syldb), which
    # _sylph_taxonomy_of_database() below reads with or without separators, so a
    # renamed copy (gtdb_r232_c200.syldb) or a release newer than sylph-tax 1.9.2
    # (gtdb-r236-...) is recognised too. The other prebuilt databases are known by a
    # fixed piece of their name, listed here with the built-in taxonomy that names
    # their genomes (as on sylph's download server and the files used to test MetaFlux).
    # A lineage table given by its path is checked the same way when its file name is
    # that of a built-in table (_sylph_taxonomy_of_table). A table of your own under
    # another name cannot be matched to a release here; for GTDB, rule sylph_tables
    # checks the release once more after profiling (38a_sylph_tables.py, section 2b).
    SYLPH_DB_NAME_HINTS = [
        ("fungi-refseq-2025-10-11", "FungiRefSeq-latest"),
        ("fungi-refseq-2024-07-25", "FungiRefSeq-2024-07-25"),
        ("smag-c", "SoilSMAG"),
        ("oceandna-c", "OceanDNA"),
        ("tara-eukmags", "TaraEukaryoticSMAG"),
    ]

    # GlobDB's own download server ships its sylph lineage table, next to
    # globdb_r232_sylph_c200.syldb, under another name than sylph-tax does. It is the
    # same GlobDB r232 taxonomy, so a path to it counts as GlobDB_r232.
    SYLPH_TABLE_ALIASES = {
        "globdb_r232_taxonomy_sylph.tsv.gz": "GlobDB_r232",
    }

    def _sylph_family(taxonomy):
        """The database family of a built-in taxonomy name: GTDB_r232 -> GTDB,
        GlobDB_r226 -> GlobDB, FungiRefSeq-latest -> FungiRefSeq. Every other name is
        its own family. Two releases of one family share genomes and accessions."""
        if taxonomy.startswith("GTDB_"):
            return "GTDB"
        if taxonomy.startswith("GlobDB_"):
            return "GlobDB"
        if taxonomy.startswith("FungiRefSeq-"):
            return "FungiRefSeq"
        return taxonomy

    def _sylph_taxonomy_of_database(file_name):
        """The sylph-tax taxonomy a database file name points to, or None when the name
        shows none (a database you named yourself).
          gtdb-r232-c200-dbv1.syldb or gtdb_r232_c200.syldb -> GTDB_r232
          globdb_r232_sylph_c200.syldb                       -> GlobDB_r232
          gtdb-r236-c200-dbv1.syldb                          -> GTDB_r236, a release
                                     that sylph-tax 1.9.2 has no lineage table for
          fungi-refseq-2025-10-11-c200-dbv1.syldb            -> FungiRefSeq-latest
        The GTDB and GlobDB release is read from 'gtdb' or 'globdb', then an optional
        '-' or '_', then 'r' and the release number."""
        name = file_name.lower()
        gtdb_release = re.search(r"gtdb[-_]?r(\d+)", name)
        if gtdb_release is not None:
            return "GTDB_r" + gtdb_release.group(1)
        globdb_release = re.search(r"globdb[-_]?r(\d+)", name)
        if globdb_release is not None:
            return "GlobDB_r" + globdb_release.group(1)
        for name_piece, taxonomy in SYLPH_DB_NAME_HINTS:
            if name_piece in name:
                return taxonomy
        return None

    def _sylph_family_of_database(file_name):
        """The family a database file name shows (GTDB, GlobDB, FungiRefSeq, ...), or
        None. Used for an entry whose taxonomy is a table of your own, which names no
        family itself. A fungal build of a date not listed above still counts as
        FungiRefSeq."""
        taxonomy = _sylph_taxonomy_of_database(file_name)
        if taxonomy is not None:
            return _sylph_family(taxonomy)
        if re.search(r"fungi[-_]refseq", file_name.lower()):
            return "FungiRefSeq"
        return None

    def _sylph_taxonomy_of_table(file_name):
        """The built-in taxonomy that a lineage table given by its path is a copy of,
        judged by its file name (gtdb_r226_metadata.tsv.gz -> GTDB_r226), or None for a
        table of your own. An unzipped copy (.tsv without .gz) counts as well."""
        if file_name.endswith(".gz"):
            zipped_name = file_name
        else:
            zipped_name = file_name + ".gz"
        for taxonomy, table_file in SYLPH_TAX_FILES.items():
            if table_file == zipped_name:
                return taxonomy
        if zipped_name in SYLPH_TABLE_ALIASES:
            return SYLPH_TABLE_ALIASES[zipped_name]
        return None

    def _sylph_c_in_name(file_name):
        """The sketch density written in a database file name, or None.
        'gtdb-r232-c200-dbv1.syldb' -> 200, 'uhgv_c100_dbv1.syldb' -> 100. The 'c' and
        its number must stand between separators (- _ .) or at either end of the name,
        so a custom name such as 'soil_c50k.syldb' is not misread as c 50."""
        match = re.search(r"(^|[-_.])c(\d+)([-_.]|$)", file_name.lower())
        if match is None:
            return None
        return int(match.group(2))

    def _sylph_prebuilt_name(file_name):
        """True when a database file name looks like one of sylph's prebuilt databases:
        it names a known database (a GTDB or GlobDB release, or a piece listed in
        SYLPH_DB_NAME_HINTS) or carries sylph's own database-version tag ('-dbv1',
        '_dbv1', 'v0.3'). Such names say truthfully what the file is; a name you chose
        yourself may not."""
        name = file_name.lower()
        if _sylph_taxonomy_of_database(name) is not None:
            return True
        if re.search(r"[-_]dbv\d", name):
            return True
        if "v0.3" in name:
            return True
        return False

    # Shown with the messages below: what a working entry looks like, and how to fill
    # the lineage-table folder.
    _SYLPH_EXAMPLE = (
        "  references:\n"
        "    sylph_databases:\n"
        "      - path: /path/to/gtdb-r232-c200-dbv1.syldb\n"
        "        taxonomy: GTDB_r232\n"
        "      - path: /path/to/fungi-refseq-2025-10-11-c200-dbv1.syldb\n"
        "        taxonomy: FungiRefSeq-latest\n"
        "    sylph_taxonomy_dir: /path/to/sylph_taxonomy\n"
        "Databases: https://sylph-docs.github.io/ (page 'Pre-built databases')"
    )
    _SYLPH_DOWNLOAD = (
        "Fill the folder once (about 80 MB) with sylph-tax 1.9.2 (bioconda):\n"
        "  sylph-tax --no-config --taxonomy-dir /path/to/sylph_taxonomy download\n"
        "or download only the tables you need; see docs/shotgun/sylph.md"
    )

    _refs = config.get("references") or {}
    _syl_dbs_cfg = _refs.get("sylph_databases") or []
    _syl_dir_cfg = _refs.get("sylph_taxonomy_dir")
    if _syl_dir_cfg:
        SYLPH_TAXONOMY_DIR = Path(str(_syl_dir_cfg)).resolve()
    else:
        SYLPH_TAXONOMY_DIR = None
    SYLPH_DB_PATHS       = []
    SYLPH_TAXONOMY_ARGS  = []
    SYLPH_TAXONOMY_FILES = []
    SYLPH_DB_LABELS      = []
    _syl_families = []      # the family of each entry whose family is known (family checks)

    if RUN_SYLPH:
        # Why: MetaFlux never downloads sylph databases, and sylph needs at least one.
        if not isinstance(_syl_dbs_cfg, list) or not _syl_dbs_cfg:
            sys.exit(f"[MetaFlux] shotgun.profiler: {PROFILER} needs at least one sylph database in "
                     f"references.sylph_databases (MetaFlux does not download them), e.g.\n"
                     + _SYLPH_EXAMPLE)
        # Why: a mistyped folder would otherwise surface only as a sylph-tax error at the end.
        if SYLPH_TAXONOMY_DIR is not None and not SYLPH_TAXONOMY_DIR.is_dir():
            sys.exit(f"[MetaFlux] references.sylph_taxonomy_dir not found: {SYLPH_TAXONOMY_DIR}\n"
                     + _SYLPH_DOWNLOAD)

        for _i, _entry in enumerate(_syl_dbs_cfg, 1):
            _where = f"references.sylph_databases entry {_i}"

            # Why: each entry must pair a database file with the taxonomy naming its genomes.
            if not isinstance(_entry, dict):
                sys.exit(f"[MetaFlux] {_where} must have a 'path' and a 'taxonomy' (got: {_entry!r}), e.g.\n"
                         + _SYLPH_EXAMPLE)
            # Why: a misspelt key (taxonmy) would otherwise read as a missing taxonomy.
            _unknown_keys = []
            for _key in _entry:
                if _key not in ("path", "taxonomy"):
                    _unknown_keys.append(str(_key))
            if _unknown_keys:
                sys.exit(f"[MetaFlux] {_where} has unknown key(s) {sorted(_unknown_keys)}; "
                         f"the allowed keys are path and taxonomy")
            # Why: there is no database to profile against without a path.
            if not _entry.get("path"):
                sys.exit(f"[MetaFlux] {_where} has no 'path'")
            _db = Path(str(_entry["path"])).resolve()

            # Why: sylph decides what a file is from its name and silently skips anything
            # else. Checked before the file's existence, so that a Kraken2 index FOLDER
            # given here by mistake gets this message.
            if _db.suffix not in (".syldb", ".syl2db"):
                sys.exit(f"[MetaFlux] {_where}: {_db.name} must end in .syldb or .syl2db. sylph decides "
                         f"what a file is from its name and silently skips anything else (a Kraken2 index "
                         f"folder does not work here)")
            # Why: a missing file would only fail inside sylph, after every sample is sketched.
            if not _db.is_file():
                sys.exit(f"[MetaFlux] {_where}: sylph database not found: {_db}")
            # Why: the same genomes loaded twice make sylph's output change from run to run.
            if _db in SYLPH_DB_PATHS:
                sys.exit(f"[MetaFlux] {_where}: {_db.name} is listed twice")

            # Why: a database denser than the sketches (c below 200) makes sylph skip every
            # sample with exit code 0. A prebuilt file name is trusted, so it stops here; a
            # name you chose may use 'c' for something else, so it only warns (rule
            # sylph_profile still stops on sylph's own error if the file really is c 100).
            _db_c = _sylph_c_in_name(_db.name)
            if _db_c is not None and _db_c < SYLPH_SKETCH_C:
                if _sylph_prebuilt_name(_db.name):
                    sys.exit(f"[MetaFlux] {_where}: {_db.name} looks like a database sketched at c {_db_c} "
                             f"(from its file name). MetaFlux sketches reads at c {SYLPH_SKETCH_C}, and sylph "
                             f"skips every sample against a denser database. Use a c 200 or c 1000 build")
                sys.stderr.write(f"[MetaFlux] warning: {_where}: the file name {_db.name} suggests a database "
                                 f"sketched at c {_db_c}. MetaFlux sketches reads at c {SYLPH_SKETCH_C}; if the "
                                 f"database really is denser, sylph_profile will stop with sylph's own error.\n")

            # Why: without a taxonomy, sylph-tax cannot name the genomes of this database.
            _tax = str(_entry.get("taxonomy") or "").strip()
            if not _tax:
                sys.exit(f"[MetaFlux] {_where} ({_db.name}) has no 'taxonomy'. Give the sylph-tax name of its "
                         f"lineage table (GTDB_r232 for gtdb-r232-*.syldb, FungiRefSeq-latest for "
                         f"fungi-refseq-2025-10-11-*.syldb, ...) or the path of your own table")
            # Why: viral ranks do not fit the seven d__..s__ slots of the sylph tables.
            if _tax in SYLPH_TAX_VIRAL:
                sys.exit(f"[MetaFlux] {_where}: taxonomy {_tax} is a viral taxonomy whose ranks (realm, viral "
                         f"kingdom, ...) do not fit the seven d__..s__ slots (domain ... species) of the "
                         f"sylph tables. Viral sylph databases are not supported yet; remove {_db.name} "
                         f"from references.sylph_databases")

            if _tax in SYLPH_TAX_FILES:
                # A built-in name: its table lives in references.sylph_taxonomy_dir.
                # Why: sylph-tax needs the folder to find the table behind a built-in name.
                if SYLPH_TAXONOMY_DIR is None:
                    sys.exit(f"[MetaFlux] {_where}: taxonomy {_tax} is read from references.sylph_taxonomy_dir, "
                             f"which is not set.\n" + _SYLPH_DOWNLOAD)
                # Why: a folder filled by an older sylph-tax can lack a newer table.
                _tax_file = SYLPH_TAXONOMY_DIR / SYLPH_TAX_FILES[_tax]
                if not _tax_file.is_file():
                    sys.exit(f"[MetaFlux] {_where}: taxonomy {_tax} needs {_tax_file.name} in "
                             f"{SYLPH_TAXONOMY_DIR}, and it is not there (folders made before sylph-tax 1.9.2 "
                             f"lack globdb_r232_taxonomy.tsv.gz).\n" + _SYLPH_DOWNLOAD)
                _tax_builtin = _tax       # the built-in taxonomy this entry uses
                _tax_arg = _tax
                _tax_label = _tax
            else:
                # Not a built-in name, so it must be the path of your own lineage table.
                # Why: built-in names are case-sensitive, and gtdb_r232 is an easy typo.
                _near = []
                for _builtin_name in SYLPH_TAX_FILES:
                    if _builtin_name.lower() == _tax.lower():
                        _near.append(_builtin_name)
                if _near:
                    sys.exit(f"[MetaFlux] {_where}: taxonomy {_tax!r} is not a sylph-tax name (names are "
                             f"case-sensitive); did you mean {_near[0]}?")
                # Why: anything else must be an existing table file.
                _tax_file = Path(_tax).resolve()
                if not _tax_file.is_file():
                    sys.exit(f"[MetaFlux] {_where}: taxonomy {_tax!r} is neither a sylph-tax name nor an existing "
                             f"file. Names (capitals matter): {', '.join(sorted(SYLPH_TAX_FILES))}")
                # A table given by its path can still be a copy of a built-in table (a
                # gtdb_r226_metadata.tsv.gz kept elsewhere, say) or GlobDB's own sylph
                # table. Recognised by its file name, it gets the release and family
                # checks below as if its built-in name had been given; sylph-tax is
                # still handed the path. None for a table of your own.
                _tax_builtin = _sylph_taxonomy_of_table(_tax_file.name)
                _tax_arg = str(_tax_file)
                _tax_label = _tax_file.name

            # Why: a taxonomy from another release names genomes wrongly without any warning.
            # Checked when both the database file name and the taxonomy show which
            # database they belong to, and they disagree.
            _db_taxonomy = _sylph_taxonomy_of_database(_db.name)
            if _tax_builtin is not None and _db_taxonomy is not None and _tax_builtin != _db_taxonomy:
                if _tax == _tax_builtin:
                    _tax_shown = _tax
                else:
                    _tax_shown = f"{_tax} (the {_tax_builtin} lineage table)"
                if _db_taxonomy in SYLPH_TAX_FILES:
                    sys.exit(f"[MetaFlux] {_where}: {_db.name} looks like a {_db_taxonomy} database (from its "
                             f"file name) but is paired with taxonomy {_tax_shown}. A taxonomy from another "
                             f"release names genomes wrongly without any warning (tested: GTDB_r226 names on an "
                             f"r232 profile). Use taxonomy: {_db_taxonomy}, or rename the file if its name is wrong")
                # Why: a release newer than the pinned sylph-tax has no built-in table, so
                # any built-in taxonomy would be from another release.
                sys.exit(f"[MetaFlux] {_where}: {_db.name} looks like a {_db_taxonomy} database (from its "
                         f"file name), a release that sylph-tax 1.9.2 (the version MetaFlux pins) has no "
                         f"lineage table for, but it is paired with taxonomy {_tax_shown}. A taxonomy from "
                         f"another release names genomes wrongly without any warning. Give the path of a "
                         f"lineage table made for that release as its taxonomy, or rename the file if its "
                         f"name is wrong")

            # The family of this database (GTDB, GlobDB, ...), for the family checks after
            # the loop: from its taxonomy when that is a built-in one or a copy of one,
            # otherwise from the database file name. A database you named yourself, with
            # a table of your own, has none and is left out of those checks.
            if _tax_builtin is not None:
                _syl_families.append(_sylph_family(_tax_builtin))
            else:
                _db_family = _sylph_family_of_database(_db.name)
                if _db_family is not None:
                    _syl_families.append(_db_family)

            SYLPH_DB_PATHS.append(_db)
            if _tax_arg not in SYLPH_TAXONOMY_ARGS:
                SYLPH_TAXONOMY_ARGS.append(_tax_arg)
                SYLPH_TAXONOMY_FILES.append(_tax_file)
            SYLPH_DB_LABELS.append(f"{_db.name} ({_tax_label})")

        # Why: two databases of one family (two GTDB releases, or the c 200 and c 1000
        # build of one release) share genomes and accessions: sylph then picks one copy of
        # each shared genome arbitrarily, and sylph-tax keeps only one release's names.
        # The families come from the loop above: from each taxonomy, or from the database
        # file name when the taxonomy is a table of your own, so a table given by its
        # path cannot slip past this check.
        _repeated_families = []
        for _family in _syl_families:
            if _syl_families.count(_family) > 1 and _family not in _repeated_families:
                _repeated_families.append(_family)
        if _repeated_families:
            sys.exit(f"[MetaFlux] references.sylph_databases lists two databases of the same family "
                     f"({', '.join(sorted(_repeated_families))}): two releases, or the c200 and c1000 (or "
                     f".syldb and .syl2db) build of one release. They share genomes and accessions, so sylph "
                     f"picks one copy of each shared genome arbitrarily and sylph-tax keeps only one release's "
                     f"names. Keep one (the family is read from each taxonomy or, for a lineage table of "
                     f"your own, from the database file name)")
        # Why: allowed, but a genome held by both databases can be split over two rows of
        # the tables, and the pair costs about twice the memory, so say so. GlobDB holds
        # every genome of the GTDB database, with the same lineage, under its accession
        # without the version suffix; the two copies tie, and which one sylph reports can
        # change from sample to sample.
        if "GTDB" in _syl_families and "GlobDB" in _syl_families:
            sys.stderr.write("[MetaFlux] warning: references.sylph_databases lists both a GTDB and a GlobDB "
                             "database. GlobDB already holds the GTDB genomes, with the same lineage, under "
                             "accessions without the version suffix (GCF_003096595 for GTDB's GCF_003096595.1). "
                             "sylph reports each such genome from one database or the other, and which one can "
                             "change from sample to sample, so one genome can end up on two rows of "
                             "sylph_table.tsv, each with zeros where the other copy was reported: genome-level "
                             "counts, prevalence and diversity are then wrong. The pair also roughly doubles the "
                             "memory (measured 40.5 -> 66.4 GB) without adding species. GlobDB alone covers "
                             "GTDB; keep one.\n")

    # What rule sylph_taxprof passes before the word taxprof: the lineage-table folder,
    # or nothing when every taxonomy is a table path of your own.
    if SYLPH_TAXONOMY_DIR is not None:
        SYLPH_TAXONOMY_DIR_OPTION = ["--taxonomy-dir", str(SYLPH_TAXONOMY_DIR)]
    else:
        SYLPH_TAXONOMY_DIR_OPTION = []

    # Memory for the ONE sylph_profile job. sylph holds every database in memory: with
    # sylph 1.0.0 the measured peak was 1.01-1.11 x the database files' size on disk
    # (26.9-29.7 GB for GTDB r232 c200 + fungi, which are 26.66 GB of files). Ask for
    # 1.2 x their size + 4 GB (34,601 MB for that pair), unless the config sets a number
    # in resources.mem_mb.sylph_profile (e.g. on a cluster). It is worked out only when
    # sylph runs (profiler sylph | both): under profiler: kraken2 there is no
    # sylph_profile rule, so a sylph-only setting must never stop a Kraken2 run.
    SYLPH_PROFILE_MEM_MB = None      # profiler: kraken2 has no sylph_profile rule
    if RUN_SYLPH:
        _resources_cfg = config.get("resources") or {}
        _mem_mb_cfg = _resources_cfg.get("mem_mb") or {}
        _syl_mem_cfg = _mem_mb_cfg.get("sylph_profile")
        if _syl_mem_cfg:
            # Why: a word such as 'auto' would otherwise end in a Python traceback.
            if not str(_syl_mem_cfg).isdigit():
                sys.exit(f"[MetaFlux] resources.mem_mb.sylph_profile must be a whole number of MB, "
                         f"or null to work it out from the database files (got: {_syl_mem_cfg!r})")
            SYLPH_PROFILE_MEM_MB = int(_syl_mem_cfg)
        else:
            _syl_db_bytes = 0
            for _db in SYLPH_DB_PATHS:
                _syl_db_bytes += _db.stat().st_size
            _syl_db_mb = _syl_db_bytes // (1024 * 1024)
            SYLPH_PROFILE_MEM_MB = int(_syl_db_mb * 1.2) + 4096

    # Why: databases listed under profiler: kraken2 are easy to mistake for a sylph run.
    if PROFILER == "kraken2" and _syl_dbs_cfg:
        sys.stderr.write("[MetaFlux] warning: references.sylph_databases is set, but shotgun.profiler is "
                         "kraken2, so sylph will not run.\n")


# ────────────────────── Sample discovery ───────────────────────
# Works out the run's sample list exactly once, here, rather than letting every rule
# file glob fastq_dir independently — SAMPLES, EXTN and FASTQ_NAMING (below) are then
# just plain variables every rule file's `expand(..., sample=SAMPLES)` calls read.
# Discovery differs by source: local files are found by globbing fastq_dir for a
# {sample}_R1 / {sample}_1 naming pattern (_discover_local); a shotgun BioProject fetch
# instead asks NCBI which SRA runs belong to that BioProject (_fetch_bioproject_runs)
# and treats each run accession as a "sample".
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

# ── Sample names sylph-tax would rewrite (profiler sylph | both) ──
# sylph-tax reads sylph's profile with pandas, which guesses each column's type. A
# sample name that looks like a number, a missing value or a true/false word is
# rewritten on the way, with exit code 0. Tested with sylph-tax 1.9.2 / pandas 3.0.6:
# samples '01' and '1' were summed into ONE file (relative abundance 200 %), '001'
# became '1', and 'NA' was dropped. Such names are refused here.
# The number test mirrors pandas, not Python's float(): float() also accepts
# underscores ('12_3' -> 123.0), but pandas keeps '12_3', '2019_05' or '1_000' as text,
# so those names are fine. pandas is not pinned by sylph-tax, so the test errs on the
# safe side (every name that pandas 3 reads as a number is refused, even when it would
# survive the round trip, like '1').
if MODE == "shotgun" and RUN_SYLPH:
    # pandas' default missing-value words (read_csv na_values). Several contain
    # characters BAD_CHARS already refuses; they are kept so the list stays complete.
    _PANDAS_NA = {"", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND",
                  "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null"}

    def _pandas_reads_as_number(name):
        """True for what pandas reads as a number: digits with an optional sign, decimal
        point and exponent ('007', '+5', '5.', '.5', '1e5'), or inf / infinity / nan in any
        capitals with an optional sign ('Infinity'). False for '12_3' or 'S001'."""
        if re.fullmatch(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?", name):
            return True
        if re.fullmatch(r"[+-]?(inf|infinity|nan)", name, flags=re.IGNORECASE):
            return True
        return False

    def _pandas_rewrites(name):
        """True when pandas would change this sample name on the way through sylph-tax."""
        if name in _PANDAS_NA:
            return True
        if name.lower() in ("true", "false"):
            return True
        return _pandas_reads_as_number(name)

    _bad_names = []
    for _s in SAMPLES:
        if _pandas_rewrites(_s):
            _bad_names.append(_s)
    if _bad_names:
        sys.exit(f"[MetaFlux] shotgun.profiler: {PROFILER} cannot use these sample names: {_bad_names}. "
                 f"sylph-tax reads sample names with pandas, which rewrites names that look like numbers "
                 f"or missing values (tested: '001' becomes '1', samples '01' and '1' end up summed in one "
                 f"file, 'NA' is dropped). Rename the FASTQ files so each name contains a letter (e.g. "
                 f"001_R1.fastq.gz -> S001_R1.fastq.gz), or use profiler: kraken2.")

# ── Output folders that already hold results of another profiler ──
# Snakemake decides what to run from the files in out_dir. The trimmed reads both
# profilers start from are temporary (deleted once read), so re-creating them re-runs
# every step downstream of them that already finished. Switching profiler in a
# finished out_dir can therefore do far more than the user expects; these checks stop
# or warn before anything runs. They look only at files in out_dir.
if MODE == "shotgun":

    def _sylph_sketch_file(sample):
        """The kept sylph sketch of one sample (rule sylph_sketch, 38_sylph.smk)."""
        return OUT / "02c.sylph" / "sketches" / f"{sample}.paired.sylsp"

    def _trimmed_r1_file(sample):
        """The temporary trimmed R1 reads of one sample (rule trim_adapters)."""
        return OUT / "01.preprocessing" / f"{sample}_trim_R1.fastq.gz"

    def _kraken2_report_file(sample):
        """The Kraken2 report of one sample (rule kraken2, 35_kraken2.smk)."""
        return OUT / "02.classification" / f"{sample}_report.txt"

    # Why: a sylph-only run in a Kraken2 folder would re-create the trimmed reads and
    # replace the MultiQC report with one that no longer shows Kraken2, next to Kraken2
    # tables from an earlier run.
    _kraken2_results_here = ((OUT / "02.classification").is_dir()
                             or (OUT / "03.abundance" / "otu_table.tsv").exists())
    if PROFILER == "sylph" and _kraken2_results_here:
        sys.exit(f"[MetaFlux] shotgun.profiler is 'sylph', but {OUT} already holds Kraken2 results "
                 f"(02.classification/ or 03.abundance/otu_table.tsv). A sylph-only run here would re-create "
                 f"the trimmed reads and replace the MultiQC report with one that no longer shows Kraken2, next "
                 f"to Kraken2 tables from an earlier run. Use a new out_dir (or profiler: both, if this folder "
                 f"was made with profiler: both).")

    if PROFILER == "both":
        # Why: adding sylph to a FINISHED Kraken2 run needs the deleted trimmed reads,
        # so fastp would run again and then Kraken2, Bracken and every Kraken2 table.
        _would_rerun_kraken2 = []
        for _s in SAMPLES:
            _kraken2_done = _kraken2_report_file(_s).exists()
            _sketched = _sylph_sketch_file(_s).exists()
            _reads_left = _trimmed_r1_file(_s).exists()
            if _kraken2_done and not _sketched and not _reads_left:
                _would_rerun_kraken2.append(_s)
        if _would_rerun_kraken2:
            _shown = ", ".join(_would_rerun_kraken2[:10])
            if len(_would_rerun_kraken2) > 10:
                _shown += ", ..."
            sys.exit(f"[MetaFlux] Kraken2 has already finished in {OUT} for {len(_would_rerun_kraken2)} "
                     f"sample(s), and their trimmed reads were deleted after use (they are temporary files). "
                     f"sylph needs exactly those reads, so Snakemake would run fastp again and then re-run "
                     f"Kraken2, Bracken and every Kraken2 table for these samples. To add sylph to this finished "
                     f"run, use profiler: sylph with a NEW out_dir. Samples: {_shown}")

        # Why: the reverse (adding Kraken2 to a finished sylph run) is allowed, but it also
        # re-runs the sylph steps, so say so.
        _sylph_only_done = []
        for _s in SAMPLES:
            _kraken2_done = _kraken2_report_file(_s).exists()
            _sketched = _sylph_sketch_file(_s).exists()
            _reads_left = _trimmed_r1_file(_s).exists()
            if _sketched and not _kraken2_done and not _reads_left:
                _sylph_only_done.append(_s)
        if _sylph_only_done:
            sys.stderr.write(f"[MetaFlux] warning: adding Kraken2 to a finished sylph run in {OUT}: the trimmed "
                             f"reads are re-created, so sylph_sketch and the sylph steps after it run again too "
                             f"(same results: sketches are made with one thread).\n")

    # Why: sylph tables left from an earlier run are no longer updated under kraken2.
    if PROFILER == "kraken2" and (OUT / "03.abundance" / "sylph_table.tsv").exists():
        sys.stderr.write(f"[MetaFlux] warning: {OUT} holds sylph results from an earlier run "
                         f"(03.abundance/sylph_table.tsv); with profiler: kraken2 they are left as they are "
                         f"and no longer updated.\n")

# The one-line summary every run prints before anything starts. In shotgun mode it also
# names the profiler, so a log shows at a glance which chain made the tables. The
# amplicon line stays word for word as before (the docs quote it).
if MODE == "shotgun":
    _profiler_part = f", profiler={PROFILER}"
else:
    _profiler_part = ""
sys.stderr.write(
    f"[MetaFlux v2.5.0] mode={MODE}{_profiler_part}, input={INPUT_MODE}, "
    f"{len(SAMPLES)} sample(s) in {FASTQ_DIR}\n"
)
# When sylph runs, a second line names its databases (with the taxonomy of each) and
# the memory the one sylph_profile job will ask for, so a log shows both at a glance.
if MODE == "shotgun" and RUN_SYLPH:
    sys.stderr.write(f"[MetaFlux] sylph databases: {', '.join(SYLPH_DB_LABELS)}; "
                     f"sylph_profile asks for {SYLPH_PROFILE_MEM_MB:,} MB\n")


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
# Snakemake builds only what it is asked to build, starting from a target file and
# working backwards through whatever rules produce it. `rule all` in the Snakefile
# (`input: all_targets()`) is what gets asked for on a plain `snakemake` invocation
# with no specific file named on the command line, so the two functions below are, in
# effect, the definition of "a complete run" for each mode: every file listed here
# (and everything upstream of it in the DAG) gets built; anything NOT listed here does
# not, even if some rule in the workflow is capable of producing it. Adding a new
# final output to the pipeline means adding its path to the relevant list here, or
# `snakemake` will happily finish without ever having produced it.
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
        # Every ASV that reached the taxonomy filter, removed ones included (same rule
        # as the three tables above, so it adds no job).
        OUT / "6.taxonomy" / "asv_table_unfiltered.txt",
        # Combined read-tracking table
        OUT / "stats" / "read_tracking.txt",
        # MultiQC aggregate
        OUT / "multiqc" / "multiqc_report.html",
    ]
    # Optional 16S phylogeny (stage 85). ONE target is registered, and it is the QC
    # report rather than the tree — which needs explaining, because it looks like the
    # tree has been left unregistered by mistake.
    #
    # The reason is the spec's requirement that a run with fewer than 4 eligible ASVs
    # produce a documented skip rather than a crash: 3 sequences cannot make a
    # meaningful tree. That decision depends on a number nobody knows until taxonomy
    # has finished, but this list is built at PARSE time, before anything has run. Ask
    # for the tree here and a legitimate skip ends the run with a missing-output error.
    #
    # So phylo_input is a Snakemake *checkpoint* — a rule after which Snakemake re-reads
    # the situation and works out the rest of the graph from what that rule actually
    # produced. phylo_qc's input list is then a function (see 85_phylogeny.smk): with 4
    # or more ASVs it asks for the alignment, tree, and params file, which is what pulls
    # that whole chain into the run; with fewer, it asks for none of them and records
    # the skip instead. Every file in the spec's 7.phylogeny/ layout is still produced
    # on any normal run — it arrives as a dependency of this one target rather than as
    # a target itself.
    #
    # Practical consequence worth knowing: because the tree is not itself a target,
    # deleting it while phylogeny_qc.json survives will not trigger a rebuild on its
    # own. Delete the QC report too, or re-run with `snakemake -R phylo_export`.
    if PHYLO_ENABLED:
        targets += [
            OUT / "stats" / "phylogeny" / "phylogeny_qc.json",
        ]
    return targets


def _shotgun_targets():
    """Default target list for mode=shotgun — the files `rule all` requests for a shotgun run."""
    # The fastp report of every sample: the read cleaning runs whichever profiler is
    # chosen.
    out = [
        *expand(OUT / "01.preprocessing" / "{s}_fastp.json",          s=SAMPLES),
    ]
    # Kraken2 + Bracken (profiler kraken2 or both). Under kraken2 the whole list comes
    # out exactly as before, in the same order, so a finished Kraken2 run sees no change.
    if RUN_KRAKEN2:
        out += [
            *expand(OUT / "02.classification" / "{s}_report.txt",         s=SAMPLES),
            *expand(OUT / "02.classification" / "{s}_output.txt",         s=SAMPLES),
            *expand(OUT / "02.classification" / "{s}_R1.fastq.gz",        s=SAMPLES),
            *expand(OUT / "02.classification" / "{s}_R2.fastq.gz",        s=SAMPLES),
            # k-mer evidence: requested unconditionally, so the evidence table and plot
            # exist on every shotgun run even when the gate is switched off. Asking for
            # these also builds the gated report (same rule, three outputs), which is
            # what makes "diff it before you trust it" possible.
            *expand(OUT / "02b.evidence" / "{s}_species_evidence.tsv",    s=SAMPLES),
            *expand(OUT / "02b.evidence" / "{s}_report_gated.txt",        s=SAMPLES),
            *expand(OUT / "stats" / "kmer_evidence" / "{s}_evidence.png", s=SAMPLES),
            *expand(OUT / "03.abundance" / "{s}_report.txt",              s=SAMPLES),
            *expand(OUT / "03.abundance" / "{s}_output.txt",              s=SAMPLES),
            OUT / "03.abundance" / "otu_table.tsv",
        ]
    # sylph (profiler sylph or both): the raw sylph profile, the per-genome tables made
    # from it, and the per-sample summary that MultiQC shows.
    if RUN_SYLPH:
        out += [
            OUT / "03.abundance" / "sylph_profile.tsv",
            OUT / "03.abundance" / "sylph_table.tsv",
            OUT / "03.abundance" / "sylph_table.genome_copies.tsv",
            OUT / "03.abundance" / "sylph_table.biom",
            OUT / "03.abundance" / "sylph_genomes.tsv",
            OUT / "stats" / "sylph" / "sylph_summary.tsv",
        ]
    # Shared by every profiler: the read-tracking table (written by
    # aggregate_read_counts in 80_stats.smk when Kraken2 runs, by sylph_read_tracking
    # in 38_sylph.smk under profiler: sylph), the MultiQC report and the sample manifest.
    out += [
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
