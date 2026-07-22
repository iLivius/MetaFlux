#!/usr/bin/env python3
"""Finalize the shotgun OTU table produced by kraken-biom.

Three things, all rendered from one filtered BIOM table so the .biom and .tsv
outputs never disagree:

  1. Strip the "_report" suffix kraken-biom copies from the Bracken report
     filenames onto the sample ids.
  2. Normalize the taxonomy lineage to the MetaFlux standard shared with the
     amplicon tables: the species slot becomes a Genus+species binomial
     (kraken-biom emits epithet-only, e.g. 's__denitrificans'), and the TSV
     joins ranks with ';' (no space).
  3. Apply the optional rank-aware taxon keep/discard filter
     (shotgun.taxonomy_filter) — the same segment-matching machinery as the
     amplicon contaminant filter. keep runs first, then discard.

Taxonomy in the BIOM is stored as a list of rank strings
(['k__Bacteria', ..., 'g__Achromobacter', 's__sp. MFA1 R4']), so filtering and
normalization operate on structured segments — no fragile string parsing.

FUTURE (known limitation): the species binomial is rebuilt from kraken-biom's
already-genus-stripped field, which is lossy — so viral strain codes (e.g.
's__TP778L' under g__Brussowvirus) are left genus-less rather than risk
mis-attaching a genus. Bacterial/archaeal binomials are recovered correctly.
For a truly complete binomial on every taxon (viruses included), source the full
NCBI names from the per-sample Bracken reports instead of kraken-biom's lineage.
"""
import biom
from biom.util import biom_open

sm = snakemake  # noqa: F821

raw_biom_path  = str(sm.input.raw_biom)
out_biom       = str(sm.output.biom)
out_tsv        = str(sm.output.tsv)
filter_enabled = bool(sm.params.filter_enabled)
filter_keep    = list(sm.params.filter_keep) or []
filter_discard = list(sm.params.filter_discard) or []

log = open(sm.log[0], "w")


def logmsg(msg: str) -> None:
    log.write(f"[finalize_otu_table] {msg}\n")
    log.flush()


table = biom.load_table(raw_biom_path)
logmsg(f"loaded {table.shape[0]} OTUs x {table.shape[1]} samples")

# ── 1. strip the "_report" suffix from sample ids ──────────────────────────
def _clean_sample(sid: str) -> str:
    s = str(sid)
    return s[: -len("_report")] if s.endswith("_report") else s

table.update_ids({sid: _clean_sample(sid) for sid in table.ids("sample")},
                 axis="sample", inplace=True)

# ── 2. normalize taxonomy (rebuild the species binomial where it is safe) ──
# kraken-biom strips the genus from the species field ONLY when the NCBI species
# name starts with the genus node, leaving a bare lowercase epithet (or "sp. …").
# When it does NOT start with the genus (reclassified synonyms, "Candidatus …",
# phage strain codes) it keeps the full name verbatim. A bacterial specific
# epithet is lowercase by nomenclatural convention, so we re-prepend the genus
# ONLY when the remaining field begins lowercase — that recovers the binomial for
# ordinary taxa (s__thermophilus -> s__Streptococcus thermophilus) while never
# doubling a genus onto an already-full or strain-coded name
# (s__Candidatus Coxiella mudrowiae, s__TP778L stay untouched).
def _normalize(tax_list) -> list:
    genus = None
    for seg in tax_list:
        seg = str(seg).strip()
        if seg.startswith("g__"):
            genus = seg[3:]
    out = []
    for seg in tax_list:
        seg = str(seg).strip()
        if seg.startswith("s__") and genus:
            epithet = seg[3:]
            if epithet[:1].islower() and not epithet.startswith(genus + " "):
                seg = "s__" + genus + " " + epithet
        out.append(seg)
    return out


norm = {}
for oid in table.ids("observation"):
    md = table.metadata(oid, "observation")
    tax = list(md["taxonomy"]) if md and md.get("taxonomy") else []
    norm[oid] = _normalize(tax)
table.add_metadata({oid: {"taxonomy": t} for oid, t in norm.items()},
                   axis="observation")

# ── 3. optional keep/discard filter (rank-aware segment matching) ──────────
def _passes(tax_list) -> bool:
    segs = {str(s).strip() for s in tax_list}
    if filter_keep and not any(tok in segs for tok in filter_keep):
        return False
    if filter_discard and any(tok in segs for tok in filter_discard):
        return False
    return True


if filter_enabled and (filter_keep or filter_discard):
    n_before = table.shape[0]
    keep_ids = [oid for oid in table.ids("observation") if _passes(norm[oid])]
    if not keep_ids:
        logmsg("WARNING: keep/discard filter removed every OTU — writing an empty table")
    table = table.filter(keep_ids, axis="observation", inplace=False)
    logmsg(f"filter ON  (keep={filter_keep or '[]'}, discard={filter_discard or '[]'}): "
           f"{table.shape[0]} / {n_before} OTUs retained")
else:
    logmsg("filter OFF (shotgun.taxonomy_filter disabled or empty)")

# ── 4. write the final BIOM ────────────────────────────────────────────────
with biom_open(out_biom, "w") as fh:
    table.to_hdf5(fh, "MetaFlux finalize_otu_table")

# ── 5. render the matching TSV (';'-joined, no space) ──────────────────────
sample_ids = [str(s) for s in table.ids("sample")]
with open(out_tsv, "w") as fh:
    fh.write("# Constructed from biom file\n")
    fh.write("#OTU ID\t" + "\t".join(sample_ids) + "\ttaxonomy\n")
    for oid in table.ids("observation"):
        counts = table.data(oid, "observation")   # aligned to sample_ids order
        md = table.metadata(oid, "observation")
        tax = md["taxonomy"] if md and md.get("taxonomy") else []
        tax_str = ";".join(str(s).strip() for s in tax)
        fh.write("\t".join([str(oid)] + [f"{v}" for v in counts] + [tax_str]) + "\n")

logmsg(f"wrote {out_biom} and {out_tsv}")
log.close()
