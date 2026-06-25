#!/usr/bin/env python3
"""Aggregate per-stage read counts into a single wide-format tracking table.

Columns (left to right):
  raw → nophix → stripped → filtered → denoised → merged → non_chimeric
  [→ post_extraction] → post_length_filter → post_taxonomy_filter

Pipeline order: dada_seqtab → target_extract (optional) → dada_length_filter
                → assign_taxonomy. When extraction is disabled, post_extraction
                is omitted from the table.

Count files produced by preprocess.smk use the format:
  /path/to/{sample}_R1.fastq[.gz] : N
  ...
  Total read count : N

DADA2 read.counts (from dada_seqtab.R) is a TSV with row-names=samples and
columns: stripped, filtered, denoised, merged, non_chimeric.

Seqtab files (TSV, row-names=samples, cols=ASV_IDs) → rowSums = reads at that stage.
asv_table.txt (TSV, row-names=ASV_IDs, cols=[samples..., taxonomy]) → colSums = reads per sample.
"""
import sys
from pathlib import Path


def log(msg: str, fh) -> None:
    print(msg, file=fh, flush=True)
    print(msg, file=sys.stderr, flush=True)


def parse_stage_counts(path: Path) -> dict[str, int]:
    """Parse preprocess.smk read-count files → {sample: count}."""
    counts: dict[str, int] = {}
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("Total") or " : " not in line:
                continue
            sample, n_str = line.split(" : ", 1)
            counts[sample.strip()] = int(n_str.strip())
    return counts


def _unq(s: str) -> str:
    """Strip surrounding double quotes that R's write.table adds to names."""
    return s.strip('"')


def parse_dada2_counts(path: Path) -> dict[str, dict[str, int]]:
    """Parse dada_seqtab.R read.counts TSV → {sample: {stage: count}}."""
    result: dict[str, dict[str, int]] = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # header[0] is empty; header[1:] are stage names (may be R-quoted)
        stages = [_unq(s) for s in header[1:]]
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            sample = _unq(cols[0])
            result[sample] = {s: int(cols[i + 1]) for i, s in enumerate(stages)}
    return result


def parse_seqtab_rowsums(path: Path) -> dict[str, int]:
    """Sum each row (sample) of a seqtab_head_names.txt → {sample: total_reads}."""
    result: dict[str, int] = {}
    with path.open() as fh:
        fh.readline()  # skip header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            sample = _unq(cols[0])
            result[sample] = sum(int(v) for v in cols[1:] if v)
    return result


def parse_asv_table_colsums(path: Path) -> dict[str, int]:
    """Sum each sample column of asv_table.txt (rows=ASVs, cols=[samples, taxonomy]).
    The last column is 'taxonomy' — excluded from the sum."""
    result: dict[str, int] = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # header[0] empty, header[1:] = sample cols + "taxonomy" (may be R-quoted)
        col_names = [_unq(h) for h in header[1:]]
        sample_cols = col_names[:-1]  # drop 'taxonomy'
        for name in sample_cols:
            result[name] = 0
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            for i, name in enumerate(sample_cols):
                result[name] += int(cols[i + 1]) if cols[i + 1] else 0
    return result


# ── Main ──────────────────────────────────────────────────────────────────────
sm = snakemake  # noqa: F821

log_path = Path(sm.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
log_fh = log_path.open("w")

samples = list(sm.params.samples)
extraction_enabled = bool(sm.params.extraction_enabled)
remove_phix = bool(sm.params.remove_phix)

raw_counts      = parse_stage_counts(Path(sm.input.raw_counts))
nophix_counts   = parse_stage_counts(Path(sm.input.nophix_counts)) if remove_phix else {}
stripped_counts = parse_stage_counts(Path(sm.input.stripped_counts))
dada2_counts    = parse_dada2_counts(Path(sm.input.dada2_counts))
lenfilt_counts  = parse_seqtab_rowsums(Path(sm.input.seqtab_lenfilt))
tax_counts      = parse_asv_table_colsums(Path(sm.input.asv_table))

if extraction_enabled and sm.input.seqtab_extracted:
    extracted_counts = parse_seqtab_rowsums(
        Path(sm.input.seqtab_extracted[0] if isinstance(sm.input.seqtab_extracted, list)
             else sm.input.seqtab_extracted)
    )
else:
    extracted_counts = {}

stages = ["raw"]
if remove_phix:
    stages.append("nophix")
stages += ["stripped", "filtered", "denoised", "merged", "non_chimeric"]
if extraction_enabled:
    stages.append("post_extraction")
stages += ["post_length_filter", "post_taxonomy_filter"]

out_path = Path(sm.output.counts)
out_path.parent.mkdir(parents=True, exist_ok=True)

with out_path.open("w") as out:
    out.write("sample\t" + "\t".join(stages) + "\n")
    for s in samples:
        d2 = dada2_counts.get(s, {})
        row = [s, str(raw_counts.get(s, "NA"))]
        if remove_phix:
            row.append(str(nophix_counts.get(s, "NA")))
        row += [
            str(stripped_counts.get(s, "NA")),
            str(d2.get("filtered", "NA")),
            str(d2.get("denoised", "NA")),
            str(d2.get("merged", "NA")),
            str(d2.get("non_chimeric", "NA")),
        ]
        if extraction_enabled:
            row.append(str(extracted_counts.get(s, "NA")))
        row.append(str(lenfilt_counts.get(s, "NA")))
        row.append(str(tax_counts.get(s, "NA")))
        out.write("\t".join(row) + "\n")

log(f"[aggregate_read_counts] Written {out_path} for {len(samples)} sample(s)", log_fh)
log_fh.close()
