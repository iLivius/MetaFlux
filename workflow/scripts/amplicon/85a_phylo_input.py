#!/usr/bin/env python3
"""Extract the eligible ASV set for the phylogeny module.

This is the first step of the optional 16S phylogeny (see 85_phylogeny.smk). Its job
is small but load-bearing: turn MetaFlux's final taxonomy tables into a plain FASTA
the aligner can read, and prove while doing so that the ASVs going into the tree are
exactly the ASVs in the abundance table the user will pair it with in R.

WHY IT READS 6.taxonomy AND NOT 5.dada2
---------------------------------------
MetaFlux's contaminant filter (chloroplast, mitochondria, wrong-domain, off-target)
runs inside stage 80, when taxonomy is assigned. The stage-60 FASTA
(5.dada2/seqs_lenfilt.fasta) therefore still contains every off-target ASV the filter
would remove. Those sequences are only distantly related to real 16S, so in a de novo
tree they hang off long branches — and long branches are precisely what distorts
Faith's PD (an absolute sum of branch length) and unweighted UniFrac. Reading the
post-filter set removes that failure mode by construction rather than trying to
detect it afterwards.

TIP LABELS
----------
Headers are bare ASV IDs (ASV_1, ASV_2, ...), never the taxonomy string. Two
practical reasons beyond tidiness: FastTree silently truncates a sequence name at the
first space unless names are quoted, and both FastTree and RAxML-NG reject the
characters ":,()" that a lineage string is full of. A tree keyed on ASV IDs joins
cleanly to the abundance table in R, where the user can attach taxonomy for display.

Note for the docs and for anyone comparing runs: ASV IDs are assigned by decreasing
abundance when the sequence table is built, so they are stable WITHIN a run but not
across runs. ASV_7 in two different runs is not the same organism.

Inputs (see rule phylo_input in 85_phylogeny.smk)
  taxon_seq_table : 6.taxonomy/taxon_seq_table.txt, from assign_taxonomy (80a or
                    80c). Rows are ASV IDs; the last column is the sequence.
  asv_table       : 6.taxonomy/asv_table.txt, same rule. Rows are ASV IDs, columns
                    are per-sample counts plus taxonomy. Used here purely as the
                    authority on which ASVs survived the contaminant filter.

Outputs
  seqs       : 7.phylogeny/asv_16s.fasta — read next by phylo_align (MAFFT).
  input_json : 7.phylogeny/input.json — the eligible count (which decides whether a
               tree is attempted at all), the ASV ID list, per-ASV lengths, and
               checksums of both source tables and the FASTA. Read by phylo_qc's
               input function to make the build-or-skip decision, by phylo_export as
               the reference ID set its tip-accounting gate compares against, and
               folded into phylogeny.params.json as provenance.

Both tables are written by R's write.table (or by 80c mimicking it), so every
character field arrives wrapped in double quotes while integer counts do not. The
parsing here strips those quotes, exactly as 80b/80c already do.
"""
import hashlib
import json
import sys
from pathlib import Path


def unquote(value: str) -> str:
    """Strip R write.table's surrounding double quotes from one field."""
    return value.strip().strip('"')


def sha256_of(path: Path) -> str:
    """Checksum a file, so phylogeny.params.json can record exactly which inputs
    produced this tree. Read in chunks — a taxon_seq_table for a large run holds
    every ASV sequence and is not something to load whole just to hash."""
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class InputError(Exception):
    """A problem with the taxonomy tables that must stop the run.

    Raised rather than calling sys.exit() directly inside the parsing helpers, so that
    main() can write the message into the rule's own log file before exiting. Snakemake
    does not populate the log of a `script:` rule — the script has to — so a bare
    sys.exit() deep in a parser puts the diagnosis on stderr only, and the log file the
    user is directed to ends up nearly empty.
    """


def read_taxon_seq_table(path: Path) -> dict[str, str]:
    """Read ASV_ID -> sequence from taxon_seq_table.txt.

    Layout written by 80a/80c: a header row naming the rank columns and then
    "sequence", followed by one row per ASV whose first field is the ASV ID and whose
    LAST field is the sequence. The sequence is located by position rather than by
    counting rank columns, because the number of ranks differs between markers (7 for
    16S, 9 for PR2/18S) — position is the stable property.
    """
    id_to_seq: dict[str, str] = {}

    with path.open() as fh:
        header = fh.readline()
        if not header:
            raise InputError(f"{path} is empty — taxonomy did not produce a table.")

        for line_number, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise InputError(
                    f"{path} line {line_number} has {len(fields)} column(s); "
                    "expected at least an ASV ID and a sequence."
                )

            asv_id   = unquote(fields[0])
            sequence = unquote(fields[-1]).upper()

            if not asv_id:
                raise InputError(f"{path} line {line_number} has an empty ASV ID.")
            if not sequence:
                raise InputError(
                    f"{path} line {line_number} ({asv_id}) has an empty sequence. The "
                    "phylogeny module cannot align a missing sequence."
                )
            if asv_id in id_to_seq:
                raise InputError(
                    f"{path} lists {asv_id} more than once. Each ASV must appear exactly "
                    "once for the tree's tips to match the abundance table."
                )

            id_to_seq[asv_id] = sequence

    return id_to_seq


def read_asv_ids(path: Path) -> list[str]:
    """Read the ASV IDs, in order, from the first column of asv_table.txt."""
    asv_ids: list[str] = []

    with path.open() as fh:
        header = fh.readline()
        if not header:
            raise InputError(f"{path} is empty — taxonomy did not produce a table.")

        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            asv_id = unquote(line.split("\t")[0])
            if asv_id:
                asv_ids.append(asv_id)

    return asv_ids


def remove_stale_outputs(phylo_dir: Path, log) -> None:
    """Clear a previous run's tree and alignment when this run is skipping the build.

    Needed because of how the skip works: with fewer than four ASVs, phylo_align,
    phylo_tree and phylo_export are never put in the DAG at all, and Snakemake only
    deletes outputs of jobs it actually runs. So an earlier successful run's tree would
    survive untouched next to a QC report saying "skipped" — a tree built from a
    different, larger ASV set, sitting in the directory as though it were current. For a
    workflow whose users publish from its output, that is the worst kind of leftover.

    Idempotent: a second consecutive skipped run finds nothing to remove.
    """
    stale_files = [
        phylo_dir / "asv_16s.aln.fasta",
        phylo_dir / "asv_16s.unrooted.nwk",
        phylo_dir / "phylogeny.params.json",
        phylo_dir / "aln_run.json",
    ]
    removed = []
    for stale in stale_files:
        if stale.exists():
            stale.unlink()
            removed.append(stale.name)

    # Backend directories too — each holds a complete tree from the earlier run.
    for backend_name in ("iqtree", "fasttree", "raxml-ng"):
        backend_dir = phylo_dir / backend_name
        if backend_dir.is_dir():
            for artifact in sorted(backend_dir.iterdir()):
                if artifact.is_file():
                    artifact.unlink()
            backend_dir.rmdir()
            removed.append(backend_name + "/")

    if removed:
        log(f"[phylo_input] Removed {len(removed)} stale output(s) from an earlier run "
            f"that DID build a tree: {', '.join(removed)}. They describe a different ASV "
            "set and would be misleading beside a skipped run.")


def main() -> int:
    sm = snakemake  # noqa: F821

    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = log_path.open("w")

    def log(message: str) -> None:
        print(message, file=log_fh, flush=True)
        print(message, file=sys.stderr, flush=True)

    taxon_seq_table = Path(sm.input.taxon_seq_table)
    asv_table       = Path(sm.input.asv_table)
    seqs_out        = Path(sm.output.seqs)
    input_json_out  = Path(sm.output.input_json)
    min_asvs        = int(sm.params.min_asvs)

    log(f"[phylo_input] Reading the post-taxonomy ASV set from {taxon_seq_table.parent}")

    # Every fatal problem below is reported through log() — which writes to both the
    # rule's log file and stderr — and the handle is closed before exiting, so the log
    # a failing run leaves behind actually contains the reason it failed.
    try:
        id_to_seq = read_taxon_seq_table(taxon_seq_table)
        table_ids = read_asv_ids(asv_table)
    except InputError as exc:
        log(f"[phylo_input] ERROR: {exc}")
        log_fh.close()
        return 1

    log(f"[phylo_input] taxon_seq_table.txt: {len(id_to_seq)} ASVs with sequences")
    log(f"[phylo_input] asv_table.txt:       {len(table_ids)} ASVs")

    # The two tables are written by the same rule from the same filtered ASV list, so
    # they must agree. If they ever do not, something upstream is inconsistent, and
    # building a tree from whichever subset happens to overlap would hand the user a
    # tree that silently disagrees with their abundance table. Stop instead.
    seq_id_set   = set(id_to_seq)
    table_id_set = set(table_ids)
    if seq_id_set != table_id_set:
        only_in_seq   = sorted(seq_id_set - table_id_set)[:10]
        only_in_table = sorted(table_id_set - seq_id_set)[:10]
        log("[phylo_input] ERROR: the two taxonomy tables describe different ASV sets.")
        if only_in_seq:
            log(f"[phylo_input]   only in taxon_seq_table.txt (up to 10): {only_in_seq}")
        if only_in_table:
            log(f"[phylo_input]   only in asv_table.txt (up to 10):       {only_in_table}")
        log("[phylo_input] Both tables come from assign_taxonomy, so this points at a "
            "corrupt or partially-written 6.taxonomy/ — re-run the taxonomy step.")
        log_fh.close()
        return 1

    if len(table_ids) != len(table_id_set):
        log(f"[phylo_input] ERROR: {asv_table} lists at least one ASV ID more than once. "
            "Each ASV must appear exactly once.")
        log_fh.close()
        return 1

    # Write the FASTA in asv_table.txt's own row order (decreasing abundance), so the
    # file reads naturally and MAFFT's --inputorder keeps that order in the alignment.
    seqs_out.parent.mkdir(parents=True, exist_ok=True)
    with seqs_out.open("w") as fh:
        for asv_id in table_ids:
            fh.write(f">{asv_id}\n{id_to_seq[asv_id]}\n")

    lengths = {asv_id: len(id_to_seq[asv_id]) for asv_id in table_ids}
    n_eligible = len(table_ids)

    log(f"[phylo_input] Wrote {n_eligible} ASVs to {seqs_out.name} "
        f"(lengths {min(lengths.values(), default=0)}-{max(lengths.values(), default=0)} bp)")

    # Below four tips there is only one possible unrooted topology, so there is nothing
    # to infer. This is recorded, not raised: phylo_qc's input function reads
    # n_eligible and simply does not ask for a tree, and 85g writes down why. See the
    # checkpoint note in 85_phylogeny.smk.
    if n_eligible < min_asvs:
        log(f"[phylo_input] NOTE: only {n_eligible} eligible ASV(s), below the minimum of "
            f"{min_asvs}. No tree will be built; stats/phylogeny/phylogeny_qc.json will "
            "record the skip.")
        remove_stale_outputs(seqs_out.parent, log)

    record = {
        "n_eligible": n_eligible,
        "min_asvs": min_asvs,
        "asv_ids": table_ids,
        "asv_lengths": lengths,
        "checksums": {
            "taxon_seq_table.txt": sha256_of(taxon_seq_table),
            "asv_table.txt": sha256_of(asv_table),
            "asv_16s.fasta": sha256_of(seqs_out),
        },
    }
    input_json_out.parent.mkdir(parents=True, exist_ok=True)
    input_json_out.write_text(json.dumps(record, indent=2) + "\n")

    log(f"[phylo_input] DONE. Recorded {n_eligible} eligible ASVs in {input_json_out.name}")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
