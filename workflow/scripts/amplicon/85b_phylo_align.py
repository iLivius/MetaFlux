#!/usr/bin/env python3
"""Align the eligible ASVs with MAFFT.

Second step of the optional 16S phylogeny (see 85_phylogeny.smk). Takes the FASTA of
post-taxonomy ASVs written by phylo_input and produces the multiple sequence alignment
every tree backend reads. The aligner is fixed — MAFFT — and not user-selectable; what
IS selectable is which of MAFFT's two relevant strategies runs.

THE TWO TIERS, AND WHY MetaFlux PICKS RATHER THAN MAFFT
-------------------------------------------------------
MAFFT ships an --auto mode, and this module deliberately does not use it. Its decision
thresholds are not published anywhere in MAFFT's documentation — they live only in the
source — and they switch on dataset size. That means adding a single ASV can cross a
boundary and silently change the alignment algorithm between two runs of the same
study, which is the opposite of what a workflow promising reproducibility should do.

Instead MetaFlux applies its own rule, at MAFFT's own *documented* ~200-sequence
ceiling for its accuracy strategies, and writes down which one it used:

    < 200 ASVs   L-INS-i   --localpair --maxiterate 1000 --threadit 0
    >= 200 ASVs  FFT-NS-2  --retree 2 --maxiterate 0

L-INS-i is the more accurate method but MAFFT documents it as applicable only up to
roughly 200 sequences; past that FFT-NS-2 is the appropriate progressive method, and
for same-region 16S reads — which are easy to align — it is entirely adequate.

--threadit 0 on the L-INS-i tier is required, not cosmetic. --maxiterate 1000 turns on
iterative refinement, which is MAFFT's one documented source of different results
between runs when using multiple threads. FFT-NS-2 never refines, so the flag is
unnecessary there.

FLAGS DELIBERATELY NOT USED
---------------------------
--adjustdirection : renames any sequence it reverse-complements with an "_R_" prefix.
                    That would break the ASV-ID match against the abundance table
                    without any error being raised. MetaFlux settles orientation
                    upstream, at primer trimming.
--thread -1       : takes every physical core on the machine, ignoring what Snakemake
                    allocated to this rule.
--auto            : see above.

Inputs (see rule phylo_align in 85_phylogeny.smk)
  seqs : 7.phylogeny/asv_16s.fasta, from phylo_input. Bare ASV_N headers.

Outputs
  aln      : 7.phylogeny/asv_16s.aln.fasta — read next by phylo_tree, and kept as a
             deliverable so it can be inspected or reused without re-running.
  run_json : 7.phylogeny/aln_run.json — resolved strategy, MAFFT version, thread
             count, and the exact command line. Folded into phylogeny.params.json by
             phylo_export; the recording is what makes a size-dependent choice
             defensible rather than hidden.
"""
import json
import shlex
import subprocess
import sys
from pathlib import Path


def parse_fasta_ids(path: Path) -> list[str]:
    """Return the sequence IDs of a FASTA, in file order.

    The ID is everything up to the first whitespace, which is also how every tree
    program reads it — so this sees the names the way the backend will.
    """
    ids: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0] if line[1:].strip() else "")
    return ids


def mafft_version(log) -> str:
    """Ask MAFFT its version, for the provenance record.

    MAFFT prints its version to stderr and exits non-zero for --version, so both
    streams are captured and the exit code is ignored. A failure here must not stop
    the run: not knowing the version is a gap in the record, not a wrong tree.
    """
    try:
        completed = subprocess.run(
            ["mafft", "--version"],
            capture_output=True, text=True, check=False,
        )
        version = (completed.stderr or completed.stdout).strip().splitlines()
        return version[0] if version else "unknown"
    except Exception as exc:                                    # noqa: BLE001
        log(f"[phylo_align] WARNING: could not read MAFFT's version ({exc})")
        return "unknown"


def resolve_strategy(requested: str, n_seqs: int, threshold: int, log) -> tuple[str, list[str]]:
    """Decide which MAFFT strategy runs, and return (name, flags).

    `requested` is phylogeny.aligner_strategy: "auto" applies the size rule described
    in the module docstring; "linsi" or "fftns2" force one regardless of size (useful
    for making two runs of different sizes directly comparable).
    """
    linsi_flags  = ["--localpair", "--maxiterate", "1000", "--threadit", "0"]
    fftns2_flags = ["--retree", "2", "--maxiterate", "0"]

    if requested == "linsi":
        log(f"[phylo_align] Strategy L-INS-i (forced by config; {n_seqs} ASVs)")
        if n_seqs >= threshold:
            log(f"[phylo_align] NOTE: MAFFT documents L-INS-i as applicable up to about "
                f"{threshold} sequences; running it on {n_seqs} may be slow.")
        return "linsi", linsi_flags

    if requested == "fftns2":
        log(f"[phylo_align] Strategy FFT-NS-2 (forced by config; {n_seqs} ASVs)")
        return "fftns2", fftns2_flags

    # auto — MetaFlux's own documented rule.
    if n_seqs < threshold:
        log(f"[phylo_align] Strategy L-INS-i ({n_seqs} ASVs, below the {threshold}-sequence "
            "ceiling MAFFT documents for its accuracy strategies)")
        return "linsi", linsi_flags

    log(f"[phylo_align] Strategy FFT-NS-2 ({n_seqs} ASVs, at or above the {threshold}-sequence "
        "ceiling for the accuracy strategies)")
    return "fftns2", fftns2_flags


def main() -> int:
    sm = snakemake  # noqa: F821

    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = log_path.open("w")

    def log(message: str) -> None:
        print(message, file=log_fh, flush=True)
        print(message, file=sys.stderr, flush=True)

    seqs_in       = Path(sm.input.seqs)
    aln_out       = Path(sm.output.aln)
    run_json_out  = Path(sm.output.run_json)

    requested_strategy = str(sm.params.strategy)
    extra_args         = str(sm.params.extra_args or "").strip()
    linsi_threshold    = int(sm.params.linsi_threshold)
    threads            = int(sm.threads)

    input_ids = parse_fasta_ids(seqs_in)
    n_seqs = len(input_ids)
    log(f"[phylo_align] {n_seqs} ASVs to align, {threads} thread(s)")

    strategy_name, strategy_flags = resolve_strategy(
        requested_strategy, n_seqs, linsi_threshold, log
    )

    # --nuc is explicit so the run never depends on MAFFT's undocumented guess at
    # whether the input is nucleotide or protein. --preservecase keeps sequences as
    # written, and --inputorder makes the alignment come back in the same order as the
    # input rather than MAFFT's internal guide-tree order, which keeps the file
    # readable next to the abundance table.
    command = [
        "mafft",
        "--thread", str(threads),
        "--nuc",
        "--preservecase",
        "--inputorder",
        *strategy_flags,
    ]
    if extra_args:
        # Advanced, unvalidated passthrough — recorded so that whatever it did is at
        # least visible in the run record afterwards.
        command += shlex.split(extra_args)
        log(f"[phylo_align] Appending extra_args.mafft: {extra_args}")
    command.append(str(seqs_in))

    log(f"[phylo_align] Running: {' '.join(command)}")

    aln_out.parent.mkdir(parents=True, exist_ok=True)
    # MAFFT writes the alignment to stdout and its progress chatter to stderr; the
    # chatter goes into this rule's log so a failure is diagnosable.
    with aln_out.open("w") as aln_fh:
        completed = subprocess.run(command, stdout=aln_fh, stderr=subprocess.PIPE, text=True)

    if completed.stderr:
        log("[phylo_align] --- MAFFT output ---")
        log(completed.stderr.rstrip())
        log("[phylo_align] --------------------")

    if completed.returncode != 0:
        log(f"[phylo_align] ERROR: MAFFT exited with status {completed.returncode}")
        log_fh.close()
        return completed.returncode

    # Verify the alignment describes exactly the ASVs that went in. MAFFT should never
    # alter a header given the flags above, but this module's whole correctness rests
    # on tips matching the abundance table, and checking after every transformation
    # means a problem is reported where it happened rather than being blamed on the
    # tree backend two rules later.
    aligned_ids = parse_fasta_ids(aln_out)
    if aligned_ids != input_ids:
        missing = sorted(set(input_ids) - set(aligned_ids))[:10]
        added   = sorted(set(aligned_ids) - set(input_ids))[:10]
        log("[phylo_align] ERROR: the alignment's sequence names do not match the input.")
        if missing:
            log(f"[phylo_align]   missing from the alignment (up to 10): {missing}")
        if added:
            log(f"[phylo_align]   unexpected in the alignment (up to 10): {added}")
        if len(aligned_ids) == len(input_ids) and set(aligned_ids) == set(input_ids):
            log("[phylo_align]   (same names, different order — check for a stray "
                "--inputorder override in extra_args.mafft)")
        log_fh.close()
        return 1

    # Every aligned sequence must be the same length; that is what makes it an
    # alignment. Cheap to check and it turns a malformed file into a clear message.
    alignment_lengths = set()
    current_length = 0
    with aln_out.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if current_length:
                    alignment_lengths.add(current_length)
                current_length = 0
            else:
                current_length += len(line.strip())
    if current_length:
        alignment_lengths.add(current_length)

    if len(alignment_lengths) != 1:
        log(f"[phylo_align] ERROR: aligned sequences have differing lengths "
            f"{sorted(alignment_lengths)} — this is not a valid alignment.")
        log_fh.close()
        return 1

    alignment_length = alignment_lengths.pop()
    log(f"[phylo_align] Alignment: {n_seqs} sequences x {alignment_length} columns")

    record = {
        "tool": "mafft",
        "version": mafft_version(log),
        "strategy_requested": requested_strategy,
        "strategy_resolved": strategy_name,
        "linsi_threshold": linsi_threshold,
        "n_sequences": n_seqs,
        "alignment_length": alignment_length,
        "threads": threads,
        "extra_args": extra_args,
        "command": " ".join(command),
    }
    run_json_out.parent.mkdir(parents=True, exist_ok=True)
    run_json_out.write_text(json.dumps(record, indent=2) + "\n")

    log(f"[phylo_align] DONE. {aln_out.name} written, all {n_seqs} ASV names preserved.")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
