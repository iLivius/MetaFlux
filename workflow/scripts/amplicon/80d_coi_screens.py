#!/usr/bin/env python3
"""COI screens: the reading-frame check and the off-target guard, per ASV.

Why COI needs them
    COI is a protein-coding gene, and two kinds of sequence reach a COI ASV table
    that are not the mitochondrial COI of a eukaryote:
      1. Pseudogenes and NUMTs (nuclear copies of mitochondrial genes) and the
         few sequencing artefacts DADA2 lets through. A functional COI amplicon
         has no stop codon in its reading frame, and its length differs from the
         expected one only by whole codons; a NUMT or an artefact with an
         insertion or deletion breaks one or the other.
      2. Bacterial COX1. Many bacteria carry COX1 homologues that COI primers
         amplify, above all Wolbachia and other Rickettsiales in insect samples.
         MIDORI2 holds eukaryotes only, so such an ASV is assigned to a eukaryote.
         The guard compares each ASV with a set of bacterial and archaeal COX1
         genes from GTDB (identity B) and with the cleaned MIDORI2 (identity E):
         an ASV closer to the bacterial set (B > E) is bacterial COX1. Why a
         competitive comparison and not an identity floor is on the manual page
         "Bacterial COX1 reference set (COI)".
    This script measures both for every ASV and writes one row per ASV. It removes
    nothing itself: rule assign_taxonomy (80c_parse_sintax.py) reads the table and
    applies the decisions inside its filter chain, so asv_table_unfiltered.txt
    keeps every ASV with the step that removed it.

The reading-frame check (amplicon.coi.frame_check)
    Every ASV starts right after the forward primer (cutadapt removed it), so all
    ASVs of a run share the same codon phase. The run's reading frame is the one
    of the three frames that is free of TAA and TAG (the stop codons of every
    animal mitochondrial code) in the largest number of reads. The reference
    length is the most frequent ASV length, by reads, among ASVs open in that
    frame. An ASV passes when
      - it has no TAA or TAG in the run's reading frame, and
      - its length differs from the reference length by a multiple of 3, by at
        most amplicon.coi.length_tolerance bases (default 9, three codons).
    TGA (tryptophan in animal and fungal mitochondria) and AGA/AGG (stops in the
    vertebrate code only) are not used, so the check works the same for
    invertebrates and vertebrates. Measured on MIDORI2 GB273 inserts of BF3/BR2,
    Leray-XT and fwhF2/fwhR2n, the frame starting at the second base of the
    insert is free of TAA/TAG in 99.3-100 % of genuine records (per group).

The off-target guard (amplicon.coi.offtarget_guard)
    VSEARCH --usearch_global of every ASV against the cleaned MIDORI2 (E) and
    against the bacterial COX1 set, whole genes and fragments (B), with
        --id 0.5 --strand plus --maxaccepts 256 --maxrejects 2048
    (eight times VSEARCH's default search depth: with distant sequences the best
    alignment is often not among the first targets VSEARCH tries). The identity
    is recounted from the aligned rows over A/C/G/T columns: identical A/C/G/T
    columns / (A/C/G/T columns + internal gap columns). VSEARCH's own identity
    counts N and other ambiguity codes as identical, and about 9 % of MIDORI2
    records contain them. A hit counts only when those columns cover at least
    80 % of the ASV. B > E: the ASV is bacterial COX1.
    The E search also gives what the species flags need (written here, decided in
    80c): the MIDORI2 records tied at the best identity, and their species.

Input  (rule coi_screens in 80_taxonomy.smk)
    seqs       05.dada2/seqs_lenfilt.fasta, the final ASV set
    seqtab     05.dada2/seqtab_lenfilt_head_names.txt, its read counts
    midori2    the cleaned MIDORI2 SINTAX file (prepare_midori2_coi_sintax)
    bacterial  the bacterial COX1 whole genes; fragments: the fragments
Output
    table      06.taxonomy/coi_screens.tsv, one row per ASV (columns below), read
               by rule assign_taxonomy (80c) only
    log        the run's frame and reference length, and counts
"""
import gzip
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sm = snakemake  # noqa: F821 (injected by Snakemake)

seqs_path = Path(sm.input.seqs)
seqtab_path = Path(sm.input.seqtab)
midori2_path = Path(sm.input.midori2)
bacterial_paths = [Path(sm.input.bacterial), Path(sm.input.fragments)]
out_path = Path(sm.output.table)
log_path = Path(sm.log[0])
threads = int(sm.threads)
frame_check = bool(sm.params.frame_check)
offtarget_guard = bool(sm.params.offtarget_guard)
length_tolerance = int(sm.params.length_tolerance)
min_query_coverage = float(sm.params.min_query_coverage)

VSEARCH_SETTINGS = ["--id", "0.5", "--strand", "plus", "--maxaccepts", "256", "--maxrejects", "2048"]
STOPS = ("TAA", "TAG")
ACGT = frozenset("ACGT")
COLUMNS = ["asv_id", "reads", "length", "open_frames", "frame_ok", "length_ok",
           "B", "best_bacterial_gene", "best_bacterial_taxonomy",
           "E", "best_midori2_record", "best_midori2_species", "verdict"]

log_path.parent.mkdir(parents=True, exist_ok=True)
log = log_path.open("w")


def logmsg(message):
    log.write(f"[coi_screens] {message}\n")
    log.flush()


def read_fasta(path):
    """{ASV id: sequence} from the ASV FASTA (one id per header, first word)."""
    sequences = {}
    current = None
    chunks = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(chunks).upper()
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        sequences[current] = "".join(chunks).upper()
    return sequences


def read_totals(path):
    """{ASV id: reads over all samples} from the seqtab (samples x ASVs)."""
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        asv_ids = [cell.strip('"') for cell in header[1:]]
        totals = Counter()
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            for asv_id, count in zip(asv_ids, cells[1:]):
                totals[asv_id] += int(count)
    return totals


def open_frames(sequence):
    """The frames (0, 1, 2) without TAA or TAG, counted from the ASV's first base."""
    frames = []
    for frame in range(3):
        has_stop = False
        for position in range(frame, len(sequence) - 2, 3):
            if sequence[position:position + 3] in STOPS:
                has_stop = True
                break
        if not has_stop:
            frames.append(frame)
    return frames


def acgt_identity(query_row, target_row):
    """Identity (%) over A/C/G/T columns and internal gaps, and the columns counted."""
    same = 0
    counted = 0
    for q, t in zip(query_row.upper(), target_row.upper()):
        if q == "-" or t == "-":
            counted += 1
        elif q in ACGT and t in ACGT:
            counted += 1
            if q == t:
                same += 1
    if counted == 0:
        return 0.0, 0
    return 100.0 * same / counted, counted


def search(database, label):
    """VSEARCH of all ASVs against one database; {ASV: [(identity, target), ...]}."""
    hits_path = out_path.parent / f"_coi_screens_{label}.tsv"
    command = ["vsearch", "--usearch_global", str(seqs_path), "--db", str(database),
               *VSEARCH_SETTINGS, "--threads", str(threads), "--userout", str(hits_path),
               "--userfields", "query+target+ql+qrow+trow"]
    logmsg("running: " + " ".join(command))
    result = subprocess.run(command, stdout=log, stderr=log)
    if result.returncode != 0:
        logmsg(f"ERROR: vsearch exited with code {result.returncode}")
        log.close()
        sys.exit(result.returncode)
    hits = defaultdict(list)
    with hits_path.open() as handle:
        for line in handle:
            query, target, query_length, query_row, target_row = line.rstrip("\n").split("\t")
            identity, counted = acgt_identity(query_row, target_row)
            if counted >= min_query_coverage * int(query_length):
                hits[query.split()[0]].append((identity, target))
    hits_path.unlink()
    return hits


def species_of(label):
    """The species name in a SINTAX label (…;tax=…,s:Genus_species;), or ''."""
    tax = label.split(";tax=", 1)[1] if ";tax=" in label else ""
    for field in tax.rstrip(";").split(","):
        if field.startswith("s:"):
            return field[2:]
    return ""


sequences = read_fasta(seqs_path)
totals = read_totals(seqtab_path)
asv_ids = list(sequences)
logmsg(f"{len(asv_ids)} ASVs, {sum(totals[a] for a in asv_ids)} reads")

# ── Reading frame and reference length of the run ──
frames_of = {asv_id: open_frames(sequences[asv_id]) for asv_id in asv_ids}
reads_open = Counter()
for asv_id in asv_ids:
    for frame in frames_of[asv_id]:
        reads_open[frame] += totals[asv_id]
run_frame = max(range(3), key=lambda frame: (reads_open[frame], -frame))
length_reads = Counter()
for asv_id in asv_ids:
    if run_frame in frames_of[asv_id]:
        length_reads[len(sequences[asv_id])] += totals[asv_id]
reference_length = max(length_reads, key=lambda length: (length_reads[length], -length)) if length_reads else 0
logmsg(f"reads open per frame: {dict(reads_open)}; run frame: {run_frame} "
       f"(0 = first base); reference length: {reference_length} bp")

# ── B and E (always measured; the switch only decides whether B > E moves an ASV) ──
# E also serves the species flags and the far_from_midori2 flag of 80c, which do
# not depend on the guard being on.
best_b = {}
best_e = {}
bacterial_hits = defaultdict(list)
for number, database in enumerate(bacterial_paths):
    for asv_id, hits in search(database, f"bacterial{number}").items():
        bacterial_hits[asv_id].extend(hits)
midori2_hits = search(midori2_path, "midori2")
for asv_id in asv_ids:
    if bacterial_hits.get(asv_id):
        best_b[asv_id] = max(bacterial_hits[asv_id])
    if midori2_hits.get(asv_id):
        top_identity = max(identity for identity, _ in midori2_hits[asv_id])
        tied_targets = [target for identity, target in midori2_hits[asv_id] if identity == top_identity]
        tied_species = sorted({species_of(target) or "no_species" for target in tied_targets})
        first_record = min(target.split(";")[0] for target in tied_targets)
        best_e[asv_id] = (top_identity, first_record, tied_species)

# ── One row per ASV ──
verdicts = Counter()
with out_path.open("w") as out:
    out.write("\t".join(COLUMNS) + "\n")
    for asv_id in asv_ids:
        length = len(sequences[asv_id])
        frame_ok = run_frame in frames_of[asv_id]
        offset = length - reference_length
        length_ok = offset % 3 == 0 and abs(offset) <= length_tolerance
        b_identity, b_gene = best_b.get(asv_id, (0.0, ""))
        e_identity, e_record, e_species = best_e.get(asv_id, (0.0, "", []))
        if frame_check and not (frame_ok and length_ok):
            verdict = "frame"
        elif offtarget_guard and b_identity > e_identity:
            verdict = "bacterial"
        else:
            verdict = "pass"
        verdicts[verdict] += 1
        out.write("\t".join([
            asv_id, str(totals[asv_id]), str(length),
            ",".join(str(f) for f in frames_of[asv_id]) or "none",
            str(int(frame_ok)), str(int(length_ok)),
            f"{b_identity:.1f}",
            b_gene.split(";")[0] if b_gene else "NA",
            b_gene.split(";tax=", 1)[1].rstrip(";") if ";tax=" in b_gene else "NA",
            f"{e_identity:.1f}", e_record or "NA",
            "|".join(e_species) if e_species else "NA", verdict,
        ]) + "\n")
logmsg(f"verdicts: {dict(verdicts)} (frame check {'on' if frame_check else 'off'}, "
       f"off-target guard {'on' if offtarget_guard else 'off'})")
log.close()
