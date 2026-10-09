#!/usr/bin/env python3
"""Build MetaFlux's sylph tables (rule sylph_tables, 38_sylph.smk).

WHAT sylph REPORTS
  sylph does not classify reads one by one, the way Kraken2 does. It compares a k-mer
  sketch of each sample with the genomes in its databases and, for every genome it
  detects (at 95 % ANI or more), estimates two shares:
    Sequence_abundance   share of the sample's reads explained by that genome.
                         MetaFlux always runs `sylph profile -u`, so these shares add
                         up to the part of the sample the detected genomes explain,
                         not to 100 %. What is left over is the unknown part.
    Taxonomic_abundance  share of genome copies (roughly, cells) among the DETECTED
                         genomes only. Adds up to 100 % in every sample with a hit.

WHAT THIS SCRIPT WRITES
  03.abundance/sylph_table.tsv     ESTIMATED read pairs per genome and sample:
                                       Sequence_abundance / 100 x read pairs after QC
                                   rounded to whole pairs, plus a last row
                                   'unassigned' = read pairs no detected genome
                                   explains. Each column adds up to the sample's read
                                   pairs after QC (the 'trimmed' column of
                                   read_tracking.txt). One exception, logged: a sample
                                   sylph explains to 100 %, where rounding can push the
                                   genome rows a few pairs past the total; 'unassigned'
                                   is then 0.
  03.abundance/sylph_table.biom    the same table in BIOM (HDF5) format, for phyloseq.
  03.abundance/sylph_table.genome_copies.tsv
                                   Taxonomic_abundance (%) per genome and sample, same
                                   rows, no 'unassigned' row.
  03.abundance/sylph_genomes.tsv   one row per genome: how its name reads, the FASTA
                                   header of its reference genome, and a flag when the
                                   two disagree in a way a reader should know about.
  stats/sylph/sylph_summary.tsv    one row per sample: how much of it the detected
                                   genomes explain, and how many genomes were found.
  stats/sylph/sylph_summary_mqc.tsv
                                   the same rows with readable titles, for MultiQC.

INPUTS (all named by rule sylph_tables)
  input.taxprof     02c.sylph/taxprof/, one {sample}.sylphmpa per sample with at least
                    one detected genome (rule sylph_taxprof, sylph-tax). Each genome row
                    carries the full lineage and the two abundances.
  input.profile     03.abundance/sylph_profile.tsv, the raw `sylph profile -u` output
                    (rule sylph_profile): one row per sample x genome.
  input.per_sample  02c.sylph/sylph_profile.per_sample.txt, sylph's own per-sample log
                    lines (rule sylph_profile).
  input.fastp_jsons the fastp reports (rule trim_adapters, 30_preprocess.smk).
  params.samples    SAMPLES, which fixes the column order.
  params.databases  'file (taxonomy)' of every database, for line 1 of the tables.

WHY SO MANY CHECKS
  sylph and sylph-tax finish with exit code 0 in several situations that would give a
  wrong table without any sign of trouble. Each one is checked below and stops the
  rule with a message that says what to fix:
    - a genome found in none of the lineage tables (NO_TAXONOMY): wrong 'taxonomy'
      for a database in references.sylph_databases;
    - a lineage that is not domain ... species (viral taxonomies);
    - a sample whose sylph-tax file is missing or renamed (sample names rewritten by
      pandas, which sylph-tax uses to read sylph's output);
    - genome counts that disagree between sylph's profile, its log and sylph-tax.
  A sample in which no genome was detected is NOT an error: it gets a column of
  zeros with all its read pairs in 'unassigned'.
"""
import json
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import biom
import numpy
from biom.util import biom_open

sm = snakemake  # noqa: F821  (the object Snakemake hands to every script: rule)

# ── Settings handed over by rule sylph_tables ────────────────────────────────
samples         = list(sm.params.samples)        # SAMPLES, sorted: the column order
databases       = list(sm.params.databases)      # 'file (taxonomy)' per database
taxprof_dir     = Path(sm.input.taxprof)
profile_path    = Path(sm.input.profile)
per_sample_path = Path(sm.input.per_sample)
fastp_jsons     = [Path(p) for p in sm.input.fastp_jsons]

log_fh = open(sm.log[0], "w")


def logmsg(msg: str) -> None:
    """One line to the rule log, flushed at once so the log is complete even if the
    script stops part-way through."""
    log_fh.write(f"[sylph_tables] {msg}\n")
    log_fh.flush()


def stop(msg: str) -> None:
    """Write the reason to the log AND to the screen, then fail the rule."""
    logmsg("ERROR: " + msg)
    sys.stderr.write(f"[MetaFlux] sylph_tables: {msg}\n")
    sys.exit(1)


# ── 1. Read pairs after QC, per sample ───────────────────────────────────────
# fastp counts the two mates separately, so its count is halved to get read PAIRS,
# exactly as 90a_aggregate_read_counts.py does for the 'trimmed' column of
# read_tracking.txt. fastp is the last cleaning step (after PhiX and host removal),
# so these are precisely the read pairs rule sylph_sketch sketched.
# Each JSON is matched to its sample by its exact file name, never by a prefix (a
# prefix would mix up samples such as S1 and S10).
fastp_by_file_name = {}
for json_path in fastp_jsons:
    fastp_by_file_name[json_path.name] = json_path

pairs_after_qc = {}
for sample in samples:
    file_name = f"{sample}_fastp.json"
    if file_name not in fastp_by_file_name:
        stop(f"no fastp report {file_name} among the inputs")
    report = json.loads(fastp_by_file_name[file_name].read_text())
    pairs_after_qc[sample] = report["filtering_result"]["passed_filter_reads"] // 2


# ── 2. The raw sylph profile ─────────────────────────────────────────────────
# Read for three things:
#   - a check that sylph ran with -u: with -u it names column 6 True_cov, without it
#     Eff_cov. Without -u the shares are stretched to 100 % and the 'unassigned' row
#     would be meaningless, so this is checked rather than assumed;
#   - how many genomes sylph reported per sample (compared with sylph-tax below);
#   - each genome's first FASTA header (Contig_name), used only for the name notes in
#     sylph_genomes.tsv. No number in any table comes from this column.
# The line is split into at most as many fields as the header has, so a FASTA header
# that itself contains a tab (possible in a custom database) stays in one piece in the
# last column, Contig_name.

# The two functions below are copied word for word from sylph-tax 1.9.2
# (sylph_tax/sylph_to_taxprof.py). They turn sylph's Genome_file into the accession
# sylph-tax looks up first, which is also the t__ label in its output. Only the genome
# label and the placeholder flag of sylph_genomes.tsv depend on them.
def genome_file_to_gcf_acc(file_name):
    if "ASM" in file_name:
        return file_name.split("/")[-1].split("_ASM")[0]
    return file_name.split("/")[-1].split("_genomic")[0]

def genome_file_to_gcf_acc_STRIP_DOT(file_name):
    return file_name.split("/")[-1].split(".")[0]


rows_in_profile = {}          # sample -> number of genome rows sylph wrote for it
for sample in samples:
    rows_in_profile[sample] = 0
contig_name_of = {}           # genome accession -> first FASTA header of that genome

with profile_path.open() as fh:
    header = fh.readline().rstrip("\n").split("\t")
    if "True_cov" not in header:
        stop(f"{profile_path.name} has no True_cov column, so it was not made with "
             f"`sylph profile -u` and the unassigned row cannot be computed. This should "
             f"not happen inside MetaFlux; delete the file and run again.")
    column = {}
    for position, name in enumerate(header):
        column[name] = position
    for line in fh:
        fields = line.rstrip("\n").split("\t", len(header) - 1)
        sample = fields[column["Sample_file"]]
        if sample not in rows_in_profile:
            stop(f"{profile_path.name} has a sample {sample!r} that is not in this run "
                 f"(samples: {', '.join(samples)})")
        rows_in_profile[sample] += 1
        genome_file = fields[column["Genome_file"]]
        contig_name = fields[column["Contig_name"]]
        # Both forms of the accession point at the same header; the first genome seen
        # under a given accession wins (they are the same genome in every sample).
        contig_name_of.setdefault(genome_file_to_gcf_acc(genome_file), contig_name)
        contig_name_of.setdefault(genome_file_to_gcf_acc_STRIP_DOT(genome_file), contig_name)


# ── 3. The genome rows of every sample, from sylph-tax ───────────────────────
# A .sylphmpa file has one row per rank level (d__, d__|p__, ...) plus one row per
# genome, whose clade name ends in '|t__<accession>'. Only the genome rows are read:
# they carry sylph's own numbers, while the higher rows are sums with floating-point
# noise (99.99990000000001). Columns: clade_name, relative_abundance (=
# Taxonomic_abundance), sequence_abundance (= Sequence_abundance), ANI (=
# Adjusted_ANI), Coverage (= True_cov).
#
# Lineage shape. MetaFlux tables have seven slots, k__ ... s__ (domain, phylum,
# class, order, family, genus, species). GTDB writes the domain as d__ and the fungal
# RefSeq taxonomy as k__; mixed in one table, phyloseq and mia would read them as two
# different ranks. So the first slot is always written k__ (it holds the domain, as
# in otu_table.tsv). Names are not touched: GTDB's _A/_B suffixes and placeholder names
# stay exactly as GTDB writes them. Two quirks are accepted: a rank the taxonomy leaves
# empty stays a bare prefix (c__), and GlobDB r226's 'p_Unclassified' (one underscore)
# becomes p__Unclassified. Anything else (an UNKNOWN level, eight levels as in the
# IMG/VR viral taxonomy, another prefix) stops the rule rather than producing a table
# phyloseq would read wrongly.
SLOT_LETTERS = ["k", "p", "c", "o", "f", "g", "s"]


def seven_slots(ranks: list, where: str) -> list:
    """['d__Bacteria', 'p__Bacillota', ..., 's__Veillonella rogosae'] ->
    ['k__Bacteria', 'p__Bacillota', ..., 's__Veillonella rogosae']."""
    lineage_text = "|".join(ranks)
    if len(ranks) != 7:
        stop(f"{where}: the lineage has {len(ranks)} ranks where MetaFlux expects 7 "
             f"(domain ... species): {lineage_text}. Only GTDB-style and fungal RefSeq "
             f"taxonomies are supported.")
    slots = []
    for position in range(7):
        rank = ranks[position]
        expected_letter = SLOT_LETTERS[position]
        # one letter, then one or two underscores, then the name (possibly empty)
        match = re.fullmatch(r"([a-z])_{1,2}(.*)", rank)
        if position == 0:
            allowed_letters = ("d", "k")       # GTDB d__ or RefSeq k__: both the domain
        else:
            allowed_letters = (expected_letter,)
        if match is None or match.group(1) not in allowed_letters:
            stop(f"{where}: expected a {expected_letter}__ rank in position {position + 1}, "
                 f"found {rank!r}, in the lineage {lineage_text}")
        slots.append(f"{expected_letter}__{match.group(2)}")
    return slots


def four_decimals(text: str) -> Decimal:
    """sylph prints every share with 4 decimals; sylph-tax writes the same number as a
    Python float (8.2320 -> '8.232'). Decimal keeps it exact, and the rounding to 4
    decimals is only a safety net (none of the values checked needed it)."""
    return Decimal(text).quantize(Decimal("0.0001"))


lineage_of = {}           # genome id -> its seven slots
seq_abund = {}            # (genome id, sample) -> Sequence_abundance, %
tax_abund = {}            # (genome id, sample) -> Taxonomic_abundance, %
ani_of = {}               # (genome id, sample) -> Adjusted_ANI, as sylph-tax wrote it
coverage_of = {}          # (genome id, sample) -> True_cov, as sylph-tax wrote it
rows_in_taxprof = {}      # sample -> genome rows sylph-tax wrote for it
for sample in samples:
    rows_in_taxprof[sample] = 0
samples_with_file = set()
no_lineage = []           # 'sample: genome' for genomes sylph-tax could not name

for mpa_path in sorted(taxprof_dir.glob("*.sylphmpa")):
    sample = mpa_path.name[: -len(".sylphmpa")]
    # sylph-tax names each file after sylph's Sample_file column, which is the MetaFlux
    # sample name (rule sylph_sketch passes -S {sample}). A file named after anything
    # else means the name was changed on the way. 00_common.smk refuses the sample
    # names known to do that, so this is a safety net.
    if sample not in rows_in_taxprof:
        stop(f"sylph-tax wrote {mpa_path.name}, which matches no sample of this run "
             f"(names rewritten on the way?). Samples: {', '.join(samples)}")
    samples_with_file.add(sample)
    with mpa_path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("clade_name"):
                continue                         # the two header lines
            fields = line.rstrip("\n").split("\t")
            clade = fields[0]
            if "|t__" not in clade:
                continue                         # a rank-sum row, not a genome
            lineage_part, genome_id = clade.split("|t__", 1)
            # sylph-tax writes 'NO_TAXONOMY|t__<genome file>:<contig>' for a genome it
            # found in none of the lineage tables, and carries on with exit code 0.
            # Collect them all, so the message below can list them.
            if lineage_part == "NO_TAXONOMY":
                no_lineage.append(f"{sample}: {genome_id}")
                continue
            slots = seven_slots(lineage_part.split("|"), f"sample {sample}, genome {genome_id}")
            if (genome_id, sample) in seq_abund:
                stop(f"sample {sample}: genome {genome_id} is reported twice; is the same "
                     f"genome in two of references.sylph_databases?")
            if genome_id in lineage_of and lineage_of[genome_id] != slots:
                stop(f"genome {genome_id} has two different lineages "
                     f"({';'.join(lineage_of[genome_id])} and {';'.join(slots)}); is the "
                     f"same genome in two databases with different taxonomies?")
            lineage_of[genome_id] = slots
            tax_abund[(genome_id, sample)] = four_decimals(fields[1])
            seq_abund[(genome_id, sample)] = four_decimals(fields[2])
            ani_of[(genome_id, sample)] = fields[3]
            coverage_of[(genome_id, sample)] = fields[4]
            rows_in_taxprof[sample] += 1

if no_lineage:
    first_ones = "; ".join(no_lineage[:5])
    stop(f"{len(no_lineage)} detected genome(s) have no lineage in the taxonomies listed "
         f"under references.sylph_databases (first ones: {first_ones}). Each database "
         f"needs its own matching taxonomy, e.g. GTDB_r232 for gtdb-r232-*.syldb and "
         f"FungiRefSeq-latest for fungi-refseq-2025-10-11-*.syldb. Add or correct the "
         f"'taxonomy:' entries and run again (only sylph_taxprof and sylph_tables run "
         f"again).")


# ── 4. Cross-check: sylph-tax wrote every genome sylph reported ──────────────
# Per sample, the genome rows in sylph-tax's file must equal the rows in sylph's
# profile (0 for a sample with no hit, which has no file at all). A sample with
# genomes in the profile but no file is the pandas renaming problem: sylph-tax read
# the sample name as a number or a missing value and wrote its file under another
# name (or merged it with another sample).
for sample in samples:
    in_profile = rows_in_profile[sample]
    in_taxprof = rows_in_taxprof[sample]
    if in_profile > 0 and sample not in samples_with_file:
        stop(f"sylph-tax wrote no {sample}.sylphmpa although sylph detected {in_profile} "
             f"genome(s) in {sample}; names rewritten by pandas? sylph-tax reads sample "
             f"names with pandas, which rewrites names that look like numbers (001 "
             f"becomes 1) or like missing values (NA). Rename such samples.")
    if in_profile != in_taxprof:
        stop(f"sample {sample}: sylph's profile has {in_profile} genome row(s) but "
             f"sylph-tax wrote {in_taxprof}; the two outputs do not match.")


# ── 5. Per-sample facts from sylph's own log lines ───────────────────────────
# rule sylph_profile keeps three kinds of line from sylph's log, per sample:
#   '<sketch> has 16 genomes passing profiling threshold.'
#   '<sketch> has 39.58% of reads detected in database by profile'
#   '<R1 file> short-read sample has high diversity compared to sequencing depth ...
#    Using 99.5% as read accuracy estimate instead of automatic detection ...'
# The first two name the sketch file (02c.sylph/sketches/<sample>.paired.sylsp); the
# sample is its file name without '.paired.sylsp'. The third, printed only for a
# sample too shallow for sylph to estimate its read accuracy, names the trimmed R1
# file the sketch was made from (<sample>_trim_R1.fastq.gz). That fallback moves the
# unknown (unassigned) estimate by several points, so the summary says which samples
# it applied to. Samples are always matched by their exact name, never by a prefix.
GENOMES_LINE  = re.compile(r"\] (.+?) has (\d+) genomes passing profiling threshold")
DETECTED_LINE = re.compile(r"\] (.+?) has ([0-9.]+)% of reads detected in database by profile")
FALLBACK_LINE = re.compile(r"\] (.+?) short-read sample has high diversity")


def sample_of_sketch(path_text: str) -> str:
    """'/out/02c.sylph/sketches/SRR29240428.paired.sylsp' -> 'SRR29240428'."""
    file_name = Path(path_text).name
    if file_name.endswith(".paired.sylsp"):
        file_name = file_name[: -len(".paired.sylsp")]
    return file_name


genomes_passing = {}      # sample -> genome count sylph printed
detected_pct = {}         # sample -> '% of reads detected' sylph printed
low_depth_samples = set() # samples that got the fixed 99.5 % read-accuracy fallback

trimmed_r1_to_sample = {}
for sample in samples:
    trimmed_r1_to_sample[f"{sample}_trim_R1.fastq.gz"] = sample

with per_sample_path.open() as fh:
    for line in fh:
        match = GENOMES_LINE.search(line)
        if match:
            sample = sample_of_sketch(match.group(1))
            if sample not in rows_in_profile:
                stop(f"sylph's log names a sketch {match.group(1)} that matches no sample")
            genomes_passing[sample] = int(match.group(2))
            continue
        match = DETECTED_LINE.search(line)
        if match:
            sample = sample_of_sketch(match.group(1))
            if sample not in rows_in_profile:
                stop(f"sylph's log names a sketch {match.group(1)} that matches no sample")
            if sample in detected_pct:
                stop(f"sylph's log has two results for sample {sample}")
            detected_pct[sample] = Decimal(match.group(2))
            continue
        match = FALLBACK_LINE.search(line)
        if match:
            r1_name = Path(match.group(1)).name
            if r1_name in trimmed_r1_to_sample:
                low_depth_samples.add(trimmed_r1_to_sample[r1_name])
            else:
                logmsg(f"warning: sylph used the low-depth read-accuracy fallback for "
                       f"{r1_name}, which matches no trimmed read file of this run")

# Every sample must have its result line, and sylph's genome count must equal the
# rows it wrote. A missing result line means sylph never profiled that sample (rule
# sylph_profile already checks the number of lines); without this check the sample
# would silently look like 'nothing detected'.
for sample in samples:
    if sample not in detected_pct:
        stop(f"sylph's log ({per_sample_path.name}) has no result for sample {sample}, so "
             f"it was not profiled. Delete 03.abundance/sylph_profile.tsv and run again.")
    if genomes_passing.get(sample) != rows_in_profile[sample]:
        stop(f"sample {sample}: sylph's log reports {genomes_passing.get(sample)} genome(s) "
             f"but its profile has {rows_in_profile[sample]} row(s)")


# ── 6. Row order ─────────────────────────────────────────────────────────────
# Sorted by the full lineage, then by genome id: related genomes sit next to each
# other, and the file is the same on every run whatever order sylph wrote.
def lineage_then_id(genome_id: str):
    return (";".join(lineage_of[genome_id]), genome_id)


genome_ids = sorted(lineage_of, key=lineage_then_id)


# ── 7. Name notes, one per genome (sylph_genomes.tsv) ────────────────────────
# GTDB names are shipped exactly as GTDB writes them. Two things help a reader:
#
# name_status, from the lineage alone:
#   'placeholder genus'    the genus name contains a digit (GTDB's code for a genus
#                          without a Latin name, e.g. G047199095);
#   'no species name'      the species slot has fewer than two words;
#   'placeholder species'  the species epithet contains a digit or a capital letter
#                          (sp018646675). GTDB's split suffixes (_A, _E) are removed
#                          first: 'Escherichia coli_E' is a real, named species;
#   'named'                everything else.
#
# flag, for one case worth knowing about: GTDB puts the genome in a placeholder genus
# while the genome's own FASTA header (sylph's Contig_name, the name the submitter
# gave it in NCBI) is a Latin binomial of ANOTHER genus. Measured case: GTDB r232 puts
# two 2025 Escherichia coli assemblies in their own placeholder genera; on the Zymo
# mock they take the E. coli reads, and no 'Escherichia coli' row exists. The flag
# never changes a number or a name in the tables; it only explains a row. Across
# 9,570 benchmark hits it fired only on those E. coli genomes.
NOT_A_SPECIES_WORD = {"sp.", "bacterium", "archaeon", "genomosp.", "cf.", "aff.",
                      "uncultured", "endosymbiont", "symbiont", "group", "cyanobacterium",
                      "proteobacterium", "actinobacterium", "clade"}


def name_status(slots: list) -> str:
    genus = slots[5][3:]
    species = slots[6][3:]
    if re.search(r"\d", genus):
        return "placeholder genus"
    words = species.split()
    if len(words) < 2:
        return "no species name"
    epithet = re.sub(r"_[A-Z]+$", "", words[1])      # GTDB split suffix: coli_E -> coli
    if re.search(r"[\dA-Z]", epithet):
        return "placeholder species"
    return "named"


def latin_binomial(label: str):
    """'Escherichia coli strain ECGI0523 chromosome' -> ('Escherichia', 'coli').
    None when the label does not start with a plain Latin genus and species: square
    brackets are ignored, a leading 'Candidatus' is skipped, and a first word that is
    a family, order or phylum name, or a second word such as 'sp.' or 'bacterium',
    does not count."""
    words = label.replace("[", "").replace("]", "").split()
    if words and words[0] == "Candidatus":
        words = words[1:]
    if len(words) < 2:
        return None
    genus, epithet = words[0], words[1]
    if not re.fullmatch(r"[A-Z][a-z]+", genus):
        return None
    if re.search(r"(aceae|ales|ota|bacteria)$", genus):
        return None
    if not re.fullmatch(r"[a-z][a-z-]+", epithet) or epithet in NOT_A_SPECIES_WORD:
        return None
    return genus, epithet


status_of = {}            # genome id -> name_status
label_of = {}             # genome id -> FASTA header without its first word (the sequence id)
flag_of = {}              # genome id -> flag text, or ''
for genome_id in genome_ids:
    slots = lineage_of[genome_id]
    status_of[genome_id] = name_status(slots)

    contig_name = contig_name_of.get(genome_id, "")
    label = " ".join(contig_name.split()[1:])        # drop the sequence id (CP178306.1)
    label_of[genome_id] = label

    flag = ""
    gtdb_genus = slots[5][3:]
    if re.search(r"\d", gtdb_genus):
        binomial = latin_binomial(label)
        if binomial is not None:
            gtdb_genus_without_suffix = re.sub(r"_[A-Z]+$", "", gtdb_genus)
            if binomial[0] != gtdb_genus_without_suffix:
                flag = f"GTDB placeholder genus; genome labelled {binomial[0]} {binomial[1]}"
    flag_of[genome_id] = flag
    if flag:
        logmsg(f"flag: {genome_id} ({slots[6][3:]}): {flag}")


# ── 8. Estimated read pairs, and the unassigned remainder ───────────────────
# estimated pairs = Sequence_abundance / 100 x read pairs after QC, rounded to a
# whole pair with a half rounded up (as rule bracken rounds). Decimal does this
# exactly on sylph's 4-decimal percentages, so a value sitting exactly on .5 is never
# tipped either way by binary floating-point error.
# 'unassigned' takes whatever is left, so each column adds up to the read pairs after
# QC. In a sample sylph explains to 100 % (mock communities, tiny samples), the
# rounded genome rows can add up to a few pairs MORE than the total; 'unassigned' is
# then set to 0 and the overshoot logged, rather than writing a negative count.
def estimated_pairs(percent: Decimal, pairs: int) -> int:
    exact = percent / Decimal(100) * Decimal(pairs)
    return int(exact.quantize(Decimal(1), rounding=ROUND_HALF_UP))


est_pairs = {}            # (genome id, sample) -> estimated read pairs
explained_pairs = {}      # sample -> sum of the genome rows
unassigned_pairs = {}     # sample -> read pairs no detected genome explains
explained_pct = {}        # sample -> sum of Sequence_abundance, %
for sample in samples:
    pairs = pairs_after_qc[sample]
    assigned = 0
    percent_sum = Decimal(0)
    for genome_id in genome_ids:
        percent = seq_abund.get((genome_id, sample), Decimal(0))
        n_pairs = estimated_pairs(percent, pairs)
        est_pairs[(genome_id, sample)] = n_pairs
        assigned += n_pairs
        percent_sum += percent
    left_over = pairs - assigned
    if left_over < 0:
        logmsg(f"{sample}: genome estimates exceed the read pairs by {-left_over} "
               f"(rounding); unassigned set to 0")
        left_over = 0
    explained_pairs[sample] = assigned
    unassigned_pairs[sample] = left_over
    explained_pct[sample] = percent_sum


# ── 9. sylph_table.tsv and sylph_table.biom, from the same numbers ───────────
# One matrix feeds both files, so the TSV and the BIOM can never disagree. The layout
# is otu_table.tsv's: line 1 a comment, line 2 '#OTU ID', the samples, 'taxonomy';
# then one row per genome and the 'unassigned' row last. Values are whole numbers.
UNASSIGNED_SLOTS = ["k__Unassigned", "p__", "c__", "o__", "f__", "g__", "s__"]
database_text = ", ".join(databases)

row_ids = []
row_slots = []
matrix = []
for genome_id in genome_ids:
    row_ids.append(genome_id)
    row_slots.append(lineage_of[genome_id])
    values = []
    for sample in samples:
        values.append(est_pairs[(genome_id, sample)])
    matrix.append(values)
row_ids.append("unassigned")
row_slots.append(UNASSIGNED_SLOTS)
unassigned_values = []
for sample in samples:
    unassigned_values.append(unassigned_pairs[sample])
matrix.append(unassigned_values)

with open(sm.output.table, "w") as fh:
    fh.write("# MetaFlux sylph table: estimated reads (coverage-derived), counted in read "
             "pairs; last row 'unassigned' = read pairs no detected genome explains "
             f"(sylph -u). Databases: {database_text}\n")
    fh.write("#OTU ID\t" + "\t".join(samples) + "\ttaxonomy\n")
    for row_id, slots, values in zip(row_ids, row_slots, matrix):
        value_texts = [str(v) for v in values]
        fh.write("\t".join([row_id] + value_texts + [";".join(slots)]) + "\n")

# BIOM twin, for phyloseq (import_biom) and other tools. biom-format stores every
# matrix as floating point, so the whole numbers come back as 896380.0. The second
# argument of to_hdf5 is a free-text "generated by" tag stored in the file. The file
# also records its creation time, so two runs give different bytes: compare BIOM files
# by content, as for otu_table.biom.
observation_metadata = []
for slots in row_slots:
    observation_metadata.append({"taxonomy": slots})
biom_table = biom.Table(numpy.array(matrix, dtype=float), row_ids, samples,
                        observation_metadata=observation_metadata, type="OTU table")
with biom_open(str(sm.output.biom), "w") as fh:
    biom_table.to_hdf5(fh, "MetaFlux sylph_tables (estimated read pairs, coverage-derived)")


# ── 10. sylph_table.genome_copies.tsv ────────────────────────────────────────
# sylph's Taxonomic_abundance, 4 decimals, same genome rows in the same order, no
# 'unassigned' row (these shares are among detected genomes only). A genome absent
# from a sample is 0.0000. Each column with a hit adds up to 100 within rounding:
# every value is rounded to 4 decimals, so the sum can be off by up to about
# 0.00005 x the number of genomes. A sample with no hit is all 0.0000.
with open(sm.output.copies, "w") as fh:
    fh.write("# MetaFlux sylph table: genome copies (sylph Taxonomic_abundance, % of "
             f"DETECTED genomes; no unassigned row). Databases: {database_text}\n")
    fh.write("#OTU ID\t" + "\t".join(samples) + "\ttaxonomy\n")
    for genome_id in genome_ids:
        value_texts = []
        for sample in samples:
            share = tax_abund.get((genome_id, sample), Decimal(0))
            value_texts.append(f"{share:.4f}")
        fh.write("\t".join([genome_id] + value_texts + [";".join(lineage_of[genome_id])]) + "\n")


# ── 11. sylph_genomes.tsv ────────────────────────────────────────────────────
# One row per genome row of sylph_table.tsv, in the same order: the GTDB species
# name, name_status and flag (section 7), the FASTA header of the reference genome
# (empty if it could not be linked), and in how many samples it was detected with
# its highest coverage (True_cov) and identity (Adjusted_ANI) over those samples.
def highest(texts: list) -> str:
    """The largest of several numbers written as text, returned as written."""
    best = texts[0]
    for text in texts[1:]:
        if Decimal(text) > Decimal(best):
            best = text
    return best


with open(sm.output.genomes, "w") as fh:
    fh.write("\t".join(["genome_id", "species", "name_status", "genome_label", "flag",
                        "samples_detected", "max_true_cov", "max_adjusted_ani"]) + "\n")
    for genome_id in genome_ids:
        coverages = []
        anis = []
        for sample in samples:
            if (genome_id, sample) in tax_abund:
                coverages.append(coverage_of[(genome_id, sample)])
                anis.append(ani_of[(genome_id, sample)])
        fh.write("\t".join([
            genome_id,
            lineage_of[genome_id][6][3:],
            status_of[genome_id],
            label_of[genome_id],
            flag_of[genome_id],
            str(len(coverages)),
            highest(coverages),
            highest(anis),
        ]) + "\n")


# ── 12. Per-sample summary (sylph_summary.tsv and its MultiQC copy) ──────────
# One row per sample, in SAMPLES order:
#   read_pairs_after_qc        fastp's count, the column total of sylph_table.tsv
#   est_pairs_explained        sum of the genome rows
#   est_pairs_unassigned       the 'unassigned' row
#   pct_explained              sum of Sequence_abundance, 2 decimals
#   explained_capped           'yes' when that rounds to 100.00: sylph explained the
#                              whole sample, so its unknown share cannot be estimated
#   read_accuracy              'automatic', or 'fixed 99.5% (low depth)' (section 5)
#   genomes_detected           genome rows sylph reported
#   pct_pairs_unnamed_species  estimated pairs on genomes whose name_status is not
#                              'named', as % of est_pairs_explained (NA if that is 0)
#   flagged_genomes            flagged genomes (section 7) detected in this sample
summary_rows = []
for sample in samples:
    pct_text = str(explained_pct[sample].quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
    if Decimal(pct_text) >= Decimal("100.00"):
        capped = "yes"
    else:
        capped = "no"
    if sample in low_depth_samples:
        read_accuracy = "fixed 99.5% (low depth)"
    else:
        read_accuracy = "automatic"

    unnamed_pairs = 0
    flagged = 0
    for genome_id in genome_ids:
        if (genome_id, sample) not in tax_abund:
            continue                                  # not detected in this sample
        if status_of[genome_id] != "named":
            unnamed_pairs += est_pairs[(genome_id, sample)]
        if flag_of[genome_id]:
            flagged += 1
    if explained_pairs[sample] > 0:
        unnamed_share = Decimal(100) * Decimal(unnamed_pairs) / Decimal(explained_pairs[sample])
        unnamed_text = str(unnamed_share.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
    else:
        unnamed_text = "NA"

    summary_rows.append([
        sample,
        str(pairs_after_qc[sample]),
        str(explained_pairs[sample]),
        str(unassigned_pairs[sample]),
        pct_text,
        capped,
        read_accuracy,
        str(rows_in_profile[sample]),
        unnamed_text,
        str(flagged),
    ])
    logmsg(f"{sample}: {pairs_after_qc[sample]} read pairs after QC, "
           f"{rows_in_profile[sample]} genome(s), {pct_text} % explained "
           f"(sylph's log: {detected_pct[sample]} %), {explained_pairs[sample]} pairs "
           f"assigned, {unassigned_pairs[sample]} unassigned, read accuracy {read_accuracy}")

with open(sm.output.summary, "w") as fh:
    fh.write("\t".join(["sample", "read_pairs_after_qc", "est_pairs_explained",
                        "est_pairs_unassigned", "pct_explained", "explained_capped",
                        "read_accuracy", "genomes_detected", "pct_pairs_unnamed_species",
                        "flagged_genomes"]) + "\n")
    for row in summary_rows:
        fh.write("\t".join(row) + "\n")

# The MultiQC copy. MultiQC picks it up because its name ends in _mqc.tsv and it sits
# in stats/sylph/, a folder rule multiqc scans only when sylph runs (90_report.smk). It
# reads the '# ' lines as the section's settings and shows the rest as one table in its
# own section. Rows keep the plain sample names: MultiQC does not add the folder
# prefix to the rows of a custom-content table (checked with MultiQC 1.35).
# The wording avoids any phrase that would suggest reads were classified one by one.
with open(sm.output.mqc, "w") as fh:
    fh.write("# id: 'metaflux_sylph'\n")
    fh.write("# section_name: 'sylph (coverage-based profile)'\n")
    fh.write("# description: 'Read pairs here are ESTIMATED from genome coverage (sylph "
             "profile -u); sylph labels no individual read. Explained = the share of the "
             "read pairs after QC that the detected genomes account for; the rest is the "
             "\"unassigned\" row of 03.abundance/sylph_table.tsv. Capped = sylph reported "
             "100 %, so the unknown share of that sample cannot be estimated. Read accuracy "
             "\"fixed 99.5% (low depth)\" = too little depth for sylph to estimate it, so "
             "the unassigned share is rougher. Flagged genomes: see "
             "03.abundance/sylph_genomes.tsv.'\n")
    fh.write("# plot_type: 'table'\n")
    fh.write("# pconfig:\n")
    fh.write("#     id: 'metaflux_sylph_table'\n")
    fh.write("#     namespace: 'sylph'\n")
    fh.write("\t".join(["Sample", "Read pairs after QC",
                        "Explained, est. pairs (coverage-derived)",
                        "Unassigned, est. pairs", "Explained (%)", "Capped at 100 %",
                        "Read accuracy", "Genomes",
                        "Unnamed species (% of explained pairs)",
                        "Flagged genomes"]) + "\n")
    for row in summary_rows:
        fh.write("\t".join(row) + "\n")

logmsg(f"wrote {len(genome_ids)} genome(s) x {len(samples)} sample(s)")
log_fh.close()
