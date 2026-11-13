#!/usr/bin/env python3
"""Decide whether primer trimming removes the last base of every read of one FASTQ file.

Rule trim_primers (workflow/rules/amplicon/30_preprocess.smk) runs this program from its
shell command, once for the R1 file and once for the R2 file of each sample, before the
first Cutadapt call. It is a small command-line program, not a Snakemake script: it runs
in the middle of the rule's shell command, in the Python of the Cutadapt environment, and
hands its answer back to that command, which passes it on to Cutadapt.

Why the last base
-----------------
An Illumina run sequences one cycle more than the read length it is set up for: 2 x 151,
2 x 251 or 2 x 301 cycles for reads of 150, 250 or 300 bases (76 and 101 for 75 and 100).
The instrument corrects the signal of each cycle for phasing (molecules of a cluster that
fall one cycle behind, or run one ahead) with the help of the cycle after it. The extra
cycle is there for that correction of the cycle before it; its own calls get none, and
they are poor. Measured on the MiSeq test sets, by comparing the last base of a read with
the base the other mate of the pair reads at the same place of the template: the last
base of R2 agrees in 60 to 86 % of the pairs (the base before it in 96 to 98 %), the last
base of R1 in 74 to 99 % (97 to 99.7 %). In the ITS2 yeast mocks, where the template
carries a T under R2's cycle 301, that cycle reads C in 79 % of the pairs.

Without this step the base reaches DADA2 whenever truncation does not remove it: always on
ITS, which is not truncated (truncLen 0, 0). A base called wrong in the same way in most
reads of a sequence survives denoising, and mergePairs, which allows no mismatch where the
two mates overlap (amplicon.dada2.merge.max_mismatch, 0 by default), then rejects the pair.
On the other markers truncLen decides what happens to the last base: a cut below the read
length removes it anyway, and a cut at the read-length limit keeps it, but removing it
there changes the read lengths pick_trunclen measures, so the cut itself can move. On the
18S test set it moved from 271 to 274 bp for R1, and 1.2 % fewer reads reached the
tables, with the same merging rate; nothing was gained. So auto leaves truncated reads
alone.

What it does
------------
Reads the first SCAN_READS reads of the file (all of them in a smaller file), takes the
length of the longest, and decides with the setting amplicon.cutadapt.trim_final_cycle:
  - auto: remove the last base when the read direction is not truncated to a fixed length
    (ITS; the other markers only with trunc_len.mode: manual and a cut of 0) and the
    longest read is one of FINAL_CYCLE_LENGTHS (76, 101, 151, 251 or 301 bases, the read
    lengths of a run with the extra cycle), keep it otherwise. The longest read, not every
    read, because a facility may have cut the adapter off the reads that ran past the end
    of a short amplicon, so that many reads are shorter than the run's read length;
  - true: always remove it;
  - false: never remove it, as MetaFlux did up to v2.6.1.
The decision is per file, so per sample and read direction: samples from different runs
in one output folder each get their own.

Input (command line, four words): the setting (auto, true or false, already checked at
start-up by _trim_final_cycle_setting in 00_common.smk), the FASTQ file (rm_phix's output,
gzipped, or the raw reads when PhiX removal is off, gzipped or not), the read direction
(R1 or R2, used in the message only), and whether DADA2 cuts that direction to a fixed
length ("truncated" or "untruncated", worked out at start-up by
_read_direction_truncation in 00_common.smk; read by auto only).
Output: one word on standard output, "remove" or "keep", which the rule turns into
Cutadapt's --cut -1 (R1) or -U -1 (R2); one line on standard error saying why, which the
rule writes into its log (logs/trim_primers/{sample}.log).

Why Python reads the file, and not zcat | head in the rule's shell command: head stops
reading after the reads it needs, zcat is then ended by the signal SIGPIPE, and the strict
mode Snakemake runs every shell command in (bash's set -euo pipefail) counts that as a
failure of the whole command. Python opens the file, reads the first reads and closes it.

What a change to this file costs: rule trim_primers lists it as an input, so editing it,
a comment included, makes primer trimming and every step after it run again in every
amplicon output folder.
"""
import gzip
import sys
from pathlib import Path

# How many reads to read from the start of the file. Every read the facility did not
# shorten has the run's read length, so the longest read is found among the first few
# reads, unless nearly every read was shortened, and then hardly any read carries the
# last cycle at all. Reading 10,000 reads takes well under a second.
SCAN_READS = 10_000

# The read lengths of a run sequenced with one extra cycle: 75, 100, 150, 250 and 300
# bases, each plus one. A run set up without the extra cycle (2 x 300, 2 x 250) gives the
# round length, and auto keeps its last base. Other set-ups with an extra cycle, such as
# 2 x 126 or an uneven split of the cycles of a kit (321 + 281), are not covered by auto;
# trim_final_cycle: true covers them. On the one uneven split among the test sets, a
# NextSeq 2000 run with 321- and 281-bp reads, the last cycle is called as well as the
# one before it, so keeping it there costs nothing.
FINAL_CYCLE_LENGTHS = (76, 101, 151, 251, 301)

# The words the setting can take, checked again here so that a call with a mistyped word
# stops instead of keeping the base without a word.
SETTINGS = ("auto", "true", "false")

# Whether DADA2 cuts the read direction to a fixed length (truncLen), as 00_common.smk
# passes it: "untruncated" for ITS (truncLen 0, 0) and for a manual cut of 0, "truncated"
# otherwise. Checked here too, for the same reason as the setting.
TRUNCATION_STATES = ("truncated", "untruncated")


def open_fastq(path):
    """Open a FASTQ file for reading as text, gzipped or plain.

    Whether a file is gzipped is read from its first two bytes (1f 8b, the mark every
    gzip file starts with), not from its name: the raw reads may be plain FASTQ files.
    """
    with open(path, "rb") as handle:
        first_two_bytes = handle.read(2)
    if first_two_bytes == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def longest_read_length(path, max_reads=SCAN_READS):
    """The length of the longest read among the first `max_reads` reads of a FASTQ file.

    Input: the path of the file, and how many reads to read at most.
    Output: (the longest length in bases, the number of reads read); (0, 0) for an empty
    file. A FASTQ record is four lines: the name, the sequence, a "+" line and the
    quality values; only the sequence line is measured.
    """
    longest = 0
    reads_read = 0
    with open_fastq(path) as handle:
        while reads_read < max_reads:
            name_line = handle.readline()
            if not name_line:
                break                                  # the end of the file
            sequence_line = handle.readline()
            handle.readline()                          # the "+" line
            handle.readline()                          # the quality line
            length = len(sequence_line.rstrip("\r\n"))
            if length > longest:
                longest = length
            reads_read += 1
    return longest, reads_read


def decide(setting, longest, truncation):
    """Whether to remove the last base of every read of the file.

    Input: the setting (auto, true or false), the length of the longest read, and
    whether DADA2 cuts this read direction to a fixed length ("truncated" or
    "untruncated").
    Output: True (remove) or False (keep). true and false do not look at the length or
    the truncation; auto removes the base only from untruncated reads of a run with one
    extra cycle.
    """
    if setting == "true":
        return True
    if setting == "false":
        return False
    if truncation == "truncated":
        return False
    return longest in FINAL_CYCLE_LENGTHS


def explanation(direction, file_name, setting, longest, reads_read, remove, truncation):
    """The line for the rule's log: the read direction, what happens to its last base,
    and why, with the longest read, the setting and, under auto, the truncation."""
    if remove:
        action = "the last base of every read is removed"
    else:
        action = "the last base is kept"

    if reads_read == 0:
        measured = f"{file_name} holds no reads"
    elif reads_read < SCAN_READS:
        measured = f"the longest of the {reads_read:,} reads of {file_name} is {longest} bp"
    else:
        measured = (f"the longest of the first {SCAN_READS:,} reads of {file_name} "
                    f"is {longest} bp")

    if setting != "auto":
        reason = f"trim_final_cycle is {setting}"
    elif truncation == "truncated":
        reason = (f"trim_final_cycle is auto, and DADA2 cuts {direction} to a fixed length "
                  f"(truncLen), which decides what happens to the last base")
    elif remove:
        reason = (f"trim_final_cycle is auto, {direction} is not truncated, and {longest} bp "
                  f"is the read length of a run with one extra cycle")
    else:
        lengths_text = ", ".join(str(length) for length in FINAL_CYCLE_LENGTHS)
        reason = (f"trim_final_cycle is auto, {direction} is not truncated, and "
                  f"{lengths_text} bp are the lengths it acts on")

    return f"[trim_primers] {direction} final cycle: {action} ({measured}; {reason})"


def main(arguments):
    """Read the four command-line words, decide, print the word and the log line."""
    if len(arguments) != 4:
        sys.exit("usage: 30a_final_cycle.py <auto|true|false> <fastq file> <R1|R2> "
                 "<truncated|untruncated>")
    setting, fastq_path, direction, truncation = arguments
    if setting not in SETTINGS:
        sys.exit(f"[trim_primers] trim_final_cycle must be auto, true or false, got {setting!r}")
    if truncation not in TRUNCATION_STATES:
        sys.exit(f"[trim_primers] the truncation of {direction} must be truncated or "
                 f"untruncated, got {truncation!r}")

    longest, reads_read = longest_read_length(fastq_path)
    remove = decide(setting, longest, truncation)
    file_name = Path(fastq_path).name

    # The line for the log first (standard error), then the word for the rule's shell
    # command (standard output), which reads nothing else.
    print(explanation(direction, file_name, setting, longest, reads_read, remove, truncation),
          file=sys.stderr)
    if remove:
        print("remove")
    else:
        print("keep")


if __name__ == "__main__":
    main(sys.argv[1:])
