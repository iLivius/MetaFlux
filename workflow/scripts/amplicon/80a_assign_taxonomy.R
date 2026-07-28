#!/usr/bin/env Rscript
# Assigns a taxonomic lineage to every ASV, using DADA2's implementation of
# the RDP Naive Bayesian classifier (Wang et al. 2007, Applied and
# Environmental Microbiology).
#
# What this is for: for each ASV sequence, assignTaxonomy compares its
# k-mers (short, fixed-length subsequences — 8 bp here, the method's own
# setting) against a reference database of sequences with known taxonomy
# (SILVA for 16S, UNITE for ITS, PR2 for 18S, ...), and works down the rank
# hierarchy (kingdom -> ... -> genus) assigning the best-supported taxon at
# each level. Confidence at each rank comes from BOOTSTRAPPING: the
# classifier repeats the classification 100 times over random subsets of the
# sequence's own k-mers, and the fraction of those 100 repeats that agree is
# the "bootstrap confidence" for that rank. A rank is only kept if its
# confidence is >= min_boot; ranks that don't clear the bar come back NA
# rather than a shaky guess. For 16S only, addSpecies then adds one more,
# final rank on top — see the addSpecies block below for how that step
# differs from assignTaxonomy.
#
# Input: the length-filtered ASV FASTA (seqs_lenfilt.fasta) and its matching
# ASV_#-keyed sequence table (seqtab_lenfilt_head_names.txt), both from rule
# dada_length_filter (60d_dada_length_filter.py) — the FASTA after any
# Metaxa2/ITSx region extraction and the post-extraction ASV length window,
# i.e. this run's final ASV set. Also reads the marker's reference database
# (refdb) and, for 16S only, the addSpecies reference (species_db).
#
# Output: three views of the same classification — asv_table.txt (rows =
# ASV_N ids, one taxonomy string per ASV), asv_table_seqs.txt (same, but rows
# = the raw ASV sequence instead of its short id), and taxon_seq_table.txt
# (rows = ASV_N ids, one column per rank instead of a single semicolon-joined
# string, plus the sequence). All three have the optional keep/discard
# contaminant filter already applied (see further down this file). Of the
# three, only asv_table.txt is read again inside this pipeline: rule
# aggregate_read_counts (80b_aggregate_read_counts.py) sums each sample's
# column as the final "taxonomy" stage of the pipeline-wide read-tracking
# table (read_tracking.txt). The other two are terminal outputs for whatever
# analysis comes after this pipeline.

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, data.table))

# sink() below redirects everything this script prints from here on —
# messages, warnings, DADA2's own progress output — into the rule's log file
# instead of the console, so the log: file Snakemake tracks has the full record.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# ── Parameters ────────────────────────────────────────────────────────────────
amp_type        <- snakemake@params[["amp_type"]]
min_boot        <- as.integer(snakemake@params[["min_boot"]])
try_rc          <- as.logical(snakemake@params[["try_rc"]])
tax_levels      <- unlist(snakemake@params[["tax_levels"]])       # ordered rank columns
rank_prefixes   <- unlist(snakemake@params[["rank_prefixes"]])    # k__/p__/… per level
prefix_style    <- snakemake@params[["prefix_style"]]             # "bare" | "embedded"
filter_enabled  <- as.logical(snakemake@params[["filter_enabled"]])
filter_keep     <- unlist(snakemake@params[["filter_keep"]])      # rank-token keep list
filter_discard  <- unlist(snakemake@params[["filter_discard"]])   # rank-token discard list
threads         <- as.integer(snakemake@threads)
seed            <- as.integer(snakemake@params[["seed"]])

# assignTaxonomy's RDP classifier draws bootstrap resamples via R's RNG to
# compute minBoot confidence; fix it so calls near the threshold don't flip.
set.seed(seed)

message("[assign_taxonomy] amp_type=", amp_type, " (", prefix_style, " ranks)")

# ── Read FASTA (simple parser — no seqinr dependency) ─────────────────────────
# Plain-text FASTA parsing, no extra package needed: find every header line
# (starts with ">"), take the id as everything up to the first space, and
# join all the lines between one header and the next into a single sequence
# string. Works whether the FASTA wraps sequences over multiple lines or
# writes one line per record.
read_fasta_simple <- function(path) {
  lines <- readLines(path)
  header_idx <- which(startsWith(lines, ">"))
  ids <- sub("^>([^ ]+).*", "\\1", lines[header_idx])
  starts <- header_idx + 1L
  ends <- c(header_idx[-1L] - 1L, length(lines))
  seqs <- vapply(seq_along(ids), function(i) {
    paste(lines[starts[i]:ends[i]], collapse = "")
  }, character(1L))
  setNames(seqs, ids)
}

fasta <- read_fasta_simple(snakemake@input[["seqs"]])
asv_ids  <- names(fasta)
seqs_vec <- unname(fasta)
names(seqs_vec) <- asv_ids
message("[assign_taxonomy] FASTA: ", length(fasta), " ASVs loaded from ", snakemake@input[["seqs"]])

# ── Read seqtab (ASV_N cols) ──────────────────────────────────────────────────
seqtab_raw <- fread(snakemake@input[["seqtab_names"]], sep = "\t",
                    data.table = FALSE, showProgress = FALSE)
rownames(seqtab_raw) <- seqtab_raw[[1L]]
seqtab_raw <- seqtab_raw[, -1L, drop = FALSE]

# Subset to ASVs present in the FASTA (length filter / extraction may have removed some)
shared_ids <- intersect(asv_ids, colnames(seqtab_raw))
seqtab_sub <- seqtab_raw[, shared_ids, drop = FALSE]
seqs_sub   <- seqs_vec[shared_ids]
message("[assign_taxonomy] Seqtab: ", ncol(seqtab_sub), " ASVs shared with FASTA")

# ── Taxonomy assignment ───────────────────────────────────────────────────────
# ITSx extraction can produce identical sequences for distinct ASVs (same ITS
# region recovered from different full-length inputs). Run assignTaxonomy on the
# unique sequences only, then expand back so every ASV ID gets its own row.
refdb_path <- snakemake@input[["refdb"]]
message("[assign_taxonomy] assignTaxonomy against ", refdb_path)
unique_seqs <- unique(unname(seqs_sub))
message("[assign_taxonomy] Classifying ", length(unique_seqs),
        " unique sequence(s) from ", length(seqs_sub), " input ASVs")
# taxLevels comes from the marker profile so the classifier matches the reference's
# actual depth. For 16S/ITS these are the Linnaean 7, identical to DADA2's own
# default (so passing them explicitly changes nothing); a 9-rank reference such as
# PR2 (domain..species) needs them declared — assignTaxonomy sizes the output from
# the reference's own depth, so undeclared levels arrive as unnamed columns AND the
# names that were declared slide onto the wrong ranks (PR2's Family would be labelled
# Species). make_tax_string below looks values up by name, so that misnaming is the
# damaging half, not the unnamed tail.
taxa_unique <- assignTaxonomy(unique_seqs, refdb_path,
                              minBoot = min_boot,
                              multithread = threads,
                              tryRC = try_rc,
                              taxLevels = tax_levels,
                              verbose = TRUE)

# addSpecies for 16S only; species_db is empty list for ITS
species_db_list <- snakemake@input[["species_db"]]
if (length(species_db_list) > 0L) {
  message("[assign_taxonomy] addSpecies against ", species_db_list[[1]])
  taxa_unique <- addSpecies(taxa_unique, species_db_list[[1]])
}

# Expand back: index by each ASV's sequence so identical-sequence ASVs share rows
taxa <- taxa_unique[unname(seqs_sub), , drop = FALSE]
taxa_df <- as.data.frame(taxa, stringsAsFactors = FALSE)
rownames(taxa_df) <- names(seqs_sub)
message("[assign_taxonomy] Taxa assigned for ", nrow(taxa_df), " ASVs")

# ── Build taxonomy string (registry-driven; see 00_common.smk MARKERS) ────────
# prefix_style == "bare"     → values carry no rank prefix: prepend rank_prefixes,
#                              and render the last (species) slot as a
#                              Genus+species binomial — but ONLY when the value
#                              looks like a bare epithet (SILVA/16S addSpecies
#                              output, e.g. "coli": lowercase, nomenclatural
#                              convention for species epithets). A reference that
#                              already hands back a complete Genus_species token
#                              (PR2/18S: "Unruhdinium_kevei") starts uppercase, so
#                              it is emitted as-is instead of gluing the genus on
#                              a second time ("Unruhdinium Unruhdinium_kevei").
#                              Same lowercase-epithet heuristic already used on
#                              the shotgun side — see 45a_finalize_otu_table.py.
# prefix_style == "embedded" → values already carry k__/p__… prefixes; emit
#                              them verbatim (UNITE/ITS behaviour).
make_tax_string <- function(row) {
  parts <- character(0)
  n <- length(tax_levels)
  for (i in seq_len(n)) {
    val <- row[[tax_levels[i]]]
    if (is.na(val) || val == "") next
    if (identical(prefix_style, "embedded")) {
      parts <- c(parts, val)
    } else if (i == n && n >= 2L) {
      # species slot → binomial; requires the genus (level n-1) to be present
      gval <- row[[tax_levels[n - 1L]]]
      if (is.na(gval) || gval == "") next
      # Glue the genus on only for a bare, lowercase epithet. A value that
      # already starts uppercase, or already begins with the genus (PR2 uses
      # both "Genus_species" binomials and, for ~18 unresolved environmental
      # clades, lowercase placeholders already prefixed with their own clade
      # name, e.g. genus "alphaI_cluster" / species "alphaI_cluster_sp." —
      # confirmed against the full shipped pr2_SSU_dada2.fasta.gz: every
      # lowercase-starting species value in it already starts with its own
      # genus), is a complete name as-is — gluing again would double it.
      needs_genus <- grepl("^[a-z]", val) && !startsWith(val, gval)
      if (needs_genus) {
        parts <- c(parts, paste0(rank_prefixes[i], gval, " ", val))
      } else {
        parts <- c(parts, paste0(rank_prefixes[i], val))
      }
    } else {
      parts <- c(parts, paste0(rank_prefixes[i], val))
    }
  }
  paste(parts, collapse = ";")
}
tax_strings <- apply(taxa_df, 1L, make_tax_string)

# ── ASV table: rows=ASV_IDs, cols=[samples..., taxonomy] ─────────────────────
# seqtab_sub has rows=samples, cols=ASV_IDs → transpose for ASV table orientation
count_t <- as.data.frame(t(seqtab_sub[, rownames(taxa_df), drop = FALSE]))
asv_table <- cbind(count_t, taxonomy = tax_strings)

# ── Contaminant filter (rank-aware keep/discard lists; config-only) ──────────
# Each token (e.g. k__Bacteria, o__Chloroplast) is matched against whole
# ';'-delimited segments of the taxonomy string, so g__Bacillus hits only a
# genus named exactly that — never a substring elsewhere. keep runs first
# (defines the universe), then discard prunes. An empty list = that direction
# is a no-op. There is NO per-marker code default: the lists come from config.
seg_match <- function(tax_string, tokens) {
  segs <- trimws(strsplit(tax_string, ";", fixed = TRUE)[[1]])
  any(tokens %in% segs)
}
n_before_filter <- nrow(asv_table)
if (isTRUE(filter_enabled)) {
  if (length(filter_keep) > 0L) {
    keep <- vapply(asv_table$taxonomy, seg_match, logical(1L),
                   tokens = filter_keep, USE.NAMES = FALSE)
    asv_table <- asv_table[keep, , drop = FALSE]
    message("[assign_taxonomy] Keep filter (", paste(filter_keep, collapse = ","), "): ",
            sum(keep), " / ", n_before_filter, " ASVs retained")
  }
  if (length(filter_discard) > 0L) {
    drop <- vapply(asv_table$taxonomy, seg_match, logical(1L),
                   tokens = filter_discard, USE.NAMES = FALSE)
    asv_table <- asv_table[!drop, , drop = FALSE]
    message("[assign_taxonomy] Discard filter (", paste(filter_discard, collapse = ","), "): ",
            sum(drop), " ASV(s) removed")
  }
}
message("[assign_taxonomy] Final ASV count: ", nrow(asv_table))

# Fail loud rather than write a header-only table: a filter that removes 100% of
# ASVs is almost always a marker/config mismatch (e.g. 16S keep tokens on an ITS run).
if (isTRUE(filter_enabled) && n_before_filter > 0L && nrow(asv_table) == 0L) {
  stop("[assign_taxonomy] the keep/discard filter removed ALL ", n_before_filter,
       " ASVs (amp_type=", amp_type, "). The keep/discard lists likely do not match this ",
       "marker: for ITS use keep: [k__Fungi]; for 16S use keep: [k__Bacteria, k__Archaea]. ",
       "Fix amplicon.taxonomy.filter, or set amplicon.taxonomy.filter.enabled: false.")
}

# ── Output (a): asv_table.txt — rows=ASV_IDs ─────────────────────────────────
dir.create(dirname(snakemake@output[["asv_table"]]), recursive = TRUE, showWarnings = FALSE)
write.table(asv_table, file = snakemake@output[["asv_table"]], col.names = NA, sep = "\t")

# ── Output (b): asv_table_seqs.txt — rows=sequences ──────────────────────────
asv_table_seqs <- asv_table
rownames(asv_table_seqs) <- seqs_sub[rownames(asv_table)]
write.table(asv_table_seqs, file = snakemake@output[["asv_table_seqs"]], col.names = NA, sep = "\t")

# ── Output (c): taxon_seq_table.txt — rows=ASV_IDs, parsed cols + sequence ───
taxon_table <- taxa_df[rownames(asv_table), , drop = FALSE]
taxon_table$sequence <- seqs_sub[rownames(asv_table)]
write.table(taxon_table, file = snakemake@output[["taxon_seq_table"]], col.names = NA, sep = "\t")

message("[assign_taxonomy] DONE. Outputs written to ", dirname(snakemake@output[["asv_table"]]))

# The main table: one row per ASV, one column per sample (read counts) plus a
# taxonomy column. Its sample columns are summed again downstream by rule
# aggregate_read_counts to build the final stage of read_tracking.txt.
sink(type = "message")
sink(type = "output")
close(log_con)
# Identical to (a), except each row is named by the ASV's actual sequence
# instead of its short ASV_N id — for matching these ASVs directly against
# sequences from another tool or database, without going through the id.
# Same ASVs and taxonomy as (a), but with the taxonomy split into one column
# per rank (Kingdom, Phylum, ...) instead of one semicolon-joined string, plus
# the ASV's sequence — convenient for filtering or grouping by a single rank
# without string-parsing the taxonomy column first. Has no per-sample counts.
