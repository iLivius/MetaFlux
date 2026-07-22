#!/usr/bin/env Rscript
# Taxonomy assignment for DADA2 ASVs via DADA2's RDP classifier.
# assignTaxonomy (SILVA/UNITE DADA2 format) + addSpecies (16S only).
# Reads the (possibly Metaxa2/ITSx-extracted) seqs FASTA and seqtab_head_names.txt.
# Writes three output tables + applies optional contaminant filter.

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org/")
suppressPackageStartupMessages(pacman::p_load(dada2, data.table))

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# ── Parameters ────────────────────────────────────────────────────────────────
amp_type        <- snakemake@params[["amp_type"]]
min_boot        <- as.integer(snakemake@params[["min_boot"]])
try_rc          <- as.logical(snakemake@params[["try_rc"]])
filter_enabled  <- as.logical(snakemake@params[["filter_enabled"]])
include_pattern <- snakemake@params[["include_pattern"]]
exclude_pattern <- snakemake@params[["exclude_pattern"]]
threads         <- as.integer(snakemake@threads)
seed            <- as.integer(snakemake@params[["seed"]])

# assignTaxonomy's RDP classifier draws bootstrap resamples via R's RNG to
# compute minBoot confidence; fix it so calls near the threshold don't flip.
set.seed(seed)

# Default contaminant filter patterns when config leaves them null
if (is.null(include_pattern) || include_pattern == "") {
  include_pattern <- if (amp_type == "16S") "k__Archaea|k__Bacteria" else "k__Fungi"
}
if (is.null(exclude_pattern) || exclude_pattern == "") {
  exclude_pattern <- if (amp_type == "16S") "chloroplast|mitochondria" else NULL
}

message("[assign_taxonomy] amp_type=", amp_type)

# ── Read FASTA (simple parser — no seqinr dependency) ─────────────────────────
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
taxa_unique <- assignTaxonomy(unique_seqs, refdb_path,
                              minBoot = min_boot,
                              multithread = threads,
                              tryRC = try_rc,
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

# ── Build taxonomy string ─────────────────────────────────────────────────────
make_tax_string_16S <- function(row) {
  parts <- character(0)
  if (!is.na(row["Kingdom"])) parts <- c(parts, paste0("k__", row["Kingdom"]))
  if (!is.na(row["Phylum"]))  parts <- c(parts, paste0("p__", row["Phylum"]))
  if (!is.na(row["Class"]))   parts <- c(parts, paste0("c__", row["Class"]))
  if (!is.na(row["Order"]))   parts <- c(parts, paste0("o__", row["Order"]))
  if (!is.na(row["Family"]))  parts <- c(parts, paste0("f__", row["Family"]))
  if (!is.na(row["Genus"]))   parts <- c(parts, paste0("g__", row["Genus"]))
  if (!is.na(row["Genus"]) && !is.na(row["Species"]))
    parts <- c(parts, paste0("s__", row["Genus"], " ", row["Species"]))
  paste(parts, collapse = ";")
}

make_tax_string_ITS <- function(row) {
  # UNITE columns already carry k__/p__/... prefixes; NAs are plain R NA
  lvls <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  parts <- row[lvls]
  parts <- parts[!is.na(parts) & parts != ""]
  s <- paste(parts, collapse = ";")
  # Strip any residual ";NA" or "; NA" artifacts
  s <- gsub(";NA|; NA", "", s)
  s
}

if (amp_type == "16S") {
  tax_strings <- apply(taxa_df, 1L, make_tax_string_16S)
} else {
  tax_strings <- apply(taxa_df, 1L, make_tax_string_ITS)
}

# ── ASV table: rows=ASV_IDs, cols=[samples..., taxonomy] ─────────────────────
# seqtab_sub has rows=samples, cols=ASV_IDs → transpose for ASV table orientation
count_t <- as.data.frame(t(seqtab_sub[, rownames(taxa_df), drop = FALSE]))
asv_table <- cbind(count_t, taxonomy = tax_strings)

# ── Contaminant filter ────────────────────────────────────────────────────────
n_before_filter <- nrow(asv_table)
if (isTRUE(filter_enabled)) {
  if (!is.null(include_pattern) && nchar(include_pattern) > 0L) {
    keep <- grepl(include_pattern, asv_table$taxonomy, ignore.case = TRUE)
    asv_table <- asv_table[keep, , drop = FALSE]
    message("[assign_taxonomy] Include filter (", include_pattern, "): ",
            sum(keep), " / ", n_before_filter, " ASVs retained")
  }
  if (!is.null(exclude_pattern) && nchar(exclude_pattern) > 0L) {
    drop <- grepl(exclude_pattern, asv_table$taxonomy, ignore.case = TRUE)
    asv_table <- asv_table[!drop, , drop = FALSE]
    message("[assign_taxonomy] Exclude filter (", exclude_pattern, "): ",
            sum(drop), " ASV(s) removed")
  }
}
message("[assign_taxonomy] Final ASV count: ", nrow(asv_table))

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

sink(type = "message")
sink(type = "output")
close(log_con)
