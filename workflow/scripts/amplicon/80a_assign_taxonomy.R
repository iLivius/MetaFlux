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
# (refdb) and, for 16S only, the addSpecies reference (species_db) and
# Metaxa2's own origin call per ASV (metaxa2_origin, 5.dada2/metaxa2_origin.tsv
# from target_extract — present only when Metaxa2 ran), which the organelle
# filter uses (see "Organelle filter" below).
#
# Output: three views of the same classification — asv_table.txt (rows =
# ASV_N ids, one taxonomy string per ASV), asv_table_seqs.txt (same, but rows
# = the raw ASV sequence instead of its short id), and taxon_seq_table.txt
# (rows = ASV_N ids, one column per rank instead of a single semicolon-joined
# string, plus the sequence). All three have the optional keep/discard
# contaminant filter already applied (see further down this file). Of the
# three, asv_table.txt is read again inside this pipeline: rule
# aggregate_read_counts (80b_aggregate_read_counts.py) sums each sample's
# column as the final "taxonomy" stage of the pipeline-wide read-tracking
# table (read_tracking.txt), and the phylogeny (85a_phylo_input.py) reads it
# together with taxon_seq_table.txt. A fourth table, asv_table_unfiltered.txt,
# holds EVERY ASV that reached the filter, with two extra columns
# (metaxa2_origin, filter_status) — for the user's own analyses, e.g. how many
# organelle reads two primer sets pick up; nothing in the pipeline reads it.

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
# Best-hit check of the discard filter (see "Discard filter, best-hit check" below):
# minimum bootstrap support, 0-100, at the token's rank; NULL = check switched off.
discard_min_support <- snakemake@params[["discard_min_support"]]
# Organelle filter (see "Organelle filter" below): TRUE when the ASVs Metaxa2 calls
# mitochondria or chloroplast are to be removed. Decided in 00_common.smk: only ever
# TRUE for 16S with extraction on, filter.enabled true and filter.organelles: metaxa2.
organelle_filter <- isTRUE(as.logical(snakemake@params[["organelle_filter"]]))
# Metaxa2's per-ASV origin table: a path when Metaxa2 ran in this run (16S with
# extraction on), NULL for every other marker — see rule assign_taxonomy.
origin_path     <- snakemake@input[["metaxa2_origin"]]
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
#
# assignTaxonomy is called with minBoot = 0 and outputBootstraps = TRUE, so it hands
# back the classifier's best-hit lineage at EVERY rank together with the bootstrap
# support of each rank. min_boot is then applied here, by apply_min_boot() below,
# with exactly the rule DADA2 itself uses, so the reported taxonomy is the same as
# a call with minBoot = min_boot (the bootstrap draws do not depend on minBoot).
# Keeping the unfiltered best hit is what lets the discard filter look one step
# further than the reported lineage (see "Discard filter, best-hit check" below).
classified <- assignTaxonomy(unique_seqs, refdb_path,
                             minBoot = 0,
                             outputBootstraps = TRUE,
                             multithread = threads,
                             tryRC = try_rc,
                             taxLevels = tax_levels,
                             verbose = TRUE)
best_hit  <- classified$tax     # best-hit name at every rank, whatever its support
best_boot <- classified$boot    # bootstrap support (0-100) of each of those names

# ── UNITE species names: put the genus back (ITS) ────────────────────────────
# What this is for: DADA2 shortens UNITE's species names. When assignTaxonomy
# recognises the UNITE general release (by "FU|refs" or "FU|reps" in its first ten
# headers), it cuts the genus out of every species field before classifying, so
# UNITE's s__Hanseniaspora_uvarum comes back as s__uvarum. The sintax path keeps
# UNITE's whole name, and so does every other marker's reference. Without this
# step an ITS species would read differently under the two methods, and a species
# token such as s__Hanseniaspora_uvarum would match nothing under rdp.
#
# What it does: for each classified sequence, the best-hit genus is put back in
# front of the best-hit species (g__Hanseniaspora + s__uvarum ->
# s__Hanseniaspora_uvarum) when all of these hold:
#   - the marker's reference values carry their own prefixes (prefix_style:
#     embedded), as UNITE's do;
#   - DADA2 did shorten this reference, decided with DADA2's own header test
#     (reference_is_unite below). So the gyrB and rpoB references are never
#     touched, although a few rpoB species also start lowercase
#     (s__endosymbiont_of_Sipalinus_gigas);
#   - both names are known: a "Genus_sp" placeholder species, and an
#     "unidentified" genus, come back from DADA2 as NA;
#   - the species is a plain epithet: it starts with a lowercase letter and has
#     no "." in it.
# It is not put back for the var./subsp./f. names, about 1 % of UNITE's named
# species. For those DADA2 keeps only the part from "var." (or "subsp.", "f.")
# onwards: s__Diaporthe_phaseolorum_var._meridionalis comes back as
# s__var._meridionalis. The species epithet in front of it is lost, so the name
# cannot be rebuilt; it is left as DADA2 wrote it. A name DADA2 did not shorten is
# whole already and starts uppercase, so it is left as it is too: DADA2's cut skips
# a hyphenated genus (s__Pseudo-nitzschia_limii).
#
# Why here, in the best-hit matrix: everything later reads this matrix. The
# reported lineage is cut from it by apply_min_boot() just below, and the best-hit
# check of the discard filter compares each discard token with it (see "Discard
# filter, best-hit check"). Fixed once here, both see s__Hanseniaspora_uvarum.
# Input: best_hit (above). Output: best_hit with its Species column rebuilt.

# TRUE when DADA2 treats the reference file at `path` as UNITE: the test in
# dada2::assignTaxonomy, "FU|refs" or "FU|reps" in each of the first ten headers.
# Only the start of the file is read (readLines also opens a gzipped file).
reference_is_unite <- function(path) {
  first_lines <- readLines(path, n = 2000L)
  headers     <- first_lines[startsWith(first_lines, ">")]
  if (length(headers) < 10L) return(FALSE)
  all(grepl("FU\\|re[pf]s", headers[1:10]))
}

# One best-hit species value with the genus put back, or the value unchanged when
# it cannot be rebuilt (see the list above). With the prefixes s__ and g__:
#   ("s__uvarum", "g__Hanseniaspora")         -> "s__Hanseniaspora_uvarum"
#   ("s__var._meridionalis", "g__Diaporthe")  -> "s__var._meridionalis" (unchanged)
put_genus_back <- function(species_val, genus_val, species_pfx, genus_pfx) {
  if (is.na(species_val) || is.na(genus_val)) return(species_val)
  if (!startsWith(species_val, species_pfx) || !startsWith(genus_val, genus_pfx)) {
    return(species_val)
  }
  epithet <- substring(species_val, nchar(species_pfx) + 1L)
  genus   <- substring(genus_val, nchar(genus_pfx) + 1L)
  is_plain_epithet <- grepl("^[a-z]", epithet) && !grepl(".", epithet, fixed = TRUE)
  if (!is_plain_epithet) return(species_val)
  paste0(species_pfx, genus, "_", epithet)
}

# Only when the matrix has exactly the pack's ranks: with any other depth the
# reference depth check below stops the run with its own message.
n_ranks <- length(tax_levels)
if (identical(prefix_style, "embedded") && ncol(best_hit) == n_ranks &&
    reference_is_unite(refdb_path)) {
  species_col <- tax_levels[n_ranks]
  genus_col   <- tax_levels[n_ranks - 1L]
  n_named     <- sum(!is.na(best_hit[, species_col]))
  n_rebuilt   <- 0L
  for (i in seq_len(nrow(best_hit))) {
    old_val <- best_hit[i, species_col]
    new_val <- put_genus_back(old_val, best_hit[i, genus_col],
                              rank_prefixes[n_ranks], rank_prefixes[n_ranks - 1L])
    if (!identical(new_val, old_val)) {
      best_hit[i, species_col] <- new_val
      n_rebuilt <- n_rebuilt + 1L
    }
  }
  n_kept_as_is <- n_named - n_rebuilt
  message("[assign_taxonomy] UNITE reference: genus put back in front of the best-hit ",
          "species epithet of ", n_rebuilt, " sequence(s); ", n_kept_as_is, " other species ",
          "name(s) kept as DADA2 returned them (var./subsp./f. names, or names DADA2 did ",
          "not shorten)")
}

# DADA2's own minBoot rule (dada2::assignTaxonomy, taxonomy.R): for each sequence,
# keep the ranks whose support reaches min_boot and write them from the first
# column onwards. Support can only fall from domain to genus, so in practice this
# keeps the lineage down to the last well-supported rank and blanks the rest.
apply_min_boot <- function(tax, boot, min_boot) {
  out <- matrix(NA_character_, nrow = nrow(tax), ncol = ncol(tax), dimnames = dimnames(tax))
  for (i in seq_len(nrow(tax))) {
    kept <- tax[i, boot[i, ] >= min_boot]
    if (length(kept) > 0L) out[i, seq_along(kept)] <- kept
  }
  out
}
taxa_unique <- apply_min_boot(best_hit, best_boot, min_boot)

# ── Reference depth check ─────────────────────────────────────────────────────
# assignTaxonomy returns one column per rank of the DEEPEST lineage in the
# reference file, and simply names them from taxLevels in order — it never checks
# that the reference's ranks are the ones the names describe. A reference one rank
# deeper than the marker pack expects therefore comes back with every name shifted
# one slot (SILVA 144's own file, with its extra Kingdom rank, would put each
# family in the Genus column and each genus in Species), and nothing errors.
# So compare the column count with the pack's rank model here, before anything is
# built from it. For 16S the last slot (Species) is filled later by addSpecies, so
# the training set itself must stop one rank short of the full model; every other
# marker's reference covers all of its ranks.
species_db_list <- snakemake@input[["species_db"]]
expected_ranks  <- length(tax_levels) - if (length(species_db_list) > 0L) 1L else 0L
if (ncol(taxa_unique) != expected_ranks) {
  stop("[assign_taxonomy] the reference ", refdb_path, " has ", ncol(taxa_unique),
       " ranks, but the marker pack's rank model (", paste(tax_levels, collapse = ","),
       ") expects ", expected_ranks, " in this file",
       if (length(species_db_list) > 0L) " (the Species slot comes from addSpecies)" else "",
       ". The names would land on the wrong ranks. For SILVA 144 use the 6-rank file ",
       "written by convert_silva_train, not the download itself.")
}

# addSpecies for 16S only; species_db is an empty list for every other marker
# (ITS, 18S, gyrB, rpoB — see 80_taxonomy.smk), so this block just does
# nothing for them. addSpecies works completely differently from
# assignTaxonomy above: it is an EXACT sequence match against a separate
# species-level reference, not the bootstrapped RDP classifier, so there is no
# confidence score at this step — a species is only assigned when the match
# is unambiguous (addSpecies' allowMultiple defaults to FALSE, so ties are
# left blank rather than guessed), and even then only kept if the genus
# assignTaxonomy already gave that ASV agrees with the matched reference's
# genus. tryRC is passed here too, using the same amplicon.taxonomy.try_rc
# setting as assignTaxonomy above, so an ASV whose stored orientation is the
# reverse complement of the reference gets the same chance at a species match
# as it did at a genus match — otherwise a sequence that only classified
# because assignTaxonomy tried its reverse complement could reach addSpecies
# and silently fail to match anything, for no reason a user could see from
# the config.
#
# SILVA 144 "--other" genera (see silva_other in config.yaml): assignTaxonomy may
# call an ASV "Pseudomonas--other" — a Pseudomonas clade that does not contain the
# type strain — while the species reference only ever writes "Pseudomonas
# fluorescens". addSpecies' genus check is an exact string comparison, so it would
# reject every species match for such an ASV. The check is therefore run on a COPY
# of the table whose Genus column has the suffix removed, and only the resulting
# Species column is taken back; the genus itself stays exactly as SILVA wrote it.
# When silva_other is "strip" (or under SILVA 138.2) no name carries the suffix and
# the copy is identical to the original.
if (length(species_db_list) > 0L) {
  message("[assign_taxonomy] addSpecies against ", species_db_list[[1]])
  genus_col <- tax_levels[length(tax_levels) - 1L]
  taxa_for_species <- taxa_unique
  taxa_for_species[, genus_col] <- sub("--other$", "", taxa_for_species[, genus_col])
  with_species <- addSpecies(taxa_for_species, species_db_list[[1]], tryRC = try_rc)
  taxa_unique <- cbind(taxa_unique, Species = with_species[, "Species"])
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
      # A species name is a binomial, so it is written with the plain genus even
      # when the genus slot reads "Genus--other" (SILVA 144, silva_other:
      # verbatim): g__Pseudomonas--other;s__Pseudomonas fluorescens.
      gval <- sub("--other$", "", gval)
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
#
# SILVA 144 "--other" names (silva_other: verbatim) count as the plain name here:
# the suffix is removed from both the segments and the tokens before comparing, so
# g__Bradyrhizobium removes "g__Bradyrhizobium--other" too. A filter token is meant
# to name a taxon, and SILVA considers the "--other" clade part of that taxon. With
# silva_other: strip, or SILVA 138.2, no name carries the suffix and this changes
# nothing.
drop_other <- function(x) sub("--other$", "", x)
seg_match <- function(tax_string, tokens) {
  segs <- drop_other(trimws(strsplit(tax_string, ";", fixed = TRUE)[[1]]))
  any(drop_other(tokens) %in% segs)
}

# Every ASV that reaches the filter, kept aside before anything is removed: this
# is what asv_table_unfiltered.txt is written from (Output (d) below).
asv_table_all <- asv_table

# ── Metaxa2 origin per ASV (16S only; empty otherwise) ──
# metaxa2_origin.tsv has one row per ASV Metaxa2 extracted: asv_id, origin
# (mitochondria, chloroplast, bacteria, archaea, eukaryota, uncertain, or NA for
# no call; read.delim turns the text NA into a missing value). Used by the
# organelle filter below and written, for every ASV, to the metaxa2_origin column
# of asv_table_unfiltered.txt. Looking up an ASV that has no entry gives NA.
metaxa2_origin <- setNames(character(0), character(0))
if (!is.null(origin_path)) {
  origin_tbl     <- read.delim(origin_path, colClasses = "character")
  metaxa2_origin <- setNames(origin_tbl$origin, origin_tbl$asv_id)
  message("[assign_taxonomy] Metaxa2 origin calls read for ", sum(!is.na(metaxa2_origin)),
          " ASVs from ", origin_path)
}
# The rule only switches the filter on when Metaxa2 ran, so the table is always
# there when it is needed; this guards against the two drifting apart.
if (organelle_filter && is.null(origin_path)) {
  stop("[assign_taxonomy] the organelle filter is on but no Metaxa2 origin table was ",
       "given (input metaxa2_origin of rule assign_taxonomy)")
}

# Total reads, over all samples, of the given ASVs (seqtab_sub: rows = samples,
# columns = ASV ids), as text for the log lines; format() keeps a round number such
# as 200000 from being printed as 2e+05.
reads_of <- function(ids) format(sum(seqtab_sub[, ids, drop = FALSE]), scientific = FALSE)

# What the filter did to each ASV, for the filter_status column of
# asv_table_unfiltered.txt. The steps below run in a fixed order — keep, discard,
# best-hit check, organelle filter — and each one only looks at the ASVs the
# earlier steps kept. So a removed ASV is labelled with the FIRST step that
# removed it: removed_keep, removed_discard, removed_best_hit or
# removed_organelle. removed_organelle therefore means "only Metaxa2 caught it".
filter_status <- setNames(rep("kept", nrow(asv_table)), rownames(asv_table))

n_before_filter <- nrow(asv_table)
if (isTRUE(filter_enabled)) {
  if (length(filter_keep) > 0L) {
    keep <- vapply(asv_table$taxonomy, seg_match, logical(1L),
                   tokens = filter_keep, USE.NAMES = FALSE)
    filter_status[rownames(asv_table)[!keep]] <- "removed_keep"
    asv_table <- asv_table[keep, , drop = FALSE]
    message("[assign_taxonomy] Keep filter (", paste(filter_keep, collapse = ","), "): ",
            sum(keep), " / ", n_before_filter, " ASVs retained")
  }
  if (length(filter_discard) > 0L) {
    drop <- vapply(asv_table$taxonomy, seg_match, logical(1L),
                   tokens = filter_discard, USE.NAMES = FALSE)
    filter_status[rownames(asv_table)[drop]] <- "removed_discard"
    asv_table <- asv_table[!drop, , drop = FALSE]
    message("[assign_taxonomy] Discard filter (", paste(filter_discard, collapse = ","), "): ",
            sum(drop), " ASV(s) removed")
  }

  # ── Discard filter, best-hit check ─────────────────────────────────────────
  # The reported taxonomy stops at the last rank whose bootstrap support reaches
  # min_boot (80 by default). An ASV can therefore have Mitochondria as its best
  # hit at family level with, say, 62% support, be reported only as
  # k__Bacteria;p__Pseudomonadota, and slip past f__Mitochondria. Organelle 16S
  # is the typical case: plant and algal mitochondria are diverse, and with
  # SILVA 144 the RDP classifier often stays below 80% at order and family for
  # them although the best hit is clearly mitochondrial. Such ASVs then end up in
  # the abundance table and, as very long branches, in the tree.
  # So, for every ASV still in the table, each discard token is also compared with
  # the best-hit name at the token's own rank (o__ = Order, f__ = Family, ...), and
  # the ASV is removed when the names match and the support at that rank is at
  # least amplicon.taxonomy.filter.discard_min_support (default 50). The keep list
  # is not affected, and a token at a rank without bootstrap support (s__, filled
  # by addSpecies) is only ever matched against the reported taxonomy.
  if (length(filter_discard) > 0L && !is.null(discard_min_support) && nrow(asv_table) > 0L) {
    asv_seqs <- unname(seqs_sub[rownames(asv_table)])
    weak_hits <- character(0)
    for (j in seq_along(asv_seqs)) {
      for (tok in filter_discard) {
        # The rank prefix this token starts with (longest match, so "sg__" is not
        # read as "s__"), and the rank column it stands for.
        hits <- rank_prefixes[nchar(rank_prefixes) > 0L & startsWith(tok, rank_prefixes)]
        if (length(hits) == 0L) next
        pfx   <- hits[which.max(nchar(hits))]
        level <- tax_levels[match(pfx, rank_prefixes)]
        if (!(level %in% colnames(best_hit))) next
        name    <- best_hit[asv_seqs[j], level]
        support <- best_boot[asv_seqs[j], level]
        if (is.na(name) || is.na(support) || support < discard_min_support) next
        # Under "embedded" style (UNITE) the reference value already carries its
        # prefix (o__...), so compare the whole token; otherwise only the name.
        wanted <- if (identical(prefix_style, "embedded")) tok else substring(tok, nchar(pfx) + 1L)
        if (drop_other(name) == drop_other(wanted)) {
          weak_hits[rownames(asv_table)[j]] <- paste0(tok, " at ", support, "% support")
          break
        }
      }
    }
    if (length(weak_hits) > 0L) {
      for (id in names(weak_hits)) {
        message("[assign_taxonomy]   ", id, ": best hit ", weak_hits[[id]],
                " (reported as ", asv_table[id, "taxonomy"], ")")
      }
      filter_status[names(weak_hits)] <- "removed_best_hit"
      asv_table <- asv_table[!(rownames(asv_table) %in% names(weak_hits)), , drop = FALSE]
    }
    message("[assign_taxonomy] Discard filter, best-hit check (support >= ",
            discard_min_support, "%): ", length(weak_hits), " more ASV(s) removed")
  }

  # ── Organelle filter (Metaxa2's own origin call; 16S) ──────────────────────
  # Removes every ASV that Metaxa2 placed in its mitochondria or chloroplast FASTA
  # (amplicon.taxonomy.filter.organelles: metaxa2, the default). It runs last, so
  # it only removes what keep, discard and the best-hit check let through. That is
  # exactly the case it is for: plant mitochondria whose classifier best hit is a
  # bacterium, so no token can see them. On the 16S test set, 8 ASVs with a
  # mitochondrion as nearest SILVA reference (97-99 % identity) were reported as
  # Enterobacterales by the RDP classifier on both SILVA releases and passed every
  # token, while Metaxa2 called all 8 mitochondrial.
  # Metaxa2's other calls (eukaryota, uncertain, the domains) remove nothing; the
  # keep tokens decide domains. The log reports how many organelle ASVs and reads
  # Metaxa2 found, how many of them the earlier steps had already removed, and
  # lists each ASV that only Metaxa2 caught.
  if (organelle_filter) {
    all_ids       <- rownames(asv_table_all)
    origin_all    <- unname(metaxa2_origin[all_ids])          # NA where no call
    organelle_ids <- all_ids[origin_all %in% c("mitochondria", "chloroplast")]
    n_mito        <- sum(metaxa2_origin[organelle_ids] == "mitochondria")
    n_chloro      <- length(organelle_ids) - n_mito
    # Organelle ASVs the earlier steps had already removed, per step.
    status_now    <- unname(filter_status[organelle_ids])
    already_keep  <- organelle_ids[status_now == "removed_keep"]
    already_disc  <- organelle_ids[status_now == "removed_discard"]
    already_best  <- organelle_ids[status_now == "removed_best_hit"]
    already_ids   <- c(already_keep, already_disc, already_best)
    # Organelle ASVs still in the table: only Metaxa2 caught these.
    only_metaxa2  <- organelle_ids[status_now == "kept"]
    filter_status[only_metaxa2] <- "removed_organelle"
    asv_table <- asv_table[!(rownames(asv_table) %in% only_metaxa2), , drop = FALSE]

    message("[assign_taxonomy] Organelle filter (Metaxa2 origin): ", length(organelle_ids),
            " ASV(s), ", reads_of(organelle_ids), " reads, called mitochondria (", n_mito,
            ") or chloroplast (", n_chloro, ")")
    message("[assign_taxonomy]   already removed above by the tokens or the best-hit check: ",
            length(already_ids), " ASV(s), ", reads_of(already_ids), " reads (keep ",
            length(already_keep), ", discard ", length(already_disc), ", best-hit check ",
            length(already_best), ")")
    message("[assign_taxonomy]   removed only by the organelle filter: ", length(only_metaxa2),
            " ASV(s), ", reads_of(only_metaxa2), " reads")
    for (id in only_metaxa2) {
      reported <- asv_table_all[id, "taxonomy"]
      if (is.na(reported) || reported == "") reported <- "unclassified"
      message("[assign_taxonomy]     ", id, ": ", metaxa2_origin[[id]], " by Metaxa2, ",
              reads_of(id), " reads (reported as ", reported, ")")
    }
  }
}
message("[assign_taxonomy] Final ASV count: ", nrow(asv_table))

# Fail loud rather than write a header-only table: a filter that removes 100% of
# ASVs is almost always a marker/config mismatch (e.g. 16S keep tokens on an ITS run).
if (isTRUE(filter_enabled) && n_before_filter > 0L && nrow(asv_table) == 0L) {
  # If the organelle filter took part, say how much, because then the lists may
  # well be right and the sample simply holds little but organelle 16S.
  n_organelle_only <- sum(filter_status == "removed_organelle")
  organelle_note <- if (n_organelle_only > 0L) {
    paste0(" ", n_organelle_only, " of them were removed only because Metaxa2 called them ",
           "mitochondria or chloroplast; amplicon.taxonomy.filter.organelles: off keeps those.")
  } else ""
  stop("[assign_taxonomy] the keep/discard filter removed ALL ", n_before_filter,
       " ASVs (amp_type=", amp_type, "). The keep/discard lists likely do not match this ",
       "marker: for ITS use keep: [k__Fungi]; for 16S use keep: [k__Bacteria, k__Archaea]. ",
       "Fix amplicon.taxonomy.filter, or set amplicon.taxonomy.filter.enabled: false.",
       organelle_note)
}

# ── Output (a): asv_table.txt — rows=ASV_IDs ─────────────────────────────────
# The main table: one row per ASV, one column per sample (read counts) plus a
# taxonomy column. Its sample columns are summed again downstream by rule
# aggregate_read_counts to build the final stage of read_tracking.txt.
dir.create(dirname(snakemake@output[["asv_table"]]), recursive = TRUE, showWarnings = FALSE)
write.table(asv_table, file = snakemake@output[["asv_table"]], col.names = NA, sep = "\t")

# ── Output (b): asv_table_seqs.txt — rows=sequences ──────────────────────────
# Identical to (a), except each row is named by the ASV's actual sequence
# instead of its short ASV_N id — for matching these ASVs directly against
# sequences from another tool or database, without going through the id.
asv_table_seqs <- asv_table
rownames(asv_table_seqs) <- seqs_sub[rownames(asv_table)]
write.table(asv_table_seqs, file = snakemake@output[["asv_table_seqs"]], col.names = NA, sep = "\t")

# ── Output (c): taxon_seq_table.txt — rows=ASV_IDs, parsed cols + sequence ───
# Same ASVs and taxonomy as (a), but with the taxonomy split into one column
# per rank (Kingdom, Phylum, ...) instead of one semicolon-joined string, plus
# the ASV's sequence — convenient for filtering or grouping by a single rank
# without string-parsing the taxonomy column first. Has no per-sample counts.
taxon_table <- taxa_df[rownames(asv_table), , drop = FALSE]
# Rank cells hold the bare name, without the rank prefix ("Fungi", "Ascomycota", ...,
# "Hanseniaspora_uvarum"): the column header already names the rank. That is what
# the sintax path writes for every marker, and what this path writes for 16S and
# 18S, whose references carry no prefixes. The references of the embedded-style
# markers (ITS, gyrB, rpoB) do carry them (k__Fungi, p__Bacillota), so here the
# prefix is removed, rank by rank. Only this table changes: the taxonomy strings in
# the other three tables keep their prefixes, as for every marker (they are what the
# keep/discard filter matched). gyrB's first rank, the gene tag (gyrB or other), has
# no prefix and stays as it is.
if (identical(prefix_style, "embedded")) {
  for (i in seq_along(tax_levels)) {
    rank_name <- tax_levels[i]
    prefix    <- rank_prefixes[i]
    if (prefix == "") next
    cells <- taxon_table[[rank_name]]
    # Only a cell that starts with this rank's prefix changes; NA (a rank the
    # classifier did not reach) stays NA.
    has_prefix <- !is.na(cells) & startsWith(cells, prefix)
    cells[has_prefix] <- substring(cells[has_prefix], nchar(prefix) + 1L)
    taxon_table[[rank_name]] <- cells
  }
}
taxon_table$sequence <- seqs_sub[rownames(asv_table)]
write.table(taxon_table, file = snakemake@output[["taxon_seq_table"]], col.names = NA, sep = "\t")

# ── Output (d): asv_table_unfiltered.txt — every ASV that reached the filter ──
# The layout of (a), but for EVERY ASV that entered the filter (organelles and the
# other removed ASVs included), in the same order, plus two columns:
#   metaxa2_origin — Metaxa2's call (mitochondria, chloroplast, bacteria, archaea,
#                    eukaryota, uncertain), NA when Metaxa2 made none or did not run
#   filter_status  — kept, or the first step that removed the ASV (see above)
# So the "kept" rows are exactly asv_table.txt. Meant for analyses that need the
# removed ASVs back, e.g. how many organelle reads each primer set amplifies;
# nothing in the pipeline reads it.
asv_table_unfiltered <- asv_table_all
asv_table_unfiltered$metaxa2_origin <- unname(metaxa2_origin[rownames(asv_table_all)])
asv_table_unfiltered$filter_status  <- unname(filter_status[rownames(asv_table_all)])
write.table(asv_table_unfiltered, file = snakemake@output[["asv_table_unfiltered"]],
            col.names = NA, sep = "\t")

message("[assign_taxonomy] DONE. Outputs written to ", dirname(snakemake@output[["asv_table"]]))

sink(type = "message")
sink(type = "output")
close(log_con)
