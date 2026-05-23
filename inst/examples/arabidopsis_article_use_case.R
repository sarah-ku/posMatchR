# Arabidopsis single-nucleotide site annotation, matching, and motif-enrichment workflow
#
# This script is intended as a clean use-case example for manuscript figures.
# It assumes a stranded single-nucleotide GRanges object, such as SAC-seq,
# nanopore m6A, GLORI-like sites, or HyperTRIBE/RNA-editing sites.

library(posMatchR)
library(GenomicRanges)
library(GenomeInfoDb)
library(S4Vectors)
library(TxDb.Athaliana.BioMart.plantsmart51)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(org.At.tair.db)
library(ggplot2)

# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------

sites_path <- "/path/to/arabidopsis_sites_GRanges.rds"  # or .Rdat/.RData containing one GRanges
outdir <- "posMatchR_arabidopsis_article_use_case"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(outdir, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

arab_chrs <- c("1", "2", "3", "4", "5")

# Keep the TxDb chromosome naming style throughout this example. The TAIR9
# BSgenome object normally uses Chr1-Chr5/ChrM/ChrC, so we rename this local
# genome object to match the TxDb and site GRanges.
txdb <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
genome <- BSgenome.Athaliana.TAIR.TAIR9::Athaliana
GenomeInfoDb::seqlevels(genome) <- c("1", "2", "3", "4", "5", "Mt", "Pt")

save_plot <- function(plot, file, width = 7, height = 4) {
  ggplot2::ggsave(file.path(figdir, file), plot, width = width, height = height)
}

# -------------------------------------------------------------------------
# 1. Load and normalise input sites
# -------------------------------------------------------------------------

sites <- load_sites(sites_path)

# Convert Chr1-Chr5 to 1-5 if required, then keep only the nuclear chromosomes
# represented by the TxDb. This deliberately drops Pt/Mt/non-standard contigs.
GenomeInfoDb::seqlevels(sites) <- sub(
  "^Chr([1-5])$",
  "\\1",
  GenomeInfoDb::seqlevels(sites)
)

sites <- GenomeInfoDb::keepSeqlevels(sites, arab_chrs, pruning.mode = "coarse")
sites <- sites[as.character(GenomicRanges::strand(sites)) %in% c("+", "-")]
sites <- prepare_sites(sites, site_id_col = "site_id", strip_mcols = FALSE)

# -------------------------------------------------------------------------
# 2. Annotate sites and add site-centred 5-mers
# -------------------------------------------------------------------------

ann <- annotate_sites(
  gr = sites,
  txdb = txdb,
  chrs = arab_chrs,
  tx_select = "longest",
  orgdb = org.At.tair.db::org.At.tair.db,
  gene_keytype = "TAIR",
  gene_symbol_col = "SYMBOL",
  gene_name_col = "GENENAME"
)

ann <- add_kmer(
  ann,
  genome = genome,
  k = 5L,
  chrs = arab_chrs
)

saveRDS(ann, file.path(outdir, "arabidopsis_annotated_sites.rds"))
write.csv(as_basic_site_table(ann), file.path(outdir, "arabidopsis_annotated_sites_basic.csv"), row.names = FALSE)
write.csv(as_site_table(ann, columns = "all"), file.path(outdir, "arabidopsis_annotated_sites_full.csv"), row.names = FALSE)

message("Annotated site classes:")
print(table(ann$location, useNA = "ifany"))
message("Top 5-mers:")
print(head(sort(table(ann$kmer), decreasing = TRUE), 20))

save_plot(plot_metagene_density(ann), "01_sites_metagene_density.pdf")
save_plot(plot_junction_distance_density(ann), "02_sites_exon_junction_distance_density.pdf")
save_plot(plot_kmer_counts(ann, top_n = 40), "03_sites_5mer_counts.pdf", width = 8, height = 5)

# -------------------------------------------------------------------------
# 3. Build a candidate background universe
# -------------------------------------------------------------------------

resources <- build_tx_resources(txdb)

# This scans the same transcript space represented in the foreground and finds
# all instances of foreground-observed 5-mers occurring at least min_count times.
# Save the universe because this is the slowest step in the example workflow.
bg_cache <- file.path(outdir, "arabidopsis_background_5mer_universe.rds")
if (file.exists(bg_cache)) {
  bg_5mer <- readRDS(bg_cache)
} else {
  bg_5mer <- make_kmer_universe(
    foreground = ann,
    txdb = txdb,
    genome = genome,
    kmer_col = "kmer",
    min_count = 10,
    scope = "transcripts",
    resources = resources
  )
  saveRDS(bg_5mer, bg_cache)
}

message("Candidate background universe:")
print(bg_5mer)
print(head(sort(table(bg_5mer$candidate_kmer), decreasing = TRUE), 20))

# -------------------------------------------------------------------------
# 4. Combine positives and background candidates, then annotate together
# -------------------------------------------------------------------------

combined <- combine_site_sets(
  positives = ann,
  background = bg_5mer,
  label_col = "label",
  site_id_col = "site_id"
)

combined <- annotate_sites(
  gr = combined,
  txdb = txdb,
  chrs = arab_chrs,
  tx_select = "longest",
  orgdb = org.At.tair.db::org.At.tair.db,
  gene_keytype = "TAIR",
  gene_symbol_col = "SYMBOL",
  gene_name_col = "GENENAME"
)

combined <- add_kmer(
  combined,
  genome = genome,
  k = 5L,
  chrs = arab_chrs
)

saveRDS(combined, file.path(outdir, "arabidopsis_combined_foreground_background.rds"))

# -------------------------------------------------------------------------
# 5. Four matching strategies
# -------------------------------------------------------------------------

match_locations <- c("fiveUTR", "coding", "threeUTR")

# Strategy A: random same-gene and same-region background.
random_region <- match_random_background(
  combined,
  group_col = "gene_id",
  match_location = TRUE,
  location_col = "location",
  locations = match_locations,
  kmer_match = FALSE,
  seed = 1L
)
random_region_sets <- subset_matched_sets(random_region)

# Strategy B: random same-gene, same-region, same-5-mer background.
random_region_kmer <- match_random_background(
  combined,
  group_col = "gene_id",
  match_location = TRUE,
  location_col = "location",
  locations = match_locations,
  kmer_match = TRUE,
  kmer_col = "kmer",
  seed = 1L
)
random_region_kmer_sets <- subset_matched_sets(random_region_kmer)

# Strategy C: covariate-aware matching using region and metagene geometry.
covariate <- match_background(
  combined,
  label_col = "label",
  kmer_match = FALSE,
  kmer_col = "kmer",
  meta_col = "metagene_split3",
  meta_tol = 0.03,
  enforce_meta_tol = TRUE,
  bin_match = TRUE,
  seed = 1L
)
covariate_sets <- subset_matched_sets(covariate)

# Strategy D: covariate-aware matching plus exact 5-mer balancing.
covariate_kmer <- match_background(
  combined,
  label_col = "label",
  kmer_match = TRUE,
  kmer_col = "kmer",
  meta_col = "metagene_split3",
  meta_tol = 0.03,
  enforce_meta_tol = TRUE,
  bin_match = TRUE,
  seed = 1L
)
covariate_kmer_sets <- subset_matched_sets(covariate_kmer)

matched_objects <- list(
  random_region = random_region_sets,
  random_region_kmer = random_region_kmer_sets,
  covariate = covariate_sets,
  covariate_kmer = covariate_kmer_sets
)

saveRDS(matched_objects, file.path(outdir, "arabidopsis_matched_sets_four_strategies.rds"))


print(list(
  random_region = S4Vectors::metadata(random_region)$match_diagnostics,
  random_region_kmer = S4Vectors::metadata(random_region_kmer)$match_diagnostics,
  covariate = S4Vectors::metadata(covariate)$match_diagnostics,
  covariate_kmer = S4Vectors::metadata(covariate_kmer)$match_diagnostics
))

# -------------------------------------------------------------------------
# 6. Diagnostic plots for each matched set
# -------------------------------------------------------------------------

for (nm in names(matched_objects)) {
  gr <- matched_objects[[nm]]
  message(nm, ":")
  print(table(gr$match_set, useNA = "ifany"))

  save_plot(
    plot_metagene_density(gr, set_col = "match_set"),
    paste0("04_", nm, "_metagene_density.pdf")
  )

  save_plot(
    plot_junction_distance_density(gr, set_col = "match_set"),
    paste0("05_", nm, "_exon_junction_distance_density.pdf")
  )

  save_plot(
    plot_kmer_counts(gr, set_col = "match_set", top_n = 40),
    paste0("06_", nm, "_5mer_counts.pdf"),
    width = 8,
    height = 5
  )
}

# -------------------------------------------------------------------------
# 7. Motif-enrichment profiles around matched positives and negatives
# -------------------------------------------------------------------------

# RNA motif UNUNU is scanned as DNA TNTNT after strand-oriented extraction.
# scan_rc=FALSE is intentional: the window sequence is already oriented as RNA.
motif_to_show <- "UNUNU"

for (nm in names(matched_objects)) {
  gr <- matched_objects[[nm]]
  save_plot(
    plot_motif_enrichment(
      gr,
      genome = genome,
      motif = motif_to_show,
      motif_name = "UNUNU / TNTNT",
      window = 250,
      set_col = "match_set",
      hit_position = "center",
      bin_size = 1,
      smooth_window = 11,
      scan_rc = FALSE,
      chrs = arab_chrs,
      drop_edge_positions = TRUE
    ),
    paste0("07_", nm, "_UNUNU_motif_enrichment_250bp.pdf"),
    width = 7,
    height = 4
  )
}

message("Finished. Outputs written to: ", normalizePath(outdir))
