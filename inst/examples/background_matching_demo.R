# Background matching demonstration after annotating foreground sites.
#
# Assumes you already have:
#   ann    annotated foreground GRanges with add_kmer() already run
#   txdb   TxDb with seqlevels compatible with genome
#   genome BSgenome object
#
# This script creates candidate universes and compares random and covariate-aware
# matchers. Diagnostic plots are made on the matched pairs only.

library(posMatchR)

resources <- build_tx_resources(txdb)

# Candidate universe 1: all transcript-oriented adenosines in the same foreground genes.
bg_A <- make_base_universe(
  foreground = ann,
  txdb = txdb,
  genome = genome,
  base = "A",
  scope = "genes",
  resources = resources
)

# Candidate universe 2: all observed foreground 5-mers in the same foreground genes.
# Use min_count to avoid scanning very rare kmers if desired.
bg_5mer <- make_kmer_universe(
  foreground = ann,
  txdb = txdb,
  genome = genome,
  kmer_col = "kmer",
  min_count = 1,
  scope = "genes",
  resources = resources
)

# IUPAC motif candidate universe. For DRACH, the modified A is the third position.
bg_drach <- make_motif_universe(
  foreground = ann,
  txdb = txdb,
  genome = genome,
  patterns = "DRACH",
  site_offset = 3,
  scope = "genes",
  resources = resources
)

# Use one candidate universe. For motif-balanced m6A-style tests, bg_5mer is usually preferable.
bg <- bg_5mer

combined <- combine_site_sets(
  positives = ann,
  background = bg,
  label_col = "label",
  site_id_col = "site_id"
)
combined <- annotate_sites(combined, txdb = txdb, tx_select = "longest", resources = resources)
combined <- add_kmer(combined, genome = genome, k = 5)

# Baseline 1: random within gene and broad transcript region.
random_gene_region <- match_random_background(
  combined,
  group_col = "gene_id",
  match_location = TRUE,
  locations = c("fiveUTR", "coding", "threeUTR"),
  kmer_match = FALSE,
  seed = 1L
)
random_gene_region_sets <- subset_matched_sets(random_gene_region)

# Baseline 2: random within gene, broad transcript region, and exact k-mer.
random_gene_kmer <- match_random_background(
  combined,
  group_col = "gene_id",
  match_location = TRUE,
  locations = c("fiveUTR", "coding", "threeUTR"),
  kmer_match = TRUE,
  seed = 1L
)
random_gene_kmer_sets <- subset_matched_sets(random_gene_kmer)

# Covariate-aware: same location/k-mer, close metagene position, and similar feature geometry.
matched_covariate <- match_background(
  combined,
  kmer_match = TRUE,
  meta_col = "metagene_split3",
  meta_tol = 0.05,
  enforce_meta_tol = FALSE,
  bin_match = TRUE,
  seed = 1L
)
matched_covariate_sets <- subset_matched_sets(matched_covariate)

# Non-kmer matched version. Use a tolerance on the same coordinate you plan to plot.
matched_covariate_no_kmer <- match_background(
  combined,
  kmer_match = FALSE,
  meta_col = "metagene_split3",
  meta_tol = 0.03,
  enforce_meta_tol = TRUE,
  bin_match = TRUE,
  seed = 1L
)
matched_covariate_no_kmer_sets <- subset_matched_sets(matched_covariate_no_kmer)

plot_metagene_density(random_gene_region_sets, set_col = "match_set")
plot_junction_distance_density(random_gene_region_sets, set_col = "match_set")

plot_metagene_density(random_gene_kmer_sets, set_col = "match_set")
plot_junction_distance_density(random_gene_kmer_sets, set_col = "match_set")
plot_kmer_counts(random_gene_kmer_sets, set_col = "match_set", top_n = 25)

plot_metagene_density(matched_covariate_sets, set_col = "match_set")
plot_junction_distance_density(matched_covariate_sets, set_col = "match_set")
plot_kmer_counts(matched_covariate_sets, set_col = "match_set", top_n = 25)

plot_metagene_density(matched_covariate_no_kmer_sets, set_col = "match_set")
plot_junction_distance_density(matched_covariate_no_kmer_sets, set_col = "match_set")

S4Vectors::metadata(random_gene_region)$match_diagnostics
S4Vectors::metadata(random_gene_kmer)$match_diagnostics
S4Vectors::metadata(matched_covariate)$match_diagnostics
S4Vectors::metadata(matched_covariate_no_kmer)$match_diagnostics
