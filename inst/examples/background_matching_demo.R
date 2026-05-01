# Background matching demonstration after annotating foreground sites.
#
# Assumes you already have:
#   ann    annotated foreground GRanges with add_kmer() already run
#   txdb   TxDb with seqlevels compatible with genome
#   genome BSgenome object
#
# This script creates two candidate universes and compares two matchers.

library(posMatchR)

# Candidate universe 1: all transcript-oriented adenosines in the same foreground genes.
bg_A <- make_base_universe(
  foreground = ann,
  txdb = txdb,
  genome = genome,
  base = "A",
  scope = "genes"
)

# Candidate universe 2: all observed foreground 5-mers in the same foreground genes.
# Use min_count to avoid scanning very rare kmers if desired.
bg_5mer <- make_kmer_universe(
  foreground = ann,
  txdb = txdb,
  genome = genome,
  kmer_col = "kmer",
  min_count = 1,
  scope = "genes"
)

# IUPAC motif candidate universe. For DRACH, the modified A is the third position.
bg_drach <- make_motif_universe(
  foreground = ann,
  txdb = txdb,
  genome = genome,
  patterns = "DRACH",
  site_offset = 3,
  scope = "genes"
)

# Use one candidate universe. For motif-balanced m6A-style tests, bg_5mer is usually preferable.
bg <- bg_5mer

combined <- combine_site_sets(
  positives = ann,
  background = bg,
  label_col = "label",
  site_id_col = "site_id"
)
combined <- annotate_sites(combined, txdb = txdb, tx_select = "longest")
combined <- add_kmer(combined, genome = genome, k = 5)

# Baseline: random within gene, with exact k-mer matching.
random_gene_kmer <- match_random_background(
  combined,
  group_col = "gene_id",
  kmer_match = TRUE,
  seed = 1L
)

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

plot_metagene_density(random_gene_kmer, set_col = "match_set")
plot_junction_distance_density(random_gene_kmer, set_col = "match_set")
plot_metagene_density(matched_covariate, set_col = "match_set")
plot_junction_distance_density(matched_covariate, set_col = "match_set")

S4Vectors::metadata(random_gene_kmer)$match_diagnostics
S4Vectors::metadata(matched_covariate)$match_diagnostics
