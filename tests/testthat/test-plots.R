test_that("plot helpers return ggplot objects on matched mock data", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_if_not_installed("ggplot2")

  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 6),
    ranges = IRanges::IRanges(start = seq(100, 600, by = 100), width = 1),
    strand = rep("+", 6)
  )
  names(gr) <- paste0("s", seq_along(gr))
  S4Vectors::mcols(gr)$label <- c(1, 1, 1, 0, 0, 0)
  S4Vectors::mcols(gr)$location <- c("coding", "coding", "threeUTR", "coding", "coding", "threeUTR")
  S4Vectors::mcols(gr)$metagene_split3 <- c(1.10, 1.55, 2.20, 1.11, 1.60, 2.18)
  S4Vectors::mcols(gr)$nearest_exon_junction_dist <- c(5, 20, 10, 6, 18, 11)
  S4Vectors::mcols(gr)$start_dist <- c(30, 80, 180, 32, 79, 182)
  S4Vectors::mcols(gr)$stop_dist <- c(200, 140, 20, 202, 138, 19)
  S4Vectors::mcols(gr)$start_dist_tx <- c(30, 80, 180, 32, 79, 182)
  S4Vectors::mcols(gr)$stop_dist_tx <- c(200, 140, 20, 202, 138, 19)
  S4Vectors::mcols(gr)$feature_width <- c(100, 100, 80, 100, 90, 80)
  S4Vectors::mcols(gr)$kmer <- c("GGACT", "AGACT", "AAACA", "GGACT", "AGACT", "AAACA")

  gr <- match_background(gr, bin_match = FALSE)

  expect_s3_class(plot_metagene_density(gr), "gg")
  expect_s3_class(plot_junction_distance_density(gr), "gg")
  expect_s3_class(plot_start_stop_distance_density(gr), "gg")
  expect_s3_class(plot_feature_width_bins(gr), "gg")
  expect_s3_class(plot_kmer_counts(gr), "gg")
})
