test_that("prepare_sites creates width-1 sites and stable IDs", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(10, 20), width = c(3, 1)),
    strand = c("+", "-")
  )

  expect_warning(out <- prepare_sites(gr), "width != 1")
  expect_true(all(GenomicRanges::width(out) == 1L))
  expect_true("site_id" %in% colnames(S4Vectors::mcols(out)))
  expect_false(anyDuplicated(names(out)) > 0L)
})

test_that("match_background returns reciprocal matched pairs", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 8),
    ranges = IRanges::IRanges(start = seq(100, 800, by = 100), width = 1),
    strand = rep("+", 8)
  )
  names(gr) <- paste0("s", seq_along(gr))

  S4Vectors::mcols(gr)$label <- c(1, 1, 1, 0, 0, 0, 0, 0)
  S4Vectors::mcols(gr)$location <- c("coding", "coding", "threeUTR", "coding", "coding", "threeUTR", "threeUTR", "fiveUTR")
  S4Vectors::mcols(gr)$metagene_split3 <- c(1.10, 1.55, 2.20, 1.11, 1.60, 2.18, 2.90, 0.50)
  S4Vectors::mcols(gr)$feature_width <- c(100, 100, 80, 100, 90, 80, 70, 30)
  S4Vectors::mcols(gr)$nearest_exon_junction_dist <- c(5, 20, 10, 6, 18, 11, 100, 4)
  S4Vectors::mcols(gr)$start_dist_tx <- c(30, 80, 180, 32, 79, 182, 210, 15)
  S4Vectors::mcols(gr)$stop_dist_tx <- c(200, 140, 20, 202, 138, 19, 10, 300)

  out <- match_background(gr, bin_match = FALSE, return_diagnostics = TRUE)

  expect_equal(sum(S4Vectors::mcols(out)$is_matched_negative), 3L)
  expect_equal(sum(!is.na(S4Vectors::mcols(out)$matched_negative_id[S4Vectors::mcols(out)$label == 1L])), 3L)

  paired <- subset_matched_pairs(out)
  expect_equal(length(paired), 6L)
  expect_true(all(S4Vectors::mcols(paired)$match_set %in% c("positive", "matched_negative")))
})

test_that("as_site_table returns scalar data frame", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10, 10), strand = "+")
  gr <- prepare_sites(gr)
  S4Vectors::mcols(gr)$location <- "coding"
  S4Vectors::mcols(gr)$feature <- "CDS"

  tab <- as_site_table(gr)
  expect_s3_class(tab, "data.frame")
  expect_true("primary_region_class" %in% colnames(tab))
  expect_equal(tab$primary_region_class[1], "CDS")
})
