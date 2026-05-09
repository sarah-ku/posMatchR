test_that("site_enrichment_profile counts same-strand point features in RNA orientation", {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(c(100, 200), width = 1),
    strand = c("+", "-")
  )
  S4Vectors::mcols(gr)$match_set <- c("positive", "matched_negative")

  feat <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(c(110, 210, 190, 120), width = 1),
    strand = c("+", "-", "-", "-")
  )

  prof <- site_enrichment_profile(
    gr,
    query = feat,
    query_name = "iCLIP",
    window = 20,
    set_col = "match_set",
    require_same_strand = TRUE,
    smooth_window = 1,
    drop_edge_positions = FALSE,
    normalise_per = 1000
  )

  pos10 <- prof[prof$group == "positive" & prof$relative_position == 10 & prof$feature == "iCLIP", ]
  neg_m10 <- prof[prof$group == "negative" & prof$relative_position == -10 & prof$feature == "iCLIP", ]
  neg_p10 <- prof[prof$group == "negative" & prof$relative_position == 10 & prof$feature == "iCLIP", ]

  expect_equal(pos10$count, 1L)
  expect_equal(neg_m10$count, 1L)
  expect_equal(neg_p10$count, 1L)
  expect_equal(pos10$hits_per_n_sites, 1000)

  # The opposite-strand feature at 120 should not be counted for the + focal site.
  pos20 <- prof[prof$group == "positive" & prof$relative_position == 20 & prof$feature == "iCLIP", ]
  expect_equal(pos20$count, 0L)
})

test_that("plot_site_enrichment returns a ggplot for point features", {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(c(100, 200), width = 1),
    strand = c("+", "-")
  )
  S4Vectors::mcols(gr)$match_set <- c("positive", "matched_negative")
  feat <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(c(110, 210), width = 1),
    strand = c("+", "-")
  )

  p <- plot_site_enrichment(
    gr,
    query = feat,
    query_name = "iCLIP",
    window = 20,
    set_col = "match_set",
    smooth_window = 1,
    drop_edge_positions = FALSE
  )
  testthat::expect_s3_class(p, "gg")
})
