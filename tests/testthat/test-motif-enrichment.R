test_that("site_enrichment_profile and plot_site_enrichment work on a tiny motif example", {
  skip_if_not_installed("Biostrings")

  seq <- paste(rep("C", 101), collapse = "")
  substr(seq, 51, 53) <- "AAA"
  genome <- Biostrings::DNAStringSet(c(chr1 = seq))

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(52, 80), width = 1),
    strand = c("+", "+")
  )
  S4Vectors::mcols(gr)$match_set <- c("positive", "matched_negative")

  prof <- site_enrichment_profile(
    gr,
    genome = genome,
    query = "AAA",
    window = 10,
    hit_position = "center",
    smooth_window = 1
  )

  expect_true(all(c("relative_position", "group", "feature", "count", "hits_per_site") %in% names(prof)))
  expect_true(any(prof$count[prof$group == "positive"] > 0))

  p <- plot_site_enrichment(
    gr,
    genome = genome,
    query = "AAA",
    window = 10,
    hit_position = "center",
    smooth_window = 1
  )
  expect_s3_class(p, "gg")
})

test_that("motif_enrichment_profile can omit edge positions", {
  skip_if_not_installed("Biostrings")

  seq <- paste(rep("A", 101), collapse = "")
  genome <- Biostrings::DNAStringSet(c(chr1 = seq))
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 51, width = 1),
    strand = "+"
  )
  S4Vectors::mcols(gr)$match_set <- "positive"

  prof <- site_enrichment_profile(
    gr,
    genome = genome,
    query = "A",
    window = 10,
    smooth_window = 1,
    drop_edge_positions = TRUE
  )

  expect_false(any(prof$relative_position %in% c(-10, 10)))

  prof2 <- site_enrichment_profile(
    gr,
    genome = genome,
    query = "A",
    window = 10,
    smooth_window = 1,
    edge_trim = 3
  )

  expect_false(any(prof2$relative_position %in% c(-10, -9, -8, 8, 9, 10)))
  expect_true(all(c(-7, 0, 7) %in% prof2$relative_position))

  prof3 <- site_enrichment_profile(
    gr,
    genome = genome,
    query = "A",
    window = 10,
    smooth_window = 1,
    center_exclude = 2
  )

  expect_false(any(prof3$relative_position %in% -2:2))
  expect_true(all(c(-3, 3) %in% prof3$relative_position))

})
