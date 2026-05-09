test_that("motif_enrichment_profile and plot_motif_enrichment work on a tiny genome", {
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

  prof <- motif_enrichment_profile(
    gr,
    genome = genome,
    motif = "AAA",
    window = 10,
    hit_position = "center",
    smooth_window = 1
  )

  expect_true(all(c("relative_position", "group", "motif", "count", "hits_per_site") %in% names(prof)))
  expect_true(any(prof$count[prof$group == "positive"] > 0))

  p <- plot_motif_enrichment(
    gr,
    genome = genome,
    motif = "AAA",
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

  prof <- motif_enrichment_profile(
    gr,
    genome = genome,
    motif = "A",
    window = 10,
    smooth_window = 1,
    drop_edge_positions = TRUE
  )

  expect_false(any(prof$relative_position %in% c(-10, 10)))

  prof2 <- motif_enrichment_profile(
    gr,
    genome = genome,
    motif = "A",
    window = 10,
    smooth_window = 1,
    edge_trim = 3
  )

  expect_false(any(prof2$relative_position %in% c(-10, -9, -8, 8, 9, 10)))
  expect_true(all(c(-7, 0, 7) %in% prof2$relative_position))

  prof3 <- motif_enrichment_profile(
    gr,
    genome = genome,
    motif = "A",
    window = 10,
    smooth_window = 1,
    center_exclude = 2
  )

  expect_false(any(prof3$relative_position %in% -2:2))
  expect_true(all(c(-3, 3) %in% prof3$relative_position))

})

test_that("motif_enrichment_profile can use spliced transcript windows and RNA/IUPAC motifs", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  seq <- paste(rep("A", 30), collapse = "")
  substr(seq, 1, 5) <- "CCCCC"
  substr(seq, 21, 25) <- "TTTTT"
  genome <- Biostrings::DNAStringSet(c(chr1 = seq))

  ex <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1, 21), width = 5),
    strand = c("+", "+")
  )
  S4Vectors::mcols(ex)$exon_rank <- c(1L, 2L)
  resources <- list(exons_by_tx = GenomicRanges::GRangesList(tx1 = ex))

  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 21, width = 1),
    strand = "+"
  )
  S4Vectors::mcols(gr)$match_set <- "positive"
  S4Vectors::mcols(gr)$tx_name <- "tx1"
  S4Vectors::mcols(gr)$tx_pos <- 6L

  prof_tx <- motif_enrichment_profile(
    gr,
    genome = genome,
    motif = "CCTTT",
    window = 4,
    resources = resources,
    window_mode = "transcript",
    smooth_window = 1
  )
  expect_gt(sum(prof_tx$count), 0L)
  expect_equal(unique(prof_tx$window_mode), "transcript")

  prof_u <- motif_enrichment_profile(
    gr,
    genome = genome,
    motif = "UUUUU",
    window = 4,
    resources = resources,
    window_mode = "transcript",
    smooth_window = 1
  )
  expect_gt(sum(prof_u$count), 0L)

  prof_genomic <- motif_enrichment_profile(
    gr,
    genome = genome,
    motif = "CCTTT",
    window = 4,
    window_mode = "genomic",
    smooth_window = 1
  )
  expect_equal(sum(prof_genomic$count), 0)

  expect_error(
    motif_enrichment_profile(
      gr,
      genome = genome,
      motif = "UZX",
      window = 4,
      resources = resources,
      window_mode = "transcript"
    ),
    "non-IUPAC"
  )
})
