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
