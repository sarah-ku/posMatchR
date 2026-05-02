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

test_that("basic site table is compact", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10, 10), strand = "+")
  gr <- prepare_sites(gr)
  S4Vectors::mcols(gr)$location <- "coding"
  S4Vectors::mcols(gr)$feature <- "CDS"
  S4Vectors::mcols(gr)$tx_len <- 1000L
  S4Vectors::mcols(gr)$nearest_exon_junction_dist <- 12L

  tab_all <- as_site_table(gr, columns = "all")
  tab_basic <- as_basic_site_table(gr)
  expect_true("tx_len" %in% colnames(tab_all))
  expect_false("tx_len" %in% colnames(tab_basic))
  expect_true("nearest_exon_junction_dist" %in% colnames(tab_basic))
})

test_that("match_random_background can optionally preserve location", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 6),
    ranges = IRanges::IRanges(start = seq(10, 60, by = 10), width = 1),
    strand = rep("+", 6)
  )
  names(gr) <- paste0("r", seq_along(gr))
  S4Vectors::mcols(gr)$label <- c(1, 1, 0, 0, 0, 0)
  S4Vectors::mcols(gr)$gene_id <- rep("g1", 6)
  S4Vectors::mcols(gr)$location <- c("coding", "threeUTR", "coding", "coding", "threeUTR", "threeUTR")

  out <- match_random_background(
    gr,
    group_col = "gene_id",
    match_location = TRUE,
    seed = 1L
  )
  paired <- subset_matched_sets(out)

  mc <- S4Vectors::mcols(out)
  pos_idx <- which(mc$label == 1L & !is.na(mc$matched_negative_id))
  neg_idx <- match(as.character(mc$matched_negative_id[pos_idx]), names(out))
  expect_true(all(as.character(mc$location[pos_idx]) == as.character(mc$location[neg_idx])))
  expect_equal(length(paired), 4L)
})

test_that("match_random_background gives exact kmer balance for reciprocal pairs", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 10),
    ranges = IRanges::IRanges(start = seq(10, 100, by = 10), width = 1),
    strand = rep("+", 10)
  )
  names(gr) <- paste0("k", seq_along(gr))
  S4Vectors::mcols(gr)$label <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  S4Vectors::mcols(gr)$gene_id <- rep("g1", 10)
  S4Vectors::mcols(gr)$location <- rep("coding", 10)
  S4Vectors::mcols(gr)$kmer <- c("AAAAA", "AAAAA", "CCCCC", "GGGGG", "AAAAA", "AAAAA", "CCCCC", "CCCCC", "TTTTT", "GGGGG")

  out <- match_random_background(
    gr,
    group_col = "gene_id",
    match_location = TRUE,
    kmer_match = TRUE,
    kmer_col = "kmer",
    seed = 1L
  )
  paired <- subset_matched_sets(out)
  bal <- summarise_matched_kmer_balance(paired, subset_first = FALSE)

  expect_equal(sum(S4Vectors::mcols(out)$is_matched_positive), 4L)
  expect_equal(sum(S4Vectors::mcols(out)$is_matched_negative), 4L)
  expect_true(all(bal$difference == 0L))
  expect_equal(S4Vectors::metadata(out)$match_diagnostics$kmer_mismatch_in_validated_pairs, 0L)
})
