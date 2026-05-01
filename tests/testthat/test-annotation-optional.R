test_that("human TxDb annotation smoke test runs when optional packages are installed", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(100000, 100000),
    strand = "+"
  )
  ann <- annotate_sites(gr, txdb = txdb, seqstyle = "UCSC", chrs = "chr1")
  expect_true(all(c("location", "feature", "metagene_prop", "nearest_exon_junction_dist") %in% colnames(S4Vectors::mcols(ann))))
})

test_that("Arabidopsis chromosome filtering smoke test runs when optional package is installed", {
  skip_if_not_installed("TxDb.Athaliana.BioMart.plantsmart51")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  txdb <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
  gr <- GenomicRanges::GRanges(
    seqnames = c("1", "Mt"),
    ranges = IRanges::IRanges(c(1000, 1000), c(1000, 1000)),
    strand = c("+", "+")
  )
  expect_warning(
    ann <- annotate_sites(gr, txdb = txdb, chrs = c("1", "2", "3", "4", "5")),
    "Dropping site seqlevels"
  )
  expect_true(length(ann) <= length(gr))
  expect_true(all(as.character(GenomicRanges::seqnames(ann)) %in% c("1", "2", "3", "4", "5")))
})
