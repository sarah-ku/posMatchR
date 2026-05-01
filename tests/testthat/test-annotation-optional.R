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

test_that("Arabidopsis add_kmer handles TxDb/genome seqlevel name differences", {
  skip_if_not_installed("TxDb.Athaliana.BioMart.plantsmart51")
  skip_if_not_installed("BSgenome.Athaliana.TAIR.TAIR9")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  genome <- BSgenome.Athaliana.TAIR.TAIR9::Athaliana
  txdb <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
  chrs <- paste0("Chr", 1:5)
  txdb <- standardize_seqlevels(txdb, target = genome, chrs = chrs, keep = TRUE)

  gr <- GenomicRanges::GRanges(
    seqnames = "Chr1",
    ranges = IRanges::IRanges(4506, 4506),
    strand = "+"
  )

  ann <- annotate_sites(gr, txdb = txdb, chrs = chrs)
  ann <- add_kmer(ann, genome = genome, k = 5L, chrs = chrs)
  expect_true("kmer" %in% colnames(S4Vectors::mcols(ann)))
  expect_equal(length(S4Vectors::mcols(ann)$kmer), length(ann))
})
