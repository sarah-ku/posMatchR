#!/usr/bin/env Rscript
# Local smoke test for human single-nucleotide GRanges, e.g. GLORI/m6A or A-to-I editing sites.
# Usage:
#   Rscript inst/examples/local_human_smoke_test.R /path/to/sites.rds
#   Rscript inst/examples/local_human_smoke_test.R /path/to/sites.Rdat output_dir

suppressPackageStartupMessages({
  library(posMatchR)
  library(GenomicRanges)
  library(S4Vectors)
  library(GenomeInfoDb)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
})

helper <- file.path("inst", "examples", "smoke_test_helpers.R")
if (!file.exists(helper)) helper <- system.file("examples", "smoke_test_helpers.R", package = "posMatchR", mustWork = TRUE)
source(helper)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript inst/examples/local_human_smoke_test.R /path/to/sites.rds [output_dir]")
}
input_path <- args[[1]]
outdir <- if (length(args) >= 2L && nzchar(args[[2]])) args[[2]] else "posmatchr_human_smoke_results"

sites <- posMatchR::load_sites(input_path)
GenomeInfoDb::seqlevelsStyle(sites) <- "UCSC"
canonical_chrs <- paste0("chr", c(1:22, "X", "Y"))
sites <- GenomeInfoDb::keepSeqlevels(sites, canonical_chrs, pruning.mode = "coarse")
if (length(sites) > 1000L) {
  message("Using the first 1000 sites for this smoke test. Edit the script to run all sites.")
  sites <- sites[seq_len(1000L)]
}

summarise_granges_for_smoke(sites, "raw human input")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

ann <- posMatchR::annotate_sites(
  gr = sites,
  txdb = txdb,
  seqstyle = "UCSC",
  chrs = canonical_chrs,
  tx_select = "longest",
  orgdb = org.Hs.eg.db::org.Hs.eg.db,
  gene_keytype = "ENTREZID",
  gene_symbol_col = "SYMBOL",
  gene_name_col = "GENENAME",
  quiet = FALSE
)

summarise_granges_for_smoke(ann, "annotated human output")
cat("\nLocation counts:\n")
print(table(S4Vectors::mcols(ann)$location, useNA = "ifany"))
cat("\nExample annotated rows:\n")
print(utils::head(posMatchR::as_site_table(ann), 10))

save_posmatchr_smoke_outputs(ann, outdir, "human")
