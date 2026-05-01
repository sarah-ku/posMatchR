#!/usr/bin/env Rscript
# Local smoke test for Arabidopsis single-nucleotide GRanges, e.g. SAC-seq or HyperTRIBE sites.
# Usage:
#   Rscript inst/examples/local_arabidopsis_smoke_test.R /path/to/sacsGR.Rdat sacsGR
#   Rscript inst/examples/local_arabidopsis_smoke_test.R /path/to/sites.rds "" output_dir

suppressPackageStartupMessages({
  library(posMatchR)
  library(GenomicRanges)
  library(S4Vectors)
  library(GenomeInfoDb)
  library(TxDb.Athaliana.BioMart.plantsmart51)
  library(BSgenome.Athaliana.TAIR.TAIR9)
  library(org.At.tair.db)
  library(AnnotationDbi)
  library(ggplot2)
})

helper <- file.path("inst", "examples", "smoke_test_helpers.R")
if (!file.exists(helper)) helper <- system.file("examples", "smoke_test_helpers.R", package = "posMatchR", mustWork = TRUE)
source(helper)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript inst/examples/local_arabidopsis_smoke_test.R /path/to/sites.rds [object_name] [output_dir]")
}
input_path <- args[[1]]
object_name <- if (length(args) >= 2L && nzchar(args[[2]])) args[[2]] else NULL
outdir <- if (length(args) >= 3L && nzchar(args[[3]])) args[[3]] else "posmatchr_arabidopsis_smoke_results"

canonical_chrs <- paste0("Chr", 1:5)
genome <- BSgenome.Athaliana.TAIR.TAIR9::Athaliana

sites <- posMatchR::load_sites(input_path, object = object_name, chrs = canonical_chrs)
if (length(sites) > 1000L) {
  message("Using the first 1000 sites for this smoke test. Edit the script to run all sites.")
  sites <- sites[seq_len(1000L)]
}

summarise_granges_for_smoke(sites, "raw Arabidopsis input")

txdb <- posMatchR::standardize_seqlevels(
  TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51,
  target = genome,
  chrs = canonical_chrs,
  keep = TRUE
)
orgdb <- org.At.tair.db::org.At.tair.db
org_cols <- AnnotationDbi::columns(orgdb)
gene_name_col <- if ("GENENAME" %in% org_cols) "GENENAME" else NULL

ann <- posMatchR::annotate_sites(
  gr = sites,
  txdb = txdb,
  chrs = canonical_chrs,
  tx_select = "longest",
  orgdb = orgdb,
  gene_keytype = "TAIR",
  gene_symbol_col = "SYMBOL",
  gene_name_col = gene_name_col,
  quiet = FALSE
)

ann <- posMatchR::add_kmer(ann, genome = genome, k = 5L, chrs = canonical_chrs)

summarise_granges_for_smoke(ann, "annotated Arabidopsis output")
cat("\nLocation counts:\n")
print(table(S4Vectors::mcols(ann)$location, useNA = "ifany"))
cat("\nNon-canonical sites are expected to be dropped when chrs = paste0('Chr', 1:5).\n")
cat("\nExample annotated rows:\n")
print(utils::head(posMatchR::as_basic_site_table(ann), 10))

save_posmatchr_smoke_outputs(ann, outdir, "arabidopsis")
