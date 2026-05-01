# Human single-nucleotide site annotation and matched background example.
#
# This script assumes:
#   observed_sites: GRanges of observed sites, width 1
#   candidate_background_sites: GRanges of eligible background sites, width 1
#
# For m6A, candidate_background_sites might be all assayed adenosines or all assayed
# DRACH-centred adenosines after expression/coverage filters.

library(posMatchR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg38

res <- posmatchr_quickstart(
  sites = observed_sites,
  background = candidate_background_sites,
  txdb = txdb,
  genome = genome,
  orgdb = org.Hs.eg.db,
  seqstyle = "UCSC",
  chrs = paste0("chr", c(1:22, "X", "Y")),
  tx_select = "longest",
  k = 5,
  match_args = list(
    kmer_match = TRUE,
    meta_col = "metagene_split3",
    bin_match = TRUE,
    bin_cols = c("nearest_exon_junction_dist", "start_dist_tx", "stop_dist_tx", "tx_len"),
    bin_mode = "soft",
    meta_tol = 0.05
  )
)

annotated_sites <- res$gr
matched_pairs <- res$matched_pairs
annotation_table <- res$table

res$plots$metagene
res$plots$junction_distance
res$plots$start_stop_distance
res$plots$kmer_counts
