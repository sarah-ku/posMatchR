# Arabidopsis single-nucleotide site annotation example.
#
# This script assumes:
#   plant_sites: GRanges of observed plant sites, width 1
#
# TxDb.Athaliana.BioMart.plantsmart51 generally supports the main chromosomes "1" to "5".
# Mitochondrial/plastid contigs in the input are dropped when chrs is restricted below.

library(posMatchR)
library(TxDb.Athaliana.BioMart.plantsmart51)
library(org.At.tair.db)

txdb <- TxDb.Athaliana.BioMart.plantsmart51

ann <- annotate_sites(
  gr = plant_sites,
  txdb = txdb,
  chrs = c("1", "2", "3", "4", "5"),
  orgdb = org.At.tair.db,
  gene_keytype = "TAIR",
  gene_symbol_col = "SYMBOL",
  gene_name_col = "GENENAME",
  tx_select = "longest"
)

plot_metagene_density(ann)
plot_junction_distance_density(ann)
as_site_table(ann)
