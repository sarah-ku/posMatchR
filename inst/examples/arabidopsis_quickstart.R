# Arabidopsis single-nucleotide site annotation example.
#
# This script assumes:
#   plant_sites: GRanges of observed plant sites, width 1
#
# The BioMart TxDb commonly uses seqlevels "1" to "5", whereas
# BSgenome.Athaliana.TAIR.TAIR9 uses "Chr1" to "Chr5". This script keeps the
# TxDb convention and renames the in-memory BSgenome/site seqlevels to match it.
# Mitochondrial and plastid contigs are dropped here because this TxDb contains
# only the five nuclear chromosomes.

library(posMatchR)
library(TxDb.Athaliana.BioMart.plantsmart51)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(org.At.tair.db)
library(GenomeInfoDb)

arab_chrs <- c("1", "2", "3", "4", "5")
genome <- BSgenome.Athaliana.TAIR.TAIR9::Athaliana
GenomeInfoDb::seqlevels(genome) <- c("1", "2", "3", "4", "5", "Mt", "Pt")
txdb <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51

GenomeInfoDb::seqlevels(plant_sites) <- sub("^Chr([1-5])$", "\\1", GenomeInfoDb::seqlevels(plant_sites))
plant_sites <- GenomeInfoDb::keepSeqlevels(plant_sites, arab_chrs, pruning.mode = "coarse")

ann <- annotate_sites(
  gr = plant_sites,
  txdb = txdb,
  chrs = arab_chrs,
  orgdb = org.At.tair.db::org.At.tair.db,
  gene_keytype = "TAIR",
  gene_symbol_col = "SYMBOL",
  gene_name_col = "GENENAME",
  tx_select = "longest"
)

ann <- add_kmer(ann, genome = genome, k = 5L, chrs = arab_chrs)

plot_metagene_density(ann)
plot_junction_distance_density(ann)
as_basic_site_table(ann)
