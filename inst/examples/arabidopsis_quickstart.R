# Arabidopsis single-nucleotide site annotation example.
#
# This script assumes:
#   plant_sites: GRanges of observed plant sites, width 1
#
# The BioMart TxDb commonly uses seqlevels "1" to "5", whereas
# BSgenome.Athaliana.TAIR.TAIR9 uses "Chr1" to "Chr5". This script keeps the
# seqlevel handling explicit rather than hiding it in the loader. Mitochondrial
# and plastid contigs are dropped here because this TxDb contains only the five
# nuclear chromosomes.

library(posMatchR)
library(TxDb.Athaliana.BioMart.plantsmart51)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(org.At.tair.db)

arab_chrs <- paste0("Chr", 1:5)
genome <- BSgenome.Athaliana.TAIR.TAIR9::Athaliana
txdb <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
GenomeInfoDb::seqlevels(txdb) <- arab_chrs

GenomeInfoDb::seqlevels(plant_sites) <- sub("^([1-5])$", "Chr\\1", GenomeInfoDb::seqlevels(plant_sites))
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
