# Neuron RNA editing-style input conversion example.
#
# Starting from a data.frame/data.table final_sites with at least:
#   chrom, pos, biological_strand, site_id

library(posMatchR)
library(GenomicRanges)
library(IRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

sites_gr <- GRanges(
  seqnames = final_sites$chrom,
  ranges = IRanges(start = as.integer(final_sites$pos), end = as.integer(final_sites$pos)),
  strand = final_sites$biological_strand
)
S4Vectors::mcols(sites_gr)$site_id <- as.character(final_sites$site_id)
names(sites_gr) <- as.character(final_sites$site_id)

ann <- annotate_sites(
  gr = sites_gr,
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  orgdb = org.Hs.eg.db,
  gene_keytype = "ENTREZID",
  gene_symbol_col = "SYMBOL",
  gene_name_col = "GENENAME",
  seqstyle = "UCSC",
  chrs = paste0("chr", c(1:22, "X", "Y"))
)

ann_table <- as_site_table(ann, compatibility_names = TRUE)

# Example join back:
# final_sites_annotated <- merge(final_sites, ann_table, by = "site_id", all.x = TRUE)
