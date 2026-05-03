# posMatchR

`posMatchR` annotates single-nucleotide `GRanges` sites with transcript context and constructs matched background site sets for post-transcriptional analyses such as m6A, A-to-I editing and precise RNA-binding maps.

The expected input is a single-nucleotide `GRanges`. Wider ranges are resized to width 1 by `prepare_sites()` and `annotate_sites()`.

## Installation

Install the package from GitHub:

```r
install.packages("remotes")
remotes::install_github("sarah-ku/posMatchR")
```

For the human hg38 examples below, install the required Bioconductor annotation packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GenomicRanges",
  "GenomeInfoDb",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "BSgenome.Hsapiens.UCSC.hg38"
))
```

For Arabidopsis examples, install the plant annotation packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GenomicRanges",
  "GenomeInfoDb",
  "TxDb.Athaliana.BioMart.plantsmart51",
  "org.At.tair.db",
  "BSgenome.Athaliana.TAIR.TAIR9"
))
```

## Basic human annotation workflow

This example assumes that you already have a single-nucleotide `GRanges` object or an `.rds`/`.RData`/`.Rdat` file containing one. The example uses human hg38/UCSC chromosome names.

```r
library(posMatchR)
library(GenomicRanges)
library(GenomeInfoDb)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Option 1: load a GRanges object from file.
sites <- load_sites("path/to/sites.rds")

# Option 2: if you already have a GRanges object, use it directly.
# sites <- my_sites_granges

human_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Standardise to width-1 sites and keep the main chromosomes.
sites <- prepare_sites(sites)
GenomeInfoDb::seqlevelsStyle(sites) <- "UCSC"
sites <- GenomeInfoDb::keepSeqlevels(
  sites,
  human_chrs,
  pruning.mode = "coarse"
)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

ann <- annotate_sites(
  gr = sites,
  txdb = txdb,
  seqstyle = "UCSC",
  chrs = human_chrs,
  tx_select = "longest",
  orgdb = org.Hs.eg.db::org.Hs.eg.db,
  gene_keytype = "ENTREZID",
  gene_symbol_col = "SYMBOL",
  gene_name_col = "GENENAME"
)

# Compact table for inspection or export.
basic <- as_basic_site_table(ann)
head(basic)

# Basic diagnostic plots.
plot_metagene_density(ann)
plot_junction_distance_density(ann)
```

The annotated `GRanges` keeps the original site coordinates and adds transcript, gene, region, metagene and local feature-geometry columns. For most reporting purposes, `as_basic_site_table()` gives a smaller table with the main annotations.
