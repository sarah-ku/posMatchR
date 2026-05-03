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

## Basic matched-background workflow

To construct a matched background, provide a second single-nucleotide `GRanges` containing candidate background sites. For m6A, this might be a set of eligible adenosines or eligible DRACH-centred adenosines. The candidates should be defined before using `posMatchR`, based on the biological question and the data quality filters appropriate for the experiment.

```r
library(posMatchR)
library(GenomeInfoDb)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

human_chrs <- paste0("chr", c(1:22, "X", "Y"))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

positives <- load_sites("path/to/positive_sites.rds")
background <- load_sites("path/to/candidate_background_sites.rds")

positives <- prepare_sites(positives)
background <- prepare_sites(background)

GenomeInfoDb::seqlevelsStyle(positives) <- "UCSC"
GenomeInfoDb::seqlevelsStyle(background) <- "UCSC"

positives <- GenomeInfoDb::keepSeqlevels(positives, human_chrs, pruning.mode = "coarse")
background <- GenomeInfoDb::keepSeqlevels(background, human_chrs, pruning.mode = "coarse")

combined <- combine_site_sets(
  positives = positives,
  background = background,
  label_col = "label",
  site_id_col = "site_id"
)

combined <- annotate_sites(
  gr = combined,
  txdb = txdb,
  seqstyle = "UCSC",
  chrs = human_chrs,
  tx_select = "longest",
  orgdb = org.Hs.eg.db::org.Hs.eg.db,
  gene_keytype = "ENTREZID",
  gene_symbol_col = "SYMBOL",
  gene_name_col = "GENENAME"
)

matched <- match_background(
  combined,
  label_col = "label",
  meta_col = "metagene_split3",
  meta_tol = 0.03,
  enforce_meta_tol = TRUE,
  bin_match = TRUE,
  kmer_match = FALSE,
  seed = 1L
)

matched_sets <- subset_matched_sets(matched)

plot_metagene_density(matched_sets, set_col = "match_set")
plot_junction_distance_density(matched_sets, set_col = "match_set")

S4Vectors::metadata(matched)$match_diagnostics
```

If exact 5-mer balancing is needed, add sequence context before matching and set `kmer_match = TRUE`:

```r
library(BSgenome.Hsapiens.UCSC.hg38)

combined <- add_kmer(
  combined,
  genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
  k = 5L,
  seqstyle = "UCSC",
  chrs = human_chrs
)

matched_kmer <- match_background(
  combined,
  label_col = "label",
  meta_col = "metagene_split3",
  meta_tol = 0.03,
  enforce_meta_tol = TRUE,
  bin_match = TRUE,
  kmer_match = TRUE,
  kmer_col = "kmer",
  seed = 1L
)

matched_kmer_sets <- subset_matched_sets(matched_kmer)
plot_kmer_counts(matched_kmer_sets, set_col = "match_set", top_n = 25)
```

## Useful output helpers

```r
# Small table with the main fields.
basic <- as_basic_site_table(ann)

# Full table with all metadata columns.
full <- as_site_table(ann, columns = "all")

# Matched rows only, for plotting or export.
matched_sets <- subset_matched_sets(matched)
```
