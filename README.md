# posMatchR

<img src="man/figures/logo.png" align="right" width="180" alt="posMatchR logo" />

`posMatchR` annotates single-nucleotide `GRanges` sites with transcript context and builds matched background sets for post-transcriptional point-site analyses.

The package expects single-nucleotide sites. Wider ranges are resized to width one by `prepare_sites()` and `annotate_sites()`.

## Installation

```r
install.packages("remotes")
remotes::install_github("sarah-ku/posMatchR")
```

For the human hg38 example, install the required Bioconductor resources:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(c(
    "GenomicRanges",
    "GenomicFeatures",
    "VariantAnnotation",
    "GenomeInfoDb",
    "BSgenome",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db",
    "BSgenome.Hsapiens.UCSC.hg38"
))
```

## Minimal GLORI example

This example uses the bundled 25,000-site GLORI HEK293T `GRanges` object in `inst/extdata/GLORI_HEK293T_25K.Rdat`, restricted to the first 5,000 sites so the quick-start remains fast. The workflow annotates the positive sites, loads a precomputed human transcript-resource object, builds a candidate background universe from foreground-observed centred 5-mers with at least 100 foreground instances, and matches positives to negatives with exact centred 5-mer matching.

```r
library(posMatchR)
library(GenomeInfoDb)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
orgdb <- org.Hs.eg.db::org.Hs.eg.db
human_chrs <- paste0("chr", c(1:22, "X", "Y"))

sites <- load_sites(system.file(
    "extdata", "GLORI_HEK293T_25K.Rdat",
    package = "posMatchR",
    mustWork = TRUE
))
sites <- sites[seq_len(min(length(sites), 5000L))]

resources <- readRDS(system.file(
    "extdata", "hg38_knownGene_tx_resources_GLORI_demo.rds",
    package = "posMatchR",
    mustWork = TRUE
))

sites <- annotate_sites(
    gr = sites,
    txdb = txdb,
    resources = resources,
    seqstyle = "UCSC",
    chrs = human_chrs,
    tx_select = "longest",
    orgdb = orgdb,
    gene_keytype = "ENTREZID",
    gene_symbol_col = "SYMBOL",
    gene_name_col = "GENENAME",
    drop_unannotated = TRUE
)

sites <- add_kmer(
    gr = sites,
    genome = genome,
    k = 5,
    out_col = "kmer",
    seqstyle = "UCSC",
    chrs = human_chrs
)

background <- make_kmer_universe(
    foreground = sites,
    txdb = txdb,
    genome = genome,
    kmer_col = "kmer",
    min_count = 100,
    scope = "transcripts",
    regions = c("fiveUTR", "coding", "threeUTR"),
    resources = resources,
    exclude_foreground = TRUE
)

combined <- combine_site_sets(sites, background)
combined <- annotate_sites(
    gr = combined,
    txdb = txdb,
    resources = resources,
    seqstyle = "UCSC",
    chrs = human_chrs,
    tx_select = "longest",
    orgdb = orgdb,
    gene_keytype = "ENTREZID",
    gene_symbol_col = "SYMBOL",
    gene_name_col = "GENENAME",
    drop_unannotated = TRUE
)
combined <- add_kmer(combined, genome = genome, k = 5, seqstyle = "UCSC", chrs = human_chrs)

matched <- match_background(
    gr = combined,
    kmer_match = TRUE,
    kmer_col = "kmer",
    bin_match = FALSE
)
matched_pairs <- subset_matched_sets(matched)
```

Basic outputs:

```r
head(as_basic_site_table(sites))
plot_metagene_density(sites)
plot_junction_distance_density(sites)
plot_kmer_counts(sites, top_n = 15)

plot_metagene_density(matched_pairs)
plot_junction_distance_density(matched_pairs)
plot_kmer_counts(matched_pairs, top_n = 15)
head(summarise_matched_kmer_balance(matched_pairs))
```

`matched_pairs` contains the reciprocal positive and matched-negative sites. When `kmer_match = TRUE`, each retained pair has the same centred 5-mer.
