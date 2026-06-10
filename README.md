# posMatchR

<img src="man/figures/logo.png" align="right" width="180" alt="posMatchR logo" />

`posMatchR` annotates single-nucleotide `GRanges` sites with transcript name, context and features, as well as sequence context. It then can be used to build matched and/or kmer-sequence balanced background sets useful for post-transcriptional point-site analyses, such as motif analysis, comparisons with local RNA features or metagene plots.

The package expects single-nucleotide sites (e.g. genomic positions for m6A sites). Wider ranges are automatically resized to width one during data prepration and annotation.

## Sequence-universe functions

In order to generate a set of background candidates (e.g. for motif analysis), `make_kmer_universe()` searches exhaustively for exact occurrences of foreground-observed k-mers (RNA `U` is converted to DNA `T` and k-mers must be same-width `A`/`C`/`G`/`T` strings).

`make_motif_universe()` searches for IUPAC motifs (e.g. such as `DRACH` for m6A analysis) using `Biostrings::matchPattern()` functionality. Note that for motifs like `DRACH`, the default site is the central `A` because `site_offset` defaults to the motif centre.

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

This example uses the bundled 25,000-site GLORI HEK293T `GRanges` object in `inst/extdata/GLORI_HEK293T_25K.Rdat`, restricted to the first 5,000 sites so the quick-start remains fast. The workflow annotates the positive sites, loads a precomputed human transcript-resource object (this can be generated in the package for the species of interest, but we supply a cached version for hg38 since it takes some time), builds a candidate background universe from foreground-observed centred 5-mers with at least 100 foreground instances, and matches positives to negatives with exact centred 5-mer matching.

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

## References

The following datasets were used for package testing and demonstration purposes:

GLORI HEK293T m6A sites (Liu et al. 2022)
https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01487-9/MediaObjects/41587_2022_1487_MOESM3_ESM.xlsx
 
m6A sites in Arabidopsis (Parker et al. 2020)
https://cdn.elifesciences.org/articles/49658/elife-49658-fig5-data1-v1.tds
 
ECT iCLIP (Arribas-Hernández et al. 2021)
https://cdn.elifesciences.org/articles/72375/elife-72375-supp2-v2.xlsx

These files are supplied in the repository https://github.com/sarah-ku/posMatchR/tree/master/inst/extdata together with processed GRanges for getting started with the R package.

