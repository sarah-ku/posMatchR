# Example data

`GLORI_HEK293T_25K.Rdat` contains one `GRanges` object named `GLORI` with 25,000 human HEK293T GLORI m6A sites for the README and vignette examples. It is subsampled from the Liu et al. 2022 dataset (see below) and saved out as a `GRanges` object for ease of package demonstration.

`hg38_knownGene_tx_resources_GLORI_demo.rds` is a reduced precomputed transcript-resource object for `TxDb.Hsapiens.UCSC.hg38.knownGene`, restricted to transcripts overlapping the GLORI sites used in the README and vignette. It is used only to keep examples and vignette builds fast. This file can be rebuilt from the package root with:

```sh
Rscript inst/scripts/make_human_glori_tx_resources_cache.R
```
# Exeternal datasets

We further include the external datasets used for testing and demonstration in their full form (supplementary tables of processed data containing site coordinates) are as follows.

GLORI HEK293T m6A sites (Liu et al. 2022) https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01487-9/MediaObjects/41587_2022_1487_MOESM3_ESM.xlsx

m6A sites in Arabidopsis (Parker et al. 2020) https://cdn.elifesciences.org/articles/49658/elife-49658-fig5-data1-v1.tds

ECT iCLIP (Arribas-Hernández et al. 2021) https://cdn.elifesciences.org/articles/72375/elife-72375-supp2-v2.xlsx
