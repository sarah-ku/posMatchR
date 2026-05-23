# Example data

`GLORI_HEK293T_25K.Rdat` contains one `GRanges` object named `GLORI` with
25,000 human HEK293T GLORI m6A sites for the README and vignette examples.

`hg38_knownGene_tx_resources_GLORI_demo.rds` is a reduced precomputed transcript-
resource object for `TxDb.Hsapiens.UCSC.hg38.knownGene`, restricted to
transcripts overlapping the GLORI sites used in the README and vignette. It is
used only to keep examples and vignette builds fast. Rebuild it from the package root with:

```sh
Rscript inst/scripts/make_human_glori_tx_resources_cache.R
```
