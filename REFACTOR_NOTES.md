# posMatchR refactor notes

This refactor narrows the package around the intended paper use case: single-nucleotide `GRanges` sites as the primary input, transcript-aware annotation, metagene and feature-distance variables, matched background construction, and simple diagnostic plots.

## Main design choices

The package now treats the input as a point-site object. `prepare_sites()` and `annotate_sites()` resize wider ranges to width 1 by default, but this can be changed to an error with `width_action = "error"` in `prepare_sites()`. Site IDs are stored in `site_id` and also used as `names(gr)`, which makes later matching and joins more stable.

Existing metadata are preserved by default. Re-running `annotate_sites()` replaces only columns that posMatchR itself generates. This avoids the problem of accumulating stale annotation columns while not discarding useful user columns such as editing level, coverage, sample labels, or experiment-specific scores. Use `preserve_mcols = FALSE` or `strip_mcols = TRUE` only when a clean object is explicitly wanted.

`annotate_sites()` uses `VariantAnnotation::locateVariants()` for region classes and then attaches transcript metrics from a precomputed resource object. The user-facing region columns are:
- `location`: raw `VariantAnnotation` class, for example `coding`, `fiveUTR`, `threeUTR`, `intron`, `spliceSite`, `promoter`, `intergenic`, or `unannotated`.
- `feature` and `region_class`: simplified labels such as `CDS`, `5UTR`, and `3UTR`.

The helper bug that required patching `queryHits()`/`subjectHits()` externally has been removed by explicitly calling `S4Vectors::queryHits()` and `S4Vectors::subjectHits()` inside package functions.

## Metagene and geometry variables

The package now stores both simple and split metagene coordinates:
- `metagene_prop`: a 0--3 axis where 5'UTR, CDS, and 3'UTR each occupy one unit.
- `metagene_split`: a 0--1 axis using transcript-length-derived 5'UTR/CDS/3'UTR widths.
- `metagene_split3`: the same split axis scaled to 0--3 for plotting.

Feature distances include `nearest_exon_junction_dist`, `dist_from_feature_start`, `dist_from_feature_end`, `feature_width`, and transcript-coordinate distances to CDS anchors (`start_dist_tx`, `stop_dist_tx`). The older genomic-coordinate `start_dist` and `stop_dist` columns are retained for compatibility.

## Matching strategy

`match_background()` uses hard strata for region, optional k-mer, and optional train/test split. Within each stratum it greedily selects one unique background site per positive site by metagene distance, with optional quantile-binned covariates such as exon-junction distance, CDS-anchor distance, transcript length, feature width, and feature-boundary distances.

This remains deliberately simple and fast. It is not a full optimal transport or propensity-score implementation. For motif-enrichment applications, the important practical point is that the candidate background set should already represent the experimental search space: for example all assayed DRACH adenosines for m6A, or all assayed adenosines passing coverage/expression filters for A-to-I editing. posMatchR then balances against annotation and geometry biases.

## Plant use

The refactor adds clearer `seqstyle` and `chrs` handling. For Arabidopsis, passing `chrs = c("1", "2", "3", "4", "5")` drops mitochondrial/plastid/non-standard contigs when they are absent from the TxDb. This is intentional and should happen early, with a warning, rather than failing later during overlap or transcript mapping.

## Tests and examples

The package contains lightweight tests for:
- preparation and stable IDs,
- mock matched-background behaviour,
- scalar table conversion,
- plot constructors,
- optional human and Arabidopsis annotation smoke tests.

Examples are under `inst/examples/`:
- `human_quickstart.R`
- `arabidopsis_quickstart.R`
- `neuron_editing_input.R`

## Validation caveat

This archive was statically reviewed and assembled in a sandbox that does not contain an R executable, so `R CMD check` could not be run here. Run the following locally after unpacking:

```r
devtools::document()
devtools::test()
devtools::check()
```

or from a shell:

```sh
R CMD build posMatchR
R CMD check posMatchR_*.tar.gz
```
