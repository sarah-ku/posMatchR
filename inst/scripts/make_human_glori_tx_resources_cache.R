#!/usr/bin/env Rscript

## Build the precomputed transcript-resource object used by the human GLORI
## vignette and README example. Run this script once from the package root
## after installing the human TxDb and posMatchR dependencies:
##
##   Rscript inst/scripts/make_human_glori_tx_resources_cache.R
##
## The script intentionally builds the full TxDb resources first, then reduces
## the object to transcripts overlapping the GLORI example sites used by the
## quick-start. This keeps the vignette fast while avoiding a large full-genome
## resource file in the software package. Set POSMATCHR_GLORI_CACHE_MAX_SITES=0
## to create a resource object for all bundled GLORI sites instead of the
## default first 5,000 sites.

pkg_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(pkg_root, "DESCRIPTION")) ||
    !dir.exists(file.path(pkg_root, "R"))) {
    stop("Run this script from the posMatchR package root.")
}

required <- c(
    "posMatchR",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "VariantAnnotation",
    "GenomeInfoDb",
    "S4Vectors"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
    stop(
        "Missing required packages: ", paste(missing, collapse = ", "),
        ". Install them before running this script."
    )
}

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
human_chrs <- paste0("chr", c(1:22, "X", "Y"))

glori_file <- file.path(pkg_root, "inst", "extdata", "GLORI_HEK293T_25K.Rdat")
if (!file.exists(glori_file)) {
    stop("Missing bundled GLORI file: ", glori_file)
}

message("Loading GLORI example sites ...")
sites <- posMatchR::load_sites(glori_file)

max_sites <- suppressWarnings(as.integer(Sys.getenv("POSMATCHR_GLORI_CACHE_MAX_SITES", "5000")))
if (length(max_sites) != 1L || is.na(max_sites)) max_sites <- 0L
if (max_sites > 0L && length(sites) > max_sites) {
    message("Restricting cache-generation sites to the first ", max_sites, " records.")
    sites <- sites[seq_len(max_sites)]
}

sites <- posMatchR::standardize_seqlevels(
    sites,
    target = txdb,
    seqstyle = "UCSC",
    chrs = human_chrs,
    keep = TRUE
)

message("Finding transcripts overlapping the GLORI example sites ...")
loc <- suppressWarnings(
    VariantAnnotation::locateVariants(
        sites,
        txdb,
        VariantAnnotation::AllVariants()
    )
)
loc_df <- as.data.frame(loc, stringsAsFactors = FALSE)
if (!nrow(loc_df) || !("TXID" %in% colnames(loc_df))) {
    stop("Could not identify overlapping transcript IDs for the GLORI sites.")
}
loc_txid <- unique(as.character(loc_df$TXID))
loc_txid <- loc_txid[!is.na(loc_txid) & nzchar(loc_txid)]

message("Building full human transcript resources. This is the slow step.")
full_resources <- posMatchR::build_tx_resources(txdb)

tx_metrics <- as.data.frame(full_resources$tx_metrics, stringsAsFactors = FALSE)
tx_key <- if (!is.null(full_resources$tx_key)) {
    as.character(full_resources$tx_key)
} else {
    as.character(tx_metrics$tx_name)
}

id_to_name <- stats::setNames(tx_key, as.character(tx_metrics$tx_id))
keep_tx <- unique(id_to_name[loc_txid])
keep_tx <- keep_tx[!is.na(keep_tx) & nzchar(keep_tx)]
keep_tx <- intersect(keep_tx, tx_key)

if (!length(keep_tx)) {
    stop("No overlapping transcript names were found in the transcript resources.")
}

message("Keeping ", length(keep_tx), " transcripts in the reduced resource object.")

subset_named <- function(x, keys) {
    if (is.null(x)) return(NULL)
    if (is.null(names(x))) return(x)
    x[intersect(keys, names(x))]
}

subset_list <- function(x, keys) {
    if (is.null(x)) return(NULL)
    if (!length(x)) return(x)
    x[intersect(keys, names(x))]
}

reduced <- full_resources
metric_keep <- tx_key %in% keep_tx
reduced$tx_metrics <- full_resources$tx_metrics[metric_keep, , drop = FALSE]
reduced$tx_key <- tx_key[metric_keep]

for (nm in c(
    "tx_len", "tx_strand", "cds_len", "five_len", "three_len",
    "intr_len", "cds_start", "cds_stop"
)) {
    reduced[[nm]] <- subset_named(full_resources[[nm]], keep_tx)
}

for (nm in c(
    "exons_by_tx", "exon_junctions", "cds_by_tx", "five_by_tx",
    "three_by_tx", "introns_by_tx"
)) {
    reduced[[nm]] <- subset_list(full_resources[[nm]], keep_tx)
}

reduced$cache_info <- list(
    resource_name = "hg38_knownGene_tx_resources_GLORI_demo",
    source_txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    source_sites = "GLORI_HEK293T_25K.Rdat",
    source_site_count = length(sites),
    max_sites_environment = Sys.getenv("POSMATCHR_GLORI_CACHE_MAX_SITES", "5000"),
    overlapping_transcript_count = length(keep_tx),
    created = as.character(Sys.time()),
    note = paste(
        "Reduced transcript-resource object for the posMatchR human GLORI",
        "vignette. Rebuild from the package root with",
        "Rscript inst/scripts/make_human_glori_tx_resources_cache.R"
    )
)

out_file <- file.path(
    pkg_root,
    "inst",
    "extdata",
    "hg38_knownGene_tx_resources_GLORI_demo.rds"
)

message("Writing ", out_file)
saveRDS(reduced, out_file, compress = "xz")

file_size <- file.info(out_file)$size
message("Wrote ", out_file, " (", round(file_size / 1024^2, 2), " MiB).")
if (is.finite(file_size) && file_size > 5 * 1024^2) {
    warning(
        "The resource file is larger than 5 MiB. Bioconductor software ",
        "packages currently warn on individual files larger than 5 MiB. ",
        "Consider reducing the demonstration resource further if BiocCheck ",
        "reports a file-size problem."
    )
}
