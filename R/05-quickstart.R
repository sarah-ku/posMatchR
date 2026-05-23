#' Quickstart: annotate sites and optionally build a matched background
#'
#' This is the simplest entry point for a vanilla single-nucleotide \code{GRanges}.
#' With only \code{sites} and \code{txdb}, it returns annotated sites. If
#' \code{background} is supplied, positives and background candidates are combined,
#' annotated, and matched.
#'
#' @param sites Positive/input sites as a \code{GRanges}.
#' @param txdb A \code{TxDb}.
#' @param background Optional candidate background \code{GRanges}.
#' @param genome Optional genome object for k-mer extraction.
#' @param orgdb Optional OrgDb for gene symbols.
#' @param seqstyle Optional sequence naming style.
#' @param chrs Optional seqlevels to keep.
#' @param tx_select Transcript selection rule.
#' @param label_col Label column.
#' @param site_id_col Site ID column.
#' @param k K-mer width when \code{genome} is supplied. Set NULL to skip.
#' @param match If TRUE, run \code{\link{match_background}}. Defaults to TRUE when
#'   \code{background} is supplied.
#' @param match_args Named list of additional arguments passed to \code{\link{match_background}}.
#' @param resources Optional precomputed transcript resources.
#' @param drop_unannotated If TRUE, drop unannotated rows.
#' @param make_plots If TRUE, include default ggplot objects in the returned list.
#' @param quiet If TRUE, suppress messages.
#'
#' @examples
#' if (interactive()) {
#'     txdb <- get(
#'         "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'         envir = asNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")
#'     )
#'     gr <- GenomicRanges::GRanges(
#'         "chr1", IRanges::IRanges(11874, width = 1), strand = "+"
#'     )
#'     posmatchr_quickstart(gr, txdb = txdb)
#' }
#'
#' @return A list with \code{gr}, \code{table}, \code{matched_pairs}, \code{diagnostics}, and optionally \code{plots}.
#' @export
posmatchr_quickstart <- function(sites,
                                 txdb,
                                 background = NULL,
                                 genome = NULL,
                                 orgdb = NULL,
                                 seqstyle = NULL,
                                 chrs = NULL,
                                 tx_select = "longest",
                                 label_col = "label",
                                 site_id_col = "site_id",
                                 k = 5L,
                                 match = !is.null(background),
                                 match_args = list(),
                                 resources = NULL,
                                 drop_unannotated = FALSE,
                                 make_plots = TRUE,
                                 quiet = FALSE) {
  if (is.null(background)) {
    gr <- prepare_sites(sites, label = NULL, label_col = label_col, site_id_col = site_id_col)
  } else {
    gr <- combine_site_sets(sites, background, label_col = label_col, site_id_col = site_id_col)
  }

  ann <- annotate_sites(
    gr = gr,
    txdb = txdb,
    tx_select = tx_select,
    seqstyle = seqstyle,
    chrs = chrs,
    resources = resources,
    drop_unannotated = drop_unannotated,
    orgdb = orgdb,
    site_id_col = site_id_col,
    preserve_mcols = TRUE,
    overwrite = TRUE,
    quiet = quiet
  )

  if (!is.null(genome) && !is.null(k)) {
    ann <- add_kmer(ann, genome = genome, k = k, out_col = "kmer", seqstyle = seqstyle, chrs = chrs)
  }

  if (isTRUE(match)) {
    if (is.null(background)) {
      warning("match=TRUE requested but no `background` was supplied; returning annotated sites only.")
    } else {
      args <- c(list(gr = ann, label_col = label_col), match_args)
      ann <- do.call(match_background, args)
    }
  }

  paired <- NULL
  if ("matched_negative_id" %in% colnames(S4Vectors::mcols(ann))) {
    paired <- tryCatch(subset_matched_pairs(ann, label_col = label_col), error = function(e) NULL)
  }

  plots <- list()
  if (isTRUE(make_plots)) {
    plots$metagene <- tryCatch(plot_metagene_density(ann), error = function(e) NULL)
    if ("nearest_exon_junction_dist" %in% colnames(S4Vectors::mcols(ann))) {
      plots$junction_distance <- tryCatch(plot_junction_distance_density(ann), error = function(e) NULL)
    }
    if ("start_dist" %in% colnames(S4Vectors::mcols(ann)) && "stop_dist" %in% colnames(S4Vectors::mcols(ann))) {
      plots$start_stop_distance <- tryCatch(plot_start_stop_distance_density(ann), error = function(e) NULL)
    }
    if ("kmer" %in% colnames(S4Vectors::mcols(ann))) {
      plots$kmer_counts <- tryCatch(plot_kmer_counts(ann), error = function(e) NULL)
    }
  }

  list(
    gr = ann,
    table = as_site_table(ann, compatibility_names = TRUE),
    matched_pairs = paired,
    diagnostics = list(
      posMatchR = S4Vectors::metadata(ann)$posMatchR,
      match = S4Vectors::metadata(ann)$match_diagnostics
    ),
    plots = plots
  )
}
