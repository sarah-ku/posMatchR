#' Standardise seqlevel names against a target object or chromosome list
#'
#' Renames seqlevels such as \code{1} to \code{Chr1} or \code{chr1} when a
#' matching target seqlevel is available. This is useful for cases where TxDb and
#' BSgenome packages for the same organism use different chromosome names.
#'
#' @param x A \code{GRanges}, \code{TxDb}, BSgenome, or another object with
#'   \code{GenomeInfoDb::seqlevels()}.
#' @param target Optional object or character vector providing the desired
#'   seqlevels. If supplied, \code{x} is renamed toward these names.
#' @param chrs Optional chromosome list to keep after renaming.
#' @param seqstyle Optional GenomeInfoDb seqlevel style to apply before direct
#'   name matching.
#' @param keep If TRUE and \code{chrs} is supplied, drop seqlevels not matching
#'   \code{chrs}.
#' @param pruning.mode Passed to \code{GenomeInfoDb::keepSeqlevels()} and
#'   seqlevel-renaming operations.
#'
#' @return \code{x} with renamed and optionally filtered seqlevels.
#' @export
standardize_seqlevels <- function(x,
                                  target = NULL,
                                  chrs = NULL,
                                  seqstyle = NULL,
                                  keep = FALSE,
                                  pruning.mode = "coarse") {
  if (!is.null(seqstyle)) {
    suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(x) <- seqstyle, silent = TRUE))
  }

  target_levels <- NULL
  if (!is.null(target)) {
    target_levels <- if (is.character(target)) target else GenomeInfoDb::seqlevels(target)
    x <- .rename_seqlevels_to_target(x, target_levels, pruning.mode = pruning.mode)
  }

  if (!is.null(chrs)) {
    if (is.null(target_levels)) {
      x <- .rename_seqlevels_to_target(x, chrs, pruning.mode = pruning.mode)
    }
    mapped <- .map_to_target_seqlevels(chrs, GenomeInfoDb::seqlevels(x))$mapped
    mapped <- unique(mapped[!is.na(mapped) & nzchar(mapped)])
    if (keep && length(mapped)) {
      x <- GenomeInfoDb::keepSeqlevels(x, mapped, pruning.mode = pruning.mode)
    }
  }

  x
}
#' Load a GRanges site set from an RDS or RData file
#'
#' A deliberately small file loader for examples and local smoke tests. For
#' \file{.rds} files, the file must contain a \code{GRanges}. For RData-style
#' files, the file must contain exactly one \code{GRanges}; if it contains more
#' than one, load it manually and pass the object to \code{prepare_sites()}.
#'
#' @param path Path to an \code{.rds}, \code{.RDS}, \code{.rda}, \code{.RData},
#'   or \code{.Rdat} file.
#'
#' @return A prepared single-nucleotide \code{GRanges}.
#' @export
load_sites <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) stop("`path` must be a single file path.")
  if (!file.exists(path)) stop("File does not exist: ", path)

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    gr <- readRDS(path)
  } else {
    env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = env)
    if (!length(loaded)) stop("No objects were loaded from: ", path)

    is_gr <- vapply(loaded, function(nm) inherits(get(nm, envir = env), "GRanges"), logical(1))
    gr_names <- loaded[is_gr]
    if (length(gr_names) != 1L) {
      stop(
        "Could not infer a single GRanges object from ", path,
        ". Load the file manually and call prepare_sites() on the object. GRanges objects found: ",
        paste(gr_names, collapse = ", ")
      )
    }
    gr <- get(gr_names[1L], envir = env)
  }

  if (!inherits(gr, "GRanges")) stop("Loaded object is not a GRanges.")
  prepare_sites(gr)
}

#' Import single-nucleotide sites from a BED file
#'
#' Uses \code{rtracklayer::import()} to read a BED-like file as \code{GRanges},
#' then applies \code{prepare_sites()}. BED coordinate conversion is handled by
#' rtracklayer.
#'
#' @param path Path to a BED file.
#'
#' @return A prepared single-nucleotide \code{GRanges}.
#' @export
import_bed_sites <- function(path) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required for BED import. Install it with BiocManager::install('rtracklayer').")
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) stop("`path` must be a single file path.")
  if (!file.exists(path)) stop("File does not exist: ", path)

  gr <- rtracklayer::import(path, format = "BED")
  if (!inherits(gr, "GRanges")) stop("rtracklayer::import() did not return a GRanges.")

  mc <- S4Vectors::mcols(gr)
  if (!("site_id" %in% colnames(mc)) && "name" %in% colnames(mc)) {
    S4Vectors::mcols(gr)[["site_id"]] <- as.character(mc$name)
  }

  prepare_sites(gr)
}

#' Return a small display table for annotated sites
#'
#' @param gr An annotated \code{GRanges}.
#' @param compatibility_names Passed to \code{as_site_table()}.
#'
#' @return A compact data.frame.
#' @export
as_basic_site_table <- function(gr, compatibility_names = TRUE) {
  as_site_table(gr, compatibility_names = compatibility_names, columns = "basic")
}

#' Nest less frequently used metadata columns inside one DataFrame column
#'
#' This is intended for display or storage after annotation/matching. Downstream
#' posMatchR plotting and matching functions expect scalar columns such as
#' \code{metagene_split3} to remain top-level columns, so compact after those steps
#' rather than before them.
#'
#' @param gr An annotated \code{GRanges}.
#' @param keep_cols Metadata columns to keep at top level.
#' @param nested_col Name of the nested metadata column.
#'
#' @return A \code{GRanges} with compacted metadata.
#' @export
compact_site_mcols <- function(gr,
                               keep_cols = c(
                                 "site_id", "label", "location",
                                 "gene_id", "gene_symbol", "gene_name", "tx_id", "tx_name",
                                 "metagene_split3", "kmer", "match_set", "matched_negative_id",
                                 "matched_positive_id"
                               ),
                               nested_col = "posmatchr_details") {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  mc <- S4Vectors::mcols(gr)
  if (!ncol(mc)) return(gr)

  keep <- intersect(keep_cols, colnames(mc))
  detail <- setdiff(colnames(mc), c(keep, nested_col))
  if (!length(detail)) return(gr)

  details <- mc[, detail, drop = FALSE]
  detail_rows <- lapply(seq_len(length(gr)), function(i) details[i, , drop = FALSE])
  new_mc <- mc[, keep, drop = FALSE]
  new_mc[[nested_col]] <- detail_rows
  S4Vectors::mcols(gr) <- new_mc

  md <- S4Vectors::metadata(gr)
  if (is.null(md$posMatchR)) md$posMatchR <- list()
  md$posMatchR$compacted_columns <- detail
  md$posMatchR$compacted_into <- nested_col
  S4Vectors::metadata(gr) <- md

  gr
}
