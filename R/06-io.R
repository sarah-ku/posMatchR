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
#' Loads a \code{GRanges} object and normalises it for posMatchR. If an RData file
#' contains exactly one \code{GRanges}, the object name can be omitted.
#'
#' @param path Path to an \code{.rds}, \code{.RDS}, \code{.rda}, \code{.RData},
#'   or \code{.Rdat} file.
#' @param object Optional object name for RData-style files.
#' @param label Optional label passed to \code{prepare_sites()}.
#' @param label_col Metadata column used for labels.
#' @param site_id_col Metadata column used for stable site IDs.
#' @param width_action Either \code{"resize"} or \code{"error"}.
#' @param strip_mcols If TRUE, remove existing metadata except the generated site ID
#'   and optional label.
#' @param seqstyle Optional GenomeInfoDb seqlevel style.
#' @param chrs Optional seqlevels to keep.
#'
#' @return A prepared single-nucleotide \code{GRanges}.
#' @export
load_sites <- function(path,
                       object = NULL,
                       label = NULL,
                       label_col = "label",
                       site_id_col = "site_id",
                       width_action = c("resize", "error"),
                       strip_mcols = FALSE,
                       seqstyle = NULL,
                       chrs = NULL) {
  width_action <- match.arg(width_action)
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) stop("`path` must be a single file path.")
  if (!file.exists(path)) stop("File does not exist: ", path)

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    gr <- readRDS(path)
  } else {
    env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = env)
    if (!length(loaded)) stop("No objects were loaded from: ", path)

    if (!is.null(object)) {
      if (!(object %in% loaded)) {
        stop("Object '", object, "' was not found in ", path, ". Available objects: ", paste(loaded, collapse = ", "))
      }
      gr <- get(object, envir = env)
    } else {
      is_gr <- vapply(loaded, function(nm) inherits(get(nm, envir = env), "GRanges"), logical(1))
      gr_names <- loaded[is_gr]
      if (length(gr_names) != 1L) {
        stop(
          "Could not infer a single GRanges object from ", path,
          ". Set `object=` explicitly. GRanges objects found: ", paste(gr_names, collapse = ", ")
        )
      }
      gr <- get(gr_names[1L], envir = env)
    }
  }

  if (!inherits(gr, "GRanges")) stop("Loaded object is not a GRanges.")

  if (!is.null(seqstyle)) suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- seqstyle, silent = TRUE))
  if (!is.null(chrs)) {
    gr <- .rename_seqlevels_to_target(gr, chrs)
    mapped <- .map_to_target_seqlevels(chrs, GenomeInfoDb::seqlevels(gr))$mapped
    mapped <- unique(mapped[!is.na(mapped) & nzchar(mapped)])
    if (length(mapped)) gr <- GenomeInfoDb::keepSeqlevels(gr, mapped, pruning.mode = "coarse")
  }

  prepare_sites(
    gr,
    label = label,
    label_col = label_col,
    site_id_col = site_id_col,
    width_action = width_action,
    strip_mcols = strip_mcols
  )
}

#' Import single-nucleotide sites from a BED file
#'
#' Uses \code{rtracklayer::import()} to read BED-like files as \code{GRanges},
#' then prepares the result for posMatchR. BED coordinates are handled by
#' rtracklayer.
#'
#' @param path Path to a BED file.
#' @param genome Optional genome argument forwarded to \code{rtracklayer::import()}.
#' @param label Optional label passed to \code{prepare_sites()}.
#' @param label_col Metadata column used for labels.
#' @param site_id_col Metadata column used for stable site IDs. If the imported BED
#'   has a \code{name} column and no site ID, it is copied to this column.
#' @param width_action Either \code{"resize"} or \code{"error"}.
#' @param strip_mcols If TRUE, remove existing metadata except site ID and label.
#' @param seqstyle Optional GenomeInfoDb seqlevel style.
#' @param chrs Optional seqlevels to keep.
#' @param ... Additional arguments passed to \code{rtracklayer::import()}.
#'
#' @return A prepared single-nucleotide \code{GRanges}.
#' @export
import_bed_sites <- function(path,
                             genome = NULL,
                             label = NULL,
                             label_col = "label",
                             site_id_col = "site_id",
                             width_action = c("resize", "error"),
                             strip_mcols = FALSE,
                             seqstyle = NULL,
                             chrs = NULL,
                             ...) {
  width_action <- match.arg(width_action)
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required for BED import. Install it with BiocManager::install('rtracklayer').")
  }
  if (!file.exists(path)) stop("File does not exist: ", path)

  args <- list(con = path, format = "BED", ...)
  if (!is.null(genome)) args$genome <- genome
  gr <- do.call(rtracklayer::import, args)

  if (!inherits(gr, "GRanges")) stop("rtracklayer::import() did not return a GRanges.")
  mc <- S4Vectors::mcols(gr)
  if (!(site_id_col %in% colnames(mc)) && "name" %in% colnames(mc)) {
    S4Vectors::mcols(gr)[[site_id_col]] <- as.character(mc$name)
  }

  if (!is.null(seqstyle)) suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- seqstyle, silent = TRUE))
  if (!is.null(chrs)) {
    gr <- .rename_seqlevels_to_target(gr, chrs)
    mapped <- .map_to_target_seqlevels(chrs, GenomeInfoDb::seqlevels(gr))$mapped
    mapped <- unique(mapped[!is.na(mapped) & nzchar(mapped)])
    if (length(mapped)) gr <- GenomeInfoDb::keepSeqlevels(gr, mapped, pruning.mode = "coarse")
  }

  prepare_sites(
    gr,
    label = label,
    label_col = label_col,
    site_id_col = site_id_col,
    width_action = width_action,
    strip_mcols = strip_mcols
  )
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
                                 "site_id", "label", "location", "feature", "region_class",
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
