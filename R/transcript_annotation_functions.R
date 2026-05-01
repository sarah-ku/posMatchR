#' @import GenomicRanges
#' @import GenomicFeatures
#' @import VariantAnnotation
#' @importFrom S4Vectors mcols DataFrame elementNROWS
#' @importFrom IRanges resize
#' @importFrom GenomeInfoDb keepSeqlevels seqlevelsStyle seqinfo
#' @importFrom eulerr euler
NULL


#' Internal helper
#' @noRd
.validate_sites_gr <- function(gr, label_col = "label") {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")

  # Enforce width=1 for point-site annotations
  if (any(GenomicRanges::width(gr) != 1L)) {
    warning("Some ranges have width != 1. Resizing to width=1 (fix='center').")
    gr <- IRanges::resize(gr, width = 1L, fix = "center")
  }

  # Stable IDs
  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_along(gr))
  }

  mc <- S4Vectors::mcols(gr)
  if (!is.null(label_col) && (label_col %in% colnames(mc))) {
    y <- mc[[label_col]]

    if (is.logical(y)) {
      mc[[label_col]] <- as.integer(y)
    } else if (is.factor(y)) {
      y_chr <- as.character(y)
      ok01 <- (!is.na(y_chr)) & (y_chr %in% c("0", "1"))
      if (all(is.na(y_chr) | ok01)) mc[[label_col]] <- as.integer(y_chr)
    } else if (is.character(y)) {
      ok01 <- (!is.na(y)) & (y %in% c("0", "1"))
      if (all(is.na(y) | ok01)) mc[[label_col]] <- as.integer(y)
    }
  }
  S4Vectors::mcols(gr) <- mc

  gr
}


#' Internal helper
#' @noRd
.validate_input_gr <- function(gr, label_col = "label", require_both = TRUE) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  if (!(label_col %in% colnames(S4Vectors::mcols(gr)))) {
    stop("`gr` must contain a metadata column named `", label_col, "`.")
  }

  y <- S4Vectors::mcols(gr)[[label_col]]
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y))  y <- as.integer(as.character(y))
  y <- as.integer(y)

  if (any(!is.na(y) & !(y %in% c(0L, 1L)))) {
    stop("`", label_col, "` must be 0/1 (or coercible to 0/1).")
  }

  if (require_both) {
    if (sum(y == 1L, na.rm = TRUE) < 1L) stop("No positives found (", label_col, "==1).")
    if (sum(y == 0L, na.rm = TRUE) < 1L) stop("No negatives found (", label_col, "==0).")
  }

  if (any(GenomicRanges::width(gr) != 1L)) {
    warning("Some ranges have width != 1. Resizing to width=1 (fix='center').")
    gr <- IRanges::resize(gr, width = 1L, fix = "center")
  }

  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_along(gr))
  }

  S4Vectors::mcols(gr)[[label_col]] <- y
  gr
}


#' Internal helper
#' @noRd
.harmonize_seqlevels <- function(gr, txdb, seqstyle = NULL, chrs = NULL) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")

  # TxDb objects can behave reference-y in some contexts (notably when you pass the
  # package-level singleton directly). Force copy-on-modify by creating a second
  # reference before any seqlevels/seqinfo mutation.
  txdb_ref <- txdb
  txdb <- txdb_ref

  txdb_seqs <- GenomeInfoDb::seqlevels(txdb)
  if (length(txdb_seqs) == 0L) {
    stop(
      "`txdb` has no seqlevels. This usually means it was previously pruned to zero ",
      "seqlevels (e.g. by keepSeqlevels). Reload/recreate the TxDb object (or restart R) ",
      "and try again."
    )
  }

  # Apply seqlevel style harmonization
  if (!is.null(seqstyle)) {
    GenomeInfoDb::seqlevelsStyle(gr) <- seqstyle
    GenomeInfoDb::seqlevelsStyle(txdb) <- seqstyle
  } else {
    st <- tryCatch(GenomeInfoDb::seqlevelsStyle(txdb), error = function(e) character(0))
    if (length(st) > 0L) {
      suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- st[1], silent = TRUE))
    }
  }

  # Helper: try to map requested seqlevels onto target seqlevels using simple
  # prefix/case heuristics (handles "1" vs "chr1"/"Chr1" style issues).
  .map_to_target_seqlevels <- function(requested, target) {
    requested <- unique(as.character(requested))
    target <- as.character(target)

    out <- requested
    ok <- out %in% target

    # Try adding "chr"
    idx <- which(!ok)
    if (length(idx)) {
      cand <- paste0("chr", out[idx])
      m <- match(cand, target)
      hit <- !is.na(m)
      if (any(hit)) {
        out[idx[hit]] <- target[m[hit]]
        ok[idx[hit]] <- TRUE
      }
    }

    # Try adding "Chr"
    idx <- which(!ok)
    if (length(idx)) {
      cand <- paste0("Chr", out[idx])
      m <- match(cand, target)
      hit <- !is.na(m)
      if (any(hit)) {
        out[idx[hit]] <- target[m[hit]]
        ok[idx[hit]] <- TRUE
      }
    }

    # Try stripping "chr"/"Chr"
    idx <- which(!ok)
    if (length(idx)) {
      cand <- sub("^chr", "", out[idx], ignore.case = TRUE)
      m <- match(cand, target)
      hit <- !is.na(m)
      if (any(hit)) {
        out[idx[hit]] <- target[m[hit]]
        ok[idx[hit]] <- TRUE
      }
    }

    # Case-insensitive exact matching
    idx <- which(!ok)
    if (length(idx)) {
      m <- match(tolower(out[idx]), tolower(target))
      hit <- !is.na(m)
      if (any(hit)) {
        out[idx[hit]] <- target[m[hit]]
        ok[idx[hit]] <- TRUE
      }
    }

    list(mapped = unique(out[ok]), dropped = unique(out[!ok]))
  }

  # Decide which chromosomes to keep
  if (is.null(chrs)) {
    chrs_final <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(txdb))
  } else {
    chrs_in <- unique(as.character(chrs))
    mapped <- .map_to_target_seqlevels(chrs_in, GenomeInfoDb::seqlevels(txdb))
    if (length(mapped$dropped)) {
      warning(
        "Dropping chromosomes not present in txdb after harmonization: ",
        paste(mapped$dropped, collapse = ", ")
      )
    }
    chrs_final <- mapped$mapped
    # Also require presence in gr (after any style conversion)
    chrs_final <- intersect(chrs_final, GenomeInfoDb::seqlevels(gr))
  }

  if (length(chrs_final) == 0L) {
    stop(
      "No overlapping seqlevels between `gr` and `txdb` after harmonization.\n",
      "seqlevels(gr):  ", paste(GenomeInfoDb::seqlevels(gr), collapse = ", "), "\n",
      "seqlevels(txdb): ", paste(GenomeInfoDb::seqlevels(txdb), collapse = ", "), "\n",
      if (!is.null(chrs)) paste0("requested chrs: ", paste(unique(as.character(chrs)), collapse = ", "), "\n") else ""
    )
  }

  gr   <- GenomeInfoDb::keepSeqlevels(gr,   chrs_final, pruning.mode = "coarse")
  txdb <- GenomeInfoDb::keepSeqlevels(txdb, chrs_final, pruning.mode = "coarse")

  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(txdb)[chrs_final]

  gr <- suppressWarnings(GenomicRanges::trim(gr))

  list(gr = gr, txdb = txdb, chrs = chrs_final)
}

#' Internal helper
#' @noRd
.first_or_na <- function(x) {
  if (is(x, "List") || is.list(x)) {
    return(vapply(x, function(z) if (length(z)) as.character(z[[1]]) else NA_character_, character(1)))
  }
  as.character(x)
}

# ----------------------------
# Transcript metrics & resources
# ----------------------------

#' Compute transcript-level length metrics
#'
#' Computes transcript lengths and (optionally) intron counts/lengths from a \code{TxDb}.
#' This is used internally by \code{\link{build_tx_resources}}, but can also be useful
#' when you want the raw per-transcript metric table.
#'
#' @param txdb A \code{TxDb} object.
#' @param include_introns Logical; if TRUE, compute intron counts and total intron length per transcript.
#' @param introns_by_tx Optional \code{GRangesList} from \code{GenomicFeatures::intronsByTranscript(txdb)}.
#'   Supplying this avoids recomputation.
#'
#' @returns A \code{data.frame} with one row per transcript (tx), including \code{tx_len},
#'   \code{cds_len}, \code{utr5_len}, \code{utr3_len}, and optionally \code{intron_len}, \code{n_intron}.
#'
#' @seealso \code{\link{build_tx_resources}}
#' @export
compute_tx_metrics <- function(txdb, include_introns = TRUE, introns_by_tx = NULL) {
  tl <- GenomicFeatures::transcriptLengths(
    txdb,
    with.cds_len  = TRUE,
    with.utr5_len = TRUE,
    with.utr3_len = TRUE
  )
  tl <- as.data.frame(tl, stringsAsFactors = FALSE)
  tl$tx_id   <- as.character(tl$tx_id)
  tl$tx_name <- as.character(tl$tx_name)

  if (include_introns) {
    if (is.null(introns_by_tx)) {
      introns_by_tx <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
    }

    n <- S4Vectors::elementNROWS(introns_by_tx)
    intr_n <- n
    names(intr_n) <- names(introns_by_tx)

    intr_len <- numeric(length(n))
    names(intr_len) <- names(introns_by_tx)

    if (sum(n) > 0L) {
      u <- unlist(introns_by_tx, use.names = FALSE)
      grp <- rep.int(seq_along(n), n)

      # grouped sum of widths; rownames are group indices
      intr_sum <- rowsum(width(u), grp, reorder = FALSE)
      intr_len[as.integer(rownames(intr_sum))] <- intr_sum[, 1]
    }

    tl$intron_len <- intr_len[match(tl$tx_name, names(intr_len))]
    tl$n_intron   <- intr_n[match(tl$tx_name, names(intr_n))]
  }

  tl
}

#' Internal helper: compute metagene region splits (median/mean-derived)
#'
#' Mimics the "median-derived region splits" idea: estimate a typical start-codon position,
#' stop-codon position, and transcript length (median or mean across transcripts), then
#' allocate the [0,1] meta-axis proportionally to 5'UTR, CDS, and 3'UTR.
#'
#' start_codon_pos ~ utr5_len
#' stop_codon_pos  ~ utr5_len + cds_len
#'
#' @noRd
.compute_metagene_splits <- function(tx_metrics,
                                     stat = c("median", "mean"),
                                     require_cds = TRUE) {
  stat <- match.arg(stat)

  default_splits <- c(fiveUTR = 1/3, coding = 1/3, threeUTR = 1/3)
  default_breaks <- c(
    fiveUTR_end = unname(default_splits["fiveUTR"]),
    coding_end  = unname(default_splits["fiveUTR"] + default_splits["coding"])
  )

  if (is.null(tx_metrics) || nrow(tx_metrics) == 0L) {
    return(list(
      splits = default_splits,
      breaks = default_breaks,
      medians = list(
        start_codon_pos = NA_real_,
        stop_codon_pos  = NA_real_,
        tx_len          = NA_real_,
        stat = stat,
        require_cds = require_cds,
        n_transcripts = 0L
      )
    ))
  }

  df <- as.data.frame(tx_metrics, stringsAsFactors = FALSE)

  needed <- c("tx_len", "utr5_len", "cds_len")
  missing <- setdiff(needed, names(df))
  if (length(missing)) {
    warning(
      "tx_metrics missing columns needed for metagene splits: ",
      paste(missing, collapse = ", "),
      ". Falling back to equal splits."
    )
    return(list(
      splits = default_splits,
      breaks = default_breaks,
      medians = list(
        start_codon_pos = NA_real_,
        stop_codon_pos  = NA_real_,
        tx_len          = NA_real_,
        stat = stat,
        require_cds = require_cds,
        n_transcripts = 0L
      )
    ))
  }

  utr5   <- as.numeric(df$utr5_len)
  cds    <- as.numeric(df$cds_len)
  tx_len <- as.numeric(df$tx_len)

  start_pos <- utr5
  stop_pos  <- utr5 + cds

  ok <- is.finite(start_pos) & is.finite(stop_pos) & is.finite(tx_len) & (tx_len > 0)
  if (require_cds) ok <- ok & is.finite(cds) & (cds > 0)

  # Basic sanity: positions should be ordered and within transcript length
  ok <- ok & (stop_pos >= start_pos) & (stop_pos <= tx_len) & (start_pos <= tx_len)

  n_ok <- sum(ok)
  if (n_ok == 0L) {
    warning(
      "No transcripts had finite/valid CDS/UTR lengths for metagene split estimation. ",
      "Falling back to equal splits."
    )
    return(list(
      splits = default_splits,
      breaks = default_breaks,
      medians = list(
        start_codon_pos = NA_real_,
        stop_codon_pos  = NA_real_,
        tx_len          = NA_real_,
        stat = stat,
        require_cds = require_cds,
        n_transcripts = 0L
      )
    ))
  }

  fun <- if (stat == "median") stats::median else base::mean

  start_est <- fun(start_pos[ok], na.rm = TRUE)
  stop_est  <- fun(stop_pos[ok],  na.rm = TRUE)
  len_est   <- fun(tx_len[ok],    na.rm = TRUE)

  L5   <- max(0, start_est)
  LCDS <- max(0, stop_est - start_est)
  L3   <- max(0, len_est - stop_est)

  total <- L5 + LCDS + L3

  if (!is.finite(total) || total <= 0) {
    warning("Metagene split estimation produced non-finite/zero total length. Falling back to equal splits.")
    splits <- default_splits
  } else {
    splits <- c(
      fiveUTR  = L5 / total,
      coding   = LCDS / total,
      threeUTR = L3 / total
    )
    splits[!is.finite(splits) | splits < 0] <- 0
    if (sum(splits) <= 0) splits <- default_splits
    splits <- splits / sum(splits)
  }

  breaks <- c(
    fiveUTR_end = unname(splits["fiveUTR"]),
    coding_end  = unname(splits["fiveUTR"] + splits["coding"])
  )

  list(
    splits = splits,
    breaks = breaks,
    medians = list(
      start_codon_pos = as.numeric(start_est),
      stop_codon_pos  = as.numeric(stop_est),
      tx_len          = as.numeric(len_est),
      stat = stat,
      require_cds = require_cds,
      n_transcripts = as.integer(n_ok)
    )
  )
}

#' Build transcript resources for fast site annotation
#'
#' Precomputes transcript-level metrics and transcript feature structures (CDS/UTR/introns)
#' from a \code{TxDb}. The returned list can be passed to annotation/geometry helpers
#' and reused across many GRanges objects to avoid repeated expensive extraction steps.
#'
#' @param txdb A \code{TxDb} object.
#' @param include_introns Logical; if TRUE, include intron GRangesLists and intron metrics.
#' @param metagene_split_stat Character; "median" (default) or "mean". Controls how the
#'   typical transcript geometry is estimated when computing metagene region splits.
#' @param metagene_require_cds Logical; if TRUE (default), compute splits using only
#'   transcripts with cds_len > 0 and consistent coordinates.
#'
#' @returns A named \code{list} containing transcript metrics and feature GRangesLists,
#'   plus length vectors and CDS anchors. Also includes metagene axis info:
#'   \code{metagene_splits}, \code{metagene_breaks}, and \code{metagene_medians}.
#'
#' @export
build_tx_resources <- function(txdb,
                               include_introns = TRUE,
                               metagene_split_stat = c("median", "mean"),
                               metagene_require_cds = TRUE) {
  metagene_split_stat <- match.arg(metagene_split_stat)

  # Transcript strand map (names are transcript names when use.names=TRUE)
  tx_gr <- GenomicFeatures::transcripts(txdb, use.names = TRUE)
  tx_strand <- as.character(strand(tx_gr))
  names(tx_strand) <- names(tx_gr)

  # Feature GRangesLists
  cds_by_tx     <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  five_by_tx    <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  three_by_tx   <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
  introns_by_tx <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)

  tx_metrics <- compute_tx_metrics(txdb, include_introns = include_introns, introns_by_tx = introns_by_tx)

  # Named length vectors
  cds_len   <- setNames(as.numeric(tx_metrics$cds_len),  tx_metrics$tx_name)
  five_len  <- setNames(as.numeric(tx_metrics$utr5_len), tx_metrics$tx_name)
  three_len <- setNames(as.numeric(tx_metrics$utr3_len), tx_metrics$tx_name)
  intr_len  <- if ("intron_len" %in% names(tx_metrics)) {
    setNames(as.numeric(tx_metrics$intron_len), tx_metrics$tx_name)
  } else {
    setNames(rep(NA_real_, nrow(tx_metrics)), tx_metrics$tx_name)
  }

  tx_names <- names(cds_by_tx)
  cds_start <- setNames(rep(NA_integer_, length(tx_names)), tx_names)
  cds_stop  <- setNames(rep(NA_integer_, length(tx_names)), tx_names)

  n_cds <- S4Vectors::elementNROWS(cds_by_tx)
  if (sum(n_cds) > 0L) {
    u <- unlist(cds_by_tx, use.names = FALSE)
    offsets <- c(0L, cumsum(n_cds))[seq_along(n_cds)]
    first_idx <- offsets + 1L
    last_idx  <- offsets + n_cds

    ok <- n_cds > 0L
    first_seg <- u[first_idx[ok]]
    last_seg  <- u[last_idx[ok]]

    st <- tx_strand[tx_names[ok]]
    # fall back to segment strand if tx_strand missing
    seg_st <- as.character(strand(first_seg))
    st[is.na(st) | st == "*"] <- seg_st[is.na(st) | st == "*"]

    plus  <- st != "-"
    minus <- !plus

    ok_names <- tx_names[ok]
    cds_start[ok_names[plus]]  <- start(first_seg[plus])
    cds_stop[ok_names[plus]]   <- end(last_seg[plus])

    cds_start[ok_names[minus]] <- end(first_seg[minus])
    cds_stop[ok_names[minus]]  <- start(last_seg[minus])
  }

  # Metagene splits (paper-like): estimate typical start/stop/len and allocate axis widths
  mg <- .compute_metagene_splits(
    tx_metrics,
    stat = metagene_split_stat,
    require_cds = metagene_require_cds
  )

  list(
    tx_metrics    = tx_metrics,
    tx_strand     = tx_strand,
    cds_by_tx     = cds_by_tx,
    five_by_tx    = five_by_tx,
    three_by_tx   = three_by_tx,
    introns_by_tx = introns_by_tx,
    cds_len       = cds_len,
    five_len      = five_len,
    three_len     = three_len,
    intr_len      = intr_len,
    cds_start     = cds_start,
    cds_stop      = cds_stop,
    metagene_splits  = mg$splits,   # named c(fiveUTR, coding, threeUTR), sums to 1
    metagene_breaks  = mg$breaks,   # c(fiveUTR_end, coding_end) on [0,1]
    metagene_medians = mg$medians   # audit trail
  )
}

# ----------------------------
# Location mapping helpers
# ----------------------------

#' Internal helper
#' @noRd
.location_priority <- function(loc) {
  # Smaller = higher priority when ties occur (e.g. promoter + intergenic)
  loc <- as.character(loc)
  pri <- rep.int(99L, length(loc))
  pri[loc == "coding"]     <- 1L
  pri[loc == "fiveUTR"]    <- 2L
  pri[loc == "threeUTR"]   <- 3L
  pri[loc == "intron"]     <- 4L
  pri[loc == "spliceSite"] <- 5L
  pri[loc == "promoter"]   <- 6L
  pri[loc == "intergenic"] <- 7L
  pri[loc == "unannotated"] <- 100L
  pri
}

#' Internal helper
#' @noRd
.map_location_to_feature <- function(location) {
  loc <- as.character(location)
  out <- rep.int("other", length(loc))

  out[loc == "coding"]     <- "CDS"
  out[loc == "fiveUTR"]    <- "5UTR"
  out[loc == "threeUTR"]   <- "3UTR"
  out[loc == "intron"]     <- "intron"
  out[loc == "intergenic"] <- "intergenic"
  out[loc == "spliceSite"] <- "spliceSite"
  out[loc == "promoter"]   <- "promoter"
  out[loc == "unannotated"] <- "unannotated"

  out
}

# ----------------------------
# Optional gene symbol mapping helper
# ----------------------------

#' Internal helper
#' @noRd
.add_gene_symbols <- function(gr,
                              orgdb = NULL,
                              gene_id_col = "gene_id",
                              keytype = NULL,
                              symbol_col = "SYMBOL",
                              name_col = NULL,
                              out_symbol_col = "gene_symbol",
                              out_name_col = "gene_name") {
  if (is.null(orgdb)) return(gr)

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("orgdb provided but package 'AnnotationDbi' is not installed/available.")
  }

  if (!(gene_id_col %in% colnames(S4Vectors::mcols(gr)))) return(gr)

  # Always create output columns with the correct length (including length-0 GRanges)
  S4Vectors::mcols(gr)[[out_symbol_col]] <- rep(NA_character_, length(gr))
  if (!is.null(name_col)) {
    S4Vectors::mcols(gr)[[out_name_col]] <- rep(NA_character_, length(gr))
  }

  if (length(gr) == 0L) return(gr)

  gids <- as.character(S4Vectors::mcols(gr)[[gene_id_col]])
  keys <- unique(gids)
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (!length(keys)) return(gr)

  if (is.null(keytype)) {
    kts <- AnnotationDbi::keytypes(orgdb)
    keytype <- if ("ENTREZID" %in% kts) "ENTREZID" else kts[1]
  }

  cols <- unique(c(symbol_col, if (!is.null(name_col)) name_col))
  sel <- AnnotationDbi::select(orgdb, keys = keys, columns = cols, keytype = keytype)

  if (is.null(sel) || nrow(sel) == 0L) return(gr)

  # Guardrails: coerce key column to character
  sel[[keytype]] <- as.character(sel[[keytype]])
  sel <- sel[!is.na(sel[[keytype]]) & nzchar(sel[[keytype]]), , drop = FALSE]
  if (nrow(sel) == 0L) return(gr)

  pick_first <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[1] else NA_character_
  }

  sym_map <- tapply(sel[[symbol_col]], sel[[keytype]], pick_first)
  S4Vectors::mcols(gr)[[out_symbol_col]] <- unname(as.character(sym_map[gids]))

  if (!is.null(name_col)) {
    nm_map <- tapply(sel[[name_col]], sel[[keytype]], pick_first)
    S4Vectors::mcols(gr)[[out_name_col]] <- unname(as.character(nm_map[gids]))
  }

  gr
}

# ----------------------------
# Feature geometry and metagene coords (per-row / per-transcript)
# ----------------------------

#' Compute within-feature geometry and metagene coordinates
#'
#' Adds feature-level geometry columns and metagene coordinates.
#'
#' Metagene coordinates added:
#'   \code{metagene_prop}   in [0,3] with equal region widths
#'   \code{metagene_split}  in [0,1] with median/mean-derived region widths
#'   \code{metagene_split3} in [0,3] scaled from \code{metagene_split}
#'
#' Isoform weighting:
#' \code{isoform_weight} is set to 1/k within each \code{record_id} group (if present),
#' otherwise within each \code{names(gr)}.
#'
#' @param gr A \code{GRanges} with \code{location} and \code{tx_name} in \code{mcols(gr)}.
#' @param resources Output from \code{\link{build_tx_resources}}.
#'
#' @returns The input \code{GRanges} with additional geometry columns added.
#'
#' @export
add_feature_geometry <- function(gr, resources) {
  if (!("location" %in% colnames(mcols(gr)))) {
    stop("`gr` must have a 'location' column (run annotate_sites or expand_sites_by_transcript first).")
  }
  if (!("tx_name" %in% colnames(mcols(gr)))) {
    stop("`gr` must have a 'tx_name' column (run annotate_sites or expand_sites_by_transcript first).")
  }

  mc <- mcols(gr)

  if ("feature_prop" %in% names(mc) && !is.numeric(mc[["feature_prop"]])) {
    warning("Existing mcols(gr)$feature_prop is not numeric; overwriting with numeric feature geometry.")
  }
  if ("metagene_prop" %in% names(mc) && !is.numeric(mc[["metagene_prop"]])) {
    warning("Existing mcols(gr)$metagene_prop is not numeric; overwriting with numeric metagene geometry.")
  }
  if ("metagene_split" %in% names(mc) && !is.numeric(mc[["metagene_split"]])) {
    warning("Existing mcols(gr)$metagene_split is not numeric; overwriting with numeric metagene geometry.")
  }
  if ("isoform_weight" %in% names(mc) && !is.numeric(mc[["isoform_weight"]])) {
    warning("Existing mcols(gr)$isoform_weight is not numeric; overwriting with numeric weights.")
  }

  mc[["feature_len"]]             <- rep(NA_real_,    length(gr))
  mc[["feature_prop"]]            <- rep(NA_real_,    length(gr))
  mc[["metagene_prop"]]           <- rep(NA_real_,    length(gr))
  mc[["metagene_split"]]          <- rep(NA_real_,    length(gr))
  mc[["metagene_split3"]]         <- rep(NA_real_,    length(gr))
  mc[["isoform_weight"]]          <- rep(1.0,         length(gr))
  mc[["feature_width"]]           <- rep(NA_integer_, length(gr))
  mc[["segment_rank"]]            <- rep(NA_integer_, length(gr))
  mc[["dist_from_feature_start"]] <- rep(NA_integer_, length(gr))
  mc[["dist_from_feature_end"]]   <- rep(NA_integer_, length(gr))
  mc[["start_dist"]]              <- rep(NA_integer_, length(gr))
  mc[["stop_dist"]]               <- rep(NA_integer_, length(gr))

  mcols(gr) <- mc

  .fill_one <- function(gr, idx, grl, flen_vec) {
    if (length(idx) == 0L) return(gr)

    sites <- gr[idx]
    txs <- unique(as.character(sites$tx_name))
    txs <- txs[!is.na(txs) & nzchar(txs)]
    if (!length(txs)) return(gr)

    txs2 <- intersect(txs, names(grl))
    if (!length(txs2)) return(gr)

    grl_sub <- grl[txs2]
    if (!length(grl_sub)) return(gr)

    gr$feature_len[idx] <- flen_vec[as.character(sites$tx_name)]

    feat_gr <- unlist(grl_sub, use.names = FALSE)
    nseg <- S4Vectors::elementNROWS(grl_sub)
    feat_gr$tx_name <- rep(names(grl_sub), nseg)

    if ("exon_rank" %in% colnames(mcols(feat_gr))) {
      feat_gr$segment_rank <- as.integer(feat_gr$exon_rank)
    } else {
      feat_gr$segment_rank <- sequence(nseg)
    }

    hits <- GenomicRanges::findOverlaps(sites, feat_gr, ignore.strand = FALSE)
    if (length(hits) > 0L) {
      qh <- queryHits(hits)
      sh <- subjectHits(hits)

      keep <- as.character(sites$tx_name[qh]) == as.character(feat_gr$tx_name[sh])
      qh <- qh[keep]; sh <- sh[keep]

      if (length(qh) > 0L) {
        o <- order(qh)
        qh <- qh[o]; sh <- sh[o]
        keep_first <- !duplicated(qh)
        qh <- qh[keep_first]; sh <- sh[keep_first]

        seg <- feat_gr[sh]
        site_sub_idx <- idx[qh]

        gr$feature_width[site_sub_idx] <- width(seg)
        gr$segment_rank[site_sub_idx]  <- seg$segment_rank

        pos <- start(gr[site_sub_idx])
        seg_start <- start(seg)
        seg_end   <- end(seg)
        st <- as.character(strand(gr[site_sub_idx]))

        plus <- st != "-"
        dist_start <- integer(length(pos))
        dist_end   <- integer(length(pos))

        dist_start[plus] <- pmax(0L, pos[plus] - seg_start[plus])
        dist_end[plus]   <- pmax(0L, seg_end[plus] - pos[plus])

        minus <- !plus
        dist_start[minus] <- pmax(0L, seg_end[minus] - pos[minus])
        dist_end[minus]   <- pmax(0L, pos[minus] - seg_start[minus])

        gr$dist_from_feature_start[site_sub_idx] <- dist_start
        gr$dist_from_feature_end[site_sub_idx]   <- dist_end
      }
    }

    mapped <- suppressWarnings(
      GenomicFeatures::mapToTranscripts(sites, transcripts = grl_sub, ignore.strand = FALSE)
    )

    if (length(mapped) > 0L) {
      q_idx <- mapped$xHits
      mapped_tx <- as.character(seqnames(mapped))
      q_tx <- as.character(sites$tx_name[q_idx])
      keep <- !is.na(q_tx) & nzchar(q_tx) & (mapped_tx == q_tx)
      mapped <- mapped[keep]
    }

    if (length(mapped) > 0L) {
      mapped_tx <- as.character(seqnames(mapped))
      fl <- as.numeric(flen_vec[mapped_tx])

      prop <- start(mapped) / fl
      ok <- is.finite(prop) & prop >= 0 & prop <= 1
      prop[!ok] <- NA_real_

      prop_by_hit <- tapply(prop, mapped$xHits, function(z) {
        if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
      })

      hit_idx <- as.integer(names(prop_by_hit))
      hit_idx <- hit_idx[!is.na(hit_idx) & hit_idx >= 1L & hit_idx <= length(idx)]

      if (length(hit_idx)) {
        gr$feature_prop[idx[hit_idx]] <- as.numeric(prop_by_hit[as.character(hit_idx)])
      }
    }

    gr
  }

  gr <- .fill_one(gr, which(gr$location == "fiveUTR"),  resources$five_by_tx,    resources$five_len)
  gr <- .fill_one(gr, which(gr$location == "coding"),   resources$cds_by_tx,     resources$cds_len)
  gr <- .fill_one(gr, which(gr$location == "threeUTR"), resources$three_by_tx,   resources$three_len)
  gr <- .fill_one(gr, which(gr$location == "intron"),   resources$introns_by_tx, resources$intr_len)

  fp <- gr$feature_prop
  i5 <- which(gr$location == "fiveUTR"  & is.finite(fp))
  ic <- which(gr$location == "coding"   & is.finite(fp))
  i3 <- which(gr$location == "threeUTR" & is.finite(fp))

  gr$metagene_prop[i5] <- fp[i5]
  gr$metagene_prop[ic] <- fp[ic] + 1
  gr$metagene_prop[i3] <- fp[i3] + 2

  splits <- resources$metagene_splits
  breaks <- resources$metagene_breaks

  if (is.null(splits) || length(splits) < 3L || any(!c("fiveUTR", "coding", "threeUTR") %in% names(splits))) {
    mg <- .compute_metagene_splits(resources$tx_metrics, stat = "median", require_cds = TRUE)
    splits <- mg$splits
    breaks <- mg$breaks
  } else {
    splits <- splits[c("fiveUTR", "coding", "threeUTR")]
    if (is.null(breaks) || length(breaks) < 2L) {
      breaks <- c(
        fiveUTR_end = unname(splits["fiveUTR"]),
        coding_end  = unname(splits["fiveUTR"] + splits["coding"])
      )
    }
  }

  if (length(i5)) gr$metagene_split[i5] <- fp[i5] * unname(splits["fiveUTR"])
  if (length(ic)) gr$metagene_split[ic] <- unname(splits["fiveUTR"]) + fp[ic] * unname(splits["coding"])
  if (length(i3)) gr$metagene_split[i3] <- unname(splits["fiveUTR"] + splits["coding"]) + fp[i3] * unname(splits["threeUTR"])

  ok_ms <- is.finite(gr$metagene_split)
  gr$metagene_split3[ok_ms] <- gr$metagene_split[ok_ms] * 3

  rid <- if ("record_id" %in% colnames(mcols(gr))) as.character(gr$record_id) else names(gr)
  w <- rep(1.0, length(gr))
  ok_r <- !is.na(rid) & nzchar(rid)
  if (any(ok_r)) {
    k <- as.numeric(ave(rid[ok_r], rid[ok_r], FUN = length))
    k[!is.finite(k) | k <= 0] <- NA_real_
    ww <- 1 / k
    ww[!is.finite(ww)] <- 1.0
    w[ok_r] <- ww
  }
  gr$isoform_weight <- w

  tx <- as.character(gr$tx_name)
  gr$start_dist <- abs(start(gr) - resources$cds_start[tx])
  gr$stop_dist  <- abs(start(gr) - resources$cds_stop[tx])

  meta <- S4Vectors::metadata(gr)
  meta$metagene_splits  <- splits
  meta$metagene_breaks  <- breaks
  meta$metagene_breaks3 <- breaks * 3
  S4Vectors::metadata(gr) <- meta

  gr
}

# ----------------------------
# Expand sites to isoforms (optional, heavier object)
# ----------------------------

#' Expand sites to all overlapping transcripts (core regions) for isoform-weighted metagenes
#'
#' Creates a multi-isoform representation of point sites by expanding each input site
#' into one row per overlapping transcript feature (5'UTR/CDS/3'UTR by default).
#'
#' @param gr A \code{GRanges} of point sites (width 1 recommended).
#' @param resources Output from \code{\link{build_tx_resources}}.
#' @param regions Regions to include (default c("fiveUTR","coding","threeUTR")).
#' @param ignore.strand Passed to \code{findOverlaps} (default FALSE).
#' @param record_id_col Optional input column to use as the record id; otherwise uses names(gr).
#' @param out_record_id_col Output column name for record ids (default "record_id").
#' @param make_unique_names If TRUE, names are set to record_id__tx_name.
#' @param drop_unannotated If TRUE, drops sites that do not overlap requested regions.
#'
#' @returns Expanded \code{GRanges} with \code{record_id}, \code{tx_name}, and \code{location}.
#'
#' @export
expand_sites_by_transcript <- function(gr,
                                       resources,
                                       regions = c("fiveUTR", "coding", "threeUTR"),
                                       ignore.strand = FALSE,
                                       record_id_col = NULL,
                                       out_record_id_col = "record_id",
                                       make_unique_names = TRUE,
                                       drop_unannotated = TRUE) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  if (any(GenomicRanges::width(gr) != 1L)) {
    warning("Some ranges have width != 1. Resizing to width=1 (fix='center').")
    gr <- IRanges::resize(gr, width = 1L, fix = "center")
  }
  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_along(gr))
  }

  rid <- NULL
  if (!is.null(record_id_col) && (record_id_col %in% colnames(S4Vectors::mcols(gr)))) {
    rid <- as.character(S4Vectors::mcols(gr)[[record_id_col]])
  } else {
    rid <- names(gr)
  }
  bad <- is.na(rid) | !nzchar(rid)
  if (any(bad)) rid[bad] <- names(gr)[bad]

  .get_grl <- function(region) {
    if (region == "fiveUTR")  return(resources$five_by_tx)
    if (region == "coding")   return(resources$cds_by_tx)
    if (region == "threeUTR") return(resources$three_by_tx)
    if (region == "intron")   return(resources$introns_by_tx)
    stop("Unknown region: ", region)
  }

  .collect_hits <- function(gr, grl, location, ignore.strand) {
    if (is.null(grl) || length(grl) == 0L) return(NULL)

    nseg <- S4Vectors::elementNROWS(grl)
    if (sum(nseg) == 0L) return(NULL)

    feat_gr <- unlist(grl, use.names = FALSE)
    feat_gr$tx_name <- rep(names(grl), nseg)

    hits <- GenomicRanges::findOverlaps(gr, feat_gr, ignore.strand = ignore.strand)
    if (length(hits) == 0L) return(NULL)

    q <- queryHits(hits)
    tx <- as.character(feat_gr$tx_name[subjectHits(hits)])

    key <- paste(q, tx, location, sep = "\t")
    keep <- !duplicated(key)
    q <- q[keep]
    tx <- tx[keep]

    data.frame(
      q = q,
      tx_name = tx,
      location = rep.int(location, length(q)),
      priority = .location_priority(rep.int(location, length(q))),
      stringsAsFactors = FALSE
    )
  }

  regions <- unique(as.character(regions))
  if (!length(regions)) stop("`regions` must have at least one entry.")

  all_hits <- NULL
  for (reg in regions) {
    grl <- .get_grl(reg)
    h <- .collect_hits(gr, grl, location = reg, ignore.strand = ignore.strand)
    if (!is.null(h) && nrow(h) > 0L) {
      all_hits <- if (is.null(all_hits)) h else rbind(all_hits, h)
    }
  }

  if (is.null(all_hits) || nrow(all_hits) == 0L) {
    if (drop_unannotated) {
      return(gr[0])
    } else {
      out <- gr
      S4Vectors::mcols(out)[[out_record_id_col]] <- rid
      S4Vectors::mcols(out)[["tx_name"]] <- NA_character_
      S4Vectors::mcols(out)[["location"]] <- NA_character_
      return(out)
    }
  }

  all_hits <- all_hits[order(all_hits$q, all_hits$tx_name, all_hits$priority), , drop = FALSE]
  keep_best <- !duplicated(paste(all_hits$q, all_hits$tx_name, sep = "\t"))
  all_hits <- all_hits[keep_best, , drop = FALSE]

  out <- gr[all_hits$q]
  S4Vectors::mcols(out)[[out_record_id_col]] <- rid[all_hits$q]
  S4Vectors::mcols(out)[["tx_name"]] <- all_hits$tx_name
  S4Vectors::mcols(out)[["location"]] <- all_hits$location

  if (make_unique_names) {
    names(out) <- paste0(S4Vectors::mcols(out)[[out_record_id_col]], "__", S4Vectors::mcols(out)[["tx_name"]])
  }

  out
}

# ----------------------------
# NEW: Add isoform-averaged metagene split to original sites (no expanded output)
# ----------------------------

#' Add an isoform-averaged metagene split coordinate to an annotated GRanges
#'
#' This keeps the input object at the same length (one row per original site) and adds
#' \code{metagene_split} (in [0,1]) plus \code{metagene_split3} (in [0,3]).
#'
#' Internally it (i) expands *only a minimal copy* of the sites (no metadata copied),
#' (ii) computes per-isoform metagene positions, and (iii) collapses back to one value
#' per original site using the 1/k isoform weighting (equivalently: mean across isoform mappings).
#'
#' @param gr A \code{GRanges}, typically the output from your main annotation function.
#' @param resources Output from \code{\link{build_tx_resources}}.
#'
#' @returns The input \code{GRanges} with \code{metagene_split} and \code{metagene_split3} added.
#'
#' @export
add_metagene_split <- function(gr, resources) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  gr <- .validate_sites_gr(gr, label_col = NULL)

  # Initialize output columns (always correct length/type)
  mc <- S4Vectors::mcols(gr)
  mc[["metagene_split"]]  <- rep(NA_real_, length(gr))
  mc[["metagene_split3"]] <- rep(NA_real_, length(gr))
  S4Vectors::mcols(gr) <- mc

  if (length(gr) == 0L) return(gr)

  splits <- resources$metagene_splits
  breaks <- resources$metagene_breaks
  if (is.null(splits) || length(splits) < 3L || any(!c("fiveUTR", "coding", "threeUTR") %in% names(splits))) {
    mg <- .compute_metagene_splits(resources$tx_metrics, stat = "median", require_cds = TRUE)
    splits <- mg$splits
    breaks <- mg$breaks
  } else {
    splits <- splits[c("fiveUTR", "coding", "threeUTR")]
    if (is.null(breaks) || length(breaks) < 2L) {
      breaks <- c(
        fiveUTR_end = unname(splits["fiveUTR"]),
        coding_end  = unname(splits["fiveUTR"] + splits["coding"])
      )
    }
  }

  # Minimal sites object to avoid copying metadata into expanded isoforms
  gr0 <- gr
  #S4Vectors::mcols(gr0) <- S4Vectors::DataFrame()

  gr_exp <- expand_sites_by_transcript(
    gr0, resources,
    regions = c("fiveUTR", "coding", "threeUTR"),
    ignore.strand = FALSE,
    record_id_col = NULL,
    out_record_id_col = "record_id",
    make_unique_names = FALSE,
    drop_unannotated = TRUE
  )

  if (length(gr_exp) == 0L) {
    meta <- S4Vectors::metadata(gr)
    meta$metagene_splits  <- splits
    meta$metagene_breaks  <- breaks
    meta$metagene_breaks3 <- breaks * 3
    S4Vectors::metadata(gr) <- meta
    return(gr)
  }

  # Compute within-feature proportions (minimal: no segment rank/edge distances)
  feature_prop <- rep(NA_real_, length(gr_exp))

  .prop_for_region <- function(region, grl, flen_vec) {
    idx <- which(gr_exp$location == region)
    if (!length(idx)) return(invisible(NULL))

    sites <- gr_exp[idx]
    txs <- unique(as.character(sites$tx_name))
    txs <- txs[!is.na(txs) & nzchar(txs)]
    if (!length(txs)) return(invisible(NULL))

    txs2 <- intersect(txs, names(grl))
    if (!length(txs2)) return(invisible(NULL))

    grl_sub <- grl[txs2]
    if (!length(grl_sub)) return(invisible(NULL))

    mapped <- suppressWarnings(
      GenomicFeatures::mapToTranscripts(sites, transcripts = grl_sub, ignore.strand = FALSE)
    )

    if (length(mapped) > 0L) {
      q_idx <- mapped$xHits
      mapped_tx <- as.character(seqnames(mapped))
      q_tx <- as.character(sites$tx_name[q_idx])
      keep <- !is.na(q_tx) & nzchar(q_tx) & (mapped_tx == q_tx)
      mapped <- mapped[keep]
    }

    if (length(mapped) > 0L) {
      fl <- as.numeric(flen_vec[as.character(seqnames(mapped))])
      prop <- start(mapped) / fl

      ok <- is.finite(prop) & prop >= 0 & prop <= 1
      prop[!ok] <- NA_real_

      prop_by_hit <- tapply(prop, mapped$xHits, function(z) {
        if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
      })

      hit_idx <- as.integer(names(prop_by_hit))
      hit_idx <- hit_idx[!is.na(hit_idx) & hit_idx >= 1L & hit_idx <= length(idx)]
      if (length(hit_idx)) {
        feature_prop[idx[hit_idx]] <<- as.numeric(prop_by_hit[as.character(hit_idx)])
      }
    }

    invisible(NULL)
  }

  .prop_for_region("fiveUTR",  resources$five_by_tx,  resources$five_len)
  .prop_for_region("coding",   resources$cds_by_tx,   resources$cds_len)
  .prop_for_region("threeUTR", resources$three_by_tx, resources$three_len)

  # Map feature_prop -> metagene_split (per isoform mapping row)
  ms <- rep(NA_real_, length(gr_exp))
  okp <- is.finite(feature_prop)

  i5 <- which(gr_exp$location == "fiveUTR"  & okp)
  ic <- which(gr_exp$location == "coding"   & okp)
  i3 <- which(gr_exp$location == "threeUTR" & okp)

  if (length(i5)) ms[i5] <- feature_prop[i5] * unname(splits["fiveUTR"])
  if (length(ic)) ms[ic] <- unname(splits["fiveUTR"]) + feature_prop[ic] * unname(splits["coding"])
  if (length(i3)) ms[i3] <- unname(splits["fiveUTR"] + splits["coding"]) + feature_prop[i3] * unname(splits["threeUTR"])

  # Isoform weights 1/k per record_id among expanded rows
  rid <- as.character(gr_exp$record_id)
  bad <- is.na(rid) | !nzchar(rid)
  if (any(bad)) rid[bad] <- names(gr_exp)[bad]

  k <- as.numeric(ave(rid, rid, FUN = length))
  k[!is.finite(k) | k <= 0] <- NA_real_
  w <- 1 / k

  # Collapse back to one value per original site (weighted mean, robust to NA)
  ok <- is.finite(ms) & is.finite(w) & (w > 0) & !is.na(rid) & nzchar(rid)
  wm <- rep(NA_real_, length(gr))
  names(wm) <- names(gr)

  if (any(ok)) {
    fac <- factor(rid[ok])
    num <- rowsum(ms[ok] * w[ok], fac, reorder = FALSE)[, 1]
    den <- rowsum(w[ok], fac, reorder = FALSE)[, 1]
    val <- num / den
    names(val) <- rownames(rowsum(w[ok], fac, reorder = FALSE))

    hit <- intersect(names(gr), names(val))
    wm[hit] <- val[hit]
  }

  S4Vectors::mcols(gr)[["metagene_split"]]  <- wm
  S4Vectors::mcols(gr)[["metagene_split3"]] <- wm * 3

  meta <- S4Vectors::metadata(gr)
  meta$metagene_splits  <- splits
  meta$metagene_breaks  <- breaks
  meta$metagene_breaks3 <- breaks * 3
  S4Vectors::metadata(gr) <- meta

  gr
}

# ----------------------------
# SIMPLE plot helper (base R)
# ----------------------------

#' Simple metagene density plot (base R)
#'
#' Uses \code{stats::density()} on \code{metagene_split3} (preferred), otherwise falls back
#' to \code{metagene_split} (scaled to 0–3) or \code{metagene_prop}. Adds 5'UTR/CDS/3'UTR
#' boundaries using metadata breakpoints when available.
#'
#' @param gr A \code{GRanges} with metagene columns.
#' @param x_col Optional column name to use (default auto).
#' @param main Plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param ... Passed to \code{graphics::plot()}.
#'
#' @returns The \code{density} object (invisibly).
#' @export
plot_metagene_density <- function(gr,
                                  x_col = NULL,
                                  main = "Metagene density",
                                  xlab = "Metagene position",
                                  ylab = "Density",
                                  ...) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  mc <- S4Vectors::mcols(gr)

  if (is.null(x_col)) {
    if ("metagene_split3" %in% colnames(mc)) {
      x_col <- "metagene_split3"
    } else if ("metagene_prop" %in% colnames(mc)) {
      x_col <- "metagene_prop"
    } else if ("metagene_split" %in% colnames(mc)) {
      x_col <- "metagene_split"
    } else {
      stop("No metagene column found. Expected one of: metagene_split3, metagene_prop, metagene_split.")
    }
  } else {
    if (!(x_col %in% colnames(mc))) stop("x_col not found in mcols(gr): ", x_col)
  }

  x <- suppressWarnings(as.numeric(mc[[x_col]]))
  x <- x[is.finite(x)]
  if (!length(x)) stop("No finite values in ", x_col)

  md <- S4Vectors::metadata(gr)

  # Always plot on a 0–3 axis for readability
  if (identical(x_col, "metagene_split")) {
    x3 <- x * 3
  } else {
    x3 <- x
  }

  breaks3 <- md$metagene_breaks3
  if (is.null(breaks3) || length(breaks3) < 2L || any(!is.finite(as.numeric(breaks3[1:2])))) {
    breaks3 <- c(1, 2)
  } else {
    breaks3 <- as.numeric(breaks3[1:2])
  }

  d <- stats::density(x3, from = 0, to = 3)

  graphics::plot(d, type = "n", main = main, xlab = paste0(xlab, " (0–3)"), ylab = ylab, xlim = c(0, 3), ...)
  usr <- graphics::par("usr")

  cols <- grDevices::adjustcolor(c("grey70", "grey60", "grey70"), alpha.f = 0.15)
  graphics::rect(
    xleft = c(0, breaks3[1], breaks3[2]),
    xright = c(breaks3[1], breaks3[2], 3),
    ybottom = usr[3], ytop = usr[4],
    col = cols, border = NA
  )

  graphics::lines(d, lwd = 1)
  graphics::abline(v = breaks3, lty = 2, lwd = 0.8)

  mids <- c((0 + breaks3[1]) / 2, (breaks3[1] + breaks3[2]) / 2, (breaks3[2] + 3) / 2)
  graphics::text(x = mids, y = usr[4] - 0.05 * (usr[4] - usr[3]),
                 labels = c("5'UTR", "CDS", "3'UTR"), cex = 0.9)

  invisible(d)
}

# ----------------------------
# Sequence context helper
# ----------------------------

#' Add sequence k-mer context around each site
#'
#' Extracts a centered k-mer from a reference genome and stores it as a metadata column.
#'
#' @param gr A \code{GRanges} of single-nucleotide sites.
#' @param genome A BSgenome object (or other supported genome object for \code{Biostrings::getSeq}).
#' @param k Integer k-mer length (default 5).
#' @param out_col Name of the metadata column to write (default \code{"kmer"}).
#' @param seqstyle Optional seqlevel style (e.g. \code{"UCSC"}).
#' @param chrs Optional character vector of seqlevels to keep.
#'
#' @returns The input \code{GRanges} with \code{out_col} added to \code{mcols(gr)}.
#' @export
add_kmer <- function(gr, genome, k = 5L, out_col = "kmer", seqstyle = NULL, chrs = NULL) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  if (k < 1L) stop("`k` must be >= 1.")
  if ((k %% 2L) == 0L) warning("Even k: window won't be perfectly symmetric around the center base.")

  if (!is.null(seqstyle)) {
    GenomeInfoDb::seqlevelsStyle(gr) <- seqstyle
    GenomeInfoDb::seqlevelsStyle(genome) <- seqstyle
  } else {
    st <- GenomeInfoDb::seqlevelsStyle(genome)
    if (length(st) > 0L) suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- st[1], silent = TRUE))
  }

  if (is.null(chrs)) chrs <- intersect(seqlevels(gr), seqlevels(genome))
  gr <- GenomeInfoDb::keepSeqlevels(gr, chrs, pruning.mode = "coarse")

  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(genome)[chrs]

  win <- IRanges::resize(gr, width = k, fix = "center")
  win <- GenomicRanges::trim(win)

  mcols(gr)[[out_col]] <- as.character(Biostrings::getSeq(genome, win))
  gr
}
