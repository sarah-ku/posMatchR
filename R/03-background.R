#' Match positive sites to a candidate background
#'
#' Greedily pairs each positive site to at most one unique negative/background site
#' within hard strata such as transcript region and optional exact k-mer. Within each
#' stratum, the primary cost is distance in a metagene coordinate; optional numeric
#' covariates can be quantile-binned and used as soft or hard secondary constraints.
#'
#' @param gr A labelled \code{GRanges} containing positive and candidate background sites.
#' @param label_col Metadata column containing 1 for positives and 0 for negatives.
#' @param location_col Region column, usually \code{"location"}.
#' @param locations Eligible regions.
#' @param return_diagnostics If TRUE, store diagnostics in \code{metadata(gr)$match_diagnostics}.
#' @param kmer_match If TRUE, require exact \code{kmer_col} agreement.
#' @param kmer_col K-mer column.
#' @param meta_col Metagene coordinate column. If \code{"metagene_split3"} is requested
#'   but absent, \code{"metagene_prop"} is used when present.
#' @param meta_tol Optional tolerance in metagene-coordinate units, not genomic bp.
#' @param enforce_meta_tol If TRUE, leave positives unmatched when no candidate lies inside tolerance.
#' @param meta_k Initial number of candidate negatives to inspect around each positive in metagene-sorted order. This is a computational candidate-pool size, not a genomic bp window.
#' @param bin_match If TRUE, use quantile-binned numeric covariates as secondary constraints.
#' @param bin_cols Candidate numeric covariates to bin; missing columns are ignored.
#' @param n_bins Number of quantile bins.
#' @param bin_within Compute quantile-bin cutoffs globally or separately within transcript region.
#' @param bin_mode \code{"soft"} adds a penalty; \code{"hard"} requires matching bins.
#' @param bin_weight Weight for soft bin mismatch. The default gives geometry-bin differences similar influence to the metagene cost.
#' @param log_bin_cols Columns log1p-transformed before binning.
#' @param seed Random seed for deterministic tie behaviour.
#'
#' @return A \code{GRanges} with canonical match columns.
#' @export
match_background <- function(gr,
                             label_col = "label",
                             location_col = "location",
                             locations = c("fiveUTR", "coding", "threeUTR"),
                             return_diagnostics = TRUE,
                             kmer_match = FALSE,
                             kmer_col = "kmer",
                             meta_col = "metagene_split3",
                             meta_tol = 0.05,
                             enforce_meta_tol = FALSE,
                             meta_k = 200L,
                             bin_match = TRUE,
                             bin_cols = c(
                               "nearest_exon_junction_dist", "start_dist_tx", "stop_dist_tx",
                               "start_dist", "stop_dist", "tx_len", "feature_width", "segment_rank",
                               "dist_from_feature_start", "dist_from_feature_end"
                             ),
                             n_bins = 5L,
                             bin_within = c("location", "global"),
                             bin_mode = c("soft", "hard"),
                             bin_weight = 1.0,
                             log_bin_cols = c(
                               "tx_len", "feature_width", "nearest_exon_junction_dist",
                               "dist_from_feature_start", "dist_from_feature_end",
                               "start_dist_tx", "stop_dist_tx", "start_dist", "stop_dist"
                             ),
                             seed = 1L) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  mc <- S4Vectors::mcols(gr)

  if (!(label_col %in% colnames(mc))) stop("Missing label_col: ", label_col)
  if (!(location_col %in% colnames(mc))) stop("Missing location_col: ", location_col)
  if (!(meta_col %in% colnames(mc))) {
    if (identical(meta_col, "metagene_split3") && "metagene_prop" %in% colnames(mc)) {
      meta_col <- "metagene_prop"
    } else {
      stop("Missing meta_col: ", meta_col)
    }
  }
  if (kmer_match && !(kmer_col %in% colnames(mc))) {
    stop("kmer_match=TRUE but missing kmer_col: ", kmer_col)
  }

  bin_within <- match.arg(bin_within)
  bin_mode <- match.arg(bin_mode)

  gr <- prepare_sites(gr, label = NULL, label_col = label_col, strip_mcols = FALSE)
  mc <- S4Vectors::mcols(gr)
  y <- .normalise_binary_label(mc[[label_col]], strict = TRUE)
  mc[[label_col]] <- y

  loc <- as.character(mc[[location_col]])
  meta <- suppressWarnings(as.numeric(mc[[meta_col]]))

  mc$loc_fiveUTR <- as.integer(loc == "fiveUTR")
  mc$loc_coding <- as.integer(loc == "coding")
  mc$loc_threeUTR <- as.integer(loc == "threeUTR")

  in_scope <- !is.na(loc) & loc %in% locations & is.finite(meta)
  pos_idx <- which(y == 1L & in_scope)
  neg_idx <- which(y == 0L & in_scope)

  if (kmer_match) {
    km <- as.character(mc[[kmer_col]])
    pos_idx <- pos_idx[!is.na(km[pos_idx]) & nzchar(km[pos_idx])]
    neg_idx <- neg_idx[!is.na(km[neg_idx]) & nzchar(km[neg_idx])]
  }

  matched_neg_for_pos <- rep(NA_character_, length(gr))
  matched_pos_for_neg <- rep(NA_character_, length(gr))
  match_meta_delta <- rep(NA_real_, length(gr))
  neg_used <- rep(FALSE, length(gr))

  quantile_bin <- function(x, n_bins) {
    x <- suppressWarnings(as.numeric(x))
    out <- rep.int(NA_integer_, length(x))
    ok <- is.finite(x)
    if (!any(ok)) return(out)
    br <- stats::quantile(x[ok], probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE, type = 7)
    br <- unique(as.numeric(br))
    if (length(br) < 2L) {
      out[ok] <- 1L
      return(out)
    }
    br[1] <- -Inf
    br[length(br)] <- Inf
    out[ok] <- as.integer(cut(x[ok], breaks = br, include.lowest = TRUE, labels = FALSE))
    out
  }

  bin_mat <- NULL
  cols_ok <- character(0)
  if (bin_match) {
    cols_ok <- intersect(bin_cols, colnames(mc))
    if (length(cols_ok)) {
      get_feat <- function(col) {
        v <- suppressWarnings(as.numeric(mc[[col]]))
        if (col %in% log_bin_cols) v <- log1p(pmax(v, 0))
        v
      }
      bin_mat <- matrix(NA_integer_, nrow = length(gr), ncol = length(cols_ok))
      colnames(bin_mat) <- cols_ok

      if (bin_within == "global") {
        for (j in seq_along(cols_ok)) bin_mat[, j] <- quantile_bin(get_feat(cols_ok[j]), n_bins)
      } else {
        for (lv in locations) {
          idx_lv <- which(loc == lv & loc %in% locations)
          if (!length(idx_lv)) next
          for (j in seq_along(cols_ok)) {
            v <- get_feat(cols_ok[j])
            b <- rep.int(NA_integer_, length(v))
            b[idx_lv] <- quantile_bin(v[idx_lv], n_bins)
            old <- bin_mat[, j]
            old[is.na(old)] <- b[is.na(old)]
            bin_mat[, j] <- old
          }
        }
      }

      for (j in seq_along(cols_ok)) mc[[paste0("bin_", cols_ok[j])]] <- bin_mat[, j]
    } else {
      bin_match <- FALSE
    }
  }

  if (length(pos_idx) && length(neg_idx)) {
    if (kmer_match) {
      km <- as.character(mc[[kmer_col]])
      key <- paste0(loc, "||", km)
    } else {
      key <- loc
    }
    pos_by_key <- split(pos_idx, key[pos_idx], drop = TRUE)
    neg_by_key <- split(neg_idx, key[neg_idx], drop = TRUE)
    keys <- intersect(names(pos_by_key), names(neg_by_key))
    keys <- keys[!is.na(keys) & nzchar(keys)]

    set.seed(seed)

    for (kkey in keys) {
      P <- pos_by_key[[kkey]]
      N <- neg_by_key[[kkey]]
      if (!length(P) || !length(N)) next

      P <- P[order(meta[P])]
      N_sorted <- N[order(meta[N])]
      N_meta <- meta[N_sorted]

      score_candidates <- function(pi, cand, x) {
        cand <- cand[!neg_used[cand]]
        if (!length(cand)) return(list(chosen = NA_integer_, n_considered = 0L))

        md <- abs(meta[cand] - x)
        if (!is.null(meta_tol) && is.finite(meta_tol) && meta_tol > 0 && enforce_meta_tol) {
          keep_tol <- md <= meta_tol
          cand <- cand[keep_tol]
          md <- md[keep_tol]
          if (!length(cand)) return(list(chosen = NA_integer_, n_considered = 0L))
        }

        meta_scale <- if (!is.null(meta_tol) && is.finite(meta_tol) && meta_tol > 0) meta_tol else 0.05
        cost_meta <- md / meta_scale
        cost <- cost_meta

        if (bin_match && !is.null(bin_mat) && ncol(bin_mat) > 0L) {
          bp <- bin_mat[pi, , drop = TRUE]
          bn <- bin_mat[cand, , drop = FALSE]

          if (bin_mode == "hard") {
            ok <- rep(TRUE, nrow(bn))
            for (jj in seq_len(ncol(bn))) {
              ok <- ok & (bn[, jj] == bp[jj] | (is.na(bn[, jj]) & is.na(bp[jj])))
            }
            cand <- cand[ok]
            md <- md[ok]
            cost_meta <- cost_meta[ok]
            bn <- bn[ok, , drop = FALSE]
            if (!length(cand)) return(list(chosen = NA_integer_, n_considered = 0L))
          }

          denom <- max(1, n_bins - 1L)
          bin_diff <- abs(sweep(bn, 2, bp, "-"))
          cost_bin <- rowMeans(bin_diff / denom, na.rm = TRUE)
          cost_bin[!is.finite(cost_bin)] <- 0
          cost <- cost_meta + bin_weight * cost_bin
        }

        ord <- order(cost, md)
        list(chosen = cand[ord[1L]], n_considered = length(cand))
      }

      for (pi in P) {
        if (!is.na(matched_neg_for_pos[pi])) next
        if (!any(!neg_used[N_sorted])) break

        x <- meta[pi]
        chosen <- NA_integer_

        # When a strict metagene tolerance is requested, score all unused
        # candidates in the stratum that fall within that tolerance. This makes
        # transcript-geometry bin penalties meaningful: a slightly more distant
        # metagene neighbour can be preferred if it is much better matched in
        # exon/transcript geometry. `meta_k` is only used as a computational
        # shortcut when the tolerance is not being enforced.
        if (!is.null(meta_tol) && is.finite(meta_tol) && meta_tol > 0 && enforce_meta_tol) {
          lo <- findInterval(x - meta_tol, N_meta) + 1L
          hi <- findInterval(x + meta_tol, N_meta)
          lo <- max(1L, lo)
          hi <- min(length(N_sorted), hi)
          if (lo <= hi) {
            cand <- N_sorted[lo:hi]
            res <- score_candidates(pi, cand, x)
            chosen <- res$chosen
          }
        } else {
          j <- findInterval(x, N_meta)
          j <- max(1L, min(j, length(N_sorted)))
          win <- max(10L, as.integer(meta_k))
          repeat {
            half <- win %/% 2L
            lo <- max(1L, j - half)
            hi <- min(length(N_sorted), j + half)
            cand <- N_sorted[lo:hi]
            res <- score_candidates(pi, cand, x)
            chosen <- res$chosen
            if (!is.na(chosen)) break
            if (win >= length(N_sorted)) break
            win <- min(length(N_sorted), win * 2L)
          }
        }

        if (!is.na(chosen)) {
          neg_used[chosen] <- TRUE
          matched_neg_for_pos[pi] <- names(gr)[chosen]
          matched_pos_for_neg[chosen] <- names(gr)[pi]
          match_meta_delta[pi] <- abs(meta[chosen] - meta[pi])
        }
      }
    }
  }

  mc$is_positive <- as.integer(y == 1L)
  mc$matched_negative_id <- matched_neg_for_pos
  mc$matched_positive_id <- matched_pos_for_neg
  mc$meta_delta <- match_meta_delta

  matched_neg_idx <- which(!is.na(matched_pos_for_neg))
  matched_pos_idx <- which(!is.na(matched_neg_for_pos))
  mc$is_matched_negative <- as.integer(seq_along(gr) %in% matched_neg_idx)
  mc$is_matched_positive <- as.integer(seq_along(gr) %in% matched_pos_idx)
  mc$match_set <- ifelse(
    mc$is_matched_positive == 1L,
    "positive",
    ifelse(
      mc$is_positive == 1L,
      "unmatched_positive",
      ifelse(mc$is_matched_negative == 1L, "matched_negative", "other")
    )
  )

  S4Vectors::mcols(gr) <- mc

  diag <- list(
    n_pos_in_scope = length(which(y == 1L & in_scope)),
    n_pos_eligible = length(pos_idx),
    n_neg_eligible = length(neg_idx),
    n_matched = sum(!is.na(matched_neg_for_pos[pos_idx])),
    n_unmatched = length(pos_idx) - sum(!is.na(matched_neg_for_pos[pos_idx])),
    meta_col = meta_col,
    kmer_match = kmer_match,
    bin_match = bin_match,
    bin_cols_used = cols_ok,
    bin_mode = bin_mode,
    bin_weight = bin_weight,
    meta_tol = meta_tol,
    enforce_meta_tol = enforce_meta_tol,
    meta_k = as.integer(meta_k)
  )

  if (return_diagnostics) {
    md <- S4Vectors::metadata(gr)
    md$match_diagnostics <- diag
    S4Vectors::metadata(gr) <- md
  }

  gr
}

#' Subset to positives and selected matched negatives
#'
#' @param gr A \code{GRanges}.
#' @param set_col Membership column.
#' @param positive_value Value marking positives.
#' @param negative_value Value marking matched negatives.
#'
#' @return A subset \code{GRanges}.
#' @export
subset_bg_set <- function(gr,
                          set_col = "match_set",
                          positive_value = "positive",
                          negative_value = "matched_negative") {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  if (!(set_col %in% colnames(S4Vectors::mcols(gr)))) stop("Missing set_col: ", set_col)
  s <- as.character(S4Vectors::mcols(gr)[[set_col]])
  gr[s %in% c(positive_value, negative_value)]
}

#' Subset to reciprocal matched positive/background pairs
#'
#' @param gr A \code{GRanges} with canonical match columns.
#' @param matched_negative_id_col Column containing the matched negative ID for positives.
#' @param matched_positive_id_col Column containing the matched positive ID for negatives.
#' @param is_positive_col Optional positive flag column.
#' @param label_col Label column used if \code{is_positive_col} is absent.
#' @param strict_reciprocal Require the negative to point back to the same positive.
#' @param drop_conflicts Drop duplicated negative assignments.
#' @param return_diagnostics Store diagnostics in metadata.
#'
#' @return A paired \code{GRanges}.
#' @export
subset_matched_pairs <- function(gr,
                                 matched_negative_id_col = "matched_negative_id",
                                 matched_positive_id_col = "matched_positive_id",
                                 is_positive_col = NULL,
                                 label_col = "label",
                                 strict_reciprocal = TRUE,
                                 drop_conflicts = TRUE,
                                 return_diagnostics = TRUE) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  gr <- prepare_sites(gr, label = NULL, label_col = label_col, strip_mcols = FALSE)
  mc <- S4Vectors::mcols(gr)
  ids <- names(gr)

  if (!(matched_negative_id_col %in% colnames(mc))) {
    stop("Missing `", matched_negative_id_col, "` in mcols(gr).")
  }
  if (strict_reciprocal && !(matched_positive_id_col %in% colnames(mc))) {
    stop("strict_reciprocal=TRUE but missing `", matched_positive_id_col, "`.")
  }

  work_idx <- seq_along(gr)

  if (!is.null(is_positive_col) && is_positive_col %in% colnames(mc)) {
    pos_flag <- .normalise_binary_label(mc[[is_positive_col]], strict = FALSE)
    pos_idx_all <- which(pos_flag == 1L)
  } else {
    if (!(label_col %in% colnames(mc))) stop("Missing label_col: ", label_col)
    y <- .normalise_binary_label(mc[[label_col]], strict = TRUE)
    pos_idx_all <- which(y == 1L)
  }
  pos_idx_all <- intersect(pos_idx_all, work_idx)

  neg_id_for_pos <- as.character(mc[[matched_negative_id_col]][pos_idx_all])
  ok <- !is.na(neg_id_for_pos) & nzchar(neg_id_for_pos)
  pos_idx <- pos_idx_all[ok]
  pos_ids <- ids[pos_idx]
  neg_ids <- neg_id_for_pos[ok]
  neg_idx <- match(neg_ids, ids)

  ok2 <- !is.na(neg_idx) & neg_idx %in% work_idx
  pos_idx <- pos_idx[ok2]
  pos_ids <- pos_ids[ok2]
  neg_ids <- neg_ids[ok2]
  neg_idx <- neg_idx[ok2]

  if (drop_conflicts && length(neg_ids)) {
    keep <- !duplicated(neg_ids)
    pos_idx <- pos_idx[keep]
    pos_ids <- pos_ids[keep]
    neg_ids <- neg_ids[keep]
    neg_idx <- neg_idx[keep]
  }

  if (strict_reciprocal && length(neg_idx)) {
    back <- as.character(mc[[matched_positive_id_col]][neg_idx])
    ok3 <- !is.na(back) & nzchar(back) & back == pos_ids
    pos_idx <- pos_idx[ok3]
    neg_idx <- neg_idx[ok3]
  }

  keep_idx <- unique(c(pos_idx, neg_idx))
  out <- gr[keep_idx]

  if (return_diagnostics) {
    md <- S4Vectors::metadata(out)
    md$subset_pairs_diagnostics <- list(
      n_total_in = length(gr),
      n_total_out = length(out),
      n_pairs_out_non_testing = length(pos_idx),
      expected_rows_out_non_testing = 2L * length(pos_idx),
      strict_reciprocal = strict_reciprocal,
      drop_conflicts = drop_conflicts
    )
    S4Vectors::metadata(out) <- md
  }

  out
}

#' Subset to matched positive/background rows
#'
#' Convenience alias for \code{\link{subset_matched_pairs}}. This is the
#' object that should usually be passed to matched-set diagnostic plots.
#'
#' @param gr A \code{GRanges} returned by \code{match_background()} or
#'   \code{match_random_background()}.
#' @param ... Passed to \code{\link{subset_matched_pairs}}.
#'
#' @return A paired \code{GRanges} containing only matched positives and their
#'   selected matched negatives.
#' @export
subset_matched_sets <- function(gr, ...) {
  subset_matched_pairs(gr, ...)
}
