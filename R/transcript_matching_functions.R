#' @import GenomicRanges
#' @import GenomicFeatures
#' @import VariantAnnotation
#' @importFrom S4Vectors mcols DataFrame elementNROWS
#' @importFrom IRanges resize
#' @importFrom GenomeInfoDb keepSeqlevels seqlevelsStyle seqinfo
#' @importFrom eulerr euler
NULL

#' Match positives to unique negatives using metagene position and optional bin constraints
#'
#' Creates a background set by pairing each positive site to a unique negative site within
#' the same region (and optionally the same k-mer), prioritizing closeness in metagene
#' coordinate and optionally refining matches using quantile-binned feature covariates.
#'
#' The function writes standard output columns used by downstream plotting helpers:
#' \code{matched_negative_id}, \code{matched_positive_id}, \code{is_matched_negative},
#' \code{is_positive}, \code{meta_delta}, and \code{match_set}.
#'
#' NEW (split-aware): If \code{split_col} is provided, positives are only matched to negatives
#' with the same split label. You can also restrict matching to specific splits via
#' \code{match_splits} (e.g. only "training").
#'
#' @param gr A \code{GRanges} containing labels and required annotation columns.
#' @param label_col Metadata label column (0/1).
#' @param location_col Region column used for matching (default \code{"location"}).
#' @param locations Which regions are eligible for matching.
#' @param return_diagnostics If TRUE attach a \code{"match_diagnostics"} attribute.
#' @param kmer_match If TRUE require exact match on \code{kmer_col}.
#' @param kmer_col Column containing k-mers.
#' @param meta_col Column containing metagene coordinate (default \code{"metagene_prop"}).
#' @param meta_tol Optional tolerance window for metagene closeness.
#' @param enforce_meta_tol If TRUE leave positives unmatched when no negative within tolerance exists.
#' @param meta_k Local search window size in meta-sorted negatives.
#' @param bin_match If TRUE use quantile bins of \code{bin_cols} as a secondary match criterion.
#' @param bin_cols Numeric feature columns to bin.
#' @param n_bins Number of quantile bins.
#' @param bin_within Compute bins globally or within each location.
#' @param bin_mode \code{"soft"} uses a penalty; \code{"hard"} requires exact bin equality.
#' @param bin_weight Weight of the bin penalty relative to metagene closeness.
#' @param log_bin_cols Columns log1p-transformed prior to binning.
#' @param seed Random seed (affects tie-breaking order only).
#'
#' @param split_col Optional column in \code{mcols(gr)} defining data splits (e.g. "dataset_split").
#'   If provided, matching is enforced strictly within split.
#' @param match_splits Optional character vector of split labels to actively match within.
#'   If provided, only rows in those splits are eligible to be matched. Other splits are left
#'   untouched (their match columns remain NA/other).
#' @param drop_na_split If TRUE (default), rows with NA/"" split are ineligible when split_col is used.
#'
#' @returns A \code{GRanges} with standard match columns added.
#'
#' @seealso \code{\link{subset_matched_pairs}}, \code{\link{plot_metagene_density}}
#' @export
match_background <- function(gr,
                             label_col = "label",
                             location_col = "location",
                             locations = c("fiveUTR","coding","threeUTR"),
                             return_diagnostics = FALSE,

                             # 2) optional kmer matching
                             kmer_match = FALSE,
                             kmer_col = "kmer",

                             # 3) metagene closeness (primary)
                             meta_col = "metagene_prop",
                             meta_tol = 0.05,
                             enforce_meta_tol = FALSE,
                             meta_k = 200L,

                             # 4) optional "match level bins" (secondary)
                             bin_match = TRUE,
                             bin_cols = c("feature_width","segment_rank","nexon","tx_len",
                                          "dist_from_feature_start","dist_from_feature_end"),
                             n_bins = 5L,
                             bin_within = c("location","global"),
                             bin_mode = c("soft","hard"),
                             bin_weight = 0.25,
                             log_bin_cols = c("tx_len","feature_width","dist_from_feature_start","dist_from_feature_end"),

                             seed = 1L,

                             # NEW: split-aware matching
                             split_col = NULL,
                             match_splits = NULL,
                             drop_na_split = TRUE) {

  stopifnot(inherits(gr, "GRanges"))
  if (!(label_col %in% colnames(S4Vectors::mcols(gr)))) stop("Missing label_col: ", label_col)
  if (!(location_col %in% colnames(S4Vectors::mcols(gr)))) stop("Missing location_col: ", location_col)
  if (!(meta_col %in% colnames(S4Vectors::mcols(gr)))) stop("Missing meta_col: ", meta_col)
  if (kmer_match && !(kmer_col %in% colnames(S4Vectors::mcols(gr)))) stop("kmer_match=TRUE but missing kmer_col: ", kmer_col)

  bin_within <- match.arg(bin_within)
  bin_mode <- match.arg(bin_mode)

  # Ensure stable names (used as IDs)
  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_along(gr))
  }

  mc <- S4Vectors::mcols(gr)

  y <- mc[[label_col]]
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y))  y <- as.integer(as.character(y))
  y <- as.integer(y)

  loc  <- as.character(mc[[location_col]])
  meta <- suppressWarnings(as.numeric(mc[[meta_col]]))

  # ---- NEW: split-aware eligibility ----
  split_vec <- rep(NA_character_, length(gr))
  split_ok  <- rep(TRUE, length(gr))

  if (!is.null(split_col)) {
    if (!(split_col %in% colnames(mc))) {
      stop("split_col='", split_col, "' not found in mcols(gr).")
    }
    split_vec <- as.character(mc[[split_col]])

    if (isTRUE(drop_na_split)) {
      split_ok <- !is.na(split_vec) & nzchar(split_vec)
    } else {
      split_vec[is.na(split_vec) | !nzchar(split_vec)] <- "__unsplit__"
      split_ok <- rep(TRUE, length(split_vec))
    }

    if (!is.null(match_splits)) {
      match_splits <- as.character(match_splits)
      split_ok <- split_ok & (split_vec %in% match_splits)
    }
  }

  # One-hot location columns (useful for modeling later)
  mc$loc_fiveUTR  <- as.integer(loc == "fiveUTR")
  mc$loc_coding   <- as.integer(loc == "coding")
  mc$loc_threeUTR <- as.integer(loc == "threeUTR")

  # In-scope sites (only these locations can be matched)
  in_scope <- !is.na(loc) & loc %in% locations & is.finite(meta) & split_ok

  pos_idx <- which(y == 1L & in_scope)
  neg_idx <- which(y == 0L & in_scope)

  # If kmer matching is on, require non-missing kmers
  if (kmer_match) {
    km <- as.character(mc[[kmer_col]])
    pos_idx <- pos_idx[!is.na(km[pos_idx]) & nzchar(km[pos_idx])]
    neg_idx <- neg_idx[!is.na(km[neg_idx]) & nzchar(km[neg_idx])]
  }

  # Edge case
  if (!length(pos_idx) || !length(neg_idx)) {
    mc$is_positive <- as.integer(y == 1L)
    mc$is_matched_negative <- 0L
    mc$matched_negative_id <- rep(NA_character_, length(gr))
    mc$matched_positive_id <- rep(NA_character_, length(gr))
    mc$meta_delta <- rep(NA_real_, length(gr))
    mc$match_set <- ifelse(mc$is_positive == 1L, "positive", "other")
    S4Vectors::mcols(gr) <- mc

    if (return_diagnostics) {
      attr(gr, "match_diagnostics") <- list(
        n_pos_in_scope = length(pos_idx),
        n_neg_in_scope = length(neg_idx),
        n_matched = 0L,
        split_col = split_col,
        match_splits = match_splits,
        drop_na_split = drop_na_split,
        reason = "No eligible positives or negatives in scope after filtering."
      )
    }
    return(gr)
  }

  # ---- binning helper (quantile bins) ----
  quantile_bin <- function(x, n_bins) {
    x <- suppressWarnings(as.numeric(x))
    out <- rep.int(NA_integer_, length(x))
    ok <- is.finite(x)
    if (!any(ok)) return(out)

    probs <- seq(0, 1, length.out = n_bins + 1L)
    br <- stats::quantile(x[ok], probs = probs, na.rm = TRUE, type = 7)
    br <- unique(as.numeric(br))
    if (length(br) < 2L) {
      out[ok] <- 1L
      return(out)
    }
    br[1] <- -Inf
    br[length(br)] <- Inf
    out[ok] <- as.integer(cut(x[ok], breaks = br, include.lowest = TRUE, right = TRUE, labels = FALSE))
    out
  }

  # ---- build bin matrix (optional) ----
  bin_mat <- NULL
  if (bin_match) {
    cols_ok <- intersect(bin_cols, colnames(mc))
    if (!length(cols_ok)) stop("bin_match=TRUE but none of bin_cols exist in mcols(gr).")

    get_feat <- function(col) {
      v <- suppressWarnings(as.numeric(mc[[col]]))
      if (col %in% log_bin_cols) v <- log1p(pmax(v, 0))
      v
    }

    bin_mat <- matrix(NA_integer_, nrow = length(gr), ncol = length(cols_ok))
    colnames(bin_mat) <- cols_ok

    if (bin_within == "global") {
      for (j in seq_along(cols_ok)) {
        bin_mat[, j] <- quantile_bin(get_feat(cols_ok[j]), n_bins = n_bins)
      }
    } else {
      for (lv in locations) {
        idx_lv <- which((loc == lv) & (!is.na(loc)) & (loc %in% locations) & split_ok)
        if (!length(idx_lv)) next
        for (j in seq_along(cols_ok)) {
          v <- get_feat(cols_ok[j])
          b <- rep.int(NA_integer_, length(v))
          b[idx_lv] <- quantile_bin(v[idx_lv], n_bins = n_bins)
          bin_mat[, j] <- ifelse(is.na(bin_mat[, j]), b, bin_mat[, j])
        }
      }
    }

    for (j in seq_along(cols_ok)) {
      mc[[paste0("bin_", cols_ok[j])]] <- bin_mat[, j]
    }
  }

  # ---- hard constraints => group keys ----
  if (kmer_match) {
    km <- as.character(mc[[kmer_col]])
    key <- paste0(loc, "||", km)
  } else {
    key <- loc
  }

  # NEW: enforce split boundary as a hard constraint
  if (!is.null(split_col)) {
    key <- paste0(split_vec, "||", key)
  }

  # Only consider keys relevant to our filtered idx
  pos_by_key <- split(pos_idx, key[pos_idx], drop = TRUE)
  neg_by_key <- split(neg_idx, key[neg_idx], drop = TRUE)
  keys <- intersect(names(pos_by_key), names(neg_by_key))
  keys <- keys[!is.na(keys) & nzchar(keys)]

  # Track unused negatives globally
  neg_used <- rep(FALSE, length(gr))
  matched_neg_for_pos <- rep(NA_character_, length(gr))  # negative id for each positive
  matched_pos_for_neg <- rep(NA_character_, length(gr))  # positive id for each negative
  match_meta_delta <- rep(NA_real_, length(gr))

  set.seed(seed)

  # ---- main matching: process each group independently ----
  for (kkey in keys) {
    P <- pos_by_key[[kkey]]
    N <- neg_by_key[[kkey]]
    if (!length(P) || !length(N)) next

    P <- P[order(meta[P])]
    N_sorted <- N[order(meta[N])]
    N_meta <- meta[N_sorted]

    for (pi in P) {
      if (!is.na(matched_neg_for_pos[pi])) next

      avail_mask <- !neg_used[N_sorted]
      if (!any(avail_mask)) break

      x <- meta[pi]
      j <- findInterval(x, N_meta)
      j <- max(1L, min(j, length(N_sorted)))

      win <- as.integer(meta_k)
      if (win < 10L) win <- 10L

      chosen <- NA_integer_

      repeat {
        half <- win %/% 2L
        lo <- max(1L, j - half)
        hi <- min(length(N_sorted), j + half)
        cand <- N_sorted[lo:hi]
        cand <- cand[!neg_used[cand]]

        if (length(cand)) {
          md <- abs(meta[cand] - x)
          if (!is.null(meta_tol) && is.finite(meta_tol) && meta_tol > 0) {
            in_tol <- which(md <= meta_tol)
            if (length(in_tol)) {
              cand <- cand[in_tol]
              md <- md[in_tol]
            } else if (enforce_meta_tol) {
              cand <- integer(0)
            }
          }

          if (length(cand)) {
            meta_scale <- if (!is.null(meta_tol) && is.finite(meta_tol) && meta_tol > 0) meta_tol else 0.05
            cost_meta <- md / meta_scale

            if (bin_match) {
              bp <- bin_mat[pi, , drop = TRUE]
              bn <- bin_mat[cand, , drop = FALSE]

              if (bin_mode == "hard") {
                ok <- rep(TRUE, nrow(bn))
                for (jj in seq_len(ncol(bn))) {
                  ok <- ok & (bn[, jj] == bp[jj] | (is.na(bn[, jj]) & is.na(bp[jj])))
                }
                cand2 <- cand[ok]
                cost_meta2 <- cost_meta[ok]
                md2 <- md[ok]

                if (!length(cand2)) {
                  cand <- integer(0)
                } else {
                  cand <- cand2
                  cost_meta <- cost_meta2
                  md <- md2
                  bn <- bn[ok, , drop = FALSE]
                }
              }

              if (length(cand)) {
                denom <- max(1, n_bins - 1L)
                bin_diff <- abs(sweep(bn, 2, bp, "-"))
                cost_bin <- rowMeans(bin_diff / denom, na.rm = TRUE)
                cost <- cost_meta + bin_weight * cost_bin
              } else {
                cost <- numeric(0)
              }
            } else {
              cost <- cost_meta
            }

            if (length(cand)) {
              o <- order(cost, md)
              chosen <- cand[o[1]]
              break
            }
          }
        }

        if (win >= length(N_sorted)) break
        win <- min(length(N_sorted), win * 2L)
      }

      if (!is.na(chosen)) {
        neg_used[chosen] <- TRUE
        matched_neg_for_pos[pi] <- names(gr)[chosen]
        matched_pos_for_neg[chosen] <- names(gr)[pi]
        match_meta_delta[pi] <- abs(meta[chosen] - meta[pi])
      }
    }
  }

  # ---- attach result columns ----
  mc$is_positive <- as.integer(y == 1L)
  mc$matched_negative_id <- matched_neg_for_pos
  mc$matched_positive_id <- matched_pos_for_neg
  mc$meta_delta <- match_meta_delta

  matched_neg_idx <- which(!is.na(matched_pos_for_neg))
  mc$is_matched_negative <- as.integer(seq_along(gr) %in% matched_neg_idx)

  mc$match_set <- ifelse(mc$is_positive == 1L, "positive",
                         ifelse(mc$is_matched_negative == 1L, "matched_negative", "other"))

  S4Vectors::mcols(gr) <- mc

  if (return_diagnostics) {
    n_pos_scope <- length(which(y == 1L & in_scope))
    n_pos_eligible <- length(pos_idx)
    n_matched <- sum(!is.na(matched_neg_for_pos[pos_idx]))

    attr(gr, "match_diagnostics") <- list(
      n_pos_in_scope = n_pos_scope,
      n_pos_eligible = n_pos_eligible,
      n_neg_eligible = length(neg_idx),
      n_matched = n_matched,
      n_unmatched = n_pos_eligible - n_matched,
      kmer_match = kmer_match,
      bin_match = bin_match,
      bin_mode = bin_mode,
      bin_weight = bin_weight,
      meta_tol = meta_tol,
      enforce_meta_tol = enforce_meta_tol,
      meta_k = as.integer(meta_k),
      split_col = split_col,
      match_splits = match_splits,
      drop_na_split = drop_na_split
    )
  }

  gr
}

#' Sample a random DRACH background site within the same transcript
#'
#' For each positive site, attempts to select a unique negative DRACH site from the same
#' transcript, optionally matching location and/or exact k-mer.
#'
#' NEW (split-aware): If \code{split_col} is provided, positives are only matched to negatives
#' with the same split label. You can also restrict matching to specific splits via
#' \code{match_splits} (e.g. only "training").
#'
#' @export
random_drach_within_transcript <- function(gr,
                                           label_col = "label",
                                           tx_col = "tx_name",
                                           location_col = "location",
                                           locations = c("fiveUTR", "coding", "threeUTR"),
                                           kmer_col = "kmer",
                                           match_location = TRUE,
                                           match_kmer = FALSE,
                                           meta_col = "metagene_prop",
                                           seed = 1L,
                                           prefix = "txrand",
                                           write_standard_cols = TRUE,
                                           overwrite_standard_cols = TRUE,
                                           assert_invariants = TRUE,
                                           return_diagnostics = TRUE,
                                           # NEW: split-aware matching
                                           split_col = NULL,
                                           match_splits = NULL,
                                           drop_na_split = TRUE) {
  stopifnot(inherits(gr, "GRanges"))

  mc <- S4Vectors::mcols(gr)
  need <- c(label_col, tx_col, location_col, kmer_col)
  miss <- setdiff(need, colnames(mc))
  if (length(miss)) {
    stop("Missing required mcols: ", paste(miss, collapse = ", "),
         ". (Need annotate_sites() + add_kmer(k=5) first.)")
  }

  n <- length(gr)

  # Stable row IDs
  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_len(n))
  }
  id <- names(gr)

  # Labels -> integer 0/1
  y <- mc[[label_col]]
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y))  y <- as.integer(as.character(y))
  y <- as.integer(y)
  if (any(!is.na(y) & !(y %in% c(0L, 1L)))) stop("`", label_col, "` must be 0/1.")

  tx  <- as.character(mc[[tx_col]])
  loc <- as.character(mc[[location_col]])
  km  <- toupper(as.character(mc[[kmer_col]]))
  km  <- chartr("U", "T", km)

  # ---- NEW: split-aware eligibility ----
  split_vec <- rep(NA_character_, n)
  split_ok  <- rep(TRUE, n)

  if (!is.null(split_col)) {
    if (!(split_col %in% colnames(mc))) {
      stop("split_col='", split_col, "' not found in mcols(gr).")
    }
    split_vec <- as.character(mc[[split_col]])

    if (isTRUE(drop_na_split)) {
      split_ok <- !is.na(split_vec) & nzchar(split_vec)
    } else {
      split_vec[is.na(split_vec) | !nzchar(split_vec)] <- "__unsplit__"
      split_ok <- rep(TRUE, length(split_vec))
    }

    if (!is.null(match_splits)) {
      match_splits <- as.character(match_splits)
      split_ok <- split_ok & (split_vec %in% match_splits)
    }
  }

  # DRACH in DNA alphabet: D=[AGT], R=[AG], A, C, H=[ACT]
  is_drach_5mer <- function(x) {
    x <- toupper(as.character(x))
    x <- chartr("U","T", x)
    ok <- !is.na(x) & nchar(x) == 5L
    out <- rep(FALSE, length(x))
    out[ok] <- grepl("^[AGT][AG]AC[ACT]$", x[ok])
    out
  }

  # In-scope rows must have tx/loc/kmer (+ split_ok if split_col used)
  in_scope <- !is.na(tx) & nzchar(tx) &
    !is.na(loc) & nzchar(loc) &
    !is.na(km) & nzchar(km) &
    split_ok

  if (!is.null(locations)) in_scope <- in_scope & (loc %in% locations)

  pos_all <- which(y == 1L & in_scope)
  neg_all <- which(y == 0L & in_scope)

  # Candidate negatives must be DRACH
  neg_all <- neg_all[is_drach_5mer(km[neg_all])]

  # Output columns (clear old prefix cols to avoid sticky values)
  neg_id_col <- paste0(prefix, "_negative_id")  # filled for positives
  pos_id_col <- paste0(prefix, "_positive_id")  # filled for selected negatives
  isneg_col  <- paste0(prefix, "_is_negative")
  set_col    <- paste0(prefix, "_set")
  dmeta_col  <- paste0(prefix, "_meta_delta")

  for (col in c(neg_id_col, pos_id_col, isneg_col, set_col, dmeta_col)) {
    if (col %in% colnames(mc)) mc[[col]] <- NULL
  }

  # Vectors to fill
  neg_id <- rep(NA_character_, n)
  pos_id <- rep(NA_character_, n)
  isneg  <- rep.int(0L, n)
  setv   <- ifelse(y == 1L, "positive", "other")
  meta_delta <- rep(NA_real_, n)

  meta <- NULL
  if (!is.null(meta_col) && meta_col %in% colnames(mc)) {
    meta <- suppressWarnings(as.numeric(mc[[meta_col]]))
  }

  # Early exit
  if (!length(pos_all) || !length(neg_all)) {
    mc[[neg_id_col]] <- neg_id
    mc[[pos_id_col]] <- pos_id
    mc[[isneg_col]]  <- isneg
    mc[[set_col]]    <- setv
    if (!is.null(meta)) mc[[dmeta_col]] <- meta_delta
    S4Vectors::mcols(gr) <- mc

    if (write_standard_cols) {
      gr <- .txrand_write_standard_cols(
        gr,
        y = y,
        matched_negative_id = neg_id,
        matched_positive_id = pos_id,
        is_matched_negative = isneg,
        meta_delta = meta_delta,
        overwrite = overwrite_standard_cols
      )
    }

    if (return_diagnostics) {
      attr(gr, paste0(prefix, "_diagnostics")) <- list(
        n_pos_total = sum(y == 1L, na.rm = TRUE),
        n_neg_total = sum(y == 0L, na.rm = TRUE),
        n_pos_in_scope = length(pos_all),
        n_neg_drach_in_scope = length(neg_all),
        split_col = split_col,
        match_splits = match_splits,
        drop_na_split = drop_na_split,
        n_pairs = 0L,
        reason = "No eligible positives or no eligible DRACH negatives in scope."
      )
    }
    return(gr)
  }

  # Matching key builder
  make_key <- function(idx) {
    k <- tx[idx]
    if (!is.null(split_col)) k <- paste0(split_vec[idx], "||", k)  # NEW: split hard constraint
    if (match_location)      k <- paste0(k, "||", loc[idx])
    if (match_kmer)          k <- paste0(k, "||", km[idx])
    k
  }

  key_pos <- make_key(pos_all)
  key_neg <- make_key(neg_all)

  okp <- !is.na(key_pos) & nzchar(key_pos)
  okn <- !is.na(key_neg) & nzchar(key_neg)

  pos <- pos_all[okp]; key_pos <- key_pos[okp]
  neg <- neg_all[okn]; key_neg <- key_neg[okn]

  if (!length(pos) || !length(neg)) {
    mc[[neg_id_col]] <- neg_id
    mc[[pos_id_col]] <- pos_id
    mc[[isneg_col]]  <- isneg
    mc[[set_col]]    <- setv
    if (!is.null(meta)) mc[[dmeta_col]] <- meta_delta
    S4Vectors::mcols(gr) <- mc

    if (write_standard_cols) {
      gr <- .txrand_write_standard_cols(
        gr,
        y = y,
        matched_negative_id = neg_id,
        matched_positive_id = pos_id,
        is_matched_negative = isneg,
        meta_delta = meta_delta,
        overwrite = overwrite_standard_cols
      )
    }

    if (return_diagnostics) {
      attr(gr, paste0(prefix, "_diagnostics")) <- list(
        n_pos_total = sum(y == 1L, na.rm = TRUE),
        n_neg_total = sum(y == 0L, na.rm = TRUE),
        n_pos_in_scope = length(pos_all),
        n_neg_drach_in_scope = length(neg_all),
        split_col = split_col,
        match_splits = match_splits,
        drop_na_split = drop_na_split,
        n_pairs = 0L,
        reason = "No overlap of matching keys between positives and DRACH negatives."
      )
    }
    return(gr)
  }

  # Sort by key and walk runs
  op <- order(key_pos)
  on <- order(key_neg)
  pos_s <- pos[op]; key_pos_s <- key_pos[op]
  neg_s <- neg[on]; key_neg_s <- key_neg[on]

  rp <- rle(key_pos_s)
  rn <- rle(key_neg_s)

  p_end <- cumsum(rp$lengths)
  p_start <- p_end - rp$lengths + 1L
  n_end <- cumsum(rn$lengths)
  n_start <- n_end - rn$lengths + 1L

  i <- 1L; j <- 1L
  pairs <- 0L

  set.seed(seed)

  while (i <= length(rp$values) && j <= length(rn$values)) {
    kp <- rp$values[i]
    kn <- rn$values[j]

    if (kp == kn) {
      P <- pos_s[p_start[i]:p_end[i]]
      N <- neg_s[n_start[j]:n_end[j]]

      if (length(P) && length(N)) {
        m <- min(length(P), length(N))
        if (m > 0L) {
          Ppick <- P[sample.int(length(P), m, replace = FALSE)]
          Npick <- N[sample.int(length(N), m, replace = FALSE)]
          Npick <- Npick[sample.int(length(Npick), length(Npick), replace = FALSE)]

          neg_id[Ppick] <- id[Npick]
          pos_id[Npick] <- id[Ppick]
          isneg[Npick]  <- 1L
          setv[Npick]   <- "matched_negative"

          if (!is.null(meta)) {
            okm <- is.finite(meta[Ppick]) & is.finite(meta[Npick])
            md <- rep(NA_real_, length(Ppick))
            md[okm] <- abs(meta[Ppick[okm]] - meta[Npick[okm]])
            meta_delta[Ppick] <- md
          }

          pairs <- pairs + m
        }
      }

      i <- i + 1L
      j <- j + 1L
    } else if (kp < kn) {
      i <- i + 1L
    } else {
      j <- j + 1L
    }
  }

  # Attach prefix results
  mc[[neg_id_col]] <- neg_id
  mc[[pos_id_col]] <- pos_id
  mc[[isneg_col]]  <- isneg
  mc[[set_col]]    <- setv
  if (!is.null(meta)) mc[[dmeta_col]] <- meta_delta
  S4Vectors::mcols(gr) <- mc

  if (assert_invariants) {
    pos_idx <- which(y == 1L)
    neg_idx <- which(y == 0L)

    bad_pos_have_posid <- sum(!is.na(pos_id[pos_idx]))
    bad_neg_have_negid <- sum(!is.na(neg_id[neg_idx]))

    if (bad_pos_have_posid > 0L || bad_neg_have_negid > 0L) {
      stop(
        "Invariant violation: positives have ", pos_id_col, " filled (", bad_pos_have_posid,
        ") and/or negatives have ", neg_id_col, " filled (", bad_neg_have_negid, ")."
      )
    }
    if (any(isneg[pos_idx] != 0L, na.rm = TRUE)) {
      stop("Invariant violation: at least one label==1 site was marked as a negative.")
    }

    n_pos_matched <- sum(!is.na(neg_id[pos_idx]))
    n_neg_matched <- sum(!is.na(pos_id[neg_idx]))
    if (n_pos_matched != n_neg_matched) {
      stop("Invariant violation: matched positives (", n_pos_matched,
           ") != matched negatives (", n_neg_matched, ").")
    }
  }

  if (write_standard_cols) {
    gr <- .txrand_write_standard_cols(
      gr,
      y = y,
      matched_negative_id = neg_id,
      matched_positive_id = pos_id,
      is_matched_negative = isneg,
      meta_delta = meta_delta,
      overwrite = overwrite_standard_cols
    )
  }

  if (return_diagnostics) {
    attr(gr, paste0(prefix, "_diagnostics")) <- list(
      n_pos_total = sum(y == 1L, na.rm = TRUE),
      n_neg_total = sum(y == 0L, na.rm = TRUE),
      n_pos_in_scope = length(pos_all),
      n_neg_drach_in_scope = length(neg_all),
      match_location = match_location,
      match_kmer = match_kmer,
      split_col = split_col,
      match_splits = match_splits,
      drop_na_split = drop_na_split,
      n_pairs = pairs,
      note = "Prefix columns written; canonical match_* columns optionally written."
    )
  }

  gr
}

#' Internal helper - standardize columns for plotting + parity with match_background_simple()
#' @noRd
.txrand_write_standard_cols <- function(gr,
                                        y,
                                        matched_negative_id,
                                        matched_positive_id,
                                        is_matched_negative,
                                        meta_delta,
                                        overwrite = TRUE) {
  mc <- S4Vectors::mcols(gr)

  standard_cols <- c("matched_negative_id", "matched_positive_id",
                     "is_matched_negative", "is_positive", "meta_delta", "match_set")

  if (!overwrite) {
    exists_any <- any(standard_cols %in% colnames(mc))
    if (exists_any) {
      warning("Standard columns already exist and overwrite=FALSE; not overwriting.")
      return(gr)
    }
  }

  mc$matched_negative_id <- matched_negative_id
  mc$matched_positive_id <- matched_positive_id
  mc$is_matched_negative <- as.integer(is_matched_negative == 1L)
  mc$is_positive <- as.integer(y == 1L)
  mc$meta_delta <- meta_delta

  mc$match_set <- ifelse(
    mc$is_positive == 1L, "positive",
    ifelse(mc$is_matched_negative == 1L, "matched_negative", "other")
  )

  S4Vectors::mcols(gr) <- mc
  gr
}



#' Subset a GRanges to positives and selected negatives
#'
#' Convenience helper to keep only rows labelled as positives and chosen negatives
#' according to a membership column (e.g. \code{"match_set"}).
#' This is a fast filter and does not enforce 1:1 pairing; use
#' \code{\link{subset_matched_pairs}} when you need exactly one matched negative per positive.
#'
#' @param gr A \code{GRanges}.
#' @param set_col Column defining membership (default \code{"match_set"}).
#' @param positive_value Value marking positives (default \code{"positive"}).
#' @param negative_value Value marking selected negatives (default \code{"matched_negative"}).
#'
#' @returns A subset \code{GRanges}.
#' @seealso \code{\link{subset_matched_pairs}}
#' @export
subset_bg_set <- function(gr,
                          set_col = "match_set",
                          positive_value = "positive",
                          negative_value = "matched_negative") {
  stopifnot(inherits(gr, "GRanges"))
  if (!(set_col %in% colnames(S4Vectors::mcols(gr)))) {
    stop("Missing set_col: ", set_col)
  }
  s <- as.character(S4Vectors::mcols(gr)[[set_col]])
  keep <- s %in% c(positive_value, negative_value)
  gr[keep]
}


#' Internal helper
#' @noRd
.bg_plot_df <- function(gr,
                        set_col,
                        positive_value,
                        negative_value,
                        cols,
                        drop_na = TRUE) {
  mc <- S4Vectors::mcols(gr)

  if (!(set_col %in% colnames(mc))) stop("Missing set_col: ", set_col)
  miss <- setdiff(cols, colnames(mc))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  s <- as.character(mc[[set_col]])
  group <- ifelse(s == positive_value, "positive",
                  ifelse(s == negative_value, "negative", NA_character_))

  df <- data.frame(group = group, stringsAsFactors = FALSE)
  for (cl in cols) df[[cl]] <- mc[[cl]]

  # Keep only pos/neg rows
  df <- df[!is.na(df$group), , drop = FALSE]
  df$group <- factor(df$group, levels = c("positive", "negative"))

  if (drop_na) {
    ok <- rep(TRUE, nrow(df))
    for (cl in cols) ok <- ok & !is.na(df[[cl]])
    df <- df[ok, , drop = FALSE]
  }

  df
}


#' @family plots
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' # p <- plot_metagene_density(gr)
#' @export
plot_metagene_density <- function(gr,
                                  set_col = "match_set",
                                  positive_value = "positive",
                                  negative_value = "matched_negative",
                                  meta_col = "metagene_prop",
                                  facet_by_location = FALSE,
                                  location_col = "location",
                                  bw_adjust = 1,
                                  xlim = c(0, 3)) {
  stopifnot(inherits(gr, "GRanges"))
  cols <- c(meta_col, if (facet_by_location) location_col else NULL)

  df <- .bg_plot_df(gr, set_col, positive_value, negative_value, cols = cols, drop_na = TRUE)
  df[[meta_col]] <- suppressWarnings(as.numeric(df[[meta_col]]))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[meta_col]], colour = group)) +
    ggplot2::geom_density(linewidth = 1, adjust = bw_adjust, na.rm = TRUE) +  # density :contentReference[oaicite:2]{index=2}
    ggplot2::geom_vline(xintercept = c(1, 2), linetype = 2) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = meta_col, y = "Density", colour = NULL) +
    ggplot2::theme_bw()

  if (facet_by_location) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", location_col))) # facets :contentReference[oaicite:3]{index=3}
  }

  p
}


#' @family plots
#' Plot distance to nearest feature boundary (splice-proximity proxy)
#'
#' Uses nearest boundary distance = min(dist_from_feature_start, dist_from_feature_end).
#' Optionally log1p-transform to handle heavy tails.
#'
#' @param gr GRanges
#' @param set_col membership column
#' @param positive_value positives label
#' @param negative_value negatives label
#' @param start_col distance-from-feature-start column
#' @param end_col distance-from-feature-end column
#' @param transform "log1p" or "identity"
#' @param facet_by_location if TRUE, facet by location_col
#' @param location_col facet column
#' @param bw_adjust density smoothing
#'
#' @return ggplot object
#' @export
plot_splice_distance_density <- function(gr,
                                         set_col = "match_set",
                                         positive_value = "positive",
                                         negative_value = "matched_negative",
                                         start_col = "dist_from_feature_start",
                                         end_col = "dist_from_feature_end",
                                         transform = c("log1p", "identity"),
                                         facet_by_location = TRUE,
                                         location_col = "location",
                                         bw_adjust = 1) {
  transform <- match.arg(transform)

  cols <- c(start_col, end_col, if (facet_by_location) location_col else NULL)
  df <- .bg_plot_df(gr, set_col, positive_value, negative_value, cols = cols, drop_na = TRUE)

  a <- suppressWarnings(as.numeric(df[[start_col]]))
  b <- suppressWarnings(as.numeric(df[[end_col]]))
  d <- pmin(a, b)

  if (transform == "log1p") {
    df$value <- log1p(pmax(d, 0))
    xlab <- paste0("log1p(min(", start_col, ", ", end_col, "))")
  } else {
    df$value <- d
    xlab <- paste0("min(", start_col, ", ", end_col, ")")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, colour = group)) +
    ggplot2::geom_density(linewidth = 1, adjust = bw_adjust, na.rm = TRUE) +
    ggplot2::labs(x = xlab, y = "Density", colour = NULL) +
    ggplot2::theme_bw()

  if (facet_by_location) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", location_col)))
  }

  p
}


#' @family plots
#' Plot distance to start and stop codons (density)
#'
#' @param gr GRanges
#' @param set_col membership column
#' @param positive_value positives label
#' @param negative_value negatives label
#' @param start_dist_col column name for start distance (default "start_dist")
#' @param stop_dist_col column name for stop distance (default "stop_dist")
#' @param transform "log1p" or "identity"
#' @param bw_adjust density smoothing
#'
#' @return ggplot object (faceted by start/stop)
#' @export
plot_start_stop_distance_density <- function(gr,
                                             set_col = "match_set",
                                             positive_value = "positive",
                                             negative_value = "matched_negative",
                                             start_dist_col = "start_dist",
                                             stop_dist_col = "stop_dist",
                                             transform = c("log1p", "identity"),
                                             bw_adjust = 1) {
  transform <- match.arg(transform)

  df <- .bg_plot_df(gr, set_col, positive_value, negative_value,
                    cols = c(start_dist_col, stop_dist_col),
                    drop_na = FALSE)

  startv <- suppressWarnings(as.numeric(df[[start_dist_col]]))
  stopv  <- suppressWarnings(as.numeric(df[[stop_dist_col]]))

  long <- rbind(
    data.frame(group = df$group, which = "start", value = startv),
    data.frame(group = df$group, which = "stop",  value = stopv)
  )
  long <- long[is.finite(long$value), , drop = FALSE]

  if (transform == "log1p") {
    long$value <- log1p(pmax(long$value, 0))
    xlab <- "log1p(distance)"
  } else {
    xlab <- "distance"
  }

  ggplot2::ggplot(long, ggplot2::aes(x = value, colour = group)) +
    ggplot2::geom_density(linewidth = 1, adjust = bw_adjust, na.rm = TRUE) +
    ggplot2::facet_wrap(~ which, scales = "free_x") +
    ggplot2::labs(x = xlab, y = "Density", colour = NULL) +
    ggplot2::theme_bw()
}

#' @family plots
#' Plot counts by exon/segment length bins
#'
#' @param gr GRanges
#' @param set_col membership column
#' @param positive_value positives label
#' @param negative_value negatives label
#' @param width_col segment width column (default "feature_width")
#' @param n_bins number of quantile bins
#' @param log1p_first if TRUE, bin log1p(width) instead of width
#' @param facet_by_location facet by location if TRUE
#' @param location_col location column
#'
#' @return ggplot object
#' @export
plot_feature_width_bins <- function(gr,
                                    set_col = "match_set",
                                    positive_value = "positive",
                                    negative_value = "matched_negative",
                                    width_col = "feature_width",
                                    n_bins = 10,
                                    log1p_first = TRUE,
                                    facet_by_location = TRUE,
                                    location_col = "location") {
  cols <- c(width_col, if (facet_by_location) location_col else NULL)
  df <- .bg_plot_df(gr, set_col, positive_value, negative_value, cols = cols, drop_na = TRUE)

  w <- suppressWarnings(as.numeric(df[[width_col]]))
  w <- w[is.finite(w)]
  if (!length(w)) stop("No finite values in ", width_col, " after filtering.")

  v <- suppressWarnings(as.numeric(df[[width_col]]))
  if (log1p_first) v <- log1p(pmax(v, 0))

  # Quantile breaks (guard against ties)
  br <- stats::quantile(v, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE, type = 7)
  br <- unique(as.numeric(br))
  if (length(br) < 2L) {
    df$bin <- factor("all")
  } else {
    br[1] <- -Inf
    br[length(br)] <- Inf
    df$bin <- cut(v, breaks = br, include.lowest = TRUE, labels = FALSE)
    df$bin <- factor(df$bin, levels = sort(unique(df$bin)))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = bin, fill = group)) +
    ggplot2::geom_bar(position = "dodge") +
    ggplot2::labs(x = if (log1p_first) paste0("Quantile bins of log1p(", width_col, ")") else paste0("Quantile bins of ", width_col),
                  y = "Count", fill = NULL) +
    ggplot2::theme_bw()

  if (facet_by_location) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", location_col)))
  }

  p
}

#' @family plots
#' Plot top-kmer composition for positive vs negative sets
#'
#' @param gr GRanges
#' @param set_col membership column
#' @param positive_value positives label
#' @param negative_value negatives label
#' @param kmer_col kmer column (default "kmer")
#' @param top_n how many kmers to show (overall)
#'
#' @return ggplot object
#' @export
plot_kmer_counts <- function(gr,
                             set_col = "match_set",
                             positive_value = "positive",
                             negative_value = "matched_negative",
                             kmer_col = "kmer",
                             top_n = 20) {
  df <- .bg_plot_df(gr, set_col, positive_value, negative_value, cols = c(kmer_col), drop_na = TRUE)
  km <- toupper(as.character(df[[kmer_col]]))
  df$kmer <- km

  # overall top kmers
  tab <- sort(table(df$kmer), decreasing = TRUE)
  keep <- names(tab)[seq_len(min(top_n, length(tab)))]
  df <- df[df$kmer %in% keep, , drop = FALSE]
  df$kmer <- factor(df$kmer, levels = keep)

  ggplot2::ggplot(df, ggplot2::aes(x = kmer, fill = group)) +
    ggplot2::geom_bar(position = "dodge") +
    ggplot2::labs(x = kmer_col, y = "Count", fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}



#' Subset a GRanges to balanced matched positive/negative pairs
#'
#' Extracts a strictly paired dataset from canonical match columns
#' (matched_negative_id, matched_positive_id) so the output contains
#' exactly one negative per positive for all retained pairs.
#'
#' Optionally, you can preserve entire split(s) in "raw" form (unbalanced),
#' e.g. keep testing (and/or validation) untouched while subsetting only training.
#'
#' @param gr A GRanges with canonical match columns produced by match_background()
#'   or random_drach_within_transcript().
#' @param matched_negative_id_col Column containing the matched negative ID for each positive.
#' @param matched_positive_id_col Column containing the matched positive ID for each negative.
#' @param is_positive_col Optional column indicating positives (1/0). If NULL, uses label_col.
#' @param label_col Label column (0/1) used when is_positive_col is NULL.
#' @param strict_reciprocal If TRUE enforce reciprocal links (negative points back to the same positive).
#' @param drop_conflicts If TRUE drop duplicated negative assignments if present.
#' @param maintain_testing Optional column name (e.g. "dataset_split") defining splits.
#' @param testing_value Split label(s) to preserve in raw form. Can be a character vector,
#'   e.g. c("testing","validation"). Default is "testing".
#' @param return_diagnostics If TRUE attach a "subset_pairs_diagnostics" attribute.
#'
#' @returns A GRanges containing matched pairs from non-heldout splits, plus all preserved split rows.
#' @export
subset_matched_pairs <- function(gr,
                                 matched_negative_id_col = "matched_negative_id",
                                 matched_positive_id_col = "matched_positive_id",
                                 is_positive_col = NULL,
                                 label_col = "label",
                                 strict_reciprocal = TRUE,
                                 drop_conflicts = TRUE,
                                 maintain_testing = NULL,
                                 testing_value = "testing",
                                 return_diagnostics = TRUE) {
  stopifnot(inherits(gr, "GRanges"))

  mc <- S4Vectors::mcols(gr)

  # Stable row names used as IDs
  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_len(length(gr)))
  }
  ids <- names(gr)

  # Required match columns
  if (!(matched_negative_id_col %in% colnames(mc))) {
    stop("Missing `", matched_negative_id_col, "` in mcols(gr). Cannot subset to matched pairs.")
  }
  if (strict_reciprocal && !(matched_positive_id_col %in% colnames(mc))) {
    stop("strict_reciprocal=TRUE but missing `", matched_positive_id_col, "` in mcols(gr).")
  }

  # -------------------------
  # Define held-out rows to preserve (raw/unbalanced)
  # -------------------------
  heldout_idx <- integer(0)
  work_idx <- seq_along(gr)

  if (!is.null(maintain_testing)) {
    if (!(maintain_testing %in% colnames(mc))) {
      stop("maintain_testing='", maintain_testing, "' but column not found in mcols(gr). ",
           "Did you run assign_split_by_chromosome(..., split_col='", maintain_testing, "') first?")
    }

    split_vec <- as.character(mc[[maintain_testing]])

    testing_value <- unique(as.character(testing_value))
    testing_value <- testing_value[!is.na(testing_value) & nzchar(testing_value)]

    if (length(testing_value)) {
      heldout_idx <- which(!is.na(split_vec) & split_vec %in% testing_value)
      work_idx <- setdiff(work_idx, heldout_idx)
    }

    # If everything is held-out, return full GRanges unchanged
    if (!length(work_idx)) {
      out <- gr
      if (return_diagnostics) {
        attr(out, "subset_pairs_diagnostics") <- list(
          n_total_in = length(gr),
          maintain_testing = maintain_testing,
          heldout_values = testing_value,
          n_heldout_in = length(heldout_idx),
          n_pairs_out_non_heldout = 0L,
          note = "All rows are in held-out split(s); returning input unchanged."
        )
      }
      return(out)
    }
  }

  # -------------------------
  # Determine positives within work_idx
  # -------------------------
  if (!is.null(is_positive_col) && (is_positive_col %in% colnames(mc))) {
    pos_flag <- mc[[is_positive_col]]
    if (is.logical(pos_flag)) pos_flag <- as.integer(pos_flag)
    if (is.factor(pos_flag))  pos_flag <- as.integer(as.character(pos_flag))
    pos_flag <- as.integer(pos_flag)
    pos_idx_all <- which(pos_flag == 1L)
  } else {
    if (!(label_col %in% colnames(mc))) {
      stop("Missing `", label_col, "` in mcols(gr).")
    }
    y <- mc[[label_col]]
    if (is.logical(y)) y <- as.integer(y)
    if (is.factor(y))  y <- as.integer(as.character(y))
    y <- as.integer(y)
    pos_idx_all <- which(y == 1L)
  }

  pos_idx_all <- intersect(pos_idx_all, work_idx)

  # Positive -> matched negative ID
  neg_id_for_pos <- as.character(mc[[matched_negative_id_col]][pos_idx_all])
  ok <- !is.na(neg_id_for_pos) & nzchar(neg_id_for_pos)

  pos_idx <- pos_idx_all[ok]
  pos_ids <- ids[pos_idx]
  neg_ids <- neg_id_for_pos[ok]

  # referenced negatives must exist AND be in work_idx
  neg_idx <- match(neg_ids, ids)
  ok2 <- !is.na(neg_idx) & (neg_idx %in% work_idx)

  pos_idx <- pos_idx[ok2]
  pos_ids <- pos_ids[ok2]
  neg_ids <- neg_ids[ok2]
  neg_idx <- neg_idx[ok2]

  if (!length(pos_idx)) {
    out <- gr[heldout_idx]
    if (return_diagnostics) {
      attr(out, "subset_pairs_diagnostics") <- list(
        n_total_in = length(gr),
        n_total_out = length(out),
        n_pos_total_in_non_heldout = length(pos_idx_all),
        n_pairs_out_non_heldout = 0L,
        maintain_testing = maintain_testing,
        heldout_values = unique(as.character(testing_value)),
        n_heldout_in = length(heldout_idx),
        n_heldout_out = length(heldout_idx),
        reason = "No positives outside held-out splits with a valid matched_negative_id pointing to a non-heldout negative."
      )
    }
    return(out)
  }

  # Optional: drop duplicate negative assignments
  if (drop_conflicts) {
    dup_neg <- duplicated(neg_ids)
    if (any(dup_neg)) {
      keep <- !dup_neg
      pos_idx <- pos_idx[keep]
      pos_ids <- pos_ids[keep]
      neg_ids <- neg_ids[keep]
      neg_idx <- neg_idx[keep]
    }
  }

  # Optional: enforce reciprocal mapping
  if (strict_reciprocal) {
    back <- as.character(mc[[matched_positive_id_col]][neg_idx])
    ok3 <- !is.na(back) & nzchar(back) & (back == pos_ids)

    pos_idx <- pos_idx[ok3]
    pos_ids <- pos_ids[ok3]
    neg_ids <- neg_ids[ok3]
    neg_idx <- neg_idx[ok3]
  }

  # Keep:
  # - matched pairs from non-heldout
  # - plus all heldout rows raw
  keep_pairs_idx <- unique(c(pos_idx, neg_idx))
  keep_idx <- unique(c(keep_pairs_idx, heldout_idx))
  out <- gr[keep_idx]

  if (return_diagnostics) {
    mc2 <- S4Vectors::mcols(out)
    y2 <- if (label_col %in% colnames(mc2)) as.integer(mc2[[label_col]]) else rep(NA_integer_, length(out))
    split2 <- if (!is.null(maintain_testing) && maintain_testing %in% colnames(mc2)) as.character(mc2[[maintain_testing]]) else NULL

    heldout_values <- unique(as.character(testing_value))

    attr(out, "subset_pairs_diagnostics") <- list(
      n_total_in = length(gr),
      n_total_out = length(out),
      n_pairs_out_non_heldout = length(pos_idx),
      expected_rows_out_non_heldout = 2L * length(pos_idx),
      strict_reciprocal = strict_reciprocal,
      drop_conflicts = drop_conflicts,
      maintain_testing = maintain_testing,
      heldout_values = heldout_values,
      n_heldout_in = length(heldout_idx),
      n_heldout_out = if (is.null(split2)) NA_integer_ else sum(!is.na(split2) & split2 %in% heldout_values),
      n_pos_out = sum(y2 == 1L, na.rm = TRUE),
      n_neg_out = sum(y2 == 0L, na.rm = TRUE)
    )
  }

  out
}




#' Euler/Venn-style overlap plot from a named list
#'
#' Builds an incidence matrix from a named list of character vectors and plots an Euler diagram.
#'
#' @param mylist A named list of character vectors (e.g. gene sets).
#'
#' @returns A plot (invisibly returns the \code{eulerr} fit object if you choose to return it).
#' @examplesIf requireNamespace("eulerr", quietly = TRUE)
#' makeVenn(list(A = c("x","y"), B = c("y","z")))
#' @export
makeVenn <- function(mylist)
{
  if (!requireNamespace("eulerr", quietly = TRUE)) {
    stop("Package 'eulerr' is required for makeVenn(). Install it first.")
  }
  all <- unique(unlist(mylist))
  all <- all[!(all=="")]
  dmat <- matrix(0,ncol=length(mylist),nrow=length(all))
  row.names(dmat) <- all
  colnames(dmat) <- names(mylist)
  for(i in 1:length(mylist))
  {
    mgenes <- mylist[[i]]
    mgenes <- mgenes[!(mgenes=="")]
    dmat[mgenes,i] <- 1
  }
  genes.venn <- euler(dmat)
  plot(genes.venn, quantities = TRUE)
}



#' Assign train/test/validation split by chromosome
#'
#' Adds a metadata column to a GRanges object indicating whether each site belongs
#' to the training, testing, or validation set based on chromosome.
#'
#' Default behavior:
#' - chr20 -> validation
#' - chr9  -> testing
#' - everything else -> training
#'
#' The function attempts to harmonize "chr" prefix style (e.g. "9" vs "chr9")
#' between user-supplied chromosomes and seqnames(gr).
#'
#' @param gr GRanges.
#' @param split_col Name of the output metadata column.
#' @param test_chrs Character vector of chromosomes to assign as testing.
#' @param val_chrs Character vector of chromosomes to assign as validation.
#' @param train_label Label used for training rows.
#' @param test_label Label used for testing rows.
#' @param val_label Label used for validation rows.
#' @param normalize_chr_prefix If TRUE, attempt to match "chr" prefix style.
#' @param warn_missing If TRUE, warn if any requested chromosomes are not present.
#' @param overwrite If FALSE and split_col already exists, error.
#'
#' @return GRanges with an added metadata column split_col (factor).
#' @export
assign_split_by_chromosome <- function(gr,
                                       split_col = "split",
                                       test_chrs = "chr9",
                                       val_chrs  = "chr20",
                                       train_label = "training",
                                       test_label  = "testing",
                                       val_label   = "validation",
                                       normalize_chr_prefix = TRUE,
                                       warn_missing = TRUE,
                                       overwrite = TRUE) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")

  mc <- S4Vectors::mcols(gr)
  if (!overwrite && (split_col %in% colnames(mc))) {
    stop("`", split_col, "` already exists in mcols(gr). Set overwrite=TRUE to replace.")
  }

  seqs <- as.character(GenomeInfoDb::seqnames(gr))
  seqs <- ifelse(is.na(seqs), NA_character_, seqs)

  # Helper: drop NA/empty and coerce to character unique
  .clean_chr_vec <- function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    x
  }

  test_chrs <- .clean_chr_vec(test_chrs)
  val_chrs  <- .clean_chr_vec(val_chrs)

  # No overlaps allowed (ambiguous assignment)
  if (length(intersect(test_chrs, val_chrs)) > 0L) {
    stop("Overlap between test_chrs and val_chrs is not allowed: ",
         paste(intersect(test_chrs, val_chrs), collapse = ", "))
  }

  # Normalize "chr" prefix style to match seqnames(gr)
  if (normalize_chr_prefix) {
    gr_has_chr <- any(grepl("^chr", seqs[!is.na(seqs)]))

    .norm_to_gr_style <- function(chrs) {
      if (!length(chrs)) return(chrs)

      # If any already match, keep as-is (common case)
      if (any(chrs %in% seqs, na.rm = TRUE)) return(chrs)

      ch_has_chr <- any(grepl("^chr", chrs))
      out <- chrs

      if (gr_has_chr && !ch_has_chr) {
        out <- paste0("chr", chrs)
      } else if (!gr_has_chr && ch_has_chr) {
        out <- sub("^chr", "", chrs)
      }

      # Small extra for mitochondrial naming differences (optional)
      # If user gives chrM but GR has MT (or vice versa), try to map.
      if (!any(out %in% seqs, na.rm = TRUE)) {
        out2 <- out
        out2[out2 == "chrM"] <- "MT"
        out2[out2 == "MT"] <- "chrM"
        # Only use the swap if it helps
        if (any(out2 %in% seqs, na.rm = TRUE)) out <- out2
      }

      out
    }

    test_chrs <- .norm_to_gr_style(test_chrs)
    val_chrs  <- .norm_to_gr_style(val_chrs)
  }

  # Warn if requested chromosomes are absent from the object
  if (warn_missing) {
    miss_test <- setdiff(test_chrs, unique(seqs))
    miss_val  <- setdiff(val_chrs,  unique(seqs))
    if (length(miss_test)) warning("These test_chrs are not present in seqnames(gr): ",
                                   paste(miss_test, collapse = ", "))
    if (length(miss_val)) warning("These val_chrs are not present in seqnames(gr): ",
                                  paste(miss_val, collapse = ", "))
  }

  split <- rep(train_label, length(gr))
  if (length(test_chrs)) split[seqs %in% test_chrs] <- test_label
  if (length(val_chrs))  split[seqs %in% val_chrs]  <- val_label

  # Factor with stable level order
  split <- factor(split, levels = c(train_label, test_label, val_label))

  mc[[split_col]] <- split
  S4Vectors::mcols(gr) <- mc
  gr
}



#' Write classification-ready HDF5 splits (train/val/test) from a single GRanges
#'
#' This function expects a GRanges that already contains both positives and negatives
#' (e.g. after subsetting to matched pairs) and already contains a split column
#' indicating train/val/test membership.
#'
#' For each site, it extracts a centered sequence window (2*flank+1 bp) from the
#' reference genome, pads out-of-bounds sequence with 'N', orients sequences by strand
#' (reverse-complement for '-' strand), and encodes bases as integers:
#' A=0, C=1, G=2, T=3, N=4. The encoded matrix is written as dataset "X_int", and labels
#' as dataset "y" in HDF5 files {train,val,test}.h5.
#'
#' @param gr A GRanges containing both classes (label 0/1) and a split column.
#' @param bsgenome A BSgenome object (or any object supported by Biostrings::getSeq()).
#' @param out_dir Output directory.
#' @param prefix Subdirectory name inside out_dir where files will be written.
#' @param flank Integer flank size; window width is 2*flank + 1 (default 500).
#' @param label_col Metadata column containing 0/1 class labels (default "label").
#' @param split_col Metadata column containing split labels (default "split").
#' @param id_col Metadata column containing unique IDs (default "id"). If missing, one is created.
#' @param chunk_n Chunk size for HDF5 writing (rows per write block).
#' @param compression_level gzip compression level for X_int (0-9; default 6).
#' @param overwrite If TRUE, overwrite existing .h5 files.
#' @param balance_train If TRUE, downsample the *training* split to balance labels 0/1.
#' @param seed Random seed (used for balancing only).
#' @param metadata_col metadata columns (from mcols(gr)) to write to CSV.
#'        If NULL, a sensible default set is attempted (only columns that exist are used).
#' @param drop_unknown_splits If TRUE, rows with split labels not mappable to train/val/test are dropped.
#' @param quiet If TRUE, suppress messages.
#'
#' @return Invisibly returns the output directory path (file.path(out_dir, prefix)).
#' @export
write_h5_classification_from_granges <- function(gr,
                                                 bsgenome,
                                                 out_dir  = "bpnet_data",
                                                 prefix   = "dataset_cls",
                                                 flank    = 500L,
                                                 label_col = "label",
                                                 split_col = "split",
                                                 id_col    = "id",
                                                 chunk_n   = 10000L,
                                                 compression_level = 6L,
                                                 overwrite = TRUE,
                                                 balance_train = FALSE,
                                                 seed = 1337L,
                                                 metadata_cols = NULL,
                                                 drop_unknown_splits = TRUE,
                                                 quiet = FALSE) {
  stopifnot(inherits(gr, "GRanges"))

  # ---- dependency guard (Suggests-friendly) ----
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' is required for HDF5 output. Install via BiocManager::install('rhdf5').")
  }

  mc <- S4Vectors::mcols(gr)

  if (!(label_col %in% colnames(mc))) stop("Missing label_col in mcols(gr): ", label_col)
  if (!(split_col %in% colnames(mc))) stop("Missing split_col in mcols(gr): ", split_col)

  # Stable rownames/IDs
  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_len(length(gr)))
  }

  # Ensure id column exists (for downstream metadata + python join keys)
  if (!(id_col %in% colnames(mc))) {
    mc[[id_col]] <- paste0(as.character(GenomicRanges::seqnames(gr)), "_",
                           GenomicRanges::start(gr), "_",
                           as.character(GenomicRanges::strand(gr)))
    S4Vectors::mcols(gr) <- mc
    mc <- S4Vectors::mcols(gr)
  }

  # Labels -> integer 0/1 (but don't require both classes per split)
  y <- mc[[label_col]]
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y))  y <- as.integer(as.character(y))
  y <- as.integer(y)
  if (any(!is.na(y) & !(y %in% c(0L, 1L)))) {
    stop("`", label_col, "` must be 0/1 (or coercible to 0/1).")
  }

  # ---- standardize split labels to {train,val,test} ----
  split_raw <- as.character(mc[[split_col]])
  split_std <- tolower(split_raw)
  split_std[split_std %in% c("training", "train")] <- "train"
  split_std[split_std %in% c("validation", "valid", "val")] <- "val"
  split_std[split_std %in% c("testing", "test")] <- "test"

  ok_split <- split_std %in% c("train", "val", "test")
  if (drop_unknown_splits) {
    if (!all(ok_split)) {
      gr <- gr[ok_split]
      split_std <- split_std[ok_split]
      y <- y[ok_split]
      mc <- S4Vectors::mcols(gr)
      if (!quiet) message("Dropped ", sum(!ok_split), " rows with unknown split labels in `", split_col, "`.")
    }
  } else {
    if (any(!ok_split)) {
      stop("Found split labels not mappable to {train,val,test} in `", split_col, "`.")
    }
  }

  if (length(gr) == 0L) stop("After split filtering, gr has 0 rows.")

  # ---- optional balancing within training split ----
  if (balance_train) {
    set.seed(seed)
    tr_idx <- which(split_std == "train")
    if (length(tr_idx)) {
      tr_pos <- tr_idx[y[tr_idx] == 1L]
      tr_neg <- tr_idx[y[tr_idx] == 0L]
      if (length(tr_pos) && length(tr_neg)) {
        n_min <- min(length(tr_pos), length(tr_neg))
        keep_train <- c(sample(tr_pos, n_min), sample(tr_neg, n_min))
        keep_all <- c(keep_train, setdiff(seq_len(length(gr)), tr_idx))
        keep_all <- sort(keep_all)
        gr <- gr[keep_all]
        split_std <- split_std[keep_all]
        y <- y[keep_all]
        mc <- S4Vectors::mcols(gr)
        if (!quiet) message("Balanced training split to ", n_min, " positives and ", n_min, " negatives.")
      }
    }
  }

  # ---- harmonize seqlevels against genome ----
  # Use genome style where possible, keep only common seqlevels, and set seqinfo.
  st <- GenomeInfoDb::seqlevelsStyle(bsgenome)
  if (length(st) > 0L) {
    suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- st[1], silent = TRUE))
  }
  gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
  common <- intersect(GenomeInfoDb::seqlevels(gr), GenomicRanges::seqnames(bsgenome))
  gr <- GenomeInfoDb::keepSeqlevels(gr, common, pruning.mode = "coarse")
  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(bsgenome)[GenomeInfoDb::seqlevels(gr)]

  # Also subset vectors after pruning
  idx_keep <- rep(TRUE, length(split_std))
  # Note: keepSeqlevels changes length(gr) but keeps order; we must realign by names
  # Safest is to re-pull split/y from mcols after pruning if split/label exist.
  mc <- S4Vectors::mcols(gr)
  y <- mc[[label_col]]
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y))  y <- as.integer(as.character(y))
  y <- as.integer(y)
  split_raw <- as.character(mc[[split_col]])
  split_std <- tolower(split_raw)
  split_std[split_std %in% c("training", "train")] <- "train"
  split_std[split_std %in% c("validation", "valid", "val")] <- "val"
  split_std[split_std %in% c("testing", "test")] <- "test"

  # ---- build windows + N padding + strand orientation ----
  flank <- as.integer(flank)
  if (!is.finite(flank) || flank < 0L) stop("`flank` must be a non-negative integer.")
  win_width <- as.integer(2L * flank + 1L)

  gr_win <- GenomicRanges::resize(gr, width = win_width, fix = "center", ignore.strand = FALSE)

  starts <- GenomicRanges::start(gr_win)
  ends   <- GenomicRanges::end(gr_win)

  seqlen <- GenomeInfoDb::seqlengths(gr_win)[as.character(GenomicRanges::seqnames(gr_win))]
  left_pad  <- pmax(1L - starts, 0L)
  right_pad <- pmax(ends - seqlen, 0L)

  gr_trim <- GenomicRanges::trim(gr_win)

  seqs_core <- Biostrings::getSeq(bsgenome, gr_trim)

  # Cache N pads up to flank for speed
  pad_cache <- vector("list", flank + 1L)
  pad_cache[[1]] <- Biostrings::DNAString("")
  if (flank > 0L) {
    for (i in seq_len(flank)) {
      pad_cache[[i + 1L]] <- Biostrings::DNAString(paste(rep("N", i), collapse = ""))
    }
  }
  padN <- function(n) pad_cache[[as.integer(n) + 1L]]

  seqs_full <- Biostrings::DNAStringSet(
    mapply(
      function(s, lp, rp) Biostrings::xscat(padN(lp), s, padN(rp)),
      seqs_core, as.integer(left_pad), as.integer(right_pad),
      SIMPLIFY = FALSE
    )
  )
  if (!all(Biostrings::width(seqs_full) == win_width)) {
    stop("Internal error: not all sequences are the expected width (", win_width, ").")
  }

  minus <- as.character(GenomicRanges::strand(gr_win)) == "-"
  if (any(minus, na.rm = TRUE)) {
    seqs_full[minus] <- Biostrings::reverseComplement(seqs_full[minus])
  }

  # ---- encoder: A/C/G/T/N -> 0/1/2/3/4 ----
  encode_block <- function(seqs_chr) {
    s <- toupper(seqs_chr)
    s <- chartr("U", "T", s)
    s <- gsub("[^ACGTN]", "N", s)     # squash any IUPAC ambiguity to N
    s <- chartr("ACGTN", "01234", s)

    B <- length(s)
    L <- nchar(s[1L])
    out <- matrix(NA_integer_, nrow = B, ncol = L)

    z <- utf8ToInt("0")  # 48
    for (i in seq_len(B)) {
      out[i, ] <- utf8ToInt(s[i]) - z
    }
    out
  }

  # ---- metadata column selection ----
  if (is.null(metadata_cols)) {
    # conservative defaults: scalar-ish columns that are often useful
    candidate <- c("motif", "kmer", "location", "feature", "tx_name", "gene_id",
                   "gene_symbol", "metagene_prop",
                   "feature_width", "segment_rank",
                   "dist_from_feature_start", "dist_from_feature_end",
                   "start_dist", "stop_dist")
    metadata_cols <- intersect(candidate, colnames(mc))
  } else {
    metadata_cols <- intersect(as.character(metadata_cols), colnames(mc))
  }

  coerce_scalar <- function(x) {
    if (inherits(x, "DNAStringSet")) return(as.character(x))
    if (!is.null(dim(x))) return(NULL)  # drop matrices/arrays
    if (inherits(x, "List") || is.list(x)) {
      # first element (or NA) per row
      return(vapply(x, function(z) if (length(z)) as.character(z[[1]]) else NA_character_, character(1)))
    }
    # atomic vector
    if (is.factor(x)) x <- as.character(x)
    x
  }

  # ---- per-split writer ----
  out_dir_split <- file.path(out_dir, prefix)
  dir.create(out_dir_split, showWarnings = FALSE, recursive = TRUE)

  write_one_split <- function(split_name) {
    idx <- which(split_std == split_name)
    if (!length(idx)) return(invisible(NULL))

    h5file <- file.path(out_dir_split, paste0(split_name, ".h5"))
    if (file.exists(h5file)) {
      if (overwrite) file.remove(h5file) else stop("File exists and overwrite=FALSE: ", h5file)
    }

    seqs_split <- seqs_full[idx]
    y_split <- as.integer(y[idx])
    ids_split <- as.character(mc[[id_col]][idx])

    N <- length(seqs_split)
    L <- Biostrings::width(seqs_split)[1L]

    rhdf5::h5createFile(h5file)
    rhdf5::h5createDataset(
      h5file, "X_int",
      dims = c(N, L),
      storage.mode = "integer",
      chunk = c(min(N, as.integer(chunk_n)), L),
      level = as.integer(compression_level)
    )
    rhdf5::h5createDataset(h5file, "y", dims = N, storage.mode = "integer")
    rhdf5::h5write(y_split, h5file, "y")

    # Chunked write
    blocks <- split(seq_len(N), ceiling(seq_len(N) / as.integer(chunk_n)))
    for (b in seq_along(blocks)) {
      rows <- blocks[[b]]
      Xblk <- encode_block(as.character(seqs_split[rows]))
      rhdf5::h5write(Xblk, h5file, "X_int", index = list(rows, seq_len(L)))
    }

    # Metadata CSV (kept simple + robust)
    meta <- data.frame(
      id       = ids_split,
      split    = split_name,
      seqnames = as.character(GenomicRanges::seqnames(gr)[idx]),
      start    = GenomicRanges::start(gr)[idx],
      end      = GenomicRanges::end(gr)[idx],
      strand   = as.character(GenomicRanges::strand(gr)[idx]),
      label    = y_split,
      stringsAsFactors = FALSE
    )

    if (length(metadata_cols)) {
      for (nm in metadata_cols) {
        v <- coerce_scalar(mc[[nm]][idx])
        if (!is.null(v)) meta[[nm]] <- v
      }
    }

    meta_path <- file.path(out_dir_split, paste0(split_name, "_metadata.csv"))
    utils::write.csv(meta, meta_path, row.names = FALSE)

    rhdf5::H5close()

    if (!quiet) {
      message(sprintf("Wrote %s: X_int[%d x %d], y[%d] -> %s",
                      split_name, N, L, N, h5file))
    }
    invisible(h5file)
  }

  write_one_split("train")
  write_one_split("val")
  write_one_split("test")

  invisible(out_dir_split)
}


#' Internal helper for setting names
#' @noRd
.site_id <- function(gr,strand=TRUE){
  if(strand){
    paste0(as.character(seqnames(gr)), ":", start(gr), ":", as.character(strand(gr)))
  }else{
    paste0(as.character(seqnames(gr)), ":", start(gr))
  }
}


#' Add site id according to genomic coordinates to a GRange
#'
#' @param strand whether or not to include the strand as part of the id
#'
#' @return the same GRange with the names set
#' @export
makeID <- function(gr,strand=TRUE) {
  names(gr) <- .site_id(gr,strand)
  gr
}




#' Balance matched positive/negative pairs by m6A ratio bins
#'
#' Downsamples *positive* sites (label==1) to equalize their distribution across
#' user-defined ratio bins (e.g. 0-0.2, 0.2-0.4, ...), while preserving matched pairs:
#' for every retained positive, its matched negative is kept as well.
#'
#' By default, this is intended to be applied to the *training* split only,
#' leaving validation/testing untouched to preserve realistic evaluation.
#'
#' @param gr A \code{GRanges} containing positives and negatives, plus match columns.
#' @param ratio_col Metadata column containing the ratio in [0,1] (e.g. "Ratio").
#' @param label_col Metadata column containing 0/1 labels (default "label").
#' @param matched_negative_id_col Column giving the matched negative ID for each positive.
#' @param matched_positive_id_col Column giving the matched positive ID for each matched negative.
#'   Only required if \code{strict_reciprocal=TRUE}.
#' @param split_col Optional split column (e.g. "dataset_split"). If provided, only rows
#'   in \code{balance_splits} are resampled; other splits are kept unchanged.
#' @param balance_splits Which split label(s) to balance (default "training").
#' @param breaks Optional numeric breaks for binning the ratio. Default makes 0.2-width bins.
#' @param bin_width Optional alternative to \code{breaks}. If provided, builds breaks as
#'   \code{seq(0, 1, by=bin_width)} (ensuring 1 is included).
#' @param clamp_ratio If TRUE, clamp ratio values to [0,1] before binning.
#' @param missing_ratio What to do with positives that have missing/non-finite ratio:
#'   "drop" removes them from the balanced subset; "keep" keeps them (unbalanced).
#' @param target_n_per_bin Optional integer. If NULL (default), uses the minimum count
#'   among non-empty bins (pure downsampling, perfectly balanced across non-empty bins).
#'   If provided and larger than that minimum, it is reduced to the minimum (no upsampling).
#' @param strict_reciprocal If TRUE, require the negative to point back to the same positive.
#' @param drop_conflicts If TRUE, enforce that a negative is used at most once (drop duplicates).
#' @param require_same_split If TRUE and \code{split_col} is provided, require that each kept
#'   positive/negative pair lies in the same split label.
#' @param seed Random seed for reproducible downsampling.
#' @param return_diagnostics If TRUE, attach a \code{"balance_ratio_diagnostics"} attribute.
#'
#' @return A subset \code{GRanges} with balanced positives (in the chosen split(s))
#'   plus their matched negatives, and all untouched rows from other splits (if any).
#' @export
balance_pairs_by_ratio <- function(gr,
                                   ratio_col = "Ratio",
                                   label_col = "label",
                                   matched_negative_id_col = "matched_negative_id",
                                   matched_positive_id_col = "matched_positive_id",
                                   split_col = NULL,
                                   balance_splits = "training",
                                   breaks = seq(0, 1, by = 0.2),
                                   bin_width = NULL,
                                   clamp_ratio = TRUE,
                                   missing_ratio = c("drop", "keep"),
                                   target_n_per_bin = NULL,
                                   strict_reciprocal = TRUE,
                                   drop_conflicts = TRUE,
                                   require_same_split = TRUE,
                                   seed = 1L,
                                   return_diagnostics = TRUE) {
  stopifnot(inherits(gr, "GRanges"))
  missing_ratio <- match.arg(missing_ratio)

  mc <- S4Vectors::mcols(gr)

  if (!(ratio_col %in% colnames(mc))) stop("Missing ratio_col in mcols(gr): ", ratio_col)
  if (!(label_col %in% colnames(mc))) stop("Missing label_col in mcols(gr): ", label_col)
  if (!(matched_negative_id_col %in% colnames(mc))) {
    stop("Missing `", matched_negative_id_col, "` in mcols(gr). Need match_background()/random_drach_within_transcript() output.")
  }
  if (strict_reciprocal && !(matched_positive_id_col %in% colnames(mc))) {
    stop("strict_reciprocal=TRUE but missing `", matched_positive_id_col, "` in mcols(gr).")
  }

  # Stable IDs
  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_len(length(gr)))
  }
  ids <- names(gr)

  # Parse splits
  all_idx <- seq_along(gr)
  balance_idx <- all_idx
  keep_untouched_idx <- integer(0)
  split_vec <- NULL

  if (!is.null(split_col)) {
    if (!(split_col %in% colnames(mc))) stop("split_col not found in mcols(gr): ", split_col)
    split_vec <- as.character(mc[[split_col]])
    balance_splits <- as.character(balance_splits)

    balance_idx <- which(!is.na(split_vec) & split_vec %in% balance_splits)
    keep_untouched_idx <- setdiff(all_idx, balance_idx)

    if (!length(balance_idx)) {
      if (return_diagnostics) {
        attr(gr, "balance_ratio_diagnostics") <- list(
          n_total_in = length(gr),
          split_col = split_col,
          balance_splits = balance_splits,
          n_in_balance_splits = 0L,
          reason = "No rows in requested balance_splits; returning input unchanged."
        )
      }
      return(gr)
    }
  }

  # Labels -> integer 0/1
  y <- mc[[label_col]]
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y))  y <- as.integer(as.character(y))
  y <- as.integer(y)

  pos_idx0 <- intersect(which(y == 1L), balance_idx)
  if (!length(pos_idx0)) {
    out <- gr[sort(unique(c(keep_untouched_idx)))]
    if (return_diagnostics) {
      attr(out, "balance_ratio_diagnostics") <- list(
        n_total_in = length(gr),
        n_total_out = length(out),
        n_pos_in_balance_splits = 0L,
        reason = "No positives in balance_splits."
      )
    }
    return(out)
  }

  # Candidate pairs: positive must have a valid matched negative that exists
  neg_ids0 <- as.character(mc[[matched_negative_id_col]][pos_idx0])
  ok <- !is.na(neg_ids0) & nzchar(neg_ids0)
  pos_idx <- pos_idx0[ok]
  neg_ids <- neg_ids0[ok]

  neg_idx <- match(neg_ids, ids)
  ok2 <- !is.na(neg_idx)
  pos_idx <- pos_idx[ok2]
  neg_idx <- neg_idx[ok2]
  neg_ids <- neg_ids[ok2]

  # Optionally require same split (prevents cross-split contamination)
  if (require_same_split && !is.null(split_vec)) {
    sp_pos <- split_vec[pos_idx]
    sp_neg <- split_vec[neg_idx]
    okS <- !is.na(sp_pos) & !is.na(sp_neg) & (sp_pos == sp_neg)
    pos_idx <- pos_idx[okS]
    neg_idx <- neg_idx[okS]
    neg_ids <- neg_ids[okS]
  }

  # Optionally enforce reciprocal link
  if (strict_reciprocal) {
    back <- as.character(mc[[matched_positive_id_col]][neg_idx])
    ok3 <- !is.na(back) & nzchar(back) & (back == ids[pos_idx])
    pos_idx <- pos_idx[ok3]
    neg_idx <- neg_idx[ok3]
    neg_ids <- neg_ids[ok3]
  }

  if (!length(pos_idx)) {
    out <- gr[sort(unique(c(keep_untouched_idx)))]
    if (return_diagnostics) {
      attr(out, "balance_ratio_diagnostics") <- list(
        n_total_in = length(gr),
        n_total_out = length(out),
        n_pos_in_balance_splits = length(pos_idx0),
        n_pairs_candidate = 0L,
        reason = "No valid matched pairs among positives in balance_splits after filtering."
      )
    }
    return(out)
  }

  # Enforce unique negatives if requested
  if (drop_conflicts) {
    set.seed(seed)
    o <- sample(seq_along(pos_idx))  # random tie-breaking among duplicates
    pos_idx <- pos_idx[o]; neg_idx <- neg_idx[o]; neg_ids <- neg_ids[o]

    dup_neg <- duplicated(neg_ids)
    if (any(dup_neg)) {
      keep <- !dup_neg
      pos_idx <- pos_idx[keep]
      neg_idx <- neg_idx[keep]
      neg_ids <- neg_ids[keep]
    }
  }

  # Ratio vector for candidate positives
  r <- suppressWarnings(as.numeric(mc[[ratio_col]][pos_idx]))
  if (clamp_ratio) r <- pmin(1, pmax(0, r))

  okR <- is.finite(r)
  pos_missing_ratio <- pos_idx[!okR]
  neg_missing_ratio <- neg_idx[!okR]

  pos_idxR <- pos_idx[okR]
  neg_idxR <- neg_idx[okR]
  rR <- r[okR]

  # Build breaks from bin_width if provided
  if (!is.null(bin_width)) {
    bw <- as.numeric(bin_width)
    if (!is.finite(bw) || bw <= 0) stop("bin_width must be > 0.")
    breaks <- seq(0, 1, by = bw)
    if (tail(breaks, 1) < 1) breaks <- c(breaks, 1)
  }

  breaks <- as.numeric(breaks)
  breaks <- sort(unique(breaks))
  if (length(breaks) < 2) stop("breaks must have at least 2 unique values.")
  if (breaks[1] > 0) breaks <- c(0, breaks)
  if (tail(breaks, 1) < 1) breaks <- c(breaks, 1)

  # Bin positives by ratio
  bin <- cut(rR, breaks = breaks, include.lowest = TRUE, right = TRUE)
  bin_chr <- as.character(bin)
  bin_counts <- table(bin_chr, useNA = "no")

  if (!length(bin_counts)) {
    # No finite ratio positives
    keep_balance_idx <- integer(0)
  } else {
    min_nonzero <- min(as.integer(bin_counts))

    target <- min_nonzero
    if (!is.null(target_n_per_bin)) {
      target_n_per_bin <- as.integer(target_n_per_bin)
      if (!is.finite(target_n_per_bin) || target_n_per_bin < 1L) {
        stop("target_n_per_bin must be a positive integer.")
      }
      if (target_n_per_bin > min_nonzero) {
        warning("target_n_per_bin > min nonzero bin count (", min_nonzero,
                "); reducing to ", min_nonzero, " (no upsampling).")
        target <- min_nonzero
      } else {
        target <- target_n_per_bin
      }
    }

    # Sample equal number per non-empty bin
    set.seed(seed)
    pos_keep <- integer(0)
    for (bn in names(bin_counts)) {
      idx_bn <- which(bin_chr == bn)
      pos_bn <- pos_idxR[idx_bn]
      if (length(pos_bn) >= target) {
        pos_keep <- c(pos_keep, sample(pos_bn, size = target, replace = FALSE))
      }
    }

    # Get corresponding negatives via mapping pos_idxR -> neg_idxR
    sel <- match(pos_keep, pos_idxR)
    neg_keep <- neg_idxR[sel]

    keep_balance_idx <- unique(c(pos_keep, neg_keep))

    # Optionally keep missing-ratio positives (and their negatives) unbalanced
    if (missing_ratio == "keep" && length(pos_missing_ratio)) {
      keep_balance_idx <- unique(c(keep_balance_idx, pos_missing_ratio, neg_missing_ratio))
    }
  }

  # Final output:
  # - keep all rows NOT in balance_splits unchanged
  # - within balance_splits keep only balanced matched pairs (and optionally missing_ratio kept)
  out_idx <- sort(unique(c(keep_untouched_idx, keep_balance_idx)))
  out <- gr[out_idx]

  if (return_diagnostics) {
    attr(out, "balance_ratio_diagnostics") <- list(
      n_total_in = length(gr),
      n_total_out = length(out),
      split_col = split_col,
      balance_splits = balance_splits,
      n_in_balance_splits_in = length(balance_idx),
      n_pos_in_balance_splits_in = length(pos_idx0),
      n_candidate_pairs = length(pos_idx),
      breaks = breaks,
      bin_counts_candidates = as.list(bin_counts),
      target_n_per_bin = if (exists("target")) target else NA_integer_,
      missing_ratio = missing_ratio,
      strict_reciprocal = strict_reciprocal,
      drop_conflicts = drop_conflicts,
      require_same_split = require_same_split
    )
  }

  out
}
