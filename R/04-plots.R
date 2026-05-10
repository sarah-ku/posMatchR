.auto_metagene_col <- function(gr, x_col = NULL) {
  mc <- S4Vectors::mcols(gr)
  if (!is.null(x_col)) {
    if (!(x_col %in% colnames(mc))) stop("x_col not found in mcols(gr): ", x_col)
    return(x_col)
  }
  for (nm in c("metagene_split3", "metagene_prop", "metagene_split")) {
    if (nm %in% colnames(mc)) return(nm)
  }
  stop("No metagene column found. Expected metagene_split3, metagene_prop, or metagene_split.")
}

.make_plot_df <- function(gr,
                          cols,
                          set_col = NULL,
                          positive_value = "positive",
                          negative_value = "matched_negative",
                          include_other = FALSE) {
  mc <- S4Vectors::mcols(gr)
  miss <- setdiff(cols, colnames(mc))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  df <- data.frame(.posmatchr_row = seq_len(length(gr)), stringsAsFactors = FALSE)
  for (cl in cols) df[[cl]] <- mc[[cl]]
  df$.posmatchr_row <- NULL

  if (is.null(set_col) && "match_set" %in% colnames(mc)) set_col <- "match_set"

  if (!is.null(set_col) && set_col %in% colnames(mc)) {
    s <- as.character(mc[[set_col]])
    df$group <- ifelse(
      s == positive_value,
      "positive",
      ifelse(s == negative_value, "negative", ifelse(include_other, "other", NA_character_))
    )
  } else {
    df$group <- "sites"
  }

  df <- df[!is.na(df$group), , drop = FALSE]
  df$group <- factor(df$group, levels = unique(c("positive", "negative", "other", "sites", df$group)))
  df
}

.metagene_axis_settings <- function(gr, x_col) {
  md <- S4Vectors::metadata(gr)
  breaks3 <- if (identical(x_col, "metagene_prop")) c(1, 2) else md$metagene_breaks3
  if (is.null(breaks3) || length(breaks3) < 2L || any(!is.finite(as.numeric(breaks3[1:2])))) {
    breaks3 <- c(1, 2)
  } else {
    breaks3 <- as.numeric(breaks3[1:2])
  }
  centres <- c(breaks3[1] / 2, mean(breaks3), breaks3[2] + (3 - breaks3[2]) / 2)
  list(
    breaks = breaks3,
    centres = centres,
    labels = c("5' UTR", "CDS", "3' UTR")
  )
}

.pretty_bp_axis <- function(max_bp, breaks_bp = c(0, 1, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)) {
  max_bp <- suppressWarnings(as.numeric(max_bp))
  if (!is.finite(max_bp) || max_bp < 0) max_bp <- 0
  br <- breaks_bp[breaks_bp <= max_bp]
  if (!length(br)) br <- breaks_bp[1]
  unique(br)
}

#' Plot metagene density
#'
#' Plots a density of annotated sites along a 0--3 metagene axis. If a match-set column
#' is present, positives and matched negatives are shown separately; otherwise all sites
#' are shown as one group.
#'
#' @param gr Annotated \code{GRanges}.
#' @param x_col Metagene column to use. Defaults to \code{metagene_split3}, then
#'   \code{metagene_prop}, then \code{metagene_split}.
#' @param set_col Optional grouping column. Defaults to \code{"match_set"} when present.
#' @param positive_value Value in \code{set_col} marking positives.
#' @param negative_value Value in \code{set_col} marking matched negatives.
#' @param facet_by_location If TRUE, facet by region.
#' @param location_col Region column.
#' @param bw_adjust Density bandwidth multiplier.
#' @param xlim X-axis limits.
#' @param x_axis Either \code{"regions"} to label the x-axis as 5'UTR/CDS/3'UTR
#'   or \code{"numeric"} to show the numeric metagene coordinate.
#' @param x_label Optional x-axis title.
#'
#' @return A ggplot object.
#' @export
plot_metagene_density <- function(gr,
                                  x_col = NULL,
                                  set_col = NULL,
                                  positive_value = "positive",
                                  negative_value = "matched_negative",
                                  facet_by_location = FALSE,
                                  location_col = "location",
                                  bw_adjust = 1,
                                  xlim = c(0, 3),
                                  x_axis = c("regions", "numeric"),
                                  x_label = NULL) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  x_axis <- match.arg(x_axis)
  x_col <- .auto_metagene_col(gr, x_col)
  cols <- c(x_col, if (facet_by_location) location_col else NULL)
  df <- .make_plot_df(gr, cols = cols, set_col = set_col,
                      positive_value = positive_value, negative_value = negative_value)
  df$x_plot <- suppressWarnings(as.numeric(df[[x_col]]))
  if (identical(x_col, "metagene_split")) df$x_plot <- df$x_plot * 3
  df <- df[is.finite(df$x_plot), , drop = FALSE]
  if (!nrow(df)) stop("No finite values in ", x_col)

  axis_info <- .metagene_axis_settings(gr, x_col)
  xlab <- x_label %||% if (x_axis == "regions") {
    "Transcript region"
  } else {
    paste0(x_col, " (0-3 metagene coordinate)")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x_plot, colour = group)) +
    ggplot2::geom_density(linewidth = 1, adjust = bw_adjust, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = axis_info$breaks, linetype = 2) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = xlab, y = "Density", colour = NULL) +
    ggplot2::theme_bw()

  if (x_axis == "regions") {
    p <- p + ggplot2::scale_x_continuous(breaks = axis_info$centres, labels = axis_info$labels, minor_breaks = NULL)
  }

  if (facet_by_location) {
    df$facet_location <- as.character(df[[location_col]])
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x_plot, colour = group)) +
      ggplot2::geom_density(linewidth = 1, adjust = bw_adjust, na.rm = TRUE) +
      ggplot2::geom_vline(xintercept = axis_info$breaks, linetype = 2) +
      ggplot2::coord_cartesian(xlim = xlim) +
      ggplot2::facet_wrap(~ facet_location) +
      ggplot2::labs(x = xlab, y = "Density", colour = NULL) +
      ggplot2::theme_bw()
    if (x_axis == "regions") {
      p <- p + ggplot2::scale_x_continuous(breaks = axis_info$centres, labels = axis_info$labels, minor_breaks = NULL)
    }
  }

  p
}

#' Plot a generic distance density
#'
#' @param gr A \code{GRanges}.
#' @param distance_col Numeric metadata column to plot.
#' @param set_col Optional grouping column.
#' @param positive_value Positive set value.
#' @param negative_value Matched-negative set value.
#' @param transform Either \code{"log1p"} or \code{"identity"}.
#' @param facet_by_location If TRUE, facet by \code{location_col}.
#' @param location_col Location column.
#' @param bw_adjust Density bandwidth multiplier.
#' @param x_label Optional x-axis title.
#' @param breaks_bp Tick marks to show on the original bp scale when \code{transform = "log1p"}.
#'
#' @return A ggplot object.
#' @export
plot_distance_density <- function(gr,
                                  distance_col,
                                  set_col = NULL,
                                  positive_value = "positive",
                                  negative_value = "matched_negative",
                                  transform = c("log1p", "identity"),
                                  facet_by_location = TRUE,
                                  location_col = "location",
                                  bw_adjust = 1,
                                  x_label = NULL,
                                  breaks_bp = c(0, 1, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000)) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  transform <- match.arg(transform)

  cols <- c(distance_col, if (facet_by_location) location_col else NULL)
  df <- .make_plot_df(gr, cols = cols, set_col = set_col,
                      positive_value = positive_value, negative_value = negative_value)
  d <- suppressWarnings(as.numeric(df[[distance_col]]))
  df$value <- d
  df <- df[is.finite(df$value), , drop = FALSE]
  if (!nrow(df)) stop("No finite values in ", distance_col)

  raw_value <- df$value
  if (transform == "log1p") {
    df$value <- log1p(pmax(df$value, 0))
    xlab <- x_label %||% paste0(distance_col, " (bp; log scale)")
    bp_breaks <- .pretty_bp_axis(max(raw_value, na.rm = TRUE), breaks_bp = breaks_bp)
    bp_labels <- format(bp_breaks, big.mark = ",", scientific = FALSE, trim = TRUE)
  } else {
    xlab <- x_label %||% paste0(distance_col, " (bp)")
    bp_breaks <- NULL
    bp_labels <- NULL
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, colour = group)) +
    ggplot2::geom_density(linewidth = 1, adjust = bw_adjust, na.rm = TRUE) +
    ggplot2::labs(x = xlab, y = "Density", colour = NULL) +
    ggplot2::theme_bw()

  if (transform == "log1p") {
    p <- p + ggplot2::scale_x_continuous(breaks = log1p(bp_breaks), labels = bp_labels)
  }

  if (facet_by_location) {
    df$facet_location <- as.character(df[[location_col]])
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, colour = group)) +
      ggplot2::geom_density(linewidth = 1, adjust = bw_adjust, na.rm = TRUE) +
      ggplot2::facet_wrap(~ facet_location) +
      ggplot2::labs(x = xlab, y = "Density", colour = NULL) +
      ggplot2::theme_bw()
    if (transform == "log1p") {
      p <- p + ggplot2::scale_x_continuous(breaks = log1p(bp_breaks), labels = bp_labels)
    }
  }

  p
}

#' Plot distance to nearest exon junction
#'
#' @param gr A \code{GRanges}.
#' @param distance_col Junction distance column.
#' @param ... Passed to \code{\link{plot_distance_density}}.
#'
#' @return A ggplot object.
#' @export
plot_junction_distance_density <- function(gr,
                                           distance_col = "nearest_exon_junction_dist",
                                           ...) {
  args <- list(...)
  if (!("x_label" %in% names(args))) {
    args$x_label <- "Distance to nearest exon junction (bp; log scale)"
  }
  do.call(plot_distance_density, c(list(gr = gr, distance_col = distance_col), args))
}

#' Plot splice-proximity density
#'
#' Uses \code{nearest_exon_junction_dist} when present. Otherwise it falls back to
#' \code{min(dist_from_feature_start, dist_from_feature_end)}.
#'
#' @param gr A \code{GRanges}.
#' @param set_col Optional grouping column.
#' @param positive_value Positive value.
#' @param negative_value Negative value.
#' @param start_col Distance-from-feature-start column.
#' @param end_col Distance-from-feature-end column.
#' @param transform Either \code{"log1p"} or \code{"identity"}.
#' @param facet_by_location If TRUE, facet by location.
#' @param location_col Location column.
#' @param bw_adjust Bandwidth multiplier.
#'
#' @return A ggplot object.
#' @export
plot_splice_distance_density <- function(gr,
                                         set_col = NULL,
                                         positive_value = "positive",
                                         negative_value = "matched_negative",
                                         start_col = "dist_from_feature_start",
                                         end_col = "dist_from_feature_end",
                                         transform = c("log1p", "identity"),
                                         facet_by_location = TRUE,
                                         location_col = "location",
                                         bw_adjust = 1) {
  if ("nearest_exon_junction_dist" %in% colnames(S4Vectors::mcols(gr))) {
    return(plot_distance_density(
      gr,
      distance_col = "nearest_exon_junction_dist",
      set_col = set_col,
      positive_value = positive_value,
      negative_value = negative_value,
      transform = transform,
      facet_by_location = facet_by_location,
      location_col = location_col,
      bw_adjust = bw_adjust
    ))
  }

  transform <- match.arg(transform)
  cols <- c(start_col, end_col, if (facet_by_location) location_col else NULL)
  df <- .make_plot_df(gr, cols = cols, set_col = set_col,
                      positive_value = positive_value, negative_value = negative_value)
  a <- suppressWarnings(as.numeric(df[[start_col]]))
  b <- suppressWarnings(as.numeric(df[[end_col]]))
  df$value <- pmin(a, b)
  df <- df[is.finite(df$value), , drop = FALSE]
  if (!nrow(df)) stop("No finite splice/boundary distances.")

  if (transform == "log1p") {
    df$value <- log1p(pmax(df$value, 0))
    xlab <- paste0("log1p(min(", start_col, ", ", end_col, "))")
  } else {
    xlab <- paste0("min(", start_col, ", ", end_col, ")")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, colour = group)) +
    ggplot2::geom_density(linewidth = 1, adjust = bw_adjust, na.rm = TRUE) +
    ggplot2::labs(x = xlab, y = "Density", colour = NULL) +
    ggplot2::theme_bw()

  if (facet_by_location) {
    df$facet_location <- as.character(df[[location_col]])
    p <- p + ggplot2::facet_wrap(~ facet_location)
  }

  p
}

#' Plot distances to CDS start and stop anchors
#'
#' @param gr A \code{GRanges}.
#' @param set_col Optional grouping column.
#' @param positive_value Positive value.
#' @param negative_value Negative value.
#' @param start_dist_col Start-anchor distance column. Defaults to transcript-coordinate distance when available.
#' @param stop_dist_col Stop-anchor distance column. Defaults to transcript-coordinate distance when available.
#' @param transform Either \code{"log1p"} or \code{"identity"}.
#' @param bw_adjust Density bandwidth multiplier.
#'
#' @return A ggplot object.
#' @export
plot_start_stop_distance_density <- function(gr,
                                             set_col = NULL,
                                             positive_value = "positive",
                                             negative_value = "matched_negative",
                                             start_dist_col = NULL,
                                             stop_dist_col = NULL,
                                             transform = c("log1p", "identity"),
                                             bw_adjust = 1) {
  transform <- match.arg(transform)
  mc_names <- colnames(S4Vectors::mcols(gr))
  if (is.null(start_dist_col)) {
    start_dist_col <- if ("start_dist_tx" %in% mc_names) "start_dist_tx" else "start_dist"
  }
  if (is.null(stop_dist_col)) {
    stop_dist_col <- if ("stop_dist_tx" %in% mc_names) "stop_dist_tx" else "stop_dist"
  }
  df <- .make_plot_df(gr, cols = c(start_dist_col, stop_dist_col), set_col = set_col,
                      positive_value = positive_value, negative_value = negative_value)

  long <- rbind(
    data.frame(group = df$group, which_distance = "start", value = suppressWarnings(as.numeric(df[[start_dist_col]]))),
    data.frame(group = df$group, which_distance = "stop", value = suppressWarnings(as.numeric(df[[stop_dist_col]])))
  )
  long <- long[is.finite(long$value), , drop = FALSE]
  if (!nrow(long)) stop("No finite start/stop distances.")

  if (transform == "log1p") {
    long$value <- log1p(pmax(long$value, 0))
    xlab <- "log1p(distance)"
  } else {
    xlab <- "distance"
  }

  ggplot2::ggplot(long, ggplot2::aes(x = value, colour = group)) +
    ggplot2::geom_density(linewidth = 1, adjust = bw_adjust, na.rm = TRUE) +
    ggplot2::facet_wrap(~ which_distance, scales = "free_x") +
    ggplot2::labs(x = xlab, y = "Density", colour = NULL) +
    ggplot2::theme_bw()
}

#' Plot feature-width quantile bins
#'
#' @param gr A \code{GRanges}.
#' @param set_col Optional grouping column.
#' @param positive_value Positive value.
#' @param negative_value Negative value.
#' @param width_col Width column.
#' @param n_bins Number of quantile bins.
#' @param log1p_first If TRUE, bin log1p(width).
#' @param facet_by_location If TRUE, facet by location.
#' @param location_col Location column.
#'
#' @return A ggplot object.
#' @export
plot_feature_width_bins <- function(gr,
                                    set_col = NULL,
                                    positive_value = "positive",
                                    negative_value = "matched_negative",
                                    width_col = "feature_width",
                                    n_bins = 10,
                                    log1p_first = TRUE,
                                    facet_by_location = TRUE,
                                    location_col = "location") {
  cols <- c(width_col, if (facet_by_location) location_col else NULL)
  df <- .make_plot_df(gr, cols = cols, set_col = set_col,
                      positive_value = positive_value, negative_value = negative_value)
  v <- suppressWarnings(as.numeric(df[[width_col]]))
  if (log1p_first) v <- log1p(pmax(v, 0))
  ok <- is.finite(v)
  df <- df[ok, , drop = FALSE]
  v <- v[ok]
  if (!length(v)) stop("No finite values in ", width_col)

  br <- stats::quantile(v, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE, type = 7)
  br <- unique(as.numeric(br))
  if (length(br) < 2L) {
    df$bin <- factor("all")
  } else {
    br[1] <- -Inf
    br[length(br)] <- Inf
    df$bin <- factor(cut(v, breaks = br, include.lowest = TRUE, labels = FALSE))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = bin, fill = group)) +
    ggplot2::geom_bar(position = "dodge") +
    ggplot2::labs(
      x = if (log1p_first) paste0("Quantile bins of log1p(", width_col, ")") else paste0("Quantile bins of ", width_col),
      y = "Count",
      fill = NULL
    ) +
    ggplot2::theme_bw()

  if (facet_by_location) {
    df$facet_location <- as.character(df[[location_col]])
    p <- ggplot2::ggplot(df, ggplot2::aes(x = bin, fill = group)) +
      ggplot2::geom_bar(position = "dodge") +
      ggplot2::facet_wrap(~ facet_location) +
      ggplot2::labs(
        x = if (log1p_first) paste0("Quantile bins of log1p(", width_col, ")") else paste0("Quantile bins of ", width_col),
        y = "Count",
        fill = NULL
      ) +
      ggplot2::theme_bw()
  }

  p
}

#' Plot k-mer counts
#'
#' @param gr A \code{GRanges}.
#' @param set_col Optional grouping column.
#' @param positive_value Positive value.
#' @param negative_value Negative value.
#' @param kmer_col K-mer column.
#' @param top_n Number of top k-mers to show.
#'
#' @return A ggplot object.
#' @export
plot_kmer_counts <- function(gr,
                             set_col = NULL,
                             positive_value = "positive",
                             negative_value = "matched_negative",
                             kmer_col = "kmer",
                             top_n = 20) {
  df <- .make_plot_df(gr, cols = c(kmer_col), set_col = set_col,
                      positive_value = positive_value, negative_value = negative_value)
  df$kmer <- toupper(as.character(df[[kmer_col]]))
  df <- df[!is.na(df$kmer) & nzchar(df$kmer), , drop = FALSE]
  if (!nrow(df)) stop("No non-missing k-mers in ", kmer_col)

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

.matching_status <- function(gr,
                             label_col = "label",
                             set_col = "match_set",
                             positive_value = "positive",
                             negative_value = "matched_negative") {
  mc <- S4Vectors::mcols(gr)
  if (!(label_col %in% colnames(mc))) stop("Missing label_col: ", label_col)
  y <- .normalise_binary_label(mc[[label_col]], strict = TRUE)
  status <- ifelse(y == 1L, "unmatched_positive", "unselected_negative")
  if (!is.null(set_col) && set_col %in% colnames(mc)) {
    s <- as.character(mc[[set_col]])
    status[s == positive_value] <- "matched_positive"
    status[s == negative_value] <- "matched_negative"
  } else {
    if ("is_matched_positive" %in% colnames(mc)) {
      f <- .normalise_binary_label(mc$is_matched_positive, strict = FALSE)
      status[!is.na(f) & f == 1L] <- "matched_positive"
    }
    if ("is_matched_negative" %in% colnames(mc)) {
      f <- .normalise_binary_label(mc$is_matched_negative, strict = FALSE)
      status[!is.na(f) & f == 1L] <- "matched_negative"
    }
  }
  factor(status, levels = c("unselected_negative", "unmatched_positive", "matched_negative", "matched_positive"))
}

.default_matching_pca_features <- function(gr, features = NULL) {
  mc <- S4Vectors::mcols(gr)
  if (!is.null(features)) return(intersect(as.character(features), colnames(mc)))

  md <- S4Vectors::metadata(gr)
  diag <- md$match_diagnostics
  if (is.null(diag) && !is.null(md$posMatchR)) diag <- md$posMatchR$match_diagnostics

  meta_col <- NULL
  if (!is.null(diag$meta_col) && diag$meta_col %in% colnames(mc)) meta_col <- diag$meta_col
  if (is.null(meta_col)) {
    meta_col <- intersect(c("metagene_split3", "metagene_prop", "metagene_split"), colnames(mc))[1]
  }

  loc_cols <- intersect(c("loc_fiveUTR", "loc_coding", "loc_threeUTR"), colnames(mc))
  bin_cols <- character(0)
  if (!is.null(diag$bin_cols_used)) bin_cols <- intersect(as.character(diag$bin_cols_used), colnames(mc))
  if (!length(bin_cols)) {
    bin_cols <- intersect(c(
      "nearest_exon_junction_dist", "start_dist_tx", "stop_dist_tx",
      "start_dist", "stop_dist", "tx_len", "feature_width", "segment_rank",
      "dist_from_feature_start", "dist_from_feature_end"
    ), colnames(mc))
  }

  unique(c(meta_col, loc_cols, bin_cols))
}

.make_matching_pca_matrix <- function(gr, features) {
  mc <- S4Vectors::mcols(gr)
  mats <- list()
  names_out <- character(0)

  add_vec <- function(name, value) {
    value <- suppressWarnings(as.numeric(value))
    mats[[length(mats) + 1L]] <<- value
    names_out[[length(names_out) + 1L]] <<- name
  }

  for (cl in features) {
    if (!(cl %in% colnames(mc))) next
    v <- mc[[cl]]
    if (is.numeric(v) || is.integer(v) || is.logical(v)) {
      add_vec(cl, v)
    } else {
      vv <- as.character(v)
      lv <- sort(unique(vv[!is.na(vv) & nzchar(vv)]))
      if (length(lv) > 1L && length(lv) <= 25L) {
        for (x in lv) add_vec(paste0(cl, "=", x), as.integer(vv == x))
      }
    }
  }

  if (!length(mats)) stop("No usable numeric PCA features were found.")
  X <- as.data.frame(mats, stringsAsFactors = FALSE)
  colnames(X) <- names_out
  X
}

#' Build PCA coordinates for matching diagnostics
#'
#' Constructs a PCA representation of the numeric covariates used by
#' \code{match_background()}, then labels rows as matched positives, matched
#' negatives, unmatched positives, or unselected negatives. For large candidate
#' universes, unselected negatives are sampled before PCA to keep the diagnostic
#' plot manageable.
#'
#' @param gr A \code{GRanges} returned by \code{match_background()} or
#'   \code{match_random_background()}.
#' @param features Optional metadata columns to use. If NULL, the function uses
#'   \code{metadata(gr)$match_diagnostics$meta_col}, location indicator columns,
#'   and \code{metadata(gr)$match_diagnostics$bin_cols_used} when available.
#' @param label_col Binary label column.
#' @param set_col Match-set column.
#' @param positive_value Value marking matched positives in \code{set_col}.
#' @param negative_value Value marking matched negatives in \code{set_col}.
#' @param max_unselected Maximum number of unselected negatives to include.
#' @param max_unmatched_positive Maximum number of unmatched positives to include.
#' @param seed Random seed for diagnostic sampling.
#' @param center,scale Passed to \code{stats::prcomp()}.
#'
#' @return A data.frame with PC coordinates, plotting status, row identifiers,
#'   and attributes containing the \code{prcomp} object and features used.
#' @export
matching_pca_data <- function(gr,
                              features = NULL,
                              label_col = "label",
                              set_col = "match_set",
                              positive_value = "positive",
                              negative_value = "matched_negative",
                              max_unselected = 10000L,
                              max_unmatched_positive = 5000L,
                              seed = 1L,
                              center = TRUE,
                              scale = TRUE) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  gr <- prepare_sites(gr, label = NULL, label_col = label_col, strip_mcols = FALSE)
  mc <- S4Vectors::mcols(gr)

  features <- .default_matching_pca_features(gr, features = features)
  if (!length(features)) stop("No PCA features available.")

  status <- .matching_status(
    gr,
    label_col = label_col,
    set_col = set_col,
    positive_value = positive_value,
    negative_value = negative_value
  )

  keep <- rep(FALSE, length(gr))
  keep[status %in% c("matched_positive", "matched_negative")] <- TRUE

  set.seed(seed)
  idx_unselected <- which(status == "unselected_negative")
  if (length(idx_unselected)) {
    n <- min(length(idx_unselected), as.integer(max_unselected))
    keep[.posmatchr_sample(idx_unselected, size = n)] <- TRUE
  }

  idx_unmatched_pos <- which(status == "unmatched_positive")
  if (length(idx_unmatched_pos)) {
    n <- min(length(idx_unmatched_pos), as.integer(max_unmatched_positive))
    keep[.posmatchr_sample(idx_unmatched_pos, size = n)] <- TRUE
  }

  idx <- which(keep)
  if (length(idx) < 3L) stop("Too few rows available for PCA diagnostic plot.")

  X <- .make_matching_pca_matrix(gr[idx], features = features)
  ok <- stats::complete.cases(X)
  if (!any(ok)) stop("No complete rows for PCA after selecting features.")
  X <- X[ok, , drop = FALSE]
  idx <- idx[ok]

  sds <- vapply(X, stats::sd, numeric(1), na.rm = TRUE)
  keep_col <- is.finite(sds) & sds > 0
  X <- X[, keep_col, drop = FALSE]
  if (ncol(X) < 2L) stop("Fewer than two non-constant PCA features remain.")

  pca <- stats::prcomp(X, center = center, scale. = scale)
  pct <- (pca$sdev^2) / sum(pca$sdev^2)

  out <- data.frame(
    row_index = idx,
    site_id = names(gr)[idx],
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    status = status[idx],
    label = .normalise_binary_label(mc[[label_col]], strict = TRUE)[idx],
    match_set = if (!is.null(set_col) && set_col %in% colnames(mc)) as.character(mc[[set_col]][idx]) else NA_character_,
    stringsAsFactors = FALSE
  )
  out$status <- factor(out$status, levels = levels(status))

  attr(out, "pca") <- pca
  attr(out, "features_used") <- colnames(X)
  attr(out, "variance_explained") <- pct
  out
}


.matching_umap_data <- function(gr,
                                features = NULL,
                                label_col = "label",
                                set_col = "match_set",
                                positive_value = "positive",
                                negative_value = "matched_negative",
                                max_unselected = 10000L,
                                max_unmatched_positive = 5000L,
                                seed = 1L,
                                center = TRUE,
                                scale = TRUE,
                                n_neighbors = 15L,
                                min_dist = 0.1,
                                metric = "euclidean",
                                n_threads = 1L) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("UMAP matching diagnostics require the optional package 'uwot'. Install it with install.packages('uwot').")
  }
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  gr <- prepare_sites(gr, label = NULL, label_col = label_col, strip_mcols = FALSE)
  mc <- S4Vectors::mcols(gr)

  features <- .default_matching_pca_features(gr, features = features)
  if (!length(features)) stop("No matching diagnostic features available.")

  status <- .matching_status(
    gr,
    label_col = label_col,
    set_col = set_col,
    positive_value = positive_value,
    negative_value = negative_value
  )

  keep <- rep(FALSE, length(gr))
  keep[status %in% c("matched_positive", "matched_negative")] <- TRUE

  set.seed(seed)
  idx_unselected <- which(status == "unselected_negative")
  if (length(idx_unselected)) {
    n <- min(length(idx_unselected), as.integer(max_unselected))
    keep[.posmatchr_sample(idx_unselected, size = n)] <- TRUE
  }

  idx_unmatched_pos <- which(status == "unmatched_positive")
  if (length(idx_unmatched_pos)) {
    n <- min(length(idx_unmatched_pos), as.integer(max_unmatched_positive))
    keep[.posmatchr_sample(idx_unmatched_pos, size = n)] <- TRUE
  }

  idx <- which(keep)
  if (length(idx) < 4L) stop("Too few rows available for UMAP diagnostic plot.")

  X <- .make_matching_pca_matrix(gr[idx], features = features)
  ok <- stats::complete.cases(X)
  if (!any(ok)) stop("No complete rows for UMAP after selecting features.")
  X <- X[ok, , drop = FALSE]
  idx <- idx[ok]

  sds <- vapply(X, stats::sd, numeric(1), na.rm = TRUE)
  keep_col <- is.finite(sds) & sds > 0
  X <- X[, keep_col, drop = FALSE]
  if (ncol(X) < 2L) stop("Fewer than two non-constant UMAP features remain.")
  if (nrow(X) < 4L) stop("Too few complete rows remain for UMAP.")

  X_mat <- as.matrix(X)
  if (isTRUE(center) || isTRUE(scale)) {
    X_mat <- base::scale(X_mat, center = center, scale = scale)
  }
  X_mat <- as.matrix(X_mat)
  X_mat[!is.finite(X_mat)] <- 0

  n_neighbors <- as.integer(n_neighbors[1L])
  if (!is.finite(n_neighbors) || n_neighbors < 2L) n_neighbors <- 2L
  n_neighbors <- min(n_neighbors, nrow(X_mat) - 1L)

  set.seed(seed)
  emb <- uwot::umap(
    X_mat,
    n_components = 2L,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    n_threads = n_threads,
    n_sgd_threads = 1L,
    verbose = FALSE
  )

  out <- data.frame(
    row_index = idx,
    site_id = names(gr)[idx],
    UMAP1 = emb[, 1L],
    UMAP2 = emb[, 2L],
    status = status[idx],
    label = .normalise_binary_label(mc[[label_col]], strict = TRUE)[idx],
    match_set = if (!is.null(set_col) && set_col %in% colnames(mc)) as.character(mc[[set_col]][idx]) else NA_character_,
    stringsAsFactors = FALSE
  )
  out$status <- factor(out$status, levels = levels(status))

  attr(out, "features_used") <- colnames(X)
  attr(out, "umap_parameters") <- list(
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    center = center,
    scale = scale,
    seed = seed
  )
  out
}

#' Plot a PCA or UMAP diagnostic of matching covariates
#'
#' Shows sites in a two-dimensional representation of the numeric covariates
#' used for matching. By default this is PCA. Set \code{reduction = "umap"}
#' to use an optional UMAP representation. Matched positives and matched
#' negatives are highlighted, while unselected negatives are shown as a light
#' grey reference cloud.
#'
#' @inheritParams matching_pca_data
#' @param point_size Point size.
#' @param alpha_unselected Alpha for unselected negatives.
#' @param alpha_highlight Alpha for matched positives/negatives.
#' @param reduction Two-dimensional reduction to plot: \code{"pca"} or \code{"umap"}.
#' @param colour_by Colour points by \code{"status"} or by a metadata column, e.g. \code{"location"}.
#' @param umap_n_neighbors,umap_min_dist,umap_metric,umap_n_threads UMAP
#'   settings used when \code{reduction = "umap"}. The optional package
#'   \pkg{uwot} must be installed.
#'
#' @return A ggplot object.
#' @export
plot_matching_pca <- function(gr,
                              features = NULL,
                              label_col = "label",
                              set_col = "match_set",
                              positive_value = "positive",
                              negative_value = "matched_negative",
                              max_unselected = 10000L,
                              max_unmatched_positive = 5000L,
                              seed = 1L,
                              center = TRUE,
                              scale = TRUE,
                              point_size = 0.8,
                              alpha_unselected = 0.18,
                              alpha_highlight = 0.75,
                              reduction = c("pca", "umap"),
                              colour_by = "status",
                              umap_n_neighbors = 15L,
                              umap_min_dist = 0.1,
                              umap_metric = "euclidean",
                              umap_n_threads = 1L) {
  reduction <- match.arg(reduction)

  if (reduction == "umap") {
    df <- .matching_umap_data(
      gr = gr,
      features = features,
      label_col = label_col,
      set_col = set_col,
      positive_value = positive_value,
      negative_value = negative_value,
      max_unselected = max_unselected,
      max_unmatched_positive = max_unmatched_positive,
      seed = seed,
      center = center,
      scale = scale,
      n_neighbors = umap_n_neighbors,
      min_dist = umap_min_dist,
      metric = umap_metric,
      n_threads = umap_n_threads
    )
    xcol <- "UMAP1"
    ycol <- "UMAP2"
    xl <- "UMAP1"
    yl <- "UMAP2"
  } else {
  df <- matching_pca_data(
    gr = gr,
    features = features,
    label_col = label_col,
    set_col = set_col,
    positive_value = positive_value,
    negative_value = negative_value,
    max_unselected = max_unselected,
    max_unmatched_positive = max_unmatched_positive,
    seed = seed,
    center = center,
    scale = scale
  )
  pct <- attr(df, "variance_explained")
  xl <- if (length(pct) >= 1L) paste0("PC1 (", round(100 * pct[1L], 1), "%)") else "PC1"
  yl <- if (length(pct) >= 2L) paste0("PC2 (", round(100 * pct[2L], 1), "%)") else "PC2"
  xcol <- "PC1"
  ycol <- "PC2"
  }

  df$x_plot <- df[[xcol]]
  df$y_plot <- df[[ycol]]

  if (!is.null(colour_by) && !identical(colour_by, "status")) {
    mc_all <- S4Vectors::mcols(gr)
    if (!(colour_by %in% colnames(mc_all))) {
      stop("colour_by='", colour_by, "' was not found in mcols(gr).")
    }
    df$colour_value <- as.character(mc_all[[colour_by]][df$row_index])
    df$colour_value[is.na(df$colour_value) | !nzchar(df$colour_value)] <- "missing"
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = x_plot, y = y_plot, colour = colour_value)) +
        ggplot2::geom_point(alpha = alpha_highlight, size = point_size) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = xl, y = yl, colour = colour_by)
    )
  }

  bg <- df[df$status == "unselected_negative", , drop = FALSE]
  um <- df[df$status == "unmatched_positive", , drop = FALSE]
  hi <- df[df$status %in% c("matched_negative", "matched_positive"), , drop = FALSE]

  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = xl, y = yl, colour = NULL)

  if (nrow(bg)) {
    p <- p + ggplot2::geom_point(
      data = bg,
      ggplot2::aes(x = x_plot, y = y_plot),
      colour = "grey80",
      alpha = alpha_unselected,
      size = point_size
    )
  }

  if (nrow(um)) {
    p <- p + ggplot2::geom_point(
      data = um,
      ggplot2::aes(x = x_plot, y = y_plot),
      colour = "grey55",
      alpha = min(alpha_highlight, 0.35),
      size = point_size
    )
  }

  if (nrow(hi)) {
    p <- p + ggplot2::geom_point(
      data = hi,
      ggplot2::aes(x = x_plot, y = y_plot, colour = status),
      alpha = alpha_highlight,
      size = point_size
    )
  }

  p
}

#' Summarise k-mer balance in matched sets
#'
#' @param gr A matched \code{GRanges}; either a full object returned by a matcher
#'   or the paired subset returned by \code{subset_matched_sets()}.
#' @param kmer_col K-mer metadata column.
#' @param set_col Match-set column.
#' @param positive_value Value marking matched positives.
#' @param negative_value Value marking matched negatives.
#' @param subset_first If TRUE, call \code{subset_matched_sets()} before counting
#'   when canonical match-ID columns are present.
#'
#' @return A data.frame with positive and matched-negative counts by k-mer.
#' @export
summarise_matched_kmer_balance <- function(gr,
                                           kmer_col = "kmer",
                                           set_col = "match_set",
                                           positive_value = "positive",
                                           negative_value = "matched_negative",
                                           subset_first = TRUE) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  if (isTRUE(subset_first) && all(c("matched_negative_id", "matched_positive_id") %in% colnames(S4Vectors::mcols(gr)))) {
    gr <- subset_matched_sets(gr)
  }
  mc <- S4Vectors::mcols(gr)
  if (!(kmer_col %in% colnames(mc))) stop("Missing kmer_col: ", kmer_col)
  if (!(set_col %in% colnames(mc))) stop("Missing set_col: ", set_col)

  km <- as.character(mc[[kmer_col]])
  ss <- as.character(mc[[set_col]])
  km[is.na(km) | !nzchar(km)] <- NA_character_
  vals <- sort(unique(km[!is.na(km)]))
  pos_mask <- !is.na(ss) & ss == positive_value & !is.na(km)
  neg_mask <- !is.na(ss) & ss == negative_value & !is.na(km)
  pos <- tabulate(match(km[pos_mask], vals), nbins = length(vals))
  neg <- tabulate(match(km[neg_mask], vals), nbins = length(vals))
  out <- data.frame(
    kmer = vals,
    positive = as.integer(pos),
    matched_negative = as.integer(neg),
    difference = as.integer(pos - neg),
    stringsAsFactors = FALSE
  )
  out[order(abs(out$difference), decreasing = TRUE), , drop = FALSE]
}
