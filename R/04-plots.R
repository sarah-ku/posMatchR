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
