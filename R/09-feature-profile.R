.posmatchr_feature_sets <- function(features, feature_name = NULL) {
  if (inherits(features, "GRanges")) {
    out <- list(features)
  } else if (is.list(features) && all(vapply(features, inherits, logical(1), what = "GRanges"))) {
    out <- features
  } else {
    stop("`features` must be a GRanges or a named list of GRanges objects.")
  }

  nms <- names(out)
  if (!is.null(feature_name)) {
    feature_name <- as.character(feature_name)
    if (length(feature_name) == length(out)) {
      nms <- feature_name
    } else if (length(feature_name) == 1L && length(out) == 1L) {
      nms <- feature_name
    }
  }
  if (is.null(nms) || length(nms) != length(out) || any(!nzchar(nms))) {
    nms <- paste0("feature_", seq_along(out))
  }
  names(out) <- nms
  out
}

.posmatchr_prepare_feature_profile_ranges <- function(gr,
                                                       features,
                                                       seqstyle = NULL,
                                                       chrs = NULL,
                                                       drop_unstranded = TRUE,
                                                       require_same_strand = TRUE) {
  gr <- prepare_sites(gr, label = NULL, strip_mcols = FALSE)
  if (any(GenomicRanges::width(features) != 1L)) {
    features <- IRanges::resize(features, width = 1L, fix = "center")
  }

  if (isTRUE(drop_unstranded)) {
    st <- as.character(GenomicRanges::strand(gr))
    keep <- st %in% c("+", "-")
    if (!any(keep)) stop("No '+' or '-' stranded focal sites remain.")
    gr <- gr[keep]
  }

  if (isTRUE(require_same_strand)) {
    stf <- as.character(GenomicRanges::strand(features))
    keepf <- stf %in% c("+", "-")
    if (!any(keepf)) stop("No '+' or '-' stranded feature sites remain.")
    features <- features[keepf]
  }

  if (!is.null(seqstyle)) {
    suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- seqstyle, silent = TRUE))
    suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(features) <- seqstyle, silent = TRUE))
  }

  common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(features))
  if (!length(common) && is.null(seqstyle)) {
    st <- tryCatch(GenomeInfoDb::seqlevelsStyle(features), error = function(e) character(0))
    if (length(st)) suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- st[1], silent = TRUE))
    common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(features))
  }
  if (!length(common)) {
    gr <- .rename_seqlevels_to_target(gr, GenomeInfoDb::seqlevels(features))
    common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(features))
  }

  if (!is.null(chrs)) {
    target_common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(features))
    common <- .map_to_target_seqlevels(chrs, target_common)$mapped
  }
  common <- unique(common[!is.na(common) & nzchar(common)])
  if (!length(common)) {
    stop(
      "No common seqlevels between focal `gr` and `features`.\n",
      "seqlevels(gr): ", paste(GenomeInfoDb::seqlevels(gr), collapse = ", "), "\n",
      "seqlevels(features): ", paste(GenomeInfoDb::seqlevels(features), collapse = ", ")
    )
  }

  gr <- GenomeInfoDb::keepSeqlevels(gr, common, pruning.mode = "coarse")
  features <- GenomeInfoDb::keepSeqlevels(features, common, pruning.mode = "coarse")

  list(gr = gr, features = features)
}

.posmatchr_feature_positions <- function(features, feature_position) {
  switch(
    feature_position,
    start = as.numeric(GenomicRanges::start(features)),
    end = as.numeric(GenomicRanges::end(features)),
    center = (as.numeric(GenomicRanges::start(features)) + as.numeric(GenomicRanges::end(features))) / 2,
    stop("Unsupported feature_position.")
  )
}

.posmatchr_smooth_profile <- function(prof, value_col, group_cols, smooth_window) {
  smooth_window <- as.integer(smooth_window)
  out <- prof[[value_col]]
  if (!is.finite(smooth_window) || smooth_window <= 1L || !nrow(prof)) return(out)

  smooth_one <- function(x) {
    as.numeric(stats::filter(x, rep(1 / smooth_window, smooth_window), sides = 2))
  }
  key <- do.call(paste, c(prof[group_cols], sep = "||"))
  for (k in unique(key)) {
    idx <- which(key == k)
    vals <- smooth_one(prof[[value_col]][idx])
    vals[is.na(vals)] <- prof[[value_col]][idx][is.na(vals)]
    out[idx] <- vals
  }
  out
}

#' Compute a point-feature occurrence profile around single-nucleotide sites
#'
#' Counts another set of point features, such as iCLIP crosslink sites, in
#' strand-oriented windows centred on focal single-nucleotide sites. This is the
#' point-feature analogue of \code{motif_enrichment_profile()}: instead of
#' scanning sequence motifs, it overlaps focal-site windows with a supplied
#' feature \code{GRanges} and reports counts or weighted signal by relative
#' position.
#'
#' Relative positions are reported in RNA orientation. For a plus-strand focal
#' site, a feature at genomic position \code{site + 10} has relative position
#' +10. For a minus-strand focal site, a feature at genomic position
#' \code{site - 10} has relative position +10. By default, only features on the
#' same strand as the focal site are counted.
#'
#' @param gr A single-nucleotide focal \code{GRanges}; for matched-set
#'   comparisons, pass \code{subset_matched_sets(gr)}.
#' @param features A point-feature \code{GRanges}, or a named list of such
#'   objects. For example, iCLIP crosslink sites.
#' @param window Number of bp on each side of the focal site.
#' @param set_col Optional grouping column. Defaults to \code{"match_set"} if present.
#' @param positive_value Value in \code{set_col} marking matched positives.
#' @param negative_value Value in \code{set_col} marking matched negatives.
#' @param include_other Include other groups from \code{set_col}.
#' @param matched_only If TRUE and canonical match columns are present, subset to
#'   reciprocal matched pairs before profiling.
#' @param feature_name Optional display name(s) for \code{features}.
#' @param feature_position Whether each feature is represented by its start, end,
#'   or centre coordinate. Width-1 features are unaffected.
#' @param require_same_strand If TRUE, count only same-strand feature/focal-site
#'   overlaps.
#' @param weight_col Optional numeric metadata column in \code{features}; if
#'   supplied, the profile additionally sums this signal rather than only counts.
#' @param bin_size Position bin size in bp.
#' @param smooth_window Optional moving-average smoothing width in bins. Use 1 for no smoothing.
#' @param seqstyle Optional seqlevel style to apply to both focal and feature ranges.
#' @param chrs Optional seqlevels to keep. Usually not needed if the two GRanges
#'   already have matching seqlevels.
#' @param drop_unstranded If TRUE, omit focal sites on '*' strand.
#' @param drop_edge_positions Logical or integer. If TRUE, omit edge positions
#'   from the profile; if FALSE, keep them. An integer is interpreted as the
#'   number of bp to trim from each edge.
#' @param edge_trim Integer number of bp to omit from each edge of the profile.
#' @param center_exclude Integer number of bp around the central site to omit
#'   from the profile. Defaults to 0.
#' @param normalise_per Scale factor for normalised rates. For example, 1000
#'   reports feature occurrences per 1000 focal sites.
#'
#' @return A data.frame with relative position, group, feature, counts, optional
#'   signal, site counts, and normalised rates.
#' @noRd
# Internal implementation used by site_enrichment_profile().
feature_occurrence_profile <- function(gr,
                                       features,
                                       window = 250L,
                                       set_col = NULL,
                                       positive_value = "positive",
                                       negative_value = "matched_negative",
                                       include_other = FALSE,
                                       matched_only = FALSE,
                                       feature_name = NULL,
                                       feature_position = c("center", "start", "end"),
                                       require_same_strand = TRUE,
                                       weight_col = NULL,
                                       bin_size = 1L,
                                       smooth_window = 1L,
                                       seqstyle = NULL,
                                       chrs = NULL,
                                       drop_unstranded = TRUE,
                                       drop_edge_positions = FALSE,
                                       edge_trim = NULL,
                                       center_exclude = 0L,
                                       normalise_per = 1000L) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  feature_position <- match.arg(feature_position)
  window <- as.integer(window)
  bin_size <- as.integer(bin_size)
  smooth_window <- as.integer(smooth_window)
  if (!is.finite(window) || window < 1L) stop("`window` must be a positive integer.")
  if (!is.finite(bin_size) || bin_size < 1L) stop("`bin_size` must be a positive integer.")
  edge_trim <- .posmatchr_edge_trim(drop_edge_positions = drop_edge_positions, edge_trim = edge_trim)
  if (edge_trim >= window) stop("`edge_trim` must be smaller than `window`.")
  center_exclude <- as.integer(center_exclude[1L])
  if (!is.finite(center_exclude) || center_exclude < 0L) stop("`center_exclude` must be a non-negative integer.")
  if (center_exclude >= window) stop("`center_exclude` must be smaller than `window`.")
  normalise_per <- as.numeric(normalise_per[1L])
  if (!is.finite(normalise_per) || normalise_per <= 0) stop("`normalise_per` must be a positive number.")

  if (matched_only && all(c("matched_negative_id", "matched_positive_id") %in% colnames(S4Vectors::mcols(gr)))) {
    gr <- subset_matched_pairs(gr, return_diagnostics = FALSE)
  }

  group <- .posmatchr_motif_groups(
    gr,
    set_col = set_col,
    positive_value = positive_value,
    negative_value = negative_value,
    include_other = include_other
  )
  keep_group <- !is.na(group) & nzchar(group)
  if (!any(keep_group)) stop("No rows remain after applying set_col/group filters.")
  gr <- gr[keep_group]
  group <- group[keep_group]

  tmp_row_col <- ".posmatchr_feature_row"
  S4Vectors::mcols(gr)[[tmp_row_col]] <- seq_along(gr)
  feature_sets <- .posmatchr_feature_sets(features, feature_name = feature_name)

  positions <- seq(-window, window, by = bin_size)
  if (!length(positions) || utils::tail(positions, 1L) != window) positions <- unique(c(positions, window))
  if (edge_trim > 0L) {
    positions <- positions[positions >= (-window + edge_trim) & positions <= (window - edge_trim)]
  }
  if (center_exclude > 0L) {
    positions <- positions[abs(positions) > center_exclude]
  }

  out_all <- list()
  out_i <- 0L

  for (feature_label in names(feature_sets)) {
    prep <- .posmatchr_prepare_feature_profile_ranges(
      gr,
      feature_sets[[feature_label]],
      seqstyle = seqstyle,
      chrs = chrs,
      drop_unstranded = drop_unstranded,
      require_same_strand = require_same_strand
    )
    gr_i <- prep$gr
    feat_i <- prep$features

    row_index <- as.integer(S4Vectors::mcols(gr_i)[[tmp_row_col]])
    group_i <- group[row_index]
    group_levels <- unique(group_i)

    site_pos <- as.numeric(GenomicRanges::start(gr_i))
    starts <- pmax(1L, as.integer(site_pos) - window)
    ends <- as.integer(site_pos) + window
    win <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(gr_i),
      ranges = IRanges::IRanges(start = starts, end = ends),
      strand = GenomicRanges::strand(gr_i)
    )
    S4Vectors::mcols(win)[[tmp_row_col]] <- row_index

    hits <- GenomicRanges::findOverlaps(win, feat_i, ignore.strand = !isTRUE(require_same_strand))

    n_sites <- as.data.frame(table(group = group_i), stringsAsFactors = FALSE)
    names(n_sites) <- c("group", "n_sites")
    n_sites$n_sites <- as.integer(n_sites$n_sites)

    grid <- expand.grid(
      relative_position = positions,
      group = group_levels,
      feature = feature_label,
      stringsAsFactors = FALSE
    )

    if (length(hits)) {
      q <- S4Vectors::queryHits(hits)
      s <- S4Vectors::subjectHits(hits)
      fpos <- .posmatchr_feature_positions(feat_i[s], feature_position = feature_position)
      spos <- site_pos[q]
      st <- as.character(GenomicRanges::strand(gr_i))[q]
      rel <- ifelse(st == "-", spos - fpos, fpos - spos)

      keep <- is.finite(rel) & rel >= -window & rel <= window
      if (edge_trim > 0L) keep <- keep & rel >= (-window + edge_trim) & rel <= (window - edge_trim)
      if (center_exclude > 0L) keep <- keep & abs(rel) > center_exclude

      if (any(keep)) {
        q <- q[keep]
        s <- s[keep]
        rel <- .posmatchr_bin_relative_positions(rel[keep], window = window, bin_size = bin_size)
        hit_df <- data.frame(
          relative_position = rel,
          group = group_i[q],
          feature = feature_label,
          count = 1L,
          stringsAsFactors = FALSE
        )
        if (!is.null(weight_col)) {
          mc_feat <- S4Vectors::mcols(feat_i)
          if (!(weight_col %in% colnames(mc_feat))) {
            stop("`weight_col` not found in feature metadata: ", weight_col)
          }
          w <- suppressWarnings(as.numeric(mc_feat[[weight_col]][s]))
          w[!is.finite(w)] <- 0
          hit_df$signal <- w
        } else {
          hit_df$signal <- 1
        }

        agg_count <- stats::aggregate(
          count ~ relative_position + group + feature,
          data = hit_df,
          FUN = sum
        )
        agg_signal <- stats::aggregate(
          signal ~ relative_position + group + feature,
          data = hit_df,
          FUN = sum
        )
        agg <- merge(agg_count, agg_signal, by = c("relative_position", "group", "feature"), all = TRUE, sort = FALSE)
      } else {
        agg <- data.frame(relative_position = numeric(0), group = character(0), feature = character(0), count = integer(0), signal = numeric(0))
      }
    } else {
      agg <- data.frame(relative_position = numeric(0), group = character(0), feature = character(0), count = integer(0), signal = numeric(0))
    }

    prof <- merge(grid, agg, by = c("relative_position", "group", "feature"), all.x = TRUE, sort = FALSE)
    prof$count[is.na(prof$count)] <- 0L
    prof$signal[is.na(prof$signal)] <- 0
    prof <- merge(prof, n_sites, by = "group", all.x = TRUE, sort = FALSE)
    prof$hits_per_site <- prof$count / prof$n_sites
    prof$hits_per_n_sites <- prof$hits_per_site * normalise_per
    prof$hits_per_1000_sites <- prof$hits_per_site * 1000
    prof$signal_per_site <- prof$signal / prof$n_sites
    prof$signal_per_n_sites <- prof$signal_per_site * normalise_per
    prof$signal_per_1000_sites <- prof$signal_per_site * 1000
    prof$normalise_per <- normalise_per
    prof$require_same_strand <- isTRUE(require_same_strand)
    prof$weight_col <- if (is.null(weight_col)) NA_character_ else as.character(weight_col)
    prof <- prof[order(prof$feature, prof$group, prof$relative_position), , drop = FALSE]

    prof$hits_per_site_smoothed <- .posmatchr_smooth_profile(prof, "hits_per_site", c("feature", "group"), smooth_window)
    prof$hits_per_n_sites_smoothed <- prof$hits_per_site_smoothed * normalise_per
    prof$hits_per_1000_sites_smoothed <- prof$hits_per_site_smoothed * 1000
    prof$signal_per_site_smoothed <- .posmatchr_smooth_profile(prof, "signal_per_site", c("feature", "group"), smooth_window)
    prof$signal_per_n_sites_smoothed <- prof$signal_per_site_smoothed * normalise_per
    prof$signal_per_1000_sites_smoothed <- prof$signal_per_site_smoothed * 1000

    out_i <- out_i + 1L
    out_all[[out_i]] <- prof
  }

  do.call(rbind, out_all)
}

#' Plot point-feature occurrences around single-nucleotide sites
#'
#' Plots the output of \code{feature_occurrence_profile()}, typically to compare
#' iCLIP-site density around positive sites and matched background sites.
#'
#' @inheritParams feature_occurrence_profile
#' @param y Value to plot. By default, plots smoothed feature occurrences per
#'   \code{normalise_per} focal sites. Use a \code{signal_*} column when
#'   \code{weight_col} is supplied and you want score-weighted signal.
#' @param x_label Optional x-axis label.
#' @param y_label Optional y-axis label.
#'
#' @return A ggplot object.
#' @noRd
# Internal implementation used by plot_site_enrichment().
plot_feature_occurrence <- function(gr,
                                    features,
                                    window = 250L,
                                    set_col = NULL,
                                    positive_value = "positive",
                                    negative_value = "matched_negative",
                                    include_other = FALSE,
                                    matched_only = FALSE,
                                    feature_name = NULL,
                                    feature_position = c("center", "start", "end"),
                                    require_same_strand = TRUE,
                                    weight_col = NULL,
                                    bin_size = 1L,
                                    smooth_window = 7L,
                                    seqstyle = NULL,
                                    chrs = NULL,
                                    drop_unstranded = TRUE,
                                    drop_edge_positions = TRUE,
                                    edge_trim = NULL,
                                    center_exclude = 0L,
                                    normalise_per = 1000L,
                                    y = c("hits_per_n_sites_smoothed", "hits_per_n_sites", "hits_per_site_smoothed", "hits_per_site", "signal_per_n_sites_smoothed", "signal_per_n_sites", "signal_per_site_smoothed", "signal_per_site", "count", "signal"),
                                    x_label = NULL,
                                    y_label = NULL) {
  y <- match.arg(y)
  feature_position <- match.arg(feature_position)
  prof <- feature_occurrence_profile(
    gr = gr,
    features = features,
    window = window,
    set_col = set_col,
    positive_value = positive_value,
    negative_value = negative_value,
    include_other = include_other,
    matched_only = matched_only,
    feature_name = feature_name,
    feature_position = feature_position,
    require_same_strand = require_same_strand,
    weight_col = weight_col,
    bin_size = bin_size,
    smooth_window = smooth_window,
    seqstyle = seqstyle,
    chrs = chrs,
    drop_unstranded = drop_unstranded,
    drop_edge_positions = drop_edge_positions,
    edge_trim = edge_trim,
    center_exclude = center_exclude,
    normalise_per = normalise_per
  )
  prof$y_plot <- prof[[y]]
  if (center_exclude > 0L) {
    prof$plot_segment <- ifelse(prof$relative_position < 0, "left", "right")
  } else {
    prof$plot_segment <- "all"
  }

  ylab <- y_label %||% if (y %in% c("count", "signal")) {
    if (y == "count") "Feature-site count" else "Feature signal"
  } else if (grepl("signal", y)) {
    paste0("Feature signal per ", format(normalise_per, scientific = FALSE, trim = TRUE), " sites")
  } else if (grepl("hits_per_n_sites", y)) {
    paste0("Feature sites per ", format(normalise_per, scientific = FALSE, trim = TRUE), " sites")
  } else if (grepl("1000", y)) {
    "Feature sites per 1,000 sites"
  } else {
    "Feature sites per site"
  }
  xlab <- x_label %||% "Distance from focal site (bp; RNA-oriented genomic)"

  p <- ggplot2::ggplot(
    prof,
    ggplot2::aes(
      x = relative_position,
      y = y_plot,
      colour = group,
      linetype = feature,
      group = interaction(group, feature, plot_segment)
    )
  ) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::scale_x_continuous(
      limits = c(-window, window),
      breaks = unique(c(-window, -round(window / 2), 0, round(window / 2), window))
    ) +
    ggplot2::labs(x = xlab, y = ylab, colour = NULL, linetype = "Feature") +
    ggplot2::theme_bw()

  if (length(unique(prof$feature)) == 1L) {
    p <- p + ggplot2::guides(linetype = "none")
  }

  p
}
