#' Compute occurrence profiles around single-nucleotide sites
#'
#' This is the main enrichment-profile interface in posMatchR. The query can be
#' either a supplied point-feature GRanges object, such as iCLIP crosslink sites,
#' or an RNA/DNA motif such as "UNUNU", "YUGUM", or "DRACH". Motif queries are
#' scanned internally in strand-oriented windows around the focal sites and are
#' treated as motif-hit occurrences. Feature queries are counted by overlap with
#' strand-oriented windows around the focal sites.
#'
#' Character motifs are interpreted as exact/IUPAC RNA or DNA patterns. RNA U is
#' converted internally to DNA T before scanning the genome sequence. Relative
#' positions are reported in RNA orientation: for minus-strand focal sites,
#' upstream/downstream positions are flipped so that +10 means 10 nt downstream in
#' the transcript direction.
#'
#' @param gr A single-nucleotide focal GRanges. For matched-set comparisons, pass
#'   subset_matched_sets(gr) or set matched_only = TRUE.
#' @param query Either a point-feature GRanges/list of GRanges objects, or a motif
#'   query. Character motifs are treated as exact/IUPAC RNA/DNA motifs. A
#'   universalmotif object can be supplied when method = "universalmotif".
#' @param genome A BSgenome-like object. Required for motif queries; ignored for
#'   supplied GRanges feature queries.
#' @param query_type "auto", "features", or "motif". In auto mode, GRanges inputs
#'   are treated as features and character/universalmotif inputs are treated as
#'   motifs.
#' @param query_name Optional display name. For multiple motifs or feature sets,
#'   this can be a vector of names.
#' @param window Number of bp/nt on each side of the focal site.
#' @param set_col Optional grouping column. Defaults to "match_set" if present.
#' @param positive_value Value in set_col marking matched positives.
#' @param negative_value Value in set_col marking matched negatives.
#' @param include_other Include other groups from set_col.
#' @param matched_only If TRUE and canonical match columns are present, subset to
#'   reciprocal matched pairs before profiling.
#' @param method Motif-scanning method: "auto", "iupac", or "universalmotif".
#'   Ignored for supplied GRanges feature queries.
#' @param scan_rc Scan reverse-complement motif instances as well. Usually FALSE
#'   for RNA-oriented analyses when site strands are known.
#' @param hit_position For motifs, whether to record the start, end, or centre of
#'   each motif hit.
#' @param feature_position For supplied features, whether to record the start,
#'   end, or centre of each feature.
#' @param require_same_strand If TRUE, supplied feature queries are counted only
#'   on the same strand as focal sites.
#' @param weight_col Optional metadata column in supplied feature GRanges used for
#'   score-weighted profiles.
#' @param threshold Threshold passed to universalmotif::scan_sequences().
#' @param threshold_type Threshold type passed to universalmotif::scan_sequences().
#' @param nthreads Number of threads for universalmotif scanning.
#' @param bin_size Position bin size.
#' @param smooth_window Moving-average smoothing width in bins. Defaults to 7.
#' @param seqstyle Optional seqlevel style to apply.
#' @param chrs Optional seqlevels to keep. Usually not needed if inputs are already harmonised.
#' @param drop_unstranded If TRUE, omit focal sites on '*' strand.
#' @param drop_edge_positions Logical or integer. If TRUE, omit edge positions.
#' @param edge_trim Number of bp/nt to omit from each edge. Overrides drop_edge_positions.
#' @param center_exclude Number of bp/nt around the central focal site to omit.
#'   Defaults to 0.
#' @param normalise_per Scale factor for normalised rates. For example, 1000
#'   reports occurrences per 1000 focal sites.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'     "chr1", IRanges::IRanges(c(20, 40), width = 1), strand = c("+", "+")
#' )
#' features <- GenomicRanges::GRanges(
#'     "chr1", IRanges::IRanges(c(22, 43), width = 1), strand = c("+", "+")
#' )
#' site_enrichment_profile(gr, query = features, window = 5, smooth_window = 1)
#'
#' @return A data.frame with relative position, group, feature/query name, counts,
#'   site counts, and normalised rates.
#' @export
site_enrichment_profile <- function(gr,
                                    query,
                                    genome = NULL,
                                    query_type = c("auto", "features", "motif"),
                                    query_name = NULL,
                                    window = 250L,
                                    set_col = NULL,
                                    positive_value = "positive",
                                    negative_value = "matched_negative",
                                    include_other = FALSE,
                                    matched_only = FALSE,
                                    method = c("auto", "iupac", "universalmotif"),
                                    scan_rc = FALSE,
                                    hit_position = c("center", "start", "end"),
                                    feature_position = c("center", "start", "end"),
                                    require_same_strand = TRUE,
                                    weight_col = NULL,
                                    threshold = 0.8,
                                    threshold_type = "logodds",
                                    nthreads = 1L,
                                    bin_size = 1L,
                                    smooth_window = 7L,
                                    seqstyle = NULL,
                                    chrs = NULL,
                                    drop_unstranded = TRUE,
                                    drop_edge_positions = FALSE,
                                    edge_trim = NULL,
                                    center_exclude = 0L,
                                    normalise_per = 1000L) {
  query_type <- match.arg(query_type)
  method <- match.arg(method)
  hit_position <- match.arg(hit_position)
  feature_position <- match.arg(feature_position)
  if (query_type == "auto") {
    if (inherits(query, "GRanges") || (is.list(query) && all(vapply(query, inherits, logical(1), what = "GRanges")))) {
      query_type <- "features"
    } else {
      query_type <- "motif"
    }
  }

  if (query_type == "features") {
    prof <- feature_occurrence_profile(
      gr = gr,
      features = query,
      window = window,
      set_col = set_col,
      positive_value = positive_value,
      negative_value = negative_value,
      include_other = include_other,
      matched_only = matched_only,
      feature_name = query_name,
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
    prof$query_type <- "features"
    return(prof)
  }

  if (is.null(genome)) {
    stop("`genome` is required when `query` is a motif.")
  }

  prof <- motif_enrichment_profile(
    gr = gr,
    genome = genome,
    motif = query,
    window = window,
    set_col = set_col,
    positive_value = positive_value,
    negative_value = negative_value,
    include_other = include_other,
    matched_only = matched_only,
    method = method,
    motif_name = query_name,
    scan_rc = scan_rc,
    hit_position = hit_position,
    threshold = threshold,
    threshold_type = threshold_type,
    nthreads = nthreads,
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
  names(prof)[names(prof) == "motif"] <- "feature"
  prof$query_type <- "motif"
  prof
}

#' Plot occurrence enrichment around single-nucleotide sites
#'
#' A single plotting interface for either motif enrichment or supplied point
#' features such as iCLIP sites. Pass a motif string, e.g. "UNUNU", or pass a
#' single-width GRanges object of feature sites.
#'
#' @inheritParams site_enrichment_profile
#' @param y Value to plot. The default plots smoothed occurrences per
#'   normalise_per focal sites.
#' @param x_label Optional x-axis label.
#' @param y_label Optional y-axis label.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'     "chr1", IRanges::IRanges(c(20, 40), width = 1), strand = c("+", "+")
#' )
#' features <- GenomicRanges::GRanges(
#'     "chr1", IRanges::IRanges(c(22, 43), width = 1), strand = c("+", "+")
#' )
#' plot_site_enrichment(gr, query = features, window = 5, smooth_window = 1)
#'
#' @return A ggplot object.
#' @export
plot_site_enrichment <- function(gr,
                                 query,
                                 genome = NULL,
                                 query_type = c("auto", "features", "motif"),
                                 query_name = NULL,
                                 window = 250L,
                                 set_col = NULL,
                                 positive_value = "positive",
                                 negative_value = "matched_negative",
                                 include_other = FALSE,
                                 matched_only = FALSE,
                                 method = c("auto", "iupac", "universalmotif"),
                                 scan_rc = FALSE,
                                 hit_position = c("center", "start", "end"),
                                 feature_position = c("center", "start", "end"),
                                 require_same_strand = TRUE,
                                 weight_col = NULL,
                                 threshold = 0.8,
                                 threshold_type = "logodds",
                                 nthreads = 1L,
                                 bin_size = 1L,
                                 smooth_window = 7L,
                                 seqstyle = NULL,
                                 chrs = NULL,
                                 drop_unstranded = TRUE,
                                 drop_edge_positions = TRUE,
                                 edge_trim = NULL,
                                 center_exclude = 0L,
                                 normalise_per = 1000L,
                                 y = c("hits_per_n_sites_smoothed", "hits_per_n_sites", "hits_per_site_smoothed", "hits_per_site", "hits_per_1000_sites_smoothed", "hits_per_1000_sites", "signal_per_n_sites_smoothed", "signal_per_n_sites", "signal_per_site_smoothed", "signal_per_site", "count", "signal"),
                                 x_label = NULL,
                                 y_label = NULL) {
  query_type <- match.arg(query_type)
  method <- match.arg(method)
  hit_position <- match.arg(hit_position)
  feature_position <- match.arg(feature_position)
  y <- match.arg(y)

  prof <- site_enrichment_profile(
    gr = gr,
    query = query,
    genome = genome,
    query_type = query_type,
    query_name = query_name,
    window = window,
    set_col = set_col,
    positive_value = positive_value,
    negative_value = negative_value,
    include_other = include_other,
    matched_only = matched_only,
    method = method,
    scan_rc = scan_rc,
    hit_position = hit_position,
    feature_position = feature_position,
    require_same_strand = require_same_strand,
    weight_col = weight_col,
    threshold = threshold,
    threshold_type = threshold_type,
    nthreads = nthreads,
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

  if (!(y %in% names(prof))) {
    stop("Requested y column is not available for this query: ", y)
  }

  prof$y_plot <- prof[[y]]
  if (center_exclude > 0L) {
    prof$plot_segment <- ifelse(prof$relative_position < 0, "left", "right")
  } else {
    prof$plot_segment <- "all"
  }

  query_label <- ifelse(unique(prof$query_type)[1L] == "motif", "Motif", "Feature")
  ylab <- y_label %||% if (y == "count") {
    paste0(query_label, " occurrence count")
  } else if (y == "signal") {
    "Feature signal"
  } else if (grepl("signal", y)) {
    paste0("Feature signal per ", format(normalise_per, scientific = FALSE, trim = TRUE), " sites")
  } else if (grepl("hits_per_n_sites", y)) {
    paste0(query_label, " occurrences per ", format(normalise_per, scientific = FALSE, trim = TRUE), " sites")
  } else if (grepl("1000", y)) {
    paste0(query_label, " occurrences per 1,000 sites")
  } else {
    paste0(query_label, " occurrences per site")
  }

  xlab <- x_label %||% "Distance from focal site (nt; RNA-oriented genomic)"

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
    ggplot2::labs(x = xlab, y = ylab, colour = NULL, linetype = query_label) +
    ggplot2::theme_bw()

  if (length(unique(prof$feature)) == 1L) {
    p <- p + ggplot2::guides(linetype = "none")
  }

  p
}
