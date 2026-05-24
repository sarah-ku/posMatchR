.posmatchr_motif_groups <- function(gr,
                                    set_col = NULL,
                                    positive_value = "positive",
                                    negative_value = "matched_negative",
                                    include_other = FALSE) {
  mc <- S4Vectors::mcols(gr)
  if (is.null(set_col) && "match_set" %in% colnames(mc)) set_col <- "match_set"

  if (!is.null(set_col) && set_col %in% colnames(mc)) {
    s <- as.character(mc[[set_col]])
    group <- ifelse(
      s == positive_value,
      "positive",
      ifelse(s == negative_value, "negative", ifelse(include_other, as.character(s), NA_character_))
    )
  } else {
    group <- rep("sites", length(gr))
  }

  group[is.na(group) | !nzchar(group)] <- NA_character_
  group
}

.posmatchr_motif_names <- function(motif, motif_name = NULL) {
  motif_len <- if (is.character(motif)) length(motif) else if (is.list(motif)) length(motif) else 1L

  if (!is.null(motif_name)) {
    motif_name <- as.character(motif_name)
    if (length(motif_name) == motif_len) return(motif_name)
    if (length(motif_name) == 1L && motif_len == 1L) return(motif_name)
    if (length(motif_name) == 1L && motif_len > 1L) return(paste0(motif_name, "_", seq_len(motif_len)))
    return(motif_name[seq_len(min(length(motif_name), motif_len))])
  }

  if (is.character(motif)) {
    nms <- names(motif)
    if (!is.null(nms) && length(nms) == length(motif) && all(nzchar(nms))) return(as.character(nms))
    return(as.character(motif))
  }

  if (inherits(motif, "universalmotif")) {
    nm <- tryCatch(as.character(motif[["name"]]), error = function(e) NA_character_)
    if (length(nm) && !is.na(nm) && nzchar(nm)) return(nm)
    return("motif")
  }

  if (is.list(motif)) {
    nms <- names(motif)
    if (!is.null(nms) && length(nms) == length(motif) && all(nzchar(nms))) return(nms)
    out <- vapply(seq_along(motif), function(i) {
      mi <- motif[[i]]
      nm <- tryCatch(as.character(mi[["name"]]), error = function(e) NA_character_)
      if (length(nm) && !is.na(nm) && nzchar(nm)) nm else paste0("motif_", i)
    }, character(1))
    return(out)
  }

  "motif"
}

.posmatchr_motif_method <- function(motif, method) {
  if (method != "auto") return(method)
  if (is.character(motif)) return("iupac")
  "universalmotif"
}

.posmatchr_normalize_iupac_patterns <- function(patterns) {
  patterns <- as.character(patterns)
  patterns <- toupper(gsub("\\s+", "", patterns))
  patterns <- chartr("U", "T", patterns)

  ok <- !is.na(patterns) & nzchar(patterns)
  if (!all(ok)) {
    patterns <- patterns[ok]
  }
  if (!length(patterns)) stop("No non-empty motif patterns supplied.")

  allowed <- c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N")
  bad <- vapply(patterns, function(p) {
    chars <- strsplit(p, "", fixed = TRUE)[[1L]]
    any(!(chars %in% allowed))
  }, logical(1))
  if (any(bad)) {
    stop(
      "Motif pattern(s) contain non-IUPAC DNA/RNA symbols after U->T conversion: ",
      paste(patterns[bad], collapse = ", "),
      ". Allowed symbols are A,C,G,T,U,R,Y,S,W,K,M,B,D,H,V,N."
    )
  }

  patterns
}

.posmatchr_make_universalmotif_from_character <- function(motif, motif_name = NULL) {
  if (!requireNamespace("universalmotif", quietly = TRUE)) {
    stop("Package 'universalmotif' is required for universalmotif scanning. Install it with BiocManager::install('universalmotif').")
  }
  pats <- .posmatchr_normalize_iupac_patterns(motif)
  nms <- .posmatchr_motif_names(motif, motif_name = motif_name)
  if (length(nms) != length(pats)) nms <- pats
  out <- lapply(seq_along(pats), function(i) {
    universalmotif::create_motif(pats[i], alphabet = "DNA", name = nms[i])
  })
  if (length(out) == 1L) out[[1L]] else out
}

.posmatchr_hit_position <- function(start, end, hit_position = c("center", "start", "end")) {
  hit_position <- match.arg(hit_position)
  start <- as.numeric(start)
  end <- as.numeric(end)
  switch(
    hit_position,
    start = start,
    end = end,
    center = floor((start + end) / 2),
    stop("Unsupported hit_position.")
  )
}

.posmatchr_prepare_genomic_motif_windows <- function(gr,
                                                     genome,
                                                     window,
                                                     seqstyle = NULL,
                                                     chrs = NULL,
                                                     drop_unstranded = TRUE) {
  gr <- prepare_sites(gr, label = NULL, strip_mcols = FALSE)
  window <- as.integer(window)
  if (!is.finite(window) || window < 1L) stop("`window` must be a positive integer.")

  if (isTRUE(drop_unstranded)) {
    st0 <- as.character(GenomicRanges::strand(gr))
    keep_st <- st0 %in% c("+", "-")
    if (!any(keep_st)) stop("No '+' or '-' stranded sites remain for strand-oriented motif scanning.")
    gr <- gr[keep_st]
  }

  if (!is.null(seqstyle)) {
    suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- seqstyle, silent = TRUE))
    suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(genome) <- seqstyle, silent = TRUE))
  }

  genome_levels <- GenomeInfoDb::seqlevels(genome)
  if (!length(genome_levels)) stop("`genome` has no seqlevels.")

  common <- intersect(GenomeInfoDb::seqlevels(gr), genome_levels)
  if (!length(common) && is.null(seqstyle)) {
    st <- tryCatch(GenomeInfoDb::seqlevelsStyle(genome), error = function(e) character(0))
    if (length(st)) suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- st[1], silent = TRUE))
    common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(genome))
  }
  if (!length(common)) {
    gr <- .rename_seqlevels_to_target(gr, GenomeInfoDb::seqlevels(genome))
    common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(genome))
  }

  if (!is.null(chrs)) {
    target_common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(genome))
    common <- .map_to_target_seqlevels(chrs, target_common)$mapped
  }
  common <- unique(common[!is.na(common) & nzchar(common)])
  if (!length(common)) {
    stop(
      "No common seqlevels between `gr` and `genome`.\n",
      "seqlevels(gr): ", paste(GenomeInfoDb::seqlevels(gr), collapse = ", "), "\n",
      "seqlevels(genome): ", paste(GenomeInfoDb::seqlevels(genome), collapse = ", ")
    )
  }

  gr <- GenomeInfoDb::keepSeqlevels(gr, common, pruning.mode = "coarse")
  try({
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(genome)[GenomeInfoDb::seqlevels(gr)]
  }, silent = TRUE)

  win <- IRanges::resize(gr, width = 2L * window + 1L, fix = "center")
  seqlen <- GenomeInfoDb::seqlengths(win)[as.character(GenomicRanges::seqnames(win))]
  valid <- GenomicRanges::start(win) >= 1L
  valid <- valid & (is.na(seqlen) | GenomicRanges::end(win) <= seqlen)

  if (!any(valid)) stop("No full-width motif windows could be extracted.")
  seqs <- Biostrings::getSeq(genome, win[valid])
  names(seqs) <- paste0("site_", which(valid))
  list(gr = gr, valid = valid, seqs = seqs, mode = "genomic")
}

.posmatchr_scan_iupac <- function(seqs, patterns, window, scan_rc = FALSE, hit_position = "center", motif_names = NULL) {
  original_patterns <- as.character(patterns)
  ok_patterns <- !is.na(original_patterns) & nzchar(gsub("\\s+", "", original_patterns))
  patterns <- .posmatchr_normalize_iupac_patterns(original_patterns)
  if (is.null(motif_names)) motif_names <- patterns
  motif_names <- as.character(motif_names)
  if (length(motif_names) == length(original_patterns)) motif_names <- motif_names[ok_patterns]
  if (length(motif_names) != length(patterns)) motif_names <- patterns

  out <- list()
  out_i <- 0L
  site_names <- names(seqs)

  for (i in seq_along(seqs)) {
    seq_i <- seqs[[i]]
    for (j in seq_along(patterns)) {
      pat <- Biostrings::DNAString(patterns[j])
      hits <- Biostrings::matchPattern(pat, seq_i, fixed = "subject")
      if (length(hits)) {
        st <- IRanges::start(hits)
        en <- IRanges::end(hits)
        pos <- .posmatchr_hit_position(st, en, hit_position)
        out_i <- out_i + 1L
        out[[out_i]] <- data.frame(
          site_name = site_names[i],
          motif = motif_names[j],
          relative_position = pos - (window + 1L),
          strand = "+",
          stringsAsFactors = FALSE
        )
      }

      if (isTRUE(scan_rc)) {
        rc <- Biostrings::reverseComplement(pat)
        hits <- Biostrings::matchPattern(rc, seq_i, fixed = "subject")
        if (length(hits)) {
          st <- IRanges::start(hits)
          en <- IRanges::end(hits)
          pos <- .posmatchr_hit_position(st, en, hit_position)
          out_i <- out_i + 1L
          out[[out_i]] <- data.frame(
            site_name = site_names[i],
            motif = motif_names[j],
            relative_position = pos - (window + 1L),
            strand = "-",
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (!length(out)) {
    return(data.frame(site_name = character(0), motif = character(0), relative_position = numeric(0), strand = character(0)))
  }
  do.call(rbind, out)
}

.posmatchr_scan_universalmotif <- function(seqs,
                                           motif,
                                           window,
                                           motif_name = NULL,
                                           scan_rc = FALSE,
                                           hit_position = "center",
                                           threshold = 0.8,
                                           threshold_type = "logodds",
                                           nthreads = 1L) {
  if (!requireNamespace("universalmotif", quietly = TRUE)) {
    stop("Package 'universalmotif' is required for universalmotif scanning. Install it with BiocManager::install('universalmotif').")
  }
  hits <- universalmotif::scan_sequences(
    motifs = motif,
    sequences = seqs,
    RC = scan_rc,
    threshold = threshold,
    threshold.type = threshold_type,
    return.granges = TRUE,
    nthreads = nthreads
  )

  if (!length(hits)) {
    return(data.frame(site_name = character(0), motif = character(0), relative_position = numeric(0), strand = character(0)))
  }

  mc <- S4Vectors::mcols(hits)
  motif_col <- intersect(c("motif", "name", "motif.name", "motif_name"), colnames(mc))
  if (length(motif_col)) {
    motif_nm <- as.character(mc[[motif_col[1]]])
  } else {
    nms <- .posmatchr_motif_names(motif, motif_name = motif_name)
    motif_nm <- rep(nms[1L], length(hits))
  }
  strand_col <- as.character(GenomicRanges::strand(hits))

  pos <- .posmatchr_hit_position(GenomicRanges::start(hits), GenomicRanges::end(hits), hit_position)
  data.frame(
    site_name = as.character(GenomicRanges::seqnames(hits)),
    motif = motif_nm,
    relative_position = pos - (window + 1L),
    strand = strand_col,
    stringsAsFactors = FALSE
  )
}

.posmatchr_bin_relative_positions <- function(x, window, bin_size) {
  bin_size <- as.integer(bin_size)
  if (!is.finite(bin_size) || bin_size < 1L) stop("`bin_size` must be a positive integer.")
  if (bin_size == 1L) return(as.numeric(round(x)))
  bins <- seq(-window, window, by = bin_size)
  if (utils::tail(bins, 1L) < window) bins <- c(bins, window)
  idx <- findInterval(x, bins, all.inside = TRUE)
  bins[idx]
}

.posmatchr_edge_trim <- function(drop_edge_positions = FALSE, edge_trim = NULL) {
  if (!is.null(edge_trim)) {
    out <- as.integer(edge_trim[1L])
  } else if (is.numeric(drop_edge_positions) && !is.logical(drop_edge_positions)) {
    out <- as.integer(drop_edge_positions[1L])
  } else if (isTRUE(drop_edge_positions)) {
    out <- 2L
  } else {
    out <- 0L
  }
  if (!is.finite(out) || out < 0L) stop("`edge_trim` must be a non-negative integer.")
  out
}

#' Compute a motif-enrichment profile around single-nucleotide sites
#'
#' Counts motif hits in strand-oriented windows centred on each site and returns
#' a per-position profile. The default uses contiguous genomic windows extracted
#' directly from the input \code{GRanges} with \code{BSgenome::getSeq()};
#' because the strand is retained, reverse-strand sites are returned in their
#' RNA/transcript orientation. Character motifs are treated as exact/IUPAC DNA
#' or RNA patterns; RNA U is converted to DNA T. IUPAC ambiguity is
#' interpreted in the motif pattern while ambiguous genome letters are treated
#' as fixed subject letters. Set \code{method = "universalmotif"} to scan
#' consensus/PWM motifs with
#' \code{universalmotif::scan_sequences()}.
#'
#' @param gr A single-nucleotide \code{GRanges}. For matched-set comparisons, pass
#'   \code{subset_matched_pairs(gr)} or use \code{matched_only = TRUE}.
#' @param genome A BSgenome-like object accepted by \code{Biostrings::getSeq()}.
#' @param motif A character vector of exact/IUPAC motifs, or a universalmotif object/list.
#' @param window Number of bp/nt on each side of the site.
#' @param set_col Optional grouping column. Defaults to \code{"match_set"} if present.
#' @param positive_value Value in \code{set_col} marking matched positives.
#' @param negative_value Value in \code{set_col} marking matched negatives.
#' @param include_other Include other groups from \code{set_col}.
#' @param matched_only If TRUE and canonical match columns are present, subset to
#'   reciprocal matched pairs before scanning.
#' @param method \code{"auto"}, \code{"iupac"}, or \code{"universalmotif"}.
#' @param motif_name Optional display name for a supplied motif.
#' @param scan_rc Scan reverse-complement motif occurrences as well. Usually FALSE
#'   for strand-oriented transcript/RNA motifs.
#' @param hit_position Whether each motif hit is represented by its start, end, or centre.
#' @param threshold Threshold passed to \code{universalmotif::scan_sequences()}.
#' @param threshold_type Threshold type passed to \code{universalmotif::scan_sequences()}.
#' @param nthreads Number of threads for universalmotif scanning.
#' @param bin_size Position bin size in bp/nt.
#' @param smooth_window Optional moving-average smoothing width in bins. Use 1 for no smoothing.
#' @param seqstyle Optional seqlevel style passed to the internal sequence-extraction step.
#' @param chrs Optional seqlevels to keep for sequence extraction. Usually not needed if `gr` and `genome` already have matching seqlevels.
#' @param drop_unstranded If TRUE, omit sites on '*' strand because RNA-oriented
#'   motif profiles require a transcript strand.
#' @param drop_edge_positions Logical or integer. If TRUE, omit edge positions
#'   from the profile; if FALSE, keep them. For backward compatibility, an integer
#'   value is interpreted as the number of bp/nt to trim from each edge.
#' @param edge_trim Integer number of bp/nt to omit from each edge of the profile.
#' @param center_exclude Integer number of bp/nt around the central site to omit
#'   from the profile. Defaults to 0.
#' @param normalise_per Scale factor for the returned normalised rates. For example, 1000 reports motif hits per 1000 sites.
#'
#' @return A data.frame with relative position, group, motif, counts, site counts,
#'   and hits per site.
#' @noRd
# Internal implementation used by site_enrichment_profile().
motif_enrichment_profile <- function(gr,
                                     genome,
                                     motif,
                                     window = 250L,
                                     set_col = NULL,
                                     positive_value = "positive",
                                     negative_value = "matched_negative",
                                     include_other = FALSE,
                                     matched_only = FALSE,
                                     method = c("auto", "iupac", "universalmotif"),
                                     motif_name = NULL,
                                     scan_rc = FALSE,
                                     hit_position = c("center", "start", "end"),
                                     threshold = 0.8,
                                     threshold_type = "logodds",
                                     nthreads = 1L,
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
  method <- match.arg(method)
  method <- .posmatchr_motif_method(motif, method)
  hit_position <- match.arg(hit_position)
  window <- as.integer(window)
  bin_size <- as.integer(bin_size)
  smooth_window <- as.integer(smooth_window)
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

  tmp_row_col <- ".posmatchr_motif_row"
  S4Vectors::mcols(gr)[[tmp_row_col]] <- seq_along(gr)
  prep <- .posmatchr_prepare_genomic_motif_windows(
    gr,
    genome = genome,
    window = window,
    seqstyle = seqstyle,
    chrs = chrs,
    drop_unstranded = drop_unstranded
  )

  valid <- prep$valid
  seqs <- prep$seqs
  row_index <- as.integer(S4Vectors::mcols(prep$gr)[[tmp_row_col]])
  group_valid <- group[row_index[valid]]
  site_lookup <- data.frame(site_name = names(seqs), group = group_valid, stringsAsFactors = FALSE)

  if (method == "iupac") {
    motif_levels <- .posmatchr_motif_names(motif, motif_name = motif_name)
    hits <- .posmatchr_scan_iupac(
      seqs,
      motif,
      window = window,
      scan_rc = scan_rc,
      hit_position = hit_position,
      motif_names = motif_levels
    )
  } else {
    if (is.character(motif)) {
      motif <- .posmatchr_make_universalmotif_from_character(motif, motif_name = motif_name)
    }
    hits <- .posmatchr_scan_universalmotif(
      seqs,
      motif = motif,
      window = window,
      motif_name = motif_name,
      scan_rc = scan_rc,
      hit_position = hit_position,
      threshold = threshold,
      threshold_type = threshold_type,
      nthreads = nthreads
    )
    motif_levels <- unique(.posmatchr_motif_names(motif, motif_name = motif_name))
    if (nrow(hits)) motif_levels <- unique(c(motif_levels, hits$motif))
  }

  group_levels <- unique(group_valid)
  positions <- seq(-window, window, by = bin_size)
  if (!length(positions) || utils::tail(positions, 1L) != window) positions <- unique(c(positions, window))
  if (edge_trim > 0L) {
    positions <- positions[positions >= (-window + edge_trim) & positions <= (window - edge_trim)]
  }
  if (center_exclude > 0L) {
    positions <- positions[abs(positions) > center_exclude]
  }

  n_sites <- as.data.frame(table(group = group_valid), stringsAsFactors = FALSE)
  names(n_sites) <- c("group", "n_sites")
  n_sites$n_sites <- as.integer(n_sites$n_sites)

  grid <- expand.grid(
    relative_position = positions,
    group = group_levels,
    motif = motif_levels,
    stringsAsFactors = FALSE
  )

  if (nrow(hits)) {
    hits <- merge(hits, site_lookup, by = "site_name", all.x = FALSE, all.y = FALSE)
    hits <- hits[is.finite(hits$relative_position) & hits$relative_position >= -window & hits$relative_position <= window, , drop = FALSE]
    if (edge_trim > 0L && nrow(hits)) {
      hits <- hits[hits$relative_position >= (-window + edge_trim) & hits$relative_position <= (window - edge_trim), , drop = FALSE]
    }
    if (center_exclude > 0L && nrow(hits)) {
      hits <- hits[abs(hits$relative_position) > center_exclude, , drop = FALSE]
    }
    if (nrow(hits)) {
      hits$relative_position <- .posmatchr_bin_relative_positions(hits$relative_position, window = window, bin_size = bin_size)
      agg <- stats::aggregate(
        count ~ relative_position + group + motif,
        data = transform(hits, count = 1L),
        FUN = sum
      )
    } else {
      agg <- data.frame(relative_position = numeric(0), group = character(0), motif = character(0), count = integer(0))
    }
  } else {
    agg <- data.frame(relative_position = numeric(0), group = character(0), motif = character(0), count = integer(0))
  }

  prof <- merge(grid, agg, by = c("relative_position", "group", "motif"), all.x = TRUE, sort = FALSE)
  prof$count[is.na(prof$count)] <- 0L
  prof <- merge(prof, n_sites, by = "group", all.x = TRUE, sort = FALSE)
  prof$hits_per_site <- prof$count / prof$n_sites
  prof$hits_per_n_sites <- prof$hits_per_site * normalise_per
  prof$hits_per_1000_sites <- prof$hits_per_site * 1000
  prof$normalise_per <- normalise_per
  prof$window_mode <- "genomic"
  prof <- prof[order(prof$motif, prof$group, prof$relative_position), , drop = FALSE]

  if (is.finite(smooth_window) && smooth_window > 1L && nrow(prof)) {
    smooth_one <- function(x) {
      as.numeric(stats::filter(x, rep(1 / smooth_window, smooth_window), sides = 2))
    }
    key <- paste(prof$motif, prof$group, sep = "||")
    sm <- prof$hits_per_site
    for (k in unique(key)) {
      idx <- which(key == k)
      vals <- smooth_one(prof$hits_per_site[idx])
      vals[is.na(vals)] <- prof$hits_per_site[idx][is.na(vals)]
      sm[idx] <- vals
    }
    prof$hits_per_site_smoothed <- sm
  } else {
    prof$hits_per_site_smoothed <- prof$hits_per_site
  }
  prof$hits_per_n_sites_smoothed <- prof$hits_per_site_smoothed * normalise_per
  prof$hits_per_1000_sites_smoothed <- prof$hits_per_site_smoothed * 1000

  prof
}

#' Plot motif enrichment around single-nucleotide sites
#'
#' @inheritParams motif_enrichment_profile
#' @param y Value to plot: \code{"hits_per_site"}, \code{"hits_per_site_smoothed"}, or \code{"count"}.
#' @param x_label Optional x-axis label.
#' @param y_label Optional y-axis label.
#'
#' @return A ggplot object.
#' @noRd
# Internal implementation used by plot_site_enrichment().
plot_motif_enrichment <- function(gr,
                                  genome,
                                  motif,
                                  window = 250L,
                                  set_col = NULL,
                                  positive_value = "positive",
                                  negative_value = "matched_negative",
                                  include_other = FALSE,
                                  matched_only = FALSE,
                                  method = c("auto", "iupac", "universalmotif"),
                                  motif_name = NULL,
                                  scan_rc = FALSE,
                                  hit_position = c("center", "start", "end"),
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
                                  y = c("hits_per_n_sites_smoothed", "hits_per_n_sites", "hits_per_site_smoothed", "hits_per_site", "hits_per_1000_sites_smoothed", "hits_per_1000_sites", "count"),
                                  x_label = NULL,
                                  y_label = NULL) {
  y <- match.arg(y)
  method <- match.arg(method)
  hit_position <- match.arg(hit_position)
  prof <- motif_enrichment_profile(
    gr = gr,
    genome = genome,
    motif = motif,
    window = window,
    set_col = set_col,
    positive_value = positive_value,
    negative_value = negative_value,
    include_other = include_other,
    matched_only = matched_only,
    method = method,
    motif_name = motif_name,
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
  prof$y_plot <- prof[[y]]
  if (center_exclude > 0L) {
    prof$plot_segment <- ifelse(prof$relative_position < 0, "left", "right")
  } else {
    prof$plot_segment <- "all"
  }
  ylab <- y_label %||% if (y == "count") {
    "Motif-hit count"
  } else if (grepl("hits_per_n_sites", y)) {
    paste0("Motif hits per ", format(normalise_per, scientific = FALSE, trim = TRUE), " sites")
  } else if (grepl("1000", y)) {
    "Motif hits per 1,000 sites"
  } else {
    "Motif hits per site"
  }
  xlab <- x_label %||% "Distance from site centre (bp; genomic)"

  p <- ggplot2::ggplot(prof, ggplot2::aes(x = relative_position, y = y_plot, colour = group, linetype = motif, group = interaction(group, motif, plot_segment))) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::scale_x_continuous(limits = c(-window, window), breaks = unique(c(-window, -round(window / 2), 0, round(window / 2), window))) +
    ggplot2::labs(x = xlab, y = ylab, colour = NULL, linetype = "Motif") +
    ggplot2::theme_bw()

  if (length(unique(prof$motif)) == 1L) {
    p <- p + ggplot2::guides(linetype = "none")
  }

  p
}
