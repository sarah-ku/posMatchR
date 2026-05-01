#' Annotate single-nucleotide sites with transcript context
#'
#' Assigns each site to a primary transcript/region using \code{VariantAnnotation::locateVariants()},
#' attaches transcript-level metrics, computes metagene coordinates and feature distances,
#' and optionally maps gene IDs to symbols/names via an OrgDb.
#'
#' Existing user metadata are preserved by default. posMatchR-generated columns are replaced
#' unless \code{overwrite = FALSE}.
#'
#' @param gr A single-nucleotide \code{GRanges}. Wider ranges are resized to width 1.
#' @param txdb A \code{TxDb}.
#' @param label_col Optional label column name to normalise if already present.
#' @param tx_select Transcript selection rule for multi-transcript sites: \code{"longest"} or \code{"first"}.
#' @param seqstyle Optional sequence naming style, e.g. \code{"UCSC"}.
#' @param chrs Optional seqlevels to keep. Use this for species-specific canonical chromosomes.
#' @param resources Optional output of \code{\link{build_tx_resources}}.
#' @param drop_unannotated If TRUE, drop sites with \code{location == "unannotated"}.
#' @param orgdb Optional OrgDb for gene symbol/name mapping.
#' @param gene_keytype Keytype for \code{orgdb}; e.g. \code{"ENTREZID"} or \code{"TAIR"}.
#' @param gene_symbol_col OrgDb column for symbols.
#' @param gene_name_col Optional OrgDb column for names/descriptions.
#' @param site_id_col Metadata column used for stable site IDs.
#' @param preserve_mcols If FALSE, strip existing metadata before annotation.
#' @param overwrite If FALSE, error when existing posMatchR annotation columns are present.
#' @param cache Optional cache passed to \code{VariantAnnotation::locateVariants()}.
#' @param quiet If TRUE, suppress non-critical progress messages.
#'
#' @return An annotated \code{GRanges}.
#' @export
annotate_sites <- function(gr,
                           txdb,
                           label_col = "label",
                           tx_select = c("longest", "first"),
                           seqstyle = NULL,
                           chrs = NULL,
                           resources = NULL,
                           drop_unannotated = FALSE,
                           orgdb = NULL,
                           gene_keytype = NULL,
                           gene_symbol_col = "SYMBOL",
                           gene_name_col = NULL,
                           site_id_col = "site_id",
                           preserve_mcols = TRUE,
                           overwrite = TRUE,
                           cache = NULL,
                           quiet = FALSE) {
  tx_select <- match.arg(tx_select)

  gr <- prepare_sites(
    gr,
    label = NULL,
    label_col = label_col,
    site_id_col = site_id_col,
    strip_mcols = !preserve_mcols
  )

  h <- .harmonize_seqlevels(gr, txdb, seqstyle = seqstyle, chrs = chrs)
  gr <- h$gr
  txdb <- h$txdb

  if (!overwrite) {
    existing <- intersect(.posmatchr_annotation_cols(), colnames(S4Vectors::mcols(gr)))
    if (length(existing)) {
      stop(
        "Existing posMatchR annotation columns found and overwrite=FALSE: ",
        paste(existing, collapse = ", ")
      )
    }
  } else {
    keep_site_id <- site_id_col
    gr <- .drop_mcols(gr, setdiff(.posmatchr_annotation_cols(), keep_site_id))
  }

  if (is.null(resources)) {
    resources <- build_tx_resources(txdb, include_introns = TRUE)
  }
  tx_metrics <- resources$tx_metrics
  tx_key <- if (!is.null(resources$tx_key)) resources$tx_key else .tx_key_from_metrics(tx_metrics)

  loc <- suppressWarnings(
    if (is.null(cache)) {
      VariantAnnotation::locateVariants(gr, txdb, VariantAnnotation::AllVariants())
    } else {
      VariantAnnotation::locateVariants(gr, txdb, VariantAnnotation::AllVariants(), cache = cache)
    }
  )

  loc_df <- as.data.frame(loc, stringsAsFactors = FALSE)

  if (!nrow(loc_df)) {
    loc_df <- data.frame(
      QUERYID = integer(0),
      LOCATION = character(0),
      TXID = character(0),
      GENEID = character(0),
      stringsAsFactors = FALSE
    )
  }

  if (!("QUERYID" %in% colnames(loc_df))) {
    stop("locateVariants() output is missing QUERYID; unexpected VariantAnnotation result.")
  }
  if (!("LOCATION" %in% colnames(loc_df))) {
    stop("locateVariants() output is missing LOCATION; unexpected VariantAnnotation result.")
  }

  loc_df$LOCATION <- as.character(loc_df$LOCATION)
  loc_df$TXID <- if ("TXID" %in% colnames(loc_df)) .first_or_na(loc_df$TXID) else NA_character_
  loc_df$GENEID <- if ("GENEID" %in% colnames(loc_df)) .first_or_na(loc_df$GENEID) else NA_character_

  txid_to_name <- stats::setNames(tx_key, as.character(tx_metrics$tx_id))
  loc_df$tx_name <- txid_to_name[as.character(loc_df$TXID)]

  txname_to_len <- stats::setNames(as.numeric(tx_metrics$tx_len), tx_key)
  loc_df$tx_len <- txname_to_len[as.character(loc_df$tx_name)]
  loc_df$.loc_pri <- .location_priority(loc_df$LOCATION)

  if (nrow(loc_df) > 0L) {
    if (tx_select == "longest") {
      ord <- order(loc_df$QUERYID, -loc_df$tx_len, loc_df$.loc_pri, loc_df$tx_name, na.last = TRUE)
    } else {
      ord <- order(loc_df$QUERYID, loc_df$.loc_pri, loc_df$tx_name, na.last = TRUE)
    }
    loc_df <- loc_df[ord, , drop = FALSE]
    loc_sel <- loc_df[!duplicated(loc_df$QUERYID), , drop = FALSE]
  } else {
    loc_sel <- loc_df
  }

  ann <- S4Vectors::DataFrame(
    location = rep(NA_character_, length(gr)),
    feature = rep(NA_character_, length(gr)),
    region_class = rep(NA_character_, length(gr)),
    gene_id = rep(NA_character_, length(gr)),
    tx_id = rep(NA_character_, length(gr)),
    tx_name = rep(NA_character_, length(gr))
  )

  if (nrow(loc_sel) > 0L) {
    m <- match(seq_along(gr), as.integer(loc_sel$QUERYID))
    ok <- !is.na(m)
    ann$location[ok] <- loc_sel$LOCATION[m[ok]]
    ann$gene_id[ok] <- loc_sel$GENEID[m[ok]]
    ann$tx_id[ok] <- loc_sel$TXID[m[ok]]
    ann$tx_name[ok] <- loc_sel$tx_name[m[ok]]
  }

  ann$location[is.na(ann$location) | !nzchar(ann$location)] <- "unannotated"
  ann$feature <- .map_location_to_feature(ann$location)
  ann$region_class <- ann$feature

  tx_match <- match(as.character(ann$tx_name), tx_key)
  attach_cols <- intersect(
    c("nexon", "tx_len", "cds_len", "utr5_len", "utr3_len", "intron_len", "n_intron"),
    colnames(tx_metrics)
  )
  tx_attached <- S4Vectors::make_zero_col_DFrame(length(gr))
  for (cl in attach_cols) {
    tx_attached[[cl]] <- tx_metrics[[cl]][tx_match]
  }

  S4Vectors::mcols(gr) <- cbind(S4Vectors::mcols(gr), ann, tx_attached)

  gr <- add_feature_geometry(gr, resources = resources)

  gr <- .add_gene_symbols(
    gr,
    orgdb = orgdb,
    gene_id_col = "gene_id",
    keytype = gene_keytype,
    symbol_col = gene_symbol_col,
    name_col = gene_name_col,
    out_symbol_col = "gene_symbol",
    out_name_col = "gene_name"
  )

  md <- S4Vectors::metadata(gr)
  md$posMatchR <- list(
    tx_select = tx_select,
    seqlevels_kept = h$chrs,
    seqlevels_dropped = h$dropped,
    annotation_columns = .posmatchr_annotation_cols()
  )
  S4Vectors::metadata(gr) <- md

  if (drop_unannotated) gr <- gr[gr$location != "unannotated"]

  gr
}

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
    stop("orgdb provided but package 'AnnotationDbi' is not installed.")
  }
  if (!(gene_id_col %in% colnames(S4Vectors::mcols(gr)))) return(gr)

  S4Vectors::mcols(gr)[[out_symbol_col]] <- rep(NA_character_, length(gr))
  if (!is.null(name_col)) S4Vectors::mcols(gr)[[out_name_col]] <- rep(NA_character_, length(gr))
  if (length(gr) == 0L) return(gr)

  gids <- as.character(S4Vectors::mcols(gr)[[gene_id_col]])
  keys <- unique(gids)
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (!length(keys)) return(gr)

  kts <- AnnotationDbi::keytypes(orgdb)
  if (is.null(keytype)) {
    keytype <- if ("ENTREZID" %in% kts) "ENTREZID" else kts[1L]
  }
  if (!(keytype %in% kts)) {
    warning("gene_keytype='", keytype, "' is not available in the supplied OrgDb; skipping gene-symbol mapping.")
    return(gr)
  }

  available_cols <- AnnotationDbi::columns(orgdb)
  cols_requested <- unique(c(symbol_col, if (!is.null(name_col)) name_col))
  cols <- intersect(cols_requested, available_cols)
  missing_cols <- setdiff(cols_requested, available_cols)
  if (length(missing_cols)) {
    warning("OrgDb column(s) not available and skipped: ", paste(missing_cols, collapse = ", "))
  }
  if (!length(cols)) return(gr)

  sel <- tryCatch(
    AnnotationDbi::select(orgdb, keys = keys, columns = cols, keytype = keytype),
    error = function(e) {
      warning("OrgDb mapping failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(sel) || !nrow(sel)) return(gr)

  sel[[keytype]] <- as.character(sel[[keytype]])
  sel <- sel[!is.na(sel[[keytype]]) & nzchar(sel[[keytype]]), , drop = FALSE]
  if (!nrow(sel)) return(gr)

  pick_first <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[1L] else NA_character_
  }

  if (symbol_col %in% colnames(sel)) {
    sym_map <- tapply(sel[[symbol_col]], sel[[keytype]], pick_first)
    S4Vectors::mcols(gr)[[out_symbol_col]] <- unname(as.character(sym_map[gids]))
  }

  if (!is.null(name_col) && name_col %in% colnames(sel)) {
    nm_map <- tapply(sel[[name_col]], sel[[keytype]], pick_first)
    S4Vectors::mcols(gr)[[out_name_col]] <- unname(as.character(nm_map[gids]))
  }

  gr
}

#' Compute within-feature geometry and metagene coordinates
#'
#' Adds feature proportions, equal-width metagene coordinates, median/mean split metagene
#' coordinates, transcript positions, distances to CDS start/stop, and distances to the
#' nearest exon junction where available.
#'
#' @param gr A \code{GRanges} with \code{location} and \code{tx_name} columns.
#' @param resources Output of \code{\link{build_tx_resources}}.
#'
#' @return The input \code{GRanges} with geometry columns added/replaced.
#' @export
add_feature_geometry <- function(gr, resources) {
  if (!("location" %in% colnames(S4Vectors::mcols(gr)))) {
    stop("`gr` must have a 'location' column; run annotate_sites() first.")
  }
  if (!("tx_name" %in% colnames(S4Vectors::mcols(gr)))) {
    stop("`gr` must have a 'tx_name' column; run annotate_sites() first.")
  }

  mc <- S4Vectors::mcols(gr)
  mc[["feature_len"]] <- rep(NA_real_, length(gr))
  mc[["feature_prop"]] <- rep(NA_real_, length(gr))
  mc[["metagene_prop"]] <- rep(NA_real_, length(gr))
  mc[["metagene_split"]] <- rep(NA_real_, length(gr))
  mc[["metagene_split3"]] <- rep(NA_real_, length(gr))
  mc[["isoform_weight"]] <- rep(1.0, length(gr))
  mc[["feature_width"]] <- rep(NA_integer_, length(gr))
  mc[["segment_rank"]] <- rep(NA_integer_, length(gr))
  mc[["dist_from_feature_start"]] <- rep(NA_integer_, length(gr))
  mc[["dist_from_feature_end"]] <- rep(NA_integer_, length(gr))
  mc[["start_dist"]] <- rep(NA_integer_, length(gr))
  mc[["stop_dist"]] <- rep(NA_integer_, length(gr))
  mc[["start_dist_tx"]] <- rep(NA_real_, length(gr))
  mc[["stop_dist_tx"]] <- rep(NA_real_, length(gr))
  mc[["tx_pos"]] <- rep(NA_real_, length(gr))
  mc[["dist_to_tx_start"]] <- rep(NA_real_, length(gr))
  mc[["dist_to_tx_end"]] <- rep(NA_real_, length(gr))
  mc[["nearest_exon_junction_dist"]] <- rep(NA_real_, length(gr))
  S4Vectors::mcols(gr) <- mc

  .fill_one <- function(gr, idx, grl, flen_vec) {
    if (!length(idx)) return(gr)
    if (is.null(grl) || !length(grl)) return(gr)

    sites <- gr[idx]
    txs <- unique(as.character(sites$tx_name))
    txs <- txs[!is.na(txs) & nzchar(txs)]
    txs <- intersect(txs, names(grl))
    if (!length(txs)) return(gr)

    grl_sub <- grl[txs]
    nseg <- S4Vectors::elementNROWS(grl_sub)
    if (!sum(nseg)) return(gr)

    gr$feature_len[idx] <- as.numeric(flen_vec[as.character(sites$tx_name)])

    feat_gr <- unlist(grl_sub, use.names = FALSE)
    feat_gr$tx_name <- rep(names(grl_sub), nseg)
    if ("exon_rank" %in% colnames(S4Vectors::mcols(feat_gr))) {
      feat_gr$segment_rank <- as.integer(feat_gr$exon_rank)
    } else {
      feat_gr$segment_rank <- sequence(nseg)
    }

    hits <- GenomicRanges::findOverlaps(sites, feat_gr, ignore.strand = FALSE)
    if (length(hits) > 0L) {
      qh <- S4Vectors::queryHits(hits)
      sh <- S4Vectors::subjectHits(hits)
      keep <- as.character(sites$tx_name[qh]) == as.character(feat_gr$tx_name[sh])
      qh <- qh[keep]
      sh <- sh[keep]

      if (length(qh) > 0L) {
        ord <- order(qh)
        qh <- qh[ord]
        sh <- sh[ord]
        keep_first <- !duplicated(qh)
        qh <- qh[keep_first]
        sh <- sh[keep_first]

        seg <- feat_gr[sh]
        site_sub_idx <- idx[qh]
        gr$feature_width[site_sub_idx] <- GenomicRanges::width(seg)
        gr$segment_rank[site_sub_idx] <- seg$segment_rank

        pos <- GenomicRanges::start(gr[site_sub_idx])
        seg_start <- GenomicRanges::start(seg)
        seg_end <- GenomicRanges::end(seg)
        st <- as.character(GenomicRanges::strand(gr[site_sub_idx]))

        plus <- st != "-"
        dist_start <- integer(length(pos))
        dist_end <- integer(length(pos))

        dist_start[plus] <- pmax(0L, pos[plus] - seg_start[plus])
        dist_end[plus] <- pmax(0L, seg_end[plus] - pos[plus])

        minus <- !plus
        dist_start[minus] <- pmax(0L, seg_end[minus] - pos[minus])
        dist_end[minus] <- pmax(0L, pos[minus] - seg_start[minus])

        gr$dist_from_feature_start[site_sub_idx] <- dist_start
        gr$dist_from_feature_end[site_sub_idx] <- dist_end
      }
    }

    mapped <- suppressWarnings(
      GenomicFeatures::mapToTranscripts(sites, transcripts = grl_sub, ignore.strand = FALSE)
    )

    if (length(mapped) > 0L) {
      q_idx <- mapped$xHits
      mapped_tx <- as.character(GenomicRanges::seqnames(mapped))
      q_tx <- as.character(sites$tx_name[q_idx])
      keep <- !is.na(q_tx) & nzchar(q_tx) & mapped_tx == q_tx
      mapped <- mapped[keep]
    }

    if (length(mapped) > 0L) {
      mapped_tx <- as.character(GenomicRanges::seqnames(mapped))
      fl <- as.numeric(flen_vec[mapped_tx])

      prop <- (GenomicRanges::start(mapped) - 0.5) / fl
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

  gr <- .fill_one(gr, which(gr$location == "fiveUTR"), resources$five_by_tx, resources$five_len)
  gr <- .fill_one(gr, which(gr$location == "coding"), resources$cds_by_tx, resources$cds_len)
  gr <- .fill_one(gr, which(gr$location == "threeUTR"), resources$three_by_tx, resources$three_len)
  gr <- .fill_one(gr, which(gr$location == "intron"), resources$introns_by_tx, resources$intr_len)

  fp <- gr$feature_prop
  i5 <- which(gr$location == "fiveUTR" & is.finite(fp))
  ic <- which(gr$location == "coding" & is.finite(fp))
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
        coding_end = unname(splits["fiveUTR"] + splits["coding"])
      )
    }
  }

  if (length(i5)) gr$metagene_split[i5] <- fp[i5] * unname(splits["fiveUTR"])
  if (length(ic)) gr$metagene_split[ic] <- unname(splits["fiveUTR"]) + fp[ic] * unname(splits["coding"])
  if (length(i3)) {
    gr$metagene_split[i3] <- unname(splits["fiveUTR"] + splits["coding"]) +
      fp[i3] * unname(splits["threeUTR"])
  }
  ok_ms <- is.finite(gr$metagene_split)
  gr$metagene_split3[ok_ms] <- gr$metagene_split[ok_ms] * 3

  gr <- .add_transcript_position_geometry(gr, resources)

  rid <- if ("record_id" %in% colnames(S4Vectors::mcols(gr))) as.character(gr$record_id) else names(gr)
  w <- rep(1.0, length(gr))
  ok_r <- !is.na(rid) & nzchar(rid)
  if (any(ok_r)) {
    k <- as.numeric(stats::ave(rid[ok_r], rid[ok_r], FUN = length))
    k[!is.finite(k) | k <= 0] <- NA_real_
    ww <- 1 / k
    ww[!is.finite(ww)] <- 1.0
    w[ok_r] <- ww
  }
  gr$isoform_weight <- w

  tx <- as.character(gr$tx_name)
  gr$start_dist <- abs(GenomicRanges::start(gr) - resources$cds_start[tx])
  gr$stop_dist <- abs(GenomicRanges::start(gr) - resources$cds_stop[tx])

  txp <- suppressWarnings(as.numeric(gr$tx_pos))
  five_len <- suppressWarnings(as.numeric(resources$five_len[tx]))
  cds_len <- suppressWarnings(as.numeric(resources$cds_len[tx]))
  cds_start_tx <- five_len + 1
  cds_stop_tx <- five_len + cds_len
  ok_tx_anchor <- is.finite(txp) & is.finite(five_len) & is.finite(cds_len) & cds_len > 0
  gr$start_dist_tx[ok_tx_anchor] <- abs(txp[ok_tx_anchor] - cds_start_tx[ok_tx_anchor])
  gr$stop_dist_tx[ok_tx_anchor] <- abs(txp[ok_tx_anchor] - cds_stop_tx[ok_tx_anchor])

  md <- S4Vectors::metadata(gr)
  md$metagene_splits <- splits
  md$metagene_breaks <- breaks
  md$metagene_breaks3 <- breaks * 3
  S4Vectors::metadata(gr) <- md

  gr
}

.add_transcript_position_geometry <- function(gr, resources) {
  if (is.null(resources$exons_by_tx) || !length(resources$exons_by_tx)) return(gr)

  tx <- as.character(gr$tx_name)
  idx <- which(!is.na(tx) & nzchar(tx) & tx %in% names(resources$exons_by_tx))
  if (!length(idx)) return(gr)

  sites <- gr[idx]
  txs <- unique(as.character(sites$tx_name))
  txs <- txs[!is.na(txs) & nzchar(txs)]
  grl_sub <- resources$exons_by_tx[intersect(txs, names(resources$exons_by_tx))]
  if (!length(grl_sub)) return(gr)

  mapped <- suppressWarnings(
    GenomicFeatures::mapToTranscripts(sites, transcripts = grl_sub, ignore.strand = FALSE)
  )
  if (!length(mapped)) return(gr)

  q_idx <- mapped$xHits
  mapped_tx <- as.character(GenomicRanges::seqnames(mapped))
  q_tx <- as.character(sites$tx_name[q_idx])
  keep <- !is.na(q_tx) & nzchar(q_tx) & mapped_tx == q_tx
  mapped <- mapped[keep]
  if (!length(mapped)) return(gr)

  q_idx <- mapped$xHits
  mapped_tx <- as.character(GenomicRanges::seqnames(mapped))
  pos <- as.numeric(GenomicRanges::start(mapped))

  ord <- order(q_idx)
  q_idx <- q_idx[ord]
  mapped_tx <- mapped_tx[ord]
  pos <- pos[ord]
  keep_first <- !duplicated(q_idx)
  q_idx <- q_idx[keep_first]
  mapped_tx <- mapped_tx[keep_first]
  pos <- pos[keep_first]

  site_idx <- idx[q_idx]
  gr$tx_pos[site_idx] <- pos

  tlen <- as.numeric(resources$tx_len[mapped_tx])
  gr$dist_to_tx_start[site_idx] <- pmax(0, pos - 1)
  gr$dist_to_tx_end[site_idx] <- pmax(0, tlen - pos)

  jd <- vapply(seq_along(pos), function(i) {
    js <- resources$exon_junctions[[mapped_tx[i]]]
    if (is.null(js) || !length(js)) return(NA_real_)
    min(abs(pos[i] - as.numeric(js)), na.rm = TRUE)
  }, numeric(1))
  gr$nearest_exon_junction_dist[site_idx] <- jd

  gr
}

#' Expand sites to all overlapping transcripts
#'
#' Creates a one-row-per-site-per-transcript representation for 5'UTR/CDS/3'UTR
#' or intron features. This is useful for isoform-weighted metagene summaries.
#'
#' @param gr A \code{GRanges}.
#' @param resources Output of \code{\link{build_tx_resources}}.
#' @param regions Feature regions to include.
#' @param ignore.strand Passed to \code{findOverlaps()}.
#' @param record_id_col Optional input metadata column to use as record ID.
#' @param out_record_id_col Output metadata column for record IDs.
#' @param make_unique_names If TRUE, make row names unique by appending transcript IDs.
#' @param drop_unannotated If TRUE, return only overlapping rows.
#'
#' @return Expanded \code{GRanges}.
#' @export
expand_sites_by_transcript <- function(gr,
                                       resources,
                                       regions = c("fiveUTR", "coding", "threeUTR"),
                                       ignore.strand = FALSE,
                                       record_id_col = NULL,
                                       out_record_id_col = "record_id",
                                       make_unique_names = TRUE,
                                       drop_unannotated = TRUE) {
  gr <- prepare_sites(gr, label = NULL, strip_mcols = FALSE)

  rid <- if (!is.null(record_id_col) && record_id_col %in% colnames(S4Vectors::mcols(gr))) {
    as.character(S4Vectors::mcols(gr)[[record_id_col]])
  } else {
    names(gr)
  }
  bad <- is.na(rid) | !nzchar(rid)
  if (any(bad)) rid[bad] <- names(gr)[bad]

  get_grl <- function(region) {
    if (region == "fiveUTR") return(resources$five_by_tx)
    if (region == "coding") return(resources$cds_by_tx)
    if (region == "threeUTR") return(resources$three_by_tx)
    if (region == "intron") return(resources$introns_by_tx)
    stop("Unknown region: ", region)
  }

  collect_hits <- function(gr, grl, location, ignore.strand) {
    if (is.null(grl) || !length(grl)) return(NULL)
    nseg <- S4Vectors::elementNROWS(grl)
    if (!sum(nseg)) return(NULL)

    feat_gr <- unlist(grl, use.names = FALSE)
    feat_gr$tx_name <- rep(names(grl), nseg)

    hits <- GenomicRanges::findOverlaps(gr, feat_gr, ignore.strand = ignore.strand)
    if (!length(hits)) return(NULL)

    q <- S4Vectors::queryHits(hits)
    tx <- as.character(feat_gr$tx_name[S4Vectors::subjectHits(hits)])
    key <- paste(q, tx, location, sep = "\t")
    keep <- !duplicated(key)

    data.frame(
      q = q[keep],
      tx_name = tx[keep],
      location = rep.int(location, sum(keep)),
      priority = .location_priority(rep.int(location, sum(keep))),
      stringsAsFactors = FALSE
    )
  }

  all_hits <- NULL
  for (reg in unique(as.character(regions))) {
    h <- collect_hits(gr, get_grl(reg), location = reg, ignore.strand = ignore.strand)
    if (!is.null(h) && nrow(h)) all_hits <- if (is.null(all_hits)) h else rbind(all_hits, h)
  }

  if (is.null(all_hits) || !nrow(all_hits)) {
    if (drop_unannotated) return(gr[0])
    out <- gr
    S4Vectors::mcols(out)[[out_record_id_col]] <- rid
    S4Vectors::mcols(out)[["tx_name"]] <- NA_character_
    S4Vectors::mcols(out)[["location"]] <- "unannotated"
    return(out)
  }

  all_hits <- all_hits[order(all_hits$q, all_hits$tx_name, all_hits$priority), , drop = FALSE]
  all_hits <- all_hits[!duplicated(paste(all_hits$q, all_hits$tx_name, sep = "\t")), , drop = FALSE]

  out <- gr[all_hits$q]
  S4Vectors::mcols(out)[[out_record_id_col]] <- rid[all_hits$q]
  S4Vectors::mcols(out)[["tx_name"]] <- all_hits$tx_name
  S4Vectors::mcols(out)[["location"]] <- all_hits$location
  S4Vectors::mcols(out)[["feature"]] <- .map_location_to_feature(all_hits$location)
  if (make_unique_names) {
    names(out) <- paste0(S4Vectors::mcols(out)[[out_record_id_col]], "__", S4Vectors::mcols(out)[["tx_name"]])
  }

  out
}

#' Add an isoform-averaged metagene split coordinate
#'
#' Expands each site to overlapping transcript isoforms, computes per-isoform metagene
#' coordinates, and collapses back to one mean coordinate per original site.
#'
#' @param gr A \code{GRanges}.
#' @param resources Output of \code{\link{build_tx_resources}}.
#'
#' @return \code{gr} with \code{metagene_split} and \code{metagene_split3} replaced.
#' @export
add_metagene_split <- function(gr, resources) {
  gr <- prepare_sites(gr, label = NULL, strip_mcols = FALSE)
  S4Vectors::mcols(gr)[["metagene_split"]] <- rep(NA_real_, length(gr))
  S4Vectors::mcols(gr)[["metagene_split3"]] <- rep(NA_real_, length(gr))
  if (!length(gr)) return(gr)

  splits <- resources$metagene_splits
  breaks <- resources$metagene_breaks
  if (is.null(splits) || length(splits) < 3L || any(!c("fiveUTR", "coding", "threeUTR") %in% names(splits))) {
    mg <- .compute_metagene_splits(resources$tx_metrics, stat = "median", require_cds = TRUE)
    splits <- mg$splits
    breaks <- mg$breaks
  } else {
    splits <- splits[c("fiveUTR", "coding", "threeUTR")]
  }

  gr_exp <- expand_sites_by_transcript(
    gr,
    resources = resources,
    regions = c("fiveUTR", "coding", "threeUTR"),
    ignore.strand = FALSE,
    record_id_col = NULL,
    out_record_id_col = "record_id",
    make_unique_names = FALSE,
    drop_unannotated = TRUE
  )
  if (!length(gr_exp)) return(gr)

  gr_exp <- add_feature_geometry(gr_exp, resources = resources)
  rid <- as.character(gr_exp$record_id)
  ms <- as.numeric(gr_exp$metagene_split)
  w <- as.numeric(gr_exp$isoform_weight)

  ok <- is.finite(ms) & is.finite(w) & w > 0 & !is.na(rid) & nzchar(rid)
  if (any(ok)) {
    fac <- factor(rid[ok])
    num <- rowsum(ms[ok] * w[ok], fac, reorder = FALSE)[, 1]
    den <- rowsum(w[ok], fac, reorder = FALSE)[, 1]
    val <- num / den
    names(val) <- rownames(rowsum(w[ok], fac, reorder = FALSE))

    out <- rep(NA_real_, length(gr))
    names(out) <- names(gr)
    hit <- intersect(names(gr), names(val))
    out[hit] <- val[hit]
    S4Vectors::mcols(gr)[["metagene_split"]] <- out
    S4Vectors::mcols(gr)[["metagene_split3"]] <- out * 3
  }

  md <- S4Vectors::metadata(gr)
  md$metagene_splits <- splits
  md$metagene_breaks <- breaks
  md$metagene_breaks3 <- breaks * 3
  S4Vectors::metadata(gr) <- md

  gr
}

#' Add centred k-mer sequence context
#'
#' Extracts a centred k-mer around each site from a genome object supported by
#' \code{Biostrings::getSeq()}. Windows that cannot be extracted at full width are set to NA.
#'
#' @param gr A \code{GRanges}.
#' @param genome A BSgenome or compatible genome object.
#' @param k K-mer width.
#' @param out_col Output metadata column.
#' @param seqstyle Optional sequence naming style.
#' @param chrs Optional seqlevels to keep.
#'
#' @return \code{gr} with \code{out_col} added.
#' @export
add_kmer <- function(gr,
                     genome,
                     k = 5L,
                     out_col = "kmer",
                     seqstyle = NULL,
                     chrs = NULL) {
  gr <- prepare_sites(gr, label = NULL, strip_mcols = FALSE)
  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) stop("`k` must be a positive integer.")
  if ((k %% 2L) == 0L) warning("Even k: the window is not perfectly symmetric around the site.")

  if (!is.null(seqstyle)) {
    suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- seqstyle, silent = TRUE))
    suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(genome) <- seqstyle, silent = TRUE))
  }

  genome_levels <- GenomeInfoDb::seqlevels(genome)
  if (!length(genome_levels)) stop("`genome` has no seqlevels.")

  # Prefer exact names that already match. This is important for TAIR packages,
  # where the TxDb may use 1..5 but BSgenome.Athaliana.TAIR.TAIR9 uses Chr1..Chr5.
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
    mapped <- .map_to_target_seqlevels(chrs, target_common)
    common <- mapped$mapped

    if (!length(common)) {
      mapped_to_genome <- .map_to_target_seqlevels(chrs, GenomeInfoDb::seqlevels(genome))$mapped
      if (length(mapped_to_genome)) {
        gr <- .rename_seqlevels_to_target(gr, mapped_to_genome)
        target_common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(genome))
        mapped <- .map_to_target_seqlevels(mapped_to_genome, target_common)
        common <- mapped$mapped
      }
    }
  } else {
    common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(genome))
  }

  common <- unique(common[!is.na(common) & nzchar(common)])
  if (!length(common)) {
    stop(
      "No common seqlevels between `gr` and `genome`.\n",
      "seqlevels(gr): ", paste(GenomeInfoDb::seqlevels(gr), collapse = ", "), "\n",
      "seqlevels(genome): ", paste(GenomeInfoDb::seqlevels(genome), collapse = ", ")
    )
  }

  observed_before <- unique(as.character(GenomicRanges::seqnames(gr)))
  dropped <- setdiff(observed_before, common)
  if (length(dropped)) {
    warning("Dropping site seqlevels not present in the selected genome/chrs: ", paste(dropped, collapse = ", "))
  }

  gr <- GenomeInfoDb::keepSeqlevels(gr, common, pruning.mode = "coarse")
  try({
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(genome)[GenomeInfoDb::seqlevels(gr)]
  }, silent = TRUE)

  win <- IRanges::resize(gr, width = k, fix = "center")
  seqlen <- GenomeInfoDb::seqlengths(win)[as.character(GenomicRanges::seqnames(win))]
  valid <- GenomicRanges::start(win) >= 1L
  valid <- valid & (is.na(seqlen) | GenomicRanges::end(win) <= seqlen)

  seqs <- rep(NA_character_, length(gr))
  if (any(valid)) {
    got <- as.character(Biostrings::getSeq(genome, win[valid]))
    got[nchar(got) != k] <- NA_character_
    seqs[valid] <- got
  }
  S4Vectors::mcols(gr)[[out_col]] <- seqs
  gr
}
#' Convert annotated sites to a scalar data.frame
#'
#' Converts an annotated \code{GRanges} to a base \code{data.frame}. The default
#' \code{columns = "all"} keeps all scalar metadata columns. \code{columns = "basic"}
#' returns a smaller table intended for quick inspection and reporting.
#'
#' @param gr A \code{GRanges}.
#' @param compatibility_names If TRUE, include column names similar to the user's
#'   earlier project wrapper, such as \code{primary_region_class}.
#' @param columns Either \code{"all"} or \code{"basic"}. \code{"basic"} keeps
#'   coordinates, site ID, labels, region/gene/transcript fields, metagene position,
#'   k-mer and matching identifiers when present.
#'
#' @return A data.frame.
#' @export
as_site_table <- function(gr, compatibility_names = TRUE, columns = c("all", "basic")) {
  columns <- match.arg(columns)
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  mc <- S4Vectors::mcols(gr)
  n <- length(gr)

  df <- data.frame(
    site_id = if ("site_id" %in% colnames(mc)) .mcol_as_character(mc$site_id, n) else as.character(names(gr)),
    seqnames = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    strand = as.character(GenomicRanges::strand(gr)),
    stringsAsFactors = FALSE
  )

  if (columns == "basic") {
    wanted <- c(
      "label", "location", "feature", "region_class", "gene_id", "gene_symbol", "gene_name",
      "tx_id", "tx_name", "transcript", "type", "metagene_prop", "metagene_split",
      "metagene_split3", "nearest_exon_junction_dist", "start_dist", "stop_dist",
      "start_dist_tx", "stop_dist_tx", "kmer", "match_set", "matched_negative_id",
      "matched_positive_id", "match_distance", "split"
    )
  } else {
    wanted <- colnames(mc)
  }

  for (nm in wanted) {
    if (!(nm %in% colnames(mc))) next
    if (nm %in% colnames(df)) next
    x <- mc[[nm]]
    if (methods::is(x, "DataFrame")) next
    if (is.numeric(x) || is.integer(x) || is.logical(x)) {
      df[[nm]] <- x
    } else {
      df[[nm]] <- .mcol_as_character(x, n)
    }
  }

  if (compatibility_names) {
    if (!("primary_region_class" %in% colnames(df))) {
      df$primary_region_class <- if ("feature" %in% colnames(df)) df$feature else NA_character_
    }
    if (!("posmatchr_location" %in% colnames(df))) {
      df$posmatchr_location <- if ("location" %in% colnames(df)) df$location else NA_character_
    }
    if (!("posmatchr_feature" %in% colnames(df))) {
      df$posmatchr_feature <- if ("feature" %in% colnames(df)) df$feature else NA_character_
    }
    if (!("primary_transcript_id" %in% colnames(df))) {
      df$primary_transcript_id <- if ("tx_name" %in% colnames(df)) df$tx_name else if ("tx_id" %in% colnames(df)) df$tx_id else NA_character_
    }
    if (!("primary_gene_id" %in% colnames(df))) {
      df$primary_gene_id <- if ("gene_id" %in% colnames(df)) df$gene_id else NA_character_
    }
    if (!("primary_gene_symbol" %in% colnames(df))) {
      df$primary_gene_symbol <- if ("gene_symbol" %in% colnames(df)) df$gene_symbol else NA_character_
    }
    if (!("primary_gene_name" %in% colnames(df))) {
      df$primary_gene_name <- if ("gene_name" %in% colnames(df)) df$gene_name else NA_character_
    }
  }

  df
}
