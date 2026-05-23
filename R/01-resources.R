#' Compute transcript-level length metrics
#'
#' Computes transcript lengths and optional intron summaries from a \code{TxDb}.
#'
#' @param txdb A \code{TxDb}.
#' @param include_introns If TRUE, compute intron counts and total intron lengths.
#' @param introns_by_tx Optional \code{GRangesList} from \code{GenomicFeatures::intronsByTranscript()}.
#'
#' @examples
#' if (interactive()) {
#'     txdb <- get(
#'         "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'         envir = asNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")
#'     )
#'     metrics <- compute_tx_metrics(txdb)
#'     head(metrics)
#' }
#'
#' @return A data.frame with one row per transcript.
#' @export
compute_tx_metrics <- function(txdb, include_introns = TRUE, introns_by_tx = NULL) {
  tl <- GenomicFeatures::transcriptLengths(
    txdb,
    with.cds_len = TRUE,
    with.utr5_len = TRUE,
    with.utr3_len = TRUE
  )
  tl <- as.data.frame(tl, stringsAsFactors = FALSE)
  tl$tx_id <- as.character(tl$tx_id)
  tl$tx_name <- as.character(tl$tx_name)

  bad_name <- is.na(tl$tx_name) | !nzchar(tl$tx_name)
  if (any(bad_name)) tl$tx_name[bad_name] <- tl$tx_id[bad_name]

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
      intr_sum <- rowsum(GenomicRanges::width(u), grp, reorder = FALSE)
      intr_len[as.integer(rownames(intr_sum))] <- intr_sum[, 1]
    }

    tl$intron_len <- intr_len[match(tl$tx_name, names(intr_len))]
    tl$n_intron <- intr_n[match(tl$tx_name, names(intr_n))]
  }

  tl
}

.compute_metagene_splits <- function(tx_metrics,
                                     stat = c("median", "mean"),
                                     require_cds = TRUE) {
  stat <- match.arg(stat)

  default_splits <- c(fiveUTR = 1 / 3, coding = 1 / 3, threeUTR = 1 / 3)
  default_breaks <- c(fiveUTR_end = 1 / 3, coding_end = 2 / 3)

  fallback <- function(n = 0L) {
    list(
      splits = default_splits,
      breaks = default_breaks,
      medians = list(
        start_codon_pos = NA_real_,
        stop_codon_pos = NA_real_,
        tx_len = NA_real_,
        stat = stat,
        require_cds = require_cds,
        n_transcripts = as.integer(n)
      )
    )
  }

  if (is.null(tx_metrics) || nrow(tx_metrics) == 0L) return(fallback(0L))

  df <- as.data.frame(tx_metrics, stringsAsFactors = FALSE)
  needed <- c("tx_len", "utr5_len", "cds_len")
  missing <- setdiff(needed, names(df))
  if (length(missing)) {
    warning(
      "tx_metrics missing columns needed for metagene split estimation: ",
      paste(missing, collapse = ", "),
      ". Falling back to equal region widths."
    )
    return(fallback(0L))
  }

  utr5 <- as.numeric(df$utr5_len)
  cds <- as.numeric(df$cds_len)
  tx_len <- as.numeric(df$tx_len)
  start_pos <- utr5
  stop_pos <- utr5 + cds

  ok <- is.finite(start_pos) & is.finite(stop_pos) & is.finite(tx_len) & tx_len > 0
  if (require_cds) ok <- ok & is.finite(cds) & cds > 0
  ok <- ok & stop_pos >= start_pos & stop_pos <= tx_len & start_pos <= tx_len

  n_ok <- sum(ok)
  if (!n_ok) {
    warning("No valid transcript length records for metagene split estimation; using equal region widths.")
    return(fallback(0L))
  }

  fun <- if (stat == "median") stats::median else base::mean
  start_est <- fun(start_pos[ok], na.rm = TRUE)
  stop_est <- fun(stop_pos[ok], na.rm = TRUE)
  len_est <- fun(tx_len[ok], na.rm = TRUE)

  L5 <- max(0, start_est)
  LCDS <- max(0, stop_est - start_est)
  L3 <- max(0, len_est - stop_est)
  total <- L5 + LCDS + L3

  if (!is.finite(total) || total <= 0) {
    splits <- default_splits
  } else {
    splits <- c(fiveUTR = L5 / total, coding = LCDS / total, threeUTR = L3 / total)
    splits[!is.finite(splits) | splits < 0] <- 0
    if (sum(splits) <= 0) splits <- default_splits
    splits <- splits / sum(splits)
  }

  breaks <- c(
    fiveUTR_end = unname(splits["fiveUTR"]),
    coding_end = unname(splits["fiveUTR"] + splits["coding"])
  )

  list(
    splits = splits,
    breaks = breaks,
    medians = list(
      start_codon_pos = as.numeric(start_est),
      stop_codon_pos = as.numeric(stop_est),
      tx_len = as.numeric(len_est),
      stat = stat,
      require_cds = require_cds,
      n_transcripts = as.integer(n_ok)
    )
  )
}

.tx_key_from_metrics <- function(tx_metrics) {
  key <- as.character(tx_metrics$tx_name)
  bad <- is.na(key) | !nzchar(key)
  if (any(bad)) key[bad] <- as.character(tx_metrics$tx_id[bad])
  key
}

.make_exon_junction_list <- function(exons_by_tx) {
  out <- vector("list", length(exons_by_tx))
  names(out) <- names(exons_by_tx)
  n <- S4Vectors::elementNROWS(exons_by_tx)
  if (!length(n)) return(out)

  for (i in seq_along(exons_by_tx)) {
    if (n[i] <= 1L) {
      out[[i]] <- numeric(0)
    } else {
      w <- GenomicRanges::width(exons_by_tx[[i]])
      out[[i]] <- as.numeric(cumsum(w)[-length(w)])
    }
  }
  out
}

#' Build transcript resources for repeated site annotation
#'
#' Precomputes transcript lengths, feature-by-transcript structures, CDS anchors,
#' exon-junction coordinates, and metagene-axis splits from a \code{TxDb}.
#'
#' @param txdb A \code{TxDb}.
#' @param include_introns If TRUE, include introns and intron metrics.
#' @param metagene_split_stat Either \code{"median"} or \code{"mean"} for estimating
#'   proportional 5'UTR/CDS/3'UTR axis widths.
#' @param metagene_require_cds If TRUE, use coding transcripts for split estimation.
#'
#' @examples
#' if (interactive()) {
#'     txdb <- get(
#'         "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'         envir = asNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")
#'     )
#'     resources <- build_tx_resources(txdb)
#'     names(resources)
#' }
#'
#' @return A named list of transcript resources.
#' @export
build_tx_resources <- function(txdb,
                               include_introns = TRUE,
                               metagene_split_stat = c("median", "mean"),
                               metagene_require_cds = TRUE) {
  metagene_split_stat <- match.arg(metagene_split_stat)

  tx_gr <- GenomicFeatures::transcripts(txdb, use.names = TRUE)
  tx_strand <- as.character(GenomicRanges::strand(tx_gr))
  names(tx_strand) <- names(tx_gr)

  exons_by_tx <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  five_by_tx <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  three_by_tx <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
  introns_by_tx <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)

  tx_metrics <- compute_tx_metrics(txdb, include_introns = include_introns, introns_by_tx = introns_by_tx)
  tx_key <- .tx_key_from_metrics(tx_metrics)

  tx_len <- stats::setNames(as.numeric(tx_metrics$tx_len), tx_key)
  cds_len <- stats::setNames(as.numeric(tx_metrics$cds_len), tx_key)
  five_len <- stats::setNames(as.numeric(tx_metrics$utr5_len), tx_key)
  three_len <- stats::setNames(as.numeric(tx_metrics$utr3_len), tx_key)
  intr_len <- if ("intron_len" %in% names(tx_metrics)) {
    stats::setNames(as.numeric(tx_metrics$intron_len), tx_key)
  } else {
    stats::setNames(rep(NA_real_, nrow(tx_metrics)), tx_key)
  }

  tx_names <- names(cds_by_tx)
  cds_start <- stats::setNames(rep(NA_integer_, length(tx_names)), tx_names)
  cds_stop <- stats::setNames(rep(NA_integer_, length(tx_names)), tx_names)

  n_cds <- S4Vectors::elementNROWS(cds_by_tx)
  if (sum(n_cds) > 0L) {
    u <- unlist(cds_by_tx, use.names = FALSE)
    offsets <- c(0L, cumsum(n_cds))[seq_along(n_cds)]
    first_idx <- offsets + 1L
    last_idx <- offsets + n_cds

    ok <- n_cds > 0L
    first_seg <- u[first_idx[ok]]
    last_seg <- u[last_idx[ok]]

    st <- tx_strand[tx_names[ok]]
    seg_st <- as.character(GenomicRanges::strand(first_seg))
    st[is.na(st) | st == "*"] <- seg_st[is.na(st) | st == "*"]

    plus <- st != "-"
    minus <- !plus
    ok_names <- tx_names[ok]

    cds_start[ok_names[plus]] <- GenomicRanges::start(first_seg[plus])
    cds_stop[ok_names[plus]] <- GenomicRanges::end(last_seg[plus])
    cds_start[ok_names[minus]] <- GenomicRanges::end(first_seg[minus])
    cds_stop[ok_names[minus]] <- GenomicRanges::start(last_seg[minus])
  }

  mg <- .compute_metagene_splits(
    tx_metrics,
    stat = metagene_split_stat,
    require_cds = metagene_require_cds
  )

  list(
    tx_metrics = tx_metrics,
    tx_key = tx_key,
    tx_len = tx_len,
    tx_strand = tx_strand,
    exons_by_tx = exons_by_tx,
    exon_junctions = .make_exon_junction_list(exons_by_tx),
    cds_by_tx = cds_by_tx,
    five_by_tx = five_by_tx,
    three_by_tx = three_by_tx,
    introns_by_tx = introns_by_tx,
    cds_len = cds_len,
    five_len = five_len,
    three_len = three_len,
    intr_len = intr_len,
    cds_start = cds_start,
    cds_stop = cds_stop,
    metagene_splits = mg$splits,
    metagene_breaks = mg$breaks,
    metagene_medians = mg$medians
  )
}
