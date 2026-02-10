#' @import GenomicRanges
#' @import GenomicFeatures
#' @import VariantAnnotation
#' @importFrom S4Vectors mcols DataFrame elementNROWS
#' @importFrom IRanges resize
#' @importFrom GenomeInfoDb keepSeqlevels seqlevelsStyle seqinfo
#' @importFrom eulerr euler
NULL

# ----------------------------
# Utilities / validation
# ----------------------------

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
  if (!is.null(seqstyle)) {
    GenomeInfoDb::seqlevelsStyle(gr) <- seqstyle
    GenomeInfoDb::seqlevelsStyle(txdb) <- seqstyle
  } else {
    st <- GenomeInfoDb::seqlevelsStyle(txdb)
    if (length(st) > 0L) {
      suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- st[1], silent = TRUE))
    }
  }

  if (is.null(chrs)) {
    chrs <- intersect(seqlevels(gr), seqlevels(txdb))
  }

  gr   <- GenomeInfoDb::keepSeqlevels(gr, chrs, pruning.mode = "coarse")
  txdb <- GenomeInfoDb::keepSeqlevels(txdb, chrs, pruning.mode = "coarse")

  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(txdb)[chrs]
  gr <- trim(gr)

  list(gr = gr, txdb = txdb, chrs = chrs)
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

#' Build transcript resources for fast site annotation
#'
#' Precomputes transcript-level metrics and transcript feature structures (CDS/UTR/introns)
#' from a \code{TxDb}. The returned list can be passed to \code{\link{annotate_sites}}
#' and reused across many GRanges objects to avoid repeated expensive extraction steps.
#'
#' @param txdb A \code{TxDb} object.
#' @param include_introns Logical; if TRUE, include intron GRangesLists and intron metrics.
#'
#' @returns A named \code{list} containing transcript metrics and feature GRangesLists,
#'   including \code{tx_metrics}, \code{cds_by_tx}, \code{five_by_tx}, \code{three_by_tx},
#'   \code{introns_by_tx} (if requested), and length vectors (\code{cds_len}, \code{five_len},
#'   \code{three_len}, \code{intr_len}) plus CDS anchor coordinates (\code{cds_start}, \code{cds_stop}).
#'
#' @examplesIf requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' res <- build_tx_resources(txdb)
#' names(res)
#'
#' @export
build_tx_resources <- function(txdb, include_introns = TRUE) {
  # Transcript strand map (names are transcript names when use.names=TRUE)
  tx_gr <- GenomicFeatures::transcripts(txdb, use.names = TRUE)
  tx_strand <- as.character(strand(tx_gr))
  names(tx_strand) <- names(tx_gr)

  # Feature GRangesLists (do NOT add segment_rank here; too costly globally)
  cds_by_tx     <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  five_by_tx    <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  three_by_tx   <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
  introns_by_tx <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)

  tx_metrics <- compute_tx_metrics(txdb, include_introns = include_introns, introns_by_tx = introns_by_tx)

  # Named length vectors (fast; no vapply(sum(width(.))) over 112k elements)
  cds_len   <- setNames(as.numeric(tx_metrics$cds_len),  tx_metrics$tx_name)
  five_len  <- setNames(as.numeric(tx_metrics$utr5_len), tx_metrics$tx_name)
  three_len <- setNames(as.numeric(tx_metrics$utr3_len), tx_metrics$tx_name)
  intr_len  <- if ("intron_len" %in% names(tx_metrics)) {
    setNames(as.numeric(tx_metrics$intron_len), tx_metrics$tx_name)
  } else {
    setNames(rep(NA_real_, nrow(tx_metrics)), tx_metrics$tx_name)
  }

  # CDS anchors per transcript without looping over 112k transcripts:
  # Use first and last CDS segments in list order. For cdsBy(by="tx"), segments are
  # ordered by ascending exon rank within transcript.
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
    cds_stop      = cds_stop
  )
}

# ----------------------------
# Annotation of sites
# ----------------------------

# gr = testGR
# txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
# label_col = "label"
# tx_select = "longest"
# seqstyle = NULL
# chrs = NULL
# resources = NULL
# drop_unannotated = FALSE
# orgdb = org.Hs.eg.db
# gene_keytype = "ENTREZID"
# gene_symbol_col = "SYMBOL"
# gene_name_col = "GENENAME"
#gene_keytype = NULL
#gene_symbol_col = "SYMBOL"
#gene_name_col = NULL
# NEW: locateVariants cache
#cache = NULL

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
  # VariantAnnotation locations include: coding, fiveUTR, threeUTR, intron,
  # intergenic, spliceSite, promoter. :contentReference[oaicite:5]{index=5}
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

  if (!(gene_id_col %in% colnames(mcols(gr)))) return(gr)

  gids <- as.character(mcols(gr)[[gene_id_col]])
  keys <- unique(gids)
  keys <- keys[!is.na(keys) & nzchar(keys)]

  # Initialize output columns
  mcols(gr)[[out_symbol_col]] <- NA_character_
  if (!is.null(name_col)) mcols(gr)[[out_name_col]] <- NA_character_

  if (!length(keys)) return(gr)

  if (is.null(keytype)) {
    kts <- AnnotationDbi::keytypes(orgdb)
    keytype <- if ("ENTREZID" %in% kts) "ENTREZID" else kts[1]
  }

  cols <- unique(c(symbol_col, if (!is.null(name_col)) name_col))
  sel <- AnnotationDbi::select(orgdb, keys = keys, columns = cols, keytype = keytype)

  # Guardrails: coerce key column to character
  sel[[keytype]] <- as.character(sel[[keytype]])
  sel <- sel[!is.na(sel[[keytype]]) & nzchar(sel[[keytype]]), , drop = FALSE]

  # If multiple rows per key exist, keep the first non-NA per key
  pick_first <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[1] else NA_character_
  }

  sym_map <- tapply(sel[[symbol_col]], sel[[keytype]], pick_first)
  mcols(gr)[[out_symbol_col]] <- unname(sym_map[gids])

  if (!is.null(name_col)) {
    nm_map <- tapply(sel[[name_col]], sel[[keytype]], pick_first)
    mcols(gr)[[out_name_col]] <- unname(nm_map[gids])
  }

  gr
}

# ----------------------------
# Updated annotation function (your version + tweaks)
# ----------------------------


#' Annotate candidate sites with transcript context and region geometry
#'
#' Assigns each site to a transcript and region (e.g. fiveUTR/coding/threeUTR/intron),
#' attaches transcript-level length metrics, computes within-feature coordinates (feature_prop)
#' and a concatenated metagene coordinate (metagene_prop), and optionally maps gene IDs
#' to symbols/names via an \code{OrgDb}.
#'
#' @param gr A \code{GRanges} of single-nucleotide sites. Must contain \code{label_col}.
#' @param txdb A \code{TxDb} object.
#' @param label_col Metadata column in \code{gr} containing 0/1 labels.
#' @param tx_select Which transcript to select if multiple overlap a site.
#'   \code{"longest"} prefers the transcript with the largest \code{tx_len}.
#' @param seqstyle Optional seqlevel style (e.g. \code{"UCSC"}).
#' @param chrs Optional character vector of chromosomes/seqlevels to keep.
#' @param resources Optional output of \code{\link{build_tx_resources}} for reuse/caching.
#' @param drop_unannotated Logical; if TRUE drop rows with \code{location == "unannotated"}.
#' @param orgdb Optional \code{OrgDb} for gene ID mapping (e.g. \code{org.Hs.eg.db}).
#' @param gene_keytype Keytype used by \code{orgdb} for \code{gene_id} (e.g. \code{"ENTREZID"}).
#' @param gene_symbol_col Column to retrieve from \code{orgdb} for gene symbols (default \code{"SYMBOL"}).
#' @param gene_name_col Optional column to retrieve for gene names/descriptions (e.g. \code{"GENENAME"}).
#' @param cache Optional cache passed to \code{VariantAnnotation::locateVariants}.
#'
#' @returns A \code{GRanges} with added metadata including \code{location}, \code{feature},
#'   \code{gene_id}, \code{tx_name}, transcript length metrics, feature geometry columns,
#'   and (if requested) gene symbol/name columns.
#'
#' @seealso \code{\link{build_tx_resources}}, \code{\link{add_kmer}},
#'   \code{\link{match_background}}, \code{\link{random_drach_within_transcript}}
#'
#' @examplesIf requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) && requireNamespace("GenomicRanges", quietly = TRUE)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100000, 100000), strand = "+")
#' S4Vectors::mcols(gr)$label <- 1L
#' ann <- annotate_sites(gr, txdb = txdb)
#'
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
                           cache = NULL) {
  tx_select <- match.arg(tx_select)
  gr <- .validate_input_gr(gr, label_col = label_col, require_both = FALSE)
  h <- .harmonize_seqlevels(gr, txdb, seqstyle = seqstyle, chrs = chrs)
  gr <- h$gr; txdb <- h$txdb

  if (is.null(resources)) {
    resources <- build_tx_resources(txdb, include_introns = TRUE)
  }
  tx_metrics <- resources$tx_metrics

  # locateVariants can be sped up by re-using cache across calls :contentReference[oaicite:6]{index=6}
  loc <- suppressWarnings(
    VariantAnnotation::locateVariants(gr, txdb, VariantAnnotation::AllVariants(), cache = cache)
  )
  names(loc) <- NULL
  loc_df <- as.data.frame(loc, stringsAsFactors = FALSE)

  if (!("QUERYID" %in% colnames(loc_df))) stop("locateVariants output missing QUERYID; unexpected.")
  if (!("LOCATION" %in% colnames(loc_df))) stop("locateVariants output missing LOCATION; unexpected.")

  loc_df$LOCATION <- as.character(loc_df$LOCATION)

  if ("TXID" %in% colnames(loc_df)) loc_df$TXID <- .first_or_na(loc_df$TXID) else loc_df$TXID <- NA_character_
  if ("GENEID" %in% colnames(loc_df)) loc_df$GENEID <- .first_or_na(loc_df$GENEID) else loc_df$GENEID <- NA_character_

  txid_to_name <- setNames(as.character(tx_metrics$tx_name), as.character(tx_metrics$tx_id))
  loc_df$tx_name <- txid_to_name[loc_df$TXID]

  txname_to_len <- setNames(tx_metrics$tx_len, tx_metrics$tx_name)
  loc_df$tx_len <- txname_to_len[loc_df$tx_name]

  # Tie-breaker for rows that differ only by LOCATION (e.g. promoter + intergenic)
  loc_df$.loc_pri <- .location_priority(loc_df$LOCATION)

  if (nrow(loc_df) > 0L) {
    if (tx_select == "longest") {
      ord <- order(loc_df$QUERYID, -loc_df$tx_len, loc_df$.loc_pri, loc_df$tx_name)
    } else {
      ord <- order(loc_df$QUERYID, loc_df$.loc_pri, loc_df$tx_name)
    }
    loc_df <- loc_df[ord, , drop = FALSE]
    loc_sel <- loc_df[!duplicated(loc_df$QUERYID), , drop = FALSE]
  } else {
    loc_sel <- loc_df
  }

  # Ann table aligned to query GRanges
  ann <- S4Vectors::DataFrame(
    location = rep(NA_character_, length(gr)),  # raw VariantAnnotation term
    feature  = rep(NA_character_, length(gr)),  # collapsed label (CDS/5UTR/3UTR/...)
    gene_id  = rep(NA_character_, length(gr)),
    tx_id    = rep(NA_character_, length(gr)),
    tx_name  = rep(NA_character_, length(gr))
  )

  if (nrow(loc_sel) > 0L) {
    m <- match(seq_along(gr), loc_sel$QUERYID)
    ok <- !is.na(m)

    ann$location[ok] <- loc_sel$LOCATION[m[ok]]
    ann$gene_id[ok]  <- loc_sel$GENEID[m[ok]]
    ann$tx_id[ok]    <- loc_sel$TXID[m[ok]]
    ann$tx_name[ok]  <- loc_sel$tx_name[m[ok]]
  }

  ann$location[is.na(ann$location)] <- "unannotated"
  ann$feature <- .map_location_to_feature(ann$location)

  # Attach transcript-level metrics (tx_len, cds_len, utr lengths, etc.)
  txm <- tx_metrics
  # Avoid duplicating gene_id in mcols (your output had gene_id twice)
  keep_cols <- setdiff(colnames(txm), c("tx_id", "tx_name", "gene_id"))
  tx_match <- match(ann$tx_name, txm$tx_name)

  tx_attached <- S4Vectors::DataFrame()
  for (cl in keep_cols) {
    tx_attached[,cl] <- txm[[cl]][tx_match]
  }

  mcols(gr) <- cbind(mcols(gr), ann, tx_attached)

  # Add feature geometry + feature_prop / metagene_prop
  gr <- add_feature_geometry(gr, resources = resources)

  # Optional gene symbol mapping (works if you provide an OrgDb)
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

  if (drop_unannotated) {
    gr <- gr[gr$location != "unannotated"]
  }

  gr
}


#' Compute within-feature geometry and metagene coordinates
#'
#' Adds feature-level geometry columns (e.g. \code{feature_len}, \code{feature_prop},
#' \code{feature_width}, segment rank, distances to feature edges) and a concatenated
#' metagene coordinate (\code{metagene_prop}) for fiveUTR/coding/threeUTR.
#'
#' @param gr A \code{GRanges} with at least \code{location} and \code{tx_name} in \code{mcols(gr)}.
#' @param resources Output from \code{\link{build_tx_resources}}.
#'
#' @returns The input \code{GRanges} with additional geometry columns added.
#'
#' @seealso \code{\link{annotate_sites}}, \code{\link{build_tx_resources}}
#' @export
add_feature_geometry <- function(gr, resources) {
  if (!("location" %in% colnames(mcols(gr)))) stop("`gr` must have a 'location' column (run annotate_sites first).")
  if (!("tx_name" %in% colnames(mcols(gr))))  stop("`gr` must have a 'tx_name' column (run annotate_sites first).")

  init_col <- function(gr, nm, val) {
    if (!(nm %in% names(mcols(gr)))) {
      mcols(gr)[[nm]] <- val
    }
    gr
  }

  gr <- init_col(gr, "feature_len",             rep(NA_real_,    length(gr)))
  gr <- init_col(gr, "feature_prop",            rep(NA_real_,    length(gr)))
  gr <- init_col(gr, "metagene_prop",           rep(NA_real_,    length(gr)))
  gr <- init_col(gr, "feature_width",           rep(NA_integer_, length(gr)))
  gr <- init_col(gr, "segment_rank",            rep(NA_integer_, length(gr)))
  gr <- init_col(gr, "dist_from_feature_start", rep(NA_integer_, length(gr)))
  gr <- init_col(gr, "dist_from_feature_end",   rep(NA_integer_, length(gr)))
  gr <- init_col(gr, "start_dist",              rep(NA_integer_, length(gr)))
  gr <- init_col(gr, "stop_dist",               rep(NA_integer_, length(gr)))


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

    # Total feature length for each site's transcript
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
      mapped_tx <- as.character(seqnames(mapped))
      fl <- flen_vec[mapped_tx]
      prop <- start(mapped) / fl

      ok <- is.finite(prop) & prop >= 0 & prop <= 1
      prop[!ok] <- NA_real_

      prop_by_hit <- tapply(prop, mapped$xHits, function(z) mean(z, na.rm = TRUE))
      hit_idx <- as.integer(names(prop_by_hit))
      gr$feature_prop[idx[hit_idx]] <- as.numeric(prop_by_hit)
    }

    gr
  }


  # Fill in geometry for core regions
  gr <- .fill_one(gr, which(gr$location == "fiveUTR"),  resources$five_by_tx,     resources$five_len)
  gr <- .fill_one(gr, which(gr$location == "coding"),   resources$cds_by_tx,      resources$cds_len)
  gr <- .fill_one(gr, which(gr$location == "threeUTR"), resources$three_by_tx,    resources$three_len)
  gr <- .fill_one(gr, which(gr$location == "intron"),   resources$introns_by_tx,  resources$intr_len)

  # Metagene coordinate (useful for classic 5UTR/CDS/3UTR plots)
  i5 <- which(gr$location == "fiveUTR"  & is.finite(gr$feature_prop))
  ic <- which(gr$location == "coding"   & is.finite(gr$feature_prop))
  i3 <- which(gr$location == "threeUTR" & is.finite(gr$feature_prop))

  gr$metagene_prop[i5] <- gr$feature_prop[i5]
  gr$metagene_prop[ic] <- gr$feature_prop[ic] + 1
  gr$metagene_prop[i3] <- gr$feature_prop[i3] + 2

  # Distances to CDS anchors in genomic space (NA for noncoding transcripts)
  tx <- as.character(gr$tx_name)
  gr$start_dist <- abs(start(gr) - resources$cds_start[tx])
  gr$stop_dist  <- abs(start(gr) - resources$cds_stop[tx])

  gr
}


# ----------------------------
# Sequence context (k-mer)
# ----------------------------

#' Add sequence k-mer context around each site
#'
#' Extracts a centered k-mer from a reference genome and stores it as a metadata column.
#' Typically used prior to strict k-mer background matching.
#'
#' @param gr A \code{GRanges} of single-nucleotide sites.
#' @param genome A BSgenome object (or other supported genome object for \code{Biostrings::getSeq}).
#' @param k Integer k-mer length (default 5).
#' @param out_col Name of the metadata column to write (default \code{"kmer"}).
#' @param seqstyle Optional seqlevel style (e.g. \code{"UCSC"}).
#' @param chrs Optional character vector of seqlevels to keep.
#'
#' @returns The input \code{GRanges} with \code{out_col} added to \code{mcols(gr)}.
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
  win <- trim(win)

  mcols(gr)[[out_col]] <- as.character(Biostrings::getSeq(genome, win))
  gr
}


# ----------------------------
# Updated matcher

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
                                    meta_tol = 0.05,              # metagene units (0..3); within a region this is ~within-feature tolerance
                                    enforce_meta_tol = FALSE,     # if TRUE: leave unmatched if no neg within tol
                                    meta_k = 200L,                # local candidate window size around each positive in meta-sorted negatives

                                    # 4) optional "match level bins" (secondary)
                                    bin_match = TRUE,
                                    bin_cols = c("feature_width","segment_rank","nexon","tx_len",
                                                 "dist_from_feature_start","dist_from_feature_end"),
                                    n_bins = 5L,
                                    bin_within = c("location","global"),
                                    bin_mode = c("soft","hard"),
                                    bin_weight = 0.25,            # penalty weight relative to meta closeness (keep small)
                                    log_bin_cols = c("tx_len","feature_width","dist_from_feature_start","dist_from_feature_end"),

                                    seed = 1L) {
  stopifnot(inherits(gr, "GRanges"))
  if (!(label_col %in% colnames(mcols(gr)))) stop("Missing label_col: ", label_col)
  if (!(location_col %in% colnames(mcols(gr)))) stop("Missing location_col: ", location_col)
  if (!(meta_col %in% colnames(mcols(gr)))) stop("Missing meta_col: ", meta_col)
  if (kmer_match && !(kmer_col %in% colnames(mcols(gr)))) stop("kmer_match=TRUE but missing kmer_col: ", kmer_col)

  bin_within <- match.arg(bin_within)
  bin_mode <- match.arg(bin_mode)

  # Ensure stable names (used as IDs)
  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_along(gr))
  }

  y <- mcols(gr)[[label_col]]
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y))  y <- as.integer(as.character(y))
  y <- as.integer(y)

  loc <- as.character(mcols(gr)[[location_col]])
  meta <- suppressWarnings(as.numeric(mcols(gr)[[meta_col]]))

  # One-hot location columns (useful for modeling later)
  mcols(gr)$loc_fiveUTR  <- as.integer(loc == "fiveUTR")
  mcols(gr)$loc_coding   <- as.integer(loc == "coding")
  mcols(gr)$loc_threeUTR <- as.integer(loc == "threeUTR")

  # In-scope sites (only these locations can be matched)
  in_scope <- !is.na(loc) & loc %in% locations

  pos_idx <- which(y == 1L & in_scope & is.finite(meta))
  neg_idx <- which(y == 0L & in_scope & is.finite(meta))

  # If kmer matching is on, require non-missing kmers
  if (kmer_match) {
    km <- as.character(mcols(gr)[[kmer_col]])
    pos_idx <- pos_idx[!is.na(km[pos_idx]) & nzchar(km[pos_idx])]
    neg_idx <- neg_idx[!is.na(km[neg_idx]) & nzchar(km[neg_idx])]
  }

  # Edge case
  if (!length(pos_idx) || !length(neg_idx)) {
    mcols(gr)$is_positive <- as.integer(y == 1L)
    mcols(gr)$is_matched_negative <- 0L
    mcols(gr)$match_set <- ifelse(mcols(gr)$is_positive == 1L, "positive", "other")
    attr(gr, "match_diagnostics") <- list(
      n_pos_in_scope = length(pos_idx),
      n_neg_in_scope = length(neg_idx),
      n_matched = 0L,
      reason = "No eligible positives or negatives in scope after filtering."
    )
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
    # Ensure cut() covers all values
    br[1] <- -Inf
    br[length(br)] <- Inf
    out[ok] <- as.integer(cut(x[ok], breaks = br, include.lowest = TRUE, right = TRUE, labels = FALSE))
    out
  }

  # ---- build bin matrix (optional) ----
  bin_mat <- NULL
  if (bin_match) {
    cols_ok <- intersect(bin_cols, colnames(mcols(gr)))
    if (!length(cols_ok)) stop("bin_match=TRUE but none of bin_cols exist in mcols(gr).")

    # Copy numeric vectors, optionally log1p-transform selected ones
    get_feat <- function(col) {
      v <- suppressWarnings(as.numeric(mcols(gr)[[col]]))
      if (col %in% log_bin_cols) v <- log1p(pmax(v, 0))
      v
    }

    # Compute bins either globally or within each location
    bin_mat <- matrix(NA_integer_, nrow = length(gr), ncol = length(cols_ok))
    colnames(bin_mat) <- cols_ok

    if (bin_within == "global") {
      for (j in seq_along(cols_ok)) {
        bin_mat[, j] <- quantile_bin(get_feat(cols_ok[j]), n_bins = n_bins)
      }
    } else {
      # within location to keep bins interpretable per region
      for (lv in locations) {
        idx_lv <- which(loc == lv & in_scope)
        if (!length(idx_lv)) next
        for (j in seq_along(cols_ok)) {
          v <- get_feat(cols_ok[j])
          b <- rep.int(NA_integer_, length(v))
          b[idx_lv] <- quantile_bin(v[idx_lv], n_bins = n_bins)
          bin_mat[, j] <- ifelse(is.na(bin_mat[, j]), b, bin_mat[, j])
        }
      }
    }

    # Store bins as columns (useful for debugging / modeling)
    for (j in seq_along(cols_ok)) {
      mcols(gr)[[paste0("bin_", cols_ok[j])]] <- bin_mat[, j]
    }
  }

  # ---- hard constraints => group keys ----
  if (kmer_match) {
    km <- as.character(mcols(gr)[[kmer_col]])
    key <- paste0(loc, "||", km)
  } else {
    key <- loc
  }

  # Only consider keys relevant to our filtered idx
  pos_by_key <- split(pos_idx, key[pos_idx], drop = TRUE)
  neg_by_key <- split(neg_idx, key[neg_idx], drop = TRUE)
  keys <- intersect(names(pos_by_key), names(neg_by_key))
  keys <- keys[!is.na(keys) & nzchar(keys)]

  # Track unused negatives globally
  neg_used <- rep(FALSE, length(gr))
  matched_neg_for_pos <- rep(NA_character_, length(gr))  # store negative name for each positive
  matched_pos_for_neg <- rep(NA_character_, length(gr))  # store positive name for each negative
  match_meta_delta <- rep(NA_real_, length(gr))          # abs(meta_pos - meta_neg) for positives

  set.seed(seed)

  # ---- main matching: process each group independently ----
  for (kkey in keys) {
    P <- pos_by_key[[kkey]]
    N <- neg_by_key[[kkey]]
    if (!length(P) || !length(N)) next

    # Sort positives by meta to preserve global metagene structure
    P <- P[order(meta[P])]

    # Sort negatives by meta; we will do local window searches in this order
    N_sorted <- N[order(meta[N])]
    N_meta <- meta[N_sorted]

    for (pi in P) {
      # Skip if already matched (shouldn't happen, but safe)
      if (!is.na(matched_neg_for_pos[pi])) next

      # Available negatives
      avail_mask <- !neg_used[N_sorted]
      if (!any(avail_mask)) break

      # Locate position in sorted negative meta
      x <- meta[pi]
      j <- findInterval(x, N_meta)
      j <- max(1L, min(j, length(N_sorted)))

      # Start with a local window around j
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
          # Optionally enforce meta tolerance
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
            # Primary cost: meta closeness (scaled so meta_tol ~ 1 unit cost)
            meta_scale <- if (!is.null(meta_tol) && is.finite(meta_tol) && meta_tol > 0) meta_tol else 0.05
            cost_meta <- md / meta_scale

            # Secondary cost: bins (either hard constraint or soft penalty)
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
                  # If hard bins eliminate all candidates, either expand window or give up
                  # (This is exactly the situation that can destroy metagene matching in strict mode)
                  cand <- integer(0)
                } else {
                  cand <- cand2
                  cost_meta <- cost_meta2
                  md <- md2
                  bn <- bn[ok, , drop = FALSE]
                }
              }

              if (length(cand)) {
                # soft bin distance in [0,1] roughly: mean(|bin_diff|/(n_bins-1))
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
              o <- order(cost, md)  # tie-break by pure meta closeness
              chosen <- cand[o[1]]
              break
            }
          }
        }

        # No candidate found: expand window
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
  mcols(gr)$is_positive <- as.integer(y == 1L)
  mcols(gr)$matched_negative_id <- matched_neg_for_pos
  mcols(gr)$matched_positive_id <- matched_pos_for_neg
  mcols(gr)$meta_delta <- match_meta_delta

  matched_neg_idx <- which(!is.na(matched_pos_for_neg))
  mcols(gr)$is_matched_negative <- as.integer(seq_along(gr) %in% matched_neg_idx)

  mcols(gr)$match_set <- ifelse(mcols(gr)$is_positive == 1L, "positive",
                                ifelse(mcols(gr)$is_matched_negative == 1L, "matched_negative", "other"))

  if (return_diagnostics) {
    n_pos_scope <- length(which(y == 1L & in_scope))
    n_pos_eligible <- length(pos_idx)
    n_matched <- sum(!is.na(matched_neg_for_pos[pos_idx]))
    diag <- list(
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
      meta_k = as.integer(meta_k)
    )
    attr(gr, "match_diagnostics") <- diag
  }

  gr
}


#' Sample a random DRACH background site within the same transcript
#'
#' For each positive site, attempts to select a unique negative DRACH site from the same
#' transcript, optionally matching location and/or exact k-mer. This yields a transcript-local
#' background set that can be compared to \code{\link{match_background}}.
#'
#' Optionally writes the same standard match columns produced by \code{\link{match_background}}
#' so downstream plotting functions work without modification.
#'
#' @param gr A \code{GRanges} containing annotation columns from \code{\link{annotate_sites}}
#'   and sequence context from \code{\link{add_kmer}}.
#' @param label_col Metadata label column (0/1).
#' @param tx_col Transcript identifier column (default \code{"tx_name"}).
#' @param location_col Region column (default \code{"location"}).
#' @param locations Allowed regions for matching.
#' @param kmer_col Column containing k-mers (default \code{"kmer"}).
#' @param match_location If TRUE require the negative to be in the same region as the positive.
#' @param match_kmer If TRUE require exact k-mer match.
#' @param meta_col Optional metagene coordinate column used to compute \code{meta_delta}.
#' @param seed Random seed controlling sampling within transcript strata.
#' @param prefix Prefix for function-specific output columns.
#' @param write_standard_cols If TRUE write canonical columns used by plotting helpers.
#' @param overwrite_standard_cols If FALSE, do not overwrite canonical columns if they already exist.
#' @param assert_invariants If TRUE perform internal consistency checks.
#' @param return_diagnostics If TRUE attach a \code{<prefix>_diagnostics} attribute.
#'
#' @returns A \code{GRanges} with prefix match columns (and optionally canonical match columns).
#'
#' @seealso \code{\link{match_background}}, \code{\link{subset_matched_pairs}}
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
                                           # NEW: write standard columns used by plot_* and match_background_simple()
                                           write_standard_cols = TRUE,
                                           overwrite_standard_cols = TRUE,
                                           assert_invariants = TRUE,
                                           return_diagnostics = TRUE) {
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

  # DRACH in DNA alphabet: D=[AGT], R=[AG], A, C, H=[ACT]
  is_drach_5mer <- function(x) {
    x <- toupper(as.character(x))
    x <- chartr("U","T", x)
    ok <- !is.na(x) & nchar(x) == 5L
    out <- rep(FALSE, length(x))
    out[ok] <- grepl("^[AGT][AG]AC[ACT]$", x[ok])
    out
  }

  # In-scope rows must have tx/loc/kmer
  in_scope <- !is.na(tx) & nzchar(tx) &
    !is.na(loc) & nzchar(loc) &
    !is.na(km) & nzchar(km)

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
  # IMPORTANT: align values with plot functions + match_background_simple
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

    # Optionally also write canonical columns so plotting works anyway
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
        n_pairs = 0L,
        reason = "No eligible positives or no eligible DRACH negatives in scope."
      )
    }
    return(gr)
  }

  # Matching key builder
  make_key <- function(idx) {
    k <- tx[idx]
    if (match_location) k <- paste0(k, "||", loc[idx])
    if (match_kmer)     k <- paste0(k, "||", km[idx])
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
        n_pairs = 0L,
        reason = "No overlap of matching keys between positives and DRACH negatives."
      )
    }
    return(gr)
  }

  # Sort by key and walk runs (memory-friendly)
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
        # match as many as possible (1:1, no replacement)
        m <- min(length(P), length(N))
        if (m > 0L) {
          # SAFE sampling: sample positions, not index values
          Ppick <- P[sample.int(length(P), m, replace = FALSE)]
          Npick <- N[sample.int(length(N), m, replace = FALSE)]

          # random pairing permutation
          Npick <- Npick[sample.int(length(Npick), length(Npick), replace = FALSE)]

          # Fill pairing columns
          neg_id[Ppick] <- id[Npick]
          pos_id[Npick] <- id[Ppick]
          isneg[Npick]  <- 1L
          setv[Npick]   <- "matched_negative"  # <- aligns with plot functions

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

  # Invariants: positives must never be treated as negatives
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

    # Matched positives should equal matched negatives in 1:1 mode
    n_pos_matched <- sum(!is.na(neg_id[pos_idx]))
    n_neg_matched <- sum(!is.na(pos_id[neg_idx]))
    if (n_pos_matched != n_neg_matched) {
      stop("Invariant violation: matched positives (", n_pos_matched,
           ") != matched negatives (", n_neg_matched, ").")
    }
  }

  # NEW: write canonical columns used by plot_* and match_background_simple()
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
    diag <- list(
      n_pos_total = sum(y == 1L, na.rm = TRUE),
      n_neg_total = sum(y == 0L, na.rm = TRUE),
      n_pos_in_scope = length(pos_all),
      n_neg_drach_in_scope = length(neg_all),
      match_location = match_location,
      match_kmer = match_kmer,
      n_pairs = pairs,
      note = "Prefix columns are written, and (optionally) canonical match_* columns are also written."
    )
    attr(gr, paste0(prefix, "_diagnostics")) <- diag
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
#' (\code{matched_negative_id}, \code{matched_positive_id}) so the output contains
#' exactly one negative per positive for all retained pairs.
#'
#' @param gr A \code{GRanges} with canonical match columns produced by
#'   \code{\link{match_background}} or \code{\link{random_drach_within_transcript}}.
#' @param matched_negative_id_col Column containing the matched negative ID for each positive.
#' @param matched_positive_id_col Column containing the matched positive ID for each negative.
#' @param is_positive_col Optional column indicating positives (1/0). If NULL, uses \code{label_col}.
#' @param label_col Label column (0/1) used when \code{is_positive_col} is NULL.
#' @param strict_reciprocal If TRUE enforce reciprocal links (negative points back to the same positive).
#' @param drop_conflicts If TRUE drop duplicated negative assignments if present.
#' @param return_diagnostics If TRUE attach a \code{"subset_pairs_diagnostics"} attribute.
#'
#' @returns A \code{GRanges} containing only matched pairs (2 rows per pair).
#' @export
subset_matched_pairs <- function(gr,
                                 # canonical columns produced by match_background_simple() and random_drach_within_transcript()
                                 matched_negative_id_col = "matched_negative_id",
                                 matched_positive_id_col = "matched_positive_id",
                                 is_positive_col = NULL,   # if NULL, will fall back to label_col
                                 label_col = "label",
                                 strict_reciprocal = TRUE, # enforce neg->pos points back correctly
                                 drop_conflicts = TRUE,    # drop duplicate neg assignments if present
                                 return_diagnostics = TRUE) {
  stopifnot(inherits(gr, "GRanges"))

  mc <- S4Vectors::mcols(gr)

  if (is.null(names(gr)) || anyDuplicated(names(gr))) {
    names(gr) <- paste0("site_", seq_len(length(gr)))
  }
  ids <- names(gr)

  if (!(matched_negative_id_col %in% colnames(mc))) {
    stop("Missing `", matched_negative_id_col, "` in mcols(gr). Cannot subset to matched pairs.")
  }
  if (strict_reciprocal && !(matched_positive_id_col %in% colnames(mc))) {
    stop("strict_reciprocal=TRUE but missing `", matched_positive_id_col, "` in mcols(gr).")
  }

  # Determine which rows are positives
  if (!is.null(is_positive_col) && (is_positive_col %in% colnames(mc))) {
    pos_flag <- mc[[is_positive_col]]
    if (is.logical(pos_flag)) pos_flag <- as.integer(pos_flag)
    if (is.factor(pos_flag))  pos_flag <- as.integer(as.character(pos_flag))
    pos_flag <- as.integer(pos_flag)
    pos_idx_all <- which(pos_flag == 1L)
  } else {
    if (!(label_col %in% colnames(mc))) stop("Missing `", label_col, "` in mcols(gr).")
    y <- mc[[label_col]]
    if (is.logical(y)) y <- as.integer(y)
    if (is.factor(y))  y <- as.integer(as.character(y))
    y <- as.integer(y)
    pos_idx_all <- which(y == 1L)
  }

  # Build (positive_id -> negative_id) table from the *positive* rows
  neg_id_for_pos <- as.character(mc[[matched_negative_id_col]][pos_idx_all])
  ok <- !is.na(neg_id_for_pos) & nzchar(neg_id_for_pos)

  pos_idx <- pos_idx_all[ok]
  pos_ids <- ids[pos_idx]
  neg_ids <- neg_id_for_pos[ok]

  # Keep only pairs where the referenced negative exists in this GRanges
  neg_idx <- match(neg_ids, ids)
  ok2 <- !is.na(neg_idx)

  pos_idx <- pos_idx[ok2]
  pos_ids <- pos_ids[ok2]
  neg_ids <- neg_ids[ok2]
  neg_idx <- neg_idx[ok2]

  if (!length(pos_idx)) {
    out <- gr[0]
    if (return_diagnostics) {
      attr(out, "subset_pairs_diagnostics") <- list(
        n_total = length(gr),
        n_pos_total = length(pos_idx_all),
        n_pos_with_match = 0L,
        n_pairs = 0L,
        reason = "No positives with a valid matched_negative_id."
      )
    }
    return(out)
  }

  # Optional: drop duplicated negative assignments (should not happen, but prevents imbalance)
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

  # Optional: enforce reciprocal mapping (negative points back to the same positive)
  if (strict_reciprocal) {
    back <- as.character(mc[[matched_positive_id_col]][neg_idx])
    ok3 <- !is.na(back) & nzchar(back) & (back == pos_ids)

    pos_idx <- pos_idx[ok3]
    pos_ids <- pos_ids[ok3]
    neg_ids <- neg_ids[ok3]
    neg_idx <- neg_idx[ok3]
  }

  # Final balanced set: exactly one positive + one negative per remaining pair
  keep_idx <- unique(c(pos_idx, neg_idx))
  out <- gr[keep_idx]

  # (Optional) diagnostics
  if (return_diagnostics) {
    # recompute counts inside output
    mc2 <- S4Vectors::mcols(out)
    # try label
    y2 <- if (label_col %in% colnames(mc2)) as.integer(mc2[[label_col]]) else rep(NA_integer_, length(out))
    attr(out, "subset_pairs_diagnostics") <- list(
      n_total_in = length(gr),
      n_total_out = length(out),
      n_pos_total_in = length(pos_idx_all),
      n_pairs_out = length(pos_idx),
      expected_rows_out = 2L * length(pos_idx),
      n_pos_out = sum(y2 == 1L, na.rm = TRUE),
      n_neg_out = sum(y2 == 0L, na.rm = TRUE),
      strict_reciprocal = strict_reciprocal,
      drop_conflicts = drop_conflicts
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





