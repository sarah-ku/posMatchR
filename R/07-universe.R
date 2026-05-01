#' Select transcript names for background-universe construction
#'
#' @keywords internal
#' @noRd
.select_universe_transcripts <- function(foreground,
                                         resources,
                                         scope = c("genes", "transcripts", "all"),
                                         gene_id_col = "gene_id",
                                         tx_name_col = "tx_name") {
  scope <- match.arg(scope)
  tx_metrics <- resources$tx_metrics
  tx_key <- if (!is.null(resources$tx_key)) resources$tx_key else .tx_key_from_metrics(tx_metrics)

  if (scope == "all") return(unique(tx_key))
  if (is.null(foreground)) stop("`foreground` is required when scope is 'genes' or 'transcripts'.")
  if (!inherits(foreground, "GRanges")) stop("`foreground` must be a GRanges.")

  mc <- S4Vectors::mcols(foreground)

  if (scope == "genes" && gene_id_col %in% colnames(mc) && "gene_id" %in% colnames(tx_metrics)) {
    genes <- unique(as.character(mc[[gene_id_col]]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    tx <- tx_key[as.character(tx_metrics$gene_id) %in% genes]
    tx <- tx[!is.na(tx) & nzchar(tx)]
    if (length(tx)) return(unique(tx))
    warning("No transcripts found for foreground gene IDs; falling back to transcript-level scope.")
  }

  if (!(tx_name_col %in% colnames(mc))) {
    stop("Could not select transcripts: missing `", tx_name_col, "` in foreground metadata.")
  }
  tx <- unique(as.character(mc[[tx_name_col]]))
  tx <- tx[!is.na(tx) & nzchar(tx)]
  tx <- intersect(tx, tx_key)
  if (!length(tx)) stop("No foreground transcripts overlap the TxDb transcript names.")
  tx
}

.universe_region_list <- function(resources, region) {
  switch(
    region,
    fiveUTR = resources$five_by_tx,
    coding = resources$cds_by_tx,
    threeUTR = resources$three_by_tx,
    exon = resources$exons_by_tx,
    transcript = resources$exons_by_tx,
    stop("Unsupported region: ", region, ". Use fiveUTR, coding, threeUTR, exon, or transcript.")
  )
}

.iupac_to_plain <- function(x) {
  x <- toupper(as.character(x))
  chartr("U", "T", x)
}

.scan_one_sequence <- function(sequence, pattern, fixed = TRUE) {
  sequence <- .iupac_to_plain(sequence)
  pattern <- .iupac_to_plain(pattern)
  if (is.na(sequence) || is.na(pattern) || !nzchar(sequence) || !nzchar(pattern)) return(integer(0))
  hits <- Biostrings::matchPattern(pattern, Biostrings::DNAString(sequence), fixed = fixed)
  as.integer(IRanges::start(hits))
}


.can_use_pdict <- function(patterns, fixed = TRUE) {
  if (!isTRUE(fixed)) return(FALSE)
  patterns <- .iupac_to_plain(patterns)
  if (!length(patterns)) return(FALSE)
  if (!all(nchar(patterns) == nchar(patterns[1L]))) return(FALSE)
  all(grepl("^[ACGT]+$", patterns))
}

.make_pdict <- function(patterns) {
  if (!.can_use_pdict(patterns, fixed = TRUE)) return(NULL)
  tryCatch(
    Biostrings::PDict(Biostrings::DNAStringSet(patterns)),
    error = function(e) NULL
  )
}

.scan_pattern_positions <- function(sequence, patterns, fixed = TRUE, pdict = NULL) {
  sequence <- .iupac_to_plain(sequence)
  if (is.na(sequence) || !nzchar(sequence)) {
    return(data.frame(pattern_i = integer(0), hit_start = integer(0)))
  }

  if (!is.null(pdict)) {
    hits <- Biostrings::matchPDict(pdict, Biostrings::DNAString(sequence))
    starts <- vector("list", length(patterns))
    for (j in seq_along(patterns)) {
      hj <- hits[[j]]
      if (length(hj)) starts[[j]] <- as.integer(IRanges::start(hj))
    }
    lens <- vapply(starts, length, integer(1))
    idx <- which(lens > 0L)
    if (!length(idx)) return(data.frame(pattern_i = integer(0), hit_start = integer(0)))
    return(data.frame(
      pattern_i = rep(idx, lens[idx]),
      hit_start = as.integer(unlist(starts[idx], use.names = FALSE))
    ))
  }

  out <- vector("list", length(patterns))
  n_out <- 0L
  for (j in seq_along(patterns)) {
    hits <- .scan_one_sequence(sequence, patterns[j], fixed = fixed)
    if (!length(hits)) next
    n_out <- n_out + 1L
    out[[n_out]] <- data.frame(
      pattern_i = rep.int(j, length(hits)),
      hit_start = as.integer(hits)
    )
  }
  if (!n_out) return(data.frame(pattern_i = integer(0), hit_start = integer(0)))
  do.call(rbind, out[seq_len(n_out)])
}

.coord_key <- function(gr) {
  paste0(as.character(GenomicRanges::seqnames(gr)), ":", GenomicRanges::start(gr), ":", as.character(GenomicRanges::strand(gr)))
}

.scan_universe_patterns <- function(foreground,
                                    txdb,
                                    genome,
                                    patterns,
                                    fixed = TRUE,
                                    resources = NULL,
                                    scope = c("genes", "transcripts", "all"),
                                    regions = c("fiveUTR", "coding", "threeUTR"),
                                    gene_id_col = "gene_id",
                                    tx_name_col = "tx_name",
                                    site_offset = NULL,
                                    exclude_foreground = TRUE,
                                    site_id_prefix = "bg_",
                                    candidate_type = "sequence") {
  scope <- match.arg(scope)
  if (missing(genome) || is.null(genome)) stop("`genome` is required.")
  if (missing(txdb) || is.null(txdb)) stop("`txdb` is required.")
  if (!length(patterns)) stop("`patterns` must contain at least one sequence pattern.")

  patterns <- unique(.iupac_to_plain(patterns))
  patterns <- patterns[!is.na(patterns) & nzchar(patterns)]
  if (!length(patterns)) stop("No valid sequence patterns supplied.")

  if (is.null(resources)) resources <- build_tx_resources(txdb)
  tx_names <- .select_universe_transcripts(
    foreground = foreground,
    resources = resources,
    scope = scope,
    gene_id_col = gene_id_col,
    tx_name_col = tx_name_col
  )

  if (is.null(site_offset)) {
    site_offset <- ceiling(nchar(patterns) / 2)
  } else if (length(site_offset) == 1L) {
    site_offset <- rep.int(as.integer(site_offset), length(patterns))
  } else if (length(site_offset) != length(patterns)) {
    stop("`site_offset` must be NULL, length 1, or length(patterns).")
  }
  site_offset <- as.integer(site_offset)
  if (any(!is.finite(site_offset) | site_offset < 1L | site_offset > nchar(patterns))) {
    stop("Each `site_offset` must be between 1 and the corresponding pattern width.")
  }

  pdict <- .make_pdict(patterns)

  tx_metrics <- as.data.frame(resources$tx_metrics, stringsAsFactors = FALSE)
  tx_key <- if (!is.null(resources$tx_key)) resources$tx_key else .tx_key_from_metrics(tx_metrics)
  gene_by_tx <- if ("gene_id" %in% colnames(tx_metrics)) {
    stats::setNames(as.character(tx_metrics$gene_id), tx_key)
  } else {
    stats::setNames(rep(NA_character_, length(tx_key)), tx_key)
  }

  out <- list()
  out_i <- 0L

  for (region in regions) {
    gl <- .universe_region_list(resources, region)
    keep_tx <- intersect(names(gl), tx_names)
    if (!length(keep_tx)) next
    gl <- gl[keep_tx]

    for (tx in names(gl)) {
      segs <- gl[[tx]]
      if (!length(segs)) next
      seqs <- as.character(Biostrings::getSeq(genome, segs))

      if (region == "transcript") {
        tx_seq <- paste0(seqs, collapse = "")
        if (!nzchar(tx_seq)) next
        widths <- as.integer(GenomicRanges::width(segs))
        cum_end <- cumsum(widths)
        cum_start0 <- c(0L, cum_end[-length(cum_end)])

        hit_df <- .scan_pattern_positions(tx_seq, patterns, fixed = fixed, pdict = pdict)
        if (!nrow(hit_df)) next

        tx_pos <- hit_df$hit_start + site_offset[hit_df$pattern_i] - 1L
        ok <- tx_pos >= 1L & tx_pos <= cum_end[length(cum_end)]
        if (!any(ok)) next
        tx_pos <- as.integer(tx_pos[ok])
        pattern_i <- hit_df$pattern_i[ok]

        exon_idx <- findInterval(tx_pos - 1L, cum_end) + 1L
        ok_ex <- exon_idx >= 1L & exon_idx <= length(segs)
        if (!any(ok_ex)) next
        tx_pos <- tx_pos[ok_ex]
        pattern_i <- pattern_i[ok_ex]
        exon_idx <- exon_idx[ok_ex]
        offset_within <- tx_pos - cum_start0[exon_idx]
        ex <- segs[exon_idx]
        st <- as.character(GenomicRanges::strand(ex))
        st[is.na(st) | !nzchar(st)] <- "*"

        pos <- ifelse(
          st == "-",
          GenomicRanges::end(ex) - offset_within + 1L,
          GenomicRanges::start(ex) + offset_within - 1L
        )
        ok_pos <- is.finite(pos) & pos >= GenomicRanges::start(ex) & pos <= GenomicRanges::end(ex)
        if (!any(ok_pos)) next

        gr <- GenomicRanges::GRanges(
          seqnames = as.character(GenomicRanges::seqnames(ex[ok_pos])),
          ranges = IRanges::IRanges(start = as.integer(pos[ok_pos]), width = 1L),
          strand = st[ok_pos]
        )
        S4Vectors::mcols(gr)$candidate_pattern <- patterns[pattern_i[ok_pos]]
        S4Vectors::mcols(gr)$candidate_kmer <- patterns[pattern_i[ok_pos]]
        S4Vectors::mcols(gr)$candidate_type <- candidate_type
        S4Vectors::mcols(gr)$candidate_region <- region
        S4Vectors::mcols(gr)$candidate_tx_name <- tx
        S4Vectors::mcols(gr)$candidate_gene_id <- gene_by_tx[tx]
        out_i <- out_i + 1L
        out[[out_i]] <- gr
        next
      }

      for (i in seq_along(segs)) {
        seg <- segs[i]
        seq_i <- seqs[[i]]
        if (!nzchar(seq_i)) next
        st <- as.character(GenomicRanges::strand(seg))
        if (is.na(st) || !nzchar(st)) st <- "*"

        hit_df <- .scan_pattern_positions(seq_i, patterns, fixed = fixed, pdict = pdict)
        if (!nrow(hit_df)) next

        offset0 <- hit_df$hit_start + site_offset[hit_df$pattern_i] - 2L
        if (st == "-") {
          pos <- GenomicRanges::end(seg) - offset0
        } else {
          pos <- GenomicRanges::start(seg) + offset0
        }

        ok <- is.finite(pos) & pos >= GenomicRanges::start(seg) & pos <= GenomicRanges::end(seg)
        if (!any(ok)) next
        pos <- as.integer(pos[ok])
        pattern_i <- hit_df$pattern_i[ok]

        gr <- GenomicRanges::GRanges(
          seqnames = as.character(GenomicRanges::seqnames(seg)),
          ranges = IRanges::IRanges(start = pos, width = 1L),
          strand = st
        )
        S4Vectors::mcols(gr)$candidate_pattern <- patterns[pattern_i]
        S4Vectors::mcols(gr)$candidate_kmer <- patterns[pattern_i]
        S4Vectors::mcols(gr)$candidate_type <- candidate_type
        S4Vectors::mcols(gr)$candidate_region <- region
        S4Vectors::mcols(gr)$candidate_tx_name <- tx
        S4Vectors::mcols(gr)$candidate_gene_id <- gene_by_tx[tx]
        out_i <- out_i + 1L
        out[[out_i]] <- gr
      }
    }
  }

  if (!length(out)) {
    empty <- GenomicRanges::GRanges()
    S4Vectors::mcols(empty)$site_id <- character(0)
    return(empty)
  }

  gr <- do.call(c, out)
  gr <- prepare_sites(gr, site_id_col = "site_id", strip_mcols = FALSE)

  key <- .coord_key(gr)
  gr <- gr[!duplicated(key)]

  if (exclude_foreground && !is.null(foreground) && inherits(foreground, "GRanges")) {
    fg_key <- .coord_key(foreground)
    gr <- gr[!(.coord_key(gr) %in% fg_key)]
  }

  ids <- paste0(site_id_prefix, .coord_key(gr))
  ids <- make.unique(ids, sep = "_")
  S4Vectors::mcols(gr)$site_id <- ids
  names(gr) <- ids
  gr
}

#' Build an eligible single-base background universe
#'
#' Finds all occurrences of one or more transcript-oriented bases, usually
#' \code{"A"}, in the selected transcript space and returns them as point-site
#' \code{GRanges}. For a negative-strand transcript, bases are interpreted in the
#' transcript/RNA orientation, not the reference plus-strand orientation.
#'
#' @param foreground Annotated foreground \code{GRanges}; used to define genes or
#'   transcripts when \code{scope} is not \code{"all"}.
#' @param txdb A \code{TxDb} with seqlevels compatible with \code{genome}.
#' @param genome A BSgenome-like object accepted by \code{Biostrings::getSeq()}.
#' @param base One or more bases, e.g. \code{"A"}.
#' @param scope Use all transcripts in foreground genes, foreground transcripts,
#'   or all TxDb transcripts.
#' @param regions Transcript regions to scan.
#' @param resources Optional \code{build_tx_resources(txdb)} output.
#' @param gene_id_col Foreground gene ID column.
#' @param tx_name_col Foreground transcript-name column.
#' @param exclude_foreground If TRUE, remove coordinates already present in foreground.
#'
#' @return Candidate background sites as \code{GRanges}.
#' @export
make_base_universe <- function(foreground,
                               txdb,
                               genome,
                               base = "A",
                               scope = c("genes", "transcripts", "all"),
                               regions = c("fiveUTR", "coding", "threeUTR"),
                               resources = NULL,
                               gene_id_col = "gene_id",
                               tx_name_col = "tx_name",
                               exclude_foreground = TRUE) {
  .scan_universe_patterns(
    foreground = foreground,
    txdb = txdb,
    genome = genome,
    patterns = base,
    fixed = TRUE,
    resources = resources,
    scope = scope,
    regions = regions,
    gene_id_col = gene_id_col,
    tx_name_col = tx_name_col,
    site_offset = 1L,
    exclude_foreground = exclude_foreground,
    site_id_prefix = "base_bg_",
    candidate_type = "base"
  )
}

#' Build a foreground k-mer background universe
#'
#' Finds occurrences of observed k-mers, usually the site-centred k-mers already
#' added by \code{add_kmer()}, within the selected transcript space. Matching later
#' with \code{kmer_match = TRUE} gives an exact foreground/background k-mer balance
#' for successfully matched pairs.
#'
#' @param foreground Annotated foreground \code{GRanges} containing \code{kmer_col}.
#' @param txdb A \code{TxDb} with seqlevels compatible with \code{genome}.
#' @param genome A BSgenome-like object accepted by \code{Biostrings::getSeq()}.
#' @param kmers Optional explicit k-mer vector. If NULL, k-mers are taken from \code{foreground}.
#' @param kmer_col Metadata column containing observed k-mers.
#' @param min_count Drop foreground k-mers observed fewer than this many times.
#' @param site_offset One-based position in the k-mer to use as the point site. NULL uses the centre.
#' @inheritParams make_base_universe
#'
#' @return Candidate background sites as \code{GRanges}.
#' @export
make_kmer_universe <- function(foreground,
                               txdb,
                               genome,
                               kmers = NULL,
                               kmer_col = "kmer",
                               min_count = 1L,
                               site_offset = NULL,
                               scope = c("genes", "transcripts", "all"),
                               regions = c("fiveUTR", "coding", "threeUTR"),
                               resources = NULL,
                               gene_id_col = "gene_id",
                               tx_name_col = "tx_name",
                               exclude_foreground = TRUE) {
  if (is.null(kmers)) {
    if (is.null(foreground) || !inherits(foreground, "GRanges")) stop("`foreground` must be a GRanges when `kmers` is NULL.")
    mc <- S4Vectors::mcols(foreground)
    if (!(kmer_col %in% colnames(mc))) stop("Missing kmer column in foreground: ", kmer_col)
    kmers <- .iupac_to_plain(as.character(mc[[kmer_col]]))
    kmers <- kmers[!is.na(kmers) & nzchar(kmers)]
    if (!length(kmers)) stop("No non-missing k-mers found in foreground column: ", kmer_col)
    tab <- table(kmers)
    kmers <- names(tab)[tab >= as.integer(min_count)]
    if (!length(kmers)) stop("No k-mers remain after min_count filtering.")
  }

  .scan_universe_patterns(
    foreground = foreground,
    txdb = txdb,
    genome = genome,
    patterns = unique(kmers),
    fixed = TRUE,
    resources = resources,
    scope = scope,
    regions = regions,
    gene_id_col = gene_id_col,
    tx_name_col = tx_name_col,
    site_offset = site_offset,
    exclude_foreground = exclude_foreground,
    site_id_prefix = "kmer_bg_",
    candidate_type = "kmer"
  )
}

#' Build an IUPAC motif background universe
#'
#' Scans transcript features for IUPAC DNA/RNA patterns such as \code{"DRACH"}.
#' The point site is normally the centre of the motif, or a user-specified
#' \code{site_offset}.
#'
#' @param patterns One or more IUPAC patterns. \code{U} is treated as \code{T}.
#' @param site_offset One-based motif position to use as the point site. NULL uses the centre.
#' @inheritParams make_base_universe
#'
#' @return Candidate background sites as \code{GRanges}.
#' @export
make_motif_universe <- function(foreground,
                                txdb,
                                genome,
                                patterns,
                                site_offset = NULL,
                                scope = c("genes", "transcripts", "all"),
                                regions = c("fiveUTR", "coding", "threeUTR"),
                                resources = NULL,
                                gene_id_col = "gene_id",
                                tx_name_col = "tx_name",
                                exclude_foreground = TRUE) {
  .scan_universe_patterns(
    foreground = foreground,
    txdb = txdb,
    genome = genome,
    patterns = patterns,
    fixed = FALSE,
    resources = resources,
    scope = scope,
    regions = regions,
    gene_id_col = gene_id_col,
    tx_name_col = tx_name_col,
    site_offset = site_offset,
    exclude_foreground = exclude_foreground,
    site_id_prefix = "motif_bg_",
    candidate_type = "motif"
  )
}

#' Randomly match positives to background candidates within groups
#'
#' This provides a simple baseline matcher: each positive is matched to one
#' negative sampled from the same group, usually \code{gene_id} or \code{tx_name}.
#' Setting \code{kmer_match = TRUE} additionally requires the same k-mer, which
#' balances motif/k-mer frequencies among successfully matched pairs.
#'
#' @param gr Combined and annotated foreground/background \code{GRanges}.
#' @param label_col Binary label column, where 1 is foreground and 0 is background.
#' @param group_col Metadata column used to constrain matching, e.g. \code{"gene_id"}.
#' @param kmer_match If TRUE, match within \code{group_col} and \code{kmer_col}.
#' @param kmer_col K-mer column used when \code{kmer_match = TRUE}.
#' @param location_col Optional location column.
#' @param locations Optional location classes eligible for matching.
#' @param seed Random seed.
#' @param return_diagnostics If TRUE, attach diagnostics to metadata.
#'
#' @return A \code{GRanges} with the standard posMatchR match columns.
#' @export
match_random_background <- function(gr,
                                    label_col = "label",
                                    group_col = "gene_id",
                                    kmer_match = FALSE,
                                    kmer_col = "kmer",
                                    location_col = "location",
                                    locations = NULL,
                                    seed = 1L,
                                    return_diagnostics = TRUE) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  gr <- prepare_sites(gr, label = NULL, label_col = label_col, strip_mcols = FALSE)
  mc <- S4Vectors::mcols(gr)
  if (!(label_col %in% colnames(mc))) stop("Missing label_col: ", label_col)
  if (!(group_col %in% colnames(mc))) stop("Missing group_col: ", group_col)
  if (kmer_match && !(kmer_col %in% colnames(mc))) stop("kmer_match=TRUE but missing kmer_col: ", kmer_col)

  y <- .normalise_binary_label(mc[[label_col]], strict = TRUE)
  mc[[label_col]] <- y

  eligible <- rep(TRUE, length(gr))
  if (!is.null(locations)) {
    if (!(location_col %in% colnames(mc))) stop("Missing location_col: ", location_col)
    eligible <- eligible & as.character(mc[[location_col]]) %in% as.character(locations)
  }

  group <- as.character(mc[[group_col]])
  group[is.na(group) | !nzchar(group)] <- NA_character_
  if (kmer_match) {
    base_group <- group
    km <- .iupac_to_plain(as.character(mc[[kmer_col]]))
    km[is.na(km) | !nzchar(km)] <- NA_character_
    group <- paste0(base_group, "||", km)
    group[is.na(base_group) | is.na(km)] <- NA_character_
  }

  eligible <- eligible & !is.na(group)
  pos_idx <- which(y == 1L & eligible)
  neg_idx <- which(y == 0L & eligible)

  matched_neg_for_pos <- rep(NA_character_, length(gr))
  matched_pos_for_neg <- rep(NA_character_, length(gr))
  neg_used <- rep(FALSE, length(gr))
  meta_delta <- rep(NA_real_, length(gr))
  meta <- if ("metagene_split3" %in% colnames(mc)) suppressWarnings(as.numeric(mc$metagene_split3)) else rep(NA_real_, length(gr))

  set.seed(seed)
  keys <- intersect(unique(group[pos_idx]), unique(group[neg_idx]))
  keys <- keys[!is.na(keys) & nzchar(keys)]

  for (key in keys) {
    P <- pos_idx[group[pos_idx] == key]
    N <- neg_idx[group[neg_idx] == key]
    if (!length(P) || !length(N)) next
    P <- sample(P)
    N <- sample(N)
    n <- min(length(P), length(N))
    if (n <= 0L) next
    P <- P[seq_len(n)]
    N <- N[seq_len(n)]
    matched_neg_for_pos[P] <- names(gr)[N]
    matched_pos_for_neg[N] <- names(gr)[P]
    neg_used[N] <- TRUE
    ok_meta <- is.finite(meta[P]) & is.finite(meta[N])
    if (any(ok_meta)) meta_delta[P[ok_meta]] <- abs(meta[P[ok_meta]] - meta[N[ok_meta]])
  }

  mc$is_positive <- as.integer(y == 1L)
  mc$matched_negative_id <- matched_neg_for_pos
  mc$matched_positive_id <- matched_pos_for_neg
  mc$meta_delta <- meta_delta
  mc$is_matched_negative <- as.integer(neg_used)
  mc$is_matched_positive <- as.integer(!is.na(matched_neg_for_pos))
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

  if (return_diagnostics) {
    diag <- list(
      method = "random_within_group",
      group_col = group_col,
      kmer_match = kmer_match,
      n_positive = length(pos_idx),
      n_negative = length(neg_idx),
      n_matched_positive = sum(!is.na(matched_neg_for_pos)),
      n_matched_negative = sum(neg_used),
      unmatched_positive = length(pos_idx) - sum(!is.na(matched_neg_for_pos))
    )
    md <- S4Vectors::metadata(gr)
    md$match_diagnostics <- diag
    if (is.null(md$posMatchR)) md$posMatchR <- list()
    md$posMatchR$match_diagnostics <- diag
    S4Vectors::metadata(gr) <- md
  }

  gr
}
