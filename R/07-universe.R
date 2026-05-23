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

.scan_exact_kmer_positions <- function(sequence, patterns) {
  sequence <- .iupac_to_plain(sequence)
  patterns <- .iupac_to_plain(patterns)
  if (is.na(sequence) || !nzchar(sequence)) {
    return(data.frame(pattern_i = integer(0), hit_start = integer(0)))
  }
  if (!length(patterns)) {
    return(data.frame(pattern_i = integer(0), hit_start = integer(0)))
  }

  widths <- nchar(patterns)
  if (!all(widths == widths[1L]) || !all(grepl("^[ACGT]+$", patterns))) {
    stop("Internal error: .scan_exact_kmer_positions requires exact same-width A/C/G/T patterns.")
  }

  k <- as.integer(widths[1L])
  n <- nchar(sequence)
  if (!is.finite(n) || n < k) {
    return(data.frame(pattern_i = integer(0), hit_start = integer(0)))
  }

  starts <- seq_len(n - k + 1L)
  words <- substring(sequence, starts, starts + k - 1L)
  pattern_i <- match(words, patterns)
  ok <- !is.na(pattern_i)
  if (!any(ok)) {
    return(data.frame(pattern_i = integer(0), hit_start = integer(0)))
  }

  data.frame(
    pattern_i = as.integer(pattern_i[ok]),
    hit_start = as.integer(starts[ok])
  )
}

.scan_pattern_positions <- function(sequence, patterns, fixed = TRUE, pdict = NULL) {
  sequence <- .iupac_to_plain(sequence)
  patterns <- .iupac_to_plain(patterns)
  if (is.na(sequence) || !nzchar(sequence)) {
    return(data.frame(pattern_i = integer(0), hit_start = integer(0)))
  }

  # Fast exact-k-mer path: enumerate all width-k windows in the subject once
  # and use vectorised membership testing. This avoids looping over every
  # k-mer for every transcript segment when the dictionary contains many
  # concrete 5-mers.
  if (.can_use_pdict(patterns, fixed = fixed)) {
    return(.scan_exact_kmer_positions(sequence, patterns))
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

.scan_exact_kmer_segments <- function(segs,
                                      seqs,
                                      tx_for_seg,
                                      region,
                                      patterns,
                                      site_offset,
                                      gene_by_tx,
                                      candidate_type) {
  patterns <- .iupac_to_plain(patterns)
  if (!.can_use_pdict(patterns, fixed = TRUE)) {
    stop("Internal error: .scan_exact_kmer_segments requires exact same-width A/C/G/T patterns.")
  }

  widths <- nchar(seqs)
  keep <- !is.na(seqs) & nzchar(seqs) & is.finite(widths) & widths > 0L
  if (!any(keep)) return(GenomicRanges::GRanges())

  segs <- segs[keep]
  seqs <- .iupac_to_plain(seqs[keep])
  tx_for_seg <- tx_for_seg[keep]
  widths <- as.integer(widths[keep])

  k_width <- as.integer(nchar(patterns[1L]))
  sep_len <- k_width
  sep <- paste(rep("N", sep_len), collapse = "")

  seg_start <- integer(length(widths))
  seg_start[1L] <- 1L
  if (length(widths) > 1L) {
    seg_start[-1L] <- 1L + cumsum(widths[-length(widths)] + sep_len)
  }
  seg_end <- seg_start + widths - 1L

  combined <- paste0(seqs, collapse = sep)
  hit_df <- .scan_exact_kmer_positions(combined, patterns)
  if (!nrow(hit_df)) return(GenomicRanges::GRanges())

  hit_start <- as.integer(hit_df$hit_start)
  hit_end <- hit_start + k_width - 1L
  seg_idx <- findInterval(hit_start, seg_start)
  ok <- seg_idx >= 1L & seg_idx <= length(segs) &
    hit_start >= seg_start[seg_idx] & hit_end <= seg_end[seg_idx]
  if (!any(ok)) return(GenomicRanges::GRanges())

  hit_start <- hit_start[ok]
  seg_idx <- seg_idx[ok]
  pattern_i <- as.integer(hit_df$pattern_i[ok])
  tx <- tx_for_seg[seg_idx]
  seg_hit <- segs[seg_idx]

  offset0 <- hit_start - seg_start[seg_idx] + site_offset[pattern_i] - 1L
  st <- as.character(GenomicRanges::strand(seg_hit))
  st[is.na(st) | !nzchar(st)] <- "*"

  pos <- ifelse(
    st == "-",
    GenomicRanges::end(seg_hit) - offset0,
    GenomicRanges::start(seg_hit) + offset0
  )
  ok_pos <- is.finite(pos) &
    pos >= GenomicRanges::start(seg_hit) &
    pos <= GenomicRanges::end(seg_hit)
  if (!any(ok_pos)) return(GenomicRanges::GRanges())

  pattern_i <- pattern_i[ok_pos]
  tx <- tx[ok_pos]
  seg_hit <- seg_hit[ok_pos]
  st <- st[ok_pos]
  pos <- as.integer(pos[ok_pos])

  gr <- GenomicRanges::GRanges(
    seqnames = as.character(GenomicRanges::seqnames(seg_hit)),
    ranges = IRanges::IRanges(start = pos, width = 1L),
    strand = st
  )
  S4Vectors::mcols(gr)$candidate_pattern <- patterns[pattern_i]
  S4Vectors::mcols(gr)$candidate_kmer <- patterns[pattern_i]
  S4Vectors::mcols(gr)$candidate_type <- candidate_type
  S4Vectors::mcols(gr)$candidate_region <- region
  S4Vectors::mcols(gr)$candidate_tx_name <- tx
  S4Vectors::mcols(gr)$candidate_gene_id <- gene_by_tx[tx]
  gr
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

  pdict <- if (.can_use_pdict(patterns, fixed = fixed)) NULL else .make_pdict(patterns)

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

    if (region != "transcript") {
      seg_n <- S4Vectors::elementNROWS(gl)
      if (!sum(seg_n)) next
      tx_for_seg <- rep(names(gl), seg_n)
      segs_all <- unlist(gl, use.names = FALSE)
      seqs_all <- as.character(Biostrings::getSeq(genome, segs_all))

      if (.can_use_pdict(patterns, fixed = fixed)) {
        gr_region <- .scan_exact_kmer_segments(
          segs = segs_all,
          seqs = seqs_all,
          tx_for_seg = tx_for_seg,
          region = region,
          patterns = patterns,
          site_offset = site_offset,
          gene_by_tx = gene_by_tx,
          candidate_type = candidate_type
        )
        if (length(gr_region)) {
          out_i <- out_i + 1L
          out[[out_i]] <- gr_region
        }
        next
      }

      for (i in seq_along(segs_all)) {
        tx <- tx_for_seg[[i]]
        seg <- segs_all[i]
        seq_i <- seqs_all[[i]]
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
      next
    }

    for (tx in names(gl)) {
      segs <- gl[[tx]]
      if (!length(segs)) next
      seqs <- as.character(Biostrings::getSeq(genome, segs))

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
#' @examples
#' if (interactive()) {
#'     txdb <- get(
#'         "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'         envir = asNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")
#'     )
#'     genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'     resources <- build_tx_resources(txdb)
#'     foreground <- GenomicRanges::GRanges(
#'         "chr1", IRanges::IRanges(11874, width = 1), strand = "+"
#'     )
#'     make_base_universe(foreground, txdb, genome, resources = resources)
#' }
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
#' @examples
#' if (interactive()) {
#'     txdb <- get(
#'         "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'         envir = asNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")
#'     )
#'     genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'     resources <- build_tx_resources(txdb)
#'     foreground <- GenomicRanges::GRanges(
#'         "chr1", IRanges::IRanges(11874, width = 1), strand = "+"
#'     )
#'     S4Vectors::mcols(foreground)$kmer <- "GGACT"
#'     make_kmer_universe(foreground, txdb, genome, resources = resources)
#' }
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
#' @examples
#' if (interactive()) {
#'     txdb <- get(
#'         "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'         envir = asNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")
#'     )
#'     genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'     resources <- build_tx_resources(txdb)
#'     foreground <- GenomicRanges::GRanges(
#'         "chr1", IRanges::IRanges(11874, width = 1), strand = "+"
#'     )
#'     make_motif_universe(foreground, txdb, genome, "DRACH", resources)
#' }
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
#' @param match_location If TRUE, additionally require positives and negatives to
#'   have the same \code{location_col}. This is useful when the simple random
#'   baseline should preserve broad transcript region, e.g. 5'UTR/CDS/3'UTR.
#' @param location_col Optional location column.
#' @param locations Optional location classes eligible for matching.
#' @param match_cols Additional metadata columns to append to the matching key,
#'   for example \code{"region_class"} or \code{"type"}.
#' @param seed Random seed.
#' @param return_diagnostics If TRUE, attach diagnostics to metadata.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'     "chr1", IRanges::IRanges(seq(10, 60, by = 10), width = 1), strand = "+"
#' )
#' gr <- prepare_sites(gr, label = c(1, 1, 1, 0, 0, 0))
#' S4Vectors::mcols(gr)$gene_id <- rep("gene1", length(gr))
#' S4Vectors::mcols(gr)$kmer <- rep("GGACT", length(gr))
#' match_random_background(gr, group_col = "gene_id", kmer_match = TRUE)
#'
#' @return A \code{GRanges} with the standard posMatchR match columns.
#' @export
match_random_background <- function(gr,
                                    label_col = "label",
                                    group_col = "gene_id",
                                    kmer_match = FALSE,
                                    kmer_col = "kmer",
                                    match_location = FALSE,
                                    location_col = "location",
                                    locations = NULL,
                                    match_cols = NULL,
                                    seed = 1L,
                                    return_diagnostics = TRUE) {
  if (!inherits(gr, "GRanges")) stop("`gr` must be a GRanges.")
  gr <- prepare_sites(gr, label = NULL, label_col = label_col, strip_mcols = FALSE)
  mc <- S4Vectors::mcols(gr)
  ids <- names(gr)

  if (!(label_col %in% colnames(mc))) stop("Missing label_col: ", label_col)
  if (!(group_col %in% colnames(mc))) stop("Missing group_col: ", group_col)
  if (kmer_match && !(kmer_col %in% colnames(mc))) stop("kmer_match=TRUE but missing kmer_col: ", kmer_col)
  if (isTRUE(match_location) && !(location_col %in% colnames(mc))) stop("match_location=TRUE but missing location_col: ", location_col)
  if (!is.null(locations) && !(location_col %in% colnames(mc))) stop("`locations` supplied but missing location_col: ", location_col)
  if (!is.null(match_cols)) {
    missing_cols <- setdiff(as.character(match_cols), colnames(mc))
    if (length(missing_cols)) stop("Missing match_cols: ", paste(missing_cols, collapse = ", "))
  }

  y <- .normalise_binary_label(mc[[label_col]], strict = TRUE)
  mc[[label_col]] <- y

  eligible <- rep(TRUE, length(gr))
  if (!is.null(locations)) {
    eligible <- eligible & as.character(mc[[location_col]]) %in% as.character(locations)
  }

  group <- as.character(mc[[group_col]])
  group[is.na(group) | !nzchar(group)] <- NA_character_

  append_to_group <- function(group, value) {
    value <- as.character(value)
    value[is.na(value) | !nzchar(value)] <- NA_character_
    out <- paste0(group, "||", value)
    out[is.na(group) | is.na(value)] <- NA_character_
    out
  }

  key_cols_used <- group_col
  if (isTRUE(match_location)) {
    group <- append_to_group(group, mc[[location_col]])
    key_cols_used <- c(key_cols_used, location_col)
  }

  if (!is.null(match_cols) && length(match_cols)) {
    for (cl in as.character(match_cols)) {
      group <- append_to_group(group, mc[[cl]])
    }
    key_cols_used <- c(key_cols_used, as.character(match_cols))
  }

  if (kmer_match) {
    km <- .iupac_to_plain(as.character(mc[[kmer_col]]))
    group <- append_to_group(group, km)
    key_cols_used <- c(key_cols_used, kmer_col)
  }

  eligible <- eligible & !is.na(group)
  pos_idx <- which(y == 1L & eligible)
  neg_idx <- which(y == 0L & eligible)

  matched_neg_for_pos <- rep(NA_character_, length(gr))
  matched_pos_for_neg <- rep(NA_character_, length(gr))
  meta_delta <- rep(NA_real_, length(gr))
  meta <- if ("metagene_split3" %in% colnames(mc)) {
    suppressWarnings(as.numeric(mc$metagene_split3))
  } else if ("metagene_prop" %in% colnames(mc)) {
    suppressWarnings(as.numeric(mc$metagene_prop))
  } else {
    rep(NA_real_, length(gr))
  }

  pair_pos <- integer(0)
  pair_neg <- integer(0)

  keys <- intersect(unique(group[pos_idx]), unique(group[neg_idx]))
  keys <- keys[!is.na(keys) & nzchar(keys)]

  withr::with_seed(seed, {
    for (key in keys) {
      P <- pos_idx[group[pos_idx] == key]
      N <- neg_idx[group[neg_idx] == key]
      if (!length(P) || !length(N)) next

      P <- .posmatchr_sample(P)
      N <- .posmatchr_sample(N)
      n <- min(length(P), length(N))
      if (n <= 0L) next

      pair_pos <- c(pair_pos, P[seq_len(n)])
      pair_neg <- c(pair_neg, N[seq_len(n)])
    }
  })

  # Rebuild all reciprocal pair columns from the selected pairs only. This is
  # deliberately done in one block so stale matching columns from a previous call
  # cannot leak into diagnostics or downstream matched-set subsetting.
  if (length(pair_pos)) {
    ord <- order(pair_pos, pair_neg)
    pair_pos <- pair_pos[ord]
    pair_neg <- pair_neg[ord]

    # Defensive conflict resolution. The loop above samples without replacement
    # within each key, but this protects against malformed inputs with duplicated
    # names or keys after user-side modification.
    keep <- !duplicated(pair_pos) & !duplicated(pair_neg)
    pair_pos <- pair_pos[keep]
    pair_neg <- pair_neg[keep]

    matched_neg_for_pos[pair_pos] <- ids[pair_neg]
    matched_pos_for_neg[pair_neg] <- ids[pair_pos]
    ok_meta <- is.finite(meta[pair_pos]) & is.finite(meta[pair_neg])
    if (any(ok_meta)) meta_delta[pair_pos[ok_meta]] <- abs(meta[pair_pos[ok_meta]] - meta[pair_neg[ok_meta]])
  }

  # Final reciprocal validation used for both flags and diagnostics.
  matched_pos_idx <- which(y == 1L & !is.na(matched_neg_for_pos) & nzchar(matched_neg_for_pos))
  matched_neg_idx <- match(matched_neg_for_pos[matched_pos_idx], ids)
  ok <- !is.na(matched_neg_idx) & y[matched_neg_idx] == 0L
  if (any(ok)) {
    back <- matched_pos_for_neg[matched_neg_idx[ok]]
    ok2 <- !is.na(back) & nzchar(back) & back == ids[matched_pos_idx[ok]]
    ok[ok] <- ok2
  }
  matched_pos_idx <- matched_pos_idx[ok]
  matched_neg_idx <- matched_neg_idx[ok]

  # If anything failed the reciprocal check, clear it rather than letting the
  # subsetter or plotting helpers see a partial/non-reciprocal pair.
  matched_neg_for_pos2 <- rep(NA_character_, length(gr))
  matched_pos_for_neg2 <- rep(NA_character_, length(gr))
  meta_delta2 <- rep(NA_real_, length(gr))
  if (length(matched_pos_idx)) {
    matched_neg_for_pos2[matched_pos_idx] <- ids[matched_neg_idx]
    matched_pos_for_neg2[matched_neg_idx] <- ids[matched_pos_idx]
    meta_delta2[matched_pos_idx] <- meta_delta[matched_pos_idx]
  }
  matched_neg_for_pos <- matched_neg_for_pos2
  matched_pos_for_neg <- matched_pos_for_neg2
  meta_delta <- meta_delta2

  kmer_mismatch <- NA_integer_
  if (kmer_match && length(matched_pos_idx)) {
    km_final <- .iupac_to_plain(as.character(mc[[kmer_col]]))
    kmer_mismatch <- sum(km_final[matched_pos_idx] != km_final[matched_neg_idx], na.rm = TRUE)
    if (kmer_mismatch > 0L) {
      warning("Internal k-mer matching validation found ", kmer_mismatch,
              " mismatched pairs. Clearing those pairs.")
      keep_km <- km_final[matched_pos_idx] == km_final[matched_neg_idx]
      keep_km[is.na(keep_km)] <- FALSE
      matched_neg_for_pos[matched_pos_idx[!keep_km]] <- NA_character_
      matched_pos_for_neg[matched_neg_idx[!keep_km]] <- NA_character_
      meta_delta[matched_pos_idx[!keep_km]] <- NA_real_
      matched_pos_idx <- matched_pos_idx[keep_km]
      matched_neg_idx <- matched_neg_idx[keep_km]
      kmer_mismatch <- 0L
    }
  }

  mc$is_positive <- as.integer(y == 1L)
  mc$matched_negative_id <- matched_neg_for_pos
  mc$matched_positive_id <- matched_pos_for_neg
  mc$meta_delta <- meta_delta
  mc$is_matched_positive <- as.integer(seq_along(gr) %in% matched_pos_idx)
  mc$is_matched_negative <- as.integer(seq_along(gr) %in% matched_neg_idx)
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
      key_cols_used = key_cols_used,
      kmer_match = isTRUE(kmer_match),
      kmer_col = if (isTRUE(kmer_match)) kmer_col else NULL,
      kmer_mismatch_in_validated_pairs = kmer_mismatch,
      match_location = isTRUE(match_location),
      location_col = location_col,
      locations = locations,
      match_cols = match_cols,
      n_positive_total = sum(y == 1L, na.rm = TRUE),
      n_negative_total = sum(y == 0L, na.rm = TRUE),
      n_positive_eligible = length(pos_idx),
      n_negative_eligible = length(neg_idx),
      n_matched_pairs = length(matched_pos_idx),
      n_matched_positive = length(matched_pos_idx),
      n_matched_negative = length(matched_neg_idx),
      unmatched_positive_eligible = length(pos_idx) - length(matched_pos_idx),
      unmatched_positive_total = sum(y == 1L, na.rm = TRUE) - length(matched_pos_idx)
    )
    md <- S4Vectors::metadata(gr)
    md$match_diagnostics <- diag
    if (is.null(md$posMatchR)) md$posMatchR <- list()
    md$posMatchR$match_diagnostics <- diag
    S4Vectors::metadata(gr) <- md
  }

  gr
}
