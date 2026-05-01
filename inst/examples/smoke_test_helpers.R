# Helper functions for local posMatchR smoke tests.
# These helpers are intended for local testing, not package unit tests.

load_granges_object <- function(path, object_name = NULL) {
  if (!file.exists(path)) stop("Input file does not exist: ", path)
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("rds", "rda", "rdata", "rdat")) {
    if (ext == "rds") {
      obj <- readRDS(path)
      if (!inherits(obj, "GRanges")) stop("readRDS(path) did not return a GRanges object.")
      return(obj)
    }

    env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = env)
    if (!is.null(object_name) && nzchar(object_name)) {
      if (!(object_name %in% loaded)) {
        stop("Object '", object_name, "' was not found in ", path,
             ". Available objects: ", paste(loaded, collapse = ", "))
      }
      obj <- get(object_name, envir = env)
    } else {
      gr_names <- loaded[vapply(loaded, function(nm) inherits(get(nm, envir = env), "GRanges"), logical(1))]
      if (length(gr_names) != 1L) {
        stop("Could not infer a unique GRanges object from ", path,
             ". Available GRanges objects: ", paste(gr_names, collapse = ", "),
             ". Pass the object name as the second command-line argument.")
      }
      obj <- get(gr_names[[1]], envir = env)
    }
    if (!inherits(obj, "GRanges")) stop("Selected object is not a GRanges.")
    return(obj)
  }

  stop("Unsupported input extension: ", ext, ". Use .rds, .RData, .Rda or .Rdat.")
}

summarise_granges_for_smoke <- function(gr, label = "input") {
  cat("\n--- ", label, " ---\n", sep = "")
  cat("n sites: ", length(gr), "\n", sep = "")
  cat("width summary:\n")
  print(summary(as.integer(GenomicRanges::width(gr))))
  cat("strand counts:\n")
  print(table(as.character(GenomicRanges::strand(gr)), useNA = "ifany"))
  cat("seqnames counts, first 20:\n")
  print(head(sort(table(as.character(GenomicRanges::seqnames(gr))), decreasing = TRUE), 20))
  invisible(NULL)
}

save_posmatchr_smoke_outputs <- function(ann, outdir, prefix) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(ann, file.path(outdir, paste0(prefix, "_annotated_granges.rds")))

  tab <- posMatchR::as_basic_site_table(ann)
  full_tab <- posMatchR::as_site_table(ann, columns = "all")
  utils::write.csv(full_tab, file.path(outdir, paste0(prefix, "_annotated_full_table.csv")), row.names = FALSE)
  utils::write.csv(tab, file.path(outdir, paste0(prefix, "_annotated_basic_table.csv")), row.names = FALSE)

  utils::write.csv(
    as.data.frame(table(location = S4Vectors::mcols(ann)$location, useNA = "ifany")),
    file.path(outdir, paste0(prefix, "_location_counts.csv")),
    row.names = FALSE
  )

  p <- posMatchR::plot_metagene_density(ann)
  ggplot2::ggsave(file.path(outdir, paste0(prefix, "_metagene_density.pdf")), p, width = 7, height = 4)

  if ("nearest_exon_junction_dist" %in% colnames(S4Vectors::mcols(ann))) {
    p2 <- posMatchR::plot_junction_distance_density(ann)
    ggplot2::ggsave(file.path(outdir, paste0(prefix, "_junction_distance_density.pdf")), p2, width = 7, height = 4)
  }

  cat("\nSaved smoke-test outputs to: ", normalizePath(outdir), "\n", sep = "")
  invisible(tab)
}
