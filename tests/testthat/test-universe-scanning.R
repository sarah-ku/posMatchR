test_that("IUPAC motif scanning treats the subject as fixed", {
  skip_if_not_installed("Biostrings")

  hits <- posMatchR:::.scan_pattern_positions(
    sequence = "GAACC",
    patterns = "DRACH",
    fixed = "subject"
  )
  expect_equal(nrow(hits), 1L)
  expect_equal(hits$hit_start, 1L)

  fixed_subject_hits <- posMatchR:::.scan_pattern_positions(
    sequence = "GAACN",
    patterns = "DRACH",
    fixed = "subject"
  )
  expect_equal(nrow(fixed_subject_hits), 0L)

  ambiguous_subject_hits <- posMatchR:::.scan_pattern_positions(
    sequence = "GAACN",
    patterns = "DRACH",
    fixed = FALSE
  )
  expect_gt(nrow(ambiguous_subject_hits), 0L)
})

test_that("make_kmer_universe rejects non-concrete k-mers", {
  expect_error(
    make_kmer_universe(
      foreground = NULL,
      txdb = NULL,
      genome = NULL,
      kmers = "AANAA"
    ),
    "A/C/G/T"
  )

  expect_error(
    make_kmer_universe(
      foreground = NULL,
      txdb = NULL,
      genome = NULL,
      kmers = c("AAAAA", "AAAA")
    ),
    "same-width"
  )
})
