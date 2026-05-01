# Motif enrichment profile around matched foreground/background sites.
#
# Assumes you already have:
#   matched_covariate_sets  GRanges from subset_matched_sets(match_background(...))
#   genome                  compatible BSgenome object
#
# Character motifs are exact/IUPAC patterns. RNA U is accepted and converted to DNA T.

library(posMatchR)

# Example exact/IUPAC-style RBP motifs.
# PUM1/PUM2 PRE consensus: UGUANAUA -> TGTANATA on DNA alphabet.
# AU-rich/HuR-like core: AUUUA -> ATTTA on DNA alphabet.
plot_motif_enrichment(
  matched_covariate_sets,
  genome = genome,
  motif = c(PUM_PRE = "TGTANATA", ARE_AUUUA = "ATTTA"),
  window = 250,
  set_col = "match_set",
  hit_position = "center",
  bin_size = 1,
  smooth_window = 11
)

# The same idea with universalmotif objects. This is most useful when you want
# PWM/log-odds motif scanning rather than exact/IUPAC matching.
if (requireNamespace("universalmotif", quietly = TRUE)) {
  pum <- universalmotif::create_motif(
    "TGTANATA",
    alphabet = "DNA",
    name = "PUM_PRE"
  )

  plot_motif_enrichment(
    matched_covariate_sets,
    genome = genome,
    motif = pum,
    method = "universalmotif",
    motif_name = "PUM_PRE",
    window = 250,
    set_col = "match_set",
    threshold = 0.8,
    threshold_type = "logodds",
    nthreads = 1L,
    hit_position = "center",
    smooth_window = 11
  )
}
