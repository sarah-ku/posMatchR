# posMatchR 0.99.0

* Initial Bioconductor submission version.
* Added transcript-context annotation, metagene summaries, sequence-context extraction, matched-background construction, and motif/point-feature enrichment plotting for single-nucleotide RNA site maps.
* Added a human GLORI quick-start vignette using bundled example sites and precomputed demonstration transcript resources.
* Simplified transcript-universe sequence searching: `make_kmer_universe()` now uses the exact same-width k-mer fast path and requires concrete A/C/G/T k-mers after U-to-T conversion, while `make_motif_universe()` uses `Biostrings::matchPattern(..., fixed = "subject")` for IUPAC motifs such as DRACH.
* Fixed the `make_motif_universe()` example to pass `resources` by name.
