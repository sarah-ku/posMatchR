\name{posMatchR}
\alias{posMatchR}
\alias{posMatchR-package}
\alias{prepare_sites}
\alias{combine_site_sets}
\alias{makeID}
\alias{compute_tx_metrics}
\alias{build_tx_resources}
\alias{annotate_sites}
\alias{add_feature_geometry}
\alias{expand_sites_by_transcript}
\alias{add_metagene_split}
\alias{add_kmer}
\alias{as_site_table}
\alias{match_background}
\alias{subset_bg_set}
\alias{subset_matched_pairs}
\alias{assign_split_by_chromosome}
\alias{plot_metagene_density}
\alias{plot_distance_density}
\alias{plot_junction_distance_density}
\alias{plot_splice_distance_density}
\alias{plot_start_stop_distance_density}
\alias{plot_feature_width_bins}
\alias{plot_kmer_counts}
\alias{posmatchr_quickstart}
\title{posMatchR: annotation and matched backgrounds for single-nucleotide RNA sites}
\description{
posMatchR annotates single-nucleotide \code{GRanges} objects with gene,
transcript, region, metagene, feature-distance, and sequence-context variables,
and constructs matched background site sets.
}
\usage{
prepare_sites(gr, ...)
combine_site_sets(positives, background, ...)
compute_tx_metrics(txdb, ...)
build_tx_resources(txdb, ...)
annotate_sites(gr, txdb, ...)
add_feature_geometry(gr, resources)
expand_sites_by_transcript(gr, resources, ...)
add_metagene_split(gr, resources)
add_kmer(gr, genome, ...)
as_site_table(gr, ...)
match_background(gr, ...)
subset_bg_set(gr, ...)
subset_matched_pairs(gr, ...)
assign_split_by_chromosome(gr, ...)
plot_metagene_density(gr, ...)
plot_distance_density(gr, distance_col, ...)
plot_junction_distance_density(gr, ...)
plot_splice_distance_density(gr, ...)
plot_start_stop_distance_density(gr, ...)
plot_feature_width_bins(gr, ...)
plot_kmer_counts(gr, ...)
posmatchr_quickstart(sites, txdb, ...)
}
\details{
The usual workflow is \code{GRanges -> annotate_sites() -> plots/tables}, with
an optional matched background branch using \code{combine_site_sets()} and
\code{match_background()}. \code{posmatchr_quickstart()} wraps these steps for
vanilla use cases.
}
\examples{
if (requireNamespace("GenomicRanges", quietly = TRUE) &&
    requireNamespace("IRanges", quietly = TRUE)) {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100), strand = "+")
  gr <- prepare_sites(gr)
}
}
\keyword{package}
