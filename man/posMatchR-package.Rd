\name{posMatchR-package}
\alias{posMatchR-package}
\alias{posMatchR}
\docType{package}
\title{posMatchR: annotation and matched backgrounds for single-nucleotide RNA sites}
\description{
posMatchR annotates single-nucleotide \code{GRanges} objects with gene,
transcript, region, metagene, feature-distance, and sequence-context variables,
and constructs matched background site sets for post-transcriptional site
analyses.
}
\details{
The usual workflow is \code{GRanges -> annotate_sites() -> plots/tables}, with
an optional matched-background branch using \code{combine_site_sets()} and
\code{match_background()}. \code{posmatchr_quickstart()} wraps these steps for
simple use cases.
}
\seealso{
\code{\link{annotate_sites}}, \code{\link{posmatchr_quickstart}},
\code{\link{match_background}}
}
\keyword{package}
