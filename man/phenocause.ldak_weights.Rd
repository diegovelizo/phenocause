% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_documentation.R
\docType{data}
\name{phenocause.ldak_weights}
\alias{phenocause.ldak_weights}
\title{SNP LD weights from the phenocause example dataset}
\format{
A list of three \code{data.frame} objects with two columns, corresponding to the *weights.short files generated by ldak \code{--calc-weights-all}:
\describe{
\item{\code{rsid}}{The SNP identifier.}
\item{\code{ldak_weight}}{The LDAK-computed LD weight.}
}
}
\usage{
phenocause.ldak_weights
}
\description{
A list of data frames containing LDAK-computed Linkage Disequilibrium (LD) weights for all SNPs in the example genetic dataset.
}
\details{
An LDAK weight for a given SNP measures the amount of local LD between that SNP and its neighbors. Within \code{phenocause}, these weights can be used to define sampling probabilities for causal sites, allowing for simulations where causal variants are enriched in either high- or low-LD regions.
}
\examples{
# Load the weights data
data(phenocause.ldak_weights)
# View chromosome names
names(phenocause.ldak_weights)
# View weights for chromosome 1
head(phenocause.ldak_weights[["chr1"]])

}
\seealso{
\code{\link[phenocause]{sample_causal_sites}}, \code{\link[phenocause]{metadata}}
}
\keyword{datasets}
