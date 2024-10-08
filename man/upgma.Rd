% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/upgma.R
\name{upgma}
\alias{upgma}
\alias{wpgma}
\alias{supgma}
\title{UPGMA, WPGMA and sUPGMA}
\usage{
upgma(D, method = "average", ...)

wpgma(D, method = "mcquitty", ...)

supgma(D, tip.dates, trace = 0)
}
\arguments{
\item{D}{A distance matrix.}

\item{method}{The agglomeration method to be used. This should be (an
unambiguous abbreviation of) one of "ward", "single", "complete", "average",
"mcquitty", "median" or "centroid". The default is "average".}

\item{\dots}{Further arguments passed to or from other methods.}

\item{tip.dates}{A named vector of sampling times associated to the tips.}

\item{trace}{Show output during optimization (see details).}
}
\value{
A phylogenetic tree of class \code{phylo}.
}
\description{
UPGMA and WPGMA clustering. UPGMA and WPGMA are a wrapper function around
\code{\link[stats]{hclust}} returning a \code{phylo} object.
\code{supgma} perform serial sampled UPGMA similar to Drummond and Rodrigo
(2000).
}
\examples{

data(Laurasiatherian)
dm <- dist.ml(Laurasiatherian)
tree <- upgma(dm)
plot(tree)

}
\references{
Sneath, P. H., & Sokal, R. R. (1973). \emph{Numerical taxonomy.
The principles and practice of numerical classification.}

Sokal, R. R., & Michener, C. D. (1958). A statistical method for evaluating
systematic relationships. \emph{University of Kansas Scientific Bulletin},
v. 38.

Drummond, A., & Rodrigo, A. G. (2000). Reconstructing genealogies of serial
samples under the assumption of a molecular clock using serial-sample UPGMA.
\emph{Molecular Biology and Evolution}, \bold{17(12)}, 1807-1815.
}
\seealso{
\code{\link{hclust}}, \code{\link{dist.hamming}}, \code{\link{NJ}},
\code{\link[ape]{as.phylo}}, \code{\link[ape]{fastme}},
\code{\link{nnls.tree}}, \code{\link[ape]{rtt}}
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\keyword{cluster}
