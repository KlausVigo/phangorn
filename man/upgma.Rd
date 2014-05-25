\name{upgma}
\alias{upgma}
\alias{wpgma}
\title{ UPGMA and WPGMA }
\description{
  UPGMA and WPGMA clustering. Just a wrapper function around \code{\link[stats]{hclust}}.
}
\usage{
upgma(D, method = "average", ...)
wpgma(D, method = "mcquitty", ...)
}
\arguments{
  \item{D}{A distance matrix.}
  \item{method}{The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". The default is "average".}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
A phylogenetic tree of class \code{phylo}. 
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}

\seealso{ \code{\link{hclust}}, \code{\link{dist.hamming}}, \code{\link{NJ}}, \code{\link{as.phylo}}, \code{\link{fastme}}, \code{\link{nnls.tree}} }
\examples{
data(Laurasiatherian)
dm = dist.ml(Laurasiatherian)
tree = upgma(dm)
plot(tree)
}
\keyword{cluster}
