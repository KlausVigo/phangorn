\name{distanceHadamard}

\alias{distanceHadamard}
\title{Distance Hadamard}
\description{
Distance Hadamard produces spectra of splits from a distance matrix.
}
\usage{
distanceHadamard(dm, eps=0.001)
}
\arguments{
  \item{dm}{A distance matrix.}
  \item{eps}{Threshold value for splits.}
}

\value{
\code{distanceHadamard} returns a matrix. The first column 
contains the distance spectra, the second one the edge-spectra.
If eps is positive an object of with all splits greater eps is returned.
}
\references{Hendy, M. D. and Penny, D. (1993). Spectral Analysis of Phylogenetic
Data. \emph{Journal of Classification}, \bold{10}, 5-24.
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}, Tim White}
\seealso{
\code{\link{hadamard}}, \code{\link{lento}}, \code{\link{plot.networx}}, \code{\link{neighborNet}}
}

\examples{
data(yeast)
dm = dist.hamming(yeast)
dm = as.matrix(dm)
fit = distanceHadamard(dm)
lento(fit)
plot(as.networx(fit), "2D")
}
\keyword{ cluster }% at least one, from doc/KEYWORDS


