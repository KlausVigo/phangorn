\name{SOWH.test}
\alias{SOWH.test}
\title{Swofford-Olsen-Waddell-Hillis Test}
\usage{
SOWH.test(x, n = 100, restricted = list(optNni=FALSE), optNni=TRUE, ...)
}
\arguments{
  \item{x}{an object of class  \code{"pml"}.}
  \item{n}{the number of bootstrap replicates.}
  \item{restricted}{list of restricted parameter settings.}
  \item{optNni}{Logical value indicating whether toplogy gets optimized (NNI).}
  \item{\dots}{Further arguments passed to \code{"optim.pml"}.}
}
\description{
  This function computes the Swofford--Olsen--Waddell--Hillis (SOWH) test, a parametric bootstrap test.
  The function is likely to be very slow. 
}
\value{
  a numeric vector with the P-value associated with each tree given in
  \code{...}.
}
\references{
  Goldman, N., Anderson, J. P., and Rodrigo, A. G. (2000) Likelihood
  -based tests of topologies in phylogenetics. \emph{Systematic Biology} \bold{49} 652-670.
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{
  \code{\link{pml}}, \code{\link{pmlPart}}, \code{\link{pmlCluster}}
}
\examples{
# in real analysis use larger n, e.g. 100
\dontrun{
data(Laurasiatherian)
dm <- dist.logDet(Laurasiatherian)
tree <- NJ(dm)
fit <- pml(tree, Laurasiatherian)
res <- SOWH.test(fit, n=10)
}
}
\keyword{models}
