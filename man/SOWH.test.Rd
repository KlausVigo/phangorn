\name{SOWH.test}
\alias{SOWH.test}
\title{Swofford-Olsen-Waddell-Hillis Test}
\usage{
SOWH.test(..., B = 10000, data=NULL)
}
\arguments{
  \item{...}{either a series of objects of class \code{"pml"}
    separated by commas, a list containing such objects
    or an object  of class \code{"pmlPart"}.}
  \item{B}{the number of bootstrap replicates.}
  \item{data}{an object of class  \code{"phyDat"}.}
}
\description{
  This function computes the Swofford--Olsen--Waddell--Hillis (SOWH) test, a parametric bootstrap test.
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
SOWH.test(fit, n=10)
}
}
\keyword{models}
