\name{SH.test}
\alias{SH.test}
\title{Shimodaira-Hasegawa Test}
\usage{
SH.test(..., B = 10000, data=NULL)
}
\arguments{
  \item{...}{either a series of objects of class \code{"pml"}
    separated by commas, a list containing such objects
    or an object  of class \code{"pmlPart"}.}
  \item{B}{the number of bootstrap replicates.}
  \item{data}{an object of class  \code{"phyDat"}.}
}
\description{
  This function computes the Shimodaira--Hasegawa test for a set of
  trees.
}
\value{
  a numeric vector with the P-value associated with each tree given in
  \code{...}.
}
\references{
  Shimodaira, H. and Hasegawa, M. (1999) Multiple comparisons of
  log-likelihoods with applications to phylogenetic
  inference. \emph{Molecular Biology and Evolution}, \bold{16},
  1114--1116.
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{
  \code{\link{pml}}, \code{\link{pmlPart}}, \code{\link{pmlCluster}}
}
\examples{
data(Laurasiatherian)
dm <- dist.logDet(Laurasiatherian)
tree1 <- NJ(dm)
tree2 <- unroot(upgma(dm))
fit1 <- pml(tree1, Laurasiatherian)
fit2 <- pml(tree2, Laurasiatherian)
fit1 <- optim.pml(fit1) # optimize edge weights
fit2 <- optim.pml(fit2)
SH.test(fit1, fit2, B=500)
# in real analysis use larger B, e.g. 10000
\dontrun{
example(pmlPart)
SH.test(sp, B=1000)
}
}
\keyword{models}
