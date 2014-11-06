\name{SOWH.test}
\alias{SOWH.test}
\title{Swofford-Olsen-Waddell-Hillis Test}
\usage{
SOWH.test(x, n = 100, restricted = list(optNni=FALSE), optNni=TRUE, trace = 1, ...)
}
\arguments{
  \item{x}{an object of class  \code{"pml"}.}
  \item{n}{the number of bootstrap replicates.}
  \item{restricted}{list of restricted parameter settings.}
  \item{optNni}{Logical value indicating whether toplogy gets optimized (NNI).}
  \item{trace}{Show output during computations.}
  \item{\dots}{Further arguments passed to \code{"optim.pml"}.}
}
\description{
This function computes the Swofford--Olsen--Waddell--Hillis (SOWH) test, a parametric bootstrap test. The function is computational very demanding and likely to be very slow. 
}
\details{
\code{SOWH.test} performs a parametric bootstrap test to compare two trees. 
It makes extensive use \code{simSeq} and \code{optim.pml} and can take quite long.  
}
\value{
  an object of class SOWH. That is a list with three elements, one is a matrix
  containing for each bootstrap replicate the (log-) likelihood of the restricted and   unrestricted estimate and two pml objetcs of the restricted and unrestricted model. 
}
\references{
Goldman, N., Anderson, J. P., and Rodrigo, A. G. (2000) Likelihood
-based tests of topologies in phylogenetics. \emph{Systematic Biology} \bold{49} 652-670.
  
Swofford, D.L., Olsen, G.J.,  Waddell, P.J. and Hillis, D.M. (1996) 
Phylogenetic Inference in Hillis, D.M., Moritz, C. and Mable, B.K. (Eds.)
\emph{Molecular Systematics} (2nd ed.) 407-514, Sunderland, MA: Sinauer
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{
  \code{\link{pml}}, \code{\link{pmlPart}}, \code{\link{pmlCluster}}, \code{\link{simSeq}}, \code{\link{SH.test}}
}
\examples{
# in real analysis use larger n, e.g. 500 preferably more
\dontrun{
data(Laurasiatherian)
dm <- dist.logDet(Laurasiatherian)
tree <- NJ(dm)
fit <- pml(tree, Laurasiatherian)
fit <- optim.pml(fit, TRUE)
set.seed(6)
tree <- rNNI(fit$tree, 1)
fit <- update(fit, tree = tree)
(res <- SOWH.test(fit, n=100))
summary(res)
}
}
\keyword{models}
