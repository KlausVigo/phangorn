\name{pmlPart}
\alias{pmlPart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Partition model. }
\description{
  Model to estimate phylogenies for partitioned data. 
}
\usage{
pmlPart(formula, object, control = pml.control(epsilon=1e-8, maxit=10, trace=1), 
    model=NULL, rooted=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ a formula object (see details).}
  \item{object}{ an object of class \code{pml} or a list of objects of class  \code{pml} . }
  \item{control}{A list of parameters for controlling the fitting process.} 
  \item{model}{A vector containing the models containing a model for each partition.} 
  \item{rooted}{Are the gene trees rooted (ultrametric) or unrooted.} 
  \item{\dots}{Further arguments passed to or from other methods. }
}
\details{
    The \code{formula} object allows to specify which parameter get     
    optimized. The formula is generally of the form \code{edge + bf + Q 
    ~ rate + shape + \dots}, on the left side are the parameters which 
    get optimized over all partitions, on the right the parameter which
    are optimized specific to each partition. The parameters available 
    are \code{"nni", "bf", "Q", "inv", "shape", "edge", "rate"}.
    Each parameters can be used only once in the formula. 
    \code{"rate"} and \code{"nni"} are only available for the right side
    of the formula.  
    
    For partitions with different edge weights, but same topology, \code{pmlPen}
    can try to find more parsimonious models (see example).             
}
\value{
  \code{kcluster} returns a list with elements
  \item{logLik}{log-likelihood of the fit}
  \item{trees}{a list of all trees during the optimization.} 
  \item{object}{an object of class \code{"pml"} or \code{"pmlPart"}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}

\seealso{  \code{\link{pml}},\code{\link{pmlCluster}},\code{\link{pmlMix}},\code{\link{SH.test}} }
\examples{
data(yeast)
dm <- dist.logDet(yeast)
tree <- NJ(dm)
fit <- pml(tree,yeast)
fits <- optim.pml(fit)

weight=xtabs(~ index+genes,attr(yeast, "index"))[,1:10]

sp <- pmlPart(edge ~ rate + inv, fits, weight=weight)
sp

\dontrun{
sp2 <- pmlPart(~ edge + inv, fits, weight=weight)
sp2
AIC(sp2)

sp3 <- pmlPen(sp2, lambda = 2) 
AIC(sp3)
}
}
\keyword{ cluster }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
