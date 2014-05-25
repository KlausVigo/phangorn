\name{pmlCluster}
\alias{pmlCluster}
%\alias{pmlCluster2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Stochastic Partitioning }
\description{
  Stochastic Partitioning of genes into p cluster. 
}
\usage{
pmlCluster(formula, fit, weight, p=1:5, part=NULL, nrep = 10, 
    control=pml.control(epsilon=1e-8, maxit=10, trace=1),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ a formula object (see details).}
  \item{fit}{ an object of class \code{pml}. }
  \item{weight}{ \code{weight}  is matrix of frequency of site patterns for all genes. }
  \item{p}{ number of clusters. }
  \item{part}{ starting partition, otherwise a random partition is generated. }
  \item{nrep}{ number of replicates for each p.  }
  \item{control}{A list of parameters for controlling the fitting process.} 
  \item{\dots}{ Further arguments passed to or from other methods. }
}
\details{
    The \code{formula} object allows to specify which parameter get     
    optimized. The formula is generally of the form \code{edge + bf + Q 
    ~ rate + shape + \dots}, on the left side are the parameters which 
    get optimized over all cluster, on the right the parameter which
    are optimized specific to each cluster. The parameters available 
    are \code{"nni", "bf", "Q", "inv", "shape", "edge", "rate"}.
    Each parameter can be used only once in the formula. 
    There are also some restriction on the combinations how parameters 
    can get used. \code{"rate"} is only available for the right side.  
    When \code{"rate"} is specified on the left hand side \code{"edge"}
    has to be specified (on either side), if  \code{"rate"} is specified on 
    the right hand side it follows directly that \code{edge} is too. 
}
\value{
  \code{pmlCluster} returns a list with elements
  \item{logLik}{log-likelihood of the fit}
  \item{trees}{a list of all trees during the optimization.} 
  \item{fits}{fits for the final partitions}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{  \code{\link{pml}},\code{\link{pmlPart}},\code{\link{pmlMix}},\code{\link{SH.test}} }
\examples{
\dontrun{
data(yeast)
dm <- dist.logDet(yeast)
tree <- NJ(dm)
fit=pml(tree,yeast)
fit = optim.pml(fit)

weight=xtabs(~ index+genes,attr(yeast, "index"))
set.seed(1)

sp <- pmlCluster(edge~rate, fit, weight, p=1:4)
sp
SH.test(sp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
