\name{ancestral.pml}
\alias{ancestral.pml}
\alias{ancestral.pars}
\alias{pace}
\alias{plotAnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ancestral character reconstruction. 
}
\description{
Marginal reconstruction of the ancestral character states.
}
\usage{
ancestral.pml(object, type = c("ml", "bayes"))
ancestral.pars(tree, data, type = c("MPR", "ACCTRAN"))
pace(tree, data, type = c("MPR", "ACCTRAN"))
plotAnc(tree, data, i, col=NULL, cex.pie=par("cex"), pos="bottomright", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class pml}
%  \item{type}{either "ml" or "bayes"}
  \item{tree}{a tree, i.e. an object of class pml}
  \item{data}{an object of class phyDat} 
  \item{type}{method used to assign characters to internal nodes, see details.}  
  \item{i}{plots the i-th character of the \code{data}.}
  \item{col}{a vector containing the colors for all possible states.}
  \item{cex.pie}{a numeric defining the size of the pie graphs}
  \item{pos}{a character string defining the positiond of the legend}
  \item{\dots}{Further arguments passed to or from other methods.}
}
%  \item{eps}{a small value to prevent rounding errors}
\details{
The argument "type" defines the criterion to assign the internal nodes. For \code{ancestral.pml} so far "ml" and
(empirical) "bayes" and for \code{ancestral.pars} "MPR" and "ACCTRAN" are possible. 

With parsimony reconstruction one has to keep in mind that there will be often no unique solution.

For further details see vignette("Ancestral").
}
\value{
%A matrix containing the the estimates character states. 
An object of class "phyDat", containing the ancestral states of all nodes. 
}
\references{
Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer Associates, Sunderland.

Swofford, D.L., Maddison, W.P. (1987) Reconstructing ancestral character states under Wagner parsimony. \emph{Math. Biosci.} \bold{87}: 199--229

Yang, Z. (2006). \emph{Computational Molecular evolution}. Oxford University Press, Oxford.
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{
\code{pml}, \code{parsimony}, \code{ace}, \code{root}
}
\examples{
example(NJ)
fit = pml(tree, Laurasiatherian)
anc.ml = ancestral.pml(fit, type = "ml")
anc.p = ancestral.pars(tree, Laurasiatherian)
\dontrun{
require(seqLogo)
seqLogo( t(subset(anc.ml, 48, 1:20)[[1]]), ic.scale=FALSE)
seqLogo( t(subset(anc.p, 48, 1:20)[[1]]), ic.scale=FALSE)
}
plotAnc(tree, anc.ml, 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
