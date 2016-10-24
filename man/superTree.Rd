\name{superTree}
\alias{superTree}
\alias{coalSpeciesTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Super Tree and Species Tree methods
}
\description{
These function \code{superTree} allows the estimation of a supertree from a set of trees using either Matrix representation parsimony, Robinson-Foulds or SPR as criterion.  \code{coalSpeciesTree} estimates species trees and can multiple individuals per species.}
\usage{
superTree(tree, method = "MRP", rooted=FALSE, trace=0, ...)
coalSpeciesTree(tree, X, sTree = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
an object of class \code{multiPhylo}
}
\item{method}{
An argument defining which algorithm is used to optimize the tree. 
Possible are "MRP", "NNI", and "SPR". 
}
\item{rooted}{
should the resulting supertrees be rooted.
}
\item{X}{
A \code{phyDat} object to define which individual belongs to which species. 
} 
\item{trace}{
defines how much information is printed during optimization.
}
\item{sTree}{
A species tree which forces the topology. 
} 
\item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
The function \code{superTree} extends the function mrp.supertree from Liam Revells, 
with artificial adding an outgroup on the root of the trees. 
This allows to root the supertree afterwards. The functions is internally used in DensiTree. The implementation for the RF- and SPR-supertree are very basic so far and assume that all trees share the same set of taxa.  

\code{coalSpeciesTree} estimates a single linkage tree as suggested by Liu et al. (2010) from the element wise minima of the cophenetic matrices of the gene trees. It extends \code{speciesTree} in ape as it allows that have several individuals per gene tree.  
}
\value{
The function returns an object of class \code{phylo}. 
}
\references{
Liu, L., Yu, L. and Pearl, D. K. (2010) Maximum tree: a consistent estimator of the species tree. \emph{Journal of Mathematical Biology}, \bold{60}, 95--106.


}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
Liam Revell 
Emmanuel Paradies
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{mrp.supertree},  \code{\link{speciesTree}}, \code{\link{densiTree}}, \code{\link{RF.dist}}, \code{\link{SPR.dist}}
}
\examples{
data(Laurasiatherian)
set.seed(1)
bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)upgma(dist.hamming(x)), bs=50)
class(bs) <- 'multiPhylo'

mrp_st <- superTree(bs, rooted=TRUE)
plot(superTree(mrp_st))
\dontrun{
rf_st <- superTree(bs, method = "RF")
spr_st <- superTree(bs, method = "SPR")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

