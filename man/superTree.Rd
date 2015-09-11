\name{superTree}
\alias{superTree}
\alias{coalSpeciesTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Super Tree and Species Tree methods
}
\description{
These function \code{superTree} allows the estimation of a rooted supertree from a set of rooted trees using Matrix representation parsimony.  \code{coalSpeciesTree} estimates species trees and can multiple individuals per species.}
\usage{
superTree(tree, method = "pratchet", rooted=TRUE, ...)
coalSpeciesTree(tree, X, sTree = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
an object of class \code{multiPhylo}
}
\item{method}{
An argument defining which algorithm is used to optimize the tree.
}
\item{rooted}{
should the resulting supertrees be rooted.
}
\item{X}{
A \code{phyDat} object to define which individual belongs to which species. 
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
This allows to root the supertree afterwards. The functions is internally used in DensiTree.

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
\code{mrp.supertree},  \code{\link{speciesTree}}, \code{\link{densiTree}}, \code{\link{hclust}}
}
\examples{
data(Laurasiatherian)
set.seed(1)
bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)upgma(dist.hamming(x)), bs=50)
class(bs) <- 'multiPhylo'
plot(superTree(bs))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

