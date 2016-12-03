\name{superTree}
\alias{superTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Super Tree methods
}
\description{
These function \code{superTree} allows the estimation of a supertree from a set of trees using either Matrix representation parsimony, Robinson-Foulds or SPR as criterion. }
\usage{
superTree(tree, method = "MRP", rooted=FALSE, trace=0, ...)
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
\item{trace}{
defines how much information is printed during optimization.
}
\item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
The function \code{superTree} extends the function mrp.supertree from Liam Revells, 
with artificial adding an outgroup on the root of the trees. 
This allows to root the supertree afterwards. The functions is internally used in DensiTree. The implementation for the RF- and SPR-supertree are very basic so far and assume that all trees share the same set of taxa.  
}
\value{
The function returns an object of class \code{phylo}. 
}
\references{
Baum, B. R., (1992) Combining trees as a way of combining data sets for phylogenetic inference, and the desirability of combining gene trees. \emph{Taxon}, \bold{41}, 3-10.

Ragan, M. A. (1992) Phylogenetic inference based on matrix representation of trees. \emph{Molecular Phylogenetics and Evolution}, \bold{1}, 53-58.
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
Liam Revell 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{mrp.supertree}, \code{\link{densiTree}}, \code{\link{RF.dist}}, \code{\link{SPR.dist}}
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

