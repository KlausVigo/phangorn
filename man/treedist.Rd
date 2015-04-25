\name{treedist}
\alias{treedist}
\alias{RF.dist}
%\alias{print.treedist}
\title{ Distances between trees }
\description{
  \code{treedist} computes different tree distance methods and \code{RF.dist} the Robinson-Foulds or symmetric distance.
}
\usage{
treedist(tree1, tree2, check.labels = TRUE)
RF.dist(tree1, tree2=NULL, check.labels=TRUE, rooted=FALSE)
}
\arguments{
  \item{tree1}{ A phylogenetic tree (class \code{phylo})
  or vector of trees (an object of class \code{multiPhylo}). See details }
  \item{tree2}{ A phylogenetic tree. }
  \item{check.labels}{compares labels of the trees.}
  \item{rooted}{take bipartitions for rooted trees into account, default is unrooting the trees.}
}
\value{
  \code{treedist} returns a vector containing the following tree distance methods
  \item{symmetric.difference}{symmetric.difference or  Robinson-Foulds distance}
  \item{branch.score.difference}{branch.score.difference}
  \item{path.difference}{path.difference}
  \item{weighted.path.difference}{weighted.path.difference}
}
\details{
The Robinson-Foulds distance is well defined only for bifurcating trees. 

RF.dist returns the Robinson-Foulds distance between either 2 trees or computes 
a matrix of all pairwise distances if a \code{multiPhylo} object is given. 
For large number of trees RF.dist can use a lot of memory!

% The function used internally is 2 * (nt - m) where nt is the number of tips and 
% m is the number of shared bipartitions. When there are multifurcations the 
% distance is therefore increasing!! This may be different to other implementations!

}
\references{Steel M. A. and Penny P. (1993) \emph{Distributions of tree comparison metrics - some new results}, Syst. Biol.,42(2), 126-141}
\author{ Klaus P. Schliep \email{klaus.schliep@gmail.com}} 
\examples{
tree1 <- rtree(100, rooted=FALSE)
tree2 <- rSPR(tree1, 3)
RF.dist(tree1, tree2)
treedist(tree1, tree2)
}
\keyword{ classif }% at least one, from doc/KEYWORDS
