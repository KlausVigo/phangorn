\name{treedist}
\alias{treedist}
\alias{RF.dist}
\alias{KF.dist}
\alias{path.dist}
%\alias{print.treedist}
\title{ Distances between trees }
\description{
  \code{treedist} computes different tree distance methods and \code{RF.dist} the Robinson-Foulds or symmetric distance. 
}
\usage{
treedist(tree1, tree2, check.labels = TRUE)
RF.dist(tree1, tree2=NULL, check.labels=TRUE, rooted=FALSE)
KF.dist(tree1, tree2=NULL, check.labels=TRUE, rooted=FALSE)
path.dist(tree1, tree2=NULL, check.labels=TRUE, use.weight=FALSE)
}
\arguments{
  \item{tree1}{ A phylogenetic tree (class \code{phylo})
  or vector of trees (an object of class \code{multiPhylo}). See details }
  \item{tree2}{ A phylogenetic tree. }
  \item{check.labels}{compares labels of the trees.}
  \item{rooted}{take bipartitions for rooted trees into account, default is unrooting the trees.}
  \item{use.weight}{use edge.length argument or just count number of edges on the path (default)}
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
RF.dist returns the Robinson-Foulds distance (Robinson & Foulds 1981) between either 2 trees or computes a matrix of all pairwise distances if a \code{multiPhylo} object is given. 
\code{KF.dist} return the branch score distance (Kuhner & Felsenstein 1994). 
\code{path.dist} computes the path difference metric as described in Steel and Penny 1993).

For large number of trees the distance functions will use a lot of memory!


% The function used internally is 2 * (nt - m) where nt is the number of tips and 
% m is the number of shared bipartitions. When there are multifurcations the 
% distance is therefore increasing!! This may be different to other implementations!

}
\references{
Steel M. A. and Penny P. (1993) \emph{Distributions of tree comparison metrics - some new results}, Syst. Biol., \bold{42(2)}, 126--141

Kuhner, M. K. and Felsenstein, J. (1994) \emph{A simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates}, Molecular Biology and Evolution, \bold{11(3)}, 459--468

D.F. Robinson and L.R. Foulds (1981) \emph{Comparison of phylogenetic trees}, Mathematical Biosciences, \bold{53(1)}, 131--147
}
\author{ Klaus P. Schliep \email{klaus.schliep@gmail.com}} 
\seealso{\code{\link[ape]{dist.topo}}}
\examples{
tree1 <- rtree(100, rooted=FALSE)
tree2 <- rSPR(tree1, 3)
RF.dist(tree1, tree2)
treedist(tree1, tree2)
}
\keyword{ classif }% at least one, from doc/KEYWORDS
