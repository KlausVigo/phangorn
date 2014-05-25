\name{designTree}
\alias{designTree}
\alias{designSplits}
\alias{nnls.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute a design matrix or non-negative LS }
\description{
\code{nnls.tree} estimates the branch length using non-negative least squares given a tree and a distance matrix. 
\code{designTree} and \code{designSplits} compute design matrices for the 
estimation of edge length of (phylogenetic) trees using linear models. 
For larger trees a sparse design matrix can save a lot of memory.
%\code{designTree} also computes a contrast matrix if the method is "rooted".
}
\usage{
designTree(tree, method = "unrooted", sparse=FALSE, ...)
designSplits(x, splits = "all", ...)
nnls.tree(dm, tree, rooted=FALSE, trace=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{an object of class \code{phylo} }
  \item{method}{design matrix for an "unrooted" or "rooted" ultrametric tree.}
  \item{sparse}{retyurn a sparse design matrix.}  
  \item{x}{number of taxa. }
  \item{splits}{one of "all", "star".}
  \item{dm}{a distance matrix.} 
  \item{rooted}{compute a "rooted" or "unrooted" tree.}
  \item{trace}{defines how much information is printed during optimisation.} 
  \item{\dots}{further arguments, passed to other methods.}
}
\value{
\code{nnls.tree} return a tree, i.e. an object of class \code{phylo}. 
\code{designTree} and \code{designSplits} a matrix, possibly sparse. 
}
\author{ Klaus Schliep \email{klaus.schliep@gmail.com} }

\seealso{ \code{\link[ape]{fastme}}, \code{\link[phangorn]{distanceHadamard}}, \code{\link[phangorn]{upgma}} }
\examples{
example(NJ)
dm <-  as.matrix(dm)
y <- dm[lower.tri(dm)]
X <- designTree(tree)
lm(y~X-1)
# avoids negative edge weights 
tree2 = nnls.tree(dm, tree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
