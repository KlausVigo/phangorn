% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distTree.R
\name{designTree}
\alias{designTree}
\alias{nnls.tree}
\alias{nnls.phylo}
\alias{nnls.splits}
\alias{nnls.networx}
\alias{designSplits}
\title{Compute a design matrix or non-negative LS}
\usage{
designTree(tree, method = "unrooted", sparse = FALSE, tip.dates = NULL,
  calibration = NULL, ...)

nnls.tree(dm, tree, method = c("unrooted", "ultrametric", "tipdated"),
  rooted = NULL, trace = 1, weight = NULL, balanced = FALSE,
  tip.dates = NULL)

nnls.phylo(x, dm, method = "unrooted", trace = 0, ...)

nnls.splits(x, dm, trace = 0, eps = 1e-08)

nnls.networx(x, dm, eps = 1e-08)

designSplits(x, splits = "all", ...)
}
\arguments{
\item{tree}{an object of class \code{phylo}}

\item{method}{compute an "unrooted", "ultrametric" or "tipdated" tree.}

\item{sparse}{return a sparse design matrix.}

\item{tip.dates}{a named vector of sampling times associated to the tips of
the tree.}

\item{calibration}{a named vector of calibration times associated to nodes of
the tree.}

\item{\dots}{further arguments, passed to other methods.}

\item{dm}{a distance matrix.}

\item{rooted}{compute a "ultrametric" or "unrooted" tree (better use method).}

\item{trace}{defines how much information is printed during optimization.}

\item{weight}{vector of weights to be used in the fitting process.
Weighted least squares is used with weights w, i.e., sum(w * e^2) is
minimized.}

\item{balanced}{use weights as in balanced fastME}

\item{x}{number of taxa.}

\item{eps}{minimum edge length (default s 1e-8).}

\item{splits}{one of "all", "star".}
}
\value{
\code{nnls.tree} return a tree, i.e. an object of class
\code{phylo}.  \code{designTree} and \code{designSplits} a matrix, possibly
sparse.
}
\description{
\code{nnls.tree} estimates the branch length using non-negative least
squares given a tree and a distance matrix.  \code{designTree} and
\code{designSplits} compute design matrices for the estimation of edge
length of (phylogenetic) trees using linear models.  For larger trees a
sparse design matrix can save a lot of memory. %\code{designTree} also
computes a contrast matrix if the method is "rooted".
}
\examples{

example(NJ)
dm <-  as.matrix(dm)
y <- dm[lower.tri(dm)]
X <- designTree(tree)
lm(y~X-1)
# avoids negative edge weights
tree2 <- nnls.tree(dm, tree)

}
\seealso{
\code{\link[ape]{fastme}}, \code{\link[ape]{rtt}},
\code{\link[phangorn]{distanceHadamard}},
\code{\link[phangorn]{splitsNetwork}}, \code{\link[phangorn]{upgma}}
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\keyword{cluster}
