\name{allTrees}
\alias{allTrees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute all trees topologies.}
\description{
\code{allTrees} computes all tree topologies for rooted or unrooted trees with up to 10 tips. \code{allTrees}  returns bifurcating trees.
}
\usage{
allTrees(n, rooted = FALSE, tip.label = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of tips (<=10).}
  \item{rooted}{Rooted or unrooted trees (default: rooted).  }
  \item{tip.label}{Tip labels.}
}
\value{
   an object of class multiPhylo. 
}
\author{ Klaus Schliep \email{klaus.schliep@gmail.com} }

\seealso{\code{\link[ape]{rtree}}, \code{\link{nni}}}
\examples{
trees <- allTrees(5)
par(mfrow = c(3,5))
for(i in 1:15)plot(trees[[i]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
