\name{midpoint}
\alias{midpoint}
\alias{pruneTree}
\alias{getRoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tree manipulation}
\description{
\code{midpoint} performs midpoint rooting of a tree. 
\code{pruneTree} produces a consensus tree.
}
\usage{
midpoint(tree, node.labels = "support")
pruneTree(tree, ..., FUN = ">=")
getRoot(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{an object of class \code{phylo} }
  \item{FUN}{a function evaluated on the nodelabels, result must be logical.}
  \item{node.labels}{are nodel labels 'support' values, 'label' or  should be 'deleted'}
  \item{\dots}{further arguments, passed to other methods.  }
}
\details{
\code{pruneTree} prunes back a tree and produces a consensus tree, for trees already containing nodelabels. 
It assumes that nodelabels are numerical or character genereated from numericals, it uses as.numeric(as.character(tree$node.labels)) 
to convert them. 
\code{midpoint} so far does not transform node.labels properly.   
}
\value{
\code{pruneTree} and \code{midpoint} a tree. \code{getRoot} returns the root node.
}
\author{ Klaus Schliep \email{klaus.schliep@gmail.com} }

\seealso{ \code{\link[ape]{consensus}}, \code{\link[ape]{root}}, \code{\link[ape]{di2multi}} }
\examples{
tree = unroot(rtree(10))
tree$node.label = c("", round(runif(tree$Nnode-1), 3))

tree2 = midpoint(tree)
tree3 = pruneTree(tree, .5)

par(mfrow = c(3,1))
plot(tree, show.node.label=TRUE)
plot(tree2, show.node.label=TRUE)
plot(tree3, show.node.label=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }


