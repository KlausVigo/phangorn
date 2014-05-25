\name{cladePar}
\alias{cladePar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility function to plot.phylo 
}
\description{
cladePar can help you coloring (choosing edge width/type) of clades. 
}
\usage{
cladePar(tree, node, edge.color = "red", tip.color = edge.color, edge.width = 1,
    edge.lty = 1, x = NULL, plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
an object of class phylo.
}
  \item{node}{
the node which is the common ancestor of the clade. 
}
  \item{edge.color}{
see plot.phylo.
}
  \item{tip.color}{
see plot.phylo.
}
  \item{edge.width}{
see plot.phylo.
}
  \item{edge.lty}{
see plot.phylo.
}
  \item{x}{
the result of a previous call to cladeInfo.
}
  \item{plot}{
logical, if TRUE the tree is plotted.   
}
  \item{\dots}{
Further arguments passed to or from other methods.
}
}
\value{
A list containing the information about the edges and tips. 
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\seealso{
\code{\link{plot.phylo}}
}
\examples{
tree = rtree(10)
plot(tree)
nodelabels()
x = cladePar(tree, 12)
cladePar(tree, 18, "blue", "blue", x=x, plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

