\name{plot.networx}
\alias{plot.networx}
\alias{as.networx}
\alias{as.networx.splits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phylogenetic networks
}
\description{
A \code{networx} object hold the information for a phylogenetic network 
and are an generalisation of \code{phylo}. Therefore some generic function for \code{phylo} objects will also work for \code{networx} objects. 
\code{as.networx} convert \code{splits} objects into a  \code{networx} object. 
\code{plot.networx} plot phylogenetic network or split graphs. 
}
\usage{
as.networx(x, ...)
\method{as.networx}{splits}(x, ...)
\method{plot}{networx}(x, type="3D", show.tip.label=TRUE, show.edge.label=FALSE, 
    show.nodes=FALSE, tip.color="blue", edge.color="grey", edge.width=3, font=3, cex=1, 
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{"splits"} (as.networx) or \code{"networx"} (plot)
}
  \item{type}{
"3D" to plot using rgl or "2D" in the normal device. 
}
  \item{show.tip.label}{
a logical indicating whether to show the tip labels on the 
graph (defaults to \code{TRUE}, i.e. the labels are shown).
}
  \item{show.edge.label}{
a logical indicating whether to show the tip labels on the graph.
}
  \item{show.nodes}{
a logical indicating whether to show the nodes (see example).
}
  \item{tip.color}{
the colours used for the tip labels.
}
  \item{edge.color}{
the colours used to draw edges.
}
  \item{edge.width}{
the width used to draw edges.
}
  \item{font}{
an integer specifying the type of font for the labels: 1 (plain text), 
2 (bold), 3 (italic, the default), or 4 (bold italic).
}
  \item{cex}{
a numeric value giving the factor scaling of the labels.}
  \item{\dots}{
Further arguments passed to or from other methods.
}
}
\note{
The internal representation is likely to change. 
}
\references{
Dress, A.W.M. and Huson, D.H. (2004) Constructicting Splits Graphs \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB)}, \bold{1(3)}, 109--115
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\seealso{
\code{\link{consensusNet}}, \code{\link{hadamard}}, \code{\link{distanceHadamard}}, \code{\link{layout.kamada.kawai}}, \code{\link[ape]{evonet}}, \code{\link[ape]{as.igraph}}
}
\examples{
set.seed(1)
tree1 = rtree(20, rooted=FALSE)
sp = as.splits(rNNI(tree1, n=10))
net = as.networx(sp)
plot(net)
\dontrun{
# also see example in consensusNet 
example(consensusNet)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

