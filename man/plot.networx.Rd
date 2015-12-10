\name{plot.networx}
\alias{plot.networx}
\alias{as.networx}
\alias{as.networx.splits}
%\alias{reorder.networx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phylogenetic networks
}
\description{
\code{as.networx} convert \code{splits} objects into a  \code{networx} object. 
\code{plot.networx} plot phylogenetic network or split graphs. 
}
\usage{
as.networx(x, ...)
\method{as.networx}{splits}(x, planar = FALSE, ...)
\method{plot}{networx}(x, type="3D", use.edge.length = TRUE, show.tip.label=TRUE, 
     show.edge.label=FALSE, edge.label = NULL, show.node.label=FALSE, 
     node.label = NULL, show.nodes=FALSE, tip.color="black", 
     edge.color="black", edge.width=3, edge.lty=1, font=3, 
     cex=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
an object of class \code{"splits"} (as.networx) or \code{"networx"} (plot)
}
\item{planar}{
logical whether to produce a planar graph from only cyclic splits 
(may excludes splits). 
}
\item{type}{
"3D" to plot using rgl or "2D" in the normal device. 
}
\item{use.edge.length}{
a logical indicating whether to use the edge weights of the network 
to draw the branches (the default) or not. 
}
\item{show.tip.label}{
a logical indicating whether to show the tip labels on the 
graph (defaults to \code{TRUE}, i.e. the labels are shown).
}
\item{show.edge.label}{
a logical indicating whether to show the tip labels on the graph.
}
\item{edge.label}{
an additional vector of edge labels (normally not needed).
}
\item{show.node.label}{
a logical indicating whether to show the node labels (see example).
}
\item{node.label}{
an additional vector of node labels (normally not needed).
}
\item{show.nodes}{
a logical indicating whether to show the nodes (see example).
}
\item{tip.color}{
the colors used for the tip labels.
}
\item{edge.color}{
the colors used to draw edges.
}
\item{edge.width}{
the width used to draw edges.
}
\item{edge.lty}{
a vector of line types.
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
\details{
A \code{networx} object hold the information for a phylogenetic network 
and extends the \code{phylo} object. Therefore some generic function for 
\code{phylo} objects will also work for \code{networx} objects. 
The argument planar = FALSE will create a planar split graph based on a 
cyclic ordering. These objects can be nicely plotted in "2D".
So far not all parameters behave the same on the the rgl "3D"
and basic graphic "2D" device. 
}
\note{
The internal representation is likely to change. 
}
\references{
Dress, A.W.M. and Huson, D.H. (2004) Constructing Splits Graphs \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB)}, \bold{1(3)}, 109--115
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\seealso{
\code{\link{consensusNet}}, \code{\link{neighborNet}}, \code{\link{splitsNetwork}}, \code{\link{hadamard}}, 
\code{\link{distanceHadamard}}, \code{\link{layout_with_kk}}, 
\code{\link[ape]{evonet}}, \code{\link[ape]{as.igraph}}, \code{\link{densiTree}}
}
\examples{
set.seed(1)
tree1 = rtree(20, rooted=FALSE)
sp = as.splits(rNNI(tree1, n=10))
net = as.networx(sp)
plot(net, "2D")
\dontrun{
# also see example in consensusNet 
example(consensusNet)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

