% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_networx.R
\name{plot.networx}
\alias{plot.networx}
\title{plot phylogenetic networks}
\usage{
\method{plot}{networx}(x, type = "equal angle", use.edge.length = TRUE,
  show.tip.label = TRUE, show.edge.label = FALSE, edge.label = NULL,
  show.node.label = FALSE, node.label = NULL, show.nodes = FALSE,
  tip.color = "black", edge.color = "black", edge.width = 3,
  edge.lty = 1, split.color = NULL, split.width = NULL,
  split.lty = NULL, font = 3, cex = par("cex"), cex.node.label = cex,
  cex.edge.label = cex, col.node.label = tip.color,
  col.edge.label = tip.color, font.node.label = font,
  font.edge.label = font, underscore = FALSE, angle = 0, digits = 3,
  ...)
}
\arguments{
\item{x}{an object of class \code{"networx"}}

\item{type}{"3D" to plot using rgl or "equal angle" and "2D" in the normal
device.}

\item{use.edge.length}{a logical indicating whether to use the edge weights
of the network to draw the branches (the default) or not.}

\item{show.tip.label}{a logical indicating whether to show the tip labels on
the graph (defaults to \code{TRUE}, i.e. the labels are shown).}

\item{show.edge.label}{a logical indicating whether to show the tip labels
on the graph.}

\item{edge.label}{an additional vector of edge labels (normally not needed).}

\item{show.node.label}{a logical indicating whether to show the node labels
(see example).}

\item{node.label}{an additional vector of node labels (normally not needed).}

\item{show.nodes}{a logical indicating whether to show the nodes (see
example).}

\item{tip.color}{the colors used for the tip labels.}

\item{edge.color}{the colors used to draw edges.}

\item{edge.width}{the width used to draw edges.}

\item{edge.lty}{a vector of line types.}

\item{split.color}{the colors used to draw edges.}

\item{split.width}{the width used to draw edges.}

\item{split.lty}{a vector of line types.}

\item{font}{an integer specifying the type of font for the labels: 1 (plain
text), 2 (bold), 3 (italic, the default), or 4 (bold italic).}

\item{cex}{a numeric value giving the factor scaling of the labels.}

\item{cex.node.label}{a numeric value giving the factor scaling of the node
labels.}

\item{cex.edge.label}{a numeric value giving the factor scaling of the edge
labels.}

\item{col.node.label}{the colors used for the node labels.}

\item{col.edge.label}{the colors used for the edge labels.}

\item{font.node.label}{the font used for the node labels.}

\item{font.edge.label}{the font used for the edge labels.}

\item{underscore}{a logical specifying whether the underscores in tip labels
should be written as spaces (the default) or left as are (if TRUE).}

\item{angle}{rotate the plot.}

\item{digits}{if edge labels are numerical a positive integer indicating how
many significant digits are to be used.}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
\code{plot.networx} returns invisibly a list with paramters of the
plot.
}
\description{
So far not all parameters behave the same on the the \code{rgl} \code{"3D"}
and basic graphic \code{"2D"} device.
}
\details{
Often it is easier and safer to supply vectors of graphical parameters for
splits (e.g. splits.color) than for edges. These overwrite values edge.color.
}
\note{
The internal representation is likely to change.
}
\examples{

set.seed(1)
tree1 <- rtree(20, rooted=FALSE)
sp <- as.splits(rNNI(tree1, n=10))
net <- as.networx(sp)
plot(net)
plot(net, direction="axial")
\dontrun{
# also see example in consensusNet
example(consensusNet)
}
}
\references{
Dress, A.W.M. and Huson, D.H. (2004) Constructing Splits Graphs
\emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics
(TCBB)}, \bold{1(3)}, 109--115

Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W. (2017),
Intertwining phylogenetic trees and networks. \emph{Methods Ecol Evol}.
\bold{8}, 1212--1220. doi:10.1111/2041-210X.12760
}
\seealso{
\code{\link{consensusNet}}, \code{\link{neighborNet}},
\code{\link{splitsNetwork}}, \code{\link{hadamard}},
\code{\link{distanceHadamard}}, \code{\link{as.networx}},
\code{\link[ape]{evonet}}, \code{\link[ape]{as.phylo}},
\code{\link{densiTree}}, \code{\link[ape]{nodelabels}}
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\keyword{plot}
