\name{plot.networx}
\alias{plot.networx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot phylogenetic networks
}
\description{
\code{plot.networx} plot phylogenetic network or split graphs. 
}
\usage{
\method{plot}{networx}(x, type="3D", show.tip.label=TRUE, show.edge.label=FALSE, 
    show.nodes=FALSE, tip.color="blue", edge.color="grey", edge.width=3, font=3, cex=1, 
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{"networx"}
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
\references{
Dress, A.W.M. and Huson, D.H. (2004) Constructicting Splits Graphs \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB)}, \bold{1(3)}, 109--115
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\seealso{
\code{\link{consensusNet}}, \code{\link{hadamard}}, \code{\link{distanceHadamard}}, \code{\link{layout.kamada.kawai}}
}
\examples{
\dontrun{
see example in consensusNet 
example(consensusNet)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

