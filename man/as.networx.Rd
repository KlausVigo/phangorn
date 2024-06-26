% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networx.R
\name{as.networx}
\alias{as.networx}
\alias{networx}
\alias{as.networx.splits}
\alias{as.networx.phylo}
\title{Conversion among phylogenetic network objects}
\usage{
as.networx(x, ...)

\method{as.networx}{splits}(x, planar = FALSE, coord = "none", ...)

\method{as.networx}{phylo}(x, ...)
}
\arguments{
\item{x}{an object of class \code{"splits"} or \code{"phylo"}}

\item{\dots}{Further arguments passed to or from other methods.}

\item{planar}{logical whether to produce a planar graph from only cyclic
splits (may excludes splits).}

\item{coord}{add coordinates of the nodes, allows to reproduce the plot.}
}
\value{
an object of class \code{networx}.
}
\description{
\code{as.networx} convert \code{splits} objects into a \code{networx}
object. And most important there exists a generic \code{plot} function to
plot phylogenetic network or split graphs.
}
\details{
A \code{networx} object hold the information for a phylogenetic
network and extends the \code{phylo} object. Therefore some generic function
for \code{phylo} objects will also work for \code{networx} objects.  The
argument \code{planar = TRUE} will create a planar split graph based on a
cyclic ordering. These objects can be nicely plotted in \code{"2D"}.
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
\dontrun{
# also see example in consensusNet
example(consensusNet)
}

}
\references{
Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W. (2017),
Intertwining phylogenetic trees and networks. \emph{Methods Ecol Evol}.
\bold{8}, 1212--1220. doi:10.1111/2041-210X.12760
}
\seealso{
\code{\link{consensusNet}}, \code{\link{neighborNet}},
\code{\link{splitsNetwork}}, \code{\link{hadamard}},
\code{\link{distanceHadamard}}, \code{\link{plot.networx}},
\code{\link[ape]{evonet}}, \code{\link[ape]{as.phylo}}
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\keyword{plot}
