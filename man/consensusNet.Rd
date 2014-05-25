\name{consensusNet}
\alias{consensusNet}
\alias{as.networx}
\alias{as.networx.splits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes a \code{networx} object from a collection of splits.
}
\description{
Computes a \code{networx} object from a collection of splits.
}
\usage{
consensusNet(obj, prob=.3, ...)
as.networx(x, ...)
\method{as.networx}{splits}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of class multiPhylo.} 
  \item{prob}{the proportion a split has to be present in all trees to be represented in the network.} 
  \item{x}{An object of class splits.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\value{
\code{as.networx} returns an object of class networx. 
This is just an intermediate to plot phylogenetic networks with igraph.   
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}

\note{
The internal representation is likely to change. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{ 
Holland B.R., Huber K.T., Moulton V., Lockhart P.J. (2004) Using consensus networks to visualize contradictory evidence for species phylogeny. \emph{Molecular Biology and Evolution}, \bold{21}, 1459--61
}
\seealso{
\code{\link{splitsNetwork}}, \code{\link{lento}}, \code{\link{distanceHadamard}}
}
\examples{
data(Laurasiatherian)
set.seed(1)
bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)nj(dist.hamming(x)), 
    bs=100, multicore=FALSE)
class(bs) <- 'multiPhylo'
cnet = consensusNet(bs, .3)
plot(cnet, show.edge.label=TRUE)
open3d()
plot(cnet, show.tip.label=FALSE, show.nodes=TRUE)
#plot(cnet, type = "2D", show.edge.label=TRUE)
set.seed(1)
tree1 = rtree(20, rooted=FALSE)
sp = as.splits(rNNI(tree1, n=10))
net = as.networx(sp)
open3d()
plot(net)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

