\name{consensusNet}
\alias{consensusNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes a consensusNetwork from a list of trees
Computes a \code{networx} object from a collection of splits.
}
\description{
Computes a consensusNetwork, i.e. an object of class \code{networx} from a list of trees,  i.e. an class of class \code{multiPhylo}.
Computes a \code{networx} object from a collection of splits.
}
\usage{
consensusNet(obj, prob=.3, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of class multiPhylo.} 
  \item{prob}{the proportion a split has to be present in all trees to be represented in the network.} 
  \item{\dots}{Further arguments passed to or from other methods.}
}
\value{
\code{consensusNet} returns an object of class networx. 
This is just an intermediate to plot phylogenetic networks with igraph.   
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{ 
Holland B.R., Huber K.T., Moulton V., Lockhart P.J. (2004) Using consensus networks to visualize contradictory evidence for species phylogeny. \emph{Molecular Biology and Evolution}, \bold{21}, 1459--61
}
\seealso{
\code{\link{splitsNetwork}}, \code{\link{neighborNet}}, \code{\link{lento}}, \code{\link{distanceHadamard}}, \code{\link{plot.networx}}, \code{\link{maxCladeCred}}
}
\examples{
data(Laurasiatherian)
set.seed(1)
bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)nj(dist.hamming(x)), 
    bs=50)
class(bs) <- 'multiPhylo'
cnet = consensusNet(bs, .3)
plot(cnet, "2D")
\dontrun{
library(rgl)
open3d()
plot(cnet, show.tip.label=FALSE, show.nodes=TRUE)
plot(cnet, type = "2D", show.edge.label=TRUE)

tmpfile <- normalizePath(system.file("extdata/trees/RAxML_bootstrap.woodmouse", package="phangorn"))
trees <- read.tree(tmpfile)
cnet_woodmouse = consensusNet(trees, .3)
plot(cnet_woodmouse, type = "2D", show.edge.label=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

