\name{maxCladeCred}
\alias{maxCladeCred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum clade credibility tree
}
\description{
\code{maxCladeCred} computes the maximum clade credibility tree from a sample of trees. 
}
\usage{
maxCladeCred(x, tree = TRUE, part = NULL, rooted = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  \code{x} is an object of class \code{multiPhylo} or \code{phylo}
}
  \item{tree}{
logical indicating whether return the tree with the clade credibility (default) or 
the clade credibility score for all trees.
}
  \item{rooted}{
logical, if FALSE the tree with highest maximum bipartition credibility is returned.
}
  \item{part}{a list of partitions as returned by \code{prop.part}}
}
\details{
So far just the best tree is returned. No annotations or transformations of edge length are performed. 

If a list of partition is provided then the clade credibility is computed for the trees in x.
}
\value{
a tree (an object of class \code{phylo}) with the highest clade credibility or a numeric vector of clade credibilities for each tree.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com} 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{consensus}}, \code{\link{consensusNet}}, \code{\link{prop.part}}
}
\examples{

data(Laurasiatherian)
set.seed(42)
bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)upgma(dist.hamming(x)), 
    bs=100)
class(bs) <- 'multiPhylo'

strict_consensus <- consensus(bs)
majority_consensus <- consensus(bs, p=.5)
max_clade_cred <- maxCladeCred(bs)
par(mfrow = c(1,3), mar = c(1,4,1,1))
plot(strict_consensus, main="Strict consensus tree")
plot(majority_consensus, main="Majority consensus tree")
plot(max_clade_cred, main="Maximum clade credibility tree")

# compute clade credibility for trees given a prop.part object
pp <- prop.part(bs)
tree <- rNNI(bs[[1]], 20)
maxCladeCred(c(tree, bs[[1]]), tree=FALSE, part = pp)
# first value likely be -Inf
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }

