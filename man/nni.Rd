\name{nni}
\alias{nni}
\alias{rNNI}
\alias{rSPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tree rearrangements.}
\description{
\code{nni} returns a list of all trees which are one nearest neighbor interchange away. \code{rNNI} and \code{rSPR} are two methods which simulate random trees which are a specified number of rearrangement apart from the input tree. Both methods assume that the input tree is bifurcating. These methods may be useful in simulation studies.
}
\usage{
nni(tree)
rSPR(tree, moves=1, n=length(moves), k=NULL)
rNNI(tree, moves=1, n=length(moves))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{A phylogenetic \code{tree}, object of class \code{phylo}.}
  \item{moves}{Number of tree rearrangements to be transformed on a tree.
  Can be a vector}
  \item{n}{Number of trees to be simulated.}
  \item{k}{If defined just SPR of distance k are performed.}
}
\value{
   an object of class multiPhylo. 
}
\author{ Klaus Schliep \email{klaus.schliep@gmail.com} }

% \seealso{ }
\examples{
tree = unroot(rtree(20))
trees1 <- nni(tree)
trees2 <- rSPR(tree, 2, 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
