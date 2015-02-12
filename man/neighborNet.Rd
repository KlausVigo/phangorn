\name{neighborNet}
\alias{neighborNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes a neighborNet from a distance matrix
}
\description{
Computes a neighborNet, i.e. an object of class \code{networx} from a distance matrix.
}
\usage{
neighborNet(x, ord = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a distance matrix.} 
  \item{ord}{a circular ordering.} 
}
\details{
\code{neighborNet} is still experimental. The cyclic ordering sometimes differ from the SplitsTree implementation, the \emph{ord} argument can be used to enforce a certain circular ordering.   
}
\value{
\code{neighborNet} returns an object of class networx.   
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{ 
Bryant, D. & Moulton, V. (2004) Neighbor-Net: An Agglomerative Method for the Construction of Phylogenetic Networks. \emph{Molecular Biology and Evolution}, 2004, \bold{21}, 255-265
}
\seealso{
\code{\link{splitsNetwork}}, \code{\link{consensusNet}}, \code{\link{plot.networx}}, \code{\link{lento}}
}
\examples{
data(yeast)
dm <- dist.ml(yeast)
nnet <- neighborNet(dm)
plot(nnet, "2D")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

