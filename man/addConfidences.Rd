\name{addConfidences}
\alias{addConfidences}
\alias{presenceAbsence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare splits and add support values to an object
}
\description{
Add support values to a \code{splits}, \code{phylo} or \code{networx} object. 
}
\usage{
addConfidences(x, y)
presenceAbsence(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{splits}, \code{phylo} or \code{networx}
}
  \item{y}{
an object of class \code{splits}, \code{phylo}, \code{multiPhylo} or \code{networx}
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The object \code{x} with added bootstrap / MCMC support values. 
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{as.splits}}, \code{\link{as.networx}}, \code{\link{RF.dist}}, \code{\link{plot.phylo}}
}
\examples{
data(woodmouse)
woodmouse <- phyDat(woodmouse)
tmpfile <- normalizePath(system.file("extdata/trees/RAxML_bootstrap.woodmouse", package="phangorn"))
boot_trees <- read.tree(tmpfile)

dm <- dist.ml(woodmouse)
tree <- upgma(dm)
nnet <- neighborNet(dm)

tree <- addConfidences(tree, boot_trees)
nnet <- addConfidences(nnet, boot_trees)

plot(tree, show.node.label=TRUE)
plot(nnet, "2D", show.edge.label=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
