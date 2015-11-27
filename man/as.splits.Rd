\name{as.splits}
\alias{as.splits}
\alias{as.prop.part.splits}
\alias{as.splits.phylo}
\alias{as.splits.multiPhylo}
\alias{as.splits.networx}
\alias{as.matrix.splits}
\alias{as.Matrix}
\alias{as.Matrix.splits}
\alias{print.splits}
\alias{write.splits}
\alias{allSplits}
\alias{compatible}
\alias{write.nexus.splits}
\alias{read.nexus.splits}
\alias{as.phylo.splits}
\alias{countCycles}
\alias{addTrivialSplits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splits representation of graphs and trees.
}
\description{
\code{as.splits} produces a list of splits or bipartitions. 
}
\usage{
as.splits(x, ...)
\method{as.splits}{phylo}(x, ...)
\method{as.splits}{multiPhylo}(x, ...)
\method{print}{splits}(x, maxp = getOption("max.print"), zero.print = ".", 
    one.print = "|", ...)
\method{as.prop.part}{splits}(x, ...)    
compatible(obj)
allSplits(k, labels = NULL)
write.nexus.splits(obj, file="", weights=NULL)
read.nexus.splits(file)
addTrivialSplits(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class phylo or multiPhylo.}
  \item{maxp}{integer, default from \code{options(max.print)}, influences how many entries of large matrices are printed at all.} 
  \item{zero.print}{character which should be printed for zeroes.} 
  \item{one.print}{character which should be printed for ones.} 
  \item{\dots}{Further arguments passed to or from other methods.}
  \item{obj}{an object of class splits.} 
  \item{k}{number of taxa.}
  \item{labels}{names of taxa.}   
  \item{file}{ a file name.}
  \item{weights}{ Edge weights.}
}
\value{
\code{as.splits} returns an object of class splits, which is mainly a list of splits and some attributes.
Often a \code{splits} object will contain attributes \code{confidences} for bootstrap or Bayesian support values and 
\code{weight} storing edge weights. \code{read.nexus.splits} and \code{write.nexus.splits} allow to exchange a
\code{splits} object with other software like Splitstree.

\code{compatible} return a lower triangular matrix where an 1 indicates that two splits are incompatible.   
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\note{
The internal representation is likely to change. 
\code{read.nexus.splits} reads in the splits block of a nexus file. It assumes that different co-variables are tab delimited and the bipartition are separated with white-space. Comments in square brackets are ignored. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{prop.part}}, \code{\link{lento}}, \code{\link{distanceHadamard}}, \code{\link{as.networx}}
}
%% as.phylo.splits
\examples{
(sp <- as.splits(rtree(5)))
write.nexus.splits(sp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }

