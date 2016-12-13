% Please edit documentation in R/networx.R
\name{as.splits}
\alias{addTrivialSplits}
\alias{allCircularSplits}
\alias{allSplits}
\alias{as.Matrix}
\alias{as.Matrix.splits}
\alias{as.bitsplits.splits}
\alias{as.matrix.splits}
\alias{as.phylo.splits}
\alias{as.prop.part.splits}
\alias{as.splits}
\alias{as.splits.multiPhylo}
\alias{as.splits.networx}
\alias{as.splits.phylo}
\alias{c.splits}
\alias{compatible}
\alias{countCycles}
\alias{distinct.splits}
\alias{matchSplits}
\alias{print.splits}
\alias{read.nexus.splits}
\alias{write.nexus.splits}
\alias{write.splits}
\title{Splits representation of graphs and trees.}
\usage{
as.splits(x, ...)

\method{print}{splits}(x, maxp = getOption("max.print"), zero.print = ".", 
    one.print = "|", ...)

\method{as.splits}{phylo}(x, ...)

\method{as.splits}{multiPhylo}(x, ...)

\method{as.prop.part}{splits}(x, ...)

compatible(obj)

allSplits(k, labels = NULL)

allCircularSplits(k, labels = NULL)

write.nexus.splits(obj, file = "", weights = NULL, taxa = TRUE, append = FALSE)

read.nexus.splits(file)
}
\arguments{
\item{x}{An object of class phylo or multiPhylo.}

\item{maxp}{integer, default from \code{options(max.print)}, influences how
many entries of large matrices are printed at all.}

\item{zero.print}{character which should be printed for zeros.}

\item{one.print}{character which should be printed for ones.}

\item{obj}{an object of class splits.}

\item{k}{number of taxa.}

\item{labels}{names of taxa.}

\item{file}{a file name.}

\item{weights}{Edge weights.}

\item{taxa}{logical. If TRUE a taxa block is added}

\item{append}{logical. If TRUE the nexus blocks will be added to a file.}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
\code{as.splits} returns an object of class splits, which is mainly
a list of splits and some attributes. Often a \code{splits} object will
contain attributes \code{confidences} for bootstrap or Bayesian support
values and \code{weight} storing edge weights. \code{read.nexus.splits} and
\code{write.nexus.splits} allow to exchange a \code{splits} object with
other software like Splitstree.

\code{compatible} return a lower triangular matrix where an 1 indicates that
two splits are incompatible.
}
\description{
\code{as.splits} produces a list of splits or bipartitions.
}
\note{
The internal representation is likely to change.
\code{read.nexus.splits} reads in the splits block of a nexus file. It
assumes that different co-variables are tab delimited and the bipartition
are separated with white-space. Comments in square brackets are ignored.
}
\examples{

(sp <- as.splits(rtree(5)))
write.nexus.splits(sp)
spl <- allCircularSplits(5)
plot(as.networx(spl), "2D")

}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\seealso{
\code{\link{prop.part}}, \code{\link{lento}},
\code{\link{distanceHadamard}}, \code{\link{as.networx}}
}
\keyword{cluster}
