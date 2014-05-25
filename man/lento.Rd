\name{lento}
\alias{lento}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lento plot
}
\description{
The lento plot represents support and conflict of splits/bipartitions. 
}
\usage{
lento(obj, xlim = NULL, ylim = NULL, main = "Lento plot", sub = NULL, xlab = NULL, 
    ylab = NULL, bipart=TRUE, trivial=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{an object of class phylo, multiPhylo or splits}
  \item{xlim}{graphical parameter}
  \item{ylim}{graphical parameter}
  \item{main}{graphical parameter}
  \item{sub}{graphical parameter}
  \item{xlab}{graphical parameter}
  \item{ylab}{graphical parameter}
  \item{bipart}{plot bipartition information.}
  \item{trivial}{logical, whether to present trivial splits (default is FALSE).}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\value{
lento returns a plot.
}
\references{
Lento, G.M., Hickson, R.E., Chambers G.K., and Penny, D. (1995) Use of spectral analysis to test hypotheses on the origin of pinninpeds. \emph{Molecular Biology and Evolution}, \bold{12}, 28-52.  
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}

\seealso{
\code{\link{as.splits}, \link{hadamard}}
}
\examples{
data(yeast)
yeast.ry = acgt2ry(yeast)
splits.h = h2st(yeast.ry)
lento(splits.h, trivial=TRUE) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ plot }% __ONLY ONE__ keyword per line
