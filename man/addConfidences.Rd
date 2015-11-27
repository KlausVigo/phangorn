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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obj, phy) 
UseMethod("addConfidences")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
