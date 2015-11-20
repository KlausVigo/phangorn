\name{delta.score}
\alias{delta.score}

\title{
Computes the \eqn{\delta} score
}
\description{
Computes the treelikeness  
}
\usage{
delta.score(x, arg = "mean", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{phyDat}
}
  \item{arg}{
Specifies the retjurn value, one of "all", "mean" or "sd"
}
  \item{...}{
     further arguments passed through \code{dist.hamming}
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A vector containg the \eqn{\delta} scores. 
}
\references{
BR Holland, KT Huber, A Dress, V Moulton (2002) \eqn{\delta} Plots: a tool for analyzing phylogenetic distance data
Russell D. Gray, David Bryant, Simon J. Greenhill (2010) On the shape and fabric of human history 
\emph{Molecular Biology and Evolution}, \bold{19(12)} 2051--2059

Russell D. Gray, David Bryant, Simon J. Greenhill (2010) On the shape and fabric of human history 
\emph{Phil. Trans. R. Soc. B}, \bold{365} 3923--3933; DOI: 10.1098/rstb.2010.0162
}
\author{
Alastair Potts and Klaus Schliep 
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{dist.hamming}}
}
\examples{
data(yeast)
hist(delta.score(yeast, "all"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
