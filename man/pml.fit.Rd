\name{pml.fit}
\alias{pml.fit}
\alias{edQt}
\alias{pml.init}
\alias{pml.free}
\alias{discrete.gamma}
\alias{lli}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal maximum likelihood functions.   
}
\description{
These functions are internally used for the likelihood computations in \code{pml} or \code{optim.pml}. 
}
\usage{
pml.fit(tree, data, bf=rep(1/length(levels), length(levels)), shape=1, k=1,
   Q=rep(1, length(levels)*(length(levels)-1)/2), levels=attr(data, "levels"),
   inv=0, rate=1, g=NULL, w=NULL, eig=NULL, INV=NULL, ll.0=NULL, llMix=NULL, 
   wMix=0, ..., site=FALSE)
pml.init(data, k)
pml.free()   
edQt(Q = c(1, 1, 1, 1, 1, 1), bf = c(0.25, 0.25, 0.25, 0.25))
lli(data, tree, ...)
discrete.gamma(alpha, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{A phylogenetic \code{tree}, object of class \code{phylo}. }
  \item{data}{An alignment, object of class \code{phyDat}.}
  \item{bf}{Base frequencies.}
  \item{shape}{Shape parameter of the gamma distribution.}
  \item{alpha}{Shape parameter of the gamma distribution.}
  \item{k}{Number of intervals of the discrete gamma distribution.}
  \item{Q}{A vector containing the lower triangular part of the rate matrix.}
  \item{levels}{
%%     ~~Describe \code{levels} here~~
}
  \item{inv}{Proportion of invariable sites.}
  \item{rate}{Rate.}
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{eig}{Eigenvalue decomposition of Q}
  \item{INV}{Sparse representation of invariant sites}
  \item{ll.0}{
%%     ~~Describe \code{ll.0} here~~  
  }
  \item{llMix}{
%%     ~~Describe \code{llMix} here~~  
  }
  \item{wMix}{
%%     ~~Describe \code{wMix} here~~  
  
  }
  \item{\dots}{Further arguments passed to or from other methods.}
  \item{site}{
%%     ~~Describe \code{site} here~~
}
}
\details{
These functions are exported to be used in different packages so far only in the package coalescentMCMC, but are not intended for end user. Most of the functions call C code and are far less forgiving if the import is not what they expect than \code{pml}. 
}
\value{
\code{pml.fit} returns the log-likelihood. 
}
\references{
Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maximum
likelihood approach. \emph{Journal of Molecular Evolution}, \bold{17}, 368--376. 
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\seealso{
\code{\link{pml}, \link{pmlPart}, \link{pmlMix}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }

