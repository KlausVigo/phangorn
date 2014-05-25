\name{dist.hamming}
\alias{dist.hamming}
\alias{dist.logDet}
\alias{dist.ml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pairwise Distances from Sequences}
\description{
\code{dist.hamming} and \code{dist.logDet} compute pairwise distances 
for an object of class \code{phyDat}. \code{dist.ml} fits distances
for nucleotide and amino acid models.
}
\usage{
dist.hamming(x, ratio = TRUE, exclude="none")
dist.logDet(x)
dist.ml(x, model="JC69", exclude="none", bf=NULL, Q=NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{phyDat}}
  \item{ratio}{Compute uncorrected ('p') distance or character difference.}
  \item{model}{One of "JC69" or one of 17 amino acid models see details.}
  \item{exclude}{One of "none", "all", "pairwise" indicating whether to delete the sites with missing data (or ambigious states). The default is handle missing data as in pml.}
  \item{bf}{A vector of base frequencies.}
  \item{Q}{A vector containing the lower triangular part of the rate matrix.} 
  \item{\dots}{Further arguments passed to or from other methods.}
}
\value{
  an object of class \code{dist}
}
\details{
So far 17 amino acid models are supported ("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb", "FLU", "Blossum62", "Dayhoff_DCMut" and "JTT_DCMut") and additional rate matrices and frequences can be supplied. 
}
\references{ 
Lockhart, P. J., Steel, M. A., Hendy, M. D. and Penny, D. (1994)
Recovering evolutionary trees under a more realistic model of sequence
evolution. \emph{Molecular Biology and Evolution}, \bold{11}, 605--602.
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}

\seealso{For more distance methods for nucleotide data see \code{\link[ape]{dist.dna}}}
\examples{
data(Laurasiatherian)
dm1 <- dist.hamming(Laurasiatherian)
tree1 <- NJ(dm1)
dm2 <- dist.logDet(Laurasiatherian)
tree2 <- NJ(dm2)
treedist(tree1,tree2)
}

\keyword{ cluster }
