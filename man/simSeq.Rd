\name{simSeq}
\alias{simSeq}
\alias{simSeq.phylo}
\alias{simSeq.pml}
\title{ Simulate sequences. }
\description{
Simulate sequences for a given evolutionary tree.
}
\usage{
simSeq(x, ...)
\method{simSeq}{phylo}(x, l=1000, Q=NULL, bf=NULL, rootseq=NULL, type="DNA",
    model="", levels=NULL, rate=1, ancestral=FALSE, ...)
\method{simSeq}{pml}(x, ancestral = FALSE, ...)    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{ a phylogenetic tree \code{tree}, i.e. an object of class \code{phylo} or and object of class \code{pml}. }
\item{l}{  length of the sequence to simulate. }
\item{Q}{ the rate matrix.  }
\item{bf}{ base frequencies. }
\item{rootseq}{a vector of length l containing the root sequence, 
other root sequence is randomly generated.}
\item{type}{Type of sequences ("DNA", "AA" or "USER").}
\item{model}{Amino acid models: one of "WAG", "JTT", "Dayhoff" or "LG"}
\item{levels}{ \code{levels} takes a character vector of the different bases,
default is for nucleotide sequences, only used when type = "USER".}
\item{rate}{rate. }
\item{ancestral}{Return ancestral sequences?}
\item{\dots}{Further arguments passed to or from other methods.}
}
\details{
\code{simSeq} is now a generic function to simulate sequence alignments. 
It is quite flexible and allows to generate DNA, RNA, amino acids or binary sequences. 
It is possible to give a \code{pml} object as input simSeq return a \code{phyDat}
from these model. 
There is also a more low level version, which lacks rate variation, but one can combine different alignments having their own rate (see example).  
}
\value{
\code{simSeq} returns an object of class phyDat.
}
\author{ Klaus Schliep \email{klaus.schliep@gmail.com}}

\seealso{  \code{\link{phyDat}}, \code{\link{pml}},  \code{\link{SOWH.test}}  
}
\examples{
\dontrun{
data(Laurasiatherian)
tree = nj(dist.ml(Laurasiatherian))
fit = pml(tree, Laurasiatherian, k=4)
fit = optim.pml(fit, optNni=TRUE, model="GTR", optGamma=TRUE)
data = simSeq(fit)
}

tree = rtree(5)
plot(tree)
nodelabels()

# Example for simple DNA alignment
data = simSeq(tree, l = 10, type="DNA", bf=c(.1,.2,.3,.4), Q=1:6)
as.character(data)

# Example to simulate discrete Gamma rate variation
rates = phangorn:::discrete.gamma(1,4)
data1 = simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[1])
data2 = simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[2])
data3 = simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[3])
data4 = simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[4])
data <- c(data1,data2, data3, data4)

write.phyDat(data, file="temp.dat", format="sequential",nbcol = -1, colsep = "")
unlink("temp.dat") 
}
\keyword{ cluster }% at least one, from doc/KEYWORDS

