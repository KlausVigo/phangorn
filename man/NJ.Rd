\name{NJ}
\alias{NJ}
\alias{UNJ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Neighbor-Joining}
\description{
This function performs the neighbor-joining tree estimation of Saitou
and Nei (1987). UNJ is the unweighted version from Gascuel (1997).
}
\usage{
NJ(x)
UNJ(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A distance matrix.}
}

\value{
 an object of class \code{"phylo"}.
}
\references{Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new
method for reconstructing phylogenetic trees. \emph{Molecular Biology
and Evolution}, \bold{4}, 406--425. 

Studier, J. A and Keppler, K. J. (1988) A Note on the Neighbor-Joining
Algorithm of Saitou and Nei. \emph{Molecular Biology and Evolution}, 
\bold{6}, 729--731. 

Gascuel, O. (1997) Concerning the NJ algorithm and its unweighted 
version, UNJ. in Birkin et. al. \emph{Mathematical Hierarchies and Biology}, 
 149--170.}    
\author{Klaus P. Schliep \email{klaus.schliep@gmail.com}}

\seealso{  \code{\link[ape]{nj}}, \code{\link[ape]{dist.dna}}, \code{\link[phangorn]{dist.hamming}}, \code{\link[phangorn]{upgma}}, \code{\link[ape]{fastme}}}
\examples{
data(Laurasiatherian)
dm <- dist.ml(Laurasiatherian)
tree <- NJ(dm)
plot(tree)
}
\keyword{ cluster }% at least one, from doc/KEYWORDS
