\name{coalSpeciesTree}
\alias{coalSpeciesTree}
\title{
Species Tree
}
\description{
\code{coalSpeciesTree} estimates species trees and can handle multiple individuals per species.
}
\usage{
coalSpeciesTree(tree, X, sTree = NULL)
}
\arguments{
\item{tree}{
an object of class \code{multiPhylo}
}
\item{X}{
A \code{phyDat} object to define which individual belongs to which species. 
} 
\item{sTree}{
A species tree which fixes the topology. 
} 
}
\details{
\code{coalSpeciesTree} estimates a single linkage tree as suggested by Liu et al. (2010) from the element wise minima of the cophenetic matrices of the gene trees. It extends \code{speciesTree} in ape as it allows that have several individuals per gene tree.  
}
\value{
The function returns an object of class \code{phylo}. 
}
\references{
Liu, L., Yu, L. and Pearl, D. K. (2010) Maximum tree: a consistent estimator of the species tree. \emph{Journal of Mathematical Biology}, \bold{60}, 95--106.
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
Emmanuel Paradies
}

\seealso{
\code{\link{speciesTree}}
}
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
