\name{parsimony}
\alias{parsimony}
\alias{optim.parsimony}
\alias{sankoff}
\alias{fitch}
\alias{PNJ}
\title{Parsimony tree.}
\alias{CI}
\alias{RI}
\alias{pratchet}
\alias{random.addition}
\alias{acctran}
\description{

\code{parsimony} returns the parsimony score of a tree using either the sankoff or the fitch algorithm.
\code{optim.parsimony} tries to find the maximum parsimony tree using either Nearest Neighbor Interchange (NNI) 
rearrangements or sub tree pruning and regrafting (SPR). \code{pratchet} implements the parsimony ratchet (Nixon, 1999) and is the prefered way to search for the best tree. 
\code{random.addition} can be used to produce starting trees. 
\code{CI} and \code{RI} computes the consistency and retention index.  
}
\usage{
parsimony(tree, data, method="fitch", ...)
optim.parsimony(tree, data, method="fitch", cost=NULL, trace=1, rearrangements="SPR", ...)
pratchet(data, start=NULL, method="fitch", maxit=1000, k=10, trace=1, all=FALSE, 
    rearrangements="SPR", ...)
fitch(tree, data, site = "pscore")
sankoff(tree, data, cost = NULL, site = "pscore")
random.addition(data, method="fitch")
CI(tree, data, cost = NULL)
RI(tree, data, cost = NULL)
acctran(tree, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A object of class phyDat containing sequences.}
  \item{tree}{ tree to start the nni search from.}
  \item{method}{one of 'fitch' or 'sankoff'.}
  \item{cost}{A cost matrix for the transitions between two states.} 
  \item{site}{return either 'pscore' or 'site' wise parsimony scores.} 
  \item{trace}{defines how much information is printed during optimisation.}  
  \item{rearrangements}{SPR or NNI rearrangements.} 
  \item{start}{a starting tree can be supplied.}
  \item{maxit}{maximum number of iterations in the ratchet.}
  \item{k}{number of rounds ratchet is stopped, when there is no improvement.}
  \item{all}{return all equally good trees or just one of them.} 
  \item{...}{Further arguments passed to or from other methods (e.g. model="sankoff" and cost matrix).} 
}  
\value{
  \code{parsimony} returns the maximum parsimony score (pscore). 
  \code{optim.parsimony} returns a tree after NNI rearrangements.
  \code{pratchet} returns a tree or list of trees containing the best tree(s) found during the search.   
  \code{acctran} returns a tree with edge length according to the ACCTRAN criterion. 
}
\details{
The "SPR" rearrangements are so far only available for the "fitch" method, "sankoff" only uses "NNI". The "fitch" algorithm only works correct for binary trees. 
}
\references{
Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer Associates, Sunderland.

Nixon, K. (1999) The Parsimony Ratchet, a New Method for Rapid Parsimony Analysis. \emph{Cladistics} \bold{15}, 407-414
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{\code{\link{bab}}, \code{\link{ancestral.pml}}, \code{\link{nni}}, \code{\link{NJ}}, \code{\link{pml}}, \code{\link{getClans}} ,\code{\link{ancestral.pars}}, \code{\link{bootstrap.pml}}}
\examples{
set.seed(3)
data(Laurasiatherian)
dm = dist.hamming(Laurasiatherian)
tree = NJ(dm)
parsimony(tree, Laurasiatherian)
treeRA <- random.addition(Laurasiatherian)
treeNNI <- optim.parsimony(tree, Laurasiatherian)
treeRatchet <- pratchet(Laurasiatherian, start=tree, maxit=100, k=5)
# assign edge length
treeRatchet <- acctran(treeRatchet, Laurasiatherian)

plot(midpoint(treeRatchet))
add.scale.bar(0,0, length=100)

parsimony(c(tree,treeNNI, treeRatchet), Laurasiatherian)
}
\keyword{cluster}
