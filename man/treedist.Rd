% Please edit documentation in R/treedist.R
\name{treedist}
\alias{KF.dist}
\alias{RF.dist}
\alias{SPR.dist}
\alias{path.dist}
\alias{sprdist}
\alias{treedist}
\alias{wRF.dist}
\title{Distances between trees}
\usage{
treedist(tree1, tree2, check.labels = TRUE)

sprdist(tree1, tree2)

SPR.dist(tree1, tree2 = NULL)

RF.dist(tree1, tree2 = NULL, normalize = FALSE, check.labels = TRUE, rooted = FALSE)

wRF.dist(tree1, tree2 = NULL, normalize = FALSE, check.labels = TRUE, rooted = FALSE)

KF.dist(tree1, tree2 = NULL, check.labels = TRUE, rooted = FALSE)

path.dist(tree1, tree2 = NULL, check.labels = TRUE, use.weight = FALSE)
}
\arguments{
\item{tree1}{A phylogenetic tree (class \code{phylo}) or vector of trees (an
object of class \code{multiPhylo}). See details}

\item{tree2}{A phylogenetic tree.}

\item{check.labels}{compares labels of the trees.}

\item{normalize}{compute normalized RF-distance, see details.}

\item{rooted}{take bipartitions for rooted trees into account, default is
unrooting the trees.}

\item{use.weight}{use edge.length argument or just count number of edges on
the path (default)}
}
\value{
\code{treedist} returns a vector containing the following tree
distance methods \item{symmetric.difference}{symmetric.difference or
Robinson-Foulds distance}
\item{branch.score.difference}{branch.score.difference}
\item{path.difference}{path.difference}
\item{weighted.path.difference}{weighted.path.difference}
}
\description{
\code{treedist} computes different tree distance methods and \code{RF.dist}
the Robinson-Foulds or symmetric distance.
}
\examples{

tree1 <- rtree(100, rooted = FALSE)
tree2 <- rSPR(tree1, 3)
RF.dist(tree1, tree2)
treedist(tree1, tree2)
sprdist(tree1, tree2)
trees <- rSPR(tree1, 1:5)
SPR.dist(tree1, trees)

}
\author{
Klaus P. Schliep \email{klaus.schliep@gmail.com},

Leonardo de Oliveira Martins
}
\references{
de Oliveira Martins L., Leal E., Kishino H. (2008)
\emph{Phylogenetic Detection of Recombination with a Bayesian Prior on the
Distance between Trees}. PLoS ONE \bold{3(7)}. e2651. doi:
10.1371/journal.pone.0002651

de Oliveira Martins L., Mallo D., Posada D. (2016) \emph{A Bayesian
Supertree Model for Genome-Wide Species Tree Reconstruction}. Syst. Biol.
\bold{65(3)}: 397-416, doi:10.1093/sysbio/syu082

Steel M. A. and Penny P. (1993) \emph{Distributions of tree comparison
metrics - some new results}, Syst. Biol., \bold{42(2)}, 126--141

Kuhner, M. K. and Felsenstein, J. (1994) \emph{A simulation comparison of
phylogeny algorithms under equal and unequal evolutionary rates}, Molecular
Biology and Evolution, \bold{11(3)}, 459--468

D.F. Robinson and L.R. Foulds (1981) \emph{Comparison of phylogenetic
trees}, Mathematical Biosciences, \bold{53(1)}, 131--147

D.F. Robinson and L.R. Foulds (1979) Comparison of weighted labelled trees.
In Horadam, A. F. and Wallis, W. D. (Eds.), \emph{Combinatorial Mathematics
VI: Proceedings of the Sixth Australian Conference on Combinatorial
Mathematics, Armidale, Australia}, 119--126
}
\seealso{
\code{\link[ape]{dist.topo}}, \code{\link{nni}},
\code{\link{superTree}}
}
\keyword{classif}
