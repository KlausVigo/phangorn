\name{treedist}
\alias{treedist}
\alias{RF.dist}
\alias{wRF.dist}
\alias{KF.dist}
\alias{path.dist}
%\alias{print.treedist}
\title{ Distances between trees }
\description{
  \code{treedist} computes different tree distance methods and \code{RF.dist} the Robinson-Foulds or symmetric distance. 
}
\usage{
treedist(tree1, tree2, check.labels = TRUE)
RF.dist(tree1, tree2=NULL, normalize=FALSE, check.labels=TRUE, rooted=FALSE)
wRF.dist(tree1, tree2=NULL, normalize=FALSE, check.labels=TRUE, rooted=FALSE)
KF.dist(tree1, tree2=NULL, check.labels=TRUE, rooted=FALSE)
path.dist(tree1, tree2=NULL, check.labels=TRUE, use.weight=FALSE)
}
\arguments{
  \item{tree1}{ A phylogenetic tree (class \code{phylo})
  or vector of trees (an object of class \code{multiPhylo}). See details }
  \item{tree2}{ A phylogenetic tree. }
  \item{normalize}{compute normalized RF-distance, see details. }
  \item{check.labels}{compares labels of the trees.}
  \item{rooted}{take bipartitions for rooted trees into account, default is unrooting the trees.}
  \item{use.weight}{use edge.length argument or just count number of edges on the path (default)}
}
\value{
  \code{treedist} returns a vector containing the following tree distance methods
  \item{symmetric.difference}{symmetric.difference or  Robinson-Foulds distance}
  \item{branch.score.difference}{branch.score.difference}
  \item{path.difference}{path.difference}
  \item{weighted.path.difference}{weighted.path.difference}
}
\details{
The Robinson-Foulds distance between two trees \eqn{T_1} and \eqn{T_2} with \eqn{n} tips is defined as (following the notation Steel & Penny 1993): 
\deqn{d(T_1, T_2) = i(T_1) + i(T_2) - 2v_s(T_1, T_2)}
where \eqn{i(T_1)} denotes the number of internal edges and \eqn{v_s(T_1, T_2)} denotes the number of internal splits shared by the two trees. The normalized Robinson-Foulds distance is derived by dividing \eqn{d(T_1, T_2)} by the maximal possible distance \eqn{i(T_1) + i(T_2)}. If both trees are unrooted and binary this value is \eqn{2n-6}. 

% The weighted Robinson-Foulds distance makes use of edge weights. Let \eqn{E_1} and \eqn{E_2} the set of edges of \eqn{d(T_1, T_2)} and \eqn{w(e), e \in E_1} the weight of an edge in \eqn{E_1}
%\deqn{d(T_1, T_2) = \sum_{e \in E \setminus E_1}w(e) + \sum_{e \in E_1\cap E_2}|w(e_1) - w(e_2) + \sum_{e \in E\setminus E_2}w(e) }
% where \eqn{E = E_1\cup E_2}

%The Kuhner-Felsenstein difference (Kuhner & Felsenstein 1994) is closely related with the weighted Robinson-Foulds distance. 

%Instead of the absolute values between edges are squared 
% 0-norm RF
% 1-norm wRF
% 2-norm KF 
 
\code{RF.dist} returns the Robinson-Foulds distance (Robinson & Foulds 1981) between either 2 trees or computes a matrix of all pairwise distances if a \code{multiPhylo} object is given. The Robinson-Foulds distance only depends on the toplogy of the trees. 
If edge weights should be considered \code{wRF.dist} calculates the weighted RF distance (Robinson & Foulds 1981). and \code{KF.dist} calculates the branch score distance (Kuhner & Felsenstein 1994). 
\code{path.dist} computes the path difference metric as described in Steel and Penny 1993).

For large number of trees the distance functions can use a lot of memory!

% The function used internally is 2 * (nt - m) where nt is the number of tips and 
% m is the number of shared bipartitions. When there are multifurcations the 
% distance is therefore increasing!! This may be different to other implementations!

}
\references{
Steel M. A. and Penny P. (1993) \emph{Distributions of tree comparison metrics - some new results}, Syst. Biol., \bold{42(2)}, 126--141

Kuhner, M. K. and Felsenstein, J. (1994) \emph{A simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates}, Molecular Biology and Evolution, \bold{11(3)}, 459--468

D.F. Robinson and L.R. Foulds (1981) \emph{Comparison of phylogenetic trees}, Mathematical Biosciences, \bold{53(1)}, 131--147

D.F. Robinson and L.R. Foulds (1979) Comparison of weighted labelled trees. In Horadam, A. F. and Wallis, W. D. (Eds.), \emph{Combinatorial Mathematics VI: Proceedings of the Sixth Australian Conference on Combinatorial Mathematics, Armidale, Australia}, 119--126
}
\author{ Klaus P. Schliep \email{klaus.schliep@gmail.com}} 
\seealso{\code{\link[ape]{dist.topo}}}
\examples{
tree1 <- rtree(100, rooted=FALSE)
tree2 <- rSPR(tree1, 3)
RF.dist(tree1, tree2)
treedist(tree1, tree2)
}
\keyword{ classif }% at least one, from doc/KEYWORDS
