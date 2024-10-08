% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{bootstrap.pml}
\alias{bootstrap.pml}
\alias{bootstrap.phyDat}
\title{Bootstrap}
\usage{
bootstrap.pml(x, bs = 100, trees = TRUE, multicore = FALSE,
  mc.cores = NULL, tip.dates = NULL, ...)

bootstrap.phyDat(x, FUN, bs = 100, multicore = FALSE, mc.cores = NULL,
  jumble = TRUE, ...)
}
\arguments{
\item{x}{an object of class \code{pml} or \code{phyDat}.}

\item{bs}{number of bootstrap samples.}

\item{trees}{return trees only (default) or whole \code{pml} objects.}

\item{multicore}{logical, whether models should estimated in parallel.}

\item{mc.cores}{The number of cores to use during bootstrap. Only supported
on UNIX-alike systems.}

\item{tip.dates}{A named vector of sampling times associated to the
tips/sequences. Leave empty if not estimating tip dated phylogenies.}

\item{\dots}{further parameters used by \code{optim.pml} or
\code{plot.phylo}.}

\item{FUN}{the function to estimate the trees.}

\item{jumble}{logical, jumble the order of the sequences.}
}
\value{
\code{bootstrap.pml} returns an object of class \code{multi.phylo}
or a list where each element is an object of class \code{pml}. \code{plotBS}
returns silently a tree, i.e. an object of class \code{phylo} with the
bootstrap values as node labels. The argument \code{BStrees} is optional and
if not supplied the tree with labels supplied in the \code{node.label} slot.
}
\description{
\code{bootstrap.pml} performs (non-parametric) bootstrap analysis and
\code{bootstrap.phyDat} produces a list of bootstrapped data sets.
\code{plotBS} plots a phylogenetic tree with the bootstrap values assigned
to the (internal) edges.
}
\details{
It is possible that the bootstrap is performed in parallel, with help of the
multicore package. Unfortunately the multicore package does not work under
windows or with GUI interfaces ("aqua" on a mac). However it will speed up
nicely from the command line ("X11").
}
\examples{

\dontrun{
data(Laurasiatherian)
dm <- dist.hamming(Laurasiatherian)
tree <- NJ(dm)
# NJ
set.seed(123)
NJtrees <- bootstrap.phyDat(Laurasiatherian,
     FUN=function(x)NJ(dist.hamming(x)), bs=100)
treeNJ <- plotBS(tree, NJtrees, "phylogram")

# Maximum likelihood
fit <- pml(tree, Laurasiatherian)
fit <- optim.pml(fit, rearrangement="NNI")
set.seed(123)
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
treeBS <- plotBS(fit$tree,bs)

# Maximum parsimony
treeMP <- pratchet(Laurasiatherian)
treeMP <- acctran(treeMP, Laurasiatherian)
set.seed(123)
BStrees <- bootstrap.phyDat(Laurasiatherian, pratchet, bs = 100)
treeMP <- plotBS(treeMP, BStrees, "phylogram")
add.scale.bar()

# export tree with bootstrap values as node labels
# write.tree(treeBS)
}

}
\references{
Felsenstein J. (1985) Confidence limits on phylogenies. An
approach using the bootstrap. \emph{Evolution} \bold{39}, 783--791

Lemoine, F., Entfellner, J. B. D., Wilkinson, E., Correia, D., Felipe, M. D.,
De Oliveira, T., & Gascuel, O. (2018). Renewing Felsenstein’s phylogenetic
bootstrap in the era of big data. \emph{Nature}, \bold{556(7702)}, 452--456.

Penny D. and Hendy M.D. (1985) Testing methods evolutionary tree
construction. \emph{Cladistics} \bold{1}, 266--278

Penny D. and Hendy M.D. (1986) Estimating the reliability of evolutionary
trees. \emph{Molecular Biology and Evolution} \bold{3}, 403--417
}
\seealso{
\code{\link{optim.pml}}, \code{\link{pml}},
\code{\link[ape]{plot.phylo}}, \code{\link{maxCladeCred}}
\code{\link[ape]{nodelabels}},\code{\link{consensusNet}} and
\code{\link{SOWH.test}} for parametric bootstrap
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\keyword{cluster}
