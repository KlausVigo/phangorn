\name{bootstrap.pml}
%\Rdversion{1.1}
\alias{bootstrap.pml}
\alias{bootstrap.phyDat}
\alias{plotBS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap }
\description{
\code{bootstrap.pml} performs (non-parametric) bootstrap analysis and \code{bootstrap.phyDat} produces a list of bootstrapped data sets. \code{plotBS} plots a phylogenetic tree with the with the bootstrap values assigned to the (internal) edges.}
\usage{
bootstrap.pml(x, bs = 100, trees = TRUE, multicore=FALSE, mc.cores = NULL, ...)
bootstrap.phyDat(x, FUN, bs = 100, multicore=FALSE, mc.cores = NULL, ...)
plotBS(tree, BStrees, type="unrooted", bs.col="black", bs.adj=NULL, p=80, ...)
}
\arguments{
  \item{x}{
an object of class \code{pml} or \code{phyDat}.
}
  \item{bs}{
number of bootstrap samples.
}
  \item{trees}{
return trees only (default) or whole \code{pml} objects.  
}
%  \item{multicore}{
%logical, if TRUE analysis is performed in parallel (see details).
%}
\item{multicore}{logical, whether models should estimated in parallel.}
  \item{mc.cores}{
The number of cores to use during bootstrap. Only supported on UNIX-alike systems.
}
  \item{\dots}{
further parameters used by \code{optim.pml} or \code{plot.phylo}.
}
  \item{FUN}{ 
the function to estimate the trees.
}
  \item{tree}{ 
The tree on which edges the bootstrap values are plotted.
}
  \item{BStrees}{ 
a list of trees (object of class "multiPhylo"). 
}
  \item{type}{ 
the type of tree to plot, so far "cladogram", "phylogram" and "unrooted" are supported.
}
  \item{bs.col}{ 
color of bootstrap support labels.
}
  \item{bs.adj}{ 
one or two numeric values specifying the horizontal and vertical justification of the bootstrap labels.
}
  \item{p}{ 
only plot support values higher than this percentage number (default is 80).
}
}
\details{
It is possible that the bootstrap is performed in parallel, with help of the multicore package.
Unfortunately the  multicore package does not work under windows or with GUI interfaces ("aqua" on a mac).
However it will speed up nicely from the command line ("X11"). 
}
\value{
\code{bootstrap.pml} returns an object of class \code{multi.phylo} or a list where each 
element is an object of class \code{pml}. \code{plotBS} returns silently a tree, i.e. an object of class \code{multi.phylo} with the bootstrap values as node labels.
}
\references{
Felsenstein J. (1985) Confidence limits on phylogenies. An approach using the bootstrap. \emph{Evolution} \bold{39}, 783--791

Penny D. and Hendy M.D. (1985) Testing methods evolutionary tree construction. \emph{Cladistics} \bold{1}, 266--278

Penny D. and Hendy M.D. (1986) Estimating the reliability of evolutionary trees. \emph{Molecular Biology and Evolution} \bold{3}, 403--417
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\seealso{
\code{\link{optim.pml}}, \code{\link{pml}}, \code{\link{plot.phylo}}, \code{\link{consensusNet}}
}
\examples{
\dontrun{
data(Laurasiatherian)
dm <- dist.logDet(Laurasiatherian)
tree <- NJ(dm)
fit=pml(tree,Laurasiatherian)
fit = optim.pml(fit,TRUE)

set.seed(123)
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
treeBS <- plotBS(fit$tree,bs)



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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

