\name{Ancestors}
\alias{Ancestors}
\alias{Children}
\alias{Descendants}
\alias{Siblings}
\alias{mrca.phylo}
%\alias{mrca}
\title{tree utility function}
\description{
Functions for describing relationships among phylogenetic nodes.
}
\usage{
Ancestors(x, node, type=c("all","parent"))
Children(x, node)
Siblings(x, node, include.self=FALSE)
Descendants(x, node, type=c("tips","children","all"))
mrca.phylo(x, node=NULL)
}
%mrca(x)
\arguments{
  \item{x}{a tree (a phylo object).}
  \item{node}{an integer or a vector of integers corresponding to a node ID}
  \item{type}{specify whether to return just direct
    children / parents or all }
  \item{include.self}{whether to include self in list of siblings}
}
\details{
These functions are inspired by \code{treewalk} in phylobase package, but work on the 
S3 \code{phylo} objects.  
The nodes are the indices as given in edge matrix of an phylo object.
From taxon labels these indices can be easily derived matching against the \code{tip.label}
argument of an phylo object, see example below. 
All the functions allow \code{node} to be either a scalar or vector.  
\code{mrca} is a faster version of the mrca in ape, in phangorn only because of dependencies. 
}
\value{
   a vector or a list containing the indices of the nodes. 
}
\seealso{\code{treewalk}, \code{\link[ape]{phylo}}, \code{\link[ape]{nodelabels}} }
\examples{
tree = rtree(10)
plot(tree, show.tip.label = FALSE)
nodelabels()
tiplabels()
Ancestors(tree, 1:3, "all")
Children(tree, 11)
Descendants(tree, 11, "tips")
Siblings(tree, 3)
mrca.phylo(tree, 1:3)
mrca.phylo(tree, match(c("t1", "t2", "t3"), tree$tip))
mrca.phylo(tree)
# same as mrca(tree), but faster for large trees
}
\keyword{misc}
