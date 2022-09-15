## include in addConfidences, plotBS etc.
#' Transfer Bootstrap
#'
#' \code{transferBootstrap} assignes transfer bootstrap (Lemoine et al. 2018)
#' values to the (internal) edges.
#'
#' @param tree The tree on which edges the bootstrap values are plotted.
#' @param BStrees a list of trees (object of class "multiPhylo").
#' @return \code{plotBS} returns silently a tree, i.e. an object of class
#' \code{phylo} with the bootstrap values as node labels. The argument
#' \code{BSTrees} is optional and if not supplied the labels supplied
#' in the \code{node.label} slot will be used.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso  \code{\link{plotBS}},
#' @references Lemoine, F., Entfellner, J. B. D., Wilkinson, E., Correia, D.,
#' Felipe, M. D., De Oliveira, T., & Gascuel, O. (2018). Renewing Felsenstein’s
#' phylogenetic bootstrap in the era of big data. \emph{Nature},
#' \bold{556(7702)}, 452--456.
#' @examples
#' fdir <- system.file("extdata/trees", package = "phangorn")
#' # RAxML best-known tree with bipartition support (from previous analysis)
#' raxml.tree <- read.tree(file.path(fdir,"RAxML_bipartitions.woodmouse"))
#' # RAxML bootstrap trees (from previous analysis)
#' raxml.bootstrap <- read.tree(file.path(fdir,"RAxML_bootstrap.woodmouse"))
#'
#' tree_tbe <- transferBootstrap(raxml.tree,  raxml.bootstrap)
#' par(mfrow=c(1,2))
#' plotBS(tree_tbe)
#' # same as
#' plotBS(raxml.tree,  raxml.bootstrap, "p", "TBE")
#' @export
transferBootstrap <- function(tree, BStrees){
  if(!inherits(BStrees, "multiPhylo"))
    stop("BSTrees needs to be of class multiPhylo!")
  BStrees <- .uncompressTipLabel(BStrees)
  BStrees <- .compressTipLabel(BStrees, tree$tip.label)
  BStrees <- reorder(BStrees, "postorder")
  l <- Ntip(tree)
  bp <- prop.part(tree)
  bp <- SHORTwise(bp)[-1]
  not_cherry <- lengths(bp) != 2
  res <- numeric(length(bp))
  for(i in seq_along(BStrees)){
     tmp <- BStrees[[i]]
     bptmp <- prop.part(tmp)
     bptmp <- SHORTwise(bptmp)[-1]
     ind <- fmatch(bp, bptmp)
     res[!is.na(ind)] <- res[!is.na(ind)] + 1
     # cherries can be check outside
     ind <- which(is.na(ind) & not_cherry)
     for(j in ind) res[j] <- res[j] + Transfer_Index(bp[[j]], tmp$edge, l)
  }
  res <- res / length(BStrees) * 100
  tree$node.label <- c(NA_real_, res)
  tree
}

