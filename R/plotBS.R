support <- function(tree, trees, method="FBP", tol=2e-8, scale=TRUE){
  trees <- keep.tip(trees, tree$tip.label)
  method <- match.arg(method, c("FBP", "TBE", "MCC"), several.ok=TRUE)
  multi <- ifelse(length(method)>1, TRUE, FALSE)
  tip2root <- ifelse(method=="MCC", TRUE, FALSE)
  if(all(sapply(trees, \(x)!is.null(x$edge.length)))){
    trees <- di2multi(trees, tol=tol) # , tip2root=tip2root)
  }
  if(multi) X <- matrix(NA, Nnode(tree), length(method),
                        dimnames = list(NULL, method))
  if("MCC" %in% method){
    trees <- .uncompressTipLabel(trees) # check if needed
    if(any(!is.rooted(trees)))
      stop("All trees need to be rooted for method 'MCC'!")
    x <- prop.clades(tree, trees)
    x <- (x / length(trees))
    if(!scale) x <- x * 100
    if(multi) X[, "MCC"] <- x
  }
  if("FBP" %in% method){
    trees <- .uncompressTipLabel(trees) # check if needed
    if (any(is.rooted(trees))) trees <- unroot(trees)
    x <- prop.clades(tree, trees)
    x <- (x / length(trees))
    if(!scale) x <- x * 100
    if(multi) X[, "FBP"] <- x
   }
  if("TBE" %in% method){
    x <- transferBootstrap(tree, trees, FALSE, scale=scale)
    if(multi) X[, "TBE"] <- x
  }
  if(multi) return(X)
  x
}



#' Plotting trees with bootstrap values
#'
#' \code{plotBS} plots a phylogenetic tree with the bootstrap values assigned
#' to the (internal) edges. It can also used to assign bootstrap values to a
#' phylogenetic tree. \code{add_support} adds support values to a plot.
#'
#' The functions can either assign the classical Felsenstein’s bootstrap
#' proportions (FBP) (Felsenstein (1985), Hendy & Penny (1985))  or the
#' transfer bootstrap expectation (TBE) of Lemoine et al. (2018). Using the
#' option \code{type=="n"} just assigns the bootstrap values and return the tree
#' without plotting it.
#'
#'
#' @param tree The tree on which edges the bootstrap values are plotted.
#' @param trees a list of trees (object of class "multiPhylo").
#' @param type the type of tree to plot, one of "phylogram", "cladogram", "fan",
#' "unrooted", "radial" or "none". If type is "none" the tree is returned with
#' the bootstrap values assigned to the node labels.
#' @param method either "FBP" the classical bootstrap (default), "TBE"
#' (transfer bootstrap) or "MCC" for assigning clade credibilities. In case of
#' "MCC" all trees need to be rooted.
#' @param bs.col color of bootstrap support labels.
#' @param bs.adj one or two numeric values specifying the horizontal and
#' vertical justification of the bootstrap labels.
#' @param digits integer indicating the number of decimal places.
#' @param p only plot support values higher than this percentage number
#' (default is 0).
#' @param sep seperator between the different methods.
#' @param \dots further parameters used by \code{plot.phylo}.
#' @param frame a character string specifying the kind of frame to be printed
#' around the bootstrap values. This must be one of "none" (the default),
#' "rect" or "circle".
#' @param tol a numeric value giving the tolerance to consider a branch length
#' significantly greater than zero.
#' @param scale return ratio or percentage.
#' @return \code{plotBS} returns silently a tree, i.e. an object of class
#' \code{phylo} with the bootstrap values as node labels. The argument
#' \code{trees} is optional and if not supplied the labels supplied
#' in the \code{node.label} slot will be used.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso  \code{\link[ape]{plot.phylo}}, \code{\link{add_ci}},
#' \code{\link[ape]{nodelabels}},
#' \code{\link[ape]{prop.clades}}, \code{\link{maxCladeCred}},
#' \code{\link{transferBootstrap}}, \code{\link[ape]{consensus}},
#' \code{\link{consensusNet}}
#' @references Felsenstein J. (1985) Confidence limits on phylogenies. An
#' approach using the bootstrap. \emph{Evolution} \bold{39}, 783--791
#'
#' Lemoine, F., Entfellner, J. B. D., Wilkinson, E., Correia, D., Felipe, M. D.,
#' De Oliveira, T., & Gascuel, O. (2018). Renewing Felsenstein’s phylogenetic
#' bootstrap in the era of big data. \emph{Nature}, \bold{556(7702)}, 452--456.
#'
#' Penny D. and Hendy M.D. (1985) Testing methods evolutionary tree
#' construction. \emph{Cladistics} \bold{1}, 266--278
#'
#' Penny D. and Hendy M.D. (1986) Estimating the reliability of evolutionary
#' trees. \emph{Molecular Biology and Evolution} \bold{3}, 403--417
#' @examples
#' fdir <- system.file("extdata/trees", package = "phangorn")
#' # RAxML best-known tree with bipartition support (from previous analysis)
#' raxml.tree <- read.tree(file.path(fdir,"RAxML_bipartitions.woodmouse"))
#' # RAxML bootstrap trees (from previous analysis)
#' raxml.bootstrap <- read.tree(file.path(fdir,"RAxML_bootstrap.woodmouse"))
#' par(mfrow=c(1,2))
#' plotBS(raxml.tree,  raxml.bootstrap, "p")
#' plotBS(raxml.tree,  raxml.bootstrap, "p", "TBE")
#' @rdname plotBS
#' @export
plotBS <- function(tree, trees, type = "phylogram", method="FBP",
                   bs.col = "black", bs.adj = NULL, digits=3, p = 0,
                   frame = "none", tol=1e-6, sep = "/", ...) {
  type <- match.arg(type, c("phylogram", "cladogram", "fan", "unrooted",
                            "radial", "none"))
  if(inherits(tree, "pml")) tree <- tree$tree
  if(!inherits(tree, "phylo")) stop("tree must be of class phylo")
#  method <- match.arg(method, c("FBP", "TBE", "MCC"), several.ok=TRUE)
# wird in support gecheckt
  if (hasArg(trees)) {
    x <-support(tree, trees, method=method, tol=tol)
    x <- round(x, digits=digits)
    if(length(method)>1) x <- apply(x, 1, paste, collapse=sep)
    tree$node.label <- x
  }
  else {
    if (is.null(tree$node.label)) stop("You need to supply 'trees' or the tree needs support-values as node.label")
    x <- tree$node.label
  }
  if(type=="none") return( tree )
  plot(tree, type = type, ...)

  label <- c(rep(0, length(tree$tip.label)), x)
  ind <- get("last_plot.phylo", envir = .PlotPhyloEnv)$edge[ ,2 ]
  if (type == "phylogram" | type == "cladogram") {
    root <- getRoot(tree)
    label <- c(rep(0, length(tree$tip.label)), x)
    label[root] <- 0
    ind <- which(label > p)
    if (is.null(bs.adj)) {
      bs.adj <- c(1, 1)
    }
    if (length(ind) > 0) {
      if(is.numeric(label)) label <- round(label, digits = digits)
      nodelabels(
        text = label[ind], node = ind,
        frame = frame, col = bs.col, adj = bs.adj, ...
      )
    }
  }
  else {
    if (is.null(bs.adj)) {
      bs.adj <- c(0.5, 0.5)
    }
    ind2 <- which(label[ind] > p)
    if (length(ind2 > 0)) {
      if(is.numeric(label)) label <- round(label, digits = digits)
      edgelabels(label[ind][ind2], ind2,
                 frame = frame,
                 col = bs.col, adj = bs.adj, ...
      )
    }
  }
  invisible(tree)
}


#' @rdname plotBS
#' @export
add_support <- function(tree, trees, method="FBP", tol=1e-8,
                        scale=TRUE, frame="none", digits=3, sep="/", ...){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
#  if(!all.equal(lastPP$edge, reorder(tree)$edge))
#     stop("tree seems to differ from the last plot!")
  x <- support(tree, trees, method=method, tol=tol, scale=scale)
  x <- round(x, digits=digits)
  if(length(method)>1) x <- apply(x, 1, paste, collapse=sep)
  drawSupportOnEdges(x, frame=frame, ...)
}



