minEdge <- function(tree, tau=1e-8, enforce_ultrametric=FALSE){
  if(tau<0) stop("tau must be >= 0!")
  if(any(tree$edge.length < tau) || enforce_ultrametric){
    rooted <- is.rooted(tree)
    if(enforce_ultrametric) rooted <- TRUE
    if(rooted){
      nTip <- Ntip(tree)
      ind <- seq_len(nTip)
      nh <- nodeHeight(tree)[ind]
      if(enforce_ultrametric) nh <- rep(0, nTip)
    }
    tree$edge.length[tree$edge.length < tau] <- tau
    if(rooted){
      el <- numeric(max(tree$edge))
      el[tree$edge[,2]] <- tree$edge.length
      nh2 <- nodeHeight(tree)[ind]
      el[ind] <- el[ind] + (nh2 - nh)
      if(any(el < tau)) el[ind] <- el[ind] + tau - min(el)
      tree$edge.length <- el[tree$edge[,2]]
    }
  }
  tree
}


#' @rdname phangorn-internal
#' @importFrom stats cor
#' @export
candidate_tree <- function(x, method=c("unrooted", "ultrametric", "tipdated"),
                           eps = 1e-8, tip.dates=NULL, ...){
  method <- match.arg(method, c("unrooted", "ultrametric", "tipdated"))
  enforce_ultrametric <- FALSE
  if(method=="ultrametric"){
    dm <- dist.ml(x, ...)
    tree <- wpgma(dm)
    enforce_ultrametric <- TRUE
  }
  if(method=="unrooted"){
    tree <- pratchet(x, maxit=5L, trace=0, perturbation = "stochastic",
                     all=FALSE)
    tree <- multi2di(tree)
    tree <- unroot(tree)
    tree <- acctran(tree, x)
    tree$edge.length <- tree$edge.length / sum(attr(x, "weight"))
  }
  if(method=="tipdated"){
    if(is.null(tip.dates)) stop("Argument tip.dates is missing!")
    if(is.null(names(tip.dates))) names(tip.dates) <- names(x)
    dm <- dist.ml(x, ...)
    tree <- supgma(dm, tip.dates)
  }
  minEdge(tree, tau=eps, enforce_ultrametric=enforce_ultrametric)
}


# like is.ultrametric
check_tip_dates <- function(tree, tip.dates){
  tip.dates <- tip.dates[tree$tip.label]
  nh <- node.depth.edgelength(tree)[seq_along(tree$tip.label)]
  isTRUE(all.equal(cor(tip.dates,nh), 1))
}


proper_tree <- function(x, tree, method=c("ultrametric", "tipdated"),
                        tip.dates=NULL, eps = 1e-8, ...){
  method <- match.arg(method, c("ultrametric", "tipdated"))
  if(!is.null(tree$edge.length)){
    if(method=="ultrametric" && is.ultrametric(tree)) return(tree)
    if(method=="tipdated" && check_tip_dates(tree, tip.dates)) return(tree)
  }
  dm <- dist.ml(x, ...)
  nnls.tree(dm, tree, method)
}
