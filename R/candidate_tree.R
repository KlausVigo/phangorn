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
    tree <- pratchet(x, maxit=10L, trace=0, perturbation = "ratchet")
    tree <- multi2di(tree)
    tree <- unroot(tree)
    tree <- acctran(tree, x)
    tree$edge.length <- tree$edge.length / sum(attr(x, "weight"))
  }
  if(method=="tipdated"){
    if(is.null(tip.dates)) stop("Argument tip.dates is missing!")
    if(is.null(names(tip.dates))) names(tip.dates) <- names(x)
    dm <- dist.ml(x, ...)
    tree <- fastme.ols(dm)
    tree <- rtt(tree, tip.dates[tree$tip.label])
    tree <- nnls.tree(dm, tree, method = "tipdated",
                      tip.dates=tip.dates[tree$tip.label])
    tree$tip.dates <- tip.dates[tree$tip.label]
  }
  minEdge(tree, tau=eps, enforce_ultrametric=enforce_ultrametric)
}
