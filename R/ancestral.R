#' Ancestral character reconstruction.
#'
#' Marginal reconstruction of the ancestral character states.
#'
#' The argument "type" defines the criterion to assign the internal nodes. For
#' \code{ancestral.pml} so far "ml and marginal (empirical) "bayes" and for
#' \code{ancestral.pars} "MPR" and "ACCTRAN" are possible.
#'
#' The function return a list containing the tree with node labels, the original
#' alignment as an \code{phyDat} object, a data.frame containing the
#' probabilities belonging to a state for all (internal nodes) and the most
#' likely state. For parsimony and nucleotide data the most likely state might
#' be ambiguous. For ML this is very unlikely to be the case.
#'
#' If the input tree does not contain unique node labels the function
#' \code{ape::MakeNodeLabel} is used to create them.
#'
#' With parsimony reconstruction one has to keep in mind that there will be
#' often no unique solution.
#'
#' The functions use the node labels of the provided tree (also if part of the
#' \code{pml} object) if these are unique. Otherwise the function
#' \code{ape::MakeNodeLabel} is used to create them.
#'
#' For further details see vignette("Ancestral").
#'
#' @param object an object of class pml
#' @param tree a tree, i.e. an object of class pml
#' @param data an object of class phyDat
#' @param type method used to assign characters to internal nodes, see details.
#' @param cost A cost matrix for the transitions between two states.
#' @param return return a \code{phyDat} object or matrix of probabilities.
##  @param x an object of class ancestral.
#' @param \dots Further arguments passed to or from other methods.
#' @return An object of class ancestral. This is a list containing the tree with
#' node labels, the original alignment as an \code{phyDat} object, a
#' \code{data.frame} containing the probabilities belonging to a state for all
#' (internal nodes) and the most likely state.
## For \code{return="phyDat"} an object  of class "phyDat", containing
## the ancestral states of all nodes. For nucleotide data this can contain
## ambiguous states. Apart from fitch parsimony the most likely states are
## returned.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}}, \code{\link{parsimony}}, \code{\link[ape]{ace}},
#' \code{\link{plotAnc}}, \code{\link{anc_heatmap}},
#' \code{\link{latag2n.phyDat}},
#' \code{\link[ape]{latag2n}}, \code{\link{gap_as_state}},
#' \code{\link[ape]{root}}, \code{\link[ape]{makeNodeLabel}}
#' @references
#' \bibshow{*, Felsenstein2004, Yang2006}
#'
#' Swofford, D.L., Maddison, W.P. (1987) Reconstructing ancestral character
#' states under Wagner parsimony. \emph{Math. Biosci.} \bold{87}: 199--229
#' @keywords cluster
#' @examples
#'
#' example(NJ)
#' # generate node labels to ensure plotting will work
#' tree <- makeNodeLabel(tree)
#' fit <- pml(tree, Laurasiatherian)
#' anc.ml <- anc_pml(fit)
#' anc.p <- anc_pars(tree, Laurasiatherian)
#' # plot ancestral sequences at the root
#' plotSeqLogo( anc.ml, 48, 1, 20)
#' plotSeqLogo( anc.p, 48, 1, 20)
#' # plot the first character
#' plotAnc(anc.ml)
#' # plot the third character
#' plotAnc(anc.ml, 3)
#' # plot joint reconstruction as heatmap
#' anc_heatmap(anc.ml, select=1:50)
#'
#' @rdname ancestral.pml
#' @export
ancestral.pml <- function(object, type = "marginal", return = "ancestral", ...) {
  assert_pml(object)
  pt <- match.arg(type, c("marginal", "ml", "bayes", "lhood")) # "joint",
  rt <- match.arg(return, c("prob", "phyDat", "ancestral"))
  tree <- object$tree
  INV <- object$INV
  inv <- object$inv
  data <- object$data
  data <- data[tree$tip.label]
  data_type <- attr(data, "type")
  attrib <- attributes(data)
  pos <- match(attrib$levels, attrib$allLevels)
  if (is.null(attr(tree, "order")) || attr(tree, "order") != "postorder") {
    tree <- reorder(tree, "postorder")
  }
  w <- object$w
  g <- object$g
  k <- length(w)
  rate <- object$rate
  eig <- object$eig
  bf <- object$bf
  assert_phylo(tree, has_edge_length=TRUE)
  nTips <- as.integer(length(tree$tip.label))
  tree <- reorder(tree, "postorder")
  if (any(tree$edge.length < 0)) tree <- minEdge(tree)
  ll.0 <- as.matrix(INV %*% (bf * inv))
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  node_label <- makeAncNodeLabel(tree, ...)
  tree$node.label <- node_label
  joint <- TRUE
  if(length(w) > 1 || object$inv > 0) joint <- FALSE
  on.exit(.Call("ll_free2"))
  .Call("ll_init2", nr, nTips, nc, as.integer(k))
  tmp <- pml.fit4(tree, data, bf, k = k, levels = attr(data, "levels"),
                  inv = inv, rate = rate, g = g, w = w,
                 eig = eig, INV = INV, ll.0 = ll.0, ...)
  ll <- .Call("get_ll", nr, nTips, nc, as.integer(k))
  scm <- .Call("get_scm", nr, nTips, as.integer(k))
  dim(ll) <- c(nr, nc, nTips, k)
  dim(scm) <- c(nr, nTips, k)

  l <- length(w)
  m <- length(tree$edge[,1]) + 1 # max(edge)
  dat <- vector(mode = "list", length = m * l)
  nNode <- Nnode(tree)
  result <- vector(mode = "list", length = nNode)
  result2 <- vector(mode = "list", length = nNode)
  dim(dat) <- c(l, m)

  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  nTips <- min(parent) - 1
  # in C with scaling
  r <- getRoot(tree)
  el <- tree$edge.length
  P <- getP(el, eig, g)

  for(i in 1:l) dat[[i, r]] <- ll[,,r-nTips,i]
  for (i in 1:l) {
    for (j in (m - 1):1) {
      if (child[j] > nTips) {
        tmp2 <- (dat[[i, parent[j]]] / (ll[ , , child[j]-nTips, i] %*% P[[i, j]]))
        dat[[i, child[j]]] <- (tmp2 %*% P[[i, j]]) * ll[ , , child[j]-nTips, i]
      }
    }
  }
  SCALE_EPS <- 1.0/4294967296.0
  SCM <- scm[,1, ]
  SCM <- matrix(SCM, nrow=dim(scm)[1])
  sc_min <- apply(SCM,1,min)
  SCM <- SCM - sc_min

  for (j in unique(parent)) {
    tmp <- matrix(0, nr, nc)
    if (inv > 0) tmp <- as.matrix(INV) * inv

    for (i in 1:l) {
      tmp2 <- dat[[i, j]] * (SCALE_EPS ** SCM[,i])
      tmp <- tmp + w[i] * tmp2
    }
    if ((pt == "bayes") || (pt == "marginal")) tmp <- tmp * rep(bf, each = nr)
    #KBH: don't normalize if partial lhoods desired
    # TODO scale lhood
    if (pt != "lhood") tmp <- tmp / rowSums(tmp)
    if (data_type == "DNA") {
      tmp_max <- p2dna(tmp)
      tmp_max  <- fitchCoding2ambiguous(tmp_max)
    }
    else {
      tmp_max  <- pos[max.col(tmp)]
    }
    result[[j - nTips]] <- tmp
    result2[[j - nTips]] <- tmp_max
  }
  attrib$names <- node_label
  attributes(result2) <- attrib
  if(rt == "phyDat") return(rbind(data, result2))
  if(rt == "prob") {
    tmp <- c(unclass(new2old.phyDat(data)), result)
    attrib$names <- c(tree$tip.label, node_label)
    attributes(tmp) <- attrib
    return(tmp)
  }
  attributes(result) <- attrib
  if(joint) result2 <- joint_pml(object)
  erg <- list(tree=tree, data=data, prob=result, state=result2)
  class(erg) <- "ancestral"
  erg
}


#' @rdname ancestral.pml
#' @export
anc_pml <- function(object, type = "marginal", ...) {
  call <- c(object$call, match.call())
  res <- ancestral.pml(object, type=type, return="ancestral")
  res$call <- call
  res
}


fitchCoding2ambiguous <- function(x, type = "DNA") {
  y <- c(1L, 2L, 4L, 8L, 8L, 3L, 5L, 9L, 6L, 10L, 12L, 7L, 11L, 13L,
    14L, 15L, 15L, 15L)
  fmatch(x, y)
}



#' @rdname ancestral.pml
#' @export
ancestral.pars <- function(tree, data, type = c("MPR", "ACCTRAN", "POSTORDER"),
                           cost = NULL, return = "prob", ...) {
  call <- match.call()
  if (hasArg(tips)) tips <- list(...)$tips
  else tips <- TRUE
  type <- match.arg(type)
  assert_phylo(tree)
  assert_phyDat(data)
  tree$nodel.label <- makeAncNodeLabel(tree)
  if (type == "ACCTRAN" || type=="POSTORDER") {
    res <- ptree(tree, data, return = return, acctran=(type == "ACCTRAN"),
                 tips=tips)
    attr(res, "call") <- call
  }
  if (type == "MPR") {
    res <- mpr(tree, data, cost = cost, return = return, tips=tips)
    attr(res, "call") <- call
  }
  res
}


#' @rdname ancestral.pml
#' @export
anc_pars <- function(tree, data, type = c("MPR", "ACCTRAN", "POSTORDER"),
                           cost = NULL, ...) {
  #
  call <- match.call()
  type <- match.arg(type)
  assert_phylo(tree)
  assert_phyDat(data)
  tree$node.label <- makeAncNodeLabel(tree, ...)
  contrast <- attr(data, "contrast")
  data <- data[tree$tip.label,]

  prob <- ancestral.pars(tree, data, type, cost, tips=FALSE)
  joint <- joint_sankoff(tree, data, cost)
  ind <- identical_sites(data)
  if(length(ind)>0){
    for(k in seq_len(Nnode(tree))){
      prob[[k]][ind,] <- contrast[data[[1]][ind],]
      joint[[k]][ind] <- data[[1]][ind]
    }
  }
  erg <- list(tree=tree, data=data, prob=prob, state=joint, call=call)
  class(erg) <- "ancestral"
  erg
}



#' @rdname ancestral.pml
#' @export
pace <- ancestral.pars


mpr.help <- function(tree, data, cost = NULL) {
  tree <- reorder(tree, "postorder")
  levels <- attr(data, "levels")
  l <- length(levels)
  if (is.null(cost)) {
    cost <- matrix(1, l, l)
    cost <- cost - diag(l)
  }
  dat <- prepareDataSankoff(data)
  datp <- pnodes(tree, dat, cost)
  nr <- attr(data, "nr")
  nc <- attr(data, "nc")
  node <- as.integer(tree$edge[, 1] - 1L)
  edge <- as.integer(tree$edge[, 2] - 1L)
  res <- .Call('sankoffMPR', datp, as.numeric(cost), as.integer(nr),
    as.integer(nc), node, edge, as.integer(Nnode(tree)))
  root <- getRoot(tree)
  res[[root]] <- datp[[root]]
  res
}


mpr <- function(tree, data, cost = NULL, return="prob", tips=FALSE, ...) {
  data <- subset(data, tree$tip.label)
  att <- attributes(data)
  if(tips) att$names <- c(tree$tip.label, tree$node.label)
  else att$names <- tree$node.label
  type <- att$type
  nr <- att$nr
  nc <- att$nc
  res <- mpr.help(tree, data, cost)
  l <- length(tree$tip.label)
  m <- length(res)
  nNode <- Nnode(tree)
  ntips <- length(tree$tip.label)
  contrast <- att$contrast
  eps <- 5e-6
  rm <- apply(res[[ntips + 1]], 1, min)
  RM <- matrix(rm, nr, nc) + eps

  fun <- function(X) {
    rs <- rowSums(X) # apply(X, 1, sum)
    X / rs
  }
#  for (i in 1:ntips) res[[i]] <- contrast[data[[i]], , drop = FALSE]
  for (i in (ntips + 1):m) res[[i]][] <- as.numeric(res[[i]] < RM)
  if(tips)for(i in seq_len(ntips)) res[[i]] <- contrast[data[[i]],,drop=FALSE]
  else  res <- res[(ntips + 1):m]
  res_prob <- lapply(res, fun)
  attributes(res_prob) <- att
  if(return=="prob") return(res_prob)
  attributes(res) <- att
  res_state <- highest_state(res)
  attributes(res_state) <- att
  if(tips) res_state[seq_len(ntips)] <- data
  res_state # list(res_prob, res_state)
}


ptree <- function(tree, data, acctran=TRUE, return = "prob", tips=FALSE, ...) {
  tree <- reorder(tree, "postorder")
  data <- subset(data, tree$tip.label)
  edge <- tree$edge
  att <- attributes(data)
  att$names <- c(tree$tip.label, tree$node.label)
  nr <- att$nr
  type <- att$type
  m <- max(edge)
  nNode <- Nnode(tree)
  nTip <- Ntip(tree)
  f <- init_fitch(data, FALSE, FALSE, m=2L)
  f$traverse(edge)
  tmp <- reorder(tree)$edge
  tmp <- tmp[tmp[,2]>Ntip(tree),]
  if(length(tmp)>0 && acctran)f$acctran_traverse(tmp)
  res <- res_state <- vector("list", nNode)
  res <- vector("list", m)
  att$names <- c(tree$tip.label, tree$node.label) #makeAncNodeLabel(tree, ...)
  fun <- function(X) {
    rs <- rowSums(X)
    X / rs
  }
  contrast <- att$contrast
  for(i in seq_len(nTip)) res[[i]] <- contrast[data[[i]], , drop=FALSE]
  for(i in seq_len(nNode)) {
    res[[i+nTip]] <- f$getAnc(i+nTip)[seq_len(nr), , drop=FALSE]
  }
  res <- lapply(res, fun)
  attributes(res) <- att
  if(!tips) res <- res[tree$node.label]

  if(return=="prob") return(res)

  res_state <- highest_state(res)
  attributes(res_state) <-  attributes(res)
  res_state
}


makeAncNodeLabel <- function(tree, ...){
  if(!is.null(tree$node.label)){
    node_label <- tree$node.label
    if(length(unique(node_label)) == Nnode(tree)) return(node_label)
    else message("Node labels are not unique, used makeNodeLabel(tree, ...) to create them!")
  }
  tree <- makeNodeLabel(tree, ...)
  tree$node.label
}


identical_sites <- function(x){
  res <- rep(TRUE, attr(x, "nr"))
  for(i in seq_along(x)) res <- res & (x[[i]] == x[[1]])
  which(res)
}


#' @srrstats {G1.0} in the lines folloing: 48
#' @srrstats {G2.3, G2.3a} in lines: 78, 79, 211, 235
NULL
