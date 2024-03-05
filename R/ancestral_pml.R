#' Ancestral character reconstruction.
#'
#' Marginal reconstruction of the ancestral character states.
#'
#' The argument "type" defines the criterion to assign the internal nodes. For
#' \code{ancestral.pml} so far "ml" and (empirical) "bayes" and for
#' \code{ancestral.pars} "MPR" and "ACCTRAN" are possible.
#'
#' With parsimony reconstruction one has to keep in mind that there will be
#' often no unique solution.
#'
#' The functions use node labels if these are present and unique. Otherwise the
#' function \code{ape::MakeNodeLabel} is used to create them.
#'
#' For further details see vignette("Ancestral").
#'
#' @param object an object of class pml
#' @param tree a tree, i.e. an object of class pml
#' @param data an object of class phyDat
#' @param type method used to assign characters to internal nodes, see details.
#' @param cost A cost matrix for the transitions between two states.
#' @param return return a \code{phyDat} object or matrix of probabilities.
#' @param x an object of class ancestral.
#' @param \dots Further arguments passed to or from other methods.
#' @return An object of class ancestral containing theestimated character
#' states.
#' For \code{return="phyDat"} an object  of class "phyDat", containing
#' the ancestral states of all nodes. For nucleotide data this can contain
#' ambiguous states. Apart from fitch parsimony the most likely states are
#' returned.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}}, \code{\link{parsimony}}, \code{\link[ape]{ace}},
#' \code{\link{plotAnc}}, \code{\link[ape]{root}},
#' \code{\link[ape]{makeNodeLabel}}
#' @references Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer
#' Associates, Sunderland.
#'
#' Swofford, D.L., Maddison, W.P. (1987) Reconstructing ancestral character
#' states under Wagner parsimony. \emph{Math. Biosci.} \bold{87}: 199--229
#'
#' Yang, Z. (2006). \emph{Computational Molecular evolution}. Oxford University
#' Press, Oxford.
#' @keywords cluster
#' @importFrom fastmatch fmatch
#' @examples
#'
#' example(NJ)
#' # generate node labels to ensure plotting will work
#' tree <- makeNodeLabel(tree)
#' fit <- pml(tree, Laurasiatherian)
#' anc.ml <- ancestral.pml(fit, type = "ml")
#' anc.p <- ancestral.pars(tree, Laurasiatherian)
#' \dontrun{
#' require(seqLogo)
#' seqLogo( t(subset(anc.ml, 48, 1:20)[[1]]), ic.scale=FALSE)
#' seqLogo( t(subset(anc.p, 48, 1:20)[[1]]), ic.scale=FALSE)
#' }
#' # plot the first site pattern
#' plotAnc(tree, anc.ml, 1)
#' # plot the third character
#' plotAnc(tree, anc.ml, attr(anc.ml, "index")[3])
#'
#' @rdname ancestral.pml
#' @export
ancestral.pml <- function(object, type = "marginal", return = "prob", ...) {
  call <- match.call()
  pt <- match.arg(type, c("marginal", "joint", "ml", "bayes"))
  tree <- object$tree
  INV <- object$INV
  inv <- object$inv
  data <- getCols(object$data, tree$tip.label)
  data_type <- attr(data, "type")
  if (is.null(attr(tree, "order")) || attr(tree, "order") != "postorder") {
    tree <- reorder(tree, "postorder")
  }
  nTips <- length(tree$tip.label)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  m <- length(edge) + 1 # max(edge)
  w <- object$w
  g <- object$g
  l <- length(w)
  nr <- attr(data, "nr")
  nc <- attr(data, "nc")
  dat <- vector(mode = "list", length = m * l)
  result <- vector(mode = "list", length = m)
  dim(dat) <- c(l, m)

  x <- attributes(data)
  label <- makeAncNodeLabel(tree, ...)
  x[["names"]] <- label
  tmp <- length(data)

  if (return != "phyDat") {
    result <- new2old.phyDat(data)
  } else {
    result[1:nTips] <- data
  }
  eig <- object$eig

  bf <- object$bf
  el <- tree$edge.length
  P <- getP(el, eig, g)
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  node <- as.integer(node - min(node))
  edge <- as.integer(edge - 1)
  nTips <- as.integer(length(tree$tip.label))
  mNodes <- as.integer(max(node) + 1)
  contrast <- attr(data, "contrast")
  # proper format
  eps <- 1.0e-5
  attr <- attributes(data)
  pos <- match(attr$levels, attr$allLevels)
  nco <- as.integer(dim(contrast)[1])
  for (i in 1:l) dat[i, (nTips + 1):m] <- .Call('LogLik2', data, P[i, ], nr, nc,
      node, edge, nTips, mNodes, contrast, nco)
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  nTips <- min(parent) - 1
  # in C with scaling
  for (i in 1:l) {
    for (j in (m - 1):1) {
      if (child[j] > nTips) {
        tmp2 <- (dat[[i, parent[j]]] / (dat[[i, child[j]]] %*% P[[i, j]]))
        dat[[i, child[j]]] <- (tmp2 %*% P[[i, j]]) * dat[[i, child[j]]]
      }
    }
  }
  for (j in unique(parent)) {
    tmp <- matrix(0, nr, nc)
    if (inv > 0) tmp <- as.matrix(INV) * inv
    for (i in 1:l) {
      # scaling!!!
      tmp <- tmp + w[i] * dat[[i, j]]
    }
    if ((pt == "bayes") || (pt == "marginal")) tmp <- tmp * rep(bf, each = nr)
    tmp <- tmp / rowSums(tmp)

    if (return == "phyDat") {
      if (data_type == "DNA") {
        tmp <- p2dna(tmp)
        tmp <- fitchCoding2ambiguous(tmp)
      }
      else {
        tmp <- pos[max.col(tmp)]
      }
    }
    result[[j]] <- tmp
  }
  attributes(result) <- x
  attr(result, "call") <- call
  if(return=="prob") class(result) <- c("ancestral", "phyDat")
  result
}



#' @rdname ancestral.pml
#' @export
as.phyDat.ancestral <- function(x, ...) {
  type <- attr(x, "type")
  fun2 <- function(x) {
    x <- p2dna(x)
    fitchCoding2ambiguous(x)
  }
  if (type == "DNA") {
    res <- lapply(x, fun2)
  }
  else {
    eps <- 1.0e-5
    contr <- attr(x, "contrast")
    attr <- attributes(x)
    pos <- match(attr$levels, attr$allLevels)
    res <- lapply(x, function(x, pos) pos[max.col(x)], pos)
  }
  attributes(res) <- attributes(x)
  class(res) <- "phyDat"
  return(res)
}


#' @rdname ancestral.pml
#' @export
as.data.frame.ancestral <- function(x, ...) {
  stopifnot(inherits(x, "ancestral"))
  l <- length(x)
  nr <- attr(x, "nr")
  nc <- attr(x, "nc")
  index <- attr(x, "index")
  nr <- length(index)
  nam <- names(x)
  X <- matrix(0, l*length(index), nc)
  j <- 0
  for(i in seq_len(l)){
    X[(j+1):(j+nr), ] <- x[[i]][index, ]
    j <- j + nr
  }
  res <- data.frame(Site=rep(seq_len(nr), l), Node=rep(nam, each=nr), X)
  colnames(res) <- c("Site", "Node", attr(x, "levels"))
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
  type <- match.arg(type)
  if (type == "ACCTRAN" || type=="POSTORDER") {
    res <- ptree(tree, data, return = return, acctran=(type == "ACCTRAN"))
    attr(res, "call") <- call
  }
  if (type == "MPR") {
    res <- mpr(tree, data, cost = cost, return = return)
    attr(res, "call") <- call
  }
  res
}


#' @rdname ancestral.pml
#' @export
pace <- ancestral.pars


mpr.help <- function(tree, data, cost = NULL) {
  tree <- reorder(tree, "postorder")
  if (!inherits(data, "phyDat")) {
    stop("data must be of class phyDat")
  }
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


mpr <- function(tree, data, cost = NULL, return = "prob", ...) {
  data <- subset(data, tree$tip.label)
  att <- attributes(data)
  type <- att$type
  nr <- att$nr
  nc <- att$nc
  res <- mpr.help(tree, data, cost)
  l <- length(tree$tip.label)
  m <- length(res)
  label <- makeAncNodeLabel(tree, ...)
  att[["names"]] <- label
  ntips <- length(tree$tip.label)
  contrast <- att$contrast
  eps <- 5e-6
  rm <- apply(res[[ntips + 1]], 1, min)
  RM <- matrix(rm, nr, nc) + eps

  fun <- function(X) {
    rs <- rowSums(X) # apply(X, 1, sum)
    X / rs
  }
  for (i in 1:ntips) res[[i]] <- contrast[data[[i]], , drop = FALSE]
  for (i in (ntips + 1):m) res[[i]][] <- as.numeric(res[[i]] < RM)
  if (return == "prob") {
    #        for(i in 1:ntips) res[[i]] <- contrast[data[[i]],,drop=FALSE]
    res <- lapply(res, fun)
    attributes(res) <- att
    class(res) <- c("ancestral", "phyDat")
  }
  #    else res[1:ntips] <- data[1:ntips]
  fun2 <- function(x) {
    x <- p2dna(x)
    fitchCoding2ambiguous(x)
  }
  if (return != "prob") {
    if (type == "DNA") {
      res <- lapply(res, fun2)
      attributes(res) <- att
    }
    else {
      attributes(res) <- att
      res <- as.phyDat.ancestral(res)
    }
    res[1:ntips] <- data
  }
  res
}


#
# ACCTRAN
#
acctran2 <- function(tree, data) {
  if(!is.binary(tree)) tree <- multi2di(tree)
  tree <- reorder(tree, "postorder")
  edge <- tree$edge
  data <- subset(data, tree$tip.label)
  f <- init_fitch(data, FALSE, FALSE, m=2L)
  psc_node <- f$pscore_node(edge)
  tmp <- reorder(tree)$edge
  tmp <- tmp[tmp[,2]>Ntip(tree), ,drop=FALSE]
  f$traverse(edge)
  if(length(tmp)>0)f$acctran_traverse(tmp)
  psc <- f$pscore_acctran(edge)
  el <- psc
  parent <- unique(edge[,1])
  desc <- Descendants(tree, parent, "children")
  for(i in seq_along(parent)){
    x <- psc_node[parent[i]] -sum(psc[desc[[i]]])
    if(x>0) el[desc[[i]] ] <- el[desc[[i]] ] + x/length(desc[[i]])
  }
  tree$edge.length <- el[edge[,2]]
  tree
}


#' @rdname parsimony
#' @export
acctran <- function(tree, data) {
  if (inherits(tree, "multiPhylo")) {
    compress <- FALSE
    if (!is.null(attr(tree, "TipLabel"))){
      compress <- TRUE
      tree <- .uncompressTipLabel(tree)
    }
    res <- lapply(tree, acctran2, data)
    class(res) <- "multiPhylo"
    if (compress) res <- .compressTipLabel(res)
    return(res)
  }
  acctran2(tree, data)
}


ptree <- function(tree, data, return = "prob", acctran=TRUE, ...) {
  tree <- reorder(tree, "postorder")
  data <- subset(data, tree$tip.label)
  edge <- tree$edge
  att <- attributes(data)
  nr <- att$nr
  type <- att$type
  m <- max(edge)
  nTip <- Ntip(tree)
  f <- init_fitch(data, FALSE, FALSE, m=2L)
  f$traverse(edge)
  tmp <- reorder(tree)$edge
  tmp <- tmp[tmp[,2]>Ntip(tree),]
  if(length(tmp)>0 && acctran==TRUE)f$acctran_traverse(tmp)
  res <- vector("list", m)
  att$names <- makeAncNodeLabel(tree, ...)
  if(type=="DNA" && return != "prob"){
    indx <- c(1, 2, 6, 3, 7, 9, 12, 4, 8, 10, 13, 11, 14, 15, 16)
    res[1:nTip] <- data[1:nTip]
    for(i in (nTip+1):m)
      res[[i]] <- indx[f$getAncAmb(i)[1:nr]]
    attributes(res) <- att
    return(res)
  }
  else {
    fun <- function(X) {
      rs <- rowSums(X)
      X / rs
    }
    contrast <- att$contrast
    for(i in seq_len(nTip)) res[[i]] <- contrast[data[[i]], , drop=FALSE]
    for(i in (nTip+1):m) res[[i]] <- f$getAnc(i)[1:nr, , drop=FALSE]
    res <- lapply(res, fun)
    attributes(res) <- att
    class(res) <- c("ancestral", "phyDat")
  }
  if(return != "prob"){
    res <- as.phyDat.ancestral(res)
    class(res) <- "phyDat"
  }
  res
}


makeAncNodeLabel <- function(tree, ...){
  if(!is.null(tree$node.label)){
    node_label <- tree$node.label
    if(length(unique(node_label)) == Nnode(tree)) return(c(tree$tip.label, node_label))
    else message("Node labels are not unique, used makeNodeLabel(tree, ...) to create them!")
  }
  tree <- makeNodeLabel(tree, ...)
  c(tree$tip.label, tree$node.label)
}


identical_sites <- function(x){
  res <- rep(TRUE, attr(x, "nr"))
  for(i in seq_along(x)) res <- res & (x[[i]] == x[[1]])
  which(res)
}

