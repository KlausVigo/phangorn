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
#' \code{\link{plotAnc}}, \code{\link{latag2n.phyDat}}, \code{\link[ape]{latag2n}},
#' \code{\link{gap_as_state}}, \code{\link[ape]{root}},
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
#' anc.ml <- anc_pml(fit)
#' anc.p <- anc_pars(tree, Laurasiatherian)
#' # plot ancestral sequences at the root
#' plotSeqLogo( anc.ml, 48, 1, 20)
#' plotSeqLogo( anc.p, 48, 1, 20)
#' # plot the first character
#' plotAnc(anc.ml)
#' # plot the third character
#' plotAnc(anc.ml, 3)
#'
#' @rdname ancestral.pml
#' @export
ancestral.pml <- function(object, type = "marginal", return = "prob", ...) {
  call <- match.call()
  pt <- match.arg(type, c("marginal", "ml", "bayes")) # "joint",
  rt <- match.arg(return, c("prob", "phyDat", "ancestral"))
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
  nNode <- Nnode(tree)
  m <- length(edge) + 1 # max(edge)
  w <- object$w
  g <- object$g
  l <- length(w)
  nr <- attr(data, "nr")
  nc <- attr(data, "nc")
  dat <- vector(mode = "list", length = m * l)
  result <- vector(mode = "list", length = nNode)
  result2 <- vector(mode = "list", length = nNode)
  dim(dat) <- c(l, m)
  node_label <- makeAncNodeLabel(tree, ...)
  tree$node.label <- node_label
  tmp <- length(data)

  joint <- TRUE
  if(length(w) > 1 || object$inv > 0) joint <- FALSE

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
  attrib <- attributes(data)
  pos <- match(attrib$levels, attrib$allLevels)
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
  ind <- identical_sites(data)
  if(length(ind)>0){
    for(k in seq_len(nNode)){
      result[[k]][ind,] <- contrast[data[[1]][ind],]
      result2[[k]][ind] <- data[[1]][ind]
    }
  }
#  browser()
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
  ancestral.pml(object, type=type, return="ancestral")
}



#' Export and convenience functions for  ancestral reconstructions
#'
#' \code{write.ancestral} allows to export ancestral reconstructions. It writes
#' out the tree, a tab delimited text file with the probabilities and the
#' alignment. \code{ancestral} generates an object of class ancestral.
#'
#' This allows also to read in reconstruction  made by iqtree to use the
#' plotting capabilities of R.
#' @param x an object of class ancestral.
#' @param file a file name. File endings are added.
#' @param ... Further arguments passed to or from other methods.
#' @returns \code{write.ancestral}  returns the input x invisibly.
#' @seealso \code{\link{ancestral.pml}}, \code{\link{plotAnc}}
#' @examples
#' data(Laurasiatherian)
#' fit <- pml_bb(Laurasiatherian[,1:100], "JC", rearrangement = "none")
#' anc_ml <- anc_pml(fit)
#' write.ancestral(anc_ml)
#' # Can be also results from iqtree
#' align <- read.phyDat("ancestral_align.fasta", format="fasta")
#' tree <- read.tree("ancestral_tree.nwk")
#' df <- read.table("ancestral.state", header=TRUE)
#' anc_ml_disc <- as.ancestral(tree, align, df)
#' plotAnc(anc_ml_disc, 20)
#' unlink(c("ancestral_align.fasta", "ancestral_tree.nwk", "ancestral.state"))
#' @rdname write.ancestral
#' @export
write.ancestral <- function(x, file="ancestral"){
  stopifnot(inherits(x, "ancestral"))
  write.phyDat(x$data, file=paste0(file, "_align.fasta"), format="fasta")
  tab <- as.data.frame(x)
  write.table(tab, file=paste0(file, ".state"), quote=FALSE, row.names=FALSE,
              sep="\t")
  write.tree(x$tree, file=paste0(file, "_tree.nwk"))
  invisible(x)
}


#' @param tree an object of class phylo.
#' @param align an object of class phyDat.
#' @param prob an data.frame containing a matrix of posterior probabilities for
#' each state and site.
#' @importFrom utils head
#' @rdname write.ancestral
#' @export
as.ancestral <- function(tree, align, prob){
  stopifnot(inherits(tree, "phylo"))
  stopifnot(inherits(align, "phyDat"))
  stopifnot(inherits(prob, "data.frame"))
  if(is.null(tree$node.label))stop("tree needs node.label")
  state <- extract_states(prob, attr(align, "type"),
                          levels=attr(align, "levels"))
  prob <- extract_prob(prob, align, tree$node.label)
  erg <- list(tree=tree, data=align[tree$tip.label], prob=prob,
              state=state[tree$node.label])
  class(erg) <- "ancestral"
  erg
}


extract_states <- function(x, type, levels=NULL){
  node_label <- unique(x$"Node")
  y <- matrix(x$"State", nrow=length(node_label), byrow=TRUE,
              dimnames = list(node_label, NULL))
  if(type %in% c("DNA", "AA")) return( phyDat(y, type=type) )
  phyDat(y, type=type, levels=levels)
}

extract_prob <- function(x, align, node_label){
  index <- attr(align, "index")
  idx <- which(!duplicated(index))
  att <- attributes(align)
  res <- vector("list", length(node_label))
  for(i in seq_along(node_label)){
     tmp <- x[x$"Node"==node_label[i], -c(1:3)] |> as.matrix()
     dimnames(tmp) <- NULL
     res[[i]] <- as.matrix(tmp)[idx,]
  }
  att$names <- node_label
  attributes(res) <- att
  res
}

#' @rdname write.ancestral
#' @export
print.ancestral <- function(x, ...){
  stopifnot(inherits(x, "ancestral"))
  print(x$tree)
  cat("\n")
  print(x$data)
#  cat("\n")
#  print(head(x$prob))
}

#' @rdname ancestral.pml
#' @param x an object of class ancestral
#' @export
as.phyDat.ancestral <- function(x, ...) {
  rbind(x$data, x$state)
}



highest_state <- function(x, ...) {
  type <- attr(x, "type")
  fun2 <- function(x) {
    x <- p2dna(x)
    fitchCoding2ambiguous(x)
  }
  if (type == "DNA" && !has_gap_state(x)) {
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
  y <- rbind(x$data, x$state)
  contrast <- attr(x$data, "contrast")
  Y <- unlist(as.data.frame(y))
  n_node <- Nnode(x$tree)
  n_tip <- Ntip(x$tree)
  l <- length(y)
  nr <- attr(x$data, "nr")
  nc <- attr(x$data, "nc")
  index <- attr(x$data, "index")
  nr <- length(index)
  nam <- names(y)
  X <- matrix(0, l*length(index), nc)
  j <- 0
  for(i in seq_len(n_tip)){
    X[(j+1):(j+nr), ] <- contrast[x$data[[i]][index], ]
    j <- j + nr
  }
  for(i in seq_len(n_node)){
    X[(j+1):(j+nr), ] <- x$prob[[i]][index, ]
    j <- j + nr
  }
  res <- data.frame(Node=rep(nam, each=nr), Site=rep(seq_len(nr), l), Y, X)
  colnames(res) <- c("Node", "Site", "State",
                     paste0("p_", attr(x$data, "levels")))
  rownames(res) <- NULL
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
  if (!inherits(data, "phyDat")) stop("data must be of class phyDat")
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
#    class(res_prob) <- c("ancestral", "phyDat")
  #    else res[1:ntips] <- data[1:ntips]
  #fun2 <- function(x) {
  #  x <- p2dna(x)
  #  fitchCoding2ambiguous(x)
  #}
#  if (return != "prob") {
#  if (type == "DNA") {
#    res_state <- lapply(res, fun2)
#    attributes(res_state) <- att
#  }
#  else {
  attributes(res) <- att
  res_state <- highest_state(res)
  attributes(res_state) <- att
  if(tips)res_state[seq_len(ntips)] <- data
#  }
#    res[1:ntips] <- data
#  }
  res_state # list(res_prob, res_state)
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

#, return = "prob"
ptree <- function(tree, data, acctran=TRUE, return = "prob", tips=FALSE, ...) {
  tree <- reorder(tree, "postorder")
  data <- subset(data, tree$tip.label)
  edge <- tree$edge
  att <- attributes(data)
  att$names <- c(tree$tip.label, tree$node.label)
  #else att$names <- tree$node.label
  nr <- att$nr
  type <- att$type
  m <- max(edge)
  nNode <- Nnode(tree)
  nTip <- Ntip(tree)
  f <- init_fitch(data, FALSE, FALSE, m=2L)
  f$traverse(edge)
  tmp <- reorder(tree)$edge
  tmp <- tmp[tmp[,2]>Ntip(tree),]
  if(length(tmp)>0 && acctran==TRUE)f$acctran_traverse(tmp)
  res <- res_state <- vector("list", nNode)
  res <- vector("list", m)
  att$names <- c(tree$tip.label, tree$node.label) #makeAncNodeLabel(tree, ...)
#  else {
    fun <- function(X) {
      rs <- rowSums(X)
      X / rs
    }
    contrast <- att$contrast
    for(i in seq_len(nTip)) res[[i]] <- contrast[data[[i]], , drop=FALSE]
#    for(i in (nTip+1):m) res[[i]] <- f$getAnc(i)[1:nr, , drop=FALSE]
    for(i in seq_len(nNode)) {
      res[[i+nTip]] <- f$getAnc(i+nTip)[seq_len(nr), , drop=FALSE]
    }
    res <- lapply(res, fun)
    attributes(res) <- att
#    class(res) <- c("ancestral", "phyDat")
#  }
  if(!tips) res <- res[tree$node.label]

  if(return=="prob"){
#    if(tips)
    return(res)
#    else return(res[tree$node.label])
  }

#  if(type=="DNA"){
#    indx <- c(1, 2, 6, 3, 7, 9, 12, 4, 8, 10, 13, 11, 14, 15, 16)
#    res_state[seq_len(nTip)] <- data
#    for(i in (nTip+1):m)
#      res_state[[i]] <- indx[f$getAncAmb(i+nTip)[1:nr]]
#    attributes(res_state) <- att
#    #    return(res)
 # } else{
    res_state <- highest_state(res)
    attributes(res_state) <-  attributes(res)
#    class(res_state) <- "phyDat"
#  }
#  if(tips) return(res_state)
#  res_state[tree$node.label]
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
