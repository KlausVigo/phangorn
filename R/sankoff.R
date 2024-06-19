prepareDataSankoff <- function(data) {
  contrast <- attr(data, "contrast")
  contrast[contrast == 0] <- 1.0e+06
  contrast[contrast == 1] <- 0.0
  attr(data, "contrast") <- contrast
  data
}


fit.sankoff <- function(tree, data, cost,
                        returnData = c("pscore", "site", "data")) {
  tree <- reorder(tree, "postorder")
  returnData <- match.arg(returnData)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  weight <- attr(data, "weight")
  nr <- attr(data, "nr")

  contr <- attr(data, "contrast")

  q <- length(tree$tip.label)
  nc <- attr(data, "nc")
  m <- length(edge) + 1L
  dat <- vector(mode = "list", length = m)
  dat[1:q] <- subset(data, tree$tip.label)
  node <- as.integer(node - 1L)
  edge <- as.integer(edge - 1L)
  nTips <- as.integer(length(tree$tip.label))
  mNodes <- as.integer(max(node) + 1)
  res <- .Call('sankoff_c', dat, as.numeric(cost), as.integer(nr),
    as.integer(nc), node, edge, mNodes, nTips, as.double(contr),
    as.integer(nrow(contr)))
  root <- getRoot(tree)
  erg <- .Call('C_rowMin', res[[root]], as.integer(nr), as.integer(nc))
  if (returnData == "site") return(erg)
  pscore <- sum(weight * erg)
  result <- pscore
  if (returnData == "data") {
    res[1:nTips] <- new2old.phyDat(data)[tree$tip.label]
    result <- list(pscore = pscore, dat = res)
  }
  result
}


#' @rdname parsimony
#' @export
sankoff <- function(tree, data, cost = NULL, site = "pscore") {
  if (!inherits(data, "phyDat")) stop("data must be of class phyDat")
  data <- prepareDataSankoff(data)
  if (is.null(cost)) {
    levels <- attr(data, "levels")
    l <- length(levels)
    cost <- matrix(1, l, l)
    cost <- cost - diag(l)
  }
  if (inherits(tree, "phylo")) return(fit.sankoff(tree, data, cost,
                                                  returnData = site))
  if (inherits(tree, "multiPhylo")) {
    if (is.null(tree$TipLabel)) tree <- unclass(tree)
    return(sapply(tree, fit.sankoff, data, cost, site))
  }
}


# traverse the tree
# for NNI or ancestral reconstruction external=FALSE is enough
# offset
pnodes2 <- function(tree, rooted=is.rooted(tree), offset=Nnode(tree),
                    external=FALSE){
  if(!rooted) tree <- unroot(tree)
  offset <- as.integer(offset)
  tree <- reorder(tree, "postorder")
  tree_2 <- reorder(tree)
  sibs <- Siblings(tree)
  anc <- Ancestors(tree, type = "parent")
  root <- tree_2$edge[1,1]
  nodes <- unique(tree_2$edge[,1])[-1]
  #  me <- max(tree$edge)
  #  edge2 <- matrix
  e2 <- vector("list", length(nodes))
  for(i in seq_along(nodes)){
    ni <- nodes[i]
    if(anc[ni]==root)e2[[i]] <- sibs[[ni]]
    else e2[[i]] <- c(sibs[[ni]], anc[ni] + offset)
  }
  cbind(rep(nodes+offset, lengths(e2)), unlist(e2))
}


pnodes <- function(tree, data, cost) {
  tree <- reorder(tree, "postorder")
  weight <- attr(data, "weight")
  nr <- attr(data, "nr")
  nc <- attr(data, "nc")
  contr <- attr(data, "contrast")
  nTips <- as.integer(length(tree$tip.label))
  m <- as.integer( Ntip(tree) + 2 * Nnode(tree) )
  data <- subset(data, tree$tip.label)
  edge_2 <- pnodes2(tree)
  EDGE <- rbind(tree$edge, edge_2)
  mNodes <- m
  res <- .Call('sankoff_c', data, as.numeric(cost),
                as.integer(nr), as.integer(nc), as.integer(EDGE[,1]-1L),
                as.integer(EDGE[,2]-1L), m, nTips,
                as.double(contr), as.integer(nrow(contr)))
  res[seq_len(Ntip(tree))] <- new2old.phyDat(data)
  res
}


sankoff_nni <- function(tree, data, cost, ...) {
  if(isSymmetric(cost)) tree <- unroot(tree)
  tree <- reorder(tree, "postorder")
  isRooted <- !(isSymmetric(cost))
  INDEX <-  indexNNI_fitch(tree, offset=Nnode(tree))
  if (!inherits(data, "phyDat")) stop("data must be of class phyDat")
  data <- data[tree$tip.label]
  levels <- attr(data, "levels")
  l <- length(levels)
  weight <- attr(data, "weight")
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  contr <- attr(data, "contrast")
  ntip <- as.integer(Ntip(tree))
  nnode <- as.integer(Nnode(tree))
  p0 <- fit.sankoff(tree, data, cost, returnData = "pscore")
  dat <- pnodes(tree, data, cost)
  pscore <- .Call('sankoff_nni_c', dat, nr, cost, nc, as.double(weight),
                  INDEX[,1:4] - 1L, nrow(INDEX), ntip, as.double(contr),
                  as.integer(nrow(contr)))
  INDEX <- rbind(INDEX[, c(1, 3, 2, 4, 5, 6)], INDEX[, c(2, 3, 1, 4, 5, 6)])
  swap <- 0

  pscore <- as.vector(pscore)
  candidates <- which(pscore < p0)

  while (length(candidates)>0) {
    ind <- which.min(pscore[candidates])
    tree2 <- changeEdge(tree, INDEX[candidates[ind], c(2, 3)])
    test <- fit.sankoff(tree2, data, cost, returnData = "pscore")
    if (test < p0) {
      p0 <- test
      swap <- swap + 1
      tree <- tree2
      indi <- which(INDEX[, 5] %in% INDEX[candidates[ind],])
      candidates <- setdiff(candidates, indi)
    }
    else candidates <- candidates[-ind]
  }
  list(tree = tree, pscore = p0, swap = swap)
}


optim.sankoff <- function(tree, data, cost = NULL, trace = 1, ...) {
  if (!inherits(tree, "phylo")) stop("tree must be of class phylo")
  if (is.rooted(tree)) tree <- unroot(tree)
  tree <- reorder(tree, "postorder")
  if (!inherits(data, "phyDat")) stop("data must be of class phyDat")
  data <- data[tree$tip.label]
  addTaxa <- FALSE
  mapping <- map_duplicates(data)
  if (!is.null(mapping)) {
    addTaxa <- TRUE
    tree2 <- drop.tip(tree, mapping[, 1])
    tree2 <- unroot(tree2)
    tree <- reorder(tree2, "postorder")
  }

  rt <- FALSE
  dat <- prepareDataSankoff(data)
  l <- attr(dat, "nc")
  if (is.null(cost)) {
    cost <- matrix(1, l, l)
    cost <- cost - diag(l)
  }

  tree$edge.length <- NULL
  swap <- 0
  iter <- TRUE
  pscore <- fit.sankoff(tree, dat, cost, returnData = "pscore")

  on.exit({
    if (rt) tree <- acctran(tree, data)
    if (addTaxa) {
      if (rt) tree <- add.tips(tree, tips = mapping[, 1], where = mapping[, 2],
                               edge.length = rep(0, nrow(mapping)))
      else tree <- add.tips(tree, tips = mapping[, 1], where = mapping[, 2])
    }
    attr(tree, "pscore") <- pscore
    return(tree)
  })
  while (iter) {
    res <- sankoff_nni(tree, dat, cost, ...)
    tree <- res$tree
    if (trace > 1) cat("optimize topology: ", pscore, "-->", res$pscore, "\n")
    pscore <- res$pscore
    swap <- swap + res$swap
    if (res$swap == 0) iter <- FALSE
  }
  if (trace > 0) cat("Final p-score", pscore, "after ", swap,
                     "nni operations \n")
}



