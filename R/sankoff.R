prepareDataSankoffNew <- function(data) {
  contrast <- attr(data, "contrast")
  contrast[contrast == 0] <- 1.0e+06
  contrast[contrast == 1] <- 0.0
  attr(data, "contrast") <- contrast
  data
}


sankoffNew <- function(tree, data, cost = NULL, site = "pscore") {
  if (!inherits(data, "phyDat"))
    stop("data must be of class phyDat")
  data <- prepareDataSankoffNew(data)
  levels <- attr(data, "levels")
  l <- length(levels)

  if (is.null(cost)) {
    cost <- matrix(1, l, l)
    cost <- cost - diag(l)
  }

  l <- length(data)
  if (inherits(tree, "phylo")) return(fit.sankoffNew(tree, data, cost,
      returnData = site))
  if (inherits(tree, "multiPhylo")) {
    if (is.null(tree$TipLabel)) tree <- unclass(tree)
    return(sapply(tree, fit.sankoffNew, data, cost, site))
  }
}


fit.sankoffNew <- function(tree, data, cost, returnData = c("pscore", "site",
                                                            "data")) {
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorder(tree, "postorder")
  returnData <- match.arg(returnData)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  weight <- attr(data, "weight")
  nr <- attr(data, "nr")

  contr <- attr(data, "contrast")

  q <- length(tree$tip.label)
  nc <- attr(data, "nc")
  m <- length(edge) + 1
  dat <- vector(mode = "list", length = m)
  dat[1:q] <- data[tree$tip.label]
  node <- as.integer(node - 1)
  edge <- as.integer(edge - 1)
  nTips <- as.integer(length(tree$tip.label))
  mNodes <- as.integer(max(node) + 1)
  res <- .Call("sankoff3B", dat, as.numeric(cost), as.integer(nr),
    as.integer(nc), node, edge, mNodes, nTips, as.double(contr),
    as.integer(nrow(contr)), PACKAGE = "phangorn")
  root <- getRoot(tree)
  erg <- .Call("C_rowMin", res[[root]], as.integer(nr), as.integer(nc),
    PACKAGE = "phangorn")
  if (returnData == "site") return(erg)
  pscore <- sum(weight * erg)
  result <- pscore
  if (returnData == "data") {
    res[1:nTips] <- new2old.phyDat(data)[tree$tip.label]
    result <- list(pscore = pscore, dat = res)
  }
  result
}

