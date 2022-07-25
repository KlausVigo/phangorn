bip <- function(x) {
  x <- reorder(x, "postorder")
  nTips <- as.integer(length(x$tip.label))
  res <- .Call("_phangorn_bipCPP", x$edge, nTips)
  attr(res, "labels") <- x$tip.label
  res
}


bipart <- function(x) {
  x <- reorder(x, "postorder")
  nTips <- as.integer(length(x$tip.label))
  res <- .Call('_phangorn_bipartCPP', x$edge, nTips)
  attr(res, "labels") <- x$tip.label
  res
}

# needed??
#bipartition <- function(tree) {
#  if (is.rooted(tree)) tree <- unroot(tree)
#  tree <- reorder(tree, "postorder")
#  bp <- bip(tree)
#  nTips <- length(tree$tip.label)
#  l <- length(bp)
#  res <- matrix(0L, l, nTips)
#  for (i in 1:l) res[i, bp[[i]]] <- 1L
#  res <- res[tree$edge[, 2], , drop = FALSE]
#  colnames(res) <- tree$tip.label
#  rownames(res) <- tree$edge[, 2]
#  res[res[, 1] == 1, ] <- 1L - res[res[, 1] == 1, ]
#  res
#}
