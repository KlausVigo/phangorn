tree2phyDat <- function(trees) {
  # some minor error checking
  if (!inherits(trees, "multiPhylo"))
    stop("trees must be object of class 'multiPhylo.'")

  labels <- lapply(trees, function(x) sort(x$tip.label))
  ulabels <- unique(labels)
  lul <- length(ulabels)
  # compute matrix representation phylogenies
  X <- vector("list", lul) # list of bipartitions
  characters <- 0 # number of characters
  weights <- NULL
  species <- trees[[1]]$tip.label

  characters <- 0 # number of characters
  weights <- NULL

  for (i in 1:lul) {
    pos <- match(labels, ulabels[i])
    ind <- which(!is.na(pos))
    temp <- prop.part(trees[ind]) # find all bipartitions
    # create matrix representation of trees[[i]] in X[[i]]
    TMP <- matrix(1L, nrow = length(temp) - 1,
      ncol = length(trees[[ind[1]]]$tip.label))
    for (j in seq_len(nrow(TMP)))  TMP[j, c(temp[[j + 1]])] <- 2L
    colnames(TMP) <- attr(temp, "labels") # label rows

    X[[i]] <- TMP

    species <- union(species, trees[[ind[1]]]$tip.label) # accumulate labels
    characters <- characters + nrow(TMP) # count characters
    weights <- c(weights, attr(temp, "number")[-1])
  }
  data <- matrix(data = 3L, nrow = characters, ncol = length(species),
    dimnames = list(NULL, species))
  j <- 1
  for (i in seq_along(X)) {
    # copy each of X into supermatrix data
    data[c(j:((j - 1) + nrow(X[[i]]))), colnames(X[[i]])] <- X[[i]]
    # [1:nrow(X[[i]]),1:ncol(X[[i]])]
    j <- j + nrow(X[[i]])
  }
  data <- as.data.frame(data)
  # compute contrast matrix
  contrast <- matrix(data = c(1, 0, 0, 1, 1, 1), 3, 2,
    dimnames = list(NULL, c("0", "1")), byrow = TRUE)

  attr(data, "row.names") <- NULL
  class(data) <- "phyDat"
  attr(data, "weight") <- weights
  attr(data, "nr") <- length(weights)
  attr(data, "nc") <- 2L
  attr(data, "levels") <- c("0", "1")
  attr(data, "allLevels") <- c("0", "1", "?")
  attr(data, "type") <- "USER"
  attr(data, "contrast") <- contrast
  class(data) <- "phyDat"
  data
}


my.supertree <- function(trees, trace = 0, ...) {
  has_edge_length <- vapply(trees, function(x) !is.null(x$edge.length), FALSE)
  if (all(has_edge_length)){
    trees <- .uncompressTipLabel(trees)
    trees <- di2multi(trees)
  }
  XX <- tree2phyDat(trees)
  supertree <- pratchet(XX, trace = trace, ...)
  # maybe add edge length supertree <- acctran(supertree, XX)
  return(supertree)
}


# Robinson-Foulds supertree
fun.rf <- function(x, tree) sum(RF.dist(x, tree))
fun.spr <- function(x, tree) sum(SPR.dist(x, tree))

dist.superTree <- function(tree, trace = 0, fun, start = NULL,
                           multicore = FALSE, mc.cores = NULL) {
  if (multicore && is.null(mc.cores)) {
    mc.cores <- detectCores()
  }
  if (is.null(start)) start <- superTree(tree, rooted = FALSE)
  if (inherits(start, "multiPhylo")) start <- start[[1]]
  best_tree <- unroot(start)
  best <- fun(best_tree, tree)
  if (trace > 0) cat("best score so far:", best, "\n")
  eps <- TRUE
  while (eps) {
    nni_trees <- nni(best_tree)
    if (multicore) {
      tmp <- mclapply(nni_trees, fun, tree, mc.cores = mc.cores)
      tmp <- unlist(tmp)
    }
    else tmp <- sapply(nni_trees, fun, tree)
    if (min(tmp) < best) {
      ind <- which.min(tmp)
      best_tree <- nni_trees[[ind]]
      best <- tmp[ind]
      if (trace > 0) cat("best score so far:", best, "\n")
    }
    else eps <- FALSE
  }
  attr(best_tree, "score") <- best
  best_tree
}



#' Super Tree methods
#'
#' These function \code{superTree} allows the estimation of a supertree from a
#' set of trees using either Matrix representation parsimony, Robinson-Foulds
#' or SPR as criterion.
#'
#' The function \code{superTree} extends the function mrp.supertree from Liam
#' Revells, with artificial adding an outgroup on the root of the trees.  This
#' allows to root the supertree afterwards. The functions is internally used in
#' DensiTree. The implementation for the RF- and SPR-supertree are very basic
#' so far and assume that all trees share the same set of taxa.
#'
#' @param tree an object of class \code{multiPhylo}
#' @param method An argument defining which algorithm is used to optimize the
#' tree.  Possible are "MRP", "RF", and "SPR".
#' @param rooted should the resulting supertrees be rooted.
#' @param trace defines how much information is printed during optimization.
#' @param start a starting tree can be supplied.
#' @param multicore logical, whether models should estimated in parallel.
#' @param mc.cores The number of cores to use, i.e. at most how many child
#' processes will be run simultaneously.
#' @param \dots further arguments passed to or from other methods.
#' @return The function returns an object of class \code{phylo}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com} Liam Revell
#' @seealso \code{mrp.supertree}, \code{\link{densiTree}},
#' \code{\link{RF.dist}}, \code{\link{SPR.dist}}
#' @references Baum, B. R., (1992) Combining trees as a way of combining data
#' sets for phylogenetic inference, and the desirability of combining gene
#' trees. \emph{Taxon}, \bold{41}, 3-10.
#'
#' Ragan, M. A. (1992) Phylogenetic inference based on matrix representation of
#' trees. \emph{Molecular Phylogenetics and Evolution}, \bold{1}, 53-58.
#' @keywords cluster
#' @examples
#'
#' data(Laurasiatherian)
#' set.seed(1)
#' bs <- bootstrap.phyDat(Laurasiatherian,
#'                        FUN = function(x) upgma(dist.hamming(x)), bs=50)
#'
#' mrp_st <- superTree(bs)
#' plot(mrp_st)
#' \dontrun{
#' rf_st <- superTree(bs, method = "RF")
#' spr_st <- superTree(bs, method = "SPR")
#' }
#'
#' @export superTree
superTree <- function(tree, method = "MRP", rooted = FALSE, trace = 0,
                      start = NULL, multicore = FALSE, mc.cores = NULL, ...) {
  fun <- function(x) {
    x <- reorder(x, "postorder")
    nTips <- length(x$tip.label)
    x$edge[x$edge > nTips] <- x$edge[x$edge > nTips] + 2L
    l <- nrow(x$edge)
    oldroot <- x$edge[l, 1L]
    x$edge <- rbind(x$edge, matrix(c(rep(nTips + 2, 2), oldroot, nTips + 1),
                                   2L, 2L))
    x$edge.length <- c(x$edge.length, 100, 100)
    x$tip.label <- c(x$tip.label, "ZZZ")
    x$Nnode <- x$Nnode + 1L
    x
  }

#  labels_start <- unique(unlist(lapply(tree, function(x)x$tip.label)))
#  TODO check for missing labels
  tmp <- Nnode(tree)
  if(any(tmp < (3 + !rooted))){
    tree_tmp <- tree[ tmp > (2 + !rooted) ]
    if (length(tree_tmp) == 0) return(consensus(tree))
    tree <- tree_tmp
  }


  if (method != "MRP") rooted <- FALSE
  if (!rooted) tree <- unroot(tree)
  if (method == "MRP" | is.null(start)) {
    if (rooted) {
      if (!is.null(attr(tree, "TipLabel"))) tree <- .uncompressTipLabel(tree)
      tree <- unclass(tree)
      if (rooted) tree <- lapply(tree, fun)
      class(tree) <- "multiPhylo"
    }
    res <- my.supertree(tree, trace = trace, ...)
    if (rooted) {
      if (inherits(res, "multiPhylo")) {
        res <- lapply(res, root, "ZZZ")
        res <- lapply(res, drop.tip, "ZZZ")
        class(res) <- "multiPhylo"
      }
      else {
        res <- root(res, "ZZZ")
        res <- drop.tip(res, "ZZZ")
      }
    }
    if (inherits(res, "multiPhylo")) {
      fun <- function(x) {
        x$edge.length <- rep(.1, nrow(x$edge))
        x
      }
      res <- lapply(res, fun)
      res <- lapply(res, reorder, "postorder")
      class(res) <- "multiPhylo"
    }
    else {
      res$edge.length <- rep(.1, nrow(res$edge))
      res <- reorder(res, "postorder")
    }
  }
  if (method == "MRP") return(res)
  if (is.null(start)) start <- res
  tree <- unroot(tree)
  tree <- reorder(tree, "postorder")

  if (method == "RF") res <- dist.superTree(tree, trace = trace, fun.rf,
      start = start, multicore = multicore, mc.cores = mc.cores)
  if (method == "SPR") res <- dist.superTree(tree, trace = trace, fun.spr,
      start = start, multicore = multicore, mc.cores = mc.cores)
  res
}
