nnin <- function(tree, n) {
  attr(tree, "order") <- NULL
  tree1 <- tree
  tree2 <- tree
  edge <- matrix(tree$edge, ncol = 2)
  parent <- edge[, 1]
  child <- tree$edge[, 2]
  k <- min(parent) - 1
  ind <- which(child > k)[n]
  if (is.na(ind)) return(NULL)
  p1 <- parent[ind]
  p2 <- child[ind]
  ind1 <- which(parent == p1)
  ind1 <- ind1[ind1 != ind][1]
  ind2 <- which(parent == p2)
  e1 <- child[ind1]
  e2 <- child[ind2[1]]
  e3 <- child[ind2[2]]
  tree1$edge[ind1, 2] <- e2
  tree1$edge[ind2[1], 2] <- e1
  tree2$edge[ind1, 2] <- e3
  tree2$edge[ind2[2], 2] <- e1
  if (!is.null(tree$edge.length)) {
    tree1$edge.length[c(ind1, ind2[1])] <- tree$edge.length[c(ind2[1], ind1)]
    tree2$edge.length[c(ind1, ind2[2])] <- tree$edge.length[c(ind2[2], ind1)]
  }
  tree1 <- reorder(tree1, "postorder")
  tree2 <- reorder(tree2, "postorder")
  result <- list(tree1, tree2)
  result
}

# faster for moves >> 1
n_nni <- function(tree, moves = 1L) {
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  nb.tip  <- length(tree$tip.label)

  if (nb.tip == 1) return(tree)

  pvector <- integer(max(edge)) # parents
  pvector[child] <- parent

  ch <- Children(tree)

  edges <- child[child %in% parent]  # these are the edges we sample from

  for (i in 1:moves) {
    p2 <- sample(edges, 1)
    p1 <- pvector[p2]
    ind1 <- ch[[p1]]
    v1 <- ind1[ind1 != p2][1]
    ind2 <- ch[[p2]]
    r2 <- sample(2, 1)
    v2 <- ind2[r2]
    ind1[ind1 == v1] <- v2
    ind2[r2] <- v1
    pvector[v1] <- p2
    pvector[v2] <- p1
    ch[[p1]] <- ind1
    ch[[p2]] <- ind2
  }
  tree$edge[, 1] <- pvector[child]
  attr(tree, "order") <- NULL
  reorder(tree, "postorder")
}

# roughly 2* fast from Martin Smith
one_nnin <- function(tree, n) {
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  lengths <- tree$edge.length
  nb.tip  <- length(tree$tip.label)
  ind     <- which(child > nb.tip)[n]
  if (is.na(ind)) return(NULL)
  nb.node <- tree$Nnode
  if (nb.node == 1) return(tree)
  p1      <- parent[ind]
  p2      <- child[ind]
  ind1    <- which(parent == p1)
  ind1    <- ind1[ind1 != ind][1]
  ind2    <- which(parent == p2)[sample(2, 1)]
  new_ind <- c(ind2, ind1)
  old_ind <- c(ind1, ind2)
  child_swap <- child[new_ind]
  edge [old_ind, 2L] <- child_swap
  child[old_ind] <- child_swap

  neworder <- reorderRcpp(edge, as.integer(nb.tip), as.integer(nb.tip + 1L), 2L)
  tree$edge <- edge[neworder, ]
  if (!is.null(tree$edge.length)) {
    lengths[old_ind] <- lengths[new_ind]
    tree$edge.length <- lengths[neworder]
  }
  attr(tree, "order") <- "postorder"
  tree
}


## @aliases nni rNNI rSPR
#' Tree rearrangements.
#'
#' \code{nni} returns a list of all trees which are one nearest neighbor
#' interchange away. \code{rNNI} and \code{rSPR} are two methods which simulate
#' random trees which are a specified number of rearrangement apart from the
#' input tree. Both methods assume that the input tree is bifurcating. These
#' methods may be useful in simulation studies.
#'
#'
#' @param tree A phylogenetic \code{tree}, object of class \code{phylo}.
#' @param moves Number of tree rearrangements to be transformed on a tree.  Can
#' be a vector
#' @param n Number of trees to be simulated.
#' @param k If defined just SPR of distance k are performed.
#' @return an object of class multiPhylo.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{allTrees}}, \code{\link{SPR.dist}}
#' @keywords cluster
#' @examples
#'
#' tree <- rtree(20, rooted = FALSE)
#' trees1 <- nni(tree)
#' trees2 <- rSPR(tree, 2, 10)
#'
#' @rdname nni
#' @export nni
nni <- function(tree) {
  tip.label <- tree$tip.label
  attr(tree, "order") <- NULL
  k <- min(tree$edge[, 1]) - 1
  n <- sum(tree$edge[, 2] > k)
  result <- vector("list", 2 * n)
  l <- 1
  for (i in 1:n) {
    tmp <- nnin(tree, i)
    tmp[[1]]$tip.label <- tmp[[2]]$tip.label <- NULL
    result[c(l, l + 1)] <- tmp
    l <- l + 2
  }
  attr(result, "TipLabel") <- tip.label
  class(result) <- "multiPhylo"
  result
}


#' @rdname nni
#' @export
rNNI <- function(tree, moves = 1, n = length(moves)) {
  k <- length(na.omit(match(tree$edge[, 2], tree$edge[, 1])))

  k_nni <- function(tree, ch, pvector, moves = 1L) {
    if (nb.tip == 1) return(tree)
    for (i in seq_len(moves)) {
      p2 <- sample(edges, 1)
      p1 <- pvector[p2]
      ind1 <- ch[[p1]]
      v1 <- ind1[ind1 != p2][1]
      ind2 <- ch[[p2]]
      r2 <- sample(2, 1)
      v2 <- ind2[r2]
      ind1[ind1 == v1] <- v2
      ind2[r2] <- v1
      pvector[v1] <- p2
      pvector[v2] <- p1
      ch[[p1]] <- ind1
      ch[[p2]] <- ind2
    }
    edge[, 1] <- pvector[child]
    neworder <- reorderRcpp(edge, nb.tip, nb.tip + 1L, 2L)
    tree$edge <- edge[neworder, ]
    if (!is.null(tree$edge.length)) {
      tree$edge.length <- tree$edge.length[neworder]
    }
    attr(tree, "order") <- "postorder"
    tree
  }

  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  nb.tip  <- as.integer(length(tree$tip.label))
  pvector <- integer(max(edge)) # parents
  pvector[child] <- parent
  ch <- Children(tree)
  edges <- child[child %in% parent]

  if (n == 1) {
    trees <- tree
    if (moves > 0) {
      trees <- k_nni(tree, ch, pvector, moves = moves)
    }
    trees$tip.label <- tree$tip.label
  }
  else {
    trees <- vector("list", n)
    if (length(moves) == 1) moves <- rep(moves, n)
    for (j in seq_len(n)) {
      tmp <- tree
      if (moves[j] > 0) {
        tmp <- k_nni(tree, ch, pvector, moves = moves[j])
      }
      tmp$tip.label <- NULL
      trees[[j]] <- tmp
    }
    attr(trees, "TipLabel") <- tree$tip.label
    class(trees) <- "multiPhylo"
  }
  trees
}


rNNI_old <- function(tree, moves = 1, n = length(moves)) {
  k <- length(na.omit(match(tree$edge[, 2], tree$edge[, 1])))

  if (n == 1) {
    trees <- tree
    if (moves > 0) {
      if (moves > 1) trees <- n_nni(trees, moves)
      else {
        trees <- one_nnin(trees, sample(k, 1))
      }
    }
    trees$tip.label <- tree$tip.label
  }
  else {
    trees <- vector("list", n)
    if (length(moves) == 1) moves <- rep(moves, n)
    for (j in 1:n) {
      tmp <- tree
      if (moves[j] > 0) {
        if (moves[j] > 1) tmp <- n_nni(tmp, moves[j])
        else {
          tmp <- one_nnin(tmp, sample(k, 1))
        }
      }
      tmp$tip.label <- NULL
      trees[[j]] <- tmp
    }
    attr(trees, "TipLabel") <- tree$tip.label
    class(trees) <- "multiPhylo"
  }
  trees
}


################################################################################
# SPR
################################################################################
dn <- function(x) {
  #  if (!is.binary(x) ) x <- multi2di(x, random = FALSE)
  if (is.null(x$edge.length)) x$edge.length <- rep(1, nrow(x$edge))
  else x$edge.length[] <- 1
  dist.nodes(x)
}


#' @rdname nni
#' @export
rSPR <- function(tree, moves = 1, n = length(moves), k = NULL) {
  if (n == 1) {
    trees <- tree
    for (i in 1:moves) trees <- kSPR(trees, k = k)
  }
  else {
    trees <- vector("list", n)
    if (length(moves) == 1) moves <- rep(moves, n)

    for (j in 1:n) {
      tmp <- tree
      if (moves[j] > 0) {
        for (i in 1:moves[j]) tmp <- kSPR(tmp, k = k)
      }
      tmp$tip.label <- NULL
      trees[[j]] <- tmp
    }
    attr(trees, "TipLabel") <- tree$tip.label
    class(trees) <- "multiPhylo"
  }
  trees
}


kSPR <- function(tree, k = NULL) {
  l <- length(tree$tip.label)
  root <- getRoot(tree)
  distN <- dn(tree)[-c(1:l), -c(1:l)]
  distN[upper.tri(distN)] <- Inf
  dN <- distN[lower.tri(distN)]
  tab <- tabulate(dN)
  tab[1] <- tab[1] * 2
  tab[-1] <- tab[-1] * 8
  if (is.null(k)) k <- seq_along(tab)
  k <- na.omit( (seq_along(tab))[k])
  if (length(k) > 1) k <- sample(seq_along(tab)[k], 1,
                                 prob = tab[k] / sum(tab[k]))
  if (k == 1) return(rNNI(tree, 1, 1))
  index <- which(distN == k, arr.ind = TRUE) + l
  m <- dim(index)[1]
  if (m == 0) stop("k is chosen too big")
  ind <- index[sample(m, 1), ]
  s1 <- sample(c(1, 2), 1)
  if (s1 == 1) res <- oneOf4(tree, ind[1], ind[2], sample(c(1, 2), 1),
      sample(c(1, 2), 1), root)
  if (s1 == 2) res <- oneOf4(tree, ind[2], ind[1], sample(c(1, 2), 1),
      sample(c(1, 2), 1), root)
#  res <- reroot2(res, root)
#  reorderPruning(res)
#  reroot(res, root, FALSE)
  res
}


oneOf4 <- function(tree, ind1, ind2, from = 1, to = 1, root) {
  if (!is.binary(tree))
    stop("Sorry, trees must be binary!")
#  tree <- reroot2(tree, ind2)
#  browser()#
  tree <- reroot(tree, ind2, FALSE)
  kids1 <- Children(tree, ind1)
  anc <- Ancestors(tree, ind1, "all")
  l <- length(anc)
  kids2 <- Children(tree, ind2)
  kids2 <- kids2[kids2 != anc[l - 1]]

  child <- tree$edge[, 2]
  tmp <- numeric(max(tree$edge))
  tmp[child] <- seq_along(child)

  edge <- tree$edge
  edge[tmp[kids1[-from]], 1] <- Ancestors(tree, ind1, "parent")
  edge[tmp[kids2[to]], 1] <- ind1
  edge[tmp[ind1]] <- ind2

  tree$edge <- edge
  attr(tree, "order") <- NULL
  tree <- reroot(tree, root, FALSE)
#  reorderPruning(tree)
  #  reorder(tree)
#  reorder(tree, "postorder")
  tree
}


# faster than kSPR
rSPR_Old <- function(tree, moves = 1, n = 1) {
  k <- length(tree$edge[, 1])
  if (n == 1) {
    trees <- tree
    for (i in 1:moves) trees <- sprMove(trees, sample(k, 1))
  }
  else {
    trees <- vector("list", n)
    for (j in 1:n) {
      tmp <- tree
      for (i in 1:moves) tmp <- sprMove(tmp, sample(k, 1))
      tmp$tip.label <- NULL
      trees[[j]] <- tmp
    }
    attr(trees, "TipLabel") <- tree$tip.label
    class(trees) <- "multiPhylo"
  }
  trees
}


sprMove <- function(tree, m) {
  if (is.rooted(tree)) tree <- unroot(tree)
  if (!is.binary(tree)) stop("Sorry, trees must be binary!")

  changeEdge <- function(tree, new, old) {
    tree$edge[tree$edge == old] <- 0L
    tree$edge[tree$edge == new] <- old
    tree$edge[tree$edge == 0L] <- new
    # needed for unrooted trees
    tree <- collapse.singles(tree)
    tree
  }

  edge <- tree$edge
  k <- max(edge)
  nTips <- length(tree$tip.label)
  nEdges <- 2 * nTips - 3
  if (m > nEdges) stop("m is chosen too big")

  parent <- edge[, 1]
  child <- edge[, 2]
  pv <- integer(k)
  pv[child] <- parent
  cv <- list()
  for (i in unique(parent)) cv[[i]] <- child[parent == i]
  bp <- bip(tree)
  root <- parent[!match(parent, child, 0)][1]

  ch <- child[m]
  pa <- parent[m]

  candidates <- !logical(k)
  candidates[root] <- FALSE
  candidates[cv[[ch]]] <- FALSE
  candidates[cv[[pa]]] <- FALSE
  candidates[pv[pa]] <- FALSE
  candidates[pa] <- FALSE

  ind <- which(candidates)
  l <- sample(ind, 1)

  cr <- FALSE

  if (!any(is.na(match(bp[[l]], bp[[ch]])))) {

    newroot <- cv[[ch]] # [ 1]
    newroot <- newroot[newroot > nTips][1]
#    tree <- reroot2(tree, newroot)
    tree <- reroot(tree, newroot, FALSE)
    edge <- tree$edge
    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    pv <- integer(k)
    pv[child] <- parent
    cv <- list()
    for (i in unique(parent)) cv[[i]] <- child[parent == i]

    tmp <- pa
    pa <- ch
    ch <- tmp
    cr <- TRUE
  }

  if (pa == root) {
    cp <- cv[[pa]]
    newroot <- cp[cp != ch]

    newroot <- newroot[newroot > nTips][1]
#    tree <- reroot2(tree, newroot)
    tree <- reroot(tree, newroot, FALSE)
    edge <- tree$edge
    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    pv <- integer(k)
    pv[child] <- parent
    cv <- list()
    for (i in unique(parent)) cv[[i]] <- child[parent == i]

    cr <- TRUE
  }

  el <- tree$edge.length
  cp <- cv[[pa]]
  sib <- cp[cp != ch]

  edge[child == l, 1] <- pa
  edge[child == pa, 1] <- pv[l]
  edge[child == sib, 1] <- pv[pa]

  el[child == sib] <- el[child == sib] + el[child == pa]
  el[child == l] <- el[child == l] / 2
  el[child == pa] <- el[child == l]

  tree$edge <- edge
  tree$edge.length <- el
  if (cr) tree <- changeEdge(tree, root, newroot)
  tree <- reorder(tree, "postorder")
  tree
}
