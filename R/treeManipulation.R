#
# tree manipulation
#


# no checks for postorder
#' @rdname midpoint
#' @export
getRoot <- function(tree) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") ==
    "postorder") {
    return(tree$edge[nrow(tree$edge), 1])
  }
  z <- logical(max(tree$edge))
  z[tree$edge[, 1]] <-  TRUE
  z[tree$edge[, 2]] <-  FALSE
  z <- which(z)
  if (length(z) == 1) return(z)
  else stop("There are apparently two root edges in your tree")
}


reroot <- function (tree, node, switch_root=TRUE) {
  root <- getRoot(tree)
  if (node == root)
    return(reorder(tree, "postorder"))
  anc <- Ancestors(tree, node, "all")
  l <- length(anc)
  ind <- match(c(node, anc[-l]), tree$edge[, 2])
  tree$edge[ind, c(1, 2)] <- tree$edge[ind, c(2, 1)]
  nb.tip <- Ntip(tree)
  neworder <- reorderRcpp(tree$edge, as.integer(nb.tip), as.integer(node), 2L)
  tree$edge <- tree$edge[neworder, ]
  if(!is.null(tree$edge.length)) tree$edge.length <- tree$edge.length[neworder]
  if(switch_root){
    tree$edge[tree$edge == root] <-  0L
    tree$edge[tree$edge == node] <-  root
    tree$edge[tree$edge == 0L] <-  node
  }
  attr(tree, "order") <- "postorder"
  if(switch_root) tree <- collapse.singles(tree)
  tree
}


changeEdge <- function(tree, swap, edge = NULL, edge.length = NULL) {
  attr(tree, "order") <- NULL
  child <- tree$edge[, 2]
  tmp <- integer(max(child))
  tmp[child] <- seq_along(child)
  tree$edge[tmp[swap[1]], 2] <- as.integer(swap[2])
  tree$edge[tmp[swap[2]], 2] <- as.integer(swap[1])
  if (!is.null(edge)) {
    tree$edge.length[tmp[edge]] <- edge.length
  }
  reorder(tree, "postorder")
}


changeEdgeLength <- function(tree, edge, edge.length) {
  tree$edge.length[match(edge, tree$edge[, 2])] <- edge.length
  tree
}


## @aliases midpoint pruneTree getRoot
#' Tree manipulation
#'
#' \code{midpoint} performs midpoint rooting of a tree.  \code{pruneTree}
#' produces a consensus tree.
#' \code{pruneTree} prunes back a tree and produces a consensus tree, for trees
#' already containing nodelabels.  It assumes that nodelabels are numerical or
#' character that allows conversion to numerical, it uses
#' as.numeric(as.character(tree$node.labels)) to convert them.
#' \code{midpoint} by default assumes that node labels contain support values.
#' This works if support values are computed from splits, but should be
#' recomputed for clades.
#' \code{keep_as_tip} takes a list of tips and/or node labels and returns a tree
#' pruned to those. If node label, then it prunes all descendants of that node
#' until that internal node becomes a tip.

#' @param tree an object of class \code{phylo}.
#' @param FUN a function evaluated on the nodelabels, result must be logical.
#' @param node.labels are node labels 'support' values (edges), 'label' or
#' should labels get 'deleted'?
#' @param \dots further arguments, passed to other methods.
#' @return \code{pruneTree} and \code{midpoint} a tree. \code{getRoot} returns
#' the root node.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{consensus}}, \code{\link[ape]{root}},
#' \code{\link[ape]{multi2di}}
#' @keywords cluster
#' @examples
#'
#' tree <- rtree(10, rooted = FALSE)
#' tree$node.label <- c("", round(runif(tree$Nnode-1), digits=3))
#'
#' tree2 <- midpoint(tree)
#' tree3 <- pruneTree(tree, .5)
#'
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow = c(3,1))
#' plot(tree, show.node.label=TRUE)
#' plot(tree2, show.node.label=TRUE)
#' plot(tree3, show.node.label=TRUE)
#' par(old.par)
#'
#' @rdname midpoint
#' @export midpoint
midpoint <- function(tree, node.labels = "support", ...)
  UseMethod("midpoint")


#' @rdname midpoint
#' @method midpoint phylo
#' @export
midpoint.phylo <- function(tree, node.labels = "support", ...) {
  # distance from node to root
  node2root <- function(x) {
    x <- reorder(x, "postorder")
    el <- numeric(max(x$edge))
    parents <- x$edge[, 1]
    child <- x$edge[, 2]
    el[child] <- x$edge.length
    l <- length(parents)
    res <- numeric(max(x$edge))
    for (i in l:1) {
      res[child[i]] <- el[child[i]]  + res[parents[i]]
    }
    res
  }
  if (is.null(tree$edge.length)) {
    warning("tree needs edge length")
    return(tree)
  }
  oldtree <- tree
  if(Ntip(tree)==1) return(tree)
  if(Ntip(tree)==2){
    tree <- collapse.singles(tree)
    el <- sum(tree$edge.length)
    tree$edge.length[] <- el / 2
    return(tree)
  }
  tree <- unroot(tree)
  nTips <- length(tree$tip.label)
  maxD1 <- node2root(tree)[1:nTips]
  ind <- which.max(maxD1)
  tmproot <- Ancestors(tree, ind, "parent")
  nTips  <- length(tree$tip.label)
  if (tmproot > nTips) tree <- root(tree, node = tmproot)
  else  tree <- root(tree, tmproot)
  el <- numeric(max(tree$edge))
  el[tree$edge[, 2]] <- tree$edge.length
  maxdm <- el[ind]
  tree$edge.length[tree$edge[, 2] == ind] <- 0
  maxD1 <- node2root(tree)[1:nTips]
  tree$edge.length[tree$edge[, 2] == ind] <- maxdm
  ind <- c(ind, which.max(maxD1))
  maxdm <- maxdm + maxD1[ind[2]]
  rn <- max(tree$edge) + 1L
  edge <- tree$edge
  el <- tree$edge.length
  children <- tree$edge[, 2]
  left <- match(ind[1], children)
  tmp <- Ancestors(tree, ind[2], "all")
  tmp <- c(ind[2], tmp[-length(tmp)])
  right <- match(tmp, children)
  if (el[left] >= (maxdm / 2)) {
    edge <- rbind(edge, c(rn, ind[1]))
    edge[left, 2] <- rn
    el[left] <- el[left] - (maxdm / 2)
    el <- c(el, maxdm / 2)
  }
  else {
    sel <- cumsum(el[right])
    i <- which(sel > (maxdm / 2))[1]
    edge <- rbind(edge, c(rn, tmp[i]))
    edge[right[i], 2] <- rn
    eltmp <-  sel[i] - (maxdm / 2)
    el <- c(el, el[right[i]] - eltmp)
    el[right[i]] <- eltmp
  }
  tree$edge.length <- el
  storage.mode(edge) <- "integer"
  tree$edge <- edge
  tree$Nnode <- tree$Nnode + 1L
  attr(tree, "order") <- NULL
  tree <- reroot(tree, rn)
  if (!is.null(tree$node.label)) {
    node.label <- tree$node.label
    tmp <- node.label[1]
    node.label[1] <- node.label[rn - nTips]
    node.label[rn - nTips] <- tmp
    node.label[is.na(node.label)] <- ""
    tree$node.label <- node.label
  }
  attr(tree, "order") <- NULL
  tree <- reorder(tree)
  if (!is.null(oldtree$node.label)) {
    type <- match.arg(node.labels, c("support", "label", "delete"))
    if (type == "support") tree <- addConfidences.phylo(tree, oldtree)
    if (type == "delete") tree$node.label <- NULL
  }
  tree
}


#' @rdname midpoint
#' @method midpoint multiPhylo
#' @export
midpoint.multiPhylo <- function(tree, node.labels = "support", ...) {
  if (!is.null(attr(tree, "TipLabel"))) compress <- TRUE
  else compress <- FALSE
  tree <- lapply(tree, midpoint.phylo, node.labels = node.labels)
  class(tree) <- "multiPhylo"
  if (compress) tree <- .compressTipLabel(tree)
  tree
}


#' @rdname midpoint
#' @export
pruneTree <- function(tree, ..., FUN = ">=") {
  if (is.null(tree$node)) stop("no node labels")
  # if (is.rooted(tree)) tree <- unroot(tree)
  has_edge.length <- !is.null(tree$edge.length)
  tree <- reorder(tree)
  if(has_edge.length){
    tree <- minEdge(tree)
    nh <- nodeHeight(tree)
  }
  else tree$edge.length <- rep(1,nrow(tree$edge))
  m <- max(tree$edge)
  nTips <- length(tree$tip.label)
  bs <- rep(TRUE, m)
  bs[ (nTips + 1):m] <- sapply(as.numeric(as.character(tree$node)), FUN, ...)
  if(has_edge.length){
    for(i in seq_len(nrow(tree$edge))){
      ei <- tree$edge[i,2]
      if(!(is.na(bs[ei])) && !bs[ei]) nh[ei] <- nh[tree$edge[i,1]]
    }
    tree$edge.length <- nh[tree$edge[,1]] -  nh[tree$edge[,2]]
  }
  else tree$edge.length[!bs[tree$edge[, 2]]] <- 0
  attr(tree, "order") <- NULL
  tree <- di2multi(tree)
  if(!has_edge.length) tree$edge.length <- NULL
  reorder(tree, "postorder")
}


# requires postorder
# for internal use in fitch.spr
# pos statt i
dropTip <- function(x, i, check.binary = FALSE, check.root = TRUE) {
  edge <- x$edge
  root <- edge[nrow(edge), 1] #getRoot(x)
  ch <- match(i, edge[,2]) #which(edge[, 2] == i)
  pa <- edge[ch, 1]
  edge <- edge[-ch, ]
  ind <- which(edge[, 1] == pa)
  if (root == pa) {
    if (length(ind) == 1) {
      edge <- edge[-ind, ]
      x$Nnode <- x$Nnode - 1L
    }
    if (length(ind) == 2) {
      n <- dim(edge)[1]
      newroot <- edge[n - 2L, 1]
      newedge <- edge[ind, 2]
      if (newedge[1] == newroot) edge[n - 1, ] <- newedge
      else edge[n - 1, ] <- newedge[2:1]
      edge <- edge[-n, ]
      x$Nnode <- x$Nnode - 1L
      edge[edge == newroot] <- root
      pa <- newroot
    }
    # todo handle unrooted trees
  }
  else {
    nind <- match(pa, edge[,2]) #which(edge[, 2] == pa)
    # normal binary case
    if (length(ind) == 1) {
      edge[nind, 2] <- edge[ind, 2]
      edge <- edge[-ind, ]
      x$Nnode <- x$Nnode - 1L
    }
  }
  #
  edge[edge > pa]  <- edge[edge > pa] - 1L
  x$edge <- edge
  x
}



# like drop tip and returns two trees,
# to be used in fitch.spr
descAll <- function(x, node, nTips, ch) {
  m <- max(x)
  isInternal <- logical(m)
  isInternal[(nTips + 1):m] <- TRUE
  desc <- function(node, isInternal) {
    if (!isInternal[node]) return(node)
    res <- NULL
    while (length(node) > 0) {
      tmp <- unlist(ch[node])
      res <- c(res, tmp)
      node <- tmp[isInternal[tmp]]
    }
    res
  }
  desc(node, isInternal)
}


dropNode <- function(x, i, check.binary = FALSE, check.root = TRUE,
                     all.ch = NULL) {
  edge <- x$edge
  root <- getRoot(x)
  ch <- match(i, edge[, 2]) # which(edge[, 2] == i)

  nTips <- length(x$tip.label)
  pa <- edge[ch, 1]
  if (i > nTips) {
    if (is.null(all.ch)) all.ch <- allChildren(x)
    kids <- descAll(edge, i, nTips, all.ch)
    ind <- match(kids, edge[, 2])
    edge2 <- edge[sort(ind), ]
    edge <- edge[-c(ch, ind), ]
  }
  else edge <- edge[-ch, ]
  if (nrow(edge) < 3) return(NULL)
  ind <- which(edge[, 1] == pa)
  sibs <- edge[ind, 2L]
  if (root == pa) {
    if (length(ind) == 1) {
      edge <- edge[-ind, ]
      x$Nnode <- x$Nnode - 1L
    }
    if (length(ind) == 2) {
      n <- dim(edge)[1]
      newroot <- edge[n - 2L, 1]
      newedge <- edge[ind, 2]
      if (newedge[1] == newroot) edge[n - 1, ] <- newedge
      else edge[n - 1, ] <- newedge[2:1]
      edge <- edge[-n, ]
      x$Nnode <- as.integer(length(unique(edge[, 1])))
      edge[edge == newroot] <- root
      pa <- newroot
    }
    # todo handle unrooted trees
  }
  else {
    nind <- match(pa, edge[,2]) # which(edge[, 2] == pa)
    # normal binary case
    if (length(ind) == 1) {
      edge[nind, 2] <- edge[ind, 2]
      edge <- edge[-ind, ]
      x$Nnode <- as.integer(length(unique(edge[, 1])))
    }
  }
  x$edge <- edge
  y <- x
  y$edge <- edge2
  y$Nnode <- as.integer(length(unique(edge2[, 1])))
  list(x, y, pa, sibs)
}

# nur mit edge matrix
# postorder remained tip in 1:nTips
addOne <- function(tree, tip, i) {
  edge <- tree$edge
  parent <- edge[, 1]
  l <- dim(edge)[1]
  m <- max(edge) + 1L
  p <- edge[i, 1]
  k <- edge[i, 2]
  edge[i, 2] <- m
  ind <- match(p, parent)
  if (ind == 1) edge <- rbind(matrix(c(m, m, k, tip), 2, 2), edge)
  else edge <- rbind(edge[1:(ind - 1), ], matrix(c(m, m, k, tip), 2, 2),
      edge[ind:l, ])
  tree$edge <- edge
  tree$Nnode <- tree$Nnode + 1L
  tree
}

# raus?
addOneTree <- function(tree, subtree, i, node) {
  edge <- tree$edge
  parent <- edge[, 1]
  l <- dim(edge)[1]
  m <- node # max(edge)+1L
  p <- edge[i, 1]
  k <- edge[i, 2]
  edge[i, 2] <- m
  edge2 <- subtree$edge
  ind <- match(p, parent)
  r2 <- edge2[nrow(edge2), 1]
  if (ind == 1) edge <- rbind(edge2, matrix(c(m, m, r2, k), 2, 2), edge)
  else edge <- rbind(edge[1:(ind - 1), ], edge2, matrix(c(m, m, r2, k), 2, 2),
      edge[ind:l, ])
  tree$edge <- edge
  tree$Nnode <- tree$Nnode + subtree$Nnode + 1L
  attr(tree, "order") <- NULL
  tips1 <- as.integer(length(tree$tip.label) + 1L)
  tmproot <- getRoot(tree)
  if (tmproot != tips1) {
    tree$edge[tree$edge == tmproot] <- 0L
    tree$edge[tree$edge == tips1] <- tmproot
    tree$edge[tree$edge == 0L] <- tips1
  }
  tree <- reorder(tree, "postorder")
  if (tmproot != tips1) tree <- unroot(tree)
  tree
}



#' Add tips to a tree
#'
#' This function binds tips to nodes of a phylogenetic trees.
#'
#'
#' @param tree an object of class "phylo".
#' @param tips a character vector containing the names of the tips.
#' @param where an integer or character vector of the same length as tips giving
#' the number of the node or tip of the tree where to add the new tips.
#' @param edge.length optional numeric vector with edge length
#' @return an object of class phylo
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{bind.tree}}
#' @keywords cluster
#' @examples
#' tree <- rcoal(10)
#' plot(tree)
#' nodelabels()
#' tiplabels()
#' tree1 <- add.tips(tree, c("A", "B", "C"), c(1,2,15))
#' plot(tree1)
#' @export
add.tips <- function(tree, tips, where, edge.length = NULL) {
  nTips <- length(tree$tip.label)
  nTips_new <- length(tips)
  if (nTips_new < 1) return(tree)
  edge <- tree$edge
  if (is.character(where)) {
    where <- match(where, c(tree$tip.label, tree$node.label))
  }
  ind <- match(where, edge[, 2])

  n_internal <- as.integer(sum(unique(where) <= nTips))
  edge[edge > nTips] <- edge[edge > nTips] + nTips_new
  p_vec <- integer(max(edge) + n_internal)
  p_vec[edge[, 2]] <- edge[, 1]
  tip_index <- (nTips + 1):(nTips + nTips_new)
  c_vec <- c(edge[, 2], tip_index)
  # first handle internal nodes (easy)
  if (any(where > nTips)) {
    ind1 <- where > nTips
    p_vec[tip_index[ind1]] <- where[ind1]  + nTips_new
  }
  # handle tips
  if (any(where <= nTips)) {
    m <- max(edge)
    tmp <- unique(where)
    tmp <- tmp[tmp <= nTips]
    new_internal <- as.integer( (m + 1L):(m + n_internal))
    # add new internal node
    p_vec[new_internal] <- edge[match(tmp, edge[, 2]), 1]
    p_vec[tmp] <- new_internal
    # add tip
    ind2 <- (where <= nTips)
    p_vec[tip_index[ind2]] <- p_vec[where[ind2]]
    ind <- match(tmp, edge[, 2])
    c_vec[ind] <- new_internal
    c_vec <- c(c_vec, tmp)
    if (!is.null(tree$node.label)) {
      tree$node.label <- c(tree$node.label, rep("", n_internal))
    }
  }
  tree$edge <- matrix(c(p_vec[c_vec], c_vec), ncol = 2)

  if (!is.null(tree$edge.length)) {
    if (is.null(edge.length)) {
      tree$edge.length <- c(tree$edge.length,
        rep(0, nTips_new + n_internal))
    }
    else {
      if (length(edge.length) < nTips_new) edge.length <- rep(edge.length,
          length.out = nTips_new)
      tree$edge.length <- c(tree$edge.length, edge.length,
        rep(0, n_internal))
    }
  }
  tree$Nnode <- tree$Nnode + n_internal
  tree$tip.label <- c(tree$tip.label, tips)
  attr(tree, "order") <- NULL
  tree <- reorder(tree)
  if (!is.null(tree$edge.length)) {
    if (is.null(edge.length)) {
      nh <- nodeHeight(tree)
      nh[tip_index] <- 0
      tree$edge.length <- nh[tree$edge[, 1]] - nh[tree$edge[, 2]]
    }
  }
  tree
}



#' Compute all trees topologies.
#'
#' \code{allTrees} computes all bifurcating tree topologies for rooted or
#' unrooted trees with up to 10 tips. The number of trees grows fast.
#'
#' @param n Number of tips (<=10).
#' @param rooted Rooted or unrooted trees (default: rooted).
#' @param tip.label Tip labels.
#' @return an object of class \code{multiPhylo}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{rtree}}, \code{\link{nni}},
#' \code{\link[ape]{howmanytrees}}, \code{\link{dfactorial}}
#' @keywords cluster
#' @examples
#'
#' trees <- allTrees(5)
#'
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow = c(3,5))
#' for(i in 1:15)plot(trees[[i]])
#' par(old.par)
#'
#' @export allTrees
allTrees <- function(n, rooted = FALSE, tip.label = NULL) {
  n <- as.integer(n)
  nt <- as.integer(round(dfactorial(2 * (n + rooted) - 5)))
  Nnode <- as.integer(n - 2L + rooted)
  if ( (n + rooted) > 10) {
    stop(gettextf("That would generate %d trees, and take up more than %d MB of memory!",
      nt, as.integer(round(nt / 1000)), domain = "R-phangorn"))
  }
  if (n < 2) {
    stop("A tree must have at least two taxa.")
  }
  if (!rooted && n == 2) {
    stop("An unrooted tree must have at least three taxa.")
  }

  if (rooted) {
    edge <- matrix(NA, 2 * n - 2, 2)
    edge[1:2, ] <- c(n + 1L, n + 1L, 1L, 2L)
  }
  else {
    edge <- matrix(NA,  2 * n - 3, 2)
    edge[1:3, ] <- c(n + 1L, n + 1L, n + 1L, 1L, 2L, 3L)
  }
  edges <- list()
  edges[[1]] <- edge

  m <- 1
  nedge <- 1
  trees <- vector("list", nt)
  if ( (n + rooted) > 3) {
    i <- 3L  + (!rooted)
    pa <- n + 2L
    nr <- 2L + (!rooted)
    while (i < (n + 1L)) {
      nedge <- nedge + 2
      m2 <- m * nedge

      newedges <- vector("list", m2)
      for (j in 1:m) {
        edge <- edges[[j]]
        l <- nr

        edgeA <- edge
        edgeB <- edge

        for (k in 1L:l) {
          edge <- edgeA
          node <- edge[k, 1]
          edge[k, 1] <- pa
          edge[l + 1, ] <- c(pa, i)
          edge[l + 2, ] <- c(node, pa)

          newedges[[(j - 1) * (l + rooted) + k]] <- edge
        }

        if (rooted) {
          edgeB[] <- as.integer(sub(n + 1L, pa, edgeB))
          edge <- edgeB
          edge[l + 1, ] <- c(n + 1L, i)
          edge[l + 2, ] <- c(n + 1L, pa)
          newedges[[j * (l + 1)]] <- edge
        }
      } # end for
      edges <- newedges
      m <- m2
      i <- i + 1L
      pa <- pa + 1L
      nr <- nr + 2L
    } # end for m
  } # end if
  for (x in 1:m) {
    edge <- edges[[x]]
    edge <- edge[reorderRcpp(edge, n, n + 1L, 2L), ]
    tree <- list(edge = edge)
    tree$Nnode <- Nnode
    attr(tree, "order") <- "postorder"
    class(tree) <- "phylo"
    trees[[x]] <- tree
  }
  attr(trees, "TipLabel") <- if (is.null(tip.label))
    paste("t", 1:n, sep = "")
  else tip.label
  class(trees) <- "multiPhylo"
  trees
}


#
# some generic tree functions
#
allAncestors <- function(x) {
  x <- reorder(x, "postorder")
  parents <- x$edge[, 1]
  child <- x$edge[, 2]
  l <- length(parents)
  res <- vector("list", max(x$edge))
  for (i in l:1) {
    pa <- parents[i]
    res[[child[i]]] <- c(pa, res[[pa]])
  }
  res
}

char2pos <- function(x, node){
  if(is.null(x$node.label)){
    tmp <- as.character(seq(Ntip(x)+1, Ntip(x)+Nnode(x)))
    labels <- c(x$tip.label, tmp)
  }
  else  labels <- c(x$tip.label, x$node.label)
  x <- match(node, labels)
  if(any(is.na(x))) stop("Can't find supplied node in the labels")
  x
}

## @aliases Ancestors Children Descendants Siblings mrca.phylo
#' tree utility function
#'
#' Functions for describing relationships among phylogenetic nodes.
#'
#' These functions are inspired by \code{treewalk} in phylobase package, but
#' work on the S3 \code{phylo} objects.  The nodes are the indices as given in
#' edge matrix of an phylo object. From taxon labels these indices can be
#' easily derived matching against the \code{tip.label} argument of an phylo
#' object, see example below.  All the functions allow \code{node} to be either
#' a scalar or vector.  \code{mrca} is a faster version of the mrca in ape, in
#' phangorn only because of dependencies.
#' If the argument node is missing the function is evaluated for all nodes.
#'
#' @param x a tree (a phylo object).
#' @param node an integer or character vector (or scalar) corresponding to a
#' node ID
#' @param type specify whether to return just direct children / parents or all
#' @param include.self whether to include self in list of siblings
#' @param full a logical indicating whether to return the MRCAs among all tips
#' and nodes (if TRUE); the default is to return only the MRCAs among tips.
#' @return a vector or a list containing the indices of the nodes.
#' @seealso \code{treewalk}, \code{\link[ape]{as.phylo}},
#' \code{\link[ape]{nodelabels}}
#' @keywords misc
#' @examples
#'
#' tree <- rtree(10)
#' plot(tree, show.tip.label = FALSE)
#' nodelabels()
#' tiplabels()
#' Ancestors(tree, 1:3, "all")
#' Children(tree, 11)
#' Descendants(tree, 11, "tips")
#' Siblings(tree, 3)
#' # Siblings of all nodes
#' Siblings(tree)
#' mrca.phylo(tree, 1:3)
#' mrca.phylo(tree, match(c("t1", "t2", "t3"), tree$tip))
#' mrca.phylo(tree)
#' # same as mrca(tree), but faster for large trees
#'
#' @export
#' @rdname Ancestors
Ancestors <- function(x, node, type = c("all", "parent")) {
  if(!missing(node) && inherits(node, "character")) node <- char2pos(x, node)
  parents <- x$edge[, 1]
  child <- x$edge[, 2]
  pvector <- integer(max(x$edge)) # parents
  pvector[child] <- parents
  type <- match.arg(type)
  if (type == "parent")
    return(pvector[node])
  anc <- function(pvector, node) {
    res <- integer(0)
    repeat {
      anc <- pvector[node]
      if (anc == 0) break
      res <- c(res, anc)
      node <- anc
    }
    res
  }
  if (!missing(node) && length(node) == 1) return(anc(pvector, node))
  else allAncestors(x)[node]
}


allChildren <- function(x) {
  l <- length(x$tip.label)
  if (l < 20) {
    parent <- x$edge[, 1]
    children <- x$edge[, 2]
    res <- vector("list", max(x$edge))
    for (i in seq_along(parent)) res[[parent[i]]] <- c(res[[parent[i]]],
                                                       children[i])
    return(res)
  }
  else {
    allChildrenCPP(x$edge)
  }
}

#' @describeIn Ancestors list all the descendant nodes of each node
#' @keywords internal
#' @export
allDescendants <- function(x) {
  if (is.null(attr(x, "order")) || attr(x, "order") != "postorder")
    x <- reorder(x, "postorder")
  nTips <- as.integer(length(x$tip.label))
  allDescCPP(x$edge, nTips)
}


#' @rdname Ancestors
#' @export
Children <- function(x, node) {
  # return allChildren if node is missing
  if(!missing(node) && inherits(node, "character")) node <- char2pos(x, node)
  if (!missing(node) && length(node) == 1)
    return(x$edge[x$edge[, 1] == node, 2])
  allChildren(x)[node]
}


#' @rdname Ancestors
#' @export
Descendants <- function(x, node, type = c("tips", "children", "all")) {
  type <- match.arg(type)
  if(!missing(node) && inherits(node, "character")) node <- char2pos(x, node)
  if (type == "children") return(Children(x, node))
  if (type == "tips") return(bip(x)[node])
  # new version using Rcpp
  if (missing(node) || length(node) > 10) return(allDescendants(x)[node])
  ch <- allChildren(x) # out of the loop
  isInternal <- logical(max(x$edge))
  isInternal[ unique(x$edge[, 1]) ] <- TRUE
  desc <- function(node, isInternal) {
    if (!isInternal[node]) return(node)
    res <- NULL
    while (length(node) > 0) {
      tmp <- unlist(ch[node])
      res <- c(res, tmp)
      node <- tmp[isInternal[tmp]]
    }
    res
  }
  if (length(node) > 1) return(lapply(node, desc, isInternal))
  desc(node, isInternal)
}


#' @rdname Ancestors
#' @export
Siblings <- function(x, node, include.self = FALSE) {
  if (missing(node)) node <- as.integer(1:max(x$edge))
  if(!missing(node) && inherits(node, "character")) node <- char2pos(x, node)
  l <- length(node)
  if (l == 1) {
    v <- Children(x, Ancestors(x, node, "parent"))
    if (!include.self)
      v <- v[v != node]
    return(v)
  }
  else {
    parents <- x$edge[, 1]
    child <- x$edge[, 2]
    pvector <- integer(max(x$edge)) # parents
    pvector[child] <- parents
    root <- as.integer(parents[!match(parents, child, 0)][1])
    res <- vector("list", l)
    ch <- allChildren(x)
    if (include.self) return(ch[ pvector[node] ])
    k <- 1
    for (i in node) {
      if (i != root) {
        tmp <- ch[[ pvector[i] ]]
        res[[k]] <- tmp[tmp != i]
      }
      k <- k + 1
    }
  }
  res
}


#' @rdname Ancestors
#' @export
mrca.phylo <- function(x, node = NULL, full = FALSE) {
  if (is.null(node)) return(mrca2(x, full = full))
  if(!missing(node) && inherits(node, "character")) node <- char2pos(x, node)
  return(getMRCA(x, node))
}


# should be in ape
mrca2 <- function(phy, full = FALSE) {
  Nnode <- phy$Nnode
  Ntips <- length(phy$tip.label)
  phy <- reorder(phy)
  nodes <- unique(phy$edge[, 1])
  if (!full) {
    res <- Descendants(phy, nodes, "tips")
    M <- matrix(nodes[1], Ntips, Ntips)
    for (i in 2:Nnode) M[res[[i]], res[[i]]] <- nodes[i]
    diag(M) <- 1:Ntips
    dimnames(M) <-  list(phy$tip.label, phy$tip.label)
  }
  else {
    res <- Descendants(phy, nodes, "all")
    M <- matrix(nodes[1], Ntips + Nnode, Ntips + Nnode)
    for (i in 2:Nnode) {
      tmp  <-  c(res[[i]], nodes[i])
      M[tmp, tmp] <- nodes[i]
    }
    diag(M) <- 1:(Ntips + Nnode)
    dimnames(M) <-  list(1:(Ntips + Nnode), 1:(Ntips + Nnode))
  }
  M
}


relabel <- function(y, ref) {
  label <- y$tip.label
  if (identical(label, ref)) return(y)
  if (length(label) != length(ref))
    stop("one tree has a different number of tips")
  ilab <- match(label, ref)
  if (any(is.na(ilab)))
    stop("one tree has different tip labels")
  ie <- match(seq_along(ref), y$edge[, 2])
  y$edge[ie, 2] <- ilab
  y$tip.label <- ref
  y
}

#' @rdname midpoint
#' @param labels tip and node labels to keep as tip labels in the tree
#' @export
keep_as_tip<- function(tree, labels){
  nodes_to_keep <- labels[!is.na(match(labels, tree$node.label))]
  tips_to_remove <-  Descendants(tree, nodes_to_keep) |> unlist() |> unique()
  tree_1 <- drop.tip(tree, tips_to_remove, subtree = TRUE)
  tree_2 <- keep.tip(tree_1, labels, collapse.singles=FALSE)
  tree_2
}
