#' @rdname as.splits
#' @export
allSplits <- function(k, labels = NULL) {
  result <- lapply(1:(2^(k - 1) - 1), dec2Bin)
  if (is.null(labels)) labels <- (as.character(1:k))
  attr(result, "labels") <- labels
  class(result) <- "splits"
  result
}


#' @rdname as.splits
#' @export
allCircularSplits <- function(k, labels = NULL) {
  fun <- function(x, y) {
    tmp <- (1L:y) + x
    tmp %% (k + 1L) + tmp %/% (k + 1L)
  }

  k <- as.integer(k)
  l <- (k - 1L) %/% 2L
  res <- vector("list", k * (k - 1L) / 2)
  res[1:k] <- 1L:k
  ind <- k
  if (k > 3) {
    if (k > 4L) {
      for (i in 2:l) {
        res[(ind + 1):(ind + k)] <- lapply(0L:(k - 1L), fun, i)
        ind <- ind + k
      }
    }
    if ( (k %% 2L) == 0) {
      m <- k %/% 2
      res[(ind + 1):(ind + m)] <- lapply(0L:(m - 1L), fun, m)
    }
  }
  res <- lapply(res, sort)
  if (is.null(labels)) labels <- as.character(1:k)
  attr(res, "labels") <- labels
  attr(res, "cycle") <- 1:k
  class(res) <- "splits"
  res
}


#degree <- function(x){
#  tabulate(x$edge)
#}


splits2design <- function(obj, weight = NULL, x=TRUE) {
  labels <- attr(obj, "labels")
  m <- length(labels)
  n <- length(obj)
  l <- 1:m
  sl <- lengths(obj)
  p0 <- sl * (m - sl)
  p <- c(0, cumsum(p0))
  i <- numeric(max(p))
  for (k in 1:n) {
    sp <- obj[[k]]
    if (p0[k] != 0) i[(p[k] + 1):p[k + 1]] <- getIndex(sp, l[-sp], m)
  }
  dims <- c(m * (m - 1) / 2, n)
  if(x) return(sparseMatrix(i = i, p = p, x=1, dims = dims) )
  sparseMatrix(i = i, p = p, dims = dims)
}


addEdge <- function(network, desc, spl) {
  edge <- network$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  nTips <- length(network$tip.label)

  desc2 <- SHORTwise(desc) #, nTips)
  split <- desc2[spl]

  index <- network$splitIndex
  ind <- which(compatible(split, desc2[index]) == 1)
  if (is.null(ind) | (length(ind) == 0)) return(network)
  add <- TRUE

  v <- sort(match(split[[1]], edge[,2]))

  fromTo <- intersect(attr(desc, "cycle"), split[[1]])
  fromTo <- parent[match(fromTo, child)]
  g <- make_graph(t(edge), directed = FALSE)
  ind <- NULL
  for (i in 2:length(fromTo)) {
    d <- all_shortest_paths(g, fromTo[i - 1], fromTo[i])$res
    sptmp <- unlist( lapply(d, \(d) E(g, path=d)) )
    ind <- c(ind, sptmp) # [-c(1, length(sptmp))])
  }
  ind <- unique(ind)

  oldNodes <- unique(as.vector(edge[ind, ]))
  mNodes <- max(network$edge)
  newNodes <- (mNodes + 1L):(mNodes + length(oldNodes))

  # duplicated splits
  dSpl <- edge[ind, ]
  edge2 <- edge[v, ]
  for (i in seq_along(oldNodes)) {
    edge2[edge2 == oldNodes[i]] <- newNodes[i]
  }
  edge[v, ] <- edge2

  # alle Splits verdoppeln
  for (i in seq_along(oldNodes)) dSpl[dSpl == oldNodes[i]] <- newNodes[i]
  edge <- rbind(edge, dSpl, deparse.level = 0) # experimental: no labels
  index <- c(index, index[ind])
  # neu zu alt verbinden
  edge <- rbind(edge, cbind(oldNodes, newNodes), deparse.level = 0)
  index <- c(index, rep(spl, length(oldNodes)))
  network$edge <- edge
  network$Nnode <- max(edge) - nTips
  network$splitIndex <- index
  network
}

## as.splits.phylo
circNetwork <- function(x, ord = NULL) {
  if (is.null(ord)) ord <- attr(x, "cycle")

  weight <- attr(x, "weights")
  if (is.null(weight)) weight <- rep(1, length(x))
  nTips <- length(ord)
  tmp <- which(ord == 1)
  if (tmp != 1) ord <- c(ord[tmp:nTips], ord[1:(tmp - 1)])
  res <- stree(nTips, tip.label = attr(x, "labels"))
  res$edge[, 2] <- ord
  res$edge.length <- NULL
  x <- SHORTwise(x) #, nTips)
  spRes <- as.splits(res)[res$edge[, 2]]
  index <- match(spRes, x)

  if (any(is.na(index))) {
    l.na <- sum(is.na(index))
    x <- c(x, spRes[is.na(index)])
    weight <- c(weight, rep(0, l.na))
    index <- match(spRes, x)
  }

  l <- lengths(ONEwise(x))
  l2 <- lengths(x)

  tmp <- countCycles(x, ord = ord)
  ind <- which(tmp == 2 & l2 > 1) # & l<nTips changed with ordering

  #    ind = ind[order(l[ind])]
  ind <- ind[order(l2[ind], decreasing = TRUE)]

  dm2 <- as.matrix(compatible(x, x[ind]))

  X <- as.matrix(x)[, ord]
  Y <- X
  rsY <- rowSums(Y)
  X <- X[ind, , drop=FALSE]

  for (k in seq_along(ind)) {
    Vstart <- ord[1]
    Vstop <- ord[nTips]
    ordStart <- 1
    ordStop <- nTips
    for (j in 2:nTips) {

      if (X[k, j - 1] < X[k, j]) {
        Vstart <- ord[j]
        ordStart <- j
      }
      if (X[k, j - 1] > X[k, j]) {
        Vstop <- ord[j - 1]
        ordStop <- j - 1
      }
    }
    fromTo <- ordStart:ordStop
    if (ordStart > ordStop) fromTo <- c(ordStart:nTips, 1:ordStop)
    fromTo <- ord[fromTo]
    g <- make_graph(t(res$edge), directed = FALSE)

    isChild <- (rsY == (Y %*% X[k, ]))[index]
    sp2 <- NULL
    sp0 <- NULL

    for (i in 2:length(fromTo)) {
      sptmp <- shortest_paths(g, fromTo[i - 1], fromTo[i],
        output = c("epath"))$epath[[1]]
      sp2 <- c(sp2, sptmp[-c(1, length(sptmp))])
      sp0 <- c(sp0, sptmp)
    }
    sp0 <- unique(sp0)

    if (length(sp2) > 0) {
      #            blub = which(dm[index[sp2], ind[k]]>0)
      TMP <- rowSums(dm2[index[sp2], 1:k, drop = FALSE])
      blub <- which(TMP > 0)
      sp2 <- sp2[blub]
    }
    if (length(sp2) == 0) {
      isChild <- (rsY == (Y %*% X[k, ]))[index]
      sp0 <- which(isChild == TRUE)
      edge1 <- unique(as.vector(res$edge[sp0, ]))
      edge2 <- as.vector(res$edge[-sp0, ])
      asdf <- edge1 %in% edge2
      sp <- edge1[asdf]
    }
    if (length(sp2) > 0) sp <- unique(as.vector(t(res$edge[sp2, ])))
    parent <- res$edge[, 1]
    child <- res$edge[, 2]

    j <- ord[which(X[k, ] == 1)]
    anc <- unique(parent[match(j, child)])

    maxVert <- max(parent)
    l <- length(sp)

    newVert <- (maxVert + 1):(maxVert + l)
    sp01 <- setdiff(sp0, sp2)
    for (i in 1:l) res$edge[sp01, ][res$edge[sp01, ] == sp[i]] <- newVert[i]

    newindex <- rep(ind[k], l)
    if (length(sp) > 1) newindex <- c(index[sp2], newindex)
    index <- c(index, newindex)
    # connect new and old vertices
    newEdge <- matrix(cbind(sp, newVert), ncol = 2)
    if (length(sp) > 1) {
      # copy edges
      qwer <- match(as.vector(res$edge[sp2, ]), sp)
      newEdge <- rbind(matrix(newVert[qwer], ncol = 2), newEdge)
    }

    res$edge <- rbind(res$edge, newEdge)
    res$Nnode <-  max(res$edge) - nTips

    res$splitIndex <- index
    res$edge.length <- rep(1, nrow(res$edge))
    class(res) <- c("networx", "phylo")
    attr(res, "order") <- NULL
  }
  res$edge.length <- weight[index]  # ausserhalb
  res$Nnode <-  max(res$edge) - nTips
  res$splitIndex <- index
  res$splits <- x
  class(res) <- c("networx", "phylo")
  attr(res, "order") <- NULL
  res
}


#' Conversion among phylogenetic network objects
#'
#' \code{as.networx} convert \code{splits} objects into a \code{networx}
#' object. And most important there exists a generic \code{plot} function to
#' plot phylogenetic network or split graphs.
#'
#' @details A \code{networx} object hold the information for a phylogenetic
#' network and extends the \code{phylo} object. Therefore some generic function
#' for \code{phylo} objects will also work for \code{networx} objects.  The
#' argument \code{planar = TRUE} will create a planar split graph based on a
#' cyclic ordering. These objects can be nicely plotted in \code{"2D"}.
#'
#' @aliases networx
#' @param x an object of class \code{"splits"} or \code{"phylo"}
#' @param planar logical whether to produce a planar graph from only cyclic
#' splits (may excludes splits).
#' @param coord add coordinates of the nodes, allows to reproduce the plot.
#' @param \dots Further arguments passed to or from other methods.
#' @returns an object of class \code{networx}.
#' @note The internal representation is likely to change.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{consensusNet}}, \code{\link{neighborNet}},
#' \code{\link{splitsNetwork}}, \code{\link{hadamard}},
#' \code{\link{distanceHadamard}}, \code{\link{plot.networx}},
#' \code{\link[ape]{evonet}}, \code{\link[ape]{as.phylo}}
#' @references
#' Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W. (2017),
#' Intertwining phylogenetic trees and networks. \emph{Methods Ecol Evol}.
#' \bold{8}, 1212--1220. doi:10.1111/2041-210X.12760
#' @keywords plot
#' @importFrom igraph make_graph
#' @examples
#'
#' set.seed(1)
#' tree1 <- rtree(20, rooted=FALSE)
#' sp <- as.splits(rNNI(tree1, n=10))
#' net <- as.networx(sp)
#' plot(net)
#' \dontrun{
#' # also see example in consensusNet
#' example(consensusNet)
#' }
#'
#' @rdname as.networx
#' @export as.networx
as.networx <- function(x, ...) {
  if (inherits(x, "networx"))
    return(x)
  UseMethod("as.networx")
}



#' @export
addTrivialSplits <- function(obj) {
  label <- attr(obj, "labels")
  nTips <- length(label)
  weight <- attr(obj, "weights")
  if (is.null(weight)) weight <- rep(1, length(obj))
  STree <- stree(nTips, tip.label = attr(obj, "labels"))
  STree$edge.length <- NULL
  spRes <- as.splits(STree)[STree$edge[, 2]]
  tmpIndex <- match(spRes, SHORTwise(obj))
  if (any(is.na(tmpIndex))) {
    l.na <- sum(is.na(tmpIndex))
    obj <- c(obj, spRes[is.na(tmpIndex)])
    weight <- c(weight, rep(0, l.na))
    attr(obj, "weights") <- weight
  }
  obj
}


#' @export
removeTrivialSplits <- function(obj) {
  nTips <- length(attr(obj, "labels"))
  l <- lengths(obj)
  ind <- which( (l == 0L) | (l == 1L) | (l == nTips) | (l == (nTips - 1L)))
  obj[-ind]
}


#' @rdname as.networx
#' @importFrom igraph shortest_paths all_shortest_paths decompose E
#' @importFrom Matrix spMatrix
#' @method as.networx splits
#' @export
as.networx.splits <- function(x, planar = FALSE, coord = "none", ...) {
  label <- attr(x, "labels")
  coord <- match.arg(coord, c("none", "equal angle", "3D", "2D"))
  x <- addTrivialSplits(x)
  nTips <- length(label)
  x <- ONEwise(x)
  l <- lengths(x)
  if (any(l == nTips)) x <- x[l != nTips] # get rid of trivial splits
  l <- lengths(x)

  weight <- attr(x, "weights")
  if (is.null(weight)) weight <- rep(1, length(x))
  attr(x, "weights") <- weight

  ext <- sum(l == 1 | l == (nTips - 1))
  if (!is.null(attr(x, "cycle"))) {
    c.ord <- attr(x, "cycle")
  }
  else c.ord <- cophenetic(x) |> getOrderingNN()

  attr(x, "cycle") <- c.ord
  # check for star tree
  if(length(x)==nTips) return(as.phylo(x))
  # which splits are in circular ordering
  circSplits <- which(countCycles(x, ord = c.ord) == 2)
  if (length(circSplits) == length(x)) planar <- TRUE
  tmp <- circNetwork(x, c.ord)
  attr(tmp, "order") <- NULL
  if (planar) {
    return(reorder(tmp))
  }

  dm <- as.matrix(compatible(x))
  ll <- lengths(x)
  ind <- tmp$splitIndex     # match(sp, x)
  ind2 <- union(ind, which(ll == 0)) # which(duplicated(x))
  ind2 <- union(ind2, which(ll == nTips))
  ord <- order(colSums(dm))
  ord <- setdiff(ord, ind2)
  if (length(ord) > 0) {
    for (i in seq_along(ord)) {
      tmp <- addEdge(tmp, x, ord[i])
      tmp$edge.length <- weight[tmp$splitIndex]
      tmp$Nnode <- max(tmp$edge) - nTips
      class(tmp) <- c("networx", "phylo")
    }
  }
  tmp$edge.length <- weight[tmp$splitIndex]
  tmp$Nnode <- max(tmp$edge) - nTips
  attr(x, "cycle") <- c.ord
  tmp$splits <- x
  class(tmp) <- c("networx", "phylo")
  tmp <- reorder(tmp)
  coord <- match.arg(coord)
  vert <- switch(coord,
    "none" = NULL,
    "equal angle" = coords.equal.angle(tmp),
    "2D" = coords(tmp, dim = "2D"),
    "3D" = coords(tmp, dim = "3D"))
  #    attr(tmp, "coords") <- coordinates
  tmp$plot <- list(vertices = vert)
  tmp
}



#' @rdname as.networx
#' @method as.networx phylo
#' @export
as.networx.phylo <- function(x, ...) {
  spl <- as.splits(x)
  spl <- spl[x$tree[, 2]]
  x$splitIndex <- seq_len( nrow(x$edge) )
  x$splits <- spl
  class(x) <- c("networx", "phylo")
  x
}


# as.igraph.networx <- function(x, directed=FALSE){
#    make_graph(t(x$edge), directed=directed)
# }


#' @export
reorder.networx <- function(x, order =  "cladewise", index.only = FALSE, ...) {
  order <- match.arg(order, c("cladewise", "postorder"))
  if (!is.null(attr(x, "order")))
    if (attr(x, "order") == order)
      return(x)
  g <- make_graph(t(x$edge))
  if (order == "cladewise") neword <- topo_sort(g, "out")
  else neword <- topo_sort(g, "in")
  neworder <- order(match(x$edge[, 1], neword))
  if (index.only) return(neworder)
  x$edge <- x$edge[neworder, ]
  if (!is.null(x$edge.length))
    x$edge.length <- x$edge.length[neworder]
  if (!is.null(x$edge.labels))
    x$edge.labels <- x$edge.labels[neworder]
  if (!is.null(x$splitIndex)) x$splitIndex <- x$splitIndex[neworder]
  attr(x, "order") <- order
  x
}
