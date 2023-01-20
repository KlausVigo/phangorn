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


degree <- function(x){
  tabulate(x$edge)
}


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
  g <- graph(t(edge), directed = FALSE)
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
  X <- X[ind, ]

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
    g <- graph(t(res$edge), directed = FALSE)

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
#' @importFrom igraph graph
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
#    graph(t(x$edge), directed=directed)
# }



#' Computes a consensusNetwork from a list of trees Computes a \code{networx}
#' object from a collection of splits.
#'
#' Computes a consensusNetwork, i.e. an object of class \code{networx} from a
#' list of trees, i.e. an class of class \code{multiPhylo}. Computes a
#' \code{networx} object from a collection of splits.
#'
#'
#' @param obj An object of class multiPhylo.
#' @param prob the proportion a split has to be present in all trees to be
#' represented in the network.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{consensusNet} returns an object of class networx.  This is
#' just an intermediate to plot phylogenetic networks with igraph.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{splitsNetwork}}, \code{\link{neighborNet}},
#' \code{\link{lento}}, \code{\link{distanceHadamard}},
#' \code{\link{plot.networx}}, \code{\link{maxCladeCred}}
#' @references Holland B.R., Huber K.T., Moulton V., Lockhart P.J. (2004) Using
#' consensus networks to visualize contradictory evidence for species
#' phylogeny. \emph{Molecular Biology and Evolution}, \bold{21}, 1459--61
#' @keywords hplot
#' @examples
#'
#' data(Laurasiatherian)
#' set.seed(1)
#' bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)nj(dist.hamming(x)),
#'     bs=50)
#' cnet <- consensusNet(bs, .3)
#' plot(cnet)
#' \dontrun{
#' library(rgl)
#' open3d()
#' plot(cnet, type = "3D", show.tip.label=FALSE, show.nodes=TRUE)
#' plot(cnet, type = "equal angle", show.edge.label=TRUE)
#'
#' tmpfile <- normalizePath(system.file(
#'               "extdata/trees/RAxML_bootstrap.woodmouse", package="phangorn"))
#' trees <- read.tree(tmpfile)
#' cnet_woodmouse <- consensusNet(trees, .3)
#' plot(cnet_woodmouse, type = "equal angle", show.edge.label=TRUE)
#' }
#'
#' @export consensusNet
consensusNet <- function(obj, prob = 0.3, ...) {
  l <- length(obj)
  spl <- as.splits(obj)
  w <- attr(spl, "weights")
  ind <- (w / l) > prob
  spl <- spl[ind]
  attr(spl, "confidences") <- (w / l)[ind]
  #    attr(spl, "weights") = w[ind]
  res <- as.networx(spl)
  res$edge.labels <- as.character(res$edge.length / l * 100)
  res$edge.labels[res$edge[, 2] <= length(res$tip.label)] <- ""
  reorder(res)
}


#' @export
reorder.networx <- function(x, order =  "cladewise", index.only = FALSE, ...) {
  order <- match.arg(order, c("cladewise", "postorder"))
  if (!is.null(attr(x, "order")))
    if (attr(x, "order") == order)
      return(x)
  g <- graph(t(x$edge))
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


# some trigonemetric functions
rad2deg <- function(rad) (rad * 180) / (pi)
deg2rad <- function(deg) (deg * pi) / (180)

# circular mean
# https://en.wikipedia.org/wiki/Mean_of_circular_quantities
circ.mean <- function(deg) {
  rad.m <- (deg * pi) / (180)
  mean.cos <- mean(cos(rad.m))
  mean.sin <- mean(sin(rad.m))

  theta <- rad2deg(atan(mean.sin / mean.cos))
  if (mean.cos < 0) theta <- theta + 180
  if ((mean.sin < 0) & (mean.cos > 0)) theta <- theta + 360
  theta
}


spl2angle <- function(x) {
  l <- length(attr(x, "labels"))
  ord <- 1:l
  if (!is.null(attr(x, "cycle"))) ord <- attr(x, "cycle")
  x <- changeOrder(x, attr(x, "labels")[ord])
  y <- lapply(x, function(x, l) (x - 1) / l * 360, l = l)
  angle <- vapply(y, circ.mean, 0) |> deg2rad()
  # angle <- ((vapply(x, sum, 0) / lengths(x) - 1) / l ) * 2*pi
  # kreis2kart(attr(x, "weight"), angle)
  angle
}


coords.equal.angle <- function(obj) {
  if (is.null(attr(obj, "order")) || (attr(obj, "order") == "postorder"))
    obj <- reorder.networx(obj)
  spl <- obj$splits
  spl <- SHORTwise(spl) #, length(obj$tip.label))
  l <- length(obj$edge.length)
#  ind1 <- which(!duplicated(obj$splitIndex))
  n <- max(obj$edge)
  angle <- spl2angle(spl)
  weight <- attr(spl, "weight")
  k <- matrix(0, max(obj$splitIndex), 2)

  res <- matrix(0, max(obj$edge), 2)
  for (i in 1:l) { # unique(obj$splitIndex)
    j <- obj$edge[i, 1]
    m <- obj$edge[i, 2]
    p <- obj$splitIndex[i]
    res[m, ] <- res[j, ] + kreis2kart(weight[p], angle[p])
  }
  res
}


#' @rdname phangorn-internal
#' @export
coords <- function(obj, dim = "3D") {
  #    if(is.null(attr(obj,"order")) || (attr(obj, "order")=="postorder") )
  #        obj = reorder.networx(obj)

  if (dim == "equal_angle") return(coords.equal.angle(obj))

  l <- length(obj$edge.length)
  ind1 <- which(!duplicated(obj$splitIndex))

  n <- max(obj$edge)
  adj <- spMatrix(n, n, i = obj$edge[, 2], j = obj$edge[, 1],
                  x = rep(1, length(obj$edge.length)))
  g <- graph_from_adjacency_matrix(adj, "undirected")
  ##########
  #    add this
  #    g2 <- graph(t(obj$edge), directed=FALSE)
  #    g2 <- set.edge.attribute(g, "weight", value=rep(1, nrow(obj$edge))
  if (dim == "3D") {
    coord <- layout_nicely(g, dim = 3)
    k <- matrix(0, max(obj$splitIndex), 3)
    for (i in ind1) {
      tmp <- coord[obj$edge[i, 2], ] - coord[obj$edge[i, 1], ]
      k[obj$splitIndex[i], ] <- kart2kugel(tmp[1], tmp[2], tmp[3])
    }
    k[obj$splitIndex[ind1], 1] <- obj$edge.length[ind1]

    res <- matrix(0, vcount(g), 3)
    for (i in 1:l) {
      j <- obj$edge[i, 1]
      m <- obj$edge[i, 2]
      p <- obj$splitIndex[i]
      res[m, ] <- res[j, ] + kugel2kart(k[p, 1], k[p, 2], k[p, 3])
    }
  }
  else {
    coord <- layout_nicely(g, dim = 2)
    k <- matrix(0, max(obj$splitIndex), 2)
    for (i in ind1) {
      tmp <- coord[obj$edge[i, 2], ] - coord[obj$edge[i, 1], ]
      k[obj$splitIndex[i], ] <- kart2kreis(tmp[1], tmp[2])
    }
    k[obj$splitIndex[ind1], 1] <- obj$edge.length[ind1]
    res <- matrix(0, vcount(g), 2)
    for (i in 1:l) {
      j <- obj$edge[i, 1]
      m <- obj$edge[i, 2]
      p <- obj$splitIndex[i]
      res[m, ] <- res[j, ] + kreis2kart(k[p, 1], k[p, 2])
    }
  }
  res
}


kart2kugel <- function(x, y, z) {
  r <- sqrt(x * x + y * y + z * z)
  alpha <- atan(sqrt(x * x + y * y) / z)
  if (z < 0) alpha <- alpha + pi
  beta <- atan(y / x)
  if (x < 0) beta <- beta + pi
  c(r, alpha, beta)
}


kart2kreis <- function(x, y) {
  r <- sqrt(x * x + y * y)
  alpha <- atan(y / x)
  if (x < 0) alpha <- alpha + pi
  c(r, alpha)
}


kreis2kart <- function(r, alpha) {
  c(r * cos(alpha), r * sin(alpha))
  #    if(length(r)>1) return(matrix(c(r*cos(alpha), r*sin(alpha)), ncol=2))
  # 	else return(c(r*cos(alpha), r*sin(alpha)))
}


kugel2kart <- function(r, alpha, beta) {
  x <- r * sin(alpha) * cos(beta)
  y <- r * sin(alpha) * sin(beta)
  z <- r * cos(alpha)
  c(x, y, z)
}


edgeLabels <- function(xx, yy, zz = NULL, edge) {
  XX <- (xx[edge[, 1]] + xx[edge[, 2]]) / 2
  YY <- (yy[edge[, 1]] + yy[edge[, 2]]) / 2
  if (!is.null(zz)) {
    ZZ <- (zz[edge[, 1]] + zz[edge[, 2]]) / 2
    return(cbind(XX, YY, ZZ))
  }
  cbind(XX, YY)
}


#' plot phylogenetic networks
#'
#' So far not all parameters behave the same on the the \code{rgl} \code{"3D"}
#' and basic graphic \code{"2D"} device.
#'
#' Often it is easier and safer to supply vectors of graphical parameters for
#' splits (e.g. splits.color) than for edges. These overwrite values edge.color.
#'
#' @param x an object of class \code{"networx"}
#' @param type "3D" to plot using rgl or "equal angle" and "2D" in the normal
#' device.
#' @param use.edge.length a logical indicating whether to use the edge weights
#' of the network to draw the branches (the default) or not.
#' @param show.tip.label a logical indicating whether to show the tip labels on
#' the graph (defaults to \code{TRUE}, i.e. the labels are shown).
#' @param show.edge.label a logical indicating whether to show the tip labels
#' on the graph.
#' @param edge.label an additional vector of edge labels (normally not needed).
#' @param show.node.label a logical indicating whether to show the node labels
#' (see example).
#' @param node.label an additional vector of node labels (normally not needed).
#' @param show.nodes a logical indicating whether to show the nodes (see
#' example).
#' @param tip.color the colors used for the tip labels.
#' @param edge.color the colors used to draw edges.
#' @param edge.width the width used to draw edges.
#' @param edge.lty a vector of line types.
#' @param split.color the colors used to draw edges.
#' @param split.width the width used to draw edges.
#' @param split.lty a vector of line types.
#' @param font an integer specifying the type of font for the labels: 1 (plain
#' text), 2 (bold), 3 (italic, the default), or 4 (bold italic).
#' @param cex a numeric value giving the factor scaling of the labels.
#' @param cex.node.label a numeric value giving the factor scaling of the node
#' labels.
#' @param cex.edge.label a numeric value giving the factor scaling of the edge
#' labels.
#' @param col.node.label the colors used for the node labels.
#' @param col.edge.label the colors used for the edge labels.
#' @param font.node.label the font used for the node labels.
#' @param font.edge.label the font used for the edge labels.
#' @param underscore a logical specifying whether the underscores in tip labels
#' should be written as spaces (the default) or left as are (if TRUE).
#' @param \dots Further arguments passed to or from other methods.
#' @rdname plot.networx
#' @note The internal representation is likely to change.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{consensusNet}}, \code{\link{neighborNet}},
#' \code{\link{splitsNetwork}}, \code{\link{hadamard}},
#' \code{\link{distanceHadamard}}, \code{\link{as.networx}},
#' \code{\link[ape]{evonet}}, \code{\link[ape]{as.phylo}},
#' \code{\link{densiTree}}, \code{\link[ape]{nodelabels}}
#' @references Dress, A.W.M. and Huson, D.H. (2004) Constructing Splits Graphs
#' \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics
#' (TCBB)}, \bold{1(3)}, 109--115
#'
#' Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W. (2017),
#' Intertwining phylogenetic trees and networks. \emph{Methods Ecol Evol}.
#' \bold{8}, 1212--1220. doi:10.1111/2041-210X.12760
#' @keywords plot
#' @importFrom igraph graph
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
#' @importFrom igraph graph_from_adjacency_matrix vcount topo_sort layout_nicely
#' @method plot networx
#' @export
plot.networx <- function(x, type = "equal angle", use.edge.length = TRUE,
                         show.tip.label = TRUE, show.edge.label = FALSE,
                         edge.label = NULL, show.node.label = FALSE,
                         node.label = NULL, show.nodes = FALSE,
                         tip.color = "black", edge.color = "black",
                         edge.width = 3, edge.lty = 1, split.color = NULL,
                         split.width = NULL, split.lty = NULL, font = 3,
                         cex = par("cex"), cex.node.label = cex,
                         cex.edge.label = cex, col.node.label = tip.color,
                         col.edge.label = tip.color, font.node.label = font,
                         font.edge.label = font, underscore = FALSE, ...) {
  type <- match.arg(type, c("equal angle", "3D", "2D"))
  if (use.edge.length == FALSE){
    x$edge.length[] <- 1
    attr(x$splits, "weight") <- rep(1, length(x$splits))
  }
  nTips <- length(x$tip.label)
  conf <- attr(x$splits, "confidences")
  index <- x$splitIndex
  if (is.null(edge.label) & !is.null(conf)) {
    conf <- conf[index]
    if (!is.null(x$translate)) conf[match(x$translate$node, x$edge[, 2])] <- ""
    else conf[x$edge[, 2] <= nTips] <- ""
    edge.label <- conf
  }
  if (is.null(node.label)) node.label <- as.character(1:max(x$edge))
  if (show.tip.label) node.label[1:nTips] <- ""
  if (show.tip.label){
    if (is.expression(x$tip.label)) underscore <- TRUE
    if (!underscore) x$tip.label <- gsub("_", " ", x$tip.label)
  }

  lspl <- max(x$splitIndex)
  if (!is.null(split.color)) {
    if (length(split.color) != lspl)
      stop("split.color must be same length as splits")
    else edge.color <- split.color[x$splitIndex]
  }
  if (!is.null(split.width)) {
    if (length(split.width) != lspl)
      stop("split.color must be same length as splits")
    else edge.width <- split.width[x$splitIndex]
  }
  if (!is.null(split.lty)) {
    if (length(split.lty) != lspl)
      stop("split.color must be same length as splits")
    else edge.lty <- split.lty[x$splitIndex]
  }

  chk <- FALSE

  if (type == "3D") chk <- requireNamespace("rgl", quietly = TRUE)
  if (!chk && type == "3D") {
    warning("type='3D' requires the package 'rgl', plotting in '2D' instead!\n")
    type <- "2D"
  }
  # use precomputed vertices when available
  coord <- NULL
  if (!is.null(x$.plot)) coord <- x$.plot$vertices

  if (type == "3D") {
    if (is.null(coord) || ncol(coord) != 3)
      coord <- coords(x, dim = "3D")
    plotRGL(coord, x, show.tip.label = show.tip.label,
      show.edge.label = show.edge.label, edge.label = edge.label,
      show.node.label = show.node.label, node.label = node.label,
      show.nodes = show.nodes, tip.color = tip.color, edge.color = edge.color,
      edge.width = edge.width, font = font, cex = cex,
      cex.node.label = cex.node.label, cex.edge.label = cex.edge.label,
      col.node.label = col.node.label, col.edge.label = col.edge.label,
      font.node.label = font.node.label, font.edge.label = font.edge.label)
  }
  else {
    if (is.null(coord) || ncol(coord) != 2) {
      if (type == "equal angle") coord <- coords.equal.angle(x)
      else coord <- coords(x, dim = "2D")
    }
    plot2D(coord, x, show.tip.label = show.tip.label,
      show.edge.label = show.edge.label, edge.label = edge.label,
      show.node.label = show.node.label, node.label = node.label,
      show.nodes = show.nodes, tip.color = tip.color, edge.color = edge.color,
      edge.width = edge.width, edge.lty = edge.lty, font = font, cex = cex,
      cex.node.label = cex.node.label, cex.edge.label = cex.edge.label,
      col.node.label = col.node.label, col.edge.label = col.edge.label,
      font.node.label = font.node.label, font.edge.label = font.edge.label,
      add = FALSE)
  }
  x$.plot <- list(vertices = coord, edge.color = edge.color,
    edge.width = edge.width, edge.lty = edge.lty)
  L <- list(Ntip = nTips, type = "networx")
  assign("last_plot.phylo", c(L, list(edge = x$edge, xx = coord[, 1],
    yy = coord[, 2])), envir = .PlotPhyloEnv)
  invisible(x)
}


plotRGL <- function(coords, net, show.tip.label = TRUE, show.edge.label = FALSE,
                    edge.label = NULL, show.node.label = FALSE,
                    node.label = NULL, show.nodes = FALSE, tip.color = "blue",
                    edge.color = "grey", edge.width = 3, font = 3,
                    cex = par("cex"), cex.node.label = cex,
                    cex.edge.label = cex, col.node.label = tip.color,
                    col.edge.label = tip.color, font.node.label = font,
                    font.edge.label = font, ...) {
  open3d <- rgl::open3d
  segments3d  <- rgl::segments3d
  spheres3d <- rgl::spheres3d
  texts3d <- rgl::texts3d

  edge <- net$edge

  x <- coords[, 1]
  y <- coords[, 2]
  z <- coords[, 3]

  nTips <- length(net$tip.label)

  segments3d(x[t(edge)], y[t(edge)], z[t(edge)],
             col = rep(edge.color, each = 2), lwd = edge.width)
  radius <- 0
  if (show.nodes) {
    radius <- sqrt( (max(x) - min(x))^2 + (max(y) - min(y))^2 +
                      (max(z) - min(z))^2) / 200
    spheres3d(x[1:nTips], y[1:nTips], z[1:nTips], radius = 2 * radius,
              color = "cyan")
    spheres3d(x[-c(1:nTips)], y[-c(1:nTips)], z[-c(1:nTips)], radius = radius,
              color = "magenta")
  }
  if (show.tip.label) {
    if (is.null(net$translate))
      texts3d(x[1:nTips] + 2.05 * radius, y[1:nTips], z[1:nTips],
                net$tip.label, color = tip.color, cex = cex, font = font)
    else
      texts3d(x[net$translate$node] + 2.05 * radius, y[net$translate$node],
        z[net$translate$node], net$tip.label, color = tip.color, cex = cex,
        font = font)
  }
  if (show.edge.label) {
    ec <- edgeLabels(x, y, z, edge)
    if (is.null(edge.label)) edge.label <- net$splitIndex
    # else edge.label = net$splitIndex
    texts3d(ec[, 1], ec[, 2], ec[, 3], edge.label, color = col.edge.label,
      cex = cex.edge.label, font = font.edge.label)
  }
  if (show.node.label) {
    texts3d(x, y, z, node.label, color = col.node.label, cex = cex.node.label,
      font = font.node.label)
  }
}


plot2D <- function(coords, net, show.tip.label = TRUE, show.edge.label = FALSE,
                   edge.label = NULL, show.node.label = FALSE,
                   node.label = NULL, tip.color = "blue", edge.color = "grey",
                   edge.width = 3, edge.lty = 1, font = 3, cex = par("cex"),
                   cex.node.label = cex,  cex.edge.label = cex,
                   col.node.label = tip.color, col.edge.label = tip.color,
                   font.node.label = font, font.edge.label = font,
                   add = FALSE, ...) {
  edge <- net$edge
  label <- net$tip.label
  xx <- coords[, 1]
  yy <- coords[, 2]
  nTips <- length(label)

  xlim <- range(xx)
  ylim <- range(yy)

  if (show.tip.label) {
    offset <- max(nchar(label)) * 0.018 * cex * diff(xlim)
    xlim <- c(xlim[1] - offset, xlim[2] + offset)
    ylim <- c(ylim[1] - 0.03 * cex * diff(ylim), ylim[2] +
                0.03 * cex * diff(ylim))
  }
  if (!add) {
    plot.new()
    plot.window(xlim, ylim, asp = 1)
  }
  cladogram.plot(edge, xx, yy, edge.color, edge.width, edge.lty)
  if (show.tip.label) {
    if (is.null(net$translate)) ind <- match(1:nTips, edge[, 2])
    else ind <- match(net$translate$node, edge[, 2])
    pos <- rep(4, nTips)
    XX <- xx[edge[ind, 1]] - xx[edge[ind, 2]]
    pos[XX > 0] <- 2
    YY <- yy[edge[ind, 1]] - yy[edge[ind, 2]]
    pos2 <- rep(3, nTips)
    pos2[YY > 0] <- 1
    # needed if tiplabels are not at internal nodes
    XX[is.na(XX)] <- 0
    YY[is.na(YY)] <- 0
    pos[abs(YY) > abs(XX)] <- pos2[abs(YY) > abs(XX)]
    if (is.null(net$translate)) text(xx[1:nTips], yy[1:nTips], labels = label,
        pos = pos, col = tip.color, cex = cex, font = font)
    else text(xx[net$translate$node], yy[net$translate$node], labels = label,
        pos = pos, col = tip.color, cex = cex, font = font)
  }
  if (show.edge.label) {
    ec <- edgeLabels(xx, yy, edge = edge)
    if (is.null(edge.label)) edge.label <- net$splitIndex

    # show only one edge label
    em <- apply(ec, 1, function(x) max(abs(x)))
    si <- net$splitIndex
    for (i in unique(si)) {
      tmp <- si == i
      if (sum(tmp) > 1) {
        w <- which(tmp)
        wm <- which.max(em[w])
        edge.label[w[-wm]] <- ""
      }
    }

    text(ec[, 1], ec[, 2], labels = edge.label, col = col.edge.label,
      cex = cex.edge.label, font = font.edge.label)
  }
  if (show.node.label) {
    text(xx, yy, labels = node.label, col = col.node.label,
         cex = cex.node.label, font = font.node.label)
  }
}


closest.edge <- function(x, y, P1, P2) {
  x1 <- P1[, 1]
  x2 <- P2[, 1]
  y1 <- P1[, 2]
  y2 <- P2[, 2]

  A <- sqrt( (x2 - x)^2 + (y2 - y)^2)    # d_BC
  B <- sqrt( (x1 - x)^2 + (y1 - y)^2)    # d_AC
  C <- sqrt( (x1 - x2)^2 + (y1 - y2)^2)  # d_AB
  # Kosinussatz
  alpha <- acos( (B^2 + C^2 - A^2) / (2 * B * C))
  beta <- acos( (A^2 + C^2 - B^2) / (2 * A * C))

  d <- abs( (y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
    sqrt( (y2 - y1)^2 + (x2 - x1)^2)
  d[alpha > (pi / 2)] <- B[alpha > (pi / 2)]
  d[beta > (pi / 2)] <- A[beta > (pi / 2)]
  d
}

closest.node <- function(x, y, P) {
  x1 <- P[, 1]
  y1 <- P[, 2]
  d <- sqrt((x1 - x)^2 + (y1 - y)^2)
  d
}


#' Identify splits in a network
#'
#' \code{identify.networx} reads the position of the graphics pointer when the
#' mouse button is pressed. It then returns the split belonging to the edge
#' closest to the pointer. The network must be plotted beforehand.
#'
#' @param x an object of class \code{networx}
#' @param quiet a logical controlling whether to print a message inviting the
#' user to click on the tree.
#' @param \dots further arguments to be passed to or from other methods.
#' @return \code{identify.networx} returns a splits object.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[phangorn]{plot.networx}},
#'   \code{\link[graphics]{identify}}
#' @examples
#' \dontrun{
#' data(yeast)
#' dm <- dist.ml(yeast)
#' nnet <- neighborNet(dm)
#' plot(nnet)
#' identify(nnet) # click close to an edge
#' }
#' @importFrom graphics identify
#' @method identify networx
#' @export
identify.networx <- function(x, quiet = FALSE, ...) {
  if (!quiet)
    cat("Click close to a node or edge of the tree...\n")
  xy <- locator(1)
  if (is.null(xy))
    return(NULL)
  if (is.null(x$.plot)) {
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    edge <- lastPP$edge
    xx <- lastPP$xx
    yy <- lastPP$yy
    vertices <- cbind(xx, yy)
  }
  else {
    lastPP <- x$.plot
    edge <- x$edge
    vertices <- lastPP$vertices
  }
  P1 <- vertices[edge[, 1], , drop = FALSE]
  P2 <- vertices[edge[, 2], , drop = FALSE]
  d <- closest.edge(xy$x, xy$y, P1, P2)
  split <- x$splitIndex[which.min(d)]
  x$splits[split]
}
