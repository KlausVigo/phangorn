#' Neighbor-Joining
#'
#' This function performs the neighbor-joining tree estimation of Saitou and
#' Nei (1987). UNJ is the unweighted version from Gascuel (1997).
#'
#'
#' @param x A distance matrix.
#' @return an object of class \code{"phylo"}.
#' @author Klaus P. Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{nj}}, \code{\link[ape]{dist.dna}},
#' \code{\link[phangorn]{dist.hamming}}, \code{\link[phangorn]{upgma}},
#' \code{\link[ape]{fastme}}
#' @references Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new
#' method for reconstructing phylogenetic trees. \emph{Molecular Biology and
#' Evolution}, \bold{4}, 406--425.
#'
#' Studier, J. A and Keppler, K. J. (1988) A Note on the Neighbor-Joining
#' Algorithm of Saitou and Nei. \emph{Molecular Biology and Evolution},
#' \bold{6}, 729--731.
#'
#' Gascuel, O. (1997) Concerning the NJ algorithm and its unweighted version,
#' UNJ. in Birkin et. al. \emph{Mathematical Hierarchies and Biology},
#' 149--170.
#' @keywords cluster
#' @examples
#'
#' data(Laurasiatherian)
#' dm <- dist.ml(Laurasiatherian)
#' tree <- NJ(dm)
#' plot(tree)
#'
#' @rdname NJ
#' @export
NJ <- function(x) reorder(nj(x), "postorder")


#' @rdname NJ
#' @export
UNJ <- function(x){
  x <- as.matrix(x)
  labels <- attr(x, "Labels")[[1]]
  edge.length <- NULL
  edge <- NULL
  d <- as.matrix(x)
  if (is.null(labels))
    labels <- colnames(d)
  l <- dim(d)[1]
  n <- l
  nam <- as.character(1:l)
  m <- l - 2
  nam <- 1:l
  k <- 2 * l - 2
  w <- rep(1, l)
  while (l > 2) {
    r <- rowSums(d) / (l - 2)
#    i <- 0
#    j <- 0
    tmp <- out_cpp(d, r, l)
    e2 <- tmp[2]
    e1 <- tmp[1]
    l1 <- d[e1, e2] / 2 + sum( (d[e1, -c(e1, e2)] - d[e2, -c(e1, e2)]) *
                                w[-c(e1, e2)]) / (2 * (n - w[e1] - w[e2]))
    l2 <- d[e1, e2] / 2 + sum( (d[e2, -c(e1, e2)] - d[e1, -c(e1, e2)]) *
                                w[-c(e1, e2)]) / (2 * (n - w[e1] - w[e2]))
    edge.length <- c(l1, l2, edge.length)
    edge <- rbind(c(k, nam[e2]), edge)
    edge <- rbind(c(k, nam[e1]), edge)
    nam <- c(nam[c(-e1, -e2)], k)

    dnew <- (w[e1] * d[e1, ] + w[e2] * d[e2, ] - w[e1] * l1 - w[e2] * l2) /
      (w[e1] + w[e2])
    d <- cbind(d, dnew)
    d <- rbind(d, c(dnew, 0))
    d <- d[-c(e1, e2), -c(e1, e2)]
    w <- c(w, w[e1] + w[e2])
    w <- w[-c(e1, e2)]
    k <- k - 1
    l <- l - 1
  }
  edge.length <- c(d[2, 1], edge.length)
  result <- list(edge = rbind(c(nam[2], nam[1]), edge),
    edge.length = edge.length, tip.label = labels, Nnode = m)
  class(result) <- "phylo"
  reorder(result, "postorder")
}


#
# Distance Matrix methods
#


#' Compute a design matrix or non-negative LS
#'
#' \code{nnls.tree} estimates the branch length using non-negative least
#' squares given a tree and a distance matrix.  \code{designTree} and
#' \code{designSplits} compute design matrices for the estimation of edge
#' length of (phylogenetic) trees using linear models.  For larger trees a
#' sparse design matrix can save a lot of memory. %\code{designTree} also
#' computes a contrast matrix if the method is "rooted".
#'
#' @param tree an object of class \code{phylo}
#' @param sparse return a sparse design matrix.
#' @param x number of taxa.
#' @param splits one of "all", "star".
#' @param dm a distance matrix.
#' @param method compute an "unrooted", "ultrametric" or "tipdated" tree.
#' @param rooted compute a "ultrametric" or "unrooted" tree (better use method).
#' @param trace defines how much information is printed during optimization.
#' @param \dots further arguments, passed to other methods.
#' @param weight vector of weights to be used in the fitting process.
#' Weighted least squares is used with weights w, i.e., sum(w * e^2) is
#' minimized.
#' @param balanced use weights as in balanced fastME
#' @param tip.dates a named vector of sampling times associated to the tips of
#' the tree.
#' @param calibration a named vector of calibration times associated to nodes of
#' the tree.
#' @param eps minimum edge length (default s 1e-8).
## @param strict strict calibration.
#' @return \code{nnls.tree} return a tree, i.e. an object of class
#' \code{phylo}.  \code{designTree} and \code{designSplits} a matrix, possibly
#' sparse.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{fastme}}, \code{\link[ape]{rtt}},
#' \code{\link[phangorn]{distanceHadamard}},
#' \code{\link[phangorn]{splitsNetwork}}, \code{\link[phangorn]{upgma}}
#' @keywords cluster
#' @importFrom Matrix Matrix sparseMatrix crossprod solve
#' @importFrom quadprog solve.QP.compact
#' @examples
#'
#' example(NJ)
#' dm <-  as.matrix(dm)
#' y <- dm[lower.tri(dm)]
#' X <- designTree(tree)
#' lm(y~X-1)
#' # avoids negative edge weights
#' tree2 <- nnls.tree(dm, tree)
#'
#' @rdname designTree
#' @export
designTree <- function(tree, method = "unrooted", sparse = FALSE,
                       tip.dates=NULL, calibration=NULL, ...) { # , strict=TRUE
  method <- match.arg(method,
                      c("unrooted", "ultrametric", "rooted", "tipdated"))
  if(method == "rooted") method <- "ultrametric"
  if(has.singles(tree)) tree <- collapse.singles(tree)
  #if (!is.na(pmatch(method, "all")))
  #  method <- "unrooted"
  #METHOD <- c("unrooted", "rooted", "tipdated")
  #method <- pmatch(method, METHOD)
  #if (is.na(method)) stop("invalid method")
  #if (method == -1) stop("ambiguous method")
  #if (!is.rooted(tree) & method == 2) stop("tree has to be rooted")
  if(method == "unrooted") X <- designUnrooted(tree, sparse = sparse, ...)
  if(method == "ultrametric") X <- designUltra(tree, sparse = sparse, ...)
  if(method == "tipdated") X <- designTipDated(tree, tip.dates=tip.dates,
                                               sparse=sparse, ...)
  X
}


designUnrooted <- function(tree, sparse=FALSE, order = NULL) {
  if (inherits(tree, "phylo")) {
    if (is.rooted(tree)) tree <- unroot(tree)
    tree <- reorder(tree, "postorder")
    p <- as.matrix(as.splits(tree)[tree$edge[,2]])
#    p <- bipartition(tree)
  }
  if (inherits(tree, "splits")) p <- as.matrix(tree)
  if (!is.null(order))
    p <- p[, order]

  m <- dim(p)[2]
  ind <- rowSums(p)
  p <- p[ind != m, ]
  n <- dim(p)[1]
  res <- matrix(0, (m - 1) * m / 2, n)
  k <- 1
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      res[k, ] <- p[, i] != p[, j]
      k <- k + 1
    }
  }
  if (inherits(tree, "phylo"))
    colnames(res) <- paste(tree$edge[, 1], tree$edge[, 2], sep = "<->")
  if(sparse) res <- Matrix(res, sparse=TRUE)
  res
}


designUltra <- function(tree, sparse = TRUE, calibration=NULL) {
  if (is.null(attr(tree, "order")) || attr(tree, "order") != "postorder")
    tree <- reorder(tree, "postorder")
#  stopifnot( !(!is.null(calibration) && is.null(tree$node.label)))
  leri <- allChildren(tree)
  bp <- bip(tree)
  n <- length(tree$tip.label)
  l <- tree$Nnode
  nodes <- integer(l)
  k <- 1L
  u <- numeric(n * (n - 1) / 2)
  v <- numeric(n * (n - 1) / 2)
  m <- 1L
  for (i in seq_along(leri)) {
    if (length(leri[[i]]) > 1) {
      if (length(leri[[i]]) == 2) ind <- getIndex(bp[[leri[[i]][1] ]],
          bp[[leri[[i]][2] ]], n)
      else {
        ind <- NULL
        le <- leri[[i]]
        nl <- length(le)
        for (j in 1:(nl - 1)) ind <- c(ind, getIndex(bp[[le[j] ]],
            unlist(bp[ le[(j + 1):nl] ]), n))
      }
      li <- length(ind)
      v[m:(m + li - 1)] <- k
      u[m:(m + li - 1)] <- ind
      nodes[k] <- i
      m <- m + li
      k <- k + 1L
    }
  }
  if (sparse) X <- sparseMatrix(i = u, j = v, x = 2L)
  else {
    X <- matrix(0L, n * (n - 1) / 2, l)
    X[cbind(u, v)] <- 2L
  }
  colnames(X) <- nodes
  attr(X, "nodes") <- nodes
  X
}


designUnrooted2 <- function(tree, sparse = TRUE) {
  if (is.null(attr(tree, "order")) || attr(tree, "order") != "postorder")
    tree <- reorder(tree, "postorder")
  leri <- allChildren(tree)
  bp <- bip(tree)
  n <- length(tree$tip.label)
  l <- tree$Nnode
  nodes <- integer(l)
  nTips <- as.integer(length(tree$tip.label))
  k <- nTips
  u <- numeric(n * (n - 1) / 2)
  v <- numeric(n * (n - 1) / 2)
  z <- numeric(n * (n - 1) / 2)
  y <- numeric(n * (n - 1) / 2)
  p <- 1L
  m <- 1L
  for (i in seq_along(leri)) {
    if (length(leri[[i]]) > 1) {
      if (length(leri[[i]]) == 2) {
        ind <-  getIndex(bp[[leri[[i]][1] ]], bp[[leri[[i]][2] ]], n)
        ytmp <- rep(bp[[leri[[i]][1] ]], each = length(bp[[leri[[i]][2] ]]))
        ztmp <- rep(bp[[leri[[i]][2] ]], length(bp[[leri[[i]][1] ]]))
      }
      else {
        ind <- NULL
        le <- leri[[i]]
        nl <- length(le)
        ytmp <- NULL
        ztmp <- NULL
        for (j in 1:(nl - 1)) {
          bp1 <- bp[[le[j] ]]
          bp2 <- unlist(bp[le[(j + 1):nl] ])
          ind <- c(ind,  getIndex(bp1, unlist(bp2), n))
          ytmp <- c(ytmp, rep(bp1, each = length(bp2)))
          ztmp <- c(ztmp, rep(bp2, length(bp1)))
        }
      }

      li <- length(ind)
      v[m:(m + li - 1)] <- k
      u[m:(m + li - 1)] <- ind
      y[m:(m + li - 1)] <- ytmp
      z[m:(m + li - 1)] <- ztmp

      nodes[p] <- i
      m <- m + li
      k <- k + 1L
      p <- p + 1L
    }
  }
  jj <- c(y, z) # [ind],v)
  ii <- c(u, u) # [ind],u)
  ind <-  (jj < nTips)
  jj <- c(jj[ind], v)
  ii <- c(ii[ind], u)
  l1 <- length(u)
  l2 <- sum(ind)
  x <- rep(c(-1L, 2L), c(l2, l1))

  X <- sparseMatrix(i = ii, j = jj, x = x)
  if (!sparse) {
    X <- as.matrix(X)
  }
  nodes <- c(1:(nTips - 1L), nodes)
  colnames(X) <- nodes
  attr(X, "nodes") <- nodes
  X
}


designTipDated <- function(tree, tip.dates, sparse=TRUE){
  #, strict=TRUE
  #if(!is.numeric(tip.dates)) browser()
  #if(!length(tip.dates) >= Ntip(tree)) browser()
  stopifnot(is.numeric(tip.dates), length(tip.dates) >= Ntip(tree))
  nTip <- Ntip(tree)
  tmp <- function(n){
    x1 <- rep(seq_len(n), each=n)
    x2 <- rep(seq_len(n), n)
    ind <- x1 < x2
    sparseMatrix(i = rep(seq_len(sum(ind)), 2), j = c(x1[ind], x2[ind]))
  }
  tip.dates <- tip.dates - max(tip.dates)
  x <- tmp(nTip) %*% tip.dates
  nodes <- integer(tree$Nnode)
  X <- designUltra(tree, sparse=sparse)
  nodes <- attr(X, "nodes")
  X <- cbind(X, x)
  colnames(X) <- c(nodes, -1)
  attr(X, "nodes") <- nodes
  X
}


designCalibrated <- function(tree, sparse=TRUE, calibration=NULL){
  #, tip.dates=NULL, strict=TRUE
  #stopifnot(is.numeric(tip.dates), length(tip.dates) >= Ntip(tree))
  nTip <- Ntip(tree)
  #if(!is.null(tree$node.label))
  cname <- tree$node.label
#  nodes <- integer(tree$Nnode)
  X <- designUltra(tree, sparse=sparse)
  #if(!is.null(tree$node.label))
  colnames(X) <- cname
  x <- X[, names(calibration), drop=FALSE] %*% calibration
  X <- cbind(X[,-match(names(calibration), cname)], rate=x)
#  nodes <- attr(X, "nodes")
#  X <- cbind(X, x)
#  colnames(X) <- c(nodes, -1)
#  attr(X, "nodes") <- nodes
  X
}


designConstrained <- function(tree, sparse=TRUE, tip.dates=NULL,
                              calibration=NULL){
  stopifnot(is.numeric(tip.dates), length(tip.dates) >= Ntip(tree))
  X <- designUltra(tree, sparse=sparse)
  nTip <- Ntip(tree)
  # designTipDated
  if(!is.null(tip.dates)){
    tmp <- function(n){
      x1 <- rep(seq_len(n), each=n)
      x2 <- rep(seq_len(n), n)
      ind <- x1 < x2
      sparseMatrix(i = rep(seq_len(sum(ind)), 2), j = c(x1[ind], x2[ind]))
    }
  }
  if(!is.null(calibration)){
    cname <- tree$node.label
    colnames(X) <- cname
    x <- X[, names(calibration), drop=FALSE] %*% calibration
    X <- cbind(X[,-match(names(calibration), cname)], rate=x)
  }


}

#' @rdname designTree
#' @export
nnls.tree <- function(dm, tree, method=c("unrooted", "ultrametric", "tipdated"),
          rooted=NULL, trace=1, weight=NULL, balanced=FALSE, tip.dates=NULL) {
  method <- match.arg(method, c("unrooted", "ultrametric", "tipdated"))
  if(has.singles(tree)) tree <- collapse.singles(tree)
  if (is.rooted(tree) && method == "unrooted") tree <- unroot(tree)
  tree <- reorder(tree, "postorder")
  if (balanced) {
    if (!is.binary(tree)) stop("tree must be binary")
    weight <- rowSums(designTree(unroot(tree)))
  }
  dm <- as.matrix(dm)
  k <- dim(dm)[1]
  labels <- tree$tip.label
  dm <- dm[labels, labels]
  y <- dm[lower.tri(dm)]
  # computing the design matrix from the tree
  if(!is.null(rooted)){
    if(isTRUE(rooted)) method <- "ultrametric"
    if(isFALSE(rooted)) method <- "unrooted"
  }
  X <- switch(method,
              unrooted=designUnrooted2(tree),
              ultrametric=designUltra(tree),
              tipdated=designTipDated(tree, tip.dates))
#  if (rooted) X <- designUltra(tree)
#  else X <- designUnrooted2(tree)

  if (!is.null(weight)) {
    y <- y * sqrt(weight)
    X <- X * sqrt(weight)
  }

  lab <- attr(X, "nodes")
  ll <- length(lab) + 1L
  # na.action
  if (any(is.na(y))) {
    ind <- which(is.na(y))
    X <- X[-ind, , drop = FALSE]
    y <- y[-ind]
  }
  # LS solution
  Dmat <- crossprod(X) # cross-product computations
  dvec <- crossprod(X, y)
  betahat <- as.vector(solve(Dmat, dvec))
  betahattmp <- betahat
  bhat <- numeric(max(tree$edge))
  if(method=="tipdated"){
    nh <-  max(tip.dates) - tip.dates
    rate <-  betahat[ll]
    bhat[seq_len(Ntip(tree))] <- nh * rate
    bhat[as.integer(lab)] <- betahat[-ll]  # * (1/rate)  [-ll]
  }
  else bhat[as.integer(lab)] <- betahat
  betahat <- bhat[tree$edge[, 1]] - bhat[tree$edge[, 2]]

  if (!any(betahat < 0)) {
    RSS <- sum((y - (X %*% betahattmp))^2)
    if (trace > 1) print(paste("RSS:", RSS))
    attr(tree, "RSS") <- RSS
    if(method=="tipdated"){
      betahat <- betahat / rate
      attr(tree, "rate") <- rate
    }
    tree$edge.length <- betahat
    return(tree)
  }

  # non-negative LS
  n <- dim(X)[2]
  l <- nrow(tree$edge)

  lab <- attr(X, "nodes")
  # vielleicht solve.QP.compact
  ind1 <- match(tree$edge[, 1], lab)
  ind2 <- match(tree$edge[, 2], lab)

  Amat <- matrix(0, 2, l)
  Amat[1, ] <- 1
  Amat[2, ] <- -1

  Aind <- matrix(0L, 3, l)
  Aind[1, ] <- 2L
  Aind[2, ] <- as.integer(ind1)
  Aind[3, ] <- as.integer(ind2)
  if(method == "tipdated"){
    ind3 <- match(seq_len(Ntip(tree)), tree$edge[,2])
    Amat[2, ind3] <- -nh  ## testen
    Aind[3, ind3] <- length(lab) + 1L
  }

  if (any(is.na(Aind))) {
    na_ind <- which(is.na(Aind), arr.ind = TRUE)
    Aind[is.na(Aind)] <- 0L
    for (i in seq_len(nrow(na_ind)) ){
      Aind[1, na_ind[i, 2]] <- Aind[1, na_ind[i, 2]] - 1L
    }
  }

  betahat <- quadprog::solve.QP.compact(as.matrix(Dmat), as.vector(dvec), Amat,
                                        Aind)$sol


  # quadratic programing solving
  RSS <- sum((y - (X %*% betahat))^2)
  if (trace > 1) print(paste("RSS:", RSS))
  attr(tree, "RSS") <- RSS

  bhat <- numeric(max(tree$edge))

  if(method=="tipdated"){
    rate <-  betahat[ll]
    bhat[seq_len(Ntip(tree))] <- nh * rate
    bhat[as.integer(lab)] <- betahat[-ll]
  }
  else  bhat[as.integer(lab)] <- betahat
  betahat <- bhat[tree$edge[, 1]] - bhat[tree$edge[, 2]]
  if(method=="tipdated") {
    betahat <- betahat / rate
    attr(tree, "rate") <- rate
  }
  tree$edge.length <- betahat
  tree
}



#' @rdname designTree
#' @export
nnls.phylo <- function(x, dm, method = "unrooted", trace = 0, ...) {
  nnls.tree(dm, x, method, trace = trace, ...)
}


#' @rdname designTree
#' @export
nnls.splits <- function(x, dm, trace = 0, eps = 1e-8) {
  labels <- attr(x, "labels")
  dm <- as.matrix(dm)
  k <- dim(dm)[1]
  dm <- dm[labels, labels]
  y <- dm[lower.tri(dm)]

  x <- SHORTwise(x) #, k) # use ape version
  l <- lengths(x)
  if (any(l == 0)) x <- x[-which(l == 0)]

  X <- splits2design(x)

  if (any(is.na(y))) {
    ind <- which(is.na(y))
    X <- X[-ind, , drop = FALSE]
    y <- y[-ind]
  }

  Dmat <- crossprod(X) # cross-product computations
  dvec <- crossprod(X, y)
  betahat <- as.vector(solve(Dmat, dvec))

#  if (!any(betahat < 0)) {
#    RSS <- sum((y - (X %*% betahat))^2)
#    if (trace > 1) print(paste("RSS:", RSS))
#    attr(x, "RSS") <- RSS
#    attr(x, "weights") <- betahat
#    return(x)
#  }
  int <- lengths(x)
  if (any(betahat < 0)) {
    n <- dim(X)[2]
    # quadratic programing
    Amat <- matrix(1, 1, n)
    Aind <- matrix(0L, 2L, n)
    Aind[1, ] <- 1L
    Aind[2, ] <- as.integer(1L:n)
    betahat <- quadprog::solve.QP.compact(as.matrix(Dmat), as.vector(dvec),
                                          Amat, Aind)$sol
  }
  RSS <- sum((y - (X %*% betahat))^2)
  ind <- (betahat > 1e-8) | int == 1
  x <- x[ind]
  attr(x, "weights") <- betahat[ind]
  if (trace > 1) print(paste("RSS:", RSS))
  attr(x, "RSS") <- RSS
  x
}


#' @rdname designTree
#' @export
nnls.networx <- function(x, dm, eps = 1e-8) {
  #    spl <- attr(x, "splits")
  spl <- x$splits
  spl2 <- nnls.splits(spl, dm, eps=eps)
  weight <- attr(spl, "weight")
  weight[] <- 0
  weight[match(spl2, spl)] <- attr(spl2, "weight")
  #    attr(attr(x, "splits"), "weight") <- weight
  attr(x$splits, "weight") <- weight
  x$edge.length <- weight[x$splitIndex]
  x
}


#' @rdname designTree
#' @export
designSplits <- function(x, splits = "all", ...)
{
  if (!is.na(pmatch(splits, "all")))
    splits <- "all"
  if (inherits(x, "splits")) return(designUnrooted(x))
  SPLITS <- c("all", "star") # ,"caterpillar")
  splits <- pmatch(splits, SPLITS)
  if (is.na(splits)) stop("invalid splits method")
  if (splits == -1) stop("ambiguous splits method")
  if (splits == 1) X <-  designAll(x)
  if (splits == 2) X <-  designStar(x, ...)
  return(X)
}

# add return splits=FALSE
designAll <- function(n, add.split = FALSE) {
  Y <- matrix(0L, n * (n - 1) / 2, n)
  k <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      Y[k, c(i, j)] <- 1L
      k <- k + 1L
    }
  }
  m <- n - 1L
  X <- matrix(0L, m + 1, 2^m)
  for (i in 1:m)
    X[i, ] <- rep(rep(c(0L, 1L), each = 2^(i - 1)), 2^(m - i))
  X <- X[, -1]
  if (!add.split) return((Y %*% X) %% 2)
  list(X = (Y %*% X) %% 2, Splits = t(X))
}


# faster sparse version
designStar <- function(n, sparse = TRUE) {
  #    res=NULL
  #    for(i in 1:(n-1)) res = rbind(res,cbind(matrix(0,(n-i),i-1),1,diag(n-i)))
  res <- stree(n) |> as.splits() |> splits2design()
  if (!sparse) return(as.matrix(res))
  res
}
