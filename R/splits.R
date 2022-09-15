#' Splits representation of graphs and trees.
#'
#' \code{as.splits} produces a list of splits or bipartitions.
#'
#' @aliases splits as.Matrix distinct.splits as.phylo.splits
#' addTrivialSplits removeTrivialSplits matchSplits
#' @param x An object of class phylo or multiPhylo.
#' @param maxp integer, default from \code{options(max.print)}, influences how
#' many entries of large matrices are printed at all.
#' @param zero.print character which should be printed for zeros.
#' @param one.print character which should be printed for ones.
#' @param incomparables	only for compatibility so far.
#' @param unrooted todo.
#' @param \dots Further arguments passed to or from other methods.
#' @param recursive	logical. If recursive = TRUE, the function recursively
#' descends through lists (and pairlists) combining all their elements into a
#' vector.
#' @param obj1,obj2 an object of class splits.
#' @param k number of taxa.
#' @param labels names of taxa.
#' @return \code{as.splits} returns an object of class splits, which is mainly
#' a list of splits and some attributes. Often a \code{splits} object will
#' contain attributes \code{confidences} for bootstrap or Bayesian support
#' values and \code{weight} storing edge weights.
#' \code{compatible} return a lower triangular matrix where an 1 indicates that
#' two splits are incompatible.
#' @note The internal representation is likely to change.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{prop.part}}, \code{\link{lento}},
#' \code{\link{as.networx}}, \code{\link{distanceHadamard}},
#' \code{\link{read.nexus.splits}}
#' @keywords cluster
#' @examples
#'
#' (sp <- as.splits(rtree(5)))
#' write.nexus.splits(sp)
#' spl <- allCircularSplits(5)
#' plot(as.networx(spl))
#'
#' @rdname as.splits
#' @export
as.splits <- function(x, ...) {
  if (inherits(x, "splits")) return(x)
  UseMethod("as.splits")
}


#' @export
as.Matrix <- function(x, ...) {
  if (inherits(x, "Matrix")) return(x)
  UseMethod("as.Matrix")
}


#' @rdname as.splits
#' @method as.matrix splits
#' @export
as.matrix.splits <- function(x, zero.print = 0L, one.print = 1L, ...) {
  m <- length(x)
  labels <- attr(x, "labels")
  n <- length(labels)
  res <- matrix(zero.print, m, n)
  for (i in 1:m) res[i, x[[i]]] <- one.print
  dimnames(res) <- list(names(x), labels)
  res
}


#' @rdname as.splits
#' @importFrom Matrix sparseMatrix
#' @method as.Matrix splits
#' @export
as.Matrix.splits <- function(x, ...) {
  labels <- attr(x, "labels")
  l <- length(x)
  j <- unlist(x)
  i <- rep(1:l, lengths(x))
  sparseMatrix(i, j, x = rep(1L, length(i)), dimnames = list(NULL, labels))
  # included x und labels
}


#' @rdname as.splits
#' @export
print.splits <- function(x, maxp = getOption("max.print"),
                         zero.print = ".", one.print = "|", ...) {
  x.orig <- x
  cx <- as.matrix(x, zero.print = zero.print, one.print = one.print)
  print(cx, quote = FALSE, right = TRUE, max = maxp)
  invisible(x.orig)
}


#' @export
"[.splits" <- function(x, i) {
  tmp <- attributes(x)
  result <- unclass(x)[i]
  if (!is.null(tmp$weights)) tmp$weights <- tmp$weights[i]
  if (!is.null(tmp$confidences)) tmp$confidences <- tmp$confidences[i]
  if (!is.null(tmp$intervals)) tmp$intervals <- tmp$intervals[i]
  if (!is.null(tmp$data)) tmp$data <- tmp$data[i, , drop = FALSE]
  attributes(result) <- tmp
  result
}


changeOrder <- function(x, labels) {
  oldL <- attr(x, "labels")
  ind <- match(oldL, labels)
  for (i in seq_along(x))
    x[[i]] <- sort(ind[x[[i]]])
  if (!is.null(attr(x, "cycle")))
    attr(x, "cycle") <- ind[attr(x, "cycle")]
  attr(x, "labels") <- labels
  x
}



## @rdname as.splits
#' @export
matchSplits <- function(x, y, as.in = TRUE) {
  tiplabel <- attr(x, "labels")
  if (any(is.na(match(tiplabel, attr(y, "labels")))))
    stop("x and y have different labels!")
  nTips <- length(tiplabel)
  y <- changeOrder(y, tiplabel)
  y <- SHORTwise(y) #, nTips)
  if (as.in) return(match(SHORTwise(x), y, nomatch = 0L) > 0L)
  match(SHORTwise(x), y)
}



countCycles <- function(splits, ord = NULL) {
  M <- as.matrix(splits)
  if(is.null(ord)){
    ord <- attr(splits, "cycle")
    if(is.null(ord)) ord <- seq_along(attr(splits, "labels"))
  }
  res <- countCycle2_cpp(M[, ord])
  res
}


#' @rdname as.splits
#' @method c splits
#' @export
c.splits <- function(..., recursive = FALSE) {
  x <- list(...)
  if (length(x) == 1 && !inherits(x[[1]], "splits")) x <- x[[1]]
  n <- length(x)
  match.names <- function(a, b) {
    if (any(!(a %in% b)))
      stop("names do not match previous names")
  }
  if (n == 1)
    return(x[[1]])

  labels <- attr(x[[1]], "labels")
  cycle <- attr(x[[1]], "cycle")
  for (i in 2:n) {
    match.names(labels, attr(x[[i]], "labels"))
    x[[i]] <- changeOrder(x[[i]], labels)
  }
  w <- as.vector(unlist(lapply(x, attr, "weights")))
  x <- lapply(x, unclass)
  res <- structure(do.call("c", x), class = c("splits", "prop.part"))
  names(res) <- NULL
  attr(res, "labels") <- labels
  attr(res, "weights") <- w
  attr(res, "cycle") <- cycle
  res
}


#' @rdname as.splits
#' @method unique splits
#' @export
unique.splits <- function(x, incomparables = FALSE, unrooted = TRUE, ...) {
  nTips <- length(attr(x, "labels"))
  x <- SHORTwise(x)
  x[!duplicated(x)]
}


#' @export
distinct.splits <- function(...) {
  tmp <- c(...)
  res <- unique(tmp)
  attributes(res) <-  attributes(tmp)
  attr(res, "weights") <- tabulate(match(tmp, res))
  res
}



# computes splits from phylo
#' @rdname as.splits
#' @method as.splits phylo
#' @export
as.splits.phylo <- function(x, ...) {
  if (hasArg(as.is))
    as.is <- list(...)$as.is
  else as.is <- TRUE
  result <- bip(x)
  if (!is.null(x$edge.length)) {
    edge.weights <- numeric(max(x$edge))
    edge.weights[x$edge[, 2]] <- x$edge.length
    attr(result, "weights") <- edge.weights
  }
  if (!is.null(x$node.label)) {
    conf <- x$node.label
    if (is.character(conf)) as.is <- TRUE
      #conf <- as.numeric(conf)
    if (!as.is) if (max(na.omit(conf)) > (1 + 1e-8)) conf <- conf / 100
    attr(result, "confidences") <- c(rep(NA_real_, length(x$tip.label)), conf)
  }
  attr(result, "labels") <- x$tip.label
  class(result) <- c("splits", "prop.part")
  result
}


# computes splits from multiPhylo object (e.g. bootstrap, MCMC etc.)
#' @rdname as.splits
#' @method as.splits multiPhylo
#' @export
as.splits.multiPhylo <- function(x, ...) {
  if (hasArg(trivial))
    trivial <- list(...)$trivial
  else trivial <- TRUE
  lx <-  length(x)
  x <- unroot(x)
  splits <- prop.part(x)
  splits <- postprocess.prop.part(splits, method="SHORTwise")
  class(splits) <- "list"
  weights <- attr(splits, "number")
  lab <- attr(splits, "labels")
  attr(splits, "labels") <- attr(splits, "number") <- NULL
  l <- length(lab)
  if(trivial){
    splitTips <- vector("list", l)
    for (i in 1:l) splitTips[[i]] <- i
    result <- c(splitTips, splits)
    attr(result, "weights") <- c(rep(lx, l), weights)
  }
  else attr(result, "weights") <- weights
  attr(result, "confidences") <- attr(result, "weights") / lx
  attr(result, "summary") <- list(confidences = "ratio", ntrees = lx,
                                  clades = FALSE)
  attr(result, "labels") <- lab
  class(result) <- c("splits", "prop.part")
  result
}


#' @export
as.splits.prop.part <- function(x, ...) {
  if (is.null(attr(x, "number")))
    attr(x, "weights") <- rep(1, length(x))
  else {
    attr(x, "weights") <- attr(x, "number")
    if( is.integer(attr(x, "number")) )
      attr(x, "confidences") <- attr(x, "number") / attr(x, "number")[1]
  }
  class(x) <- c("splits", "prop.part")
  x
}


#' @rdname as.splits
#' @method as.splits networx
#' @export
as.splits.networx <- function(x, ...) {
  if (!is.null(x$splits)) x$splits
  else warning("No split object included!")
}


#' @rdname as.splits
#' @method as.prop.part splits
#' @export
as.prop.part.splits <- function(x, ...) {
  attr(x, "number") <- attr(x, "weights")
  attr(x, "weights") <- NULL
  attr(x, "confidences") <- NULL
  class(x) <- c("prop.part")
  x
}

## as.splits.phylo
## @rdname as.splits
## @method as.phylo splits
#' @export
as.phylo.splits <- function(x, check=TRUE,...){
  if(check) x <- compatibleSplits(x)
  phy <- splits2phylo(x)
  spl <- as.splits(phy)[phy$edge[,2]]
  ind <- matchSplits(spl, x, FALSE)
  phy$edge.length <- attr(x, "weights")[ind]
  phy
}


splits2phylo <- function(x){
  labels <- attr(x, "labels")
  nTips <- length(labels)
  x <- SHORTwise(x)
  l <- lengths(x)
  x <- x[order(l)]
  x <- x[lengths(x) > 1]
  nNodes <- length(x) + 1L
  node_i <- as.integer( nNodes + nTips )
  edge <- matrix(0L, node_i - 1L, 2L)
  y <- seq_len(nTips)
  x <- unclass(x)
  m <- 0
  for(i in seq_along(x)){
    tmp <- x[[i]]
    kids <- unique( y[ tmp ] )
    k <- length(kids)
    y[ tmp ] <- node_i
    edge[m + 1:k ,1] <- node_i
    edge[m + 1:k ,2] <- kids
    m <- m + k
    node_i <- node_i - 1L
  }
  kids <- unique( y )
  k <- length(kids)
  edge[m + 1:k ,1] <- node_i
  edge[m + 1:k ,2] <- kids
  phy <- structure(list(edge, labels, nNodes),
                   .Names = c("edge", "tip.label", "Nnode"),
                   class = "phylo", order = "postorder")
  phy
}

compatibleSplits <- function(x) {
  x <- postprocess.splits(x)
  labels <- attr(x, "labels")
  nTips <- length(labels)
  x <- SHORTwise(x)
#  x <- x[lengths(x)>1]
  dm <- as.matrix(compatible(x))
  rs <- rowSums(dm)
  ind <- which(rs == 0)
  if (any(rs > 0)) {
    tmp <- which(rs > 0)
    candidates <- tmp[order(rs[tmp])]
    for (i in candidates) {
      if (sum(dm[ind, i]) == 0)
        ind <- c(ind, i)
    }
  }
  x[ind]
}


postprocess.splits <- function (x)
{
  #  w <- attr(x, "number")
  tmp <- attributes(x)
  labels <- attr(x, "labels")
  x <- SHORTwise(x)
  drop <- duplicated(x)
  if (any(drop)) {
    W <- ifelse (is.null(tmp$weights), FALSE, TRUE)
    CONF <- ifelse (is.null(tmp$confidences), FALSE, TRUE)
    if(W) w <- tmp$weights
    if (CONF) conf <- tmp$confidences
    class(x) <- NULL
    attributes(x) <- NULL
    y <- x[drop]
    ind1 <- match(y, x)
    ind2 <- which(drop)
    for (i in seq_along(ind2)) {
      if(W) w[ind1[i]] <- w[ind1[i]] + w[ind2[i]]
      if(CONF) conf[ind1[i]] <- conf[ind1[i]] + conf[ind2[i]]
    }
    x <- x[!drop]
    w <- w[!drop]
    if(CONF) conf <- conf[!drop]
    attr(x, "weights") <- w
    if(CONF) attr(x, "confidences") <- conf
    attr(x, "labels") <- labels
    class(x) <- c("splits", "prop.part")
  }
  x
}


#' @rdname as.splits
#' @method as.bitsplits splits
#' @export
as.bitsplits.splits <- function(x) {
  foo <- function(vect, RAWVECT) {
    res <- RAWVECT
    for (y in vect) {
      i <- ceiling(y / 8)
      res[i] <- res[i] | as.raw(2^(8 - ((y - 1) %% 8) - 1))
    }
    res
  }
  N <- length(x)
  n <- length(attr(x, "labels"))
  nr <- ceiling(n / 8)
  mat <- raw(N * nr)
  dim(mat) <- c(nr, N)
  RAWVECT <- raw(nr)
  for (i in 1:N) mat[, i] <- foo(x[[i]], RAWVECT)
  freq <- attr(x, "weights")
  if (is.null(freq)) freq <- rep(1, N)
  structure(list(matsplit = mat, labels = attr(x, "labels"),
    freq = freq), class = "bitsplits")
}


#' @rdname as.splits
#' @method as.splits bitsplits
#' @export
as.splits.bitsplits <- function(x, ...){
  as.splits(as.prop.part(x))
}


# computes compatible splits
#' @rdname as.splits
#' @export
compatible <- function(obj1, obj2 = NULL) {
  if (!inherits(obj1, "splits"))
    stop("obj needs to be of class splits")
  labels <- attr(obj1, "labels")
  l <- length(labels)
  n <- length(obj1)
  bp1 <- as.matrix(obj1)
  bp1[bp1[, 1] == 0L, ] <- 1L - bp1[bp1[, 1] == 0L, ]
  if (!is.null(obj2)) {
    m <- length(obj2)
    bp2 <- as.matrix(obj2)
    labels2 <- attr(obj2, "labels")
    bp2 <- bp2[, match(labels2, labels), drop = FALSE]
    bp2[bp2[, 1] == 0L, ] <- 1L - bp2[bp2[, 1] == 0L, ]
  }
  else bp2 <- bp1

  if (is.null(obj2)) res <- matrix(0L, n, n)
  else res <- matrix(0L, n, m)

  tmp1 <- tcrossprod(bp1, bp2)
  tmp2 <- tcrossprod(1L - bp1, 1L - bp2)
  tmp3 <- tcrossprod(bp1, 1L - bp2)
  tmp4 <- tcrossprod(1L - bp1, bp2)
  res[(tmp1 * tmp2 * tmp3 * tmp4) > 0] <- 1L
  if (is.null(obj2)) {
    res <- res[lower.tri(res)]
    attr(res, "Size") <- n
    attr(res, "Diag") <- FALSE
    attr(res, "Upper") <- FALSE
    class(res) <- "dist"
  }
  return(res)
}

# in clanistic.R ??
compatible3 <- function(x, y = NULL) {
  if (!inherits(x, "splits"))
    stop("x needs to be of class splits")
  if (is.null(y)) y <- x

  if (!inherits(y, "splits"))
    stop("y needs to be of class splits")
  xlabels <- attr(x, "labels")
  ylabels <- attr(y, "labels")
  if (identical(xlabels, ylabels)) labels <- xlabels
  else labels <- intersect(xlabels, ylabels)
  nx <- length(x)
  ny <- length(y)
  bp1 <- as.matrix(x)[, labels, drop = FALSE]
  bp2 <- as.matrix(y)[, labels, drop = FALSE]
  rs1 <- rowSums(bp1)
  rs2 <- rowSums(bp2)
  res <- matrix(0L, nx, ny)
  tmp1 <- tcrossprod(bp1, bp2)
  res <- matrix(0L, nx, ny)
  for (i in 1:nx) {
    for (j in 1:ny) {
      if (tmp1[i, j] == rs1[i]) res[i, j] <- 1
      if (tmp1[i, j] == rs2[j]) res[i, j] <- 2
      if (tmp1[i, j] == rs1[i] & tmp1[i, j] == rs2[j]) res[i, j] <- 3
    }
  }
  if (is.null(y)) {
    res <- res[lower.tri(res)]
    attr(res, "Size") <- length(x)
    attr(res, "Diag") <- FALSE
    attr(res, "Upper") <- FALSE
    class(res) <- "dist"
  }
  return(res)
}



compatible_2 <- function(obj1, obj2) {
  ntaxa <- length(obj1$labels)
  m1 <- obj1$matsplit
  m2 <- obj2$matsplit
  n1 <- ncol(m1)
  n2 <- ncol(m2)
  res <- rep(TRUE, n1)
  for (i in 1:n1) {
    j <- 1
    while (j <= n2) {
      if (!ape::arecompatible(m1[, i], m2[, j], ntaxa)) {
        res[i] <- FALSE
        break()
      }
      j <- j + 1L
    }
  }
  res
}
