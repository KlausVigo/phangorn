dec2Bin <- function(x) {
  res <- NULL
  i <- 1L
  while (x > 0) {
    if (x %% 2L)
      res <- c(res, i)
    x <- x %/% 2L
    i <- i + 1L
  }
  res
}


# returns binary (0, 1) vector of length k
dec2bin <- function(x, k = ceiling(log2(x))) {
  i <- 1L
  res <- integer(k)
  while (x > 0) {
    if (x %% 2L)
      res[i] <- 1L
    x <- x %/% 2L
    i <- i + 1L
  }
  res
}

# double factorial: log version
#' @rdname dfactorial
#' @export
ldfactorial <- function(x) {
  x <- (x + 1) / 2
  res <- lgamma(2 * x) - (lgamma(x) + (x - 1) * log(2))
  res
}

# double factorial


#' Arithmetic Operators
#'
#' double factorial function
#'
#'
#' @param x a numeric scalar or vector
#' @return \code{dfactorial(x)} returns the double factorial, that is
#' \eqn{x\!\! = 1 * 3 * 5 * \ldots * x } and \code{ldfactorial(x)} is the
#' natural logarithm of it.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[base:Special]{factorial}},
#' \code{\link[ape]{howmanytrees}}
#' @keywords classif
#' @examples
#'
#' dfactorial(1:10)
#'
#' @rdname dfactorial
#' @export dfactorial
"dfactorial" <- function(x) exp(ldfactorial(x))


#
# Hadamard Conjugation
#


### @aliases hadamard fhm h4st h2st

#' Hadamard Matrices and Fast Hadamard Multiplication
#'
#' A collection of functions to perform Hadamard conjugation.  %Hv of a
#' Hadamard matrix H with a vector v using fast Hadamard multiplication.
#'
#' \code{h2st} and \code{h4st} perform Hadamard conjugation for 2-state
#' (binary, RY-coded) or 4-state (DNA/RNA) data. \code{write.nexus.splits}
#' writes splits returned from \code{h2st} or
#' \code{\link[phangorn]{distanceHadamard}} to a nexus file, which can be
#' processed by Spectronet or Splitstree.
#'
#' @param x a vector of length \eqn{2^n}, where n is an integer.
#' @param v a vector of length \eqn{2^n}, where n is an integer.
#' @param obj a data.frame or character matrix, typical a sequence alignment.
#' @param eps Threshold value for splits.
#' @param levels levels of the sequences.
#' @return \code{hadamard} returns a Hadamard matrix. \code{fhm} returns the
#' fast Hadamard multiplication.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{distanceHadamard}}, \code{\link{lento}},
#' \code{\link{plot.networx}}
#' @references Hendy, M.D. (1989). The relationship between simple evolutionary
#' tree models and observable sequence data. \emph{Systematic Zoology},
#' \bold{38} 310--321.
#'
#' Hendy, M. D. and Penny, D. (1993). Spectral Analysis of Phylogenetic Data.
#' \emph{Journal of Classification}, \bold{10}, 5--24.
#'
#' Hendy, M. D. (2005). Hadamard conjugation: an analytical tool for
#' phylogenetics. In O. Gascuel, editor, \emph{Mathematics of evolution and
#' phylogeny}, Oxford University Press, Oxford
#'
#' Waddell P. J. (1995). Statistical methods of phylogenetic analysis:
#' Including hadamard conjugation, LogDet transforms, and maximum likelihood.
#' \emph{PhD thesis}.
#' @keywords cluster
#' @examples
#'
#' H <- hadamard(3)
#' v <- 1:8
#' H %*% v
#' fhm(v)
#'
#' data(yeast)
#'
#' # RY-coding
#' dat_ry <- acgt2ry(yeast)
#' fit2 <- h2st(dat_ry)
#' lento(fit2)
#'
#' # write.nexus.splits(fit2, file = "test.nxs")
#' # read this file into Spectronet or Splitstree to show the network
#' \dontrun{
#' dat <- as.character(yeast)
#' dat4 <- phyDat(dat, type="USER", levels=c("a","c", "g", "t"), ambiguity=NULL)
#' fit4 <- h4st(dat4)
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow=c(3,1))
#' lento(fit4[[1]], main="Transversion")
#' lento(fit4[[2]], main="Transition 1")
#' lento(fit4[[3]], main="Transition 2")
#' par(old.par)
#' }
#'
#' @rdname hadamard
#' @export hadamard
hadamard <- function(x) {
  res <- 1
  while (x > 0) {
    res <- rbind(cbind(res, res), cbind(res, -res))
    x <- x - 1
  }
  res
}


#' @rdname hadamard
#' @export
fhm <- function(v) {
  n <- length(v)
  n <- log2(n)
  res <- .Call("_phangorn_fhm_new", v = as.double(v), n = as.integer(n))
  res
}


seq2split <- function(s) {
  n <- length(s)
  res <- fhm(log(fhm(s))) / n
  res
}


split2seq <- function(q) {
  n <- length(q)
  res <- fhm(exp(fhm(q))) / n
  res
}


#' Distance Hadamard
#'
#' Distance Hadamard produces spectra of splits from a distance matrix.
#'
#'
#' @param dm A distance matrix.
#' @param eps Threshold value for splits.
#' @return \code{distanceHadamard} returns a matrix. The first column contains
#' the distance spectra, the second one the edge-spectra. If eps is positive an
#' object of with all splits greater eps is returned.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}, Tim White
#' @seealso \code{\link{hadamard}}, \code{\link{lento}},
#' \code{\link{plot.networx}}, \code{\link{neighborNet}}
#' @references Hendy, M. D. and Penny, D. (1993). Spectral Analysis of
#' Phylogenetic Data. \emph{Journal of Classification}, \bold{10}, 5-24.
#' @keywords cluster
#' @examples
#'
#' data(yeast)
#' dm <- dist.hamming(yeast)
#' dm <- as.matrix(dm)
#' fit <- distanceHadamard(dm)
#' lento(fit)
#' plot(as.networx(fit), "2D")
#'
#' @export distanceHadamard
distanceHadamard <- function(dm, eps = 0.001) {
  if (inherits(dm, "dist")) {
    n <- attr(dm, "Size")
    Labels <- attr(dm, "Labels")
  }
  if (inherits(dm, "matrix")) {
    n <- dim(dm)[1]
    Labels <- colnames(dm)
    dm <- dm[lower.tri(dm)]
  }
  ns <- 2^(n - 1)
  if (n > 23)
    stop("Hadamard conjugation works only efficient for n < 24")
  result <- .Call("dist2spectra", dm, as.integer(n), as.integer(ns),
    PACKAGE = "phangorn")
  weights <- -fhm(result) / 2^(n - 2)

  if (eps > 0) {
    weights <- weights[-1]
    ind2 <- which(weights > eps)
    n2 <- length(ind2)
    splits <- vector("list", n2)
    for (i in 1:n2) splits[[i]] <- dec2Bin(ind2[i])
    attr(splits, "weights") <- weights[ind2]
    attr(splits, "labels") <- Labels
    attr(splits, "dm") <- dm
    class(splits) <- "splits"
    return(splits)
  }
  res <- data.frame(distance = result, edges = weights, index = 0:(ns - 1))
  attr(res, "Labels") <- Labels
  res
}


#' @rdname hadamard
#' @export
h4st <- function(obj, levels = c("a", "c", "g", "t")) {
  if (is.matrix(obj))
    obj <- as.data.frame(t(obj))
  if (inherits(obj, "phyDat"))
    obj <- as.data.frame(t(as.character(obj)))

  n <- dim(obj)[1]
  p <- dim(obj)[2]

  if (p > 11)
    stop("4-state Hadamard conjugation works only efficient for n < 12")

  DNAX <- matrix(0, n, p)
  DNAY <- matrix(0, n, p)

  DNAX[obj == levels[1]] <- 0
  DNAX[obj == levels[2]] <- 1
  DNAX[obj == levels[3]] <- 1
  DNAX[obj == levels[4]] <- 0

  DNAY[obj == levels[1]] <- 0
  DNAY[obj == levels[2]] <- 1
  DNAY[obj == levels[3]] <- 0
  DNAY[obj == levels[4]] <- 1

  DNAY <- DNAY - DNAY[, p]
  DNAX <- DNAX - DNAX[, p]

  DNAY <- abs(DNAY[, -p])
  DNAX <- abs(DNAX[, -p])
  dy <- DNAY %*% (2^(0:(p - 2)))
  dx <- DNAX %*% (2^(0:(p - 2)))

  INDEX <- dx + 2^(p - 1) * dy
  blub <- table(INDEX)
  index <- as.numeric(rownames(blub)) + 1
  sv <- numeric(4^(p - 1))
  sv[index] <- blub
  qv <- matrix(seq2split(sv), 2^(p - 1), 2^(p - 1))
  sv <- matrix(sv, 2^(p - 1), 2^(p - 1))
  transversion <- transition.1 <- transition.2 <- allSplits(p, colnames(obj))
  attr(transversion, "weights") <- qv[-1, 1]
  attr(transition.1, "weights") <- diag(qv)[-1]
  attr(transition.2, "weights") <- qv[1, -1]
  result <- list(transversion = transversion, transition.1 = transition.1,
    transition.2 = transition.2, qv = qv, sv = sv, n = sum(sv),
    names = names(obj))
  result
}


#' @rdname hadamard
#' @export
h2st <- function(obj, eps = 0.001) {
  if (!inherits(obj, "phyDat")) stop("Error")
  if (attr(obj, "nc") != 2) stop("Error")
  nr <- attr(obj, "nr") # n
  p <- length(obj) # p
  weight <- attr(obj, "weight")
  if (p > 23)
    stop("Hadamard conjugation works only efficient for n < 24")
  DNAX <- matrix(0, nr, p - 1)
  for (i in 1:(p - 1)) DNAX[, i] <- obj[[i]] - 1
  DNAX[obj[[p]] == 2, ] <- 1 - DNAX[obj[[p]] == 2, ]

  index <- DNAX %*% (2^(0:(p - 2))) + 1
  sv <- numeric(2^(p - 1))
  for (i in 1:nr) sv[index[i]] <- sv[index[i]] + weight[i]
  qv <- seq2split(sv)

  if (eps > 0) {
    qv <- qv[-1]
    ind2 <- which(qv > eps)
    indT <- c(2L^(0:(p - 2)), 2L^(p - 1) - 1)
    ind2 <- union(ind2, indT)
    n2 <- length(ind2)
    splits <- vector("list", n2)
    for (i in 1:n2) splits[[i]] <- dec2Bin(ind2[i])
    attr(splits, "weights") <- qv[ind2]
    attr(splits, "labels") <- names(obj)
    class(splits) <- "splits"
    return(splits)
  }
  result <- data.frame(edges = qv, splits = sv, index = 0:(2^(p - 1) - 1))
  attr(result, "Labels") <- names(obj)
  result
}
