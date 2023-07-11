#' Phylogenetic Network
#'
#' \code{splitsNetwork} estimates weights for a splits graph from a distance
#' matrix.
#'
#' \code{splitsNetwork} fits non-negative least-squares phylogenetic networks
#' using L1 (LASSO), L2(ridge regression) constraints.  The function minimizes
#' the penalized least squares
#' \deqn{\beta = min \sum(dm - X\beta)^2 + \lambda \|\beta \|^2_2 }{ beta = sum(dm - X*beta)^2 + lambda |beta|^2_2 }
#' with respect to \deqn{\|\beta \|_1 <= \gamma, \beta >= 0}{ |beta|_1 = gamma, beta >= 0}
#' where \eqn{X} is a design matrix constructed with \code{designSplits}.
#' External edges are fitted without L1 or L2 constraints.
#'
#' @param dm A distance matrix.
#' @param splits a splits object, containing all splits to consider, otherwise
#' all possible splits are used
#' @param gamma penalty value for the L1 constraint.
#' @param lambda penalty value for the L2 constraint.
#' @param weight a vector of weights.
#' @return \code{splitsNetwork} returns a splits object with a matrix added.
#' The first column contains the indices of the splits, the second column an
#' unconstrained fit without penalty terms and the third column the constrained
#' fit.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[phangorn]{distanceHadamard}},
#' \code{\link[phangorn]{designTree}} \code{\link[phangorn]{consensusNet}},
#' \code{\link[phangorn]{plot.networx}}
#' @references Efron, Hastie, Johnstone and Tibshirani (2004) Least Angle
#' Regression (with discussion) \emph{Annals of Statistics} \bold{32(2)},
#' 407--499
#'
#' K. P. Schliep (2009). Some Applications of statistical phylogenetics (PhD
#' Thesis)
#' @keywords cluster
#' @importFrom Matrix sparseMatrix
#' @importFrom quadprog solve.QP
#' @examples
#'
#' data(yeast)
#' dm <- dist.ml(yeast)
#' fit <- splitsNetwork(dm)
#' net <- as.networx(fit)
#' plot(net)
#' write.nexus.splits(fit)
#'
#' @export splitsNetwork
splitsNetwork <- function(dm, splits = NULL, gamma = .1, lambda = 1e-6,
                          weight = NULL) {
  dm <- as.matrix(dm)
  k <- dim(dm)[1]

  if (!is.null(splits)) {
    tmp <- which(lengths(splits) == k)
    if(length(tmp)>0) splits <- splits[-tmp]
    lab <- attr(splits, "labels")
    dm <- dm[lab, lab]
  }

  if (is.null(splits)) {
    X2 <- designAll(k, TRUE)
    X <- X2[[1]]
  }
  else X <- as.matrix(splits2design(splits))

  y <- dm[lower.tri(dm)]
  if (is.null(splits)) ind <- c(2^(0:(k - 2)), 2^(k - 1) - 1)
  else ind <- which(lengths(splits) == 1)
  #  else ind = which(sapply(splits, length)==1)
  #   y2 = lm(y~X[,ind]-1)$res
  n <- dim(X)[2]

  ridge <- lambda * diag(n)
  ridge[ind, ind] <- 0
  if (!is.null(weight)) Dmat <- crossprod(X * sqrt(weight)) + ridge
  else Dmat <- crossprod(X) + ridge
  if (!is.null(weight)) dvec <- crossprod(X * sqrt(weight), y * sqrt(weight))
  else dvec <- crossprod(X, y)

  ind1       <- rep(1, n)
  ind1[ind]  <- 0

  Amat       <- cbind(ind1, diag(n))
  bvec       <- c(gamma, rep(0, n))

  # needs quadprog::solve.QP.compact
  solution <- quadprog::solve.QP(Dmat, dvec, Amat, bvec = bvec, meq = 1)$sol

  ind2 <- which(solution > 1e-8)
  n2 <- length(ind2)

  ind3 <- which(duplicated(c(ind2, ind), fromLast = TRUE)[1:n2])
  ridge2 <- lambda * diag(n2)
  ridge2[ind3, ind3] <- 0

  if (!is.null(weight)) Dmat <- crossprod(X[, ind2] * sqrt(weight)) + ridge2
  else Dmat <- crossprod(X[, ind2]) + ridge2
  if (!is.null(weight)) dvec <- crossprod(X[, ind2] * sqrt(weight),
                                          y * sqrt(weight))
  else dvec <- crossprod(X[, ind2], y)

  Amat2 <- diag(n2)
  bvec2 <- rep(0, n2)
  # needs quadprog::solve.QP.compact
  # bvec2 not used
  solution2  <- quadprog::solve.QP(Dmat, dvec, Amat2)$sol

  RSS1 <- sum( (y - X[, ind2] %*% solution[ind2])^2)
  RSS2 <- sum( (y - X[, ind2] %*% solution2)^2)

  if (is.null(splits)) {
    splits <- vector("list", length(ind2))
    for (i in seq_along(ind2)) splits[[i]] <- which(X2[[2]][ind2[i], ] == 1)
  }
  else splits <- splits[ind2]
  attr(splits, "weights") <- solution[ind2]
  attr(splits, "unrestricted") <- solution2
  attr(splits, "stats") <- c(df = n2, RSS_p = RSS1, RSS_u = RSS2)
  attr(splits, "labels") <- dimnames(dm)[[1]]
  class(splits) <- "splits"
  return(splits)
}
