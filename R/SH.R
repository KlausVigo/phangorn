#' Shimodaira-Hasegawa Test
#'
#' This function computes the Shimodaira--Hasegawa test for a set of trees.
#'
#'
#' @param ... either a series of objects of class \code{"pml"} separated by
#' commas, a list containing such objects or an object of class \code{"pmlPart"}
#' or a matrix containing the site-wise likelihoods in columns.
#' @param B the number of bootstrap replicates.
#' @param data an object of class \code{"phyDat"}.
#' @param weight if a matrix with site (log-)likelihoods is is supplied an
#' optional vector containing the number of occurrences of each site pattern.
#' @return a numeric vector with the P-value associated with each tree given in
#' \code{...}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}}, \code{\link{pmlPart}}, \code{\link{pmlCluster}},
#' \code{\link{SOWH.test}}
#' @references Shimodaira, H. and Hasegawa, M. (1999) Multiple comparisons of
#' log-likelihoods with applications to phylogenetic inference. \emph{Molecular
#' Biology and Evolution}, \bold{16}, 1114--1116.
#' @keywords models
#' @examples
#'
#' data(Laurasiatherian)
#' dm <- dist.logDet(Laurasiatherian)
#' tree1 <- NJ(dm)
#' tree2 <- unroot(upgma(dm))
#' fit1 <- pml(tree1, Laurasiatherian)
#' fit2 <- pml(tree2, Laurasiatherian)
#' fit1 <- optim.pml(fit1) # optimize edge weights
#' fit2 <- optim.pml(fit2)
#' # with pml objects as input
#' SH.test(fit1, fit2, B=1000)
#' # in real analysis use larger B, e.g. 10000
#'
#' # with matrix as input
#' X <- matrix(c(fit1$siteLik, fit2$siteLik), ncol=2)
#' SH.test(X, weight=attr(Laurasiatherian, "weight"), B=1000)
#' \dontrun{
#' example(pmlPart)
#' SH.test(sp, B=1000)
#' }
#' @export
SH.test <- function(..., B = 10000, data = NULL, weight = NULL) {
  fits <- list(...)
  if (inherits(fits[[1]], "matrix") || inherits(fits[[1]], "data.frame"))
    return(SH.tmp(fits[[1]], weight = weight, B = B))
  p <- 1
  if (inherits(fits[[1]], "pmlPart")) {
    fits <- fits[[1]]$fits
    p <- length(fits)
  }
  if (inherits(fits[[1]], "list")) {
    fits <- fits[[1]]
  }
  k <- length(fits)
  if (is.null(data))
    data <- fits[[1]]$data
  res <- NULL
  for (h in 1:p) {
    if (p > 1)
      data <- fits[[h]]$data
    weight <- attr(data, "weight")
    lw <- length(weight)
    siteLik <- matrix(0, lw, k)
    for (i in 1:k) siteLik[, i] <- update(fits[[i]], data = data)$siteLik
    ntree <- k
    Lalpha <- drop(crossprod(siteLik, weight))
    Talpha <- max(Lalpha) - Lalpha
    M <- matrix(NA, k, B)
    wvec <- rep(1L:lw, weight)
    for (i in 1:B) {
      boot <- tabulate(sample(wvec, replace = TRUE), nbins = lw)
      M[, i] <- crossprod(siteLik, boot)
    }
    M <- M - rowMeans(M)
    #   for (i in 1:B) for (j in 1:ntree) S[j, i] <- max(M[j, i] - M[, i])
    S <- matrix(apply(M, 2, min), k, B, byrow = TRUE)
    S <- M - S
    count <- numeric(ntree)
    for (j in 1:ntree) count[j] <- sum(S[j, ] > Talpha[j])
    count <- count / B
    trees <- 1:k
    if (p == 1)
      res <- cbind(trees, Lalpha, Talpha, count)
    else res <- rbind(res, cbind(h, trees[-h], Lalpha[-h],
        Talpha[-h], count[-h]))
  }
  if (p == 1)
    colnames(res) <- c("Trees", "ln L", "Diff ln L", "p-value")
  else colnames(res) <- c("Partition", "Trees", "ln L", "Diff ln L",
      "p-value")
  res
}


SH.tmp <- function(siteLik, weight = NULL, B = 10000, boot=NULL) {
  lw <- nrow(siteLik)
  if (is.null(weight)) weight <- rep(1, lw)
  ntree <- k <- ncol(siteLik)
  Lalpha <- drop(crossprod(siteLik, weight))
  Talpha <- max(Lalpha) - Lalpha
  if(is.null(boot)){
    M <- matrix(NA, k, B)
    wvec <- rep(1L:lw, weight)
    for (i in 1:B) {
      boot <- tabulate(sample(wvec, replace = TRUE), nbins = lw)
      M[, i] <- crossprod(siteLik, boot)
    }
  } else {
    M <- crossprod(siteLik, boot)
  }
  M <- M - rowMeans(M)
  S <- matrix(apply(M, 2, min), k, B, byrow = TRUE)
  S <- M - S
  count <- numeric(ntree)
  for (j in 1:ntree) count[j] <- sum(S[j, ] > Talpha[j])
  count <- count / B
  trees <- 1:k
  res <- cbind(trees, Lalpha, Talpha, count)
  colnames(res) <- c("Trees", "ln L", "Diff ln L", "p-value")
  res
}
