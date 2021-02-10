#
# pmlPart + pmlCluster
#
optimPartQGeneral <- function(object, Q = c(1, 1, 1, 1, 1, 1),
                              subs = rep(1, length(Q)), ...) {
  m <- length(Q)
  n <- max(subs)
  ab <- numeric(n)
  for (i in 1:n) ab[i] <- log(Q[which(subs == i)[1]])
  fn <- function(ab, object, m, n, subs, ...) {
    Q <- numeric(m)
    for (i in 1:n) Q[subs == i] <- ab[i]
    Q <- exp(Q)
    result <- 0
    for (i in seq_along(object)) result <- result + update(object[[i]],
        Q = Q, ...)$logLik
    result
  }
  res <- optim(par = ab, fn = fn, gr = NULL, method = "L-BFGS-B",
    lower = -Inf, upper = Inf, control = list(fnscale = -1,
      maxit = 25), object = object, m = m, n = n, subs = subs, ...)
  Q <- rep(1, m)
  for (i in 1:n) Q[subs == i] <- exp(res[[1]][i])
  res[[1]] <- Q
  res
}


optimPartBf <- function(object, bf = c(0.25, 0.25, 0.25, 0.25), ...) {
  l <- length(bf)
  nenner <- 1 / bf[l]
  lbf <- log(bf * nenner)
  lbf <- lbf[-l]
  fn <- function(lbf, object, ...) {
    result <- 0
    bf <- exp(c(lbf, 0))
    bf <- bf / sum(bf)
    n <- length(object)
    for (i in 1:n) result <- result + update(object[[i]],
        bf = bf, ...)$logLik
    result
  }
  res <- optim(par = lbf, fn = fn, gr = NULL, method = "Nelder-Mead",
    control = list(fnscale = -1, maxit = 500), object, ...)
  bf <- exp(c(res[[1]], 0))
  bf <- bf / sum(bf)
}


optimPartInv <- function(object, inv = 0.01, ...) {
  fn <- function(inv, object, ...) {
    result <- 0
    n <- length(object)
    for (i in 1:n) result <- result + update(object[[i]], inv = inv,
        ...)$logLik
    result
  }
  res <- optimize(f = fn, interval = c(0, 1), lower = 0, upper = 1,
    maximum = TRUE, tol = 1e-04, object, ...)
  res[[1]]
}


optimPartGamma <- function(object, shape = 1, ...) {
  fn <- function(shape, object, ...) {
    result <- 0
    n <- length(object)
    for (i in 1:n) result <- result + update(object[[i]], shape = shape,
        ...)$logLik
    result
  }
  res <- optimize(f = fn, interval = c(0, 100), lower = 0, upper = 100,
    maximum = TRUE, tol = 0.01, object, ...)
  res
}


dltmp <- function(fit, i = 1, transform = transform) {
  tree <- fit$tree
  data <- getCols(fit$data, tree$tip.label)
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorder(tree, "postorder")
  q <- length(tree$tip.label)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  m <- max(edge)
  dat <- vector(mode = "list", length = m)
  eig <- fit$eig
  w <- fit$w[i]
  g <- fit$g[i]
  bf <- fit$bf
  el <- tree$edge.length
  P <- getP(el, eig, g)
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  node <- as.integer(node - min(node))
  edge <- as.integer(edge - 1)
  nTips <- as.integer(length(tree$tip.label))
  mNodes <- as.integer(max(node) + 1)
  contrast <- attr(data, "contrast")
  nco <- as.integer(dim(contrast)[1])
  dat[(q + 1):m] <- .Call("LogLik2", data, P, nr, nc, node, edge, nTips,
    mNodes, contrast, nco)

  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  nTips <- min(parent) - 1
  datp <- vector("list", m)
  el <- tree$edge.length
  if (transform) dP <- getdP(tree$edge.length, eig, g)
  else dP <- getdP2(tree$edge.length, eig, g)

  datp[(nTips + 1)] <- dat[(nTips + 1)]
  l <- length(child)
  dl <- matrix(0, nr, l)
  for (j in (m - 1):1) {
    # tips have factor format, internal edges are matrices
    if (child[j] > nTips) {
      tmp2 <- (datp[[parent[j]]] / (dat[[child[j]]] %*% P[[j]]))
      dl[, j] <- (tmp2 * (dat[[child[j]]] %*% dP[[j]])) %*% (w * bf)
      datp[[child[j]]] <- (tmp2 %*% P[[j]]) * dat[[child[j]]]
    }
    else {
      tmp2 <- datp[[parent[j]]] / ((contrast %*% P[[j]])[data[[child[j]]], ])
      dl[, j] <- (tmp2 * ((contrast %*% dP[[j]])[data[[child[j]]], ])) %*%
        (w * bf)
    }
  }
  dl
}


dl <- function(x, transform = TRUE) {
  l <- length(x$w)
  dl <- dltmp(x, 1, transform)
  i <- 2
  while (i < (l + 1)) {
    dl <- dl + dltmp(x, i, transform)
    i <- i + 1
  }
  dl
}


# add control and change edge
optimPartEdge <- function(object, ...) {
  tree <- object[[1]]$tree
  theta <- tree$edge.length
  theta <- pmax(theta, 1e-8)
  tree$edge.length <- theta
  tmptree <- tree
  n <- length(object)
  l <- length(theta)
  nrv <- numeric(n)
  for (i in 1:n) nrv[i] <- attr(object[[i]]$data, "nr")
  cnr <- cumsum(c(0, nrv))
  weight <- numeric(sum(nrv))
  dl <- matrix(NA, sum(nrv), l)
  for (i in 1:n) weight[(cnr[i] + 1):cnr[i + 1]] <- attr(object[[i]]$data,
      "weight")
  ll0 <- 0
  for (i in 1:n) object[[i]] <- update(object[[i]], tree = tree)
  for (i in 1:n) ll0 <- ll0 + object[[i]]$logLik
  eps <- 1
  scalep <- 1
  k <- 1
  while (eps > 0.001 & k < 50) {
    if (scalep == 1) {
      for (i in 1:n) {
        lv <- drop(exp(object[[i]]$siteLik))
        dl[(cnr[i] + 1):cnr[i + 1], ] <- dl(object[[i]], TRUE) / lv
      }
      sc <- colSums(weight * dl)
      F <- crossprod(dl * weight, dl) + diag(l) * 1e-10
      # add small ridge penalty for numerical stability
    }
    thetaNew <- log(theta) + scalep * solve(F, sc)
    thetaNew <- pmax(thetaNew, log(1e-8))
    tmptree$edge.length <- as.numeric(exp(thetaNew))
    for (i in 1:n) object[[i]] <- update(object[[i]], tree = tmptree)
    ll1 <- 0
    for (i in 1:n) ll1 <- ll1 + object[[i]]$logLik
    eps <- ll1 - ll0
    if (eps < 0 || is.nan(eps)) {
      scalep <- scalep / 2
      eps <- 1
      thetaNew <- log(theta)
      ll1 <- ll0
    }
    else {
      scalep <- 1
      tree <- tmptree
    }
    theta <- exp(thetaNew)
    theta <- pmax(theta, 1e-8)
    ll0 <- ll1
    k <- k + 1
  }
  for (i in 1:n) object[[i]] <- update(object[[i]], tree = tree)
  object
}


makePart <- function(fit, rooted, weight = ~index + genes) {
  if (inherits(fit, "phyDat")) {
    x <- fit
    dm <- dist.ml(x)
    if (!rooted) tree <- NJ(dm)
    else tree <- upgma(dm)
    fit <- pml(tree, x, k = 4)
  }
  dat <- fit$data
  if (class(weight)[1] == "formula")
    weight <- xtabs(weight, data = attr(dat, "index"))
  fits <- NULL
  for (i in 1:dim(weight)[2]) {
    ind <- which(weight[, i] > 0)
    dat2 <- getRows(dat, ind)
    attr(dat2, "weight") <- weight[ind, i]
    fits[[i]] <- update(fit, data = dat2)
  }
  names(fits) <- colnames(fits)
  fits
}


#' @rdname pmlPart
#' @export
multiphyDat2pmlPart <- function(x, rooted = FALSE,  ...) {
  shared_tree <- TRUE
  if (shared_tree) {
    concatenate_x <- do.call(cbind.phyDat, x@seq)
    dm <- dist.ml(concatenate_x)
    if (!rooted) tree <- NJ(dm)
    else tree <- upgma(dm)
  }
  else tree <- NULL
  fun <-  function(x, rooted = FALSE, tree, ...) {
    if (is.null(tree)) {
      dm <- dist.ml(x)
      if (!rooted) tree <- NJ(dm)
      else tree <- upgma(dm)
    }
    pml(tree, x, ...)
  }
  fits <- lapply(x@seq, fun, tree = tree, rooted = rooted, ...)
  fits
}


#' @rdname pmlPart
#' @export
pmlPart2multiPhylo <- function(x) {
  res <- lapply(x$fits, FUN = function(x) x$tree)
  class(res) <- "multiPhylo"
  res
}

#' @export
plot.pmlPart <- function(x, ...) {
  plot(pmlPart2multiPhylo(x), ...)
}



#' Partition model.
#'
#' Model to estimate phylogenies for partitioned data.
#'
#' The \code{formula} object allows to specify which parameter get optimized.
#' The formula is generally of the form \code{edge + bf + Q ~ rate + shape +
#' \dots{}}, on the left side are the parameters which get optimized over all
#' partitions, on the right the parameter which are optimized specific to each
#' partition. The parameters available are \code{"nni", "bf", "Q", "inv",
#' "shape", "edge", "rate"}.  Each parameters can be used only once in the
#' formula.  \code{"rate"} is only available for the right side of the formula.
#'
#' For partitions with different edge weights, but same topology, \code{pmlPen}
#' can try to find more parsimonious models (see example).
#'
#' \code{pmlPart2multiPhylo} is a convenience function to extract the trees out
#' of a \code{pmlPart} object.
#'
#' @aliases pmlPart
#' @param formula a formula object (see details).
#' @param object an object of class \code{pml} or a list of objects of class
#' \code{pml} .
#' @param control A list of parameters for controlling the fitting process.
#' @param model A vector containing the models containing a model for each
#' partition.
#' @param rooted Are the gene trees rooted (ultrametric) or unrooted.
#' @param \dots Further arguments passed to or from other methods.
#' @param x an object of class \code{pmlPart}
#' @return \code{kcluster} returns a list with elements
#' \item{logLik}{log-likelihood of the fit} \item{trees}{a list of all trees
#' during the optimization.} \item{object}{an object of class \code{"pml"} or
#' \code{"pmlPart"}}
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso
#' \code{\link{pml}},\code{\link{pmlCluster}},\code{\link{pmlMix}},\code{\link{SH.test}}
#' @keywords cluster
#' @examples
#'
#' data(yeast)
#' dm <- dist.logDet(yeast)
#' tree <- NJ(dm)
#' fit <- pml(tree,yeast)
#' fits <- optim.pml(fit)
#'
#' weight=xtabs(~ index+genes,attr(yeast, "index"))[,1:10]
#'
#' sp <- pmlPart(edge ~ rate + inv, fits, weight=weight)
#' sp
#'
#' \dontrun{
#' sp2 <- pmlPart(~ edge + inv, fits, weight=weight)
#' sp2
#' AIC(sp2)
#'
#' sp3 <- pmlPen(sp2, lambda = 2)
#' AIC(sp3)
#' }
#'
#' @rdname pmlPart
#' @export pmlPart
pmlPart <- function(formula, object, control = pml.control(epsilon = 1e-8,
                    maxit = 10, trace = 1), model = NULL, rooted = FALSE, ...) {
  call <- match.call()
  form <- phangornParseFormula(formula)
  opt <- c("nni", "bf", "Q", "inv", "shape", "edge", "rate")
  optAll <- match(opt, form$left)
  optPart <- match(opt, form$right)
  AllNNI <- !is.na(optAll[1])
  AllBf <- !is.na(optAll[2])
  AllQ <- !is.na(optAll[3])
  AllInv <- !is.na(optAll[4])
  AllGamma <- !is.na(optAll[5])
  AllEdge <- !is.na(optAll[6])
  PartNni <- !is.na(optPart[1])
  PartBf <- !is.na(optPart[2])
  PartQ <- !is.na(optPart[3])
  PartInv <- !is.na(optPart[4])
  PartGamma <- !is.na(optPart[5])
  PartEdge <- !is.na(optPart[6])
  PartRate <- !is.na(optPart[7])

  if (PartNni) PartEdge <- TRUE
  if(AllNNI) AllEdge <- TRUE
  if (inherits(object, "multiphyDat")) {
    if (AllNNI || AllEdge) object <- do.call(cbind.phyDat, object@seq)
    else fits <- multiphyDat2pmlPart(object, rooted = rooted, ...)
  }
  if (inherits(object, "pml")) fits <- makePart(object, rooted = rooted, ...)
  if (inherits(object, "phyDat")) fits <- makePart(object, rooted = rooted, ...)
  if (inherits(object, "pmlPart")) fits <- object$fits
  if (inherits(object, "list")) fits <- object
  if(AllNNI) if(Ntip(fits[[1]]$tree) <  (3 + !rooted)) AllNNI <- FALSE

  trace <- control$trace
  epsilon <- control$epsilon
  maxit <- control$maxit

  p <- length(fits)
  #   if(length(model)<p) model = rep(model, length = p)

  if (AllQ) {
    Q <- fits[[1]]$Q
    for (i in 1:p) fits[[i]] <- update(fits[[i]], Q = Q)
  }
  if (AllBf) {
    bf <- fits[[1]]$bf
    for (i in 1:p) fits[[i]] <- update(fits[[i]], bf = bf)
  }
  if (AllInv) {
    inv <- fits[[1]]$inv
    for (i in 1:p) fits[[i]] <- update(fits[[i]], inv = inv)
  }
  if (AllGamma) {
    shape <- fits[[1]]$shape
    for (i in 1:p) fits[[i]] <- update(fits[[i]], shape = shape)
  }
  if (AllEdge || AllNNI) {
    tree <- fits[[1]]$tree
    tree$edge.length <- pmax(tree$edge.length, 1e-6)
    for (i in 1:p) fits[[i]] <- update(fits[[i]], tree = tree)
  }


  m <- 1
  logLik <- 0
  for (i in 1:p) logLik <- logLik + fits[[i]]$log
  eps <- 10

  on.exit({
    df <- matrix(1, 6, 2)
    colnames(df) <- c("#df", "group")
    rownames(df) <- c("Edge", "Shape", "Inv", "Bf", "Q", "Rate")
    df[1, 1] <- length(fits[[1]]$tree$edge.length)
    df[2, 1] <- fits[[1]]$k > 1
    df[3, 1] <- fits[[1]]$inv > 0
    df[4, 1] <- length(unique(fits[[1]]$bf)) - 1
    df[5, 1] <- length(unique(fits[[1]]$Q)) - 1
    df[6, 1] <- 0 # rates
    if (PartEdge) df[1, 2] <- p
    if (PartGamma) df[2, 2] <- p
    if (PartInv) df[3, 2] <- p
    if (PartBf) df[4, 2] <- p
    if (PartQ) df[5, 2] <- p
    if (PartRate) df[6, 1] <- p - 1
    attr(logLik, "df") <- sum(df[, 1] * df[, 2])
    object <- list(logLik = logLik, fits = fits, call = call, df = df)
    class(object) <- "pmlPart"
    return(object)
  })

  while (eps > epsilon & m < maxit) {
    loli <- 0
    if (any(c(PartNni, PartBf, PartInv, PartQ, PartGamma, PartEdge, PartRate))){
      for (i in 1:p) {
        fits[[i]] <- optim.pml(fits[[i]], optNni = PartNni, optBf = PartBf,
          optQ = PartQ, optInv = PartInv, optGamma = PartGamma,
          optEdge = PartEdge, optRate = PartRate, optRooted = rooted,
          control = pml.control(maxit = 3, epsilon = 1e-8, trace - 1),
          model = model[i])
      }
    }
    if (AllQ) {
      Q <- fits[[1]]$Q
      subs <- c(1:(length(Q) - 1), 0)
      newQ <- optimPartQGeneral(fits, Q = Q, subs = subs)
      for (i in 1:p) fits[[i]] <- update(fits[[i]], Q = newQ[[1]])
    }
    if (AllBf) {
      bf <- fits[[1]]$bf
      newBf <- optimPartBf(fits, bf = bf)
      for (i in 1:p) fits[[i]] <- update(fits[[i]], bf = newBf)
    }
    if (AllInv) {
      inv <- fits[[1]]$inv
      newInv <- optimPartInv(fits, inv = inv)
      for (i in 1:p) fits[[i]] <- update(fits[[i]], inv = newInv)
    }
    if (AllGamma) {
      shape <- fits[[1]]$shape
      newGamma <- optimPartGamma(fits, shape = shape)[[1]]
      for (i in 1:p) fits[[i]] <- update(fits[[i]], shape = newGamma)
    }
    if (AllNNI) {
      fits <- optimPartNNI(fits, AllEdge)
      if (trace > 0) cat(attr(fits, "swap"), " NNI operations performed")
    }
    if (AllEdge)
      fits <- optimPartEdge(fits)
    if (PartRate) {
      tree <- fits[[1]]$tree
      rate <- numeric(p)
      wp <- numeric(p)
      for (i in 1:p) {
        wp[i] <- sum(fits[[i]]$weight)
        rate[i] <- fits[[i]]$rate
      }
      ratemult <- sum(wp) / sum(wp * rate)
      tree$edge.length <- tree$edge.length / ratemult
      for (i in 1:p) fits[[i]] <- update(fits[[i]], tree = tree,
                                         rate = rate[i] * ratemult)
    }
    loli <- 0
    for (i in 1:p) loli <- loli + fits[[i]]$log
    eps <- (logLik - loli) / loli
    if (trace > 0) cat("loglik:", logLik, "-->", loli, "\n")
    logLik <- loli
    m <- m + 1
  }
}


#
# pmlCluster
#


pmlCluster.fit <- function(formula, fit, weight, p = 4, part = NULL,
                           control = pml.control(epsilon = 1e-8, maxit = 10,
                                                 trace = 1), ...) {
  call <- match.call()
  form <- phangornParseFormula(formula)
  opt <- c("nni", "bf", "Q", "inv", "shape", "edge", "rate")
  optAll <- match(opt, form$left)
  optPart <- match(opt, form$right)
  AllNNI <- !is.na(optAll[1])
  AllBf <- !is.na(optAll[2])
  AllQ <- !is.na(optAll[3])
  AllInv <- !is.na(optAll[4])
  AllGamma <- !is.na(optAll[5])
  AllEdge <- !is.na(optAll[6])
  PartNni <- !is.na(optPart[1])
  PartBf <- !is.na(optPart[2])
  PartQ <- !is.na(optPart[3])
  PartInv <- !is.na(optPart[4])
  PartGamma <- !is.na(optPart[5])
  PartEdge <- !is.na(optPart[6])
  PartRate <- !is.na(optPart[7])
  if (PartNni) PartEdge <- TRUE
  nrw <- dim(weight)[1]
  ncw <- dim(weight)[2]
  if (is.null(part)) {
    part <- rep(1:p, length = ncw)
    part <- sample(part)
  }
  Part <- part
  Gtrees <- vector("list", p)
  dat <- fit$data
  attr(fit$orig.data, "index") <- attr(dat, "index") <- NULL
  for (i in 1:p) Gtrees[[i]] <- fit$tree
  fits <- vector("list", p)
  for (i in 1:p) fits[[i]] <- fit
  trace <- control$trace
  eps <- 0
  m <- 1
  logLik <- fit$log
  trees <- list()
  weights <- matrix(0, nrw, p)
  lls <- matrix(0, nrw, p)
  loli <- fit$log
  oldpart <- part
  eps2 <- 1
  iter <- 0
  swap <- 1

  on.exit({
    df <- matrix(1, 6, 2)
    colnames(df) <- c("#df", "group")
    rownames(df) <- c("Edge", "Shape", "Inv", "Bf", "Q", "Rate")
    df[1, 1] <- length(fits[[1]]$tree$edge.length)
    df[2, 1] <- fits[[1]]$k - 1
    df[3, 1] <- fits[[1]]$inv > 0
    df[4, 1] <- length(unique(fits[[1]]$bf)) - 1
    df[5, 1] <- length(unique(fits[[1]]$Q)) - 1
    df[6, 1] <- 0
    if (PartEdge)
      df[1, 2] <- p
    if (PartGamma)
      df[2, 2] <- p
    if (PartInv)
      df[3, 2] <- p
    if (PartBf)
      df[4, 2] <- p
    if (PartQ)
      df[5, 2] <- p
    if (PartRate)
      df[6, 1] <- p - 1
    attr(logLik, "df") <- sum(df[, 1] * df[, 2])
    res <- list(logLik = logLik, Partition = Part, trees = trees)
    result <- list(logLik = loli, fits = fits, Partition = part, df = df,
                   res = res, call = call)
    class(result) <- c("pmlPart")
    return(result)
  })

  while (eps < ncw || abs(eps2) > control$eps) {
    df2 <- 0
    if (any(c(PartNni, PartBf, PartInv, PartQ, PartGamma, PartEdge, PartRate))){
      for (i in 1:p) {
        weights[, i] <- rowSums(weight[, which(part == i),
          drop = FALSE])
        ind <- which(weights[, i] > 0)
        dat2 <- getRows(dat, ind)
        attr(dat2, "weight") <- weights[ind, i]
        fits[[i]] <- update(fits[[i]], data = dat2)
        fits[[i]] <- optim.pml(fits[[i]], PartNni, PartBf,
          PartQ, PartInv, PartGamma, PartEdge, PartRate,
          control = pml.control(epsilon = 1e-8, maxit = 3, trace - 1))
        lls[, i] <- update(fits[[i]], data = dat)$siteLik
        Gtrees[[i]] <- fits[[i]]$tree
      }
    }
    if (AllQ) {
      Q <- fits[[1]]$Q
      subs <- c(1:(length(Q) - 1), 0)
      newQ <- optimPartQGeneral(fits, Q = Q, subs = subs)[[1]]
      for (i in 1:p) fits[[i]] <- update(fits[[i]], Q = newQ)
      df2 <- df2 + length(unique(newQ)) - 1
    }
    if (AllBf) {
      bf <- fits[[1]]$bf
      newBf <- optimPartBf(fits, bf = bf)
      for (i in 1:p) fits[[i]] <- update(fits[[i]], bf = newBf)
      df2 <- df2 + length(unique(newBf)) - 1
    }
    if (AllInv) {
      inv <- fits[[1]]$inv
      newInv <- optimPartInv(fits, inv = inv)
      for (i in 1:p) fits[[i]] <- update(fits[[i]], inv = newInv)
      # there was an Error
      df2 <- df2 + 1
    }
    if (AllGamma) {
      shape <- fits[[1]]$shape
      newGamma <- optimPartGamma(fits, shape = shape)[[1]]
      for (i in 1:p) fits[[i]] <- update(fits[[i]], shape = newGamma)
      df2 <- df2 + 1
    }
    if (AllNNI) {
      fits <- optimPartNNI(fits, AllEdge)
      if (trace > 0) cat(attr(fits, "swap"), " NNI operations performed")
      swap <- attr(fits, "swap")
    }
    if (AllEdge) {
      fits <- optimPartEdge(fits)
      df2 <- df2 + length(fits[[1]]$tree$edge.length)
    }
    if (PartRate) {
      tree <- fits[[1]]$tree
      rate <- numeric(p)
      wp <- numeric(p)
      for (i in 1:p) {
        wp[i] <- sum(fits[[i]]$weight)
        rate[i] <- fits[[i]]$rate
      }
      ratemult <- sum(wp) / sum(wp * rate)
      tree$edge.length <- tree$edge.length / ratemult
      for (i in 1:p) fits[[i]] <- update(fits[[i]], tree = tree,
          rate = rate[i] * ratemult)
    }
    for (i in 1:p) lls[, i] <- update(fits[[i]], data = dat)$siteLik
    trees[[m]] <- Gtrees
    LL <- t(weight) %*% lls
    # choose partitions which change
    tmp <- (LL[cbind(1:ncw, part)] - apply(LL, 1, max)) / colSums(weight)
    fixi <- numeric(p)
    for (i in 1:p) {
      tmpi <- which(part == i)
      fixi[i] <- tmpi[which.max(tmp[tmpi])]
    }
    oldpart <- part
    # restrict the number of elements changing groups
    # If more than 25% would change, only the 25% with the highest increase per
    # site change
    if (sum(tmp == 0) / length(tmp) < .75) {
      medtmp <- quantile(tmp, .25)
      medind <- which(tmp <= medtmp)
      part[medind] <- max.col(LL[medind, ])
    }
    else part <- max.col(LL)
    #        else part <- apply(LL, 1, which.max)
    # force groups to have at least one member
    part[fixi] <- 1:p
    Part <- cbind(Part, part)
    eps <- sum(diag(table(part, oldpart)))
    eps2 <- loli
    loli <- sum(apply(LL, 1, max))
    eps2 <- (eps2 - loli) / loli
    logLik <- c(logLik, loli)
    if (trace > 0) print(loli)
    Part <- cbind(Part, part)
    df2 <- df2 + df2
    if (eps == ncw & swap == 0)
      AllNNI <- FALSE
    m <- m + 1
    if (eps == ncw)
      iter <- iter + 1
    if (iter == 3)
      break
  }
}



#' Stochastic Partitioning
#'
#' Stochastic Partitioning of genes into p cluster.
#'
#' The \code{formula} object allows to specify which parameter get optimized.
#' The formula is generally of the form \code{edge + bf + Q ~ rate + shape +
#' \dots{}}, on the left side are the parameters which get optimized over all
#' cluster, on the right the parameter which are optimized specific to each
#' cluster. The parameters available are \code{"nni", "bf", "Q", "inv",
#' "shape", "edge", "rate"}.  Each parameter can be used only once in the
#' formula.  There are also some restriction on the combinations how parameters
#' can get used. \code{"rate"} is only available for the right side.  When
#' \code{"rate"} is specified on the left hand side \code{"edge"} has to be
#' specified (on either side), if \code{"rate"} is specified on the right hand
#' side it follows directly that \code{edge} is too.
#'
#' @param formula a formula object (see details).
#' @param fit an object of class \code{pml}.
#' @param weight \code{weight} is matrix of frequency of site patterns for all
#' genes.
#' @param p number of clusters.
#' @param part starting partition, otherwise a random partition is generated.
#' @param nrep number of replicates for each p.
#' @param control A list of parameters for controlling the fitting process.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{pmlCluster} returns a list with elements
#' \item{logLik}{log-likelihood of the fit} \item{trees}{a list of all trees
#' during the optimization.} \item{fits}{fits for the final partitions}
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso
#' \code{\link{pml}},\code{\link{pmlPart}},\code{\link{pmlMix}},\code{\link{SH.test}}
#' @references K. P. Schliep (2009). Some Applications of statistical
#' phylogenetics (PhD Thesis)
#'
#' Lanfear, R., Calcott, B., Ho, S.Y.W. and Guindon, S. (2012) PartitionFinder:
#' Combined Selection of Partitioning Schemes and Substitution Models for
#' Phylogenetic Analyses. \emph{Molecular Biology and Evolution}, \bold{29(6)},
#' 1695-1701
#' @keywords cluster
#' @examples
#'
#' \dontrun{
#' data(yeast)
#' dm <- dist.logDet(yeast)
#' tree <- NJ(dm)
#' fit <- pml(tree,yeast)
#' fit <- optim.pml(fit)
#'
#' weight <- xtabs(~ index+genes,attr(yeast, "index"))
#' set.seed(1)
#'
#' sp <- pmlCluster(edge~rate, fit, weight, p=1:4)
#' sp
#' SH.test(sp)
#' }
#'
#' @export pmlCluster
pmlCluster <- function(formula, fit, weight, p = 1:5, part = NULL, nrep = 10,
                       control = pml.control(epsilon = 1e-08, maxit = 10,
                                             trace = 1), ...) {
  call <- match.call()
  form <- phangornParseFormula(formula)
  if (any(p == 1)) {
    opt2 <- c("nni", "bf", "Q", "inv", "shape", "edge")
    tmp1 <- opt2 %in% form$left
    tmp1 <- tmp1 | (opt2 %in% form$right)
    fit <- optim.pml(fit, tmp1[1], tmp1[2], tmp1[3], tmp1[4],
      tmp1[5], tmp1[6], control = control)
  }

  p <- p[p != 1]
  if (length(p) == 0) return(fit)
  n <- sum(weight)
  k <- 2

  BIC <- matrix(0, length(p) + 1, nrep)
  BIC[1, ] <- AIC(fit, k = log(n))
  LL <- matrix(NA, length(p) + 1, nrep)
  LL[1, ] <- logLik(fit)

  P <- array(dim = c(length(p) + 1, nrep, dim(weight)[2]))
  tmpBIC <- Inf
  choice <- c(1, 1)
  for (j in p) {
    tmp <- NULL
    for (i in 1:nrep) {
      tmp <- pmlCluster.fit(formula, fit, weight, p = j, part = part,
        control = control, ...)
      P[k, i, ] <- tmp$Partition
      BIC[k, i] <- AIC(tmp, k = log(n))
      LL[k, i] <- logLik(tmp)
      if (BIC[k, i] < tmpBIC) {
        tmpBIC <- BIC[k, i]
        result <- tmp
        choice <- c(k, i)
      }
    }
    k <- k + 1
  }

  p <- c(1, p)
  result$choice <- choice
  result$BIC <- BIC
  result$AllPartitions <- P
  result$AllLL <- LL
  result$p <- p
  class(result) <- c("pmlCluster", "pmlPart")
  result
}


#' @export
plot.pmlCluster <- function(x, which = c(1L:3L), caption =
                            list("BIC", "log-likelihood", "Partitions"), ...) {
  show <- rep(FALSE, 3)
  show[which] <- TRUE
  choice <- x$choice
  if (show[1]) {
    X <- x$AllPartitions[choice[1], , ]
    d <- dim(X)
    ind <- order(X[choice[2], ])
    im <- matrix(0, d[2], d[2])
    for (j in 1:d[1]) {
      for (i in 1:d[2]) im[i, ] <- im[i, ] + (X[j, ] == X[j, i])
    }
    image(im[ind, ind], ...)
  }

  if (show[1]) matplot(x$p, x$BIC, ylab = "BIC", xlab = "number of clusters")
  if (show[1]) matplot(x$p, x$AllLL, ylab = "log-likelihood",
                       xlab = "number of clusters")
}


#' @export
print.pmlPart <- function(x, ...) {
  nc <- attr(x$fits[[1]]$data, "nc")
  levels <- attr(x$fits[[1]]$data, "levels")
  r <- length(x$fits)
  nc <- attr(x$fits[[1]]$data, "nc")
  k <- x$fits[[1]]$k

  lbf <- x$df["Bf", 2]
  bf <- matrix(0, lbf, nc)
  if (lbf > 1) dimnames(bf) <- list(1:r, levels)
  lQ <- x$df["Q", 2]
  Q <- matrix(0, lQ, nc * (nc - 1) / 2)
  if (lQ > 1) dimnames(Q) <- list(1:r, NULL)
  type <- attr(x$fits[[1]]$data, "type")

  loli <- numeric(r)
  rate <- numeric(r)
  shape <- numeric(r)
  sizes <- numeric(r)
  inv <- numeric(r)
  for (i in 1:r) {
    loli[i] <- x$fits[[i]]$logLik
    if (i <= lbf) bf[i, ] <- x$fits[[i]]$bf
    if (i <= lQ) Q[i, ] <- x$fits[[i]]$Q
    rate[i] <- x$fits[[i]]$rate
    shape[i] <- x$fits[[i]]$shape
    inv[i] <- x$fits[[i]]$inv
    sizes[i] <- sum(attr(x$fits[[i]]$data, "weight"))
  }
  cat("\nloglikelihood:", x$logLik, "\n")
  cat("\nloglikelihood of partitions:\n ", loli, "\n")
  cat("AIC: ", AIC(x), " BIC: ", AIC(x, k = log(sum(sizes))), "\n\n")
  cat("Proportion of invariant sites:", inv, "\n")
  cat("\nRates:\n")
  cat(rate, "\n")
  if (k > 1) {
    cat("\nShape parameter:\n")
    cat(shape, "\n")
  }
  if (type == "AA") cat("Rate matrix:", x$fits[[1]]$model, "\n")
  else {
    cat("\nBase frequencies:  \n")
    print(bf)
    cat("\nRate matrix:\n")
    print(Q)
  }
}


#' @export
logLik.pmlPart <- function(object, ...) {
  res <- object$logLik
  attr(res, "df") <- sum(object$df[, 1] * object$df[, 2])
  class(res) <- "logLik"
  res
}



optNNI <- function(fit, INDEX) {
  tree <- fit$tree
  ll.0 <- fit$ll.0
  bf <- fit$bf
  eig <- fit$eig
#  k <- fit$k
  w <- fit$w
  g <- fit$g
  rootEdges <- attr(INDEX, "root")
  .dat <- NULL
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]

  data <- getCols(fit$data, tree$tip.label)
  datp <- rnodes(tree, data, w, g, eig, bf)
  tmp <- length(tree$tip.label)
  for (i in seq_along(w)) .dat[i, 1:tmp] <- new2old.phyDat(data)

  evector <- numeric(max(parent))
  evector[child] <- tree$edge.length
  m <- dim(INDEX)[1]
  loglik <- numeric(2 * m)
  edgeMatrix <- matrix(0, 2 * m, 5)
  for (i in 1:m) {
    ei <- INDEX[i, ]
    el0 <- evector[INDEX[i, ]]
    l <- length(datp[, 1])
    weight <- fit$weight
    datn <- vector("list", 4 * l)
    attr(datn, "dim") <- c(l, 4)
    datn <- .dat[, ei[1:4], drop = FALSE]
    if (!(ei[5] %in% rootEdges))
      datn[, 1] <- datp[, ei[1], drop = FALSE]
    new1 <- optim.quartet(el0[c(1, 3, 2, 4, 5)],
      eig, bf, datn[, c(1, 3, 2, 4), drop = FALSE], g,
      w, weight, ll.0, llcomp = fit$log)
    new2 <- optim.quartet(el0[c(1, 4, 3, 2, 5)],
      eig, bf, datn[, c(1, 4, 3, 2), drop = FALSE], g,
      w, weight, ll.0, llcomp = fit$log)
    loglik[(2 * i) - 1] <- new1[[2]]
    loglik[(2 * i)] <- new2[[2]]
    edgeMatrix[(2 * i) - 1, ] <- new1[[1]]
    edgeMatrix[(2 * i), ] <- new2[[1]]
  }
  list(loglik = loglik, edges = edgeMatrix)
}


optimPartNNI <- function(object, AllEdge = TRUE, ...) {
  tree <- object[[1]]$tree
  INDEX <- indexNNI(tree)
  l <- length(object)
  loglik0 <- 0
  for (i in 1:l) loglik0 <- loglik0 + logLik(object[[i]])

  l <- length(object)
  TMP <- vector("list", l)
  for (i in 1:l) {
    TMP[[i]] <- optNNI(object[[i]], INDEX)
  }
  loglik <- TMP[[1]][[1]]
  for (i in 2:l) loglik <- loglik + TMP[[i]][[1]]

  swap <- 0
  candidates <- loglik > loglik0

  while (any(candidates)) {
    ind <- which.max(loglik)
    loglik[ind] <- -Inf
    if (ind %% 2)
      swap.edge <- c(2, 3)
    else swap.edge <- c(2, 4)
    tree2 <- changeEdge(tree, INDEX[(ind + 1) %/% 2, swap.edge],
      INDEX[(ind + 1) %/% 2, ], TMP[[1]][[2]][ind, ])
    tmpll <- 0
    for (i in 1:l) {
      if (!AllEdge) tree2 <- changeEdge(object[[i]]$tree, INDEX[(ind + 1) %/% 2,
                    swap.edge], INDEX[(ind + 1) %/% 2, ], TMP[[i]][[2]][ind, ])
      tmpll <- tmpll + update(object[[i]], tree = tree2)$logLik
    }

    if (tmpll < loglik0)
      candidates[ind] <- FALSE
    if (tmpll > loglik0) {

      swap <- swap + 1
      tree <- tree2
      indi <- which(rep(colSums(apply(INDEX, 1, match,
        INDEX[(ind + 1) %/% 2, ], nomatch = 0)) > 0, each = 2))
      candidates[indi] <- FALSE
      loglik[indi] <- -Inf

      for (i in 1:l) {
        if (!AllEdge) tree2 <- changeEdge(object[[i]]$tree,
            INDEX[(ind + 1) %/% 2, swap.edge],
            INDEX[(ind + 1) %/% 2, ], TMP[[i]][[2]][ind, ])
        object[[i]] <- update(object[[i]], tree = tree2)
      }
      loglik0 <- 0
      for (i in 1:l) loglik0 <- loglik0 + logLik(object[[i]])
      cat(loglik0, "\n")
    }
  }
  if (AllEdge) object <- optimPartEdge(object)
  attr(object, "swap") <- swap
  object
}
