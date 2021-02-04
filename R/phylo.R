# allow transition probs of zero (added o)
optimQ <- function(tree, data, Q = rep(1, 6), subs = rep(1, length(Q)),
                   trace = 0, ...) {
  m <- length(Q)
  n <- max(subs)
  o <- min(subs)
  ab <- numeric(n)
  #    ab = log(Q[match(1:n, subs)])
  for (i in 1:n) ab[i] <- log(Q[which(subs == i)[1]])
  fn <- function(ab, tree, data, m, n, o, subs, ...) {
    Q <- numeric(m)
    for (i in 1:n) Q[subs == i] <- ab[i]
    if (o < 0) Q[subs < 0] <- -Inf
    pml.fit(tree, data, Q = exp(Q), ...) # Q^2, ...)
  }
  res <- optim(par = ab, fn = fn, gr = NULL, method = "L-BFGS-B", lower = -Inf,
               upper = 10, control = list(fnscale = -1, maxit = 25,
                                          trace = trace), tree = tree,
               data = data, m = m, n = n, o = o, subs = subs, ...)
  Q <- rep(1, m)
  for (i in 1:n) Q[subs == i] <- exp(res[[1]][i])
  if (o < 0) Q[subs < 0] <- 0
  res[[1]] <- Q
  res
}


CodonQ <- function(subs = .sub, syn = .syn, tstv = 1, dnds = 1) {
  Q <- numeric(1830)
  Q[subs == 1] <- 1 # transversion
  Q[subs == 2] <- tstv # transition
  Q[syn == 1] <- Q[syn == 1] * dnds
  Q[syn < 0] <- 0
  Q
}


# needs no Q
optimCodon <- function(tree, data, Q = rep(1, 1830), subs = rep(1, length(Q)),
                       syn = rep(0, length(Q)), trace = 0L, ab = c(0, 0),
                       optK = TRUE, optW = TRUE, ...) {
  m <- length(Q)
  n <- 1L
  fn <- function(ab, tree, data, m, n, subs, syn, optK, optW, ...) {
    Q <- numeric(m)
    Q[subs == 1] <- 0 # transversion
    if (optK) Q[subs == 2] <- ab[1] # transition
    else Q[subs == 2] <- 0
    if (optW) Q[syn == 1] <- Q[syn == 1] + ab[2] # ab[n+1] dnds
    Q[syn < 0] <- -Inf
    pml.fit(tree, data, Q = exp(Q), ...) # Q^2, ...)
  }
  res <- optim(par = ab, fn = fn, gr = NULL, method = "L-BFGS-B",
    lower = -Inf, upper = Inf, control = list(fnscale = -1,
      maxit = 25, trace = trace), tree = tree, data = data, m = m, n = n,
    subs = subs, syn = syn, optK = optK, optW = optW, ...)
  ab <- exp(res[[1]])
  Q[subs == 1] <- 1 # transversion
  if (optK) Q[subs == 2] <- ab[1] # transition
  else {
    Q[subs == 2] <- 1
    ab[1] <- 1
  }
  if (optW) Q[syn == 1] <- Q[syn == 1] * ab[2] # dnds
  else ab[2] <- 1
  Q[syn < 0] <- 0
  res[[5]] <- ab
  res[[1]] <- Q
  res
}


subsChoice <- function(type = c("JC", "F81", "K80", "HKY", "TrNe", "TrN",
                                "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3",
                                "TPM3u", "TIM1e", "TIM1", "TIM2e", "TIM2",
                                "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR")) {
  type <- match.arg(type)
  switch(type,
    JC = list(optQ = FALSE, optBf = FALSE,   subs = c(0, 0, 0, 0, 0, 0)),
    F81 = list(optQ = FALSE, optBf = TRUE,   subs = c(0, 0, 0, 0, 0, 0)),
    K80 = list(optQ = TRUE, optBf = FALSE,   subs = c(0, 1, 0, 0, 1, 0)),
    HKY = list(optQ = TRUE, optBf = TRUE,    subs = c(0, 1, 0, 0, 1, 0)),
    TrNe = list(optQ = TRUE, optBf = FALSE,  subs = c(0, 1, 0, 0, 2, 0)),
    TrN = list(optQ = TRUE, optBf = TRUE,    subs = c(0, 1, 0, 0, 2, 0)),
    TPM1 = list(optQ = TRUE, optBf = FALSE,  subs = c(0, 1, 2, 2, 1, 0)),
    K81 = list(optQ = TRUE, optBf = FALSE,   subs = c(0, 1, 2, 2, 1, 0)),
    TPM1u = list(optQ = TRUE, optBf = TRUE,  subs = c(0, 1, 2, 2, 1, 0)),
    TPM2 = list(optQ = TRUE, optBf = FALSE,  subs = c(1, 2, 1, 0, 2, 0)),
    TPM2u = list(optQ = TRUE, optBf = TRUE,  subs = c(1, 2, 1, 0, 2, 0)),
    TPM3 = list(optQ = TRUE, optBf = FALSE,  subs = c(1, 2, 0, 1, 2, 0)),
    TPM3u = list(optQ = TRUE, optBf = TRUE,  subs = c(1, 2, 0, 1, 2, 0)),
    TIM1e = list(optQ = TRUE, optBf = FALSE, subs = c(0, 1, 2, 2, 3, 0)),
    TIM1 = list(optQ = TRUE, optBf = TRUE,   subs = c(0, 1, 2, 2, 3, 0)),
    TIM2e = list(optQ = TRUE, optBf = FALSE, subs = c(1, 2, 1, 0, 3, 0)),
    TIM2 = list(optQ = TRUE, optBf = TRUE,   subs = c(1, 2, 1, 0, 3, 0)),
    TIM3e = list(optQ = TRUE, optBf = FALSE, subs = c(1, 2, 0, 1, 3, 0)),
    TIM3 = list(optQ = TRUE, optBf = TRUE,   subs = c(1, 2, 0, 1, 3, 0)),
    TVMe = list(optQ = TRUE, optBf = FALSE,  subs = c(1, 2, 3, 4, 2, 0)),
    TVM = list(optQ = TRUE, optBf = TRUE,    subs = c(1, 2, 3, 4, 2, 0)),
    SYM = list(optQ = TRUE, optBf = FALSE,   subs = c(1, 2, 3, 4, 5, 0)),
    GTR = list(optQ = TRUE, optBf = TRUE,    subs = c(1, 2, 3, 4, 5, 0))
  )
}



optimGamma <- function(tree, data, shape = 1, k = 4, ...) {
  fn <- function(shape, tree, data, k, ...) pml.fit(tree, data, shape = shape,
      k = k, ...)
  res <- optimize(f = fn, interval = c(0.1, 100), lower = 0.1, upper = 1000,
    maximum = TRUE,  tol = .01, tree = tree, data = data, k = k, ...)
  res
}


optimInv <- function(tree, data, inv = 0.01, INV = NULL, ll.0 = NULL, ...) {
  fn <- function(inv, tree, data, ...) pml.fit(tree, data, inv = inv, INV = INV,
      ll.0 = NULL, ...)
  res <- optimize(f = fn, interval = c(0, 1), lower = 0, upper = 1,
    maximum = TRUE, tol = .0001, tree = tree, data = data, ...)
  res
}


# changed to c(-10,10) from c(-5,5)
optimRate <- function(tree, data, rate = 1, ...) {
  fn <- function(rate, tree, data, ...)
    pml.fit(tree, data, rate = exp(rate), ...)
  res <- optimize(f = fn, interval = c(-10, 10), tree = tree, data = data,
                  ..., maximum = TRUE)
  res[[1]] <- exp(res[[1]])
  res
}


optimBf <- function(tree, data, bf = c(.25, .25, .25, .25), trace = 0, ...) {
  l <- length(bf)
  nenner <- 1 / bf[l]
  lbf <- log(bf * nenner)
  lbf <- lbf[-l]
  fn <- function(lbf, tree, data, ...) {
    bf <- exp(c(lbf, 0))
    bf <- bf / sum(bf)
    pml.fit(tree, data, bf = bf, ...)
  }
  res <- optim(par = lbf, fn = fn, gr = NULL, method = "Nelder-Mead",
               control = list(fnscale = -1, maxit = 500, trace = trace),
               tree = tree, data = data, ...)
  bf <- exp(c(res[[1]], 0))
  bf <- bf / sum(bf)
  result <- list(bf = bf, loglik = res[[2]])
  result
}


# MLF3x4 model
optimF3x4 <- function(tree, data, bf_codon = matrix(.25, 4, 3), trace = 0, ...) {
  l <- nrow(bf_codon)
  nenner <- 1 / bf_codon[l, ]
  lbf <- log(bf_codon * rep(nenner, each = 4))
  lbf <- lbf[-l, ]
  fn <- function(lbf, tree, data, ...) {
    dim(lbf) <- c(3, 3)
    bf_codon <- rbind(exp(lbf), c(1, 1, 1))
    bf_codon <- bf_codon / rep(colSums(bf_codon), each = 4)
    bf <- F3x4_freq(bf_codon)
    pml.fit(tree, data, bf = bf, ...)
  }
  res <- optim(par = lbf, fn = fn, gr = NULL, method = "Nelder-Mead", control =
    list(fnscale = -1, maxit = 500, trace = trace), tree = tree,
  data = data, ...)
  bf_codon <- rbind(exp(res[[1]]), c(1, 1, 1))
  bf_codon <- bf_codon / rep(colSums(bf_codon), each = 4)
  bf <- F3x4_freq(bf_codon)
  result <- list(bf = bf, loglik = res[[2]], bf_codon = bf_codon)
  result
}



optimW <- function(fit, ...) {
  w <- fit$w
  g <- fit$g
  siteLik <- fit$siteLik
  k <- length(w)
  l <- dim(siteLik[[1]])[1]
  x <- matrix(0, l, k)
  for (i in 1:k) x[, i] <- rowSums(siteLik[[i]])
  weight <- fit$weight
  nenner <- 1 / w[k]
  eta <- log(w * nenner)
  eta <- eta[-k]
  fn <- function(eta, x, g, weight) {
    eta <- c(eta, 0)
    p <- exp(eta) / sum(exp(eta))
    res <- x %*% p
    res <- sum(weight * log(res))  * (1 + abs(sum(p * g) - 1))
    res
  }
  res <- optim(eta, fn = fn, method = "Nelder-Mead",
               control = list(fnscale = -1, reltol = 1e-12), gr = NULL, x = x,
               g = g, weight = weight)
  p <- exp(c(res$par, 0))
  p <- p / sum(p)
  result <- list(par = p, value = res$value)
  result
}


# predict.pml <- function(object, newdata,...) sum(object$siteLik * newdata)


#' @export
logLik.pml <- function(object, ...) {
  res <- object$logLik
  attr(res, "df") <- object$df
  class(res) <- "logLik"
  res
}


#' @export
AICc <- function(object, ...)
  UseMethod("AICc")


#' @export
AICc.pml <- function(object, ...) {
  n <- sum(object$weight)
  k <- object$df
  if (k >= (n - 1)) return(NaN)
  res <- AIC(object)
  res +   (2 * k * (k + 1)) / (n - k - 1)
}


#' @export
BIC.pml <- function(object, ...) {
  res <- AIC(object, k = log(sum(object$weight)))
  res
}


#' @export
anova.pml <- function(object, ...) {
  X <- c(list(object), list(...))
  fun <- function(x) {
    tmp <- logLik(x)
    c(tmp[1], attr(tmp, "df"))
  }
  DF <- t(sapply(X, fun))
  dev <- c(NA, 2 * diff(DF[, 1]))
  ddf <- c(NA, diff(DF[, 2]))
  table <- data.frame(DF, ddf, dev, pchisq(dev, ddf, lower.tail = FALSE))
  dimnames(table) <- list(seq_along(X), c("Log lik.", "Df",
    "Df change", "Diff log lik.", "Pr(>|Chi|)"))
  structure(table, heading = "Likelihood Ratio Test Table",
    class = c("anova", "data.frame"))
}

# vcov.pml <- function(object, obs=FALSE,...){
#    if(obs) FI = score4(object)[[2]]
#    else FI = score(object,FALSE)[[2]]
#    l = dim(FI)[1]
#    res = try(solve(FI))
#    if(class(res) == "try-error"){
#        cat("Covariance is ill-conditioned !! \n")
#        res = solve(FI + diag(l)* 1e-8)
#        }
#    res
# }

#' @export
vcov.pml <- function(object, ...) {
  FI <- score(object, FALSE)[[2]]
  l <- dim(FI)[1]
  res <- try(solve(FI))
  if (inherits(res, "try-error")) {
    cat("Covariance is ill-conditioned !! \n")
    res <- solve(FI + diag(l) * 1e-8)
  }
  res
}


getd2P <- function(el, eig = edQt(), g = 1.0) {
  n <- length(eig$values)
  res <- .Call("getd2PM", eig, as.integer(n), as.double(el), as.double(g))
  attr(res, "dim") <- c(length(g), length(el))
  res
}


getdP <- function(el, eig = edQt(), g = 1.0) {
  n <- length(eig$values)
  res <- .Call("getdPM", eig, as.integer(n), as.double(el), as.double(g))
  attr(res, "dim") <- c(length(g), length(el))
  res
}


# version without transformation (used for vcov)
getdP2 <- function(el, eig = edQt(), g = 1.0) {
  n <- length(eig$values)
  res <- .Call("getdPM2", eig, as.integer(n), as.double(el), as.double(g))
  attr(res, "dim") <- c(length(g), length(el))
  res
}


# version without transformation
getd2P2 <- function(el, eig = edQt(), g = 1.0) {
  n <- length(eig$values)
  res <- .Call("getd2PM2", eig, as.integer(n), as.double(el), as.double(g))
  attr(res, "dim") <- c(length(g), length(el))
  res
}


getP <- function(el, eig = edQt(), g = 1.0) {
  n <- length(eig$values)
  res <- .Call("getPM", eig, as.integer(n), as.double(el), as.double(g))
  attr(res, "dim") <- c(length(g), length(el))
  res
}


#' @rdname pml.fit
#' @export
lli <- function(data, tree = NULL, ...)  {
  contrast <- attr(data, "contrast")
  nr <- attr(data, "nr")
  nc <- attr(data, "nc")
  nco <- as.integer(dim(contrast)[1])
  if (!is.null(tree)) data <- subset(data, tree$tip.label)
  .Call("invSites", data, as.integer(nr), as.integer(nc), contrast,
    as.integer(nco))
}


#' @rdname pml.fit
#' @export
edQt <- function(Q = c(1, 1, 1, 1, 1, 1), bf = c(0.25, 0.25, 0.25, 0.25)) {
  l <- length(bf)
  res <- matrix(0, l, l)
  res[lower.tri(res)] <- Q
  res <- res + t(res)
  res <- res * bf
  res2 <- res * rep(bf, each = l)
  diag(res) <- -colSums(res)
  res <- res / sum(res2)
  e <- eigen(res, FALSE)
  e$inv <- solve.default(e$vec)
  e
}


# some functions to get codon models to work
getScaler <- function(Q = c(1, 1, 1, 1, 1, 1), bf = c(0.25, 0.25, 0.25, 0.25)) {
  l <- length(bf)
  res <- matrix(0, l, l)
  res[lower.tri(res)] <- Q
  res <- res + t(res)
  res <- res * bf
  res2 <- res * rep(bf, each = l)
  diag(res) <- -colSums(res)
  sum(res2)
}

# eigen value decomposition without scaling
edQt2 <- function(Q = c(1, 1, 1, 1, 1, 1), bf = c(0.25, 0.25, 0.25, 0.25),
                  scale = 1) {
  l <- length(bf)
  res <- matrix(0, l, l)
  res[lower.tri(res)] <- Q
  res <- res + t(res)
  res <- res * bf
  res2 <- res * rep(bf, each = l)
  diag(res) <- -colSums(res)
  res <- res / scale
  e <- eigen(res, FALSE)
  e$inv <- solve.default(e$vec)
  e
}


df_freq_codon <- function(bf) {
  switch(bf,
         equal = 0,
         empirical = 60,
         F61 = 60,
         F3x4 = 9,
         F1x4 = 3)
}


#' @rdname pml.fit
#' @export
pml.free <- function() {
  .Call("ll_free2")
  #    rm(.INV, .iind, envir = parent.frame())
}


#' @rdname pml.fit
#' @export
pml.init <- function(data, k = 1L) {
  nTips <- length(data)
  nr <- attr(data, "nr")
  nc <- attr(data, "nc")
  .Call("ll_init2", as.integer(nr), as.integer(nTips), as.integer(nc),
    as.integer(k))
}


fn.quartet <- function(old.el, eig, bf, dat,  g = 1, w = 1, weight, ll.0) {
  l <- length(dat[, 1])
  ll <- ll.0
  res <- vector("list", 2 * l)
  tmp1 <- NULL
  tmp2 <- NULL
  attr(res, "dim") <- c(l, 2)
  for (j in 1:l) {
    P <- getP(old.el, eig, g[j])
    tmp1 <- (dat[[j, 1]] %*% P[[1]]) * (dat[[j, 2]] %*% P[[2]])
    tmp2 <- (dat[[j, 3]] %*% P[[3]]) * (dat[[j, 4]] %*% P[[4]])
    res[[j, 1]] <- tmp1 * (tmp2 %*% P[[5]])
    res[[j, 2]] <- tmp2
    ll <- ll +  res[[j, 1]] %*% (w[j] * bf)
  }
  l0 <- sum(weight * log(ll))
  list(ll = l0, res = res)
}


rnodes <- function(tree, data, w, g, eig, bf) {
  if (is.null(attr(tree, "order")) || attr(tree, "order") ==
    "cladewise")
    tree <- reorder(tree, "postorder")
  data <- getCols(data, tree$tip.label)
  q <- length(tree$tip.label)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  m <- length(edge) + 1  # max(edge)
  l <- length(w)
  dat <- vector(mode = "list", length = m * l)
  dim(dat) <- c(l, m)
  tmp <- length(data)
  #    for(i in seq_along(w))dat[i,1:tmp]=new2old.phyDat(data) #
  #    dat[1,1:tmp] <- data  vielleicht gebraucht
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
  for (i in 1:l) dat[i, (q + 1):m] <- .Call("LogLik2", data, P[i, ], nr, nc,
                                            node, edge, nTips, mNodes, contrast,
                                            nco)
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  nTips <- min(parent) - 1
  datp <- vector("list", m)
  dat2 <- vector("list", m * l)
  dim(dat2) <- c(l, m)
  for (i in 1:l) {
    datp[(nTips + 1)] <- dat[i, (nTips + 1)]
    for (j in (m - 1):1) {
      if (child[j] > nTips) {
        tmp2 <- (datp[[parent[j]]] / (dat[[i, child[j]]] %*% P[[i, j]]))
        datp[[child[j]]] <- (tmp2 %*% P[[i, j]]) * dat[[i, child[j]]]
        dat2[[i, child[j]]] <- tmp2
      }
    }
  }
  assign(".dat", dat, envir = parent.frame(n = 1))
  dat2
}


score <- function(fit, transform = TRUE) {
  tree <- fit$tree
  child <- tree$edge[, 2]
  l <- length(child)
  sc <- numeric(l)
  weight <- as.numeric(fit$weight)
  f <- drop(exp(fit$siteLik))
  dl <- dl(fit, transform)
  dl <- dl / f
  sc <- colSums(weight * dl)
  F <- crossprod(dl * weight, dl)
  names(sc) <- child
  dimnames(F) <- list(child, child)
  result <- list(sc = sc, F = F)
  result
}


# wird noch in partition models verwendet
optim.quartet <- function(old.el, eig, bf, dat, g = 1, w = 1, weight,
                          ll.0 = weight * 0, control = list(eps = 1e-08,
                                        maxit = 5, trace = 0), llcomp = -Inf) {
  eps <- 1
  iter <- 0
  evi <- (t(eig[[3]]) * bf)
  while (eps > control$eps && iter < control$maxit) {
    tmp <- fn.quartet(old.el = old.el, eig = eig, bf = bf, dat = dat,
      g = g, w = w, weight = weight, ll.0 = ll.0)
    old.ll <- tmp$ll
    el1 <- fs(old.el[1], eig, tmp$res[, 1], dat[, 1], weight,
      g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA = TRUE, getB = FALSE)
    el2 <- fs(old.el[2], eig, el1[[2]], dat[, 2], weight,
      g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA = TRUE, getB = FALSE)
    el5 <- fs(old.el[5], eig, el2[[2]], tmp$res[, 2], weight,
      g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA = FALSE, getB = TRUE)
    el3 <- fs(old.el[3], eig, el5[[3]], dat[, 3], weight,
      g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA = TRUE, getB = FALSE)
    el4 <- fs(old.el[4], eig, el3[[2]], dat[, 4], weight,
      g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA = FALSE, getB = FALSE)
    old.el[1] <- el1[[1]]
    old.el[2] <- el2[[1]]
    old.el[3] <- el3[[1]]
    old.el[4] <- el4[[1]]
    old.el[5] <- el5[[1]]
    iter <- iter + 1
    ll <- el4[[4]]
    eps <- (old.ll - ll) / ll
    if (ll < llcomp) return(list(old.el, ll))
    old.ll <- ll
  }
  list(old.el, ll)
}


#' @export
plot.pml <- function(x, ...) plot.phylo(x$tree, ...)


phangornParseFormula <- function(model) {

  parseSide <- function(model) {
    model.vars <- list()
    while (length(model) == 3 && model[[1]] == as.name("+")) {
      model.vars <- c(model.vars, model[[3]])
      model <- model[[2]]
    }
    unlist(rev(c(model.vars, model)))

  }

  if (!inherits(model, "formula"))
    stop("model must be a formula object")
  l <- length(model)
  varsLHS <- NULL
  if (l == 3) {
    modelLHS <- model[[2]]
    modelRHS <- model[[3]]
    varsRHS <- parseSide(modelRHS)
    varsRHS <- unlist(lapply(varsRHS, as.character))
    varsLHS <- parseSide(modelLHS)
    varsLHS <- unlist(lapply(varsLHS, as.character))
  }
  if (l == 2) {
    modelRHS <- model[[2]]
    varsRHS <- parseSide(modelRHS)
    varsRHS <- unlist(lapply(varsRHS, as.character))
  }
  list(left = varsLHS, right = varsRHS)
}


#' @rdname pml
#' @export
pml.control <- function(epsilon = 1e-08, maxit = 10, trace = 1) {
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  list(epsilon = epsilon, maxit = maxit, trace = trace)
}


likelihoodRatchet <- function(obj, maxit = 100, k = 10,
                              control = pml.control(epsilon = 1e-08, maxit = 10,
                                                    trace = 1L)) {
  tree <- obj$tree
  nTips <- length(tree$tip.label)
  trace <- control$trace
  control$trace <- trace - 1L
  kmax <- 1
  for (i in 1:maxit) {
    tree <- rNNI(obj$tree, moves = nTips / 3, n = 1)
    # tree <- rSPR(tree, moves=10, k=3, n=1)
    obj2 <- update(obj, tree = tree)
    obj2 <- optim.pml(obj2, TRUE, control = control)
    if (logLik(obj2) > logLik(obj)) {
      obj <- obj2
      kmax <- 1
    }
    else kmax <- kmax + 1
    if (trace > 0) print(paste("Iteration ", i, ", best pscore so far:",
        logLik(obj)))
    if (kmax == k) break()
  }
  obj
}


fs <- function(old.el, eig, parent.dat, child.dat, weight, g = g,
               w = w, bf = bf, ll.0 = ll.0, evi, getA = TRUE, getB = TRUE) {
  if (old.el < 1e-8) old.el <- 1e-8
  lg <- length(parent.dat)
  P <- getP(old.el, eig, g)
  nr <- as.integer(length(weight))
  nc <- as.integer(length(bf))
  eve <- eig[[2]]
  dad <- .Call("getDAD", parent.dat, child.dat, P, nr, nc)
  X <- .Call("getPrep", dad, child.dat, eig[[2]], evi, nr, nc)
  .Call("FS4", eig, as.integer(length(bf)), as.double(old.el),
    as.double(w), as.double(g), X, child.dat, dad, as.integer(length(w)),
    as.integer(length(weight)), as.double(weight),
    as.double(ll.0), as.integer(getA), as.integer(getB))
}


# rate not used internally
optimEdge <- function(tree, data, eig = eig, w = w, g = g, bf = bf, rate = rate,
                      ll.0 = ll.0, control = pml.control(epsilon = 1e-08,
                        maxit = 10, trace = 0), ...) {
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorder(tree, "postorder")
  nTips <- length(tree$tip.label)
  el <- tree$edge.length
  tree$edge.length[el < 1e-08] <- 1e-08
  oldtree <- tree
  k <- length(w)
  data <- subset(data, tree$tip.label)
  loglik <- pml.fit4(tree, data, bf = bf, g = g, w = w, eig = eig, ll.0 = ll.0,
                     k = k)
  start.ll <- old.ll <- loglik
  contrast <- attr(data, "contrast")
  contrast2 <- contrast %*% eig[[2]]
  evi <- (t(eig[[3]]) * bf)
  weight <- attr(data, "weight")
  eps <- 1
  iter <- 0

  treeP <- tree
  tree <- reorder(tree)

  child <- tree$edge[, 2]
  parent <- tree$edge[, 1]
  m <- max(tree$edge)
  pvec <- integer(m)
  pvec[child] <- parent

  EL <- numeric(m)
  EL[child] <- tree$edge.length

  n <- length(tree$edge.length)

  nr <- as.integer(length(weight))
  nc <- as.integer(length(bf))
  nco <- as.integer(nrow(contrast))
  eve <- eig[[2]]
  lg <- k
  ScaleEPS <- 1.0 / 4294967296.0
  anc <- Ancestors(tree, 1:m, "parent")
  anc0 <- as.integer(c(0L, anc))

  while (eps > control$eps && iter < control$maxit) {
    EL <- .Call("optE", as.integer(parent), as.integer(child),
      as.integer(anc0), eig, evi, EL, w, g, as.integer(nr),
      as.integer(nc), as.integer(nTips), as.double(contrast),
      as.double(contrast2), nco, data, as.double(weight),
      as.double(ll.0))
    iter <- iter + 1
    treeP$edge.length <- EL[treeP$edge[, 2]]
    newll <- pml.fit4(treeP, data, bf = bf, g = g, w = w, eig = eig,
                      ll.0 = ll.0, k = k)
    eps <- (old.ll - newll) / newll
    if (eps < 0) return(list(oldtree, old.ll))
    oldtree <- treeP
    if (control$trace > 1) cat(old.ll, " -> ", newll, "\n")
    old.ll <- newll
  }
  if (control$trace > 0)
    cat("optimize edge weights: ", start.ll, "-->", newll, "\n")
#    cat(start.ll, " -> ", newll, "\n")
  list(tree = treeP, logLik = newll, c(eps, iter))
}

# only use of PML3
pml.move <- function(EDGE, el, data, g, w, eig, k, nTips, bf) {
  node <- EDGE[, 1]
  edge <- EDGE[, 2]
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  node <- as.integer(node - nTips - 1L)
  edge <- as.integer(edge - 1L)
  contrast <- attr(data, "contrast")
  nco <- as.integer(dim(contrast)[1])
  tmp <- .Call("PML3", dlist = data, as.double(el), as.double(g), nr, nc, k,
    eig, as.double(bf), node, edge, nTips, nco, contrast,
    N = as.integer(length(edge)))
  # as.double(w),
  return(NULL)
}


bip <- function(x) {
  x <- reorder(x, "postorder")
  nTips <- as.integer(length(x$tip.label))
  .Call("_phangorn_bipCPP", PACKAGE = "phangorn", x$edge, nTips)
}


bipart <- function(x) {
  x <- reorder(x, "postorder")
  nTips <- as.integer(length(x$tip.label))
  .Call("_phangorn_bipartCPP", PACKAGE = "phangorn", x$edge, nTips)
}


bipartition <- function(tree) {
  if (is.rooted(tree)) tree <- unroot(tree)
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorder(tree, "postorder")
  bp <- bip(tree)
  nTips <- length(tree$tip.label)
  l <- length(bp)
  res <- matrix(0L, l, nTips)
  for (i in 1:l) res[i, bp[[i]]] <- 1L
  res <- res[tree$edge[, 2], , drop = FALSE]
  colnames(res) <- tree$tip.label
  rownames(res) <- tree$edge[, 2]
  res[res[, 1] == 1, ] <- 1L - res[res[, 1] == 1, ]
  res
}


readAArate <- function(file) {
  tmp <- read.table(system.file(file.path("extdata", file)), col.names = 1:20,
    fill = TRUE)
  Q <- tmp[1:19, 1:19]
  names <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m",
    "f", "p", "s", "t", "w",  "y", "v")
  Q <- as.numeric(Q[lower.tri(Q, TRUE)])
  bf <- as.numeric(as.character(unlist(tmp[20, ])))
  names(bf) <- names
  list(Q = Q, bf = bf)
}


# .LG <- readAArate("lg.dat")
# .WAG <- readAArate("wag.dat")
# .Dayhoff <- readAArate("dayhoff-dcmut.dat")
# .JTT <- readAArate("jtt-dcmut.dat")
# .cpREV <- readAArate("cpREV.dat")
# .mtmam <- readAArate("mtmam.dat")
# .mtArt <- readAArate("mtArt.dat")
# save(.LG,.WAG,.Dayhoff,.JTT,.cpREV,.mtmam,.mtArt, file = "sysdata2.rda")


getModelAA <- function(model, bf = TRUE, Q = TRUE) {
  model <- match.arg(eval(model), .aamodels)
  tmp <- get(paste(".", model, sep = ""), environment(pml))
  if (Q) assign("Q", tmp$Q, envir = parent.frame())
  if (bf) assign("bf", tmp$bf, envir = parent.frame())
}


#' @export
print.pml <- function(x, ...) {
  cat("\n loglikelihood:", x$logLik, "\n")
  w <- x$weight
  w <- w[w > 0]
  type <- attr(x$data, "type")
  levels <- attr(x$data, "levels")
  nc <- attr(x$data, "nc")
  ll0 <- sum(w * log(w / sum(w)))
  cat("\nunconstrained loglikelihood:", ll0, "\n")
  if (x$inv > 0) cat("Proportion of invariant sites:", x$inv, "\n")
  if (x$k > 1) {
    cat("Discrete gamma model\n")
    cat("Number of rate categories:", x$k, "\n")
    cat("Shape parameter:", x$shape, "\n")
  }
  if (type == "AA") cat("Rate matrix:", x$model, "\n")
  if (type == "DNA") {
    cat("\nRate matrix:\n")
    QM <- matrix(0, nc, nc, dimnames = list(levels, levels))
    QM[lower.tri(QM)] <- x$Q
    QM <- QM + t(QM)
    print(QM)
    cat("\nBase frequencies:  \n")
    bf <- x$bf
    names(bf) <- levels
    cat(bf, "\n")
  }
  if (type == "CODON") {
    cat("dn/ds:", x$dnds, "\n")
    cat("ts/tv:", x$tstv, "\n")
    cat("Freq:", x$frequencies, "\n")
  }
  if (type == "USER" & length(x$bf) < 11) {
    cat("\nRate matrix:\n")
    QM <- matrix(0, nc, nc, dimnames = list(levels, levels))
    QM[lower.tri(QM)] <- x$Q
    QM <- QM + t(QM)
    print(QM)
    cat("\nBase frequencies:  \n")
    bf <- x$bf
    names(bf) <- levels
    cat(bf, "\n")
  }
}


optEdgeMulti <- function(object, control = pml.control(epsilon = 1e-8,
                           maxit = 10, trace = 1), ...) {
  tree <- object$tree
  theta <- object$tree$edge.length
  weight <- attr(object$data, "weight")
  ll0 <- object$logLik
  eps <- 1
  iter <- 0
  iter2 <- 0
  scale <- 1
  # l <- length(theta)
  while (abs(eps) > control$eps && iter < control$maxit) {
    dl <- score(object)
    thetaNew <- log(theta) + scale * solve(dl[[2]], dl[[1]]) # + diag(l) * 1e-10
    newtheta <- exp(thetaNew)
    tree$edge.length <- as.numeric(newtheta)
    object <- update(object, tree = tree)
    ll1 <- object$logLik
    eps <- (ll0 - ll1) / ll1
    if (eps < 0) {
      newtheta <- theta
      scale <- scale / 2
      tree$edge.length <- as.numeric(theta)
      ll1 <- ll0
      iter2 <- iter2 + 1
    }
    else {
      scale <- 1
      iter2 <- 0
    }
    theta <- newtheta
    if (iter2 == 0 && control$trace > 0) cat("loglik: ", ll1, "\n")
    ll0 <- ll1
    if (iter2 == 10) iter2 <- 0
    if (iter2 == 0) iter <- iter + 1
  }
  object <- update(object, tree = tree)
  object
}


# add data for internal use parent.frame(n) for higher nestings
update.pmlNew <- function(object, ..., evaluate = TRUE) {
  call <- object$call
  if (is.null(call))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate)
    eval(call, object, parent.frame())
  else call
}


#' @export
update.pml <- function(object, ...) {
  extras <- match.call(expand.dots = FALSE)$...
  pmla <- c("tree", "data", "bf", "Q", "inv", "k", "shape",
    "rate", "model", "wMix", "llMix", "dnds", "tstv", "scaleQ","...")
  names(extras) <- pmla[pmatch(names(extras), pmla[-length(pmla)])]
  call <- object$call
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  existing <- match(pmla, names(extras))
  updateEig <- FALSE
  updateRates <- FALSE
  Mkv <- object$Mkv
  site.rate <- object$site.rate
  type <- attr(object$data, "type")
  if(type=="CODON"){
    bf_choice <- object$frequencies
  }
  if (is.na(existing[1])) tree <- object$tree
  else tree <- eval(extras[[existing[1]]], parent.frame())
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorder(tree, "postorder")
  if (is.na(existing[2])) {
    data <- object$data
    INV <- object$INV
  }
  else {
    data <- eval(extras[[existing[2]]], parent.frame())
    ll.0 <- numeric(attr(data, "nr"))
    INV <- Matrix(lli(data, tree), sparse = TRUE)
  }
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))

  if (is.na(existing[3])) bf <- object$bf
  else {
    bf <- eval(extras[[existing[3]]], parent.frame())
    # check for "character"
    if (is.character(bf)) {
      bf_choice <- match.arg(bf, c("equal", "empirical", "F1x4", "F3x4", "F61"))
      # if(bf_choice %in% c("F1x4", "F3x4", "F61") & type != "CODON")
      if (bf_choice == "F3x4" & type != "CODON")
        stop("F3x4 not available for this data type")
      if (bf_choice == "F1x4" & type != "CODON")
        stop("F1x4 not available for this data type")
      if (bf_choice == "F61" & type != "CODON")
        stop("F61 not available for this data type")
      bf <- switch(bf_choice,
                   equal = rep(1 / nc, nc),
                   empirical = baseFreq(data),
                   F61 = baseFreq(data),
                   F3x4 = F3x4(data),
                   F1x4 = F1x4(data))
      freq_df <- df_freq_codon(bf_choice)
      names(bf) <- NULL
    }
    updateEig <- TRUE
  }
  if (is.na(existing[4])) Q <- object$Q
  else {
    Q <- eval(extras[[existing[4]]], parent.frame())
    updateEig <- TRUE
  }
  type <- attr(object$data, "type")
  model <- NULL
  if (type == "AA") {
    if (!is.na(existing[9])) {
      model <- match.arg(eval(extras[[existing[9]]], parent.frame()),
        .aamodels)
      getModelAA(model, bf = is.na(existing[3]), Q = is.na(existing[4]))
      updateEig <- TRUE
    }
    #        else model <- object$model
  }
  scaleQ <- FALSE
  if (type == "CODON") {
    if (is.na(existing[12])) dnds <- object$dnds
    else {
      dnds <- eval(extras[[existing[12]]], parent.frame())
      updateEig <- TRUE
    }
    if (is.na(existing[13])) tstv <- object$tstv
    else {
      tstv <- eval(extras[[existing[13]]], parent.frame())
      updateEig <- TRUE
    }
    if (!is.na(existing[14])) {
      scaleQ <- eval(extras[[existing[14]]], parent.frame())
      updateEig <- TRUE
    }
    if (updateEig){
      .syn <- synonymous_subs(code=attr(data, "code"))
      .sub <- tstv_subs(code=attr(data, "code"))
      Q <- CodonQ(subs = .sub, syn = .syn, tstv = tstv,
                               dnds = dnds)
    }
  }
  if (is.na(existing[5])) inv <- object$inv
  else {
    inv <- eval(extras[[existing[5]]], parent.frame())
    updateRates <- TRUE
  }
  if (is.na(existing[6])) k <- object$k
  else {
    k <- eval(extras[[existing[6]]], parent.frame())
    updateRates <- TRUE
  }
  if (is.na(existing[7])) shape <- object$shape
  else {
    shape <- eval(extras[[existing[7]]], parent.frame())
    updateRates <- TRUE
  }
  rate <- ifelse(is.na(existing[8]), object$rate,
    eval(extras[[existing[8]]], parent.frame()))
  wMix <- ifelse(is.na(existing[10]), object$wMix,
    eval(extras[[existing[10]]], parent.frame()))
  if (is.na(existing[11])) llMix <- object$llMix
  else llMix <- eval(extras[[existing[11]]], parent.frame())
  levels <- attr(data, "levels")
  weight <- attr(data, "weight")
  if (updateEig){
    if(scaleQ) eig <- edQt2(Q = Q, bf = bf, scale = scaleQ)
    else eig <- edQt(Q = Q, bf = bf)
  }
  else {
    eig <- object$eig
    model <- object$model
  }

  rw <- rates_n_weights(shape, k, site.rate)
  g <- rw[, 1]
  w <- rw[, 2]

  if (inv > 0){
    w <- (1 - inv) * w
    g <- g / (1 - inv)
  }
  if (wMix > 0) w <- (1 - wMix) * w
  g <- g * rate

  ll.0 <- as.matrix(INV %*% (bf * inv))
  if (wMix > 0) ll.0 <- ll.0 + llMix

  m <- 1
  ### play save
  kmax <- k
  if (any(g < .gEps)) {
    for (i in seq_along(g)) {
      if (g[i] < .gEps) {
        inv <- inv + w[i]
      }
    }
    w <- w[g > .gEps]
    g <- g[g > .gEps]
    k <- length(w)
  }
  ####

  resll <- matrix(0, nr, k)
  nTips <- as.integer(length(tree$tip.label))

  data <- subset(data, tree$tip.label)

  on.exit(.Call("ll_free2"))
  .Call("ll_init2", nr, nTips, nc, as.integer(k))
  tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q,
    levels = attr(data, "levels"), inv = inv, rate = rate, g = g, w = w,
    eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix, wMix = wMix,
    site = TRUE)

  df <- ifelse(is.ultrametric(tree), tree$Nnode, length(tree$edge.length))
  df <- switch(type,
    DNA = df + (k > 1) + (inv > 0) + length(unique(bf)) - 1 +
      length(unique(Q)) - 1,
    AA = df + (k > 1) + (inv > 0),
    CODON = df + (k > 1) + (inv > 0) + length(unique(bf)) - 1 + (dnds != 1)
      + (tstv != 1),
    USER = df + (k > 1) + (inv > 0) + length(unique(bf)) - 1 +
      length(unique(Q)) - 1)

  result <- list(logLik = tmp$loglik, inv = inv, k = kmax, shape = shape,
    Q = Q, bf = bf, rate = rate, siteLik = tmp$siteLik,
    weight = weight, g = g, w = w, eig = eig, data = data,
    model = model, INV = INV, ll.0 = ll.0, tree = tree,
    lv = tmp$resll, call = call, df = df, wMix = wMix,
    llMix = llMix, Mkv=Mkv, site.rate=site.rate)
  if (type == "CODON") {
    result$dnds <- dnds
    result$tstv <- tstv
    result$frequencies <- bf_choice
  }
  class(result) <- "pml"
  result
}


### this is the version we want to optimise
pml.fit4 <- function(tree, data, bf = rep(1 / length(levels), length(levels)),
                     shape = 1, k = 1, Q = rep(1, length(levels) *
                                                 (length(levels) - 1) / 2),
                     levels = attr(data, "levels"), inv = 0, rate = 1, g = NULL,
                     w = NULL, eig = NULL, INV = NULL, ll.0 = NULL,
                     llMix = NULL, wMix = 0, ..., site = FALSE,
                     site.rate = "gamma") {
  weight <- as.double(attr(data, "weight"))
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  nTips <- as.integer(length(tree$tip.label))
  k <- as.integer(k)
  m <- 1
  if (is.null(eig))
    eig <- edQt(bf = bf, Q = Q)

  if(is.null(g) | is.null(w)){
    rw <- rates_n_weights(shape, k, site.rate)
    g <- rw[, 1]
    w <- rw[, 2]
    if (inv > 0){
      w <- (1 - inv) * w
      g <- g / (1 - inv)
    }
    if (wMix > 0) w <- (1 - wMix) * w
    g <- g * rate
  }

  if (any(g < .gEps)) {
    for (i in seq_along(g)) {
      if (g[i] < .gEps) {
        inv <- inv + w[i]
      }
    }
    w <- w[g > .gEps]
    g <- g[g > .gEps]
    #        kmax <- k
    k <- length(w)
  }
  #    .iind <- get(".iind", parent.frame())
  #    .INV <- get(".INV", parent.frame())
  # if(is.null(ll.0))

  if (is.null(ll.0)) {
    ll.0 <- numeric(attr(data, "nr"))
  }
  if (inv > 0)
    ll.0 <- as.matrix(INV %*% (bf * inv))
  #    if(inv>0)
  #        ll.0 <- as.matrix(.INV %*% (bf * inv))

  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  #    root <- as.integer(node[length(node)])
  el <- as.double(tree$edge.length)
  node <- as.integer(node - nTips - 1L) #    min(node))
  edge <- as.integer(edge - 1L)

  contrast <- attr(data, "contrast")
  nco <- as.integer(dim(contrast)[1])

  siteLik <- .Call("PML4", dlist = data, el, as.double(w), as.double(g), nr, nc,
    k, eig, as.double(bf), node, edge, nTips, nco, contrast,
    N = as.integer(length(edge)))
  ind <- which(ll.0 > 0)
  if (!is.null(ll.0)) siteLik[ind] <- log(exp(siteLik[ind]) + ll.0[ind])
  if (wMix > 0) siteLik <- log(exp(siteLik) * (1 - wMix) + llMix)
  loglik <- sum(weight * siteLik)
  if (!site) return(loglik)
  return(list(loglik = loglik, siteLik = siteLik)) # , resll=resll
}



#' Internal maximum likelihood functions.
#'
#' These functions are internally used for the likelihood computations in
#' \code{pml} or \code{optim.pml}.
#'
#' These functions are exported to be used in different packages so far only in
#' the package coalescentMCMC, but are not intended for end user. Most of the
#' functions call C code and are far less forgiving if the import is not what
#' they expect than \code{pml}.
#'
#' @param tree A phylogenetic \code{tree}, object of class \code{phylo}.
#' @param data An alignment, object of class \code{phyDat}.
#' @param bf Base frequencies.
#' @param shape Shape parameter of the gamma distribution.
#' @param k Number of intervals of the discrete gamma distribution.
#' @param Q A vector containing the lower triangular part of the rate matrix.
#' @param levels The alphabet used e.g. c("a", "c", "g", "t") for DNA
#' @param inv Proportion of invariable sites.
#' @param rate Rate.
#' @param g vector of quantiles (default is NULL)
#' @param w vector of probabilities (default is NULL)
#' @param eig Eigenvalue decomposition of Q
#' @param INV Sparse representation of invariant sites
#' @param ll.0 default is NULL
#' @param llMix default is NULL
#' @param wMix default is NULL
#' @param \dots Further arguments passed to or from other methods.
#' @param site return the log-likelihood or vector of sitewise likelihood
#' values
#' @param Mkv indicate if Lewis' Mkv should be estimated.
#' @param site.rate Indicates what type of gamma distribution to use. Options
#' are "gamma" approach of Yang 1994 (default), "quadrature" after the Laguerre
#' quadrature approach of Felsenstein 2001.
## or "lognormal" after a lognormal quadrature approach.
#' @return \code{pml.fit} returns the log-likelihood.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}, \link{pmlPart}, \link{pmlMix}}
#' @references Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a
#' maximum likelihood approach. \emph{Journal of Molecular Evolution},
#' \bold{17}, 368--376.
#' @keywords cluster
#'
#' @rdname pml.fit
#' @export pml.fit
pml.fit <- function(tree, data, bf = rep(1 / length(levels), length(levels)),
                    shape = 1, k = 1, Q = rep(1, length(levels) *
                                                (length(levels) - 1) / 2),
                    levels = attr(data, "levels"), inv = 0, rate = 1, g = NULL,
                    w = NULL, eig = NULL, INV = NULL, ll.0 = NULL, llMix = NULL,
                    wMix = 0, ..., site = FALSE, Mkv = FALSE,
                    site.rate = "gamma") {
  weight <- as.double(attr(data, "weight"))
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  nTips <- as.integer(length(tree$tip.label))
  k <- as.integer(k)
  m <- 1
  if (is.null(eig))
    eig <- edQt(bf = bf, Q = Q)
  if(is.null(g) | is.null(w)){
    rw <- rates_n_weights(shape, k, site.rate)
    g <- rw[, 1]
    w <- rw[, 2]

    if (inv > 0){
      w <- (1 - inv) * w
      g <- g / (1 - inv)
    }
    if (wMix > 0)
      w <- (1 - wMix) * w
    g <- g * rate
  }
  if (any(g < .gEps)) {
    for (i in seq_along(g)) {
      if (g[i] < .gEps) {
        inv <- inv + w[i]
      }
    }
    w <- w[g > .gEps]
    g <- g[g > .gEps]
    #        kmax <- k
    k <- length(w)
  }
  if (is.null(INV))
    INV <- Matrix(lli(data, tree), sparse = TRUE)
  if (is.null(ll.0)) {
    ll.0 <- numeric(attr(data, "nr"))
  }
  if (inv > 0)
    ll.0 <- as.matrix(INV %*% (bf * inv))
  if (Mkv)
    ll.0 <- as.matrix(INV %*% bf)
  if (wMix > 0)
    ll.0 <- ll.0 + llMix

  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  #    root <- as.integer(node[length(node)])
  el <- as.double(tree$edge.length)
  node <- as.integer(node - nTips - 1L) #    min(node))
  edge <- as.integer(edge - 1L)

  contrast <- attr(data, "contrast")
  nco <- as.integer(dim(contrast)[1])
  # dlist=data, nr, nc, weight, k ausserhalb definieren
  # pmlPart einbeziehen
  # as.double(w),
  resll <- .Call("PML0", dlist = data, el, as.double(g), nr, nc, k, eig,
    as.double(bf), node, edge, nTips, nco, contrast,
    N = as.integer(length(edge)))
  # sort(INV@i)+1L
  ind <- which(ll.0 > 0) # automatic in INV gespeichert

  sca <- .Call("rowMax", resll, length(weight), as.integer(k)) + 1
  # nr statt length(weight)
  lll <- resll - sca
  lll <- exp(lll)
  lll <- (lll %*% w)
  if (Mkv) p0 <- sum(exp(log(lll[ind]) + sca[ind]))
  if (inv > 0) lll[ind] <- lll[ind] + exp(log(ll.0[ind]) - sca[ind])
  siteLik <- lll
  siteLik <- log(siteLik) + sca
  # needs to change
  if (wMix > 0) siteLik <- log(exp(siteLik) * (1 - wMix) + llMix)
#  if (Mkv) siteLik <- siteLik - log(1 - p0)
  loglik <- sum(weight * siteLik)
  if (Mkv) loglik <- loglik - sum(weight) * log(1 - p0)
  if (!site) return(loglik)
  resll <- exp(resll)
  return(list(loglik = loglik, siteLik = siteLik, resll = resll))
}

### @param optF3x4 Logical value indicating if codon frequencies are estimated
### for the F3x4 model


#' Likelihood of a tree.
#'
#' \code{pml} computes the likelihood of a phylogenetic tree given a sequence
#' alignment and a model. \code{optim.pml} optimizes the different model
#' parameters.
#'
#' Base frequencies in \code{pml} can be supplied in different ways.
#' For amino acid they are usually defined through specifying a model, so the
#' argument bf does not need to be specified. Otherwise if \code{bf=NULL},
#' each state is given equal probabilty. It can be a numeric vector given the
#' frequencies. Last but not least \code{bf} can be string "equal", "empirical"
#' and for codon models additionally "F3x4".
#'
#' The topology search uses a nearest neighbor interchange (NNI) and the
#' implementation is similar to phyML.  The option model in pml is only used
#' for amino acid models.  The option model defines the nucleotide model which
#' is getting optimised, all models which are included in modeltest can be
#' chosen. Setting this option (e.g. "K81" or "GTR") overrules options optBf
#' and optQ.  Here is a overview how to estimate different phylogenetic models
#' with \code{pml}: \tabular{lll}{ model \tab optBf \tab optQ \cr Jukes-Cantor
#' \tab FALSE \tab FALSE \cr F81 \tab TRUE \tab FALSE \cr symmetric \tab FALSE
#' \tab TRUE \cr GTR \tab TRUE \tab TRUE } Via model in optim.pml the following
#' nucleotide models can be specified: JC, F81, K80, HKY, TrNe, TrN, TPM1, K81,
#' TPM1u, TPM2, TPM2u, TPM3, TPM3u, TIM1e, TIM1, TIM2e, TIM2, TIM3e, TIM3,
#' TVMe, TVM, SYM and GTR.  These models are specified as in Posada (2008).
#'
#' So far 17 amino acid models are supported ("WAG", "JTT", "LG", "Dayhoff",
#' "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb",
#' "FLU", "Blossum62", "Dayhoff_DCMut" and "JTT_DCMut") and additionally rate
#' matrices and amino acid frequencies can be supplied.
#'
#' It is also possible to estimate codon models (e.g. YN98), for details see
#' also the chapter in vignette("phangorn-specials").
#'
#' If the option 'optRooted' is set to TRUE than the edge lengths of rooted
#' tree are optimized. The tree has to be rooted and by now ultrametric!
#' Optimising rooted trees is generally much slower.
#'
#' \code{pml.control} controls the fitting process. \code{epsilon} and
#' \code{maxit} are only defined for the most outer loop, this affects
#' \code{pmlCluster}, \code{pmlPart} and \code{pmlMix}.  \code{epsilon} is
#' defined as (logLik(k)-logLik(k+1))/logLik(k+1), this seems to be a good
#' heuristics which works reasonably for small and large trees or alignments.
#' If \code{trace} is set to zero than no out put is shown, if functions are
#' called internally than the trace is decreased by one, so a higher of trace
#' produces more feedback.
#'
#' If \code{rearrangement} is set to \code{stochastic} a stochastic search
#' algorithm similar to Nguyen et al. (2015). and for \code{ratchet} the
#' likelihood ratchet as in Vos (2003).  This should helps often to find better
#' tree topologies, especially for larger trees.
#'
#' @aliases pml
#' @param tree A phylogenetic \code{tree}, object of class \code{phylo}.
#' @param data An alignment, object of class \code{phyDat}.
#' @param bf Base frequencies (see details).
#' @param Q A vector containing the lower triangular part of the rate matrix.
#' @param inv Proportion of invariable sites.
#' @param k Number of intervals of the discrete gamma distribution.
#' @param shape Shape parameter of the gamma distribution.
#' @param rate Rate.
#' @param model allows to choose an amino acid models or nucleotide model, see
#' details.
#' @param site.rate Indicates what type of gamma distribution to use. Options
#' are "gamma" approach of Yang 1994 (default), "quadrature" after the Laguerre
#' quadrature approach of Felsenstein 2001.
## or "lognormal" after a lognormal
#' @param object An object of class \code{pml}.
#' @param optNni Logical value indicating whether toplogy gets optimized (NNI).
#' @param optBf Logical value indicating whether base frequencies gets
#' optimized.
#' @param optQ Logical value indicating whether rate matrix gets optimized.
#' @param optInv Logical value indicating whether proportion of variable size
#' gets optimized.
#' @param optGamma Logical value indicating whether gamma rate parameter gets
#' optimized.
#' @param optEdge Logical value indicating the edge lengths gets optimized.
#' @param optRate Logical value indicating the overall rate gets optimized.
#' @param optRooted Logical value indicating if the edge lengths of a rooted
#' tree get optimized.
#' @param ratchet.par search parameter for stochastic search
#' @param rearrangement type of tree tree rearrangements to perform, one of
#' "none", "NNI", "stochastic" or "ratchet"
#' @param control A list of parameters for controlling the fitting process.
#' @param subs A (integer) vector same length as Q to specify the optimization
#' of Q
#' @param \dots Further arguments passed to or from other methods.
#' @param epsilon Stop criterion for optimisation (see details).
#' @param maxit Maximum number of iterations (see details).
#' @param trace Show output during optimization (see details).
#' @return \code{pml} or \code{optim.pml} return a list of class \code{pml},
#' some are useful for further computations like \item{tree}{the phylogenetic
#' tree.} \item{data}{the alignment.} \item{logLik}{Log-likelihood of the
#' tree.} \item{siteLik}{Site log-likelihoods.} \item{weight}{Weight of the
#' site patterns.}
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{bootstrap.pml}}, \code{\link{modelTest}},
#' \code{\link{pmlPart}}, \code{\link{pmlMix}}, \code{\link{plot.phylo}},
#' \code{\link{SH.test}}, \code{\link{ancestral.pml}}
#' @references Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a
#' maximum likelihood approach. \emph{Journal of Molecular Evolution},
#' \bold{17}, 368--376.
#'
#' Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer Associates,
#' Sunderland.
#'
#' Yang, Z. (2006). \emph{Computational Molecular evolution}. Oxford University
#' Press, Oxford.
#'
#' Adachi, J., P. J. Waddell, W. Martin, and M. Hasegawa (2000) Plastid genome
#' phylogeny and a model of amino acid substitution for proteins encoded by
#' chloroplast DNA.  \emph{Journal of Molecular Evolution}, \bold{50}, 348--358
#'
#' Rota-Stabelli, O., Z. Yang, and M. Telford. (2009) MtZoa: a general
#' mitochondrial amino acid substitutions model for animal evolutionary
#' studies. \emph{Mol. Phyl. Evol}, \bold{52(1)}, 268--72
#'
#' Whelan, S. and Goldman, N. (2001) A general empirical model of protein
#' evolution derived from multiple protein families using a maximum-likelihood
#' approach. \emph{Molecular Biology and Evolution}, \bold{18}, 691--699
#'
#' Le, S.Q. and Gascuel, O. (2008) LG: An Improved, General Amino-Acid
#' Replacement Matrix \emph{Molecular Biology and Evolution}, \bold{25(7)},
#' 1307--1320
#'
#' Yang, Z., R. Nielsen, and M. Hasegawa (1998) Models of amino acid
#' substitution and applications to Mitochondrial protein evolution.
#' \emph{Molecular Biology and Evolution}, \bold{15}, 1600--1611
#'
#' Abascal, F., D. Posada, and R. Zardoya (2007) MtArt: A new Model of amino
#' acid replacement for Arthropoda. \emph{Molecular Biology and Evolution},
#' \bold{24}, 1--5
#'
#' Kosiol, C, and Goldman, N (2005) Different versions of the Dayhoff rate
#' matrix - \emph{Molecular Biology and Evolution}, \bold{22}, 193--199
#'
#' L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, and B.Q. Minh (2015) IQ-TREE: A
#' fast and effective stochastic algorithm for estimating maximum likelihood
#' phylogenies. \emph{Molecular Biology and Evolution}, \bold{32}, 268--274.
#'
#' Vos, R. A. (2003) Accelerated Likelihood Surface Exploration: The Likelihood
#' Ratchet. \emph{Systematic Biology}, \bold{52(3)}, 368--373
#'
#' Yang, Z., and R. Nielsen (1998) Synonymous and nonsynonymous rate variation
#' in nuclear genes of mammals. \emph{Journal of Molecular Evolution},
#' \bold{46}, 409-418.
#'
#' Lewis, P.O. (2001) A likelihood approach to estimating phylogeny from
#' discrete morphological character data. \emph{Systematic Biology} \bold{50},
#' 913--925.
#' @keywords cluster
#' @importFrom Matrix Matrix
#' @examples
#'
#'   example(NJ)
#' # Jukes-Cantor (starting tree from NJ)
#'   fitJC <- pml(tree, Laurasiatherian)
#' # optimize edge length parameter
#'   fitJC <- optim.pml(fitJC)
#'   fitJC
#'
#' \dontrun{
#' # search for a better tree using NNI rearrangements
#'   fitJC <- optim.pml(fitJC, optNni=TRUE)
#'   fitJC
#'   plot(fitJC$tree)
#'
#' # JC + Gamma + I - model
#'   fitJC_GI <- update(fitJC, k=4, inv=.2)
#' # optimize shape parameter + proportion of invariant sites
#'   fitJC_GI <- optim.pml(fitJC_GI, optGamma=TRUE, optInv=TRUE)
#' # GTR + Gamma + I - model
#'   fitGTR <- optim.pml(fitJC_GI, rearrangement = "stochastic",
#'       optGamma=TRUE, optInv=TRUE, model="GTR")
#' }
#'
#'
#' # 2-state data (RY-coded)
#'   dat <- acgt2ry(Laurasiatherian)
#'   fit2ST <- pml(tree, dat)
#'   fit2ST <- optim.pml(fit2ST,optNni=TRUE)
#'   fit2ST
#' # show some of the methods available for class pml
#'   methods(class="pml")
#'
#' @rdname pml
#' @export pml
pml <- function(tree, data, bf = NULL, Q = NULL, inv = 0, k = 1, shape = 1,
                rate = 1, model = NULL, site.rate = "gamma", ...) {
  Mkv <- FALSE
  if (!is.null(model) && model == "Mkv") Mkv <- TRUE
  call <- match.call()
  extras <- match.call(expand.dots = FALSE)$...
  pmla <- c("wMix", "llMix", "dnds", "tstv")
  existing <- match(pmla, names(extras))
  wMix <- ifelse(is.na(existing[1]), 0,
    eval(extras[[existing[1]]], parent.frame()))
  llMix <- ifelse(is.na(existing[2]), 0,
    eval(extras[[existing[2]]], parent.frame()))
  # allow
  dnds <- tstv <- 1
  dnds <- ifelse(is.na(existing[3]), 1,
    eval(extras[[existing[3]]], parent.frame()))
  tstv <- ifelse(is.na(existing[4]), 1,
    eval(extras[[existing[4]]], parent.frame()))

  if (!inherits(tree, "phylo")) stop("tree must be of class phylo")
  #    if(is.null(tree$edge.length)){
  #        warning("tree has no edge length, used nnls.phylo to assign them")
  #        tree <- nnls.phylo(tree, dist.ml(data))
  #    }
  if (any(duplicated(tree$tip.label))) stop("tree must have unique labels!")
  nTips <- as.integer(length(tree$tip.label))
  if (is.null(attr(tree, "order")) || attr(tree, "order") ==
    "cladewise")
     tree <- reorder(tree, "postorder")
  if (any(tree$edge.length < 0)) {
    if(is.rooted(tree)) nh <- nodeHeight(tree)[1:nTips]
    tree$edge.length[tree$edge.length < 0] <- 1e-08
    message("negative edges length changed to 0!")
    if(is.rooted(tree)){
      ind <- match(as.integer(1:nTips), tree$edge[, 2])
      tree$edge.length[ind] <- tree$edge.length[ind] +
        (nodeHeight(tree)[1:nTips] - nh)
    }
  }
#  if ( class(data)[1] != "phyDat") stop("data must be of class phyDat")
  if (!inherits(data, "phyDat")) stop("data must be of class phyDat")
  if (is.null(tree$edge.length)) stop("tree must have edge weights")
  if (any(is.na(match(tree$tip.label, attr(data, "names")))))
    stop("tip labels are not in data")
  data <- subset(data, tree$tip.label) # needed
  levels <- attr(data, "levels")
  if (Mkv) {
    data <- cbind(constSitePattern(length(tree$tip.label), levels,
      tree$tip.label), data)
    attr(data, "weight")[1:attr(data, "nc")] <- 0.0
  }
  weight <- attr(data, "weight")
  nr <- attr(data, "nr")
  type <- attr(data, "type")
  if (type == "AA" & !is.null(model)) {
    model <- match.arg(model, .aamodels)
    getModelAA(model, bf = is.null(bf), Q = is.null(Q))
  }
  if (type == "CODON") {
    .syn <- synonymous_subs(code=attr(data, "code"))
    .sub <- tstv_subs(code=attr(data, "code"))
    Q <- CodonQ(subs = .sub, syn = .syn, tstv = tstv, dnds = dnds)
    if(is.null(bf)) bf_choice <- "equal"
  }
  if (is.null(bf))
    bf <- rep(1 / length(levels), length(levels))
  if (is.character(bf)) {
    bf_choice <- match.arg(bf, c("equal", "empirical", "F1x4", "F3x4", "F61"))
    if (bf_choice == "F3x4" & type != "CODON")
      stop("F3x4 not available for this data type")
    if (bf_choice == "F1x4" & type != "CODON")
      stop("F1x4 not available for this data type")
    if (bf_choice == "F61" & type != "CODON")
      stop("F61 not available for this data type")
    bf <- switch(bf_choice,
      equal = rep(1 / length(levels), length(levels)),
      empirical = baseFreq(data),
      F61 = baseFreq(data),
      F3x4 = F3x4(data),
      F1x4 = F1x4(data))
    names(bf) <- NULL
  }
  if (type == "CODON") freq_df  <- df_freq_codon(bf_choice)
  if (is.null(Q))
    Q <- rep(1, length(levels) * (length(levels) - 1) / 2)
  m <- 1
  eig <- edQt(bf = bf, Q = Q)

  rw <- rates_n_weights(shape, k, site.rate)
  g <- rw[, 1]
  w <- rw[, 2]
  if (inv > 0){
    w <- (1 - inv) * w
    g <- g / (1 - inv)
  }
  if (wMix > 0)
    w <- (1 - wMix) * w
  g <- g * rate

  inv0 <- inv
  kmax <- k
  if (any(g < .gEps)) {
    for (i in seq_along(g)) {
      if (g[i] < .gEps) {
        inv <- inv + w[i]
      }
    }
    w <- w[g > .gEps]
    g <- g[g > .gEps]
    k <- length(w)
  }

  INV <- Matrix(lli(data, tree), sparse = TRUE)
  ll.0 <- as.matrix(INV %*% (bf * inv))
  if(Mkv) ll.0 <- as.matrix(INV %*% rep(1, length(bf)))

  if (wMix > 0) ll.0 <- ll.0 + llMix

  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))

  on.exit(.Call("ll_free2"))
  .Call("ll_init2", nr, nTips, nc, as.integer(k))
  tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q,
    levels = attr(data, "levels"), inv = inv, rate = rate, g = g, w = w,
    eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix, wMix = wMix, site = TRUE,
    Mkv=Mkv)
  df <- ifelse(is.ultrametric(tree), tree$Nnode, length(tree$edge.length))

  df <- switch(type,
    DNA =  df + (kmax > 1) + (inv0 > 0) + length(unique(bf)) - 1 +
      length(unique(Q)) - 1,
    AA = df + (kmax > 1) + (inv0 > 0),
    CODON = df + (kmax > 1) + (inv0 > 0) + freq_df + (dnds!=1) + (tstv!=1),
    USER = df + (kmax > 1) + (inv0 > 0) + length(unique(bf)) - 1 +
      length(unique(Q)) - 1)

  result <- list(logLik = tmp$loglik, inv = inv, k = kmax, shape = shape,
    Q = Q, bf = bf, rate = rate, siteLik = tmp$siteLik, weight = weight,
    g = g, w = w, eig = eig, data = data, model = model, INV = INV,
    ll.0 = ll.0, tree = tree, lv = tmp$resll, call = call, df = df, wMix = wMix,
    llMix = llMix, Mkv=Mkv, site.rate=site.rate) #
  if (type == "CODON") {
    result$dnds <- dnds
    result$tstv <- tstv
    result$frequencies <- bf_choice
  }
  class(result) <- "pml"
  result
}


optimRooted <- function(tree, data, eig = eig, w = w, g = g, bf = bf,
                        rate = rate, ll.0 = ll.0, INV = INV,
                        control = pml.control(epsilon = 1e-08, maxit = 25,
                                              trace = 0), ...) {
  tree$edge.length[tree$edge.length < 1e-08] <- 1e-08 # nicht richtig
  nTips <- as.integer(length(tree$tip.label))
  k <- length(w)

  # optimising rooted triplets
  optRoot0 <- function(t, tree, data, g, w, eig, bf, ll.0, k) {
    l <- length(tree$edge.length)
    tree$edge.length[1:(l - 1)] <- tree$edge.length[1:(l - 1)] + t
    tree$edge.length[l] <- tree$edge.length[l] - t
    loglik <- pml.fit4(tree, data, bf = bf, g = g, w = w, eig = eig, INV = INV,
      ll.0 = ll.0, k = k) #
    loglik
  }
  # optim edges leading to the root
  optRoot2 <- function(t, tree, data, g, w, eig, bf, ll.0, k) {
    tree$edge.length <- tree$edge.length + t   # c(el1+t, el2-t)
    loglik <- pml.fit4(tree, data, bf = bf, g = g, w = w, eig = eig, INV = INV,
      ll.0 = ll.0, k = k) # , INV=INV
    loglik
  }
  # scale whole tree
  scaleEdges <- function(t = 1, trace = 0, tree, data, ...) {
    fn <- function(t, tree, data, ...) {
      tree$edge.length <- tree$edge.length * t
      pml.fit4(tree, data, ...)
    }
    optimize(f = fn, interval = c(0.25, 4), tree = tree, data = data, ...,
      maximum = TRUE, tol = .00001)
  }
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]

  anc <- Ancestors(tree, 1:max(tree$edge), "parent")
  sibs <- Siblings(tree, 1:max(tree$edge))
  allKids <- cvector <- allChildren(tree)
  rootNode <- getRoot(tree)

  child2 <- orderNNI(tree, nTips)   # (cvector, rootNode, nTips, TRUE)

  lengthList <- edgeList <- vector("list", max(tree$edge))
  for (i in tree$edge[, 2]) {
    pa <- anc[i]
    kids <- sibs[[i]]
    if (pa != rootNode) {
      edgeList[[i]] <- cbind(pa, c(anc[pa], kids))
      lengthList[[i]] <- c(pa, kids)
    }
    else {
      edgeList[[i]] <- cbind(pa, kids)
      lengthList[[i]] <- kids
    }
  }

  treeList <- vector("list", max(child2))
  for (i in child2) {
    pa <- anc[i]
    kids <- allKids[[i]]
    treeList[[i]] <- cbind(i, c(kids, pa))
  }

  ll <- pml.fit4(tree, data, bf = bf,  k = k, eig = eig, ll.0 = ll.0, INV = INV,
                 w = w, g = g)
  #    if(control$trace>2)cat("ll", ll, "\n")
  eps <- 10
  iter <- 1

  EL <- numeric(max(tree$edge))
  EL[tree$edge[, 2]] <- tree$edge.length

  eps0 <- 1e-8

  tmp <- scaleEdges(t, trace = 0, tree, data, bf = bf, k = k, ll.0 = ll.0,
    eig = eig, w = w, g = g)
  #    if(control$trace>2)cat("scale", tmp[[2]], "\n")
  t <- tmp[[1]]
  tree$edge.length <- tree$edge.length * t
  el <- tree$edge.length
  EL[tree$edge[, 2]] <- tree$edge.length
  ll2 <- pml.fit4(tree, data, bf = bf,  k = k, eig = eig, INV = INV,
                  ll.0 = ll.0, w = w, g = g)
  tmptree <- tree

  while (eps > control$eps && iter < control$maxit) {
    ll2 <- pml.fit4(tree, data, bf = bf,  k = k, eig = eig, INV = INV,
                    ll.0 = ll.0, w = w, g = g)
    loli <- rootNode

    children <- allKids[[rootNode]]
    kidsEl <- EL[children]
    minEl <- min(kidsEl)
    kidsEl <- kidsEl - minEl

    tmptree$edge <- cbind(rootNode, children)
    tmptree$edge.length <- kidsEl

    t <- optimize(f = optRoot2, interval = c(1e-8, 3), tmptree, data = data,
                  k = k, g = g, w = w, eig = eig, bf = bf, ll.0 = ll.0,
                  maximum = TRUE)
    optRoot2(t[[1]], tmptree, data = data, k = k, g = g, w = w, eig = eig,
             bf = bf, ll.0 = ll.0)
    # if(control$trace>2)cat("optRoot", t[[2]], "\n")
    ll3 <- t[[2]]
    EL[children] <- kidsEl + t[[1]]

    tree$edge.length <- EL[tree$edge[, 2]]
    ll2 <- pml.fit4(tree, data, bf = bf, k = k, eig = eig, INV = INV,
                    ll.0 = ll.0, w = w, g = g)
    for (i in seq_along(child2)) {
      dad <- child2[i]
      pa <- anc[dad]
      while (loli != pa) {
        tmpKids <- cvector[[loli]]
        tmpEdge <- cbind(loli, tmpKids)
        pml.move(tmpEdge, EL[tmpKids], data, g, w, eig, k, nTips, bf)
        loli <- anc[loli]
      }
      pml.move(edgeList[[dad]], EL[lengthList[[dad]]], data, g, w, eig, k,
        nTips, bf)

      children <- allKids[[dad]]
      kidsEl <- EL[children]
      minEl <- min(kidsEl)
      maxEl <- EL[dad]
      EDGE <- treeList[[dad]]

      tmptree$edge <- EDGE
      tmptree$edge.length <- c(kidsEl, maxEl)

      t0 <- optRoot0(0, tmptree, data, g, w, eig, bf, ll.0, k)

      t <- optimize(f = optRoot0, interval = c(-minEl + eps0, maxEl - eps0),
        tmptree, data = data, g = g, w = w, eig = eig, bf = bf,
        ll.0 = ll.0, k = k, maximum = TRUE)
      # if(control$trace>2) cat("edge", t[[2]], "\n")
      if (!is.nan(t[[2]]) & t[[2]] > ll3) {
        optRoot0(t[[1]], tmptree, data = data, g = g, w = w, eig = eig, bf = bf,
          ll.0 = ll.0, k = k)
        EL[children] <- kidsEl + t[[1]]
        EL[dad] <- maxEl - t[[1]]
        ll3 <- t[[2]]
      }
      else optRoot0(0, tmptree, data, g, w, eig, bf, ll.0, k)
      loli <- dad
      #            }
    }
    tree$edge.length <- EL[tree$edge[, 2]]
    ll2 <- pml.fit4(tree, data, bf = bf, k = k, eig = eig, ll.0 = ll.0,
                   INV = INV, w = w, g = g)
    eps <- (ll - ll2) / ll2

    if (control$trace > 1) cat("optimRooted: ", ll, " -> ", ll2, "\n")
    ll <- ll2
    iter <- iter + 1
  }
  list(tree = tree, logLik = ll, c(eps = eps, iter = iter))
}


index.nni <- function(ch, cvector, pvector, root) {
  p1 <- pvector[ch]
  k12 <- cvector[[ch]]
  k3 <- cvector[[p1]]
  k3 <- k3[k3 != ch]
  kids <- c(k12, k3, ch)
  parents <- c(ch, ch, p1, p1)
  if (p1 != root) {
    k4 <- pvector[p1]
    kids <- c(kids, k4)
    parents <- c(parents, p1)
  }
  cbind(parents, kids)
}


orderNNI <- function(tree, nTips) {
  res <- reorder(tree)$edge[, 2]
  res <- res[res > nTips]
  res
}


rooted.nni <- function(tree, data, eig, w, g, bf, rate, ll.0, INV,
                       control = pml.control(epsilon = 1e-08, maxit = 25,
                                             trace = 0), ...) {
  ind0 <- which(ll.0 > 0)
  contrast <- attr(data, "contrast")
  tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
  nTips <- as.integer(length(tree$tip.label))
  k <- length(w)
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorder.phylo(tree, "postorder")
  if (!is.rooted(tree)) stop("tree must be rooted")

  attr(tree, "order") <- NULL
  weight <- attr(data, "weight")
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))

  getEL1 <- function(t, nh) {
    el <- numeric(4)
    if (nh[1] > nh[2]) {
      el[2] <- nh[1] - nh[2]
      tnh <- nh[1] + t[1]
    }
    else {
      el[1] <- nh[2] - nh[1]
      tnh <- nh[2] + t[1]
    }
    el[1:2] <- el[1:2] + t[1]
    if (tnh > nh[3]) el[3] <- el[3] + tnh - nh[3]
    else el[4] <- el[4] - tnh + nh[3]
    el[3:4] <- el[3:4] + t[2]
    el
  }

  optRootU <- function(t, tree, data, bf, g, w, eig, ll.0, k, INV, nh) {
    tree$edge.length <- getEL1(t, nh)
    pml.fit4(tree, data, bf = bf,  k = k, eig = eig, ll.0 = ll.0, INV = INV,
             w = w, g = g)
  }

  getEL2 <- function(t, nh) {
    el <- numeric(5)
    eps <- 1e-6
    nh12.min <- max(nh[1:2]) + eps
    nh123.min <- max(nh12.min, nh[3]) + eps
    l1 <- nh[5] - nh123.min -  eps
    el[5] <- l1 * t[1] + eps
    nh123 <- nh[5] - el[5]
    l2 <- nh123 - nh12.min - eps
    nh12 <- nh12.min + l2 * t[2]
    el[1] <- nh12 - nh[1]
    el[2] <- nh12 - nh[2]
    el[3] <- nh123 - nh[3]
    el[4] <- nh123 - nh12
    el
  }

  optEdgeU <- function(t, tree, data, bf, g, w, eig, ll.0, k, INV, nh) {
    tree$edge.length <- getEL2(t, nh)
    pml.fit4(tree, data, bf = bf,  k = k, eig = eig, ll.0 = ll.0, INV = INV,
             w = w, g = g)
  }

  child <- tree$edge[, 2]
  parent <- tree$edge[, 1]
  ll <- pml.fit4(tree, data, bf = bf,  k = k, eig = eig, ll.0 = ll.0, INV = INV,
                 w = w, g = g)
  llstart <- ll
  eps <- .00001
  iter <- 1
  EL <- numeric(max(tree$edge))
  EL[tree$edge[, 2]] <- tree$edge.length
  change <- numeric(length(parent)) + 1
  rootNode <- getRoot(tree)
  anc <- Ancestors(tree, 1:max(tree$edge), "parent")
  cvector <- allChildren(tree)
  sibs <- Siblings(tree, 1:max(tree$edge))

  child2 <- orderNNI(tree, nTips)

  while (iter < 2) {
    ll2 <-  pml.fit(tree, data, bf = bf, k = k, eig = eig, ll.0 = ll.0,
                    INV = INV, w = w, g = g)
    nh <- nodeHeight(tree)
    loli <- rootNode
    pa <- rootNode
    nchanges <- 0
    ind <- 1
    i <- 1

    tree1 <- tree2 <- tree3 <- tree
    for (i in seq_along(child2)) {
      ch <- child2[i]
      dad <- anc[ch]
      if (ch > nTips) {

        EL[tree$edge[, 2]] <- tree$edge.length

        pa <- ifelse(dad == rootNode, rootNode, anc[dad])
        # should avoid unnecessary movements
        while (loli != dad && loli != rootNode) {
          if (loli == pa) {
            tmpKids <- sibs[[dad]]
            tmpEdge <- cbind(pa, c(anc[pa], tmpKids))
            pml.move(tmpEdge, EL[c(pa, tmpKids)], data, g, w, eig,
              k, nTips, bf)
            # cat("move from pa to dad \n")
            loli <- dad
          }
          else {
            #   cat("move loli up", loli, "dad", dad, "pa", pa, "ch", ch, "\n")
            tmpKids <- cvector[[loli]]
            tmpEdge <- cbind(loli, tmpKids)
            pml.move(tmpEdge, EL[tmpKids], data, g, w, eig, k,
              nTips, bf)
            loli <- anc[loli]
          }

        }

        if (loli == rootNode && dad != loli) {
          # update all nodes
          pml.fit(tree, data, bf = bf, k = k, eig = eig, ll.0 = ll.0, INV = INV,
            w = w, g = g)
          #   cat("move down loli", loli, "dad", dad, "pa", pa, "ch", ch, "\n")
          gd <- rev(Ancestors(tree, ch, "all"))

          tmpKids <- sibs[[gd[2]]]
          tmpEdge <- cbind(rootNode, tmpKids)
          pml.move(tmpEdge, EL[tmpKids], data, g, w, eig, k, nTips, bf)
          gd <- gd[-1]
          while (length(gd) > 1) {
            tmpKids <- sibs[[gd[2]]]
            tmpEdge <- cbind(gd[1], c(anc[gd[1]], tmpKids))
            pml.move(tmpEdge, EL[c(gd[1], tmpKids)], data, g, w, eig,
              k, nTips, bf)
            gd <- gd[-1]
          }
          loli <- dad
        }

        X1 <- index.nni(ch, cvector, anc, rootNode)

        if (loli != rootNode) {
          tree1$edge <- X1
          tree1$edge.length <- abs(nh[X1[, 1]] - nh[X1[, 2]])
          ll0 <- pml.fit4(tree1, data, bf = bf,  k = k, eig = eig,
            ll.0 = ll.0, INV = INV, w = w, g = g)
          #      cat("quartet", ll0, ch, dad, "\n")
        }


        if (dad == rootNode) {

          ll0 <- pml.fit(tree, data, bf = bf, g = g, w = w, eig = eig,
            ll.0 = ll.0, k = k, INV = INV)
          #   cat("at root", ll0, ch, dad, "\n")
          ind2 <- c(1, 3, 2, 4)
          ind3 <- c(3, 2, 1, 4)
          X2 <- X3 <- X1
          X2[, 2] <- X1[ind2, 2]
          X3[, 2] <- X1[ind3, 2]

          tree1$edge <- X1
          tree2$edge <- X2
          tree3$edge <- X3
          edge1 <- X1[, 2]
          edge1[4] <- dad
          res1 <- optim(par = c(.1, .1), optRootU, gr = NULL, tree = tree1,
            data = data, nh = nh[X1[, 2]], g = g, w = w, eig = eig, bf = bf,
            ll.0 = ll.0, INV = INV, k = k, method = "L-BFGS-B",
            lower = 1e-8, upper = 5, control = list(fnscale = -1))
          res2 <- optim(par = c(.1, .1), optRootU, gr = NULL, tree = tree2,
            data = data, nh = nh[X2[, 2]], g = g, w = w, eig = eig, bf = bf,
            ll.0 = ll.0, INV = INV, k = k, method = "L-BFGS-B",
            lower = 1e-8, upper = 5, control = list(fnscale = -1))
          res3 <- optim(par = c(.1, .1), optRootU, gr = NULL, tree = tree3,
            data = data,  nh = nh[X3[, 2]], g = g, w = w, eig = eig, bf = bf,
            ll.0 = ll.0, INV = INV, k = k, method = "L-BFGS-B",
            lower = 1e-8, upper = 5, control = list(fnscale = -1))
          ind <- which.max(c(res1[[2]], res2[[2]], res3[[2]]))
          if (control$trace > 2) cat("root", c(res1[[2]], res2[[2]],
              res3[[2]]), "\n")

          if (ind == 1) {
            ll2 <- res1[[2]]
            optRootU(t = res1[[1]], tree = tree1, data = data,
              nh = nh[X1[, 2]], g = g, w = w, eig = eig, bf = bf,
              ll.0 = ll.0, INV = INV, k = k)
            tmpEL <- getEL1(res1[[1]], nh[X1[, 2]])
            tree <- changeEdgeLength(tree, X1[, 2], tmpEL)
          }
          if (ind == 2) {
            ll2 <- res2[[2]]
            optRootU(t = res2[[1]], tree = tree2, data = data,
              nh = nh[X2[, 2]], g = g, w = w, eig = eig, bf = bf,
              ll.0 = ll.0, INV = INV, k = k)
            tmpEL <- getEL1(res2[[1]], nh[X2[, 2]])
            tree <- changeEdge(tree, X1[c(2, 3), 2])
            tree <- changeEdgeLength(tree, X2[, 2], tmpEL)
          }
          if (ind == 2) {
            ll2 <- res3[[2]]
            optRootU(t = res3[[1]], tree = tree3, data = data,
              nh = nh[X3[, 2]], g = g, w = w, eig = eig, bf = bf,
              ll.0 = ll.0, INV = INV, k = k)
            tmpEL <- getEL1(res3[[1]], nh[X3[, 2]])
            tree <- changeEdge(tree, X1[c(1, 3), 2])
            tree <- changeEdgeLength(tree, X3[, 2], tmpEL)
          }
        }
        else {
          loli <- dad
          ind2 <- c(1, 3, 2, 4, 5)
          ind3 <- c(3, 2, 1, 4, 5)
          X2 <- X3 <- X1
          X2[, 2] <- X1[ind2, 2]
          X3[, 2] <- X1[ind3, 2]
          tree1$edge <- X1
          tree2$edge <- X2
          tree3$edge <- X3
          tt <- c(.3, .5)

          res1 <- optim(par = tt, optEdgeU, gr = NULL, tree = tree1, data,
            nh = nh[X1[, 2]], g = g, w = w, eig = eig, bf = bf, ll.0 = ll.0,
            INV = INV, k = k, method = "L-BFGS-B", lower = 1e-4,
            upper = 1 - 1e-4, control = list(fnscale = -1))
          res2 <- optim(par = tt, optEdgeU, gr = NULL, tree = tree2, data,
            nh = nh[X2[, 2]], g = g, w = w, eig = eig, bf = bf, ll.0 = ll.0,
            INV = INV, k = k, method = "L-BFGS-B", lower = 1e-4,
            upper = 1 - 1e-4, control = list(fnscale = -1))
          res3 <- optim(par = tt, optEdgeU, gr = NULL, tree = tree3, data,
            nh = nh[X3[, 2]], g = g, w = w, eig = eig, bf = bf, ll.0 = ll.0,
            INV = INV, k = k, method = "L-BFGS-B", lower = 1e-4,
            upper = 1 - 1e-4, control = list(fnscale = -1))

          ind <- which.max(c(res1[[2]], res2[[2]], res3[[2]]))
          if (control$trace > 2) cat("edge", ch, ":", c(res1[[2]], res2[[2]],
              res3[[2]]), "\n")
          ll3 <- max(c(res1[[2]], res2[[2]], res3[[2]]))

          if ( (ll3 - 1e-5 * ll3) < ll2) {
            loli <- rootNode
            ll2 <- pml.fit4(tree, data, bf = bf,  k = k, eig = eig, ll.0 = ll.0,
                            INV = INV, w = w, g = g)
            nh <- nodeHeight(tree)
            EL[tree$edge[, 2]] <- tree$edge.length
            ind <- 0
          }
          else {
            if (ind == 1) {
              ll2 <- res1[[2]]
              optEdgeU(res1[[1]], tree = tree1, data, nh = nh[X1[, 2]], g = g,
                w = w, eig = eig, bf = bf, ll.0 = ll.0, INV = INV, k = k)
              tmpEL <- getEL2(res1[[1]], nh[X1[, 2]])
              tmpE <- X1[, 2]
              tmpE[5] <- X1[5, 1]
              tree <- changeEdgeLength(tree, tmpE, tmpEL)
            }
            if (ind == 2) {
              ll2 <- res2[[2]]
              optEdgeU(res2[[1]], tree = tree2, data, nh = nh[X2[, 2]], g = g,
                w = w, eig = eig, bf = bf, ll.0 = ll.0, INV = INV, k = k)
              tmpEL <- getEL2(res2[[1]], nh[X2[, 2]])
              tmpE <- X2[, 2]
              tmpE[5] <- X1[5, 1]
              tree <- changeEdge(tree, X1[c(2, 3), 2])
              tree <- changeEdgeLength(tree, tmpE, tmpEL)
            }
            if (ind == 3) {
              ll2 <- res3[[2]]
              optEdgeU(res3[[1]], tree = tree3, data, nh = nh[X3[, 2]], g = g,
                w = w, eig = eig, bf = bf, ll.0 = ll.0, INV = INV, k = k)
              tmpEL <- getEL2(res3[[1]], nh[X3[, 2]])
              tmpE <- X3[, 2]
              tmpE[5] <- X1[5, 1]
              tree <- changeEdge(tree, X1[c(1, 3), 2])
              tree <- changeEdgeLength(tree, tmpE, tmpEL)
            }

          }
        }
        nh <- nodeHeight(tree)
        EL[tree$edge[, 2]] <- tree$edge.length
        loli <- dad

        if (ind > 1) {
          # print("NNI swap")
          nchanges <- nchanges + 1
          anc <- Ancestors(tree, 1:max(tree$edge), "parent")
          cvector <- allChildren(tree)
          sibs <- Siblings(tree, 1:max(tree$edge))
        }
      }

    }
    ll2 <- pml.fit4(tree, data, bf = bf,  k = k, eig = eig, ll.0 = ll.0,
                    INV = INV, w = w, g = g)
    eps <- (ll - ll2) / ll2
    if (control$trace > 1) cat(ll, " -> ", ll2, "\n")
    if (control$trace > 1) cat("swap:", nchanges)
    ll <- ll2
    iter <- iter + 1
  }
  list(tree = tree, logLik = ll, iter = iter, swap = nchanges)
}


updateRates <- function(res, ll, rate, shape, k, inv, wMix, update="rate",
                        site.rate = "gamma"){
  if( is.infinite(res[[2]]) || is.nan(res[[2]])) return(NULL)
  if(res[[2]] < ll) return(NULL)
  update <- match.arg(update, c("rate", "shape", "inv"))
  if(update=="rate") rate <- res[[1]]
  if(update=="shape") shape <- res[[1]]
  if(update=="inv") inv <- res[[1]]

  rw <- rates_n_weights(shape, k, site.rate)
  g <- rw[, 1]
  w <- rw[, 2]

  if (inv > 0){
    w <- (1 - inv) * w
    g <- g / (1 - inv)
  }
  if (wMix > 0)
    w <- (1 - wMix) * w
  g <- g * rate

  assign("g", g, envir = parent.frame(n = 1))
  assign("w", w, envir = parent.frame(n = 1))
  assign("inv", inv, envir = parent.frame(n = 1))
  assign("rate", rate, envir = parent.frame(n = 1))
  assign("shape", shape, envir = parent.frame(n = 1))
  assign("ll", res[[2]], envir = parent.frame(n = 1))
}


#' @rdname pml
#' @aliases optim.pml
#' @export
optim.pml <- function(object, optNni = FALSE, optBf = FALSE, optQ = FALSE,
                      optInv = FALSE, optGamma = FALSE, optEdge = TRUE,
                      optRate = FALSE,  optRooted = FALSE, #optF3x4 = FALSE,
                      control = pml.control(epsilon = 1e-8, maxit = 10,
                                            trace = 1L), model = NULL,
                      rearrangement = ifelse(optNni, "NNI", "none"),
                      subs = NULL, ratchet.par = list(iter = 20L, maxit = 100L,
                                                      prop = 1 / 2), ...) {
  optRatchet <- FALSE
  optRatchet2 <- FALSE
  optF3x4 <- FALSE
  if (rearrangement ==  "none") {
    optNni <- FALSE
    optRatchet <- FALSE
    optRatchet2 <- FALSE
  }
  if (rearrangement ==  "NNI") optNni <- TRUE
  if (rearrangement ==  "stochastic") optRatchet <- TRUE
  if (rearrangement ==  "ratchet") optRatchet2 <- TRUE
  extras <- match.call(expand.dots = FALSE)$...
  pmla <- c("wMix", "llMix")
  wMix <- object$wMix
  llMix <- object$llMix
  Mkv <- object$Mkv
#  Mkv <- FALSE
  site.rate <- object$site.rate
  if (is.null(llMix)) llMix <- 0
  if (!is.null(extras)) {
    names(extras) <- pmla[pmatch(names(extras), pmla)]
    existing <- match(pmla, names(extras))
    if (!is.na(existing[1]))
      wMix <- eval(extras[[existing[1]]], parent.frame())
    if (!is.na(existing[2]))
      llMix <- eval(extras[[existing[2]]], parent.frame())
  }
  tree <- object$tree
  call <- object$call
  ratchet <- FALSE
  ratchet2 <- FALSE
  if (optRatchet == TRUE) {
    if (optRooted == FALSE) {
      optNni <- TRUE
      optEdge <- TRUE
      ratchet <- TRUE
    }
  }
  if (optRatchet2 == TRUE) {
    optNni <- TRUE
    optEdge <- TRUE
    ratchet2 <- TRUE
  }

  data <- object$data
  addTaxa <- FALSE

  if (optNni) {
    mapping <- map_duplicates(data)
    if (!is.null(mapping)) {
      orig.data <- data
      addTaxa <- TRUE
      tree2 <- drop.tip(tree, mapping[, 1])
      tree <- reorder(tree2, "postorder")
    }
    if (!is.binary(tree))
      tree <- multi2di(tree)
    optEdge <- TRUE
  }
  if (length(tree$tip.label) < (3 + !optRooted)) {
    optNni <- FALSE
    ratchet <- FALSE
    ratchet2 <- FALSE
  }
  if (length(tree$tip.label) < (2 + !optRooted)) {
    stop("rooted / unrooted tree needs at least 2 / 3 tips")
  }
  if (is.rooted(tree)) {
    if (optRooted == FALSE && optEdge == TRUE) {
      tree <- unroot(tree)
      attr(tree, "order") <- NULL
      tree <- reorder(tree, "postorder")
      warning("I unrooted the tree", call. = FALSE)
    }
  }
  if (is.null(attr(tree, "order")) || attr(tree, "order") ==
    "cladewise")
    tree <- reorder(tree, "postorder")
  if (any(tree$edge.length < 1e-08)) {
    tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
    if (optRooted) {
      # ensure tree is ultrametric
      nTips <- as.integer(length(tree$tip.label))
      ind <- match(as.integer(1:nTips), tree$edge[, 2])
      tree$edge.length[ind] <- tree$edge.length[ind] +
        nodeHeight(tree)[1:nTips]
    }
    object <- update.pml(object, tree = tree)
  }
  if (optEdge & optRate) {
    warning("You can't optimise edges and rates at the same time, only edges are optimised!", call. = FALSE)
    optRate <- FALSE
  }
  if (optRooted) {
    optEdge <- FALSE
    if (!is.rooted(tree)) stop("tree must be rooted")
    if (!is.ultrametric(tree)) stop("Tree must be ultrametric!")
  }
  trace <- control$trace
  data <- subset(data, tree$tip.label)
  type <- attr(data, "type")
  if (type == "AA") {
    if(!is.null(model)) object <- update(object, model = model)
    model <- object$model
  }
  if (type == "CODON") {
    .syn <- synonymous_subs(code=attr(data, "code"))
    .sub <- tstv_subs(code=attr(data, "code"))
    if (is.null(model)) model <- "codon1"
    model <- match.arg(model, c("codon0", "codon1", "codon2", "codon3", "YN98"))
    dnds <- object$dnds
    tstv <- object$tstv
    if (!is.null(model)) {
      if (model == "codon0") optQ <- FALSE
      else  optQ <- TRUE
      if (model == "YN98") optBf <- TRUE
    }
    bf_choice <- object$frequencies
    freq_df <- df_freq_codon(bf_choice)
    if(bf_choice=="F3x4" & optBf) optF3x4 <- TRUE
  }
  if (optF3x4) {
    if (type != "CODON") stop("optF3x4 needs codon data")
    optBf <- FALSE
    bf <- F3x4(data)
    bf_codon <- bf_by_codon(data)
    object <- update.pml(object, bf = bf)
  }
  Q <- object$Q
  if (is.null(subs)) subs <- c(1:(length(Q) - 1), 0)
  bf <- object$bf
  eig <- object$eig
  inv <- object$inv
  k <- object$k
  if (k == 1 & optGamma) {
    optGamma <- FALSE
    message("only one rate class, ignored optGamma")
  }
  if(Mkv==TRUE & optInv==TRUE){
    optInv <- FALSE
    message('cannot estimate invariant sites and Mkv model, ignored optInv')
  }
  shape <- object$shape
  w <- object$w
  g <- object$g
  if (type == "DNA" & !is.null(model)) {
    tmp <- subsChoice(model)
    optQ <- tmp$optQ
    if (!optQ)
      Q <- rep(1, 6)
    optBf <- tmp$optBf
    if (!optBf)
      bf <- c(0.25, 0.25, 0.25, 0.25)
    subs <- tmp$subs
  }
  ll0 <- object$logLik
  INV <- object$INV
  ll.0 <- object$ll.0
  rate <- object$rate
  ll <- ll0
  ll1 <- ll0
  opti <- TRUE

  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  nTips <- as.integer(length(tree$tip.label))
  on.exit({
    tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q,
      levels = attr(data, "levels"), inv = inv, rate = rate,
      g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix,
      wMix = wMix, site = TRUE)
    if (addTaxa) {
      tree <- add.tips(tree, tips = mapping[, 1], where = mapping[, 2],
        edge.length = rep(0, nrow(mapping)))
      data <- orig.data
    }
    df <- ifelse(optRooted, tree$Nnode, length(tree$edge.length))
    df <- switch(type,
      DNA = df + (k > 1) + (inv > 0) + length(unique(bf)) - 1 +
        length(unique(Q)) - 1,
      AA = df + (k > 1) + (inv > 0) +  optBf * (length(unique(bf)) - 1),
      CODON = df + (k > 1) + (inv > 0) + freq_df + (dnds != 1) + (tstv != 1),
      USER = df + (k > 1) + (inv > 0) + length(unique(bf)) - 1 +
        length(unique(Q)) - 1)

    object <- list(logLik = tmp$loglik, inv = inv, k = k, shape = shape,
      Q = Q, bf = bf, rate = rate, siteLik = tmp$siteLik,
      weight = attr(data, "weight"),
      g = g, w = w, eig = eig, data = data, model = model,
      INV = INV, ll.0 = ll.0, tree = tree, lv = tmp$resll,
      call = call, df = df, wMix = wMix, llMix = llMix, Mkv=Mkv,
      site.rate=site.rate)
    if (type == "CODON") {
      object$dnds <- dnds
      object$tstv <- tstv
      object$frequencies <- bf_choice
    }
    class(object) <- "pml"

    extras <- pairlist(bf = bf, Q = Q, inv = inv, shape = shape,
      rate = rate)[c(optBf, optQ, optInv, optGamma, optRate)]
    if (length(extras)) {
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if (any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }
    object$call <- call

    pml.free()
    return(object)
  })
  pml.init(data, k)

  if (optEdge) {
    # check if non-negative least-squares is better for start of
    # optimisation
#    treetmp <- nnls.phylo(tree, dist.ml(data))
#    treetmp$edge.length[treetmp$edge.length < 1e-8] <- 1e-8
#    tmplogLik <- pml.fit(treetmp, data, bf, k = k, inv = inv, g = g, w = w,
#      eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix, wMix = wMix)
#    if (tmplogLik > ll) tree <- treetmp
    res <- optimEdge(tree, data, eig = eig, w = w, g = g, bf = bf, rate = rate,
      ll.0 = ll.0, INV = INV,
      control <- pml.control(epsilon = 1e-07, maxit = 5, trace = trace))
#    if (trace > 0)
#      cat("optimize edge weights: ", ll, "-->", max(res[[2]], ll), "\n")
    if (res[[2]] > ll) {
      ll <- res[[2]]
      tree <- res[[1]]
    }
  }
  if (optRooted) {
    res <- optimRooted(tree, data, eig = eig, w = w, g = g, bf = bf,
                       rate = rate, ll.0 = ll.0, INV = INV,
                       control = pml.control(epsilon = 1e-07, maxit = 10,
                                             trace = trace - 1))
    if (trace > 0)
      cat("optimize edge weights: ", ll, "-->", max(res[[2]], ll), "\n")
    if (res[[2]] > ll) {
      ll <- res[[2]]
      tree <- res[[1]]
    }
  }
  rounds <- 1
  while (opti) {
    if (optBf) {
      res <- optimBf(tree, data, bf = bf, inv = inv, Q = Q, w = w, g = g,
        INV = INV, rate = rate, k = k, llMix = llMix)
      if (trace > 0)
        cat("optimize base frequencies: ", ll, "-->",
          max(res[[2]], ll), "\n")
      if (res[[2]] > ll) {
        bf <- res[[1]]
        eig <- edQt(Q = Q, bf = bf)
        if (inv > 0)
          ll.0 <- as.matrix(INV %*% (bf * inv))
        if (wMix > 0)
          ll.0 <- ll.0 + llMix
        ll <- res[[2]]
      }
    }
    if (optF3x4) {
      res <- optimF3x4(tree, data, bf_codon = bf_codon, inv = inv, Q = Q, w = w,
        g = g, INV = INV, rate = rate, k = k, llMix = llMix)
      if (trace > 0)
        cat("optimize base frequencies: ", ll, "-->", max(res[[2]], ll), "\n")
      if (res[[2]] > ll) {
        bf <- res[[1]]
        eig <- edQt(Q = Q, bf = bf)
        if (inv > 0)
          ll.0 <- as.matrix(INV %*% (bf * inv))
        if (wMix > 0)
          ll.0 <- ll.0 + llMix
        ll <- res[[2]]
      }
    }
    if (optQ) {
      if (type == "CODON") {
        ab <- c(tstv, dnds)
        res <- switch(model,
          YN98 = optimCodon(tree, data, Q = rep(1, 1830), subs = .sub,
            syn = .syn, bf = bf, w = w, g = g, inv = inv,
            INV = INV, ll.0 = ll.0, rate = rate, k = k,
            ab = log(ab), optK = TRUE, optW = TRUE),
          codon1 = optimCodon(tree, data, Q = rep(1, 1830), subs = .sub,
            syn = .syn, bf = bf, w = w, g = g, inv = inv,
            INV = INV, ll.0 = ll.0, rate = rate, k = k,
            ab = log(ab), optK = TRUE, optW = TRUE),
          codon2 = optimCodon(tree, data, Q = rep(1, 1830), subs = .sub,
            syn = .syn, bf = bf, w = w, g = g, inv = inv,
            INV = INV, ll.0 = ll.0, rate = rate, k = k,
            ab = log(ab), optK = FALSE, optW = TRUE),
          codon3 = optimCodon(tree, data, Q = rep(1, 1830), subs = .sub,
            syn = .syn, bf = bf, w = w, g = g, inv = inv,
            INV = INV, ll.0 = ll.0, rate = rate, k = k,
            ab = log(ab), optK = TRUE, optW = FALSE))
        tmp <- res[[5]]
        m <- length(tmp)
        dnds <- tmp[m]

        if (m > 1) tstv <- tmp[1]
      }
      else{
        res <- optimQ(tree, data, Q = Q, subs = subs, bf = bf, w = w, g = g,
                      inv = inv, INV = INV, ll.0 = ll.0, rate = rate, k = k)
      }
      Q <- res[[1]]
      eig <- edQt(Q = Q, bf = bf)
      if (trace > 0) cat("optimize rate matrix: ", ll, "-->", res[[2]], "\n")
      ll <- res[[2]]
    }
    if (optInv) {
      res <- optimInv(tree, data, inv = inv, INV = INV, Q = Q,
        bf = bf, eig = eig, k = k, shape = shape, rate = rate)
      if (trace > 0)
        cat("optimize invariant sites: ", ll, "-->", max(res[[2]], ll), "\n")
      updateRates(res, ll, rate, shape, k, inv, wMix, update="inv",
                  site.rate=site.rate)
      ll.0 <- as.matrix(INV %*% (bf * inv))
      if (wMix > 0) ll.0 <- ll.0 + llMix
    }
    if (optGamma) {
      res <- optimGamma(tree, data, shape = shape, k = k, inv = inv, INV = INV,
                        Q = Q, bf = bf, eig = eig, ll.0 = ll.0, rate = rate)
      if (trace > 0)
        cat("optimize shape parameter: ", ll, "-->", max(res[[2]], ll), "\n")
      updateRates(res, ll, rate, shape, k, inv, wMix, update="shape",
                  site.rate=site.rate)
    }
    if (optRate) {
      res <- optimRate(tree, data, rate = rate, inv = inv,
        INV = INV, Q = Q, bf = bf, eig = eig, k = k,
        shape = shape, w = w, ll.0 = ll.0)
      if (trace > 0)
        cat("optimize rate: ", ll, "-->", max(res[[2]], ll), "\n")
      updateRates(res, ll, rate, shape, k, inv, wMix, update="rate",
                  site.rate=site.rate)
    }
    if (optEdge) {
      res <- optimEdge(tree, data, eig = eig, w = w, g = g, bf = bf,
                       rate = rate, ll.0 = ll.0,
                       control = pml.control(epsilon = 1e-08, maxit = 5,
                                             trace = trace))
#      if (trace > 0)
#        cat("optimize edge weights: ", ll, "-->", max(res[[2]], ll),
#          "\n")
      if (res[[2]] > ll) {
        ll <- res[[2]]
        tree <- res[[1]]
      }
    }
    if (optRooted) {
      res <- optimRooted(tree, data, eig = eig, w = w, g = g, bf = bf,
                         rate = rate, ll.0 = ll.0, INV = INV,
                         control = pml.control(epsilon = 1e-07, maxit = 10,
                                               trace = trace - 1))
#      if (trace > 0)
#        cat("optimize edge weights: ", ll, "-->", max(res[[2]], ll), "\n")
      if (res[[2]] > ll) {
        ll <- res[[2]]
        tree <- res[[1]]
      }
    }
    if (optNni) {
      swap <- 0
      iter <- 1
      while (iter < 4) {
        if (optEdge) {
          tmp <- pml.nni(tree, data, w = w, g = g, eig = eig, bf = bf,
            ll.0 = ll.0, ll = ll, INV = INV, ...)
          swap <- swap + tmp$swap
          res <- optimEdge(tmp$tree, data, eig = eig, w = w, g = g, bf = bf,
            rate = rate, ll.0 = ll.0, control = pml.control(
              epsilon = 1e-08, maxit = 3, trace = 0))
          ll2 <- res[[2]]
          tree <- res[[1]]
        }
        else {
          tmp <- rooted.nni(tree, data, eig = eig, w = w, g = g, bf = bf,
            rate = rate, ll.0 = ll.0, INV = INV, ...)
          swap <- swap + tmp$swap
          res <- optimRooted(tmp$tree, data, eig = eig, w = w, g = g, bf = bf,
            rate = rate, ll.0 = ll.0, INV = INV, control = pml.control(
              epsilon = 1e-08, maxit = 5, trace = trace - 1))
          tree <- res$tree
          ll2 <- res$logLik
        }
        if (trace > 0)
          cat("optimize topology: ", ll, "-->", ll2, "\n")
        ll <- ll2
        iter <- iter + 1
        if (tmp$swap == 0) {
          iter <- 4
        }
      }
      if (trace > 0)
        cat(swap, "\n")
      if (swap > 0)
        rounds <- 1
      if (swap == 0)
        optNni <- FALSE
    }
    epsR <- 1e-8
    if ( (ratchet == TRUE) && (optNni == FALSE)) {
      maxR <- ratchet.par$iter
      maxit <- ratchet.par$maxit
      kmax <- 1
      i <- 1
      while (i < maxit) {
        tree2 <- rNNI(tree, moves = round(nTips * ratchet.par$prop), n = 1)
        swap <- 1
        ll2 <- pml.fit(tree2, data, bf, shape = shape, k = k, Q = Q,
          levels = attr(data, "levels"), inv = inv, rate = rate,
          g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0,
          llMix = llMix, wMix = wMix, site = FALSE)

        while (swap > 0) {
          tmp <- pml.nni(tree2, data, w = w, g = g, eig = eig, bf = bf,
            ll.0 = ll.0, ll = ll2, INV = INV, ...)
          swap <- tmp$swap
          res <- optimEdge(tmp$tree, data, eig = eig, w = w, g = g, bf = bf,
            rate = rate, ll.0 = ll.0, control = pml.control(
              epsilon = 1e-08, maxit = 3, trace = 0))
          if (trace > 1)
            cat("optimize topology: ", ll2, "-->", res[[2]], "\n",
              "swaps:", tmp$swap, "\n")
          ll2 <- res[[2]]
          tree2 <- res[[1]]
        }

        if (ll2 > (ll + epsR)) {
          tree <- tree2
          ll <- ll2
          kmax <- 1
        }
        else kmax <- kmax + 1
        if (trace > 0) print(paste("Ratchet iteration ", i,
                                   ", best pscore so far:", ll))
        i <- i + 1
        if (kmax > maxR) i <- maxit
      }
      optNni <- TRUE
      ratchet <- FALSE
      rounds <- 1
    }

    if ((ratchet2 == TRUE) && (optNni == FALSE)) {
      if (optRooted) {
        FUN <- function(x, bf, Q, k, shape) {
          dm <- dist.ml(x, bf = bf, Q = Q, k = k, shape = shape)
          tr <- wpgma(dm)
          tr
        }
      }
      else {
        FUN <- function(x, bf, Q, k, shape) {
          dm <- dist.ml(x, bf = bf, Q = Q)
          tr <- fastme.bal(dm, TRUE, FALSE, FALSE)
          tr$edge.length[tr$edge.length < 1e-8] <- 1e-8
          tr
        }
      }
      maxR <- ratchet.par$iter
      maxit <- ratchet.par$maxit
      kmax <- 1
      i <- 1
      while (i < maxit) {
        tree2 <- bootstrap.phyDat(data, FUN, bs = 1, bf = bf, Q = Q, k = k,
          shape = shape)[[1]]
        tree2 <- checkLabels(tree2, tree$tip.label)
        tree2 <- reorder(tree2, "postorder")
        swap <- 1

        ll2 <- pml.fit(tree2, data, bf, shape = shape, k = k, Q = Q,
          levels = attr(data, "levels"), inv = inv, rate = rate,
          g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0,
          llMix = llMix, wMix = wMix, site = FALSE)
        while (swap > 0) {
          tmp <- pml.nni(tree2, data, w = w, g = g, eig = eig, bf = bf,
                         ll.0 = ll.0, ll = ll2, INV = INV, ...)
          swap <- tmp$swap
          res <- optimEdge(tmp$tree, data, eig = eig, w = w, g = g, bf = bf,
            rate = rate, ll.0 = ll.0, control = pml.control(
              epsilon = 1e-08, maxit = 3, trace = 0))
          if (trace > 1)
            cat("optimize topology: ", ll2, "-->", res[[2]], "\n",
              "swaps:", tmp$swap, "\n")
          ll2 <- res[[2]]
          tree2 <- res[[1]]
        }
        if (ll2 > (ll + epsR)) {
          tree <- tree2
          ll <- ll2
          kmax <- 1
        }
        else kmax <- kmax + 1
        if (trace > 0) print(paste("Ratchet iteration ", i, ", best pscore
                                          so far:", ll))
        i <- i + 1
        if (kmax > maxR) i <- maxit
      }
      optNni <- TRUE
      ratchet2 <- FALSE
      rounds <- 1
    }

    if (rounds > control$maxit) opti <- FALSE
    if ( ( (ll1 - ll) / ll  < control$eps) && rounds > 2) # abs(ll1 - ll)
      opti <- FALSE
    rounds <- rounds + 1
    ll1 <- ll
  }
}


indexNNI3 <- function(tree) {
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  ind <- reorder(tree)$edge[, 2]
  nTips <- length(tree$tip.label)
  ind <- ind[ind > nTips]
  #    ind <- which(child %in% parent)
  Nnode <- tree$Nnode
  #     a         d
  #      \       /
  #       e-----f       c is closest to root, f is root from subtree
  #      /       \
  #     b         c     c(a,b,c,d,e,f)
  edgeMatrix <- matrix(0L, (Nnode - 1), 6)

  pvector <- integer(max(parent))
  pvector[child] <- parent
  tips  <- !logical(max(parent))
  tips[parent] <-  FALSE
  #    cvector <- allCildren(tree)
  cvector <- vector("list", max(parent))
  for (i in seq_along(parent)) cvector[[parent[i]]] <- c(cvector[[parent[i]]],
      child[i])
  k <- 1L
  for (i in ind) {
    f <- pvector[i] # f
    ab <- cvector[[i]] # a,b
    ind1 <- cvector[[f]] # c,d
    cd <- ind1[ind1 != i]
    if (pvector[f]) cd <- c(pvector[f], cd) # cd
    edgeMatrix[k, 1:6] <- c(ab, cd, i, f)
    k <- k + 1L
  }
  edgeMatrix
}

# EL ausserhalb
index2tree <- function(x, tree, root = length(tree$tip.label) + 1L) {
  EL <- numeric(max(tree$edge))
  EL[tree$edge[, 2]] <- tree$edge.length
  pa <- c(5L, 5L, 6L, 6L, 6L)
  ch <- c(1L, 2L, 5L, 4L, 3L)
  elind <- c(1L, 2L, 5L, 4L, 6L)                # raus
  if (x[6L] == root) el <- EL[x[ch]]
  else   el <- EL[x[elind]]
  structure(list(edge = structure(c(x[pa], x[ch]), .Dim = c(5L, 2L)),
    edge.length = el, Nnode = 2L), .Names = c("edge", "edge.length",
    "Nnode"), class = "phylo", order = "postorder")
}


index2tree2 <- function(x, tree, root = length(tree$tip.label) + 1L) {
  EL <- numeric(max(tree$edge))
  EL[tree$edge[, 2]] <- tree$edge.length
  pa <- c(6L, 6L, 5L, 5L, 5L)
  ch <- c(3L, 4L, 6L, 1L, 2L)
  elr <- c(3L, 4L, 5L, 1L, 2L)
  eln <- c(6L, 4L, 5L, 1L, 2L)
  if (x[6L] == root) el <- EL[x[elr]]
  else   el <- EL[x[eln]]
  structure(list(edge = structure(c(x[pa], x[ch]), .Dim = c(5L, 2L)),
    edge.length = el, Nnode = 2L),
  .Names = c("edge", "edge.length", "Nnode"), class = "phylo",
  order = "postorder")
}


# weight, nr, nc, contrast, nco (Reihenfolge beibehalten)
# INV raus
# evi, eve, contrast2 ausserhalb definieren
optimQuartet <- function(tree, data, eig, w, g, bf, rate, ll.0, nTips,
                         weight, nr, nc, contrast, nco, llcomp = -Inf,
                         control = pml.control(epsilon = 1e-08, maxit = 5,
                                               trace = 0), ...) {
  #    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") {
  #        tree <- reorder(tree, "postorder")
  #     print("reorder")
  #    }
  el <- tree$edge.length
  tree$edge.length[el < 1e-08] <- 1e-08
  oldtree <- tree
  k <- length(w)
  loglik <- pml.quartet(tree, data, bf = bf, g = g, w = w, eig = eig,
                        ll.0 = ll.0, k = k, nTips = nTips, weight = weight,
                        nr = nr, nc = nc, contrast = contrast, nco = nco)
  start.ll <- old.ll <- new.ll <- loglik
  #    contrast <- attr(data, "contrast")
  contrast2 <- contrast %*% eig[[2]]
  evi <- (t(eig[[3]]) * bf)
  #    weight <- attr(data, "weight")
  eps <- 1
  iter <- 0

  child <- tree$edge[, 2]
  parent <- tree$edge[, 1]
  m <- max(tree$edge)

  EL <- tree$edge.length
  n <- length(tree$edge.length)

  ind.inv <- which(ll.0 > 0)

  #    nr <- as.integer(length(weight))
  #    nc <- as.integer(length(bf))
  #    nco <- as.integer(nrow(contrast))
  eve <- eig[[2]]
  lg <- k
  ScaleEPS <- 1.0 / 4294967296.0
  #    anc <- Ancestors(tree, 1:m, "parent")
  #    anc0 <- as.integer(c(0L, anc))

  while (eps > control$eps && iter < control$maxit) {
    EL <- .Call("optQrtt", as.integer(parent), as.integer(child), eig, evi,
      EL, w, g, as.integer(nr), as.integer(nc), as.integer(nTips),
      as.double(contrast), as.double(contrast2), nco, data,
      as.double(weight), as.double(ll.0))
    iter <- iter + 1
    #        tree$edge.length <- EL[tree$edge[,2]]
    tree$edge.length <- EL  # [treeP$edge[,2]] # vormals treeP
    newll <- pml.quartet(tree, data, bf = bf, g = g, w = w, eig = eig,
                         ll.0 = ll.0, k = k, nTips = nTips, weight = weight,
                         nr = nr, nc = nc, contrast = contrast, nco = nco)
    eps <- (old.ll - newll) / newll
    if ( (eps < 0) || (newll < llcomp))
      return(list(tree = oldtree, logLik = old.ll, c(eps, iter)))
    oldtree <- tree # vormals treeP
    if (control$trace > 1) cat(old.ll, " -> ", newll, "\n")
    old.ll <- newll
    #        loli <- parent[1]
  }
  if (control$trace > 0) cat(start.ll, " -> ", newll, "\n")
  list(tree = tree, logLik = newll, c(eps, iter)) # vormals treeP
}


# weight, nr, nc, contrast, nco rein
# inv, INV, site, ... raus
# wMix, rate last
pml.quartet <- function(tree, data, bf = rep(.25, 4), k = 1, rate = 1, g, w,
                        eig, ll.0 = NULL, ind.ll0 = NULL, llMix = NULL,
                        wMix = 0, nTips, weight, nr, nc, contrast, nco, ...,
                        site = FALSE) {
  #    k <- as.integer(k)
  #    m = 1

  #    if(any(g<.gEps)){
  #        for(i in seq_along(g)){
  #            if(g[i]<.gEps){
  #                inv <- inv+w[i]
  #            }
  #        }
  #        w <- w[g>.gEps]
  #        g <- g[g>.gEps]
  #        k <- length(w)
  #    }
  #    if (is.null(INV))
  #        INV <- Matrix(lli(data, tree), sparse=TRUE)

  # in C
  if (is.null(ll.0)) {
    ll.0 <- numeric(nr)
  }
  if (is.null(ind.ll0)) {
    ind <- which(ll.0 > 0)
  }
  else ind <- ind.ll0
  #    if(inv>0)
  #        ll.0 <- as.matrix(INV %*% (bf * inv))
  #    if (wMix > 0)
  #        ll.0 <- ll.0 + llMix

  #    node <- tree$edge[, 1]
  #    edge <- tree$edge[, 2]
  #    root <- as.integer(node[length(node)])
  #    el <- as.double(tree$edge.length)
  node <- as.integer(tree$edge[, 1] - nTips - 1L) #    min(node))
  edge <- as.integer(tree$edge[, 2] - 1L)

  #    contrast = attr(data, "contrast")
  #    nco = as.integer(dim(contrast)[1])
  #    dlist=data, nr, nc, weight, k ausserhalb definieren
  #    pmlPart einbeziehen
  siteLik <- .Call("PML4", dlist = data, as.double(tree$edge.length),
    as.double(w), as.double(g), nr, nc, as.integer(k), eig,
    as.double(bf), node, edge, nTips, nco, contrast,
    N = as.integer(length(edge)))
  # in C
  if (!is.null(ll.0)) siteLik[ind] <- log(exp(siteLik[ind]) + ll.0[ind])
  if (wMix > 0) siteLik <- log(exp(siteLik) * (1 - wMix) + llMix)
  loglik <- sum(weight * siteLik)
  return(loglik)
}


index2edge <- function(x, root) {
  ch <- c(1L, 2L, 5L, 4L, 3L)
  elind <- c(1L, 2L, 5L, 4L, 6L)
  if (x[6L] == root) el <- x[ch]
  else   el <- x[elind]
  el
}


pml.nni <- function(tree, data, w, g, eig, bf, ll.0, ll, INV = INV, ...) {
  k <- length(w)
  INDEX <-  indexNNI3(tree)
  tmpl <- pml.fit4(tree, data, bf = bf, g = g, w = w, eig = eig, INV = INV,
                   ll.0 = ll.0, k = k, ...)
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  weight <- as.numeric(attr(data, "weight"))
  contrast <- attr(data, "contrast")
  nco <- as.integer(dim(contrast)[1])
  contrast2 <- contrast %*% eig[[2]]
  evi <- (t(eig[[3]]) * bf)
  nTips <- as.integer(length(tree$tip.label))

  m <- dim(INDEX)[1]
  loglik <- numeric(2 * m)
  edgeMatrix <- matrix(0L, 2 * m, 5)

  anc <- Ancestors(tree, 1:max(tree$edge), "parent")
  loli <- getRoot(tree)

  ind1 <- c(1L, 4L, 3L, 2L, 5L) #
  ind2 <- c(4L, 2L, 3L, 1L, 5L) #

  for (i in 1:m) {
    ei <- INDEX[i, ]
    tree0 <- index2tree(INDEX[i, ], tree, nTips + 1L)
    ch <- ei[5]
    pa <- ei[6]

    # move up
    while (pa != loli) {
      tmpr <- match(loli, INDEX[, 5])
      treetmp <- index2tree(INDEX[tmpr, ], tree, nTips + 1L)
      tmpl <- pml.quartet(treetmp, data, bf = bf, g = g, w = w, eig = eig,
        ll.0 = ll.0, k = k, nTips = nTips, weight = weight,
        nr = nr, nc = nc, contrast = contrast, nco = nco)
      loli <- anc[loli]
    }
    llt0 <- pml.quartet(tree0, data, bf = bf, g = g, w = w, eig = eig,
                        ll.0 = ll.0, k = k, nTips = nTips, weight = weight,
                        nr = nr, nc = nc, contrast = contrast, nco = nco)
    #        new0 <- optimQuartet(tree0, data, eig=eig, w=w, g=g, bf=bf,
    #                rate=rate, ll.0=ll.0, nTips=nTips,
    #                weight=weight, nr=nr, nc=nc, contrast=contrast, nco=nco,
    #                control = list(epsilon = 1e-08, maxit = 3, trace=0))
    tree2 <- tree1 <- tree0
    tree1$edge[, 2] <- tree1$edge[ind1, 2]
    tree1$edge.length <- tree1$edge.length[ind1]
    tree2$edge[, 2] <- tree2$edge[ind2, 2]
    tree2$edge.length <- tree2$edge.length[ind2]

    new1 <- optimQuartet(tree1, data, eig = eig, w = w, g = g, bf = bf,
                         ll.0 = ll.0, nTips = nTips, weight = weight, nr = nr,
                         nc = nc, contrast = contrast, nco = nco,
                         llcomp = ll + 1e-8, ...)
    # new0$logLik+1e-8)
    new2 <- optimQuartet(tree2, data, eig = eig, w = w, g = g, bf = bf,
                         ll.0 = ll.0, nTips = nTips, weight = weight, nr = nr,
                         nc = nc, contrast = contrast, nco = nco,
                         llcomp = ll + 1e-8, ...)
    # new0$logLik+1e-8)
    loglik[(2 * i) - 1] <- new1$logLik
    loglik[(2 * i)] <- new2$logLik
    edgeMatrix[(2 * i) - 1, ] <- new1$tree$edge.length
    edgeMatrix[(2 * i), ] <- new2$tree$edge.length

    # godown or recompute
    if (any (INDEX[i, c(1, 2)] > nTips)) {
      tree00 <- index2tree2(INDEX[i, ], tree, nTips + 1L)
      tmp3 <- pml.quartet(tree00, data, bf = bf, g = g, w = w, eig = eig,
        ll.0 = ll.0, k = k, nTips = nTips, weight = weight,
        nr = nr, nc = nc, contrast = contrast, nco = nco)
      loli <- getRoot(tree00)
    }
    else tmp3 <- pml.quartet(tree0, data, bf = bf, g = g, w = w, eig = eig,
        ll.0 = ll.0, k = k, nTips = nTips, weight = weight,
        nr = nr, nc = nc, contrast = contrast, nco = nco)
  }
  swap <- 0
  eps0 <- 1e-6
  candidates <- loglik > ll + eps0
  #    cat("candidates", sum(candidates), "\n")
  INDEX2 <- t(apply(INDEX, 1, index2edge, root = getRoot(tree)))
  while (any(candidates)) {
    ind <- which.max(loglik)
    loglik[ind] <- -Inf
    if (ind %% 2) swap.edge <- c(2, 4)
    else swap.edge <- c(1, 4)

    IND <- index2edge(INDEX[(ind + 1) %/% 2, ], nTips + 1L)
    treeT <- changeEdge(tree, IND[swap.edge], IND, edgeMatrix[ind, ])
    test <- pml.fit4(treeT, data, bf = bf, k = k, g = g, w = w, eig = eig,
      ll.0 = ll.0, INV = INV, ...)

    if (test <= ll + eps0) candidates[ind] <- FALSE
    if (test > ll + eps0) {
      ll <- test
      swap <- swap + 1
      tree <- treeT
      indi <- which(rep(colSums(apply(INDEX, 1, match, INDEX[(ind + 1) %/% 2, ],
        nomatch = 0)) > 0, each = 2))
      candidates[indi] <- FALSE
      loglik[indi] <- -Inf
    }
  }
  #    trees <- vector("list", 2*m)
  #    for(i in seq_along(loglik)){
  #        ind = i
  #        if( ind %% 2 ) swap.edge = c(2,4)
  #        else swap.edge = c(1,4)
  #        IND = index2edge(INDEX[(ind+1)%/%2,])
  #        tree2 <- try(changeEdge(tree, IND[swap.edge], IND,
  #                     edgeMatrix[ind,]) )
  #        trees[[i]] <- tree2
  #    }
  #    class(trees) <- "multiPhylo"
  #    trees <- .compressTipLabel(trees)
  # , all=loglik
  list(tree = tree, loglik = ll, swap = swap, candidates = candidates)
}
