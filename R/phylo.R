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


CodonQ <- function(subs, syn, tstv = 1, dnds = 1) {
  Q <- numeric(length(subs))
  Q[subs == 1] <- 1 # transversion
  Q[subs == 2] <- tstv # transition
  Q[syn == 1] <- Q[syn == 1] * dnds
  Q[syn < 0] <- 0
  Q
}


# needs no Q
optimCodon <- function(tree, data, Q, subs, syn, trace = 0L, ab = c(0, 0),
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


subsChoice <- function(type = .dnamodels) {
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


subsChoice_USER <- function(type = .usermodels, nstates) {
  nstates <- as.integer(nstates)
  n <- (nstates * (nstates-1L)/2)
  subs0 <- numeric(n)
  subsSym <- c(seq_len(n-1), 0)
  subsOrd <- rep(-1, n)
  ind <- cumsum(c(1, (nstates-1):2))
  subsOrd[ind] <- 0
  Q_ord <- rep(0, n)
  Q_ord[ind] <- 1
  type <- match.arg(type)
  switch(type,
         ER = list(optQ = FALSE, optBf = FALSE, subs = subs0),
         FREQ = list(optQ = FALSE, optBf = TRUE, subs = subs0),
         SYM = list(optQ = TRUE, optBf = FALSE, subs = subsSym),
         ORDERED = list(optQ = FALSE, optBf = FALSE, subs = subsOrd,
                        Q=Q_ord),
         GTR = list(optQ = TRUE, optBf = TRUE, subs = subsSym)
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


optimFreeRate <- function(tree, data, g = c(.25, .75, 1, 2), k=4, w=w, ...) {
  g0 <- c(g[1], diff(g))
  g0[g0 < 1e-8] <- 1e-8 # required by constrOptim
  R <- matrix(0, k, k)
  R[lower.tri(R, TRUE)] <- 1
  fn <- function(g0, tree, data, R, w, k, ...) {
    g_new <- as.vector(R %*% g0)
    pml.fit(tree, data, g=g_new, w=w, k=k, ...)
  }
  ui <- rbind(R, diag(k))
  ci <- rep(0, 2 * k)
  # Maybe constrain rates * omega
  res <- constrOptim(g0, fn, grad = NULL, ui = ui, ci = ci, mu = 1e-04,
                     control = list(fnscale = -1), method = "Nelder-Mead",
                     outer.iterations = 100, outer.eps = 1e-08, tree = tree,
                     data = data, R = R, w = w, k=k, ...)
  rate <- res[[1]]
  res[[1]] <- as.vector(R %*% rate)
  res
}


optimWs <- function(tree, data, w = c(.25, .25, .25, .25), g=g, ...) {
  k <- length(w)
  nenner <- 1 / w[1]
  eta <- log(w * nenner)
  eta <- eta[-1]
  fn <- function(eta, tree, data, g, k, ...) {
    eta <- c(0, eta)
    w_new <- exp(eta) / sum(exp(eta))
    pml.fit(tree, data, g=g, w=w_new, k=k, ...)
  }
  if (k == 2) res <- optimize(f = fn, interval = c(-3, 3), lower = -3,
                              upper = 3, maximum = TRUE,
                              tol = .Machine$double.eps^0.25, tree = tree,
                              data = data, g = g, k=k, ...)
  else res <- optim(eta, fn = fn, method = "L-BFGS-B", lower = -5, upper = 5,
                    control = list(fnscale = -1, maxit = 25), gr = NULL,
                    tree = tree, data = data, g = g, k = k, ...)
  w <- exp(c(0, res[[1]]))
  w <- w / sum(w)
  result <- list(par = w, value = res[[2]])
  result
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


# ML F3x4 model
optimF3x4 <- function(tree, data, bf_codon = matrix(.25, 4, 3), trace = 0, ...){
  l <- nrow(bf_codon)
  nenner <- 1 / bf_codon[l, ]
  lbf <- log(bf_codon * rep(nenner, each = 4))
  lbf <- lbf[-l, ]
  codon_abc <- attr(data, "levels")
  fn <- function(lbf, tree, data, ...) {
    dim(lbf) <- c(3, 3)
    bf_codon <- rbind(exp(lbf), c(1, 1, 1))
    bf_codon <- bf_codon / rep(colSums(bf_codon), each = 4)
    bf <- F3x4_freq(bf_codon, CodonAlphabet = codon_abc)
    pml.fit(tree, data, bf = bf, ...)
  }
  res <- optim(par = lbf, fn = fn, gr = NULL, method = "Nelder-Mead", control =
    list(fnscale = -1, maxit = 500, trace = trace), tree = tree,
  data = data, ...)
  bf_codon <- rbind(exp(res[[1]]), c(1, 1, 1))
  bf_codon <- bf_codon / rep(colSums(bf_codon), each = 4)
  bf <- F3x4_freq(bf_codon, CodonAlphabet = codon_abc)
  result <- list(bf = bf, loglik = res[[2]], bf_codon = bf_codon)
  result
}


# predict.pml <- function(object, newdata,...) sum(object$siteLik * newdata)


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
                          maxit = 5, trace = 0, tau=1e-8), llcomp = -Inf) {
  eps <- 1
  iter <- 0
  evi <- (t(eig[[3]]) * bf)
  tau <- control$tau
  while (eps > control$eps && iter < control$maxit) {
    tmp <- fn.quartet(old.el = old.el, eig = eig, bf = bf, dat = dat,
      g = g, w = w, weight = weight, ll.0 = ll.0)
    old.ll <- tmp$ll
    el1 <- fs(old.el[1], eig, tmp$res[, 1], dat[, 1], weight, g = g, w = w,
              bf = bf, ll.0 = ll.0, evi, tau = tau, getA = TRUE, getB = FALSE)
    el2 <- fs(old.el[2], eig, el1[[2]], dat[, 2], weight, g = g, w = w,
              bf = bf, ll.0 = ll.0, evi, tau = tau, getA = TRUE, getB = FALSE)
    el5 <- fs(old.el[5], eig, el2[[2]], tmp$res[, 2], weight, g = g, w = w,
              bf = bf, ll.0 = ll.0, evi, tau = tau, getA = FALSE, getB = TRUE)
    el3 <- fs(old.el[3], eig, el5[[3]], dat[, 3], weight, g = g, w = w,
              bf = bf, ll.0 = ll.0, evi, tau = tau, getA = TRUE, getB = FALSE)
    el4 <- fs(old.el[4], eig, el3[[2]], dat[, 4], weight, g = g, w = w,
              bf = bf, ll.0 = ll.0, evi, tau = tau, getA = FALSE, getB = FALSE)
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



fs <- function(old.el, eig, parent.dat, child.dat, weight, g = g, w = w,
               bf = bf, ll.0 = ll.0, evi, tau=1e-8, getA = TRUE, getB = TRUE) {
  if (old.el < 1e-8) old.el <- 1e-8
  lg <- length(parent.dat)
  P <- getP(old.el, eig, g)
  nr <- as.integer(length(weight))
  nc <- as.integer(length(bf))
  dad <- .Call("getDAD", parent.dat, child.dat, P, nr, nc)
  X <- .Call("getPrep", dad, child.dat, eig[[2]], evi, nr, nc)
  .Call("FS4", eig, as.integer(length(bf)), as.double(old.el), as.double(w),
    as.double(g), unlist(X), child.dat, dad, as.integer(length(w)), nr,
    as.double(weight), as.double(ll.0), as.double(tau),
    as.integer(getA), as.integer(getB))
}


optimEdge <- function(tree, data, eig = eig, w = w, g = g, bf = bf, rate = rate,
                      ll.0 = ll.0, control = pml.control(epsilon = 1e-08,
                        maxit = 10, trace = 0, tau=1e-8), ...) {
  tree <- reorder(tree, "postorder")
  nTips <- length(tree$tip.label)
  el <- tree$edge.length
  tree$edge.length[el < 1e-08] <- 1e-08
  oldtree <- tree
  k <- length(w)
  data <- subset(data, tree$tip.label)
  loglik <- pml.fit4(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, ...)
  start.ll <- old.ll <- loglik
  contrast <- attr(data, "contrast")
  contrast2 <- contrast %*% eig[[2]]
  evi <- (t(eig[[3]]) * bf)
  weight <- attr(data, "weight")
  eps <- 1
  iter <- 0
  treeP <- tree
  tree <- reorder(tree)
  tau <- control$tau
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
  lg <- k
  ScaleEPS <- 1.0 / 4294967296.0
  anc <- Ancestors(tree, 1:m, "parent")
  anc0 <- as.integer(c(0L, anc))

  while (eps > control$eps && iter < control$maxit) {
    EL <- .Call("optE", as.integer(parent), as.integer(child),
      as.integer(anc0), eig, evi, EL, w, g, as.integer(nr),
      as.integer(nc), as.integer(nTips), as.double(contrast),
      as.double(contrast2), nco, data, as.double(weight),
      as.double(ll.0), as.double(tau))
    iter <- iter + 1
    treeP$edge.length <- EL[treeP$edge[, 2]]
    newll <- pml.fit4(treeP, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, ...)
    eps <- (old.ll - newll) / newll
    if (eps < 0) return(list(tree=oldtree, logLik=old.ll))
    oldtree <- treeP
    if (control$trace > 1) cat(old.ll, " -> ", newll, "\n")
    old.ll <- newll
  }
  if (control$trace > 0)
    cat("optimize edge weights: ", start.ll, "-->", newll, "\n")
  list(tree = treeP, logLik = newll, c(eps, iter))
}


pml.move <- function(EDGE, el, data, g=1, w=1, eig=edQt(), k=1, nTips=NULL,
                     bf=length(levels)) {
  node <- EDGE[, 1]
  edge <- EDGE[, 2]
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  node <- as.integer(node - nTips - 1L)
  edge <- as.integer(edge - 1L)
  contrast <- attr(data, "contrast")
  nco <- as.integer(dim(contrast)[1])
  tmp2 <- .Call("PML4", dlist = data, as.double(el), as.double(w), as.double(g),
                nr, nc, k, eig, as.double(bf), node, edge, nTips, nco, contrast,
                N = as.integer(length(edge)))
  return(NULL)
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


guess_model <- function(x){
  model <- x$model
  type <- attr(x$data, "type")
  if(is.null(model)){
    bf <- sd(x$bf)>0
    Q <- FALSE
    if(length(Q)>1) Q <- sd(x$Q)>0
    if(type=="DNA"){
      if(bf==FALSE && Q==FALSE) model <- "JC"
      if(bf==TRUE && Q==FALSE) model <- "F81"
      if(bf==TRUE && Q==TRUE) model <- "GTR"
    }
    else if(bf==FALSE && Q==FALSE){ model <- "Mk"}
    else {model <- "UNKNOWN"}
  }
  if(x$k > 1) {
    if(x$site.rate=="free_rate") model <- paste0(model, "+R(", x$k, ")")
    else model <- paste0(model, "+G(", x$k, ")")
  }
  if(x$inv>0) model <- paste0(model, "+I")
  model
}



optEdgeMulti <- function(object, control = pml.control(epsilon = 1e-8,
                           maxit = 10, trace = 1, tau = 1e-8), ...) {
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
  pmla <- c("tree", "data", "bf", "Q", "inv", "k", "shape", "rate", "model",
            "wMix", "llMix", "dnds", "tstv", "scaleQ", "...")
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
  ASC <- object$ASC
  site.rate <- object$site.rate
  type <- attr(object$data, "type")
  if(type=="CODON"){
    bf_choice <- object$frequencies
  }
  if (is.na(existing[1])) tree <- object$tree
  else tree <- eval(extras[[existing[1]]], parent.frame())
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
    else model <- object$model
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
    model <- object$model
  }
  if (type == "DNA") {
    if (!is.na(existing[9])) {
      model <- match.arg(eval(extras[[existing[9]]], parent.frame()),
                         .dnamodels)
    }
    else model <- object$model
  }
  if (type == "USER") {
    if (!is.na(existing[9])) {
      model <- match.arg(eval(extras[[existing[9]]], parent.frame()),
                         .usermodels)
      if(model=="ORDERED") Q <- subsChoice_USER(model, nc)$Q
      updateEig <- TRUE
    }
    else model <- object$model
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
#    model <- object$model
  }
  rw <- rates_n_weights(shape, k, site.rate)
  g <- rw[, 1]
  w <- rw[, 2]

  if (inv > 0){
    w <- (1 - inv) * w
    g <- g / (1 - inv)
  }
#  if (wMix > 0) w <- (1 - wMix) * w
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

#  resll <- matrix(0, nr, k)
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
    llMix = llMix, ASC=ASC, site.rate=site.rate)
  if (type == "CODON") {
    result$dnds <- dnds
    result$tstv <- tstv
    result$frequencies <- bf_choice
  }
  class(result) <- "pml"
  result
}


### for edge length optimization
pml.fit4 <- function(tree, data, bf = rep(1 / length(levels), length(levels)),
                     inv = 0, g=1, w=1, eig=edQt(), ll.0=NULL,
                     llMix = NULL, wMix = 0, ..., site = FALSE, ASC = FALSE,
                     site.rate = "gamma") {
  weight <- as.double(attr(data, "weight"))
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  nTips <- as.integer(length(tree$tip.label))
  k <- length(w)
  m <- 1

#  if (is.null(ll.0)) {
#    ll.0 <- numeric(attr(data, "nr"))
#  }
#  if (is.null(INV))
#    INV <- Matrix(lli(data, tree), sparse = TRUE)
#  if (inv > 0)
#    ll.0 <- as.matrix(INV %*% (bf * inv))
#  if (ASC)
#    ll.0 <- as.matrix(INV %*% bf)
#  if (wMix > 0)
#    ll.0 <- ll.0 + llMix
#  node <- tree$edge[, 1]
#  edge <- tree$edge[, 2]
  el <- as.double(tree$edge.length)
  node <- as.integer(tree$edge[, 1] - nTips - 1L) #    min(node))
  edge <- as.integer(tree$edge[, 2] - 1L)

  contrast <- attr(data, "contrast")
  nco <- as.integer(dim(contrast)[1])

  siteLik <- .Call("PML4", dlist = data, el, as.double(w), as.double(g), nr, nc,
    k, eig, as.double(bf), node, edge, nTips, nco, contrast,
    N = as.integer(length(edge)))
  if (inv > 0){
    ind <- which(ll.0 > 0)
    siteLik[ind] <- log(exp(siteLik[ind]) + ll.0[ind])
  }
#  if (!is.null(ll.0)) siteLik[ind] <- log(exp(siteLik[ind]) + ll.0[ind])
#  if (wMix > 0) siteLik <- log(exp(siteLik) * (1 - wMix) + llMix)
  resll <- exp(siteLik)
  if (wMix > 0) siteLik <- log(exp(siteLik) * wMix + llMix)
  loglik <- sum(weight * siteLik)
  if (ASC) {
    ind <- seq_len(nc)
    p0 <- sum(exp(siteLik[ind]))
    if(is.nan(log(1 - p0))) browser()
    loglik <- loglik - sum(weight) * log(1 - p0)
  }
  if (!site) return(loglik)
  return(list(loglik = loglik, siteLik = siteLik, resll=resll))
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
#' @param ASC ascertainment bias correction (ASC), allows to estimate models
#' like Lewis' Mkv.
#' @param site.rate Indicates what type of gamma distribution to use. Options
#' are "gamma" approach of Yang 1994 (default), "gamma_quadrature" after the
#' Laguerre quadrature approach of Felsenstein 2001 and "freerate".
## or "lognormal" after a lognormal quadrature approach.
#' @return \code{pml.fit} returns the log-likelihood.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}, \link{pml_bb}, \link{pmlPart}, \link{pmlMix}}
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
                    wMix = 0, ..., site = FALSE, ASC = FALSE,
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
#    if (wMix > 0) w <- (1 - wMix) * w
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
  if (inv > 0) ll.0 <- as.matrix(INV %*% (bf * inv))
# Notwendig??
  if (ASC) ll.0 <- as.matrix(INV %*% bf)
#  if (wMix > 0) ll.0 <- ll.0 + llMix

  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  #    root <- as.integer(node[length(node)])
  el <- as.double(tree$edge.length)
  node <- as.integer(node - nTips - 1L) #    min(node))
  edge <- as.integer(edge - 1L)

  contrast <- attr(data, "contrast")
  nco <- as.integer(dim(contrast)[1])

  resll <- .Call("PML0", dlist = data, el, as.double(g), nr, nc, k, eig,
    as.double(bf), node, edge, nTips, nco, contrast,
    N = as.integer(length(edge)))
  sca <- .Call("rowMax", resll, length(weight), as.integer(k)) + 1
  # nr statt length(weight)
  lll <- resll - sca
  lll <- exp(lll)
  lll <- as.vector(lll %*% w)
  if (inv > 0){
    ind <- which(ll.0 > 0) # automatic in INV gespeichert
    lll[ind] <- lll[ind] + exp(log(ll.0[ind]) - sca[ind])
  }
  siteLik <- log(lll) + sca
  resll2 <- exp(siteLik)
  # needs to change
#  if (wMix > 0) siteLik <- log(exp(siteLik) * (1 - wMix) + llMix)
  if (wMix > 0) siteLik <- log(exp(siteLik) * wMix + llMix)
  loglik <- sum(weight * siteLik)
  if (ASC) {
    ind <- seq_len(nc)
    p0 <- sum(exp(log(lll[ind]) + sca[ind]))
    loglik <- loglik - sum(weight) * log(1 - p0)
  }
  if (!site) return(loglik)
  resll <- exp(resll)
  return(list(loglik=loglik, siteLik=siteLik, resll=resll2, resll2=resll))
}

### @param optF3x4 Logical value indicating if codon frequencies are estimated
### for the F3x4 model


#' Likelihood of a tree.
#'
#' \code{pml} computes the likelihood of a phylogenetic tree given a sequence
#' alignment and a model. \code{optim.pml} optimizes the different model
#' parameters. For a more user-friendly interface see \code{\link{pml_bb}}.
#'
#' Base frequencies in \code{pml} can be supplied in different ways.
#' For amino acid they are usually defined through specifying a model, so the
#' argument bf does not need to be specified. Otherwise if \code{bf=NULL},
#' each state is given equal probability. It can be a numeric vector given the
#' frequencies. Last but not least \code{bf} can be string "equal", "empirical"
#' and for codon models additionally "F3x4".
#'
#' The topology search uses a nearest neighbor interchange (NNI) and the
#' implementation is similar to phyML.  The option model in pml is only used
#' for amino acid models.  The option model defines the nucleotide model which
#' is getting optimized, all models which are included in modeltest can be
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
#' "FLU", "Blosum62", "Dayhoff_DCMut" and "JTT_DCMut") and additionally rate
#' matrices and amino acid frequencies can be supplied.
#'
#' It is also possible to estimate codon models (e.g. YN98), for details see
#' also the chapter in vignette("phangorn-specials").
#'
#' If the option 'optRooted' is set to TRUE than the edge lengths of rooted
#' tree are optimized. The tree has to be rooted and by now ultrametric!
#' Optimising rooted trees is generally much slower.
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
#' are "gamma" approach of Yang 1994 (default), ""gamma_quadrature"" after the
#' Laguerre quadrature approach of Felsenstein 2001 or "freerate".
## or "lognormal" after a lognormal
#' @param object An object of class \code{pml}.
#' @param optNni Logical value indicating whether topology gets optimized (NNI).
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
#' @param x So far only an object of class \code{modelTest}.
#' @param ASC ascertainment bias correction (ASC), allows to estimate models
#' like Lewis' Mkv.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{pml} or \code{optim.pml} return a list of class \code{pml},
#' some are useful for further computations like \item{tree}{the phylogenetic
#' tree.} \item{data}{the alignment.} \item{logLik}{Log-likelihood of the
#' tree.} \item{siteLik}{Site log-likelihoods.} \item{weight}{Weight of the
#' site patterns.}
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml_bb}}, \code{\link{bootstrap.pml}},
#' \code{\link{modelTest}}, \code{\link{pmlPart}}, \code{\link{pmlMix}},
#' \code{\link{plot.phylo}}, \code{\link{SH.test}}, \code{\link{ancestral.pml}}
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
                rate = 1, model = NULL, site.rate = "gamma", ASC = FALSE, ...) {
  call <- match.call()
  extras <- match.call(expand.dots = FALSE)$...
  pmla <- c("wMix", "llMix", "dnds", "tstv")
  existing <- match(pmla, names(extras))
  wMix <- ifelse(is.na(existing[1]), 0,
    eval(extras[[existing[1]]], parent.frame()))
#  llMix <- ifelse(is.na(existing[2]), 0,
#    eval(extras[[existing[2]]], parent.frame()))
  llMix <- 0
  if(!is.na(existing[2])) llMix <- eval(extras[[existing[2]]], parent.frame())

  # allow
  dnds <- tstv <- 1
  dnds <- ifelse(is.na(existing[3]), 1,
    eval(extras[[existing[3]]], parent.frame()))
  tstv <- ifelse(is.na(existing[4]), 1,
    eval(extras[[existing[4]]], parent.frame()))

  if (!inherits(tree, "phylo")) stop("tree must be of class phylo")
  if (is.null(tree$edge.length)) stop("tree must have edge weights")
  if (any(duplicated(tree$tip.label))) stop("tree must have unique labels!")
  nTips <- as.integer(length(tree$tip.label))
  tree <- reorder(tree, "postorder")
  if (any(tree$edge.length < 0)) tree <- minEdge(tree)
  if (!inherits(data, "phyDat")) stop("data must be of class phyDat")
  if (any(is.na(match(tree$tip.label, attr(data, "names")))))
    stop("tip labels are not in data")
  data <- subset(data, tree$tip.label) # needed
  levels <- attr(data, "levels")
  type <- attr(data, "type")
  if (ASC) {
    data <- cbind(constSitePattern(length(tree$tip.label),
                  names=tree$tip.label, levels=levels, type=type), data)
                ##  compress=FALSE)
    attr(data, "weight")[1:attr(data, "nc")] <- 0.0
  }
  weight <- attr(data, "weight")
  nr <- attr(data, "nr")
  nc <- as.integer(attr(data, "nc"))
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
  if(type=="USER" & !is.null(model)){
    model <- match.arg(model, .usermodels)
    if(model=="ORDERED") Q <- subsChoice_USER("ORDERED", nc)$Q
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
  if(site.rate=="free_rate"){
    w <- rep(1/k, k)
    g <- rep(1, k)
  }
  else{
    rw <- rates_n_weights(shape, k, site.rate)
    g <- rw[, 1]
    w <- rw[, 2]
  }
  if (inv > 0){
    w <- (1 - inv) * w
    g <- g / (1 - inv)
  }
#  if (wMix > 0) w <- (1 - wMix) * w
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
  if(ASC) ll.0 <- as.matrix(INV %*% rep(1, length(bf)))

#  if (wMix > 0) ll.0 <- ll.0 + llMix
  nr <- as.integer(attr(data, "nr"))

  on.exit(.Call("ll_free2"))
  .Call("ll_init2", nr, nTips, nc, as.integer(k))
  tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q,
    levels = attr(data, "levels"), inv = inv, rate = rate, g = g, w = w,
    eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix, wMix = wMix, site = TRUE,
    ASC=ASC)
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
    llMix = llMix, ASC=ASC, site.rate=site.rate)
  if (type == "CODON") {
    result$dnds <- dnds
    result$tstv <- tstv
    result$frequencies <- bf_choice
  }
  class(result) <- "pml"
  result
}

# removed INV=NULL
optimRooted <- function(tree, data, bf, g, w, eig, ll.0,
                        control = pml.control(epsilon = 1e-08, maxit = 25,
                                              trace = 0), ...) {
  nTips <- as.integer(length(tree$tip.label))
  k <- length(w)
  tau <- control$tau / 2

  # optimising rooted triplets
  optRoot0 <- function(t, tree, data, g, w, eig, bf, ll.0, ...) {
    l <- length(tree$edge.length)
    tree$edge.length[1:(l - 1)] <- tree$edge.length[1:(l - 1)] + t
    tree$edge.length[l] <- tree$edge.length[l] - t
    loglik <- pml.fit4(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, ...)
    loglik
  }
  # optim edges leading to the root
  # add tau == t ???
  optRoot2 <- function(t, tree, data, g, w, eig, bf, ll.0, ...) {
    tree$edge.length <- tree$edge.length + t   # c(el1+t, el2-t)
    loglik <- pml.fit4(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, ...)
    loglik
  }
  # scale whole tree
  scaleEdges <- function(tree, data, tau = 1e-8, ...) { #t = 1, trace = 0,
    fn <- function(t, tree, data, ...) {
      tree$edge.length <- tree$edge.length * t
      pml.fit4(tree, data, ...)
    }
    min_scaler <- max(.25, tau / min(tree$edge.length) )
    min_scaler <- min(min_scaler, 1)
    if(min_scaler>1) browser()
    optimize(f = fn, interval = c(min_scaler, 4), tree = tree, data = data, ...,
      maximum = TRUE, tol = .00001)
  }
  # ensure that each edge is at least tau long
  # tips have the same height
  tree <- minEdge(tree, 2*tau)

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

  ll <- pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
  start.ll <- ll
  eps <- 10
  iter <- 1

  EL <- numeric(max(tree$edge))
  EL[tree$edge[, 2]] <- tree$edge.length

  if(is.ultrametric(tree)){
    tmp <- scaleEdges(tree, data, tau=tau, bf=bf, ll.0=ll.0,
      eig=eig, w=w, g=g, ...)
    #    if(control$trace>2)cat("scale", tmp[[2]], "\n")
    if(tmp[[2]]>ll){
      t <- tmp[[1]]
      tree$edge.length <- tree$edge.length * t
    }
  }
  el <- tree$edge.length
  EL[tree$edge[, 2]] <- tree$edge.length
  ll2 <- pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
  tmptree <- tree

  while (eps > control$eps && iter < control$maxit) {
    ll2 <- pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
    loli <- rootNode

    children <- allKids[[rootNode]]
    kidsEl <- EL[children]
    minEl <- min(kidsEl)
    kidsEl <- kidsEl - minEl

    tmptree$edge <- cbind(rootNode, children)
    tmptree$edge.length <- kidsEl
    t <- optimize(f = optRoot2, interval = c(tau, 3), tmptree, data = data,
                  g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, maximum = TRUE, ...)
    optRoot2(t[[1]], tmptree, data = data, g=g, w=w, eig=eig, bf=bf,
             ll.0=ll.0, ...)
    # if(control$trace>2)cat("optRoot", t[[2]], "\n")
    if(t[[2]] > ll2){
      ll2 <- t[[2]]
      EL[children] <- kidsEl + t[[1]]
      tree$edge.length <- EL[tree$edge[, 2]]
    }
    ll2 <- pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
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

      t0 <- optRoot0(0, tmptree, data, g, w, eig, bf, ll.0, ...)

## Additional check
      if(c(-minEl + tau < maxEl - tau)) {
        t <- optimize(f = optRoot0, interval = c(-minEl + tau, maxEl - tau),
        tmptree, data = data, g = g, w = w, eig = eig, bf = bf,
        ll.0 = ll.0, maximum = TRUE, ...)
      }
      # if(control$trace>2) cat("edge", t[[2]], "\n")
      if (!is.nan(t[[2]]) & t[[2]] > ll2) {
        optRoot0(t[[1]], tmptree, data = data, g = g, w = w, eig = eig, bf = bf,
          ll.0 = ll.0, k = k, ...)
        EL[children] <- kidsEl + t[[1]]
        EL[dad] <- maxEl - t[[1]]
        ll2 <- t[[2]]
      }
      else optRoot0(0, tmptree, data, g, w, eig, bf, ll.0, k, ...)
      loli <- dad
    }
    tree$edge.length <- EL[tree$edge[, 2]]
    ll2 <- pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
    eps <- (ll - ll2) / ll2

    if (control$trace > 1) cat("optimRooted: ", ll, " -> ", ll2, "\n")
    ll <- ll2
    iter <- iter + 1
  }
  if (control$trace > 0)
    cat("optimize edge weights: ", start.ll, "-->", ll, "\n")
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


rooted.nni <- function(tree, data, eig, w, g, bf, rate, ll.0, INV, RELL=NULL,
                       control = pml.control(epsilon = 1e-08, maxit = 25,
                                             trace = 0), ...) {
#  ind0 <- which(ll.0 > 0)
  contrast <- attr(data, "contrast")
  tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
  nTips <- as.integer(length(tree$tip.label))
  k <- length(w)
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

  optRootU <- function(t, tree, data, bf, g, w, eig, ll.0, nh, ...) {
    tree$edge.length <- getEL1(t, nh)
    pml.fit4(tree, data, bf = bf, eig = eig, ll.0 = ll.0, w = w, g = g, ...)
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

  optEdgeU <- function(t, tree, data, bf, g, w, eig, ll.0, nh, ...) {
    tree$edge.length <- getEL2(t, nh)
    pml.fit4(tree, data, bf = bf, eig = eig, ll.0 = ll.0, w = w, g = g, ...)
  }

  child <- tree$edge[, 2]
  parent <- tree$edge[, 1]
  ll <- pml.fit4(tree, data, bf = bf, eig = eig, ll.0 = ll.0, w = w, g = g, ...)
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
    ll <- pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
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
            pml.move(tmpEdge, EL[c(pa, tmpKids)], data, g, w, eig, k, nTips, bf)
            # cat("move from pa to dad \n")
            loli <- dad
          }
          else {
            #   cat("move loli up", loli, "dad", dad, "pa", pa, "ch", ch, "\n")
            tmpKids <- cvector[[loli]]
            tmpEdge <- cbind(loli, tmpKids)
            pml.move(tmpEdge, EL[tmpKids], data, g, w, eig, k, nTips, bf)
            loli <- anc[loli]
          }

        }

        if (loli == rootNode && dad != loli) {
          # update all nodes
          pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
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
          ll0 <- pml.fit4(tree1, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
          #      cat("quartet", ll0, ch, dad, "\n")
        }

        if (dad == rootNode) {
          ll0 <- pml.fit4(tree1, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
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
            ll.0 = ll.0, ..., method = "L-BFGS-B",
            lower = 1e-8, upper = 5, control = list(fnscale = -1))
          res2 <- optim(par = c(.1, .1), optRootU, gr = NULL, tree = tree2,
            data = data, nh = nh[X2[, 2]], g = g, w = w, eig = eig, bf = bf,
            ll.0 = ll.0, ..., method = "L-BFGS-B",
            lower = 1e-8, upper = 5, control = list(fnscale = -1))
          res3 <- optim(par = c(.1, .1), optRootU, gr = NULL, tree = tree3,
            data = data,  nh = nh[X3[, 2]], g = g, w = w, eig = eig, bf = bf,
            ll.0 = ll.0, ..., method = "L-BFGS-B",
            lower = 1e-8, upper = 5, control = list(fnscale = -1))
          ind <- which.max(c(res1[[2]], res2[[2]], res3[[2]]))
          if (control$trace > 2) cat("root", c(res1[[2]], res2[[2]],
              res3[[2]]), "\n")
#          optRootU <- function(t, tree, data, bf, g, w, eig, ll.0, nh, ...) {
          if (ind == 1) {
            ll2 <- res1[[2]]
            optRootU(t = res1[[1]], tree = tree1, data = data,
              nh = nh[X1[, 2]], g = g, w = w, eig = eig, bf = bf,
              ll.0 = ll.0, ...)
            tmpEL <- getEL1(res1[[1]], nh[X1[, 2]])
            tree <- changeEdgeLength(tree, X1[, 2], tmpEL)
          }
          if (ind == 2) {
            ll2 <- res2[[2]]
            optRootU(t = res2[[1]], tree = tree2, data = data,
              nh = nh[X2[, 2]], g = g, w = w, eig = eig, bf = bf,
              ll.0 = ll.0, ...)
            tmpEL <- getEL1(res2[[1]], nh[X2[, 2]])
            tree <- changeEdge(tree, X1[c(2, 3), 2])
            tree <- changeEdgeLength(tree, X2[, 2], tmpEL)
          }
          if (ind == 3) {
            ll2 <- res3[[2]]
            optRootU(t = res3[[1]], tree = tree3, data = data,
              nh = nh[X3[, 2]], g = g, w = w, eig = eig, bf = bf,
              ll.0 = ll.0, ...)
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
            ..., method = "L-BFGS-B", lower = 1e-4,
            upper = 1 - 1e-4, control = list(fnscale = -1))
          res2 <- optim(par = tt, optEdgeU, gr = NULL, tree = tree2, data,
            nh = nh[X2[, 2]], g = g, w = w, eig = eig, bf = bf, ll.0 = ll.0,
            ..., method = "L-BFGS-B", lower = 1e-4,
            upper = 1 - 1e-4, control = list(fnscale = -1))
          res3 <- optim(par = tt, optEdgeU, gr = NULL, tree = tree3, data,
            nh = nh[X3[, 2]], g = g, w = w, eig = eig, bf = bf, ll.0 = ll.0,
            ..., method = "L-BFGS-B", lower = 1e-4,
            upper = 1 - 1e-4, control = list(fnscale = -1))

          ind <- which.max(c(res1[[2]], res2[[2]], res3[[2]]))
          if (control$trace > 2) cat("edge", ch, ":", c(res1[[2]], res2[[2]],
              res3[[2]]), "\n")
          ll3 <- max(c(res1[[2]], res2[[2]], res3[[2]]))

          if ( (ll3 - 1e-5 * ll3) < ll2) {
            loli <- rootNode
            ll2 <- pml.fit4(tree, data, bf = bf,  k = k, eig = eig, ll.0 = ll.0,
                            INV = INV, w = w, g = g, ...)
            nh <- nodeHeight(tree)
            EL[tree$edge[, 2]] <- tree$edge.length
            ind <- 0
          }
          else {
            if (ind == 1) {
              ll2 <- res1[[2]]
              optEdgeU(res1[[1]], tree = tree1, data, nh = nh[X1[, 2]], g = g,
                w = w, eig = eig, bf = bf, ll.0 = ll.0, ...)
              tmpEL <- getEL2(res1[[1]], nh[X1[, 2]])
              tmpE <- X1[, 2]
              tmpE[5] <- X1[5, 1]
              tree <- changeEdgeLength(tree, tmpE, tmpEL)
            }
            if (ind == 2) {
              ll2 <- res2[[2]]
              optEdgeU(res2[[1]], tree = tree2, data, nh = nh[X2[, 2]], g = g,
                w = w, eig = eig, bf = bf, ll.0 = ll.0, ...)
              tmpEL <- getEL2(res2[[1]], nh[X2[, 2]])
              tmpE <- X2[, 2]
              tmpE[5] <- X1[5, 1]
              tree <- changeEdge(tree, X1[c(2, 3), 2])
              tree <- changeEdgeLength(tree, tmpE, tmpEL)
            }
            if (ind == 3) {
              ll2 <- res3[[2]]
              optEdgeU(res3[[1]], tree = tree3, data, nh = nh[X3[, 2]], g = g,
                w = w, eig = eig, bf = bf, ll.0 = ll.0, ...)
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
          if(!is.null(RELL)){
            siteLik <- pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w,
                                g=g, site=TRUE, ...)$siteLik
            RELL <- update_rell(RELL, siteLik, tree)
          }
        }
      }

    }
    ll2 <- pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g, ...)
    eps <- (ll - ll2) / ll2
    if (control$trace > 1) cat(ll, " -> ", ll2, "\n")
    if (control$trace > 1) cat("swap:", nchanges)
    ll <- ll2
    iter <- iter + 1
  }
  list(tree = tree, logLik = ll, iter = iter, swap = nchanges, RELL=RELL)
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
  # if (wMix > 0) w <- (1 - wMix) * w
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
                      optRate = FALSE, optRooted = FALSE, #optF3x4 = FALSE,
                      control = pml.control(), model = NULL,
                      rearrangement = ifelse(optNni, "NNI", "none"),
                      subs = NULL, ratchet.par = ratchet.control(), ...) {
  rearrangement <- match.arg(rearrangement,
                        c("none", "NNI", "ratchet", "stochastic", "multi2di"))
  optNni <- ifelse(rearrangement ==  "none", FALSE, TRUE)
  perturbation <- ifelse(rearrangement %in%
                        c("ratchet", "stochastic", "multi2di"), TRUE, FALSE)
  extras <- match.call(expand.dots = FALSE)$...
  pmla <- c("wMix", "llMix")
  wMix <- object$wMix
  llMix <- object$llMix
  ASC <- object$ASC
  site.rate <- object$site.rate
  optFreeRate <- FALSE
  if(site.rate=="free_rate"){
    if(optGamma){
      optFreeRate <- TRUE
      optGamma <- FALSE
    }
  }
  optModel <- FALSE
  if (is.null(model)) model <- object$model
  else optModel <- TRUE
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
  data <- object$data
  addTaxa <- FALSE
  trace <- control$trace
  tau <- control$tau
# mit Zeile 2000 vereinheitlichen
  method <- "unrooted"
  is_ultrametric <- FALSE
  timetree <- FALSE
  if (is.rooted(tree)) {
    if (optRooted == FALSE && optEdge == TRUE) {
      tree <- unroot(tree)
      attr(tree, "order") <- NULL
      tree <- reorder(tree, "postorder")
      warning("I unrooted the tree", call. = FALSE)
    }
    else{
      is_ultrametric <- is.ultrametric(tree, option=2)
      if(!is_ultrametric) {
        timetree <- TRUE
        method <- "tipdated"
        tip.dates <- node.depth.edgelength(tree)[seq_len(Ntip(tree))]
        tip.dates <- tip.dates - min(tip.dates)
      } else {
        method <- "ultrametric"
      }
    }
  }
  if (optNni) {
    if(!timetree){
      mapping <- map_duplicates(data)
      if (!is.null(mapping)) {
        orig.data <- data
        addTaxa <- TRUE
        tree <- drop.tip(tree, mapping[, 1])
        tree <- reorder(tree, "postorder")
      }
    }
    if (!is.binary(tree)) tree <- multi2di(tree)
    optEdge <- TRUE
  }
  if (length(tree$tip.label) < (3 + !optRooted)) {
    optNni <- FALSE
    perturbation <- FALSE
  }
  if (length(tree$tip.label) < (2 + !optRooted)) {
    stop("rooted / unrooted tree needs at least 2 / 3 tips")
  }
  tree <- reorder(tree, "postorder")
  if (any(tree$edge.length < tau)) {
    tree <- minEdge(tree, tau, is_ultrametric)
    object <- update.pml(object, tree = tree)
  }
  if (optEdge & optRate & !timetree) {
    warning("You can't optimize edges and rates at the same time, only edges are optimized!", call. = FALSE)
    optRate <- FALSE
  }
  if (optRooted) {
    optEdge <- TRUE
    if (!is.rooted(tree)) stop("tree must be rooted")
  }
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
    if (bf_choice=="F3x4") {
      bf <- F3x4(data)
      bf_codon <- bf_by_codon(data)
      object <- update.pml(object, bf = bf)
    }
  }
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  if (type == "DNA" & optModel) {
    tmp <- subsChoice(model)
    optQ <- tmp$optQ
    if (!optQ) {
      Q <- rep(1, 6)
      object <- update.pml(object, Q = Q)
    }
    optBf <- tmp$optBf
    if (!optBf){
      bf <- c(0.25, 0.25, 0.25, 0.25)
    } else bf <- baseFreq(data)
    object <- update.pml(object, bf = bf)
    subs <- tmp$subs
  }
  if (type == "USER" & optModel) {
    tmp <- subsChoice_USER(model, nc)
    optQ <- tmp$optQ
    if (!optQ){
      Q <- rep(1, (nc*(nc-1L))/2)
      object <- update.pml(object, Q = Q)
    }
    optBf <- tmp$optBf
    if (!optBf){
      bf <- rep(1 / nc, nc)
      object <- update.pml(object, bf = bf)
    }
    subs <- tmp$subs
    if(model=="ORDERED") {
      Q <- tmp$Q
      object <- update.pml(object, Q = Q)
    }
  }
  Q <- object$Q
  if (is.null(subs) & optQ) subs <- c(seq_len(length(Q) - 1), 0)
  bf <- object$bf
  eig <- object$eig
  inv <- object$inv
  k <- object$k
  if (k == 1 & optGamma) {
    optGamma <- FALSE
    message("only one rate class, ignored optGamma")
  }
  if(ASC==TRUE & optInv==TRUE){
    optInv <- FALSE
    message('cannot estimate invariant sites and Mkv model, ignored optInv')
  }
  shape <- object$shape
  w <- object$w
  g <- object$g
  ll0 <- object$logLik
  INV <- object$INV
  ll.0 <- object$ll.0
  rate <- object$rate
  ll <- ll0
  ll1 <- ll0
  opti <- TRUE
  RELL <- NULL
  if(ratchet.par$rell && perturbation){
    RELL <- init_rell(data, B=ratchet.par$bs)
  }
  nTips <- as.integer(length(tree$tip.label))
  on.exit({
    tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q,
      levels = attr(data, "levels"), inv = inv, rate = rate,
      g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix,
      wMix = wMix, site = TRUE, ASC=ASC)
    if (addTaxa) {
      tree <- add.tips(tree, tips = mapping[, 1], where = mapping[, 2],
        edge.length = rep(0, nrow(mapping)))
      data <- orig.data
      if (!is.null(RELL)){
        bs <- RELL$bs
        for(i in seq_along(bs)){
          bs[[i]] <- add.tips(bs[[i]], tips = mapping[, 1],
                              where = mapping[, 2], edge.length = rep(0, nrow(mapping)))
        }
        RELL$bs <- bs
      }
    }
    df <- ifelse(optRooted, tree$Nnode, length(tree$edge.length))
    dfQ <- ifelse(is.null(subs), length(unique(Q)) - 1, max(subs))
    df <- switch(type,
      DNA = df + (k > 1) + (optInv | (inv > 0)) + length(unique(bf)) - 1 + dfQ,
      AA = df + (k > 1) + (optInv | (inv > 0)) +
           optBf * (length(unique(bf)) - 1),
      CODON = df + (k > 1) + (optInv | (inv > 0)) + freq_df + (dnds != 1) +
           (tstv != 1),
      USER = df + (k > 1) + (optInv | (inv > 0)) + length(unique(bf)) - 1 + dfQ)

    object <- list(logLik = tmp$loglik, inv = inv, k = k, shape = shape,
      Q = Q, bf = bf, rate = rate, siteLik = tmp$siteLik,
      weight = attr(data, "weight"),
      g = g, w = w, eig = eig, data = data, model = model,
      INV = INV, ll.0 = ll.0, tree = tree, lv = tmp$resll,
      call = call, df = df, wMix = wMix, llMix = llMix, ASC=ASC,
      site.rate=site.rate)
    if (type == "CODON") {
      object$dnds <- dnds
      object$tstv <- tstv
      object$frequencies <- bf_choice
    }
    class(object) <- "pml"

    extras <- pairlist(bf = bf, Q = Q, inv = inv, shape = shape, rate = rate,
               model=model)[c(optBf, optQ, optInv, optGamma, optRate, optModel)]
    if (length(extras)) {
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if (any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }
    object$call <- call
    if(!is.null(RELL)){
      bs <- RELL$bs
      class(bs) <- "multiPhylo"
      bs <- .compressTipLabel(bs)
      object$bs <- bs
      spl <- as.splits(bs)
      object$tree <- addConfidences(object$tree, spl)
    }
    pml.free()
    return(object)
  })
  pml.init(data, k)

  if (optEdge) {
    res <- opt_Edge(tree, data, rooted = optRooted, eig = eig, w = w, g = g,
      bf = bf, inv=inv, rate = rate, ll.0 = ll.0, INV = INV,
      llMix = llMix, wMix=wMix, ASC=ASC,
      control = pml.control(epsilon = 1e-07, maxit = 10, trace = trace,
                             tau = tau))
    if (res[[2]] > ll) {
      ll <- res[[2]]
      tree <- res[[1]]
    }
  }
  rounds <- 1
  while (opti) {
    if (optRate) {
      res <- optimRate(tree, data, rate = rate, inv = inv,
                       INV = INV, Q = Q, bf = bf, eig = eig, k = k,
                       shape = shape, w = w, ll.0 = ll.0)
      if (trace > 0)
        cat("optimize rate: ", ll, "-->", max(res[[2]], ll), "\n")
      updateRates(res, ll, rate, shape, k, inv, wMix, update="rate",
                  site.rate=site.rate)
    }
    if (optBf) {
      if (type=="CODON" && bf_choice=="F3x4") res <- optimF3x4(tree, data,
            bf_codon = bf_codon, inv = inv, Q = Q, w = w, g = g, INV = INV,
            rate = rate, k = k, llMix = llMix, ASC=ASC)
      else res <- optimBf(tree, data, bf = bf, inv = inv, Q = Q, w = w, g = g,
        INV = INV, rate = rate, k = k, llMix = llMix, wMix=wMix, ASC=ASC)
      if (trace > 0)
        cat("optimize base frequencies: ", ll, "-->", max(res[[2]], ll), "\n")
      if (res[[2]] > ll) {
        bf <- res[[1]]
        eig <- edQt(Q = Q, bf = bf)
        if (inv > 0) ll.0 <- as.matrix(INV %*% (bf * inv))
        if (wMix > 0) ll.0 <- ll.0 + llMix
        ll <- res[[2]]
      }
    }
    if (optQ) {
      if (type == "CODON") {
        ab <- c(tstv, dnds)
        res <- switch(model,
          YN98 = optimCodon(tree, data, Q = length(.sub), subs = .sub,
            syn = .syn, bf = bf, w = w, g = g, inv = inv,
            INV = INV, ll.0 = ll.0, rate = rate, k = k,
            ab = log(ab), optK = TRUE, optW = TRUE),
          codon1 = optimCodon(tree, data, Q = length(.sub), subs = .sub,
            syn = .syn, bf = bf, w = w, g = g, inv = inv,
            INV = INV, ll.0 = ll.0, rate = rate, k = k,
            ab = log(ab), optK = TRUE, optW = TRUE),
          codon2 = optimCodon(tree, data, Q = length(.sub), subs = .sub,
            syn = .syn, bf = bf, w = w, g = g, inv = inv,
            INV = INV, ll.0 = ll.0, rate = rate, k = k,
            ab = log(ab), optK = FALSE, optW = TRUE),
          codon3 = optimCodon(tree, data, Q = length(.sub), subs = .sub,
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
                      inv = inv, INV = INV, ll.0 = ll.0, rate = rate, k = k,
                      llMix = llMix, wMix=wMix, ASC=ASC)
      }
      Q <- res[[1]]
      eig <- edQt(Q = Q, bf = bf)
      if (trace > 0) cat("optimize rate matrix: ", ll, "-->", res[[2]], "\n")
      ll <- res[[2]]
    }
    ### start sitewise
    if (optInv) {
      res <- optimInv(tree, data, inv = inv, INV = INV, Q = Q,
        bf = bf, eig = eig, k = k, shape = shape, rate = rate,
        llMix = llMix, wMix=wMix)
      if (trace > 0)
        cat("optimize invariant sites: ", ll, "-->", max(res[[2]], ll), "\n")
      updateRates(res, ll, rate, shape, k, inv, wMix, update="inv",
                  site.rate=site.rate)
      ll.0 <- as.matrix(INV %*% (bf * inv))
      if (wMix > 0) ll.0 <- ll.0 + llMix
    }
    if (optGamma) {
      res <- optimGamma(tree, data, shape = shape, k = k, inv = inv, INV = INV,
                        Q = Q, bf = bf, eig = eig, ll.0 = ll.0, rate = rate,
                        llMix = llMix, wMix=wMix, ASC=ASC)
      if (trace > 0)
        cat("optimize shape parameter: ", ll, "-->", max(res[[2]], ll), "\n")
      updateRates(res, ll, rate, shape, k, inv, wMix, update="shape",
                  site.rate=site.rate)
    }
    if (optFreeRate) {
      # bis jetzt w nicht optimiert!
      tmp_ll <- ll
      res <- optimFreeRate(tree, data, g = g, k = k, w = w, inv = inv,
                           INV = INV, bf = bf, eig = eig,
                           ll.0 = ll.0, rate = rate)
      scale <- function(tree, g, w){
        blub <- sum(g * w)
        g <- g / blub
        tree$edge.length <- tree$edge.length * blub
        list(tree=tree, g=g)
      }
      if(res[[2]] > ll){
        tmp_sc <- scale(tree, res[[1]], w)
        g0 <- res[[1]]
        blub <- sum(g0 * w)
        g <- g0 / blub
        tree$edge.length <- tree$edge.length * blub
##        if (trace > 0) cat("optimize free rate parameters: ", ll, "-->",
##                           max(res[[2]], ll), "\n")
        ll <- res[[2]]
      }
      res2 <- optimWs(tree, data, w = w, g=g, inv = inv,
                    INV = INV, bf = bf, eig = eig,
                    ll.0 = ll.0, rate = rate)
      if(res2[[2]] > ll){
        w <- res2[[1]]
        blub <- sum(g * w)
        g <- g / blub
        tree$edge.length <- tree$edge.length * blub
        ll <- res2[[2]]
      }
      if (trace > 0) cat("optimize free rate parameters: ", tmp_ll, "-->",
                         ll, "\n")
    }
    ### end sitewise
    if (optEdge) {
      res <- opt_Edge(tree, data, rooted = optRooted, eig = eig, w = w, g = g,
                       bf = bf, inv=inv, rate = rate, ll.0 = ll.0,
                      llMix = llMix, wMix=wMix, ASC=ASC,
                       control = pml.control(epsilon = 1e-08, maxit = 10,
                                             trace = trace, tau = tau))
      if (res[[2]] > ll) {
        ll <- res[[2]]
        tree <- res[[1]]
      }
    }
    epsR <- 1e-8
    if (optNni) {
      res <- opt_nni(tree, data, rooted=optRooted, iter_max=5, trace=trace,
                     ll=ll, w = w, g = g, eig = eig, bf = bf, inv=inv,
                     ll.0 = ll.0, INV = INV, llMix = llMix, wMix=wMix, ASC=ASC,
                     RELL=NULL,
                     control = list(eps=1e-08, maxit=3, trace=trace-1, tau=tau))
      ll <- res$logLik
      tree <- res$tree
      swap <- res$swap
      rounds <- 1
      if (swap == 0) optNni <- FALSE
    }
    if ( (perturbation == TRUE) && (optNni == FALSE)) {
      maxR <- ratchet.par$iter
      maxit <- ratchet.par$maxit
      minit <- ratchet.par$minit
      kmax <- 1
      i <- 1
      if((rearrangement == "stochastic" || rearrangement == "ratchet") && optRooted){
        dm <- dist.ml(data, bf=bf, Q=Q, exclude = "pairwise")
      }
      for(i in seq_len(maxit)){
        if(rearrangement == "stochastic"){
          tree2 <- di2multi(tree, tol = 10 * tau, tip2root = TRUE)
          if (!is.binary(tree2)) {
            tree2 <- multi2di(tree2)
            if(!optRooted) tree2 <- unroot(tree2)
            tree2 <- minEdge(tree2, tau)
            tree2 <- reorder(tree2, "postorder")
          }
          tree2 <- rNNI(tree2, moves = round(nTips * ratchet.par$prop), n = 1)
          if(optRooted){
             tree2 <- nnls.tree(dm, tree2, method = method,
                                tip.dates=tip.dates)
             tree2 <- minEdge(tree2, 10*tau)
          }
        } else if(rearrangement == "ratchet"){
          tree2 <- bootstrap.phyDat(data, candidate_tree, bs = 1, method=method,
                        eps = tau, bf = bf, Q = Q, k = k, shape = shape,
                        tip.dates=tip.dates)[[1]]
          tree2 <- checkLabels(tree2, tree$tip.label)
          tree2 <- reorder(tree2, "postorder")
        } else if(rearrangement == "multi2di"){
          tree2 <- di2multi(tree, tol=10*tau, tip2root=TRUE)
          if(any(degree(tree2)>4)){
            tree2 <- multi2di(tree2)
            if(!optRooted) tree2 <- unroot(tree2)
            tree2 <- minEdge(tree2, tau)
            tree2 <- reorder(tree2, "postorder")
          }
          #else return(list(tree, ll))
          else{
            i <- maxit
            break()
          }
        }
        res <- opt_Edge(tree2, data, rooted=optRooted, eig=eig, w=w, g=g, bf=bf,
                        inv=inv, rate=rate, ll.0=ll.0, llMix = llMix, wMix=wMix,
                        ASC=ASC,
                        control = pml.control(epsilon = 1e-08, maxit = 10,
                                              trace = trace-1L, tau = tau))
        ll2 <- res[[2]]
        tree2 <- res[[1]]
        swap <- 1
#        ll2 <- pml.fit(tree2, data, bf, shape = shape, k = k, Q = Q,
#          levels = attr(data, "levels"), inv = inv, rate = rate,
#          g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0,
#          llMix = llMix, wMix = wMix, site = FALSE, Mkv=Mkv)
#        browser()
        res <- opt_nni(tree2, data, rooted=optRooted, iter_max=25, trace=trace,
                       ll=ll2, w = w, g = g, eig = eig, bf = bf, inv=inv,
                       rate=rate, ll.0 = ll.0, INV = INV, llMix = llMix,
                       wMix=wMix, ASC=ASC, RELL=RELL,
                       control=list(eps=1e-08, maxit=5, trace=trace-1, tau=tau))
        if (res$logLik > (ll + epsR)) {
          tree <- res$tree
          ll <- res$logLik
          kmax <- 1
        }
        else kmax <- kmax + 1
        if(!is.null(RELL)) RELL <- res$RELL
        if (trace > 0) print(paste("Ratchet iteration ", i,
                                   ", best pscore so far:", ll))
        # i <- i + 1
        if ( (kmax >= maxR) && (i >= minit)) break()
      }
      optNni <- TRUE
      perturbation <- FALSE
      rounds <- 1
    }
    if ( perturbation==FALSE && ((abs((ll1 - ll) / ll)  < control$eps) || rounds > control$maxit))
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
# INV raus, inv rein
# evi, eve, contrast2 ausserhalb definieren
optimQuartet <- function(tree, data, eig, w, g, bf, rate, ll.0, nTips,
                         weight, nr, nc, contrast, nco, inv=0, llcomp = -Inf,
                         control = pml.control(epsilon = 1e-08, maxit = 5,
                                               trace = 0, tau = 1e-8), ...) {
  el <- tree$edge.length
  tree$edge.length[el < 1e-08] <- 1e-08
  oldtree <- tree
  k <- length(w)
  loglik <- pml.quartet(tree, data, bf = bf, g = g, w = w, eig = eig,
                        ll.0 = ll.0, k = k, nTips = nTips, weight = weight,
                        inv = inv, nr = nr, nc = nc, contrast = contrast,
                        nco = nco, ...)
  start.ll <- old.ll <- new.ll <- loglik
  contrast2 <- contrast %*% eig[[2]]
  evi <- (t(eig[[3]]) * bf)
  eps <- 1
  iter <- 0

  child <- tree$edge[, 2]
  parent <- tree$edge[, 1]
  m <- max(tree$edge)

  EL <- tree$edge.length
  n <- length(tree$edge.length)

  ind.inv <- which(ll.0 > 0)
  tau <- control$tau
  lg <- k
  ScaleEPS <- 1.0 / 4294967296.0
  #    anc <- Ancestors(tree, 1:m, "parent")
  #    anc0 <- as.integer(c(0L, anc))

  while (eps > control$eps && iter < control$maxit) {
    EL <- .Call("optQrtt", as.integer(parent), as.integer(child), eig, evi,
      EL, w, g, as.integer(nr), as.integer(nc), as.integer(nTips),
      as.double(contrast), as.double(contrast2), nco, data,
      as.double(weight),  as.double(ll.0), as.double(tau))
    iter <- iter + 1
    tree$edge.length <- EL  # [treeP$edge[,2]]
    newll <- pml.quartet(tree, data, bf = bf, g = g, w = w, eig = eig,
                         ll.0 = ll.0, k = k, nTips = nTips, weight = weight,
                         inv = inv, nr = nr, nc = nc, contrast = contrast,
                         nco = nco)
    eps <- (old.ll - newll) / newll
    if ( (eps < 0) || (newll < llcomp))
      return(list(tree = oldtree, logLik = old.ll, c(eps, iter)))
    oldtree <- tree # vormals treeP
#    if (control$trace > 1) cat(old.ll, " -> ", newll, "\n")
    old.ll <- newll
  }
  if (control$trace > 0) cat(start.ll, " -> ", newll, "\n")
  list(tree = tree, logLik = newll, c(eps, iter))
}


pml.quartet <- function(tree, data, bf = rep(.25, 4), k = 1, rate = 1, g, w,
                        eig, ll.0 = NULL, #ind.ll0 = NULL,
                        inv=0, llMix = NULL,
                        wMix = 0, nTips, weight, nr, nc, contrast, nco, ...,
                        site = FALSE, ASC=FALSE) {
  # raus pos_ll.0
  if (is.null(ll.0)) {
    ll.0 <- numeric(nr)
  }
#  if (is.null(ind.ll0)) {
#    ind <- which(ll.0 > 0)
#  }
#  else ind <- ind.ll0
  node <- as.integer(tree$edge[, 1] - nTips - 1L) #    min(node))
  edge <- as.integer(tree$edge[, 2] - 1L)

  siteLik <- .Call("PML4", dlist = data, as.double(tree$edge.length),
    as.double(w), as.double(g), nr, nc, as.integer(k), eig,
    as.double(bf), node, edge, nTips, nco, contrast,
    N = as.integer(length(edge)))
  # in C 1st line out
  if (inv > 0){
    ind <- which(ll.0 > 0) # define outside
    siteLik[ind] <- log(exp(siteLik[ind]) + ll.0[ind])
  }
#  if (!is.null(ll.0)) siteLik[ind] <- log(exp(siteLik[ind]) + ll.0[ind])
#  if (wMix > 0) siteLik <- log(exp(siteLik) * (1 - wMix) + llMix)
  if (wMix > 0) siteLik <- log(exp(siteLik) * wMix + llMix)
  loglik <- sum(weight * siteLik)
  if (ASC) {
    ind <- seq_len(nc)
    p0 <- sum(exp(siteLik[ind]))
    loglik <- loglik - sum(weight) * log(1 - p0)
  }
  return(loglik)
}


index2edge <- function(x, root) {
  ch <- c(1L, 2L, 5L, 4L, 3L)
  elind <- c(1L, 2L, 5L, 4L, 6L)
  if (x[6L] == root) el <- x[ch]
  else   el <- x[elind]
  el
}


pml.nni <- function(tree, data, w, g, eig, bf, ll.0, ll, inv, wMix, llMix,
                    RELL=NULL, ...) {
  k <- length(w)
  INDEX <-  indexNNI3(tree)
  tmpl <- pml.fit4(tree, data, bf=bf, g=g, w=w, eig=eig, inv=inv,
                   ll.0=ll.0, k=k, wMix=wMix, llMix=llMix, ...)
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
        ll.0 = ll.0, k = k, nTips = nTips, weight = weight, inv = inv,
        nr = nr, nc = nc, contrast = contrast, nco = nco, wMix=wMix,
        llMix=llMix, ...)
      loli <- anc[loli]
    }
    llt0 <- pml.quartet(tree0, data, bf = bf, g = g, w = w, eig = eig,
                        ll.0 = ll.0, k = k, nTips = nTips, weight = weight,
                        inv = inv, nr = nr, nc = nc, contrast = contrast,
                        nco = nco, wMix=wMix, llMix=llMix, ...)
    #        new0 <- optimQuartet(tree0, data, eig=eig, w=w, g=g, bf=bf,
    #                rate=rate, ll.0=ll.0, nTips=nTips, weight=weight,
    #                nr=nr, nc=nc, contrast=contrast, nco=nco, inv=0,
    #                control = list(epsilon = 1e-08, maxit = 3, trace=0))
    tree2 <- tree1 <- tree0
    tree1$edge[, 2] <- tree1$edge[ind1, 2]
    tree1$edge.length <- tree1$edge.length[ind1]
    tree2$edge[, 2] <- tree2$edge[ind2, 2]
    tree2$edge.length <- tree2$edge.length[ind2]

    new1 <- optimQuartet(tree1, data, eig = eig, w = w, g = g, bf = bf,
                         ll.0 = ll.0, nTips = nTips, weight = weight, nr = nr,
                         nc = nc, contrast = contrast, nco = nco, inv=inv,
                         llcomp = ll + 1e-8, wMix=wMix, llMix=llMix, ...)
    # new0$logLik+1e-8)
    new2 <- optimQuartet(tree2, data, eig = eig, w = w, g = g, bf = bf,
                         ll.0 = ll.0, nTips = nTips, weight = weight, nr = nr,
                         nc = nc, contrast = contrast, nco = nco, inv=inv,
                         llcomp = ll + 1e-8, wMix=wMix, llMix=llMix, ...)
    # new0$logLik+1e-8)
    loglik[(2 * i) - 1] <- new1$logLik
    loglik[(2 * i)] <- new2$logLik
    edgeMatrix[(2 * i) - 1, ] <- new1$tree$edge.length
    edgeMatrix[(2 * i), ] <- new2$tree$edge.length

    # godown or recompute
    if (any (INDEX[i, c(1, 2)] > nTips)) {
      tree00 <- index2tree2(INDEX[i, ], tree, nTips + 1L)
      tmp3 <- pml.quartet(tree00, data, bf = bf, g = g, w = w, eig = eig,
        ll.0 = ll.0, k = k, nTips = nTips, weight = weight, inv = inv,
        nr = nr, nc = nc, contrast = contrast, nco = nco, wMix=wMix,
        llMix=llMix, ...)
      loli <- getRoot(tree00)
    }
    else tmp3 <- pml.quartet(tree0, data, bf = bf, g = g, w = w, eig = eig,
        ll.0 = ll.0, k = k, nTips = nTips, weight = weight, inv = inv,
        nr = nr, nc = nc, contrast = contrast, nco = nco, wMix=wMix,
        llMix=llMix, ...)
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
      ll.0 = ll.0, inv = inv, wMix=wMix, llMix=llMix, ...)

    if (test <= ll + eps0) candidates[ind] <- FALSE
    if (test > ll + eps0) {
      ll <- test
      swap <- swap + 1
      tree <- treeT
      indi <- which(rep(colSums(apply(INDEX, 1, match, INDEX[(ind + 1) %/% 2, ],
        nomatch = 0)) > 0, each = 2))
      candidates[indi] <- FALSE
      loglik[indi] <- -Inf
      if(!is.null(RELL)){
        siteLik <- pml.fit4(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w,
                            g=g, site=TRUE, wMix=wMix, llMix=llMix, ...)$siteLik
        RELL <- update_rell(RELL, siteLik, tree)
      }
    }
  }
  list(tree = tree, loglik = ll, swap = swap, candidates = candidates,
       RELL=RELL)
}


opt_nni <- function(tree, data, rooted, iter_max, trace, ll, RELL=NULL, ...){
  swap <- 0
  iter <- 0
  llstart <- ll
  while (iter < iter_max) {
    if (!rooted) {
      tmp <- pml.nni(tree, data, ll=ll, RELL=RELL, ...)
      res <- optimEdge(tmp$tree, data, ...)
    }
    else {
      tmp <- rooted.nni(tree, data, RELL=RELL, ...)
      res <- optimRooted(tmp$tree, data, ...)
    }
    if(!is.null(RELL)) RELL <- tmp$RELL
    ll2 <- res$logLik
    if(length(ll2)==0) browser()
    if(ll2 > (ll + 1e-8))  # epsR
      tree <- res$tree
    else {
      res$logLik <- ll
      res$tree <- tree
      tmp$swap <- 0
      ll2 <- ll
    }
    swap <- swap + tmp$swap
    if (trace > 1) cat("optimize topology: ", ll, "-->", ll2,
                       " NNI moves: ", tmp$swap, "\n")
    ll <- ll2
    iter <- iter + 1
    if (tmp$swap == 0) {
      break()
      iter <- iter_max
    }
  }
  if (trace > 0) cat("optimize topology: ", llstart, "-->", ll2,
                     " NNI moves: ", swap, "\n")
  res$iter <- iter
  res$swap <- swap
  res$RELL <- RELL
  res
}


opt_Edge <- function(tree, data, rooted, ...){
  if(rooted){
    res <- optimRooted(tree, data, ...)
  }
  else{
    res <- optimEdge(tree, data, ...)
  }
  res
}


init_rell <- function(x, B = 100L){
  weight <- as.integer(attr(x, "weight"))
  lw <- attr(x, "nr")
  X <- matrix(NA_integer_, B, lw)
  wvec <- rep( seq_len(lw), weight)
  for (i in 1:B) X[i,] <- tabulate(sample(wvec, replace = TRUE), nbins = lw)
  bs <- vector("list", B)
  logLik <- rep(-Inf, B)
  list(X=X, bs=bs, logLik=logLik)
}


update_rell <- function(obj, siteLik, tree){
  rell_tmp <- obj$X %*% siteLik
  rell_ind <- rell_tmp > obj$logLik
  if(any(rell_ind)){
    obj$logLik[rell_ind] <- rell_tmp[rell_ind]
    obj$bs[rell_ind] <- c(tree)
  }
  obj
}
