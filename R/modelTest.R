aic.weights <- function(aic) {
  diff.aic <- aic - min(aic)
  exp(-0.5 * diff.aic) / sum(exp(-0.5 * diff.aic))
}


#' ModelTest
#'
#' Comparison of different nucleotide or amino acid substitution models
#'
#' \code{modelTest} estimates all the specified models for a given tree and
#' data.  When the mclapply is available, the computations are done in
#' parallel. \code{modelTest} runs each model in one thread.  This is may not
#' work within a GUI interface and will not work under Windows.
#'
#' @aliases modelTest AICc
#' @param object an object of class phyDat or pml
#' @param tree a phylogenetic tree.
#' @param model a vector containing the substitution models to compare with
#' each other or "all" to test all available models
#' @param G logical, TRUE (default) if (discrete) Gamma model should be tested
#' @param I logical, TRUE (default) if invariant sites should be tested
#' @param FREQ logical, FALSE (default) if TRUE amino acid frequencies will be
#' estimated.
#' @param k number of rate classes
#' @param control A list of parameters for controlling the fitting process.
#' @param multicore logical, whether models should estimated in parallel.
#' @param mc.cores The number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. Must be at least one, and
#' parallelization requires at least two cores.
#' @return A data.frame containing the log-likelihood, number of estimated
#' parameters, AIC, AICc and BIC all tested models.  The data.frame has an
#' attributes "env" which is an environment which contains all the trees, the
#' data and the calls to allow get the estimated models, e.g. as a starting
#' point for further analysis (see example).
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}}, \code{\link{anova}}, \code{\link[stats]{AIC}},
#' \code{\link{codonTest}}
#' @references Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection
#' and multimodel inference: a practical information-theoretic approach}. 2nd
#' ed. Springer, New York
#'
#' Posada, D. and Crandall, K.A. (1998) MODELTEST: testing the model of DNA
#' substitution. \emph{Bioinformatics} \bold{14(9)}: 817-818
#'
#' Posada, D. (2008) jModelTest: Phylogenetic Model Averaging. \emph{Molecular
#' Biology and Evolution} \bold{25}: 1253-1256
#'
#' Darriba D., Taboada G.L., Doallo R and Posada D. (2011) ProtTest 3: fast
#' selection of best-fit models of protein evolution. . \emph{Bioinformatics}
#' \bold{27}: 1164-1165
#' @keywords cluster
#' @examples
#'
#' \dontrun{
#' example(NJ)
#' (mT <- modelTest(Laurasiatherian, tree, model = c("JC", "F81", "K80", "HKY",
#'                  "SYM", "GTR")))
#'
#' # extract best model
#' (best_model <- as.pml(mT))
#'
#'
#' data(chloroplast)
#' (mTAA <- modelTest(chloroplast, model=c("JTT", "WAG", "LG")))
#'
#' # test all available amino acid models
#' (mTAA_all <- modelTest(chloroplast, model="all", multicore=TRUE, mc.cores=2))
#' }
#'
#' @export modelTest
modelTest <- function(object, tree = NULL, model = NULL, G = TRUE, I = TRUE,
                      FREQ = FALSE, k = 4,
                      control = pml.control(epsilon = 1e-08, maxit = 10,
                      trace = 1), multicore = FALSE, mc.cores = NULL) {
  if(.Platform$OS.type=="windows") multicore <- FALSE
  if (multicore && is.null(mc.cores)) mc.cores <- detectCores()
  if(inherits(object, "DNAbin") || inherits(object, "AAbin"))
    object <- as.phyDat(object)
  if (inherits(object, "phyDat")) data <- object
  if (inherits(object, "pml")) {
    data <- object$data
    if (is.null(tree)) tree <- object$tree
  }
  if (attr(data, "type") == "DNA") type <- .dnamodels
  if (attr(data, "type") == "AA") type <- .aamodels
  if (attr(data, "type") == "USER") type <- .usermodels
  if ( is.null(model) || (length(model)==1 && model == "all") ) model <- type
  model <- match.arg(model, type, TRUE)

  env <- new.env()
  assign("data", data, envir = env)

  if (is.null(tree)) tree <- candidate_tree(data)
  if(!identical(sort(names(data)), sort(tree$tip.label))){
    stop("Labels in tree and data differ!")
  }
  if (is.null(tree$edge.length)){
      tree <- acctran(tree, data)
      tree$edge.length <- tree$edge.length / sum(attr(data, "weight"))
      tree <- minEdge(tree, tau=1e-8)
  }
  trace <- control$trace
  control$trace <- trace - 1
  fit <- pml(tree, data)
  fit <- optim.pml(fit, control = control)
  l <- length(model)
  if (attr(fit$data, "type") == "DNA") FREQ <- FALSE
  n <- 1L + sum(I + G + (G & I) + FREQ + (FREQ & I) + (FREQ & G) +
    (FREQ & G & I))
  nseq <- sum(attr(data, "weight"))

  get_pars <- function(x, nseq){
    c(x$df, x$logLik, AIC(x), AICc(x), AIC(x, k = log(nseq)))
  }

  fitPar <- function(model, fit, G, I, k, FREQ) {
    m <- 1
    res <- matrix(NA, n, 6)
    res <- as.data.frame(res)
    colnames(res) <- c("Model", "df", "logLik", "AIC", "AICc", "BIC")
    data.frame(c("Model", "df", "logLik", "AIC", "AICc", "BIC"))
    calls <- vector("list", n)
    trees <- vector("list", n)
    fittmp <- optim.pml(fit, model = model, control = control)
    res[m, 1] <- model
    pars <- get_pars(fittmp, nseq)
    res[m, 2:6] <- pars
    if (trace > 0) cat(formatC(res[m,1], width=12),
                       prettyNum(pars[-4], preserve.width="individual"), "\n")
    calls[[m]] <- fittmp$call
    trees[[m]] <- fittmp$tree
    m <- m + 1
    if (I) {
      fitI <- optim.pml(fittmp, model = model, optInv = TRUE,
        control = control)
      res[m, 1] <- paste0(model, "+I")
      pars <- get_pars(fitI, nseq)
      res[m, 2:6] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[-4],
                         preserve.width="individual"), "\n")
      calls[[m]] <- fitI$call
      trees[[m]] <- fitI$tree
      m <- m + 1
    }
    if (G) {
#      if (trace > 0) print(paste0(model, "+G"))
      fitG <- update(fittmp, k = k)
      fitG <- optim.pml(fitG, model = model, optGamma = TRUE,
        control = control)
      res[m, 1] <- paste0(model, "+G(", k, ")")
      pars <- get_pars(fitG, nseq)
      res[m, 2:6] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[-4],
                         preserve.width="individual"), "\n")
      calls[[m]] <- fitG$call
      trees[[m]] <- fitG$tree
      m <- m + 1
    }
    if (G & I) {
      fitGI <- update(fitI, k = k)
      fitGI <- optim.pml(fitGI, model = model, optGamma = TRUE,
        optInv = TRUE, control = control)
      res[m, 1] <- paste0(model, "+G(", k, ")+I")
      pars <- get_pars(fitGI, nseq)
      res[m, 2:6] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[-4],
                         preserve.width="individual"), "\n")
      calls[[m]] <- fitGI$call
      trees[[m]] <- fitGI$tree
      m <- m + 1
    }
    if (FREQ) {
      fitF <- optim.pml(fittmp, model = model, optBf = TRUE,
        control = control)
      res[m, 1] <- paste0(model, "+F")
      pars <- get_pars(fitF, nseq)
      res[m, 2:6] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[-4],
                         preserve.width="individual"), "\n")
      calls[[m]] <- fitF$call
      trees[[m]] <- fitF$tree
      m <- m + 1
    }
    if (FREQ & I) {
      fitIF <- update(fitF, inv = fitI$inv)
      fitIF <- optim.pml(fitIF, model = model, optBf = TRUE, optInv = TRUE,
        control = control)
      res[m, 1] <- paste0(model, "+I+F")
      pars <- get_pars(fitIF, nseq)
      res[m, 2:6] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[-4],
                                 preserve.width="individual"), "\n")
      calls[[m]] <- fitIF$call
      trees[[m]] <- fitIF$tree
      m <- m + 1
    }
    if (FREQ & G) {
      fitGF <- update(fitF, k = k, shape = fitG$shape)
      fitGF <- optim.pml(fitGF, model = model, optBf = TRUE,
        optGamma = TRUE, control = control)
      res[m, 1] <- paste0(model, "+G(", k, ")+F")
      pars <- get_pars(fitGF, nseq)
      res[m, 2:6] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[-4],
                         preserve.width="individual"), "\n")
      calls[[m]] <- fitGF$call
      trees[[m]] <- fitGF$tree
      m <- m + 1
    }
    if (FREQ & G & I) {
      fitGIF <- update(fitIF, k = k)
      fitGIF <- optim.pml(fitGIF, model = model, optBf = TRUE,
        optInv = TRUE, optGamma = TRUE, control = control)
      res[m, 1] <- paste0(model, "+G(", k, ")+I+F")
      pars <- get_pars(fitGIF, nseq)
      res[m, 2:6] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[-4],
                         preserve.width="individual"), "\n")
      calls[[m]] <- fitGIF$call
      trees[[m]] <- fitGIF$tree
      m <- m + 1
    }
    list(res, trees, calls)
  }
  if(trace & !multicore) cat("Model        df  logLik   AIC      BIC\n")
  eval.success <- FALSE
  if (!eval.success & multicore) {
    RES <- mclapply(model, fitPar, fit, G, I, k, FREQ, mc.cores = mc.cores)
    eval.success <- TRUE
  }
  if (!eval.success)
    RES <- lapply(model, fitPar, fit, G, I, k, FREQ)
  RESULT <- matrix(NA, n * l, 8)
  RESULT <- as.data.frame(RESULT)
  colnames(RESULT) <- c("Model", "df", "logLik", "AIC", "AICw", "AICc",
    "AICcw", "BIC")

  for (i in 1:l) RESULT[( (i - 1) * n + 1):(n * i), c(1, 2, 3, 4, 6, 8)] <-
    RES[[i]][[1]]
  RESULT[, 5] <- aic.weights(RESULT[, 4])
  RESULT[, 7] <- aic.weights(RESULT[, 6])
  for (i in 1:l) {
    for (j in 1:n) {
      mo <- RES[[i]][[1]][j, 1]
      tname <- paste0("tree_", mo)
      tmpmod <- RES[[i]][[3]][[j]]
      tmpmod["tree"] <- call(tname)
      if (!is.null(tmpmod[["k"]])) tmpmod["k"] <- k
      if (attr(data, "type") == "AA") tmpmod["model"] <- RES[[i]][[1]][1, 1]
      assign(tname, RES[[i]][[2]][[j]], envir = env)
      assign(mo, tmpmod, envir = env)
    }
  }
  attr(RESULT, "env") <- env
  class(RESULT) <- c("modelTest", "data.frame")
  RESULT
}


#' @importFrom generics tidy
#' @export
generics::tidy


#' @export
tidy.modelTest <- function(x, ...) {
  env <- attr(x, "env")
  l <- nrow(x)
  k <- rep(1L, l)
  shape <- rep(NA_real_, l)
  inv <- rep(0, l)
  for (i in seq_len(l)) {
    tmp <- get(x$Model[i], env)
    if (!is.null(tmp[["k"]])) k[i] <- tmp[["k"]]
    if (!is.null(tmp[["shape"]])) shape[i] <- tmp[["shape"]]
    if (!is.null(tmp[["inv"]])) inv[i] <- tmp[["inv"]]
  }
  data.frame(Model = x$Model, k = k, shape = shape, inv = inv)
}


#' @rdname pml
#' @export
as.pml <- function (x, ...)
{
  if (inherits(x, "pml"))
    return(x)
  UseMethod("as.pml")
}


#' @export
as.pml.modelTest <- function(x, model="BIC", ...){
  model <- match.arg(model, c("AIC", "AICc", "BIC", x$Model))
  if(model %in% c("AIC", "AICc", "BIC")){
    model <- x$Model[which.min(x[, model])]
  }
  env <- attr(x, "env")
  best_model <- get(model, env)
  fit <- eval(best_model, env)
  fit$model <- strsplit(model, "\\+")[[1]][1]
  fit
}


