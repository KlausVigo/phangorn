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
#' model="3Di" tests the three models for structural alignments mention in
#' Garg(2025).
#'
#' @aliases modelTest AICc
#' @param object an object of class phyDat or pml.
#' @param tree a phylogenetic tree.
#' @param model a vector containing the substitution models to compare with
#' each other or "all" to test all available models.
#' @param G logical, TRUE (default) if (discrete) Gamma model should be tested.
#' @param I logical, TRUE (default) if invariant sites should be tested.
#' @param FREQ logical, FALSE (default) if TRUE amino acid frequencies will be
#' estimated.
#' @param R logical, TRUE (default) if free rate model should be tested.
#' @param k number of rate classes.
#' @param control A list of parameters for controlling the fitting process.
#' @param multicore Currently not used.
## logical, whether models should estimated in parallel.
#' @param mc.cores Currently not used.
## The number of cores to use, i.e. at most how many child
## processes will be run simultaneously. Must be at least one, and
## parallelization requires at least two cores.
#' @return A data.frame containing the log-likelihood, number of estimated
#' parameters, AIC, AICc and BIC all tested models.  The data.frame has an
#' attributes "env" which is an environment which contains all the trees, the
#' data and the calls to allow get the estimated models, e.g. as a starting
#' point for further analysis (see example).
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}}, \code{\link{anova}}, \code{\link[stats]{AIC}},
#' \code{\link{codonTest}}, \code{\link[future]{plan}}
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
#'
#' Puente-Lelievre, Caroline, et al. (2023) "Tertiary-interaction characters
#' enable fast, model-based structural phylogenetics beyond the twilight zone."
#' \emph{bioRxiv}: 2023-12
#'
#' Garg, S.G., Hochberg, G.K.A. (2025) A General Substitution Matrix for
#' Structural Phylogenetics, \emph{Molecular Biology and Evolution},
#' \bold{42(6)}, https://doi.org/10.1093/molbev/msaf124
#' @keywords cluster
#' @examples
#'
#' \dontrun{
#' data(Laurasiatherian)
#' mT <- modelTest(Laurasiatherian, model = c("JC", "K80", "HKY", "GTR"),
#'                 R=TRUE)
#'
#' # Some exploratory data analysis
#' plot(mT$TL, mT$logLik, xlim=c(3,6.5))
#' text(mT$TL, mT$logLik, labels=mT$Model, pos=4)
#'
#' # extract best model
#' (best_model <- as.pml(mT))
#'
#'
#' data(chloroplast)
#' (mTAA <- modelTest(chloroplast, model=c("JTT", "WAG", "LG")))
#'
#' # test all available amino acid models
#' plan(multisession, workers = 2)
#' (mTAA_all <- modelTest(chloroplast, model="all"))
#' plan(sequential)
#' }
#'
#' @export modelTest
modelTest <- function(object, tree = NULL, model = NULL, G = TRUE, I = TRUE,
                      FREQ = FALSE, R=FALSE, k = 4, control = pml.control(),
                      multicore = FALSE, mc.cores = NULL) {
  if(inherits(object, "DNAbin") || inherits(object, "AAbin"))
    object <- as.phyDat(object)
  if (inherits(object, "phyDat")) data <- object
  if (inherits(object, "pml")) {
    data <- object$data
    if (is.null(tree)) tree <- object$tree
  }
  #  assert_phyDat(data)
  if (attr(data, "type") == "DNA") type <- .dnamodels
  if (attr(data, "type") == "AA") type <- .aa_3Di_models
  if (attr(data, "type") == "USER") type <- .usermodels
  if ( is.null(model) ) model <- "all"
  if ( length(model)==1 ) {
    if(model == "all"){
      model <- type
      if (attr(data, "type") == "AA") type <- .aamodels
    } else if(model=="3Di"){
      model <- ._3Di_models
    } else if(model=="AA_3Di") model <- .aa_3Di_models
  }
  # clean up
  if(attr(data, "type") == "DNA" &&
     any(model %in% c("K2P", "K3P", "TPM1", "JC69", "TN93"))){
    # check for models: TPM1 == K81 == K3P & K80 == K2P
    model[model=="K2P"] <- "K80"
    model[model=="K3P"] <- "K81"
    model[model=="TPM1"] <- "K81"
    model[model=="JC69"] <- "JC"
    model[model=="TN93"] <- "TrN"
    model <- unique(model)
  }

  model <- match.arg(model, type, TRUE)


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
                  (FREQ & G & I) + R + (FREQ & R))

  nseq <- sum(attr(data, "weight"))

  get_pars <- function(x, nseq){
    shape <- NA_real_
    if(x$k > 1 && x$site.rate != "free_rate") shape <- x$shape
    c(x$df, x$logLik, AIC(x), AICc(x), AIC(x, k = log(nseq)), x$k, shape, x$inv,
      sum(x$tree$edge.length))
  }

  clean_call <- function(x, k){
    if (!is.null(x[["k"]])) x["k"] <- k
    x
  }

  fitPar <- function(model, fit, G, I, k, FREQ) {
    m <- 1
    res <- matrix(NA, n, 10)
    res <- as.data.frame(res)
    colnames(res) <- c("Model", "df", "logLik", "AIC", "AICc", "BIC", "k",
                       "shape", "inv", "TL")
    calls <- vector("list", n)
    trees <- vector("list", n)
    fittmp <- optim.pml(fit, model = model, control = control)
    res[m, 1] <- model
    pars <- get_pars(fittmp, nseq)
    res[m, -1] <- pars
    if (trace > 0) cat(formatC(res[m,1], width=12),
                       prettyNum(pars[c(1,2,3,5)], preserve.width="individual"), "\n")
    calls[[m]] <- fittmp$call
    trees[[m]] <- fittmp$tree
    m <- m + 1
    if (I) {
      fitI <- optim.pml(fittmp, model = model, optInv = TRUE,
                        control = control)
      res[m, 1] <- paste0(model, "+I")
      pars <- get_pars(fitI, nseq)
      res[m, -1] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[c(1,2,3,5)],
                                                                preserve.width="individual"), "\n")
      calls[[m]] <- fitI$call
      trees[[m]] <- fitI$tree
      m <- m + 1
    }
    if (G) {
      fitG <- update(fittmp, k = k)
      fitG <- optim.pml(fitG, model = model, optGamma = TRUE,
                        control = control)
      res[m, 1] <- paste0(model, "+G(", k, ")")
      pars <- get_pars(fitG, nseq)
      res[m, -1] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[c(1,2,3,5)],
                                                                preserve.width="individual"), "\n")
      calls[[m]] <- clean_call(fitG$call, k=k)
      trees[[m]] <- fitG$tree
      m <- m + 1
    }
    if (G & I) {
      fitGI <- update(fitI, k = k)
      fitGI <- optim.pml(fitGI, model = model, optGamma = TRUE,
                         optInv = TRUE, control = control)
      res[m, 1] <- paste0(model, "+G(", k, ")+I")
      pars <- get_pars(fitGI, nseq)
      res[m, -1] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[c(1,2,3,5)],
                                                                preserve.width="individual"), "\n")
      calls[[m]] <- clean_call(fitGI$call, k=k)
      trees[[m]] <- fitGI$tree
      m <- m + 1
    }
    if (R) {
      fitR <- update(fittmp, k = k, site.rate="free_rate")
      fitR <- optim.pml(fitR, model = model, optGamma = TRUE,
                        control = control)
      res[m, 1] <- paste0(model, "+R(", k, ")")
      pars <- get_pars(fitR, nseq)
      res[m, -1] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[c(1,2,3,5)],
                                                                preserve.width="individual"), "\n")
      calls[[m]] <- clean_call(fitR$call, k=k)
      trees[[m]] <- fitR$tree
      m <- m + 1
    }
    if (FREQ) {
      fitF <- optim.pml(fittmp, model = model, optBf = TRUE,
                        control = control)
      res[m, 1] <- paste0(model, "+F")
      pars <- get_pars(fitF, nseq)
      res[m, -1] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12), prettyNum(pars[c(1,2,3,5)],
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
      res[m, -1] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12),
                         prettyNum(pars[c(1,2,3,5)],
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
      res[m, -1] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12),
                         prettyNum(pars[c(1,2,3,5)],
                                   preserve.width="individual"), "\n")
      calls[[m]] <- clean_call(fitGF$call, k=k)
      trees[[m]] <- fitGF$tree
      m <- m + 1
    }
    if (FREQ & G & I) {
      fitGIF <- update(fitIF, k = k)
      fitGIF <- optim.pml(fitGIF, model = model, optBf = TRUE,
                          optInv = TRUE, optGamma = TRUE, control = control)
      res[m, 1] <- paste0(model, "+G(", k, ")+I+F")
      pars <- get_pars(fitGIF, nseq)
      res[m, -1] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12),
                         prettyNum(pars[c(1,2,3,5)],
                         preserve.width="individual"), "\n")
      calls[[m]] <- clean_call(fitGIF$call, k=k)
      trees[[m]] <- fitGIF$tree
      m <- m + 1
    }
    if (FREQ & R) {
      fitRF <- update(fittmp, k = k, site.rate="free_rate")
      fitRF <- optim.pml(fitRF, model = model, optGamma = TRUE, optBf = TRUE,
                         control = control)
      res[m, 1] <- paste0(model, "+R(", k, ")+F")
      pars <- get_pars(fitRF, nseq)
      res[m, -1] <- pars
      if (trace > 0) cat(formatC(res[m,1], width=12),
                         prettyNum(pars[c(1,2,3,5)],
                        preserve.width="individual"), "\n")
      calls[[m]] <- clean_call(fitRF$call, k=k)
      trees[[m]] <- fitRF$tree
      m <- m + 1
    }
    list(res, trees, calls)
  }
  if(trace & !multicore) cat("Model        df  logLik   AIC      BIC\n")

  RES <- future_lapply(model, fitPar, fit, G, I, k, FREQ, future.seed = TRUE)
  RESULT <- matrix(NA, n * l, 12)
  RESULT <- as.data.frame(RESULT)
  colnames(RESULT) <- c("Model", "df", "logLik", "AIC", "AICw", "AICc",
                        "AICcw", "BIC", "k", "shape", "inv", "TL")

  trees <- vector("list", n * l)
  calls <- vector("list", n * l)

  for (i in 1:l){
    ind <- ((i - 1) * n + 1):(n * i)
    RESULT[ind, -c(5, 7)] <- RES[[i]][[1]]
    trees[ind] <- RES[[i]][[2]]
    calls[ind] <- RES[[i]][[3]]


  }
  RESULT[, 5] <- aic.weights(RESULT[, 4])
  RESULT[, 7] <- aic.weights(RESULT[, 6])

  names(trees) <- RESULT[, 1]
  names(calls) <- RESULT[, 1]
  class(trees) <- "multiPhylo"
  trees <- .compressTipLabel(trees)

  class(RESULT) <- c("modelTest", "data.frame")
  attr(RESULT, "data") <- data
  attr(RESULT, "trees") <- trees
  attr(RESULT, "calls") <- calls
  if (trace > 0) {
    tmp <- RESULT$Model[which.min(RESULT$BIC)]
    cat("\nBest fitting model according to BIC:", tmp, "\n\n")
  }
  RESULT
}



## @importFrom generics tidy
## @export
#generics::tidy


## @export
#tidy.modelTest <- function(x, ...) {
#  env <- attr(x, "env")
#  l <- nrow(x)
#  k <- rep(1L, l)
#  shape <- rep(NA_real_, l)
#  inv <- rep(0, l)
#  for (i in seq_len(l)) {
#    tmp <- get(x$Model[i], env)
#    if (!is.null(tmp[["k"]])) k[i] <- tmp[["k"]]
#    if (!is.null(tmp[["shape"]])) shape[i] <- tmp[["shape"]]
#    if (!is.null(tmp[["inv"]])) inv[i] <- tmp[["inv"]]
#  }
#  data.frame(Model = x$Model, k = k, shape = shape, inv = inv)
#}


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
  mods <- strsplit(x$Model, "\\+")
  names(mods) <- x$Model

  best_model <- attr(x, "calls")[[model]]
  tree <- attr(x, "trees")[[model]]
  data <- attr(x, "data")
  fit <- eval(best_model)
  fit$model <- strsplit(model, "\\+")[[1]][1]
  fit
}


#' @param x an object of class modelTest.
#' @param digits	default is 10, i.e. edge length for the bootstrap trees are
#' exported. For digits larger smaller than zero no edge length are exported.
#' @param file a file name. File endings are added.
#' @rdname modelTest
#' @export
write.modelTest <- function(x, file="modelTest", digits=10){
  zzfil <- paste0(file, ".nex.gz")
  zz <- gzfile(zzfil, "w")
  write.nexus(attr(x, "trees"), file=zz, digits=digits)
  close(zz)
  write.phyDat(attr(x, "data"), file=paste0(file, ".fasta"), format="fasta")
  write.table(cbind(x$Model, unlist(attr(x, "calls"))), row.names=FALSE,
              sep = " ", file=paste0(file, ".tsv"),
              col.names = c("Model", "Call"))
}



#' @srrstats {G1.0} in the lines folloing: 39
#' @srrstats {G2.3, G2.3a} in lines: 89, 295
