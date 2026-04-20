aic.weights <- function(aic) {
  diff.aic <- aic - min(aic)
  exp(-0.5 * diff.aic) / sum(exp(-0.5 * diff.aic))
}


as_par <- function(models, k=4, inv=0.2, site_rate="gamma",
                   I = TRUE, G = TRUE, FREQ=FALSE){
  # create_model
  l <- length(site_rate)
  if(!is.list(k)) k <- as.list(rep(k, l))
  l_r <- sum(lengths(k))
  l_m <- length(models)
  l_g <- length(site_rate)
  nb_m <- l_m + I * l_m + G * l_r * l_m + (G && I) * l_r * l_m
  res <- as.data.frame(matrix(NA, nb_m, 5, dimnames = list(NULL,
                          c("model_term", "model", "k", "site_rate", "inv"))))
  pos <- 1
  for(m in models){
    m_term <- model_term(model=m, k=1, site.rate=site_rate[1], inv=0)
    res[pos,] <- list(m_term=m_term, model=m, k=1, site_rate=site_rate[1], inv=0)
    pos <- pos + 1
  }
  if(I){
    for(m in models){
      m_term <- model_term(model=m, k=1, site.rate=site_rate[1], inv=inv)
      res[pos,] <- list(m_term=m_term, model=m, k=1, site_rate=site_rate[1], inv=inv)
      pos <- pos + 1
    }
  }
  if(G){
    for(m in models){
      for(j in seq_along(site_rate)){
        for(kj in k[[j]]){
          m_term <- model_term(model=m, k=kj, site.rate=site_rate[j], inv=0)
          res[pos,] <- list(m_term=m_term, model=m, k=kj, site_rate=site_rate[j], inv=0)
          pos <- pos + 1
        }
      }
    }
  }
  if(G && I){
    for(m in models){
      for(j in seq_along(site_rate)){
        for(kj in k[[j]]){
          m_term <- model_term(model=m, k=kj, site.rate=site_rate[j], inv=inv)
          res[pos,] <- list(m_term=m_term, model=m, k=kj, site_rate=site_rate[j], inv=inv)
          pos <- pos + 1
        }
      }
    }
  }
  res$anc <- res$model
  if(FREQ){
    res2 <- res
    res2$model_term <- paste0(res2$model_term, "+F")
    res2$anc <- paste0(res2$anc, "+F")
    res <- rbind(res, res2)
  }
  res
}


fitPar <- function(par, fit, trees=NULL, calls=NULL, ...) {
  if(!is.null(trees) && !is.null(calls)){
    data <- fit$data
    tree <- trees[[par$anc]]
    fit <- eval(calls[[par$anc]])
#    fit$model <- strsplit(par$model, "\\+")[[1]][1]
#     tmp <- update(fit, k = par$k, site.rate = par$site_rate,
#                inv = par$inv)
  }
  tmp <- update(fit, model = par$model, k = par$k, site.rate = par$site_rate,
                inv = par$inv)
  tmp <- pml_bb(tmp, model = par$model_term, rearrangement = "none", ...)
  pars <- glance(tmp)
  list(model_term = par$model_term, model = par$model, call = tmp$call,
       tree = tmp$tree, pars = pars)
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
#' @param RHAS a character vector specifying the rate heterogeneity among sites
#' models. Option are "gamma" for discrete gamma model with equal weights,
#' "gamma_weighted" for discrete gamma model with estimated weights,
#' "free_rate" and "gamma_quadrature".
#' @param k number of rate classes. Can be a list with a vector for each RHAS
#' term.
#' @param mt_control a list with some options.
#' @param control A list of parameters for controlling the fitting process.
#' @param ... Further arguments passed to or from other methods.
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
#' mT <- modelTest(Laurasiatherian, model = c("JC", "K80", "HKY", "GTR"))
#'
#' # Some exploratory data analysis
#' plot(mT$TL, mT$logLik, xlim=c(3,6.5))
#' text(mT$TL, mT$logLik, labels=mT$Model, pos=4)
#'
#' fit_GTR_G <- as.pml(mt, "GTR+G(4)")
#' fit_GTR_GW <- as.pml(mt, "GTR+GW(4)")
#' fit_GTR_R <- as.pml(mt, "GTR+R(4)")
#' plotRates(fit_GTR_G)
#' plotRates(fit_GTR_GW, append=TRUE, col="red")
#' plotRates(fit_GTR_R, append=TRUE, col="green")
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
#' @export
modelTest <- function(object, tree = NULL, model = NULL, G = TRUE, I = TRUE,
                       FREQ = FALSE, k = 4, control = pml.control(),
                       RHAS = "gamma", ...,
                       mt_control=list(crit="BIC", n_model=100, n_rhas=7)) {
  crit <- mt_control$crit

  if(inherits(object, "DNAbin") || inherits(object, "AAbin"))
    object <- as.phyDat(object)
  if (inherits(object, "phyDat")) data <- object
  if (inherits(object, "pml")) {
    data <- object$data
    if (is.null(tree)) tree <- object$tree
  }
  RHAS <- match.arg(RHAS, choices = c("gamma", "gamma_weighted",
                    "gamma_quadtrature", "free_rate"), several.ok = TRUE)
  gld <- glance(data)
  inv0 <- max(0, 0.9 * (gld$const_sites / gld$nchar))
  #inv0 <- 0
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

  if (is.null(tree)){
 #   p("Compute starting tree")
    tree <- candidate_tree(data)

  }
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
  fit <- pml(tree, data, model=model[1], ...)
  fit <- optim.pml(fit, control = control)
  l <- length(model)
  if (attr(fit$data, "type") == "DNA") FREQ <- FALSE
  nseq <- sum(attr(data, "weight"))

  clean_call <- function(x, k){
    if (!is.null(x[["k"]])) x["k"] <- k
    x
  }
  pars <- as_par(model, k=k, inv=inv0, site_rate=RHAS, I = I, G = G,
                 FREQ = FREQ)

  ind <- which(pars$k==1 & pars$inv == 0)
  spl_pars <- split(pars[ind, ], ind)
  names(spl_pars) <- pars$model_term[ind]

  p <- progressor(along = spl_pars)
  RES <- lapply(spl_pars, function(x, fit, control=control, ...){
    p()
    fitPar(x, fit=fit, control=control, ...)
  }, fit, control=control, future.seed = TRUE)

#  RES <- future_lapply(spl_pars, fitPar, fit, control = control,
#                       future.seed = TRUE)


  fun <- function(x){
    trees <- lapply(x, \(x)x$tree)
    class(trees) <- "multiPhylo"
    calls <- sapply(x, \(x)x$call)
    df <- as.data.frame(matrix(NA, length(x), length(x[[1]]$pars),
            dimnames = list(names(x), names(x[[1]]$pars))))
    for(i in seq_along(x)) df[i,] <- x[[i]]$pars
    list(df=df, calls=calls, trees=trees)
  }

  res1 <- fun(RES)
  df <- res1$df
  calls <- res1$calls
  trees <- res1$trees

  pars2 <- pars[-ind, , drop=FALSE]

  RESULT <- df

  if( G || I){
    best_models <- df$Substitution[order(df[, crit])] |> unique()
    #  sort(setNames(df[, crit], rownames(df)))
    #best_models <- gsub(best_models, "+F") |> unique()
    best_model <- best_models[1] #names(best_models)[1]
    best_models <- best_models[-1]#names(best_models)[-1]
    best_models <- best_models[ seq_len(min(length(best_models), mt_control$n_model)) ]
    ind2 <- which(pars2$model==best_model)
    pars3 <- pars2[-ind2,]
    pars2 <- pars2[ind2, ]
    spl_pars2 <- split(pars2, seq_len(nrow(pars2)))
    names(spl_pars2) <- pars2$model_term

#    p <- progressor(along = spl_pars2)
#    RES2 <- future_lapply(spl_pars2, fitPar, fit, trees = trees, calls = calls,
#                          control = control, future.seed = TRUE)

    p <- progressor(along = spl_pars2)
    RES2 <- lapply(spl_pars2, function(x, fit, trees = trees, calls = calls,
                                      control=control,...){
      p()
      fitPar(x, fit=fit, trees = trees, calls = calls, control=control, ...)
    }, fit,  trees = trees, calls = calls, control=control, future.seed = TRUE)
    res2 <- fun(RES2)
    df2 <- res2$df
    best_rhas <- sort(setNames(df2[, crit], rownames(df2)))
    best_rhas <- best_rhas[seq_len(min(length(best_rhas), mt_control$n_rhas))]
    best_rhas <- names(best_rhas)

    rest <- NULL
    for(i in best_models) rest <- c(rest, gsub(best_model, i, best_rhas))

    trees3 <- NULL
    calls3 <- NULL
    ind3 <- match(rest, pars3$model_term)
    if(length(ind3)>0) {
      pars3 <- pars3[ind3, ]
      spl_pars3 <- split(pars3, seq_len(nrow(pars3)))
      names(spl_pars3) <- pars3$model_term
#      p <- progressor(along = spl_pars3)
#      RES3 <- future_lapply(spl_pars3, fitPar, fit, trees=trees, calls=calls,
#                            control=control, future.seed = TRUE)

      p <- progressor(along = spl_pars3)
      RES3 <- lapply(spl_pars3, function(x, fit, trees = trees, calls = calls,
                                        control=control,...){
        p()
        fitPar(x, fit=fit, trees = trees, calls = calls, control=control, ...)
      }, fit,  trees = trees, calls = calls, control=control, future.seed = TRUE)
      res3 <- fun(RES3)
      calls3 <- res3$calls
      tree3 <- res3$trees
    }
    calls <- c(calls, res2$calls, calls3)
    trees <- c(trees, res2$trees, trees3)
    trees <- .compressTipLabel(trees)
    RESULT <- rbind(df, df2, res3$df)
  }

  RES_1 <-  cbind(Model=rownames(RESULT), AICw = aic.weights(RESULT[, "AIC"]),
                  AICcw = aic.weights(RESULT[, "AICc"]), RESULT[, 1:5])
  RES_1 <- RES_1[, c("Model", "Substitution", "df", "logLik", "AIC", "AICw", "AICc",
                     "AICcw")]
  RESULT <- cbind(RES_1, RESULT[, -(1:5)])
  RESULT <- RESULT[order(RESULT[,crit]), ]

#  RESULT <- cbind(Model=rownames(RESULT), AICw = aic.weights(RESULT[, "AIC"]),
#                  AICcw = aic.weights(RESULT[, "AICc"]), RESULT)
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
  best_model <- attr(x, "calls")[[model]]
  tree <- attr(x, "trees")[[model]]
  data <- attr(x, "data")
  fit <- eval(best_model)
  fit$model <- strsplit(model, "\\+")[[1]][1]
  fit
}

# c.modelTest <- function(...){}


#' @param x an object of class modelTest.
#' @param digits	default is 10, i.e. edge length for the bootstrap trees are
#' exported. For digits larger smaller than zero no edge length are exported.
#' @param file a file name. File endings are added.
#' @importFrom progressr progressor
#' @rdname modelTest
#' @export
write.modelTest <- function(x, file="modelTest", digits=10){
  zzfil <- paste0(file, ".nex.gz")
  zz <- gzfile(zzfil, "w")
  write.nexus(attr(x, "trees"), file=zz, digits=digits)
  close(zz)
  write.phyDat(attr(x, "data"), file=paste0(file, ".fasta"), format="fasta")

  calls <- as.character(gsub("\"", "\'", attr(x, "calls")))

  write.table(cbind(x, Call=calls), row.names=FALSE,
              sep = "\t", file=paste0(file, ".tsv"))
}


#summary.modelTest <- function(object, max.print=10, sort#="BIC", plot=TRUE, ...){
#  browser()
#  tmp <- object[order(object[, sort]),
#                c("Model", "df", "logLik", "AIC", "AICc", #"BIC", "TL")]
#  print(tmp, row.names="Model")
#  invisible(object)
#}


#plot.modelTest <- function(object, ...){
#  library(ggrepel)
#  plot(logLik ~ TL, data=object)
#  text(y=object$logLik, x=object$TL, labels=object$Model)
#  if(requireNamespace(ggplot2) && requireNamespace(ggrepel)){
#  p <- ggplot(object, aes(TL, logLik,  label = Model)) +
#    geom_point() + geom_text_repel(cex=2.25) + theme_bw() +
#    xlim(3.1, 6.2) +
#    xlab("total tree length") + ylab("log-Likelihood")
#  return(p)
#}
#  invisible(object)
#}

