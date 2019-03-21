#' codonTest
#'
#' Models for detecting positive selection
#'
#' \code{codonTest} allows to test for positive selection similar to programs
#' like PAML (Yang ) or HyPhy (Kosakovsky Pond et al. 2005).
#'
#' There are several options for deriving the codon frequencies.
#' Frequencies can be "equal" (1/61), derived from nucleotide frequencies "F1x4"
#' and "F3x4" or "empirical" codon frequencies. The frequencies taken using
#' the empirical frequencies or estimated via maximum likelihood.
#'
#' So far the M0 model (Goldman and Yang 2002), M1a and M2a are
#' implemented. The M0 model is always computed the other are optional.
#' The convergence may be very slow and sometimes fails.
#'
#' @aliases codonTest
#' @param object an object of class phyDat.
#' @param tree a phylogenetic tree.
#' @param model a vector containing the substitution models to compare with
#' each other or "all" to test all available models.
#' @param frequencies a character string or vector defining how to compute
#' the codon frequencies
#' @param opt_freq optimize frequencies (so far ignored)
#' @param codonstart an integer giving where to start the translation. This
#' should be 1, 2, or 3, but larger values are accepted and have for effect to
#' start the translation further within the sequence.
#' @param control a list of parameters for controlling the fitting process.
#' @param ... further arguments passed to or from other methods.
#' @return A list whith an element called summary containing a data.frame with
#' the log-likelihood, number of estimated parameters, etc. of all tested
#' models. An object called posterior which contains the posterior probability
#' for the rate class for each sites and the estimates of the defined models.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}}, \code{\link{pmlMix}}, \code{\link{modelTest}},
#' \code{\link[stats]{AIC}}
#' @references Ziheng Yang (2014). \emph{Molecular Evolution: A Statistical
#' Approach}. Oxford University Press, Oxford
#'
#' Sergei L. Kosakovsky Pond, Simon D. W. Frost, Spencer V. Muse (2005) HyPhy:
#' hypothesis testing using phylogenies, \emph{Bioinformatics}, \bold{21(5)}:
#' 676--679, https://doi.org/10.1093/bioinformatics/bti079
#'
#' Nielsen, R., and Z. Yang. (1998) Likelihood models for detecting positively
#' selected amino acid sites and applications to the HIV-1 envelope gene.
#' \emph{Genetics}, \bold{148}: 929--936
#'
#' @examples
#'
#' \dontrun{
#' # load woodmouse data from ape
#' data(woodmouse)
#' dat_codon <- dna2codon(as.phyDat(woodmouse))
#' tree <- NJ(dist.ml(dat_codon))
#' # optimise the model the old way
#' fit <- pml(tree, dat_codon, bf="F3x4")
#' M0 <- optim.pml(fit, model="codon1")
#' # Now using the codonTest function
#' fit_codon <- codonTest(tree, dat_codon)
#' fit_codon
#' plot(fit_codon, "M1a")
#' }
#'
#' @keywords cluster
# @param control A list of parameters for controlling the fitting process.
codonTest <- function(tree, object, model = c("M0", "M1a", "M2a"),
                      frequencies = "F3x4", opt_freq=FALSE, codonstart = 1,
                      control=pml.control(maxit = 20), ...){
  if (attr(object, "type") == "DNA")
    object <- dna2codon(object, codonstart = codonstart)
  if (is.null(tree$edge.length)) tree <- nnls.phylo(tree, dist.ml(object))
  if (!("M0" %in% model)) model <- c("M0", model)
  trace <- control$trace
  control$trace <- trace - 1

  if (trace > 2) print("optimize model M0")
  fit <- pml(tree, object, bf = frequencies)
  M0 <- optim.pml(fit, model="codon1", control = control)
  result <- cbind(model = "M0", Frequencies = frequencies,
    estimate = "empirical", glance.pml(M0), dnds_0 = M0$dnds,
    dnds_1 = NA, dnds_2 = NA, p_0 = 1, p_1 = NA, p_2 = NA,
    tstv = M0$tstv)

  choices <- c("M0", "M1a", "M2a")
  model <- match.arg(choices, model, TRUE)

  M1a <- NULL
  M2a <- NULL

  estimates <- vector("list", length(model))
  names(estimates) <- model
  estimates[["M0"]] <- M0
  prob <- list()

  if ("M1a" %in% model) {
    if (trace > 2) print("optimize model M1a")
    M1a_start <- list(update(M0, dnds = 0.1, scaleQ = 1),
      update(M0, dnds = 1, scaleQ = 1))
    M1a <- pmlMix(rate ~ M1a, M1a_start, m = 2, control = control)
    result <- rbind(result, cbind(model = "M1a", Frequencies = frequencies,
      estimate = "empirical", glance.pmlMix(M1a),
      dnds_0 = M1a$dnds[1], dnds_1 = 1, dnds_2 = NA,
      p_0 = M1a$omega[1], p_1 = M1a$omega[2], p_2 = NA,
      tstv = M1a$fits[[1]]$tstv))
    prob[["M1a"]] <- neb(M1a)
    estimates[["M1a"]] <- M1a
  }
  if ("M2a" %in% model) {
    if (trace > 2) print("optimize model M2a")
    M2a_start <- list(update(M0, dnds = 0.1, scaleQ = 1),
      update(M0, dnds = 1, scaleQ = 1),
      update(M0, dnds = 3, scaleQ = 1))
    M2a <- pmlMix(rate ~ M2a, M2a_start, m = 3, control = control)
    result <- rbind(result, cbind(model = "M2a", Frequencies = frequencies,
      estimate = "empirical", glance.pmlMix(M2a),
      dnds_0 = M2a$dnds[1], dnds_1 = 1,
      dnds_2 = M2a$dnds[3],
      p_0 = M2a$omega[1], p_1 = M2a$omega[2],
      p_2 = M2a$omega[3],
      tstv = M2a$fits[[1]]$tstv))
    prob[["M2a"]] <- neb(M2a)
    estimates[["M2a"]] <- M2a
  }

  # attr(result, "estimates") <- estimates
  # attr(result, "prob") <- prob
  res <- list(summary = result, posterior = prob, estimates = estimates)
  class(res) <- "codonTest"
  res
}


# tidy codon
glance.pml <- function(x, ...) {
  res <- data.frame(logLik = x$logLik,
    df = x$df,
    AIC = AIC(x),
    BIC = BIC(x))
  res
}


glance.pmlMix <- function(x, ...) {
  nr <- attr(x$fits[[1]]$data, "nr")
  res <- data.frame(logLik = x$logLik,
    df = attr(logLik(x), "df"),
    AIC = AIC(x),
    BIC = AIC(x, k = log(nr)))
  res
}


#' @export
print.codonTest <- function(x, ...) print(x$summary)


#' @export
plot.codonTest <- function(x, model = "M1a", col = c(2, 5, 6), ...) {
  dat <- t(x$posterior[[model]])
  colnames(dat) <- seq_len(ncol(dat))
  barplot(dat, col = col, space = 0, border = NA,
          legend.text = prettyNum(x$estimates[[model]]$dnds), ...)
}


# compute Naive Empirical Bayes (NEB) probabilities
neb <- function(x) {
  p <- x$omega
  l <- length(p)
  index <- attr(x$fits[[1]]$data, "index")
  res <- matrix(0, max(index), l)
  for (i in seq_len(l)) res[, i] <- p[i] * x$fits[[i]]$lv
  res <- res / rowSums(res)
  res[index, ]
}
