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
#' M1a can be used as the null model to test positive selection.
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
#' @param ... further arguments passed to or from other methods.
#' @return A data.frame containing the log-likelihood, number of estimated
#' parameters, AIC, AICc and BIC all tested models.  The data.frame has an
#' attributes "env" which is an environment which contains all the trees, the
#' data and the calls to allow get the estimated models, e.g. as a starting
#' point for further analysis (see example).
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
#' @examples
#'
#' \dontrun{
#' dat_dna <- read.phyDat("/home/klaus/Dropbox/Klaus/R/Paml_files/lysin.nuc")
#' dat_codon <- dna2codon(dat_dna)
#' tree <- read.tree("/home/klaus/Dropbox/Klaus/R/Paml_files/lysin.trees")
#' tree <- nnls.phylo(tree, dist.ml(dat_codon))
#' fit <- pml(tree, dat_codon, bf="F3x4")
#' M0 <- optim.pml(fit, model="codon1")
#'
#' M0_M1_M2 <- codonTest(tree, dat_codon)
#'
#'
#' }
#'
#' @keywords cluster
# @param control A list of parameters for controlling the fitting process.


codonTest <- function(tree, object, model=c("M0", "M1a", "M2a"),
                      frequencies="F3x4",
                      opt_freq=FALSE, codonstart = 1, trace=0, ...){
    if(attr(object, "type")=="DNA")
        object <- dna2codon(object, codonstart=codonstart)
    if(is.null(tree$edge.length)) tree <- nnls.phylo(tree, dist.ml(object))
    if( !("M0" %in% model)) model <- c("M0", model)


    if(trace>2) print("optimize model M0")
    fit <- pml(tree, object, bf=frequencies)
    M0 <- optim.pml(fit, model="codon1")

    choices <- c("M0", "M1a", "M2a")
    model <- match.arg(choices, model, TRUE)

#    M1_start <- list(update(M0, dnds=0), update(M0, dnds=1))
#    M1 <- pmlMix(edge ~ ., M1_start, m=2)

    M1a <- NULL
    M2a <- NULL

#    scaler_Q <- getScaler()

    if("M1a" %in% model){
      if(trace>2) print("optimize model M1a")
      M1a_start <- list(update(M0, dnds=0.1), update(M0, dnds=1))
      M1a <- pmlMix(edge ~ M1a, M1a_start, m=2)
    }
    if("M2a" %in% model){
      if(trace>2) print("optimize model M2a")
      M2a_start <- list(update(M0, dnds=0.1), update(M0, dnds=1),
                        update(M0, dnds=3))
      M2a <- pmlMix(edge ~ M2a, M2a_start, m=3)
    }

    result <- list(M0, M1a, M2a)

    M0_ <- c(model="M0", Frequencies=frequencies, "empirical", glance.pml(M0))

    print(M0_)

    class(result) <- c("codonTest") #, "data.frame")
    result
}


#tidy codon
glance.pml <- function(x, ...){
    res <- data.frame(logLik = x$logLik,
                      df = x$df,
                      AIC = AIC(x),
                      BIC = BIC(x))
    res
}


glance.pmlPart <- function(x, ...){
  res <- data.frame(logLik = x$logLik,
#                    df = x$df,
                    AIC = AIC(x),
                    BIC = BIC(x))
  res
}


plot.codonTest <- function(x, model="M1a"){
  return(NULL)
}


