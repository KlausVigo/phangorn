#' codonTest
#'
#' Models for detecting positive selection
#'
#' Currently \code{codonTest} estimates all the specified models for a given
#' tree and data.
#' There are several options for deriving the codon frequencies.
#' Frequencies can be "equal" (1/61), derived from nucleotide frequencies "F1x4"
#' and "F3x4" or "empirical" codon frequencies. The frequencies taken using
#' the empirical frequencies or estimated via maximum likelihood.
#'
#' So far the M0 model (Goldman and Yang 2002), the M1, M1a and M2a are
#' implemented. The M0 model is always computed the other a re optional.
#' M1a can be used as the null model to test positive selection
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
#' @param control A list of parameters for controlling the fitting process.
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
#' @keywords cluster

codonTest <- function(tree, object, model=c("M0", "M1a", "M2a"),
                      frequencies="F3x4",
                      opt_freq=FALSE, codonstart = 1, ...){
    if(attr(object, "type")=="DNA")
        object <- dna2codon(object, codonstart=codonstart)
    fit <- pml(tree, object, bf=frequencies)
    M0 <- optim.pml(fit, model="codon1")

    choices <- c("M0", "M1", "M1a", "M2a")
    model <- match.arg(choices, model, TRUE)

    M1_start <- list(update(M0, dnds=0), update(M0, dnds=1))
    M1 <- pmlMix(edge ~ ., M1_start, m=2)

    M1a_start <- list(update(M0, dnds=0.1), update(M0, dnds=1))
    M1a <- pmlMix(edge ~ M1a, M1a_start, m=2)

    M2a_start <- list(update(M0, dnds=0.1), update(M0, dnds=1),
                      update(M0, dnds=3))
    M2a <- pmlMix(edge ~ M2a, M2a_start, m=3)

    #    attr(RESULT, "env") <- env

    result <- list(M0, M1a, M2a)

    M0_ <- c(model="M0", Frequencies=frequencies, "empirical", glance.pml(M0))

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



