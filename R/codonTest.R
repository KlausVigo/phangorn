#' codonTest
#' 
#' Comparison of different codon substitution models
#' 
#' \code{codonTest} estimates all the specified models for a given tree and
#' data. 
#' 
#' 
#' Currently the codon frequencies, options are "equal", "F1x4", "F3x4" 
#' and "empirical" ("F61"). The frequencies can be taken from the empirical 
#' or estimated. 
#' 
#' 
#' @aliases mcodonTest 
#' @param object an object of class phyDat.
#' @param tree a phylogenetic tree.
#' @param model a vector containing the substitution models to compare with
#' each other or "all" to test all available models. 
#' @param frequencies a character string or vector defining how to compute 
#' the codon frequencies 
#' @param k number of rate classes
#' @param control A list of parameters for controlling the fitting process.
#' @return A data.frame containing the log-likelihood, number of estimated
#' parameters, AIC, AICc and BIC all tested models.  The data.frame has an
#' attributes "env" which is an environment which contains all the trees, the
#' data and the calls to allow get the estimated models, e.g. as a starting
#' point for further analysis (see example).
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}}, \code{\link{pmlMix}}, \code{\link{modelTest}},
#' \code{\link[stats]{AIC}}
#' @keywords cluster

codonTest <- function(tree, object, model=c("M0", "M1a", "M2a"), 
                      frequencies="F3x4", 
                      empirical=TRUE, codonstart = 1, ...){
    if(attr(object, "type")=="DNA") 
        object <- dna2codon(object, codonstart=codonstart)
    fit <- pml(tree, object, bf=frequencies)
    M0 <- optim.pml(fit, model="codon1")
    
#    M1_start <- list(update(M0, dnds=0), update(M0, dnds=1))
#    M1 <- pmlMix(edge ~ ., M1_start, m=2)
    
    M1a_start <- list(update(M0, dnds=0.1), update(M0, dnds=1))
    M1a <- pmlMix(edge ~ M1a, M1a_start, m=2)
    
    M2a_start <- list(update(M0, dnds=0.1), update(M0, dnds=1), 
                      update(M0, dnds=3))
    M2a <- pmlMix(edge ~., M2a_start, m=3)
    
    #    attr(RESULT, "env") <- env 
    #    class(RESULT) <- c("modelTest", "data.frame")
    #    RESULT
    
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
    
    


