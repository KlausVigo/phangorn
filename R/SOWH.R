#' Swofford-Olsen-Waddell-Hillis Test
#' 
#' This function computes the Swofford--Olsen--Waddell--Hillis (SOWH) test, a
#' parametric bootstrap test. The function is computational very demanding and
#' likely to be very slow.
#' 
#' \code{SOWH.test} performs a parametric bootstrap test to compare two trees.
#' It makes extensive use \code{simSeq} and \code{optim.pml} and can take quite
#' long.
#' 
#' @param x an object of class \code{"pml"}.
#' @param n the number of bootstrap replicates.
#' @param restricted list of restricted parameter settings.
#' @param optNni Logical value indicating whether topology gets optimized
#' (NNI).
#' @param trace Show output during computations.
#' @param \dots Further arguments passed to \code{"optim.pml"}.
#' @return an object of class SOWH. That is a list with three elements, one is
#' a matrix containing for each bootstrap replicate the (log-) likelihood of
#' the restricted and unrestricted estimate and two pml objects of the
#' restricted and unrestricted model.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}}, \code{\link{pmlPart}}, \code{\link{pmlCluster}},
#' \code{\link{simSeq}}, \code{\link{SH.test}}
#' @references Goldman, N., Anderson, J. P., and Rodrigo, A. G. (2000)
#' Likelihood -based tests of topologies in phylogenetics. \emph{Systematic
#' Biology} \bold{49} 652-670.
#' 
#' Swofford, D.L., Olsen, G.J., Waddell, P.J. and Hillis, D.M. (1996)
#' Phylogenetic Inference in Hillis, D.M., Moritz, C. and Mable, B.K. (Eds.)
#' \emph{Molecular Systematics} (2nd ed.) 407-514, Sunderland, MA: Sinauer
#' @keywords models
#' @examples
#' 
#' # in real analysis use larger n, e.g. 500 preferably more
#' \dontrun{
#' data(Laurasiatherian)
#' dm <- dist.logDet(Laurasiatherian)
#' tree <- NJ(dm)
#' fit <- pml(tree, Laurasiatherian)
#' fit <- optim.pml(fit, TRUE)
#' set.seed(6)
#' tree <- rNNI(fit$tree, 1)
#' fit <- update(fit, tree = tree)
#' (res <- SOWH.test(fit, n=100))
#' summary(res)
#' }
#' 
#' @export SOWH.test
SOWH.test <- function(x, n=100, restricted=list(optNni=FALSE), optNni=TRUE, trace = 1, ...){
  
  res <- matrix(NA, n, 2)
  extras <- match.call(expand.dots = FALSE)$...
  
  optU <- list (optNni = optNni, optBf = FALSE, optQ = FALSE, 
          optInv = FALSE, optGamma = FALSE, optEdge = TRUE, optRate = FALSE, 
          optRooted = FALSE, model = NULL)

  if(!is.null(extras)){
      namAll <-  names(extras)
      for(i in 1: length(extras))optU[[namAll[i]]] <- extras[[i]]    
  }   
  optR <- optU
  namR <- names(restricted)   
  for(i in 1: length(namR))optR[[namR[i]]] <- restricted[[i]]
  restr <- optim.pml(x, optNni = optR$optNni, optBf = optR$optBf, optQ = optR$optQ, 
        optInv = optR$optInv, optGamma = optR$optGamma, optEdge = optR$optEdge, 
        optRate = optR$optRate, optRooted = optR$optRooted, model = optR$model, 
        pml.control(trace = trace-1L))    
  unrestr <- optim.pml(restr, optNni = optU$optNni, optBf = optU$optBf, optQ = optU$optQ, 
        optInv = optU$optInv, optGamma = optU$optGamma, optEdge = optU$optEdge, 
        optRate = optU$optRate, optRooted = optU$optRooted, model = optU$model, 
        pml.control(trace = trace-1L)) 
  
  for(i in 1:n){
    if(trace>0) cat("iteration: ", i, "\n")  
    newData <- simSeq(restr)
    restrTmp <- update(restr, data=newData)
    unrestrTmp <- restrTmp # update(unrestr, data=newData)
    restrTmp <- optim.pml(restrTmp, optNni = optR$optNni, optBf = optR$optBf, optQ = optR$optQ, 
        optInv = optR$optInv, optGamma = optR$optGamma, optEdge = optR$optEdge, 
        optRate = optR$optRate, optRooted = optR$optRooted, model = optR$model, 
        pml.control(trace = trace-1L))  
    unrestrTmp <- optim.pml(unrestrTmp, optNni = optU$optNni, optBf = optU$optBf, optQ = optU$optQ, 
        optInv = optU$optInv, optGamma = optU$optGamma, optEdge = optU$optEdge, 
        optRate = optU$optRate, optRooted = optU$optRooted, model = optU$model, 
        pml.control(trace = trace-1L)) 
    res[i, 1] <- logLik(restrTmp)
    res[i, 2] <- logLik(unrestrTmp)  
  }
  result <- list("LL"=res, "restr" = restr, "unrestr" = unrestr)
  class(result) <- "SOWH" 
  result
}


print.SOWH <- function(x, digits = 4L, ...){
    resLL <- logLik(x$restr)  
    unresLL <- logLik(x$unrestr) 
    diffLL <- unresLL - resLL
    pval <- sum( (x$LL[,2] - x$LL[,1]) > diffLL) / nrow(x$LL)
    res <- c(resLL, unresLL, diffLL, pval)
    names(res) <- c("ln L restr", "ln L unrestr", "Diff ln L", "p-value")
    print(res, digits=digits)
    invisible(x)
}


summary.SOWH <- function(object, digits = 4L, plot=TRUE, ...){
    resLL <- logLik(object$restr)  
    unresLL <- logLik(object$unrestr) 
    diffLL <- unresLL - resLL
    pval <- sum( (object$LL[,2] - object$LL[,1]) > diffLL) / nrow(object$LL)
    res <- c(resLL, unresLL, diffLL, pval)
    names(res) <- c("ln L restr", "ln L unrestr", "Diff ln L", "p-value")
    print(res, digits=digits)
    if(plot){
        d <- object$LL[,2] - object$LL[,1]
        hist( d, freq=FALSE, xlim=c(0, 1.2 * max(d, diffLL)))
        abline(v=diffLL, col="red")
    } 
    invisible(object)
}



