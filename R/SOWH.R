SOWH.test <- function(x, n=100, restricted=list(optNni=FALSE), optNni=TRUE, trace = 1, ...){
  
  res = matrix(NA, n, 2)
  extras <- match.call(expand.dots = FALSE)$...
  
  optU = list (optNni = optNni, optBf = FALSE, optQ = FALSE, 
          optInv = FALSE, optGamma = FALSE, optEdge = TRUE, optRate = FALSE, 
          optRooted = FALSE, model = NULL)

  if(!is.null(extras)){
      namAll =  names(extras)
      for(i in 1: length(extras))optU[[namAll[i]]] = extras[[i]]    
  }   
  optR = optU
  namR = names(restricted)   
  for(i in 1: length(namR))optR[[namR[i]]] = restricted[[i]]
#  x <- optim.pml()
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
  result <- res 
  attr(result, "restr") <- restr
  attr(result, "unrestr") <- unrestr 
  class(result) = "SOWH" # list(restr, unrestr, res)
  result
}

print.SOWH <- function(x, digits = 4L, ...){
    resLL = logLik(attr(x, "restr"))
    unresLL = logLik(attr(x, "unrestr"))
    diffLL = unresLL - resLL
    pval <- sum( (x[,2] - x[,1]) > diffLL) / nrow(x)
    res = c(resLL, unresLL, diffLL, pval)
    names(res) = c("ln L restr", "ln L unrestr", "Diff ln L", "p-value")
    print(res, digits=digits)
    invisible(x)
}

#summary.SOWH <- function(x, digits = 4L, ...){
#    resLL = logLik(attr(x, "restr"))
#    unresLL = logLik(attr(x, "unrestr"))
#    diffLL = resll - unresll
#    pval <- sum( (x[,1] - x[,2]) > diffLL) / nrow(x)
#    res = c(resLL, unresLL, diffLL, pval)
#    print(res, digits=digits)
#    hist(x[,1] - x[,2])
#    return(NULL)
#}



