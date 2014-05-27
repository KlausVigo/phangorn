SOWH.test <- function(x, n=100, restricted=list(optNni=FALSE), optNni=TRUE, ...){
  
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
        optRate = optR$optRate, optRooted = optR$optRooted, model = optR$model)    
  unrestr <- optim.pml(restr, optNni = optU$optNni, optBf = optU$optBf, optQ = optU$optQ, 
        optInv = optU$optInv, optGamma = optU$optGamma, optEdge = optU$optEdge, 
        optRate = optU$optRate, optRooted = optU$optRooted, model = optU$model) 
  
  for(i in 1:n){
    newData <- simSeq(restr)
    restrTmp <- update(restr, data=newData)
    unrestrTmp <- restrTmp # update(unrestr, data=newData)
    restrTmp <- optim.pml(restrTmp, optNni = optR$optNni, optBf = optR$optBf, optQ = optR$optQ, 
        optInv = optR$optInv, optGamma = optR$optGamma, optEdge = optR$optEdge, 
        optRate = optR$optRate, optRooted = optR$optRooted, model = optR$model)  
    unrestrTmp <- optim.pml(unrestrTmp, optNni = optU$optNni, optBf = optU$optBf, optQ = optU$optQ, 
        optInv = optU$optInv, optGamma = optU$optGamma, optEdge = optU$optEdge, 
        optRate = optU$optRate, optRooted = optU$optRooted, model = optU$model) 
    res[i, 1] <- logLik(restrTmp)
    res[i, 2] <- logLik(unrestrTmp)
    
  }
  list(restr, unrestr, res)
}
