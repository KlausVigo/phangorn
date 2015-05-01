#  computes all n(n-1)/2 cyclic splits
cyclicSplits <- function(k, labels=NULL){
    k = as.integer(k)
    l = (k-1L) %/% 2L
    res <- vector("list", k*(k-1L)/2)
    res[1:k] = 1L:k
    ind = k
    if(k>3){
        fun = function(x,y){
            tmp = (1L:y)+x
            tmp %% (k+1L) + tmp %/% (k+1L)
        }
        if(k>4){
            for(i in 2:l){
                res[(ind+1):(ind+k)] <- lapply(0L:(k-1L), fun, i)
                ind <- ind+k
            }
        }
        if((k%%2L)==0){
            m <- k%/%2
            res[(ind+1):(ind+m)] <- lapply(0L:(m-1L), fun, m)
        }        
    }   
    if(is.null(labels)) labels=(as.character(1:k))
    attr(res, 'labels') =labels
    attr(res, "cycle") = 1:k
    class(res)="splits"
    res   
}


distC <- function(d, CL){
  l=length(CL)
  res = matrix(0, l, l)
  for(i in 1:(l-1)){
    for(j in (i+1):l)
      res[i,j] = mean.default(d[CL[[i]], CL[[j]]])
  }
  res + t(res)
}


reduc <- function(d, x, y, z){
  u <- 2/3* d[x, ] + d[y,]/3
  v <- 2/3* d[z, ] + d[y,]/3
  uv <- (d[x,y] + d[x,z] + d[y,z])/3 
  d[x, ] <- u
  d[, x] <- u
  d[z, ] <- v
  d[, z] <- v
  
  d[y, ] <- 0
  d[, y] <- 0
  
  d[x, z] <- d[z, x] <- uv
  diag(d) <- 0
  d 
}


# computes ordering
getOrderingNN <- function (x) 
{
  x = as.matrix(x)
  labels <- attr(x, "Labels")
  if (is.null(labels)) 
    labels = colnames(x)
  d = x #as.matrix(x)
  l = dim(d)[1]
  CL = vector("list", l)  
  CL[1:l] <- ORD <- 1:l
  lCL <- length(CL)
  ord <- CL   
  while (lCL>1){
    i = 0
    j = 0
 #   browser()
    DM = distC(d, CL)
    l = nrow(DM)
    if(l>2){
    r = rowSums(DM)/(l - 2)
    tmp <- .C("out", as.double(DM), as.double(r), as.integer(l), 
              as.integer(i), as.integer(j), PACKAGE = "phangorn")
    e1 = tmp[[4]]
    e2 = tmp[[5]]
    }
    else {e1 = 1
     e2=2}
    n1 <- length(CL[[e1]])
    n2 <- length(CL[[e2]])
    if(n1==1 & n2==1){
      newCL <- c(CL[[e1]], CL[[e2]])
      newOrd = newCL
      CL = c(CL[-c(e1,e2)], list(newCL))
      ord <- c(ord[-c(e1,e2)], list(newCL))
      lCL <- lCL - 1L
    }
    else{
      CLtmp = c(as.list(CL[[e1]]), as.list(CL[[e2]]), CL[-c(e1,e2)])
      ltmp =length(CLtmp)
      DM2 = distC(d, CLtmp)
      if(ltmp>2) rtmp = rowSums(DM2)/(ltmp - 2)
      DM2 = DM2 - outer(rtmp, rtmp, "+")
      
      TMP = DM2[1:n1, (n1+1):(n1+n2)]
#browser()
#      dtmp = d[CL[[e1]], CL[[e2]]]
#      rtmp = numeric(n1+n2)
#      for(ii in 1:(n1+n2)){
#          for(jj in 1:ltmp){if(ii!=jj) rtmp[ii]=rtmp[ii] + mean.default(d[CLtmp[[ii]], CLtmp[[jj]]])
#        }
#      }
#browser()      
#      rtmp = rtmp/(ltmp-2)
#      TMP2  = dtmp + rep(rtmp[1:n1],n2) + rep(rtmp[(n1+1):(n1+n2)], each=n1) 

#browser()

      blub = which.min(TMP)
#      print(blub)
#print("blub")      
      if(n1==2 & n2==1){
        if(blub == 2){
          newCL <- c(CL[[e1]][1], CL[[e2]])
          newOrd <-  c(CL[[e1]], ord[[e2]]) 
          d <- reduc(d, CL[[e1]][1], CL[[e1]][2], CL[[e2]]) 
        }
        else{
          newCL <- c(CL[[e2]], CL[[e1]][2])
          newOrd <- c(ord[[e2]], ord[[e1]])
          d <- reduc(d, CL[[e2]], CL[[e1]][1], CL[[e1]][2]) 
        }

       
      }
      if(n1==1 & n2==2){
        if(blub==1){
          newCL <- c(CL[[e1]], CL[[e2]][2])
          newOrd <-  c(CL[[e1]], ord[[e2]])
          d <- reduc(d, CL[[e1]], CL[[e2]][1], CL[[e2]][2])
        }
        else{
          newCL <- c(CL[[e2]][1], CL[[e1]])
          newOrd <- c(ord[[e2]], ord[[e1]])
          d <- reduc(d, CL[[e2]][1], CL[[e2]][2], CL[[e1]])
        }
        }
      if(n1==2 & n2==2){
        if(blub==1){
          newCL <- c(CL[[e1]][2], CL[[e2]][2])
          newOrd <-  c(rev(ord[[e1]]), ord[[e2]])
          d <- reduc(d, CL[[e1]][2], CL[[e1]][1], CL[[e2]][1]) 
          d <- reduc(d, CL[[e1]][2], CL[[e2]][1], CL[[e2]][2]) 
        }
        if(blub==2){
          newCL <- c(CL[[e1]][1], CL[[e2]][2])
          newOrd <-  c(ord[[e1]], ord[[e2]])      
          d <- reduc(d, CL[[e1]][1], CL[[e1]][2], CL[[e2]][1]) 
          d <- reduc(d, CL[[e1]][1], CL[[e2]][1], CL[[e2]][2]) 
          
        }
        if(blub==3){
          newCL <- c(CL[[e1]][2], CL[[e2]][1])
          newOrd <-  c(rev(ord[[e1]]), rev(ord[[e2]]))
          d <- reduc(d, CL[[e1]][2], CL[[e1]][1], CL[[e2]][2]) 
          d <- reduc(d, CL[[e1]][2], CL[[e2]][2], CL[[e2]][1]) 
        }
        if(blub==4){
            newCL <- c(CL[[e1]][1], CL[[e2]][1])
            newOrd <-  c(ord[[e1]], rev(ord[[e2]])) 
            d <- reduc(d, CL[[e1]][1], CL[[e1]][2], CL[[e2]][2]) 
            d <- reduc(d, CL[[e1]][1], CL[[e2]][2], CL[[e2]][1]) 
            }
        }
        CL <- c(CL[-c(e1,e2)], list(newCL))
        ord <- c(ord[-c(e1,e2)], list(newOrd))
        lCL <- lCL - 1L
        }
    }
    newOrd
} 

#
neighborNet <-  function(x, ord=NULL){
    x = as.matrix(x)
    labels <- attr(x, "Labels")[[1]]
    if (is.null(labels)) 
        labels = colnames(x)
    l <- length(labels)    
#browser()    
    if(is.null(ord))ord <- getOrderingNN(x)
    spl <- cyclicSplits(l, labels[ord])
    spl <- nnls.splits(spl, x)
    # nnls.split mit nnls statt quadprog
    attr(spl, "cycle") <- 1:l
    as.networx(spl)
} 


