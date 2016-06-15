nodeHeight <- function(tree) 
{
    if(is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
        tree <- reorder(tree, "postorder")
    edge = tree$edge[, 2]
    node = tree$edge[, 1]
    m <- max(tree$edge)
    el = double(m)
    el[edge] = tree$edge.length
    res = .C("nodeH", as.integer(edge), as.integer(node), el, as.integer(length(edge)), double(m))[[5]] 
    max(res) - res
}



ancstat = function(phy, x){
  contrast= attr(x, "contrast")
  storage.mode(contrast) = "integer"
  phy=reorder(phy, "postorder")
  res=matrix(0L, max(phy$edge), ncol(contrast))
  colnames(res) = attr(x,"levels")
  nTips=length(phy$tip.label)
  pa=phy$edge[,1]
  ch=phy$edge[,2]
  res[1:nTips, ] = contrast[as.numeric(x)[match(phy$tip.label, names(x))],, drop=FALSE]
  for(i in 1:length(pa)){
    res[pa[i],] = res[pa[i],] | res[ch[i],]    
  }
  res
}


comp <- function(x, y){
  tmp1 = matrix(rowSums(x), nrow(x), nrow(y))
  res = matrix(rowSums(y), nrow(x), nrow(y), byrow=TRUE)
  tmp3 = tcrossprod(x, 1-y)  
  tmp0 = tcrossprod(x, y)
  tmp0[tmp3>0]=0L
  res[!(tmp0>(tmp1 - 1e-8))] = 10000000L 
  apply(res, 1, which.min)
}


comp2 <- function(x, y){
  res = matrix(rowSums(x), nrow(x), nrow(y))
  tmp1 = matrix(rowSums(y), nrow(x), nrow(y), byrow=TRUE)
  tmp3 = tcrossprod(1-x, y)  
  tmp0 = tcrossprod(x, y)
  tmp0[tmp3>0]=0L
  res[tmp0<2] = Inf
  apply(res, 2, which.min)
}

# single linkage of minimal coalescent times
# extends speciesTree fom ape

coalSpeciesTree <- function(tree, X, sTree=NULL){
  
  if(is.null(X))return(speciesTree(tree))  
  trees = unclass(tree)
  States = lapply(tree, ancstat, X)
  NH = lapply(tree, nodeHeight)
  if(is.null(sTree)){
    l <- attr(X, "nc")
    m <- choose(l, 2)
    SST <- matrix(0L, m, l)
    k <- 1
    for(i in 1:(l-1)){
      for(j in (i+1):l){
        SST[k, i] <- SST[k,j] <- 1L
        k <- k+1
      }
    }
    Y=matrix(Inf, length(NH), nrow(SST)) 
    dm = rep(Inf, m)
    for(i in 1:length(NH)){
      ind = comp2(States[[i]],SST)  
      dm = pmin(dm, NH[[i]][ind])
      #       for(j in 1:length(ind))Y[i, ind[j]] = min(Y[i, ind[j]], NH[[i]][j])
    }
    dm = structure(2*dm, Labels = attr(X, "levels"), Size = l, class = "dist", Diag = FALSE, Upper = FALSE)
    
    sTree <- upgma(dm, "single")   
    # dm of pairwise states
  }   
  else{ 
    SST = ancstat(sTree, X)
    Y=matrix(Inf, length(NH), nrow(SST)) 
    for(i in 1:length(NH)){
      ind = comp(States[[i]],SST) 
      for(j in 1:length(ind))Y[i, ind[j]] = min(Y[i, ind[j]], NH[[i]][j])
    }
    STH = apply(Y, 2, min)
    sTree$edge.length = STH[sTree$edge[,1]] - STH[sTree$edge[,2]]
  }
  sTree
}

