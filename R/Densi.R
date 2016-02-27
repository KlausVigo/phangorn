getAges <- function(x){  
  fun=function(x) max(node.depth.edgelength(x))  
  height=NULL
  if(inherits(x,"phylo")) height <- fun(x)
  if(inherits(x,"multiPhylo")){
    if(!is.null(attr(x, "TipLabel"))){
      x = unclass(x)
      x = .uncompressTipLabel(x)  
      x = unclass(x)  
      height = sapply(x, fun)
    }
    else{
      x = unclass(x)
      height = sapply(x, fun) 
    }
  }
  height
}


# now more memoryefficient
# from phytools code by Liam Revell with a few changes
my.supertree<-function(trees,method=c("pratchet","optim.parsimony"), trace=0, ...){
  # set method
  method<-method[1]
  # some minor error checking
  if(!inherits(trees,"multiPhylo")) stop("trees must be object of class 'multiPhylo.'")
  
  labels <- lapply(trees, function(x)sort(x$tip.label))
  ulabels <- unique(labels)
  lul <- length(ulabels)
  # compute matrix representation phylogenies
  X<-vector("list", lul) # list of bipartitions
  characters<-0 # number of characters
  weights <- NULL
  species<-trees[[1]]$tip.label
  for(i in 1:lul){
    pos <- match(labels, ulabels[i])
    ind <- which(!is.na(pos))  
    temp<-prop.part(trees[ind]) # find all bipartitions
    # create matrix representation of trees[[i]] in X[[i]]
    X[[i]]<-matrix(0,nrow=length(trees[[ind[1]]]$tip),ncol=length(temp)-1)
    for(j in 1:ncol(X[[i]])) X[[i]][c(temp[[j+1]]),j]<-1
    rownames(X[[i]])<-attr(temp,"labels") # label rows
#    if(i==1) species<-trees[[ind[1]]]$tip.label
#    else 
        species<-union(species,trees[[ind[1]]]$tip.label) # accumulate labels
    characters<-characters+ncol(X[[i]]) # count characters
    weights <- c(weights, attr(temp, "number")[-1])
  }
  XX<-matrix(data="?",nrow=length(species),ncol=characters,dimnames=list(species))
  j<-1
  for(i in 1:length(X)){
    # copy each of X into supermatrix XX
    XX[rownames(X[[i]]),c(j:((j-1)+ncol(X[[i]])))]<-X[[i]][1:nrow(X[[i]]),1:ncol(X[[i]])]
    j<-j+ncol(X[[i]])
  }
  # compute contrast matrix for phangorn
  contrast<-matrix(data=c(1,0,0,1,1,1),3,2,dimnames=list(c("0","1","?"),c("0","1")),byrow=TRUE)
  # convert XX to phyDat object
  XX<-phyDat(XX,type="USER",contrast=contrast, compress=FALSE) 
  attr(XX, "weight") <- weights 
  # estimate supertree
  if(method=="pratchet"){
    if(hasArg(start)){
      start<-list(...)$start
      if(inherits(start,"phylo")){
        supertree<-pratchet(XX,all=TRUE, trace=0, ...)
      } else {
        if(start=="NJ") start<-NJ(dist.hamming(XX))
        else if(start=="random") start<-rtree(n=length(XX),tip.label=names(XX))
        else {
          warning("do not recognize that option for start; using random starting tree")
          tree<-rtree(n=length(XX),tip.label=names(XX))
        }
        args<-list(...)
        args$start<-start
        args$data<-XX
        args$all<-TRUE
        supertree<-do.call(pratchet,args)
      }
    } else supertree<-pratchet(XX,all=TRUE, trace=0, ...)
    if(inherits(supertree,"phylo"))
      if(trace>0)message(paste("The MRP supertree, optimized via pratchet(),\nhas a parsimony score of ",
                    attr(supertree,"pscore")," (minimum ",characters,")",sep=""))
    else if(inherits(supertree,"multiPhylo"))
      if(trace>0)message(paste("pratchet() found ",length(supertree)," supertrees\nwith a parsimony score of ",
                    attr(supertree[[1]],"pscore")," (minimum ",characters,")",sep=""))
  } else if(method=="optim.parsimony"){
    if(hasArg(start)){
      start<-list(...)$start
      if(inherits(start,"phylo")){
        supertree<-optim.parsimony(tree=start,data=XX, trace=0, ...)
      } else {
        if(start=="NJ") start<-NJ(dist.hamming(XX))
        else if(start=="random") start<-rtree(n=length(XX),tip.label=names(XX))
        else {
          warning("do not recognize that option for tree; using random starting tree")
          start<-rtree(n=length(XX),tip.label=names(XX))
        }
        supertree<-optim.parsimony(tree=start,data=XX,...)
      }			
    } else {
      if(trace>0)message("no input starting tree or option for optim.parsimony; using random addition tree")
      start<-random.addition(XX) # rtree(n=length(XX),tip.label=names(XX))
      supertree<-optim.parsimony(tree=start,data=XX, trace=0, ...)
    }
    if(inherits(supertree,"phylo"))
      if(trace>0)message(paste("The MRP supertree, optimized via optim.parsimony(),\nhas a parsimony score of ",
                    attr(supertree,"pscore")," (minimum ",characters,")",sep=""))
    else if(inherits(supertree,"multiPhylo"))
      if(trace>0)message(paste("optim.parsimony() found ",length(supertree)," supertrees\nwith a parsimony score of ",
                    attr(supertree[[1]],"pscore")," (minimum ",characters,")",sep=""))
  }
  return(supertree)
}


my.supertree.Old<-function(trees,method=c("pratchet","optim.parsimony"), trace=0, ...){
    # set method
    method<-method[1]
    # some minor error checking
    if(!inherits(trees,"multiPhylo")) stop("trees must be object of class 'multiPhylo.'")
    # compute matrix representation phylogenies
    X<-list() # list of bipartitions
    characters<-0 # number of characters
    for(i in 1:length(trees)){
        temp<-prop.part(trees[[i]]) # find all bipartitions
        # create matrix representation of trees[[i]] in X[[i]]
        X[[i]]<-matrix(0,nrow=length(trees[[i]]$tip),ncol=length(temp)-1)
        for(j in 1:ncol(X[[i]])) X[[i]][c(temp[[j+1]]),j]<-1
        rownames(X[[i]])<-attr(temp,"labels") # label rows
        if(i==1) species<-trees[[i]]$tip.label
        else species<-union(species,trees[[i]]$tip.label) # accumulate labels
        characters<-characters+ncol(X[[i]]) # count characters
    }
    XX<-matrix(data="?",nrow=length(species),ncol=characters,dimnames=list(species))
    j<-1
    for(i in 1:length(X)){
        # copy each of X into supermatrix XX
        XX[rownames(X[[i]]),c(j:((j-1)+ncol(X[[i]])))]<-X[[i]][1:nrow(X[[i]]),1:ncol(X[[i]])]
        j<-j+ncol(X[[i]])
    }
    # compute contrast matrix for phangorn
    contrast<-matrix(data=c(1,0,0,1,1,1),3,2,dimnames=list(c("0","1","?"),c("0","1")),byrow=TRUE)
    # convert XX to phyDat object
    XX<-phyDat(XX,type="USER",contrast=contrast) 
    # estimate supertree
    if(method=="pratchet"){
        if(hasArg(start)){
            start<-list(...)$start
            if(inherits(start,"phylo")){
                supertree<-pratchet(XX,all=TRUE, trace=0, ...)
            } else {
                if(start=="NJ") start<-NJ(dist.hamming(XX))
                else if(start=="random") start<-rtree(n=length(XX),tip.label=names(XX))
                else {
                    warning("do not recognize that option for start; using random starting tree")
                    tree<-rtree(n=length(XX),tip.label=names(XX))
                }
                args<-list(...)
                args$start<-start
                args$data<-XX
                args$all<-TRUE
                supertree<-do.call(pratchet,args)
            }
        } else supertree<-pratchet(XX,all=TRUE, trace=0, ...)
        if(inherits(supertree,"phylo"))
            if(trace>0)message(paste("The MRP supertree, optimized via pratchet(),\nhas a parsimony score of ",
                                     attr(supertree,"pscore")," (minimum ",characters,")",sep=""))
        else if(inherits(supertree,"multiPhylo"))
            if(trace>0)message(paste("pratchet() found ",length(supertree)," supertrees\nwith a parsimony score of ",
                                     attr(supertree[[1]],"pscore")," (minimum ",characters,")",sep=""))
    } else if(method=="optim.parsimony"){
        if(hasArg(start)){
            start<-list(...)$start
            if(inherits(start,"phylo")){
                supertree<-optim.parsimony(tree=start,data=XX, trace=0, ...)
            } else {
                if(start=="NJ") start<-NJ(dist.hamming(XX))
                else if(start=="random") start<-rtree(n=length(XX),tip.label=names(XX))
                else {
                    warning("do not recognize that option for tree; using random starting tree")
                    start<-rtree(n=length(XX),tip.label=names(XX))
                }
                supertree<-optim.parsimony(tree=start,data=XX,...)
            }			
        } else {
            if(trace>0)message("no input starting tree or option for optim.parsimony; using random addition tree")
            start<-random.addition(XX) # rtree(n=length(XX),tip.label=names(XX))
            supertree<-optim.parsimony(tree=start,data=XX, trace=0, ...)
        }
        if(inherits(supertree,"phylo"))
            if(trace>0)message(paste("The MRP supertree, optimized via optim.parsimony(),\nhas a parsimony score of ",
                                     attr(supertree,"pscore")," (minimum ",characters,")",sep=""))
        else if(inherits(supertree,"multiPhylo"))
            if(trace>0)message(paste("optim.parsimony() found ",length(supertree)," supertrees\nwith a parsimony score of ",
                                     attr(supertree[[1]],"pscore")," (minimum ",characters,")",sep=""))
    }
    return(supertree)
}


# we want a rooted supertree
superTree = function(tree, method="pratchet", rooted=TRUE, ...){
  fun = function(x){
    x=reorder(x, "postorder")
    nTips = length(x$tip)
    x$edge[x$edge>nTips] = x$edge[x$edge>nTips] + 2L
    l=nrow(x$edge)
    oldroot = x$edge[l,1L]
    x$edge=rbind(x$edge,matrix(c(rep(nTips+2,2),oldroot,nTips+1),2L,2L))
    x$edge.length=c(x$edge.length, 100, 100)
    x$tip.label=c(x$tip.label, "ZZZ")
    x$Nnode=x$Nnode+1L
    x
  }
  if(!is.null(attr(tree, "TipLabel")))tree = .uncompressTipLabel(tree)
  tree = unclass(tree)
  if(rooted) tree = lapply(tree, fun)    
  class(tree)="multiPhylo"
  res = my.supertree(tree, method=method, ...)
  if(rooted){
    if(inherits(res,"multiPhylo")){
      res = lapply(res, root, "ZZZ")
      res = lapply(res, drop.tip, "ZZZ")  
      class(res) = "multiPhylo"
    }
    else{
      res = root(res, "ZZZ")
      res = drop.tip(res, "ZZZ")  
    }
  }
  if(inherits(res,"multiPhylo")){
    fun = function(x){
      x$edge.length <- rep(.1, nrow(x$edge)) 
      x
    }
    res <- lapply(res, fun)
    res <- lapply(res, reorder, "postorder")
    class(res) = "multiPhylo"
  }       
  else{ 
    res$edge.length = rep(.1, nrow(res$edge))
    res <- reorder(res, "postorder")
  }
  res
}



densiTree <- function(x, type="cladogram", alpha=1/length(x), consensus=NULL, optim=FALSE, scaleX=FALSE, col=1, width=1, cex=.8, ...) {
  if(!inherits(x,"multiPhylo"))stop("x must be of class multiPhylo")
  compressed <- ifelse(is.null(attr(x, "TipLabel")), FALSE, TRUE)
  if(is.null(consensus))consensus <- superTree(x)
  if(inherits(consensus,"multiPhylo")) consensus = consensus[[1]]
  type <- match.arg(type, c("phylogram", "cladogram"))
  consensus = reorder(consensus, "postorder")
  e2 = reorder(consensus)$edge[,2]
  nTip = as.integer(length(consensus$tip))
  tiporder = e2[e2<=nTip]   
  maxBT = max(getAges(x))
  if(scaleX) maxBT=1.0
  label = rev(pretty(c(maxBT,0)))
  maxBT = max(label)
  xy = plotPhyloCoor(consensus, ...)
  yy = xy[,2]
  plot.new() 
  tl = which.max(nchar(consensus$tip.label))
  sw <- strwidth(consensus$tip.label[tl],cex=cex) * 1.1
  plot.window(xlim=c(0, 1.0+sw), ylim=c(0, nTip+1))
  axis(side=1,at=seq(0,1.0, length.out=length(label)), labels=label)
  text(x=rep(1.0,Ntip(consensus)),y=yy[1:nTip],labels=consensus$tip.label,pos=4,cex=cex)  
  tip.order = yy[1:nTip]
  for (treeindex in 1:length(x)) {
    tmp <- reorder(x[[treeindex]], "postorder")
    if(!compressed) tip.order <- match(tmp$tip.label, consensus$tip.label)
    xy <- plotPhyloCoor(tmp, tip.order=tiporder, ...)
    xx = xy[,1]
    yy = xy[,2]
    if(scaleX) xx <- xx/max(xx)
    else xx <- xx/maxBT 
    xx <- xx + (1.0 - max(xx))
    e1=tmp$edge[,1]
    e2=tmp$edge[,2]
    if(type=="cladogram") cladogram.plot(tmp$edge, xx, yy, edge.color=adjustcolor(col, alpha.f=alpha), edge.width=width, edge.lty=1)
    if(type=="phylogram"){
      Ntip <- min(e1)-1L 
      Nnode <- tmp$Nnode 
      phylogram.plot(tmp$edge, Ntip, Nnode, xx, yy, TRUE, edge.color=adjustcolor(col, alpha.f=alpha), edge.width=width, 1) 
    }
  }  
}






