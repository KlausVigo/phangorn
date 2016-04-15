#
# splits format, networx, Matrix, lento plot 
#
as.splits <- function (x, ...){
    if(inherits(x, "splits")) return(x)
    UseMethod("as.splits")
}


as.Matrix <- function (x, ...){
    if (inherits(x,"Matrix")) return(x)
    UseMethod("as.Matrix")
}


as.matrix.splits <- function(x, zero.print = 0L, one.print=1L, ...){
   m = length(x)
   labels = attr(x, "labels")
   n = length(labels)    
   res = matrix(zero.print, m, n)
   for(i in 1:m)res[i,x[[i]]]=one.print
   dimnames(res) = list(names(x), labels)
   res
}


as.Matrix.splits <- function(x, ...){
    labels = attr(x, "labels")
    l = length(x)
    j = unlist(x)
    i = rep(1:l, sapply(x, length))
    sparseMatrix(i,j, x = rep(1L, length(i)), dimnames = list(NULL, labels)) # included x und labels
}


print.splits <- function (x, maxp = getOption("max.print"), 
    zero.print = ".", one.print="|", ...)
{
    x.orig <- x
    cx <- as.matrix(x, zero.print = zero.print, one.print=one.print)
    print(cx, quote = FALSE, right = TRUE, max = maxp)
    invisible(x.orig)
}


"[.splits" = function(x, i){
    tmp = attributes(x)
    result = unclass(x)[i]
    if(!is.null(tmp$weights)) tmp$weights = tmp$weights[i] 
    if(!is.null(tmp$confidences)) tmp$confidences = tmp$confidences[i]
    if(!is.null(tmp$intervals)) tmp$intervals = tmp$intervals[i] 
    if(!is.null(tmp$data)) tmp$data = tmp$data[i,, drop=FALSE] 
    attributes(result) = tmp
    result
}

changeOrder <- function(x, labels){
    oldL <- attr(x, "labels")
    ind <- match(oldL,labels)
    for(i in 1:length(x))
        x[[i]] <- sort(ind[x[[i]]])
    if(!is.null(attr(x, "cycle")))
        attr(x, "cycle") <- ind[attr(x, "cycle")]
    attr(x, "labels") <- labels
    x    
}


#orderSplitLabel = function(x, order){
#    label = attr(x, "labels")
#    nTips = length(label)
#    ord = match(label, order)
#    for(i in 1:length(x))
#        x[[i]] = sort(ord[x[[i]]])
#    attr(x, "labels") = order
#    x
#}


# returns order of x$edge
presenceAbsenceOld <- function(x, y){
    X <- as.splits(x)
    Y <- as.splits(y)
    labels <- attr(X, "labels") 
    #    if(inherits(x,"phylo")) X <- X[x$edge[,2]]
    #    if(inherits(y,"phylo")) Y <- Y[y$edge[,2]]
    Y <- changeOrder(Y, labels) # orderSplitLabel
    nTips <- length(labels)
    X <- oneWise(X, nTips)
    Y <- oneWise(Y, nTips)
    res <- match(X, Y)    
    res <- !is.na(res)
    if(inherits(x, "networx")){
        res <- res[x$splitIndex]    
    }    
    if(class(x)[1]=="phylo"){
        # res <- res[x$edge[,2]]
        x$node.label = res[-c(1:length(labels))]
        return(x)
    }
    res            
}


presenceAbsence <- function(x,y){
    spl <- as.splits(y)
    l <- length(spl)
    attr(spl, "confidences") <- rep(1, l)
    addConfidences(x, y)
}


matchSplits <- function(x, y, as.in=TRUE){
    tiplabel <- attr(x, "label")
    if(any(is.na(match(tiplabel, attr(y, "label"))))) stop("x and y have different labels!")
    nTips <- length(tiplabel)
    y <- changeOrder(y, tiplabel)
    y <- SHORTwise(y, nTips)
    if(as.in) return(match(SHORTwise(x, nTips), y, nomatch = 0L) > 0L)
    match(SHORTwise(x, nTips), y)
}


optCycle <- function(splits, tree){
    tips = tree$tip.label
    tree = reorder(tree)
    nodes = sort(unique(tree$edge[,1]))
    
    M = as.matrix(splits)
    
    l = as.integer(nrow(M))
    m = as.integer(ncol(M))

    tmp = tree$edge[,2]
    tmp = tmp[tmp<=m]

    start <- .C("countCycle", M[, tmp], l, m, integer(1))[[4]]
    best = start
    eps = 1
    if(eps>0){
        for(i in 1:length(nodes)){
           tmptree = rotate(tree, nodes[i])
           tmp = tmptree$edge[,2]
           tmp = tmp[tmp<=m]
           tmpC <- .C("countCycle", M[, tmp], l, m, integer(1))[[4]]
           if(tmpC < best){
              best <- tmpC
              tree = tmptree
           }
        }
        eps = start - best
    }
    tree # list(best, tree)
}


countCycles <- function(splits, tree=NULL, ord=NULL){
  M = as.matrix(splits)
  l = as.integer(nrow(M))
  m = as.integer(ncol(M))
  if(!is.null(tree))  ord  = getOrdering(tree)
  res <- .C("countCycle2", M[, ord], l, m, integer(l))[[4]]
  res
}

  
c.splits <- function (..., recursive=FALSE) 
{
    x <- list(...)
    n <- length(x)
    match.names <- function(a, b) {
        if (any(!(a %in% b))) 
            stop("names do not match previous names")
    }
    if (n == 1) 
        return(x[[1]])

    labels <- attr(x[[1]], "labels")
    cycle <- attr(x[[1]], "cycle")
    for (i in 2:n) {
        match.names(labels, attr(x[[i]], "labels"))
    }
    res = structure(NextMethod("c"), class=c("splits", "prop.part"))
    attr(res, "labels") = labels
    attr(res, "weight") = as.vector(sapply(x, attr, "weight"))
    attr(res, "cycle") = cycle
    res
}


## as.splits.phylo
# computes splits from phylo
as.splits.phylo <- function(x, ...){
    result <- bip(x)
    if(!is.null(x$edge.length)){
        edge.weights <- numeric(max(x$edge))
        edge.weights[x$edge[,2]] <- x$edge.length
        attr(result, "weights") <- edge.weights
    }
    if(!is.null(x$node.label)){
        conf <- x$node.label
        if(is.character(conf)) conf <- as.numeric(conf)
        if(max(na.omit(conf)) > (1 + 1e-8))conf <- conf / 100
        #if(!is.null(scale)) conf <- conf / scale
        attr(result, "confidences") <- c(rep(1, length(x$tip.label)), conf)
#        attr(result, "confidences") <- c(rep("", length(x$tip.label)), x$node.label)
    }    
    attr(result, "labels") <- x$tip
    class(result) <- c('splits', 'prop.part')
    result 
}


# computes splits from multiPhylo object (e.g. bootstrap, MCMC etc.)
# unrooted trees
as.splits.multiPhylo <- function(x, ...){
#    if(inherits(x,"multiPhylo"))x = .uncompressTipLabel(x)
#    if(inherits(x,"multiPhylo"))class(x)='list'  # prop.part allows not yet multiPhylo
#    firstTip = x[[1]]$tip[1]
#    x = lapply(x, root, firstTip) # old trick 
    lx <-  length(x)
    x <- lapply(x, unroot)
    class(x) <- "multiPhylo"
    splits <- prop.part(x)
    class(splits)='list'
    weights = attr(splits, 'number')    
    lab = attr(splits,'labels')
    attr(splits,'labels') <- attr(splits, 'number') <- NULL
    l = length(lab)
    splitTips = vector('list', l)
    for(i in 1:l) splitTips[[i]] = i
    result = c(splitTips,splits)
    attr(result, "weights") = c(rep(lx, l), weights)
    attr(result, "confidences") <- attr(result, "weights") / lx
    attr(result, "summary") <- list(confidences="ratio", ntrees=l, clades=FALSE) 
    attr(result, "labels") <- lab
    class(result) = c('splits', 'prop.part')
    result  
}


as.splits.prop.part <- function(x, ...){
    if(is.null(attr(x, "number")))  
        attr(x, "weights") = rep(1, length(x)) 
	 else{ 
        attr(x, "weights") = attr(x, "number")
        attr(x, "confidences") = attr(x, "number") / attr(x, "number")[1] 
   	}    
    class(x) = c('splits', 'prop.part')	
    x
}


as.splits.networx <- function(x, ...){
#    if(!is.null(attr(x, "splits")))attr(x, "splits")
    if(!is.null(x$splits)) x$splits
    else warning("No split object included!")    
}


as.prop.part.splits <- function(x, ...){
    attr(x, "number") = attr(x, "weights")
    attr(x, "weights") = NULL
    class(x) = c('prop.part')	
    x
}

## as.splits.phylo
as.phylo.splits <- function (x, result = "phylo", ...) 
{
    result <- match.arg(result, c("phylo", "all"))
    labels = attr(x, "labels")
    nTips = length(labels)
    weights = attr(x, "weights")
    nTips = length(labels)
    x = SHORTwise(x, nTips)
    dm = as.matrix(compatible(x))
    rs = rowSums(dm)
    ind = which(rs == 0)
    if (any(rs > 0)) {
        tmp = which(rs > 0)
        candidates = tmp[order(rs[tmp])]
        for (i in candidates) {
            if (sum(dm[ind, i]) == 0) 
                ind = c(ind, i)
        }
    }
    splits = x[ind]
    weights = weights[ind]
    l = length(ind)
    res = matrix(0L, l, nTips)
    for (i in 1:l) res[i, splits[[i]]] = 1L
    dm2 = (crossprod(res * weights, 1 - res))
    dm2 = dm2 + t(dm2)
    dimnames(dm2) = list(labels, labels)
    tree <- di2multi(NJ(dm2), tol = 1e-08)
    attr(tree, "order") = NULL
    tree <- reorder(tree)    
    tree <- optCycle(x, tree)
    tree <- reorder(tree, "postorder")
    if (result == "phylo") 
        return(tree)  
#    tree = reroot(tree, Ancestors(tree, 1, "parent")) 
    spl = as.splits(tree)
    spl = SHORTwise(spl, nTips)
    spl <- spl[tree$edge[,2]]
    list(tree = tree, index = tree$edge[, 2], split = spl, rest = x[-ind])
}


# computes compatible splits
compatible <- function(obj){
    labels = attr(obj, "labels")
    if(!inherits(obj, "splits"))stop("obj needs to be of class splits")
    
    l = length(labels)
    n = length(obj)
    
    bp = matrix(0L, n, l)
    for(i in 1:n)bp[i,obj[[i]]] = 1L
    bp[bp[, 1] == 0L, ] = 1L - bp[bp[, 1] == 0L, ]
    k=1
    res = matrix(0L, n, n) 
            
    tmp1 = tcrossprod(bp) #sum(bp[i,]* bp[j,])
    tmp2 = tcrossprod(1L - bp) #sum((1L - bp[i,])*(1L - bp[j,]))
    tmp3 = tcrossprod(bp, 1L - bp) #sum(bp[i,]*(1L - bp[j,]))
    tmp4 = tcrossprod(1L - bp, bp) #sum((1L - bp[i,])*bp[j,]) 
    res[(tmp1 * tmp2 * tmp3 * tmp4)>0]=1L
    k = k+1
    
    res = res[lower.tri(res)]
    attr(res, "Size") <- n
    attr(res, "Diag") <- FALSE
    attr(res, "Upper") <- FALSE
    class(res) <- "dist"
    return(res)
}

    
compatible2 <- function (obj1, obj2=NULL) 
{   
    if (!inherits(obj1, "splits")) 
        stop("obj needs to be of class splits")
    labels = attr(obj1, "labels")    
    l = length(labels)
    n = length(obj1)
    bp1 = as.matrix(obj1)
    bp1[bp1[, 1] == 0L, ] = 1L - bp1[bp1[, 1] == 0L, ] 
    if(!is.null(obj2)){
        m = length(obj2) 
        bp2 = as.matrix(obj2)
        labels2 = attr(obj2, "labels")
        bp2 = bp2[, match(labels2, labels), drop=FALSE]
        bp2[bp2[, 1] == 0L, ] = 1L - bp2[bp2[, 1] == 0L, ]
    }
    else bp2 = bp1

    if(is.null(obj2)) res = matrix(0L, n, n)
    else res = matrix(0L, n, m)

    tmp1 = tcrossprod(bp1, bp2)
    tmp2 = tcrossprod(1L - bp1, 1L - bp2)
    tmp3 = tcrossprod(bp1, 1L - bp2)
    tmp4 = tcrossprod(1L - bp1, bp2)
    res[(tmp1 * tmp2 * tmp3 * tmp4) > 0] = 1L
    if(is.null(obj2)){
        res = res[lower.tri(res)]
        attr(res, "Size") <- n
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        class(res) <- "dist"
    }
    return(res)
}


compatible3 <- function(x, y=NULL) 
{
    if (!inherits(x, "splits")) 
        stop("x needs to be of class splits")
    if(is.null(y)) y <- x
        
    if (!inherits(y, "splits")) 
        stop("y needs to be of class splits")
    xlabels = attr(x, "labels")
    ylabels = attr(y, "labels")
    if(identical(xlabels, ylabels)) labels = xlabels 
    else labels = intersect(xlabels, ylabels)
    nx = length(x)
    ny = length(y)   
    bp1 = as.matrix(x)[,labels, drop=FALSE]
    bp2 = as.matrix(y)[,labels, drop=FALSE]
    rs1 = rowSums(bp1)
    rs2 = rowSums(bp2)
    res = matrix(0L, nx, ny)
    tmp1 = tcrossprod(bp1, bp2)
    res = matrix(0L, nx, ny)
    for(i in 1:nx){
        for(j in 1:ny){            
            if(tmp1[i, j]==rs1[i]) res[i,j] = 1
            if(tmp1[i, j]==rs2[j]) res[i,j] = 2
            if(tmp1[i, j]==rs1[i] & tmp1[i, j]==rs2[j])res[i,j] = 3
        }
    }      
    if(is.null(y)){
        res = res[lower.tri(res)]
        attr(res, "Size") <- length(x)
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        class(res) <- "dist"
    }
    return(res)
}
    

#
# splits
#
splitsNetwork <- function(dm, splits=NULL, gamma=.1, lambda=1e-6, weight=NULL){
  dm = as.matrix(dm)
  k = dim(dm)[1]
  
  if(!is.null(splits)){
    tmp = which(sapply(splits, length)==k)
    splits = splits[-tmp]
    lab = attr(splits, "labels")
    dm = dm[lab, lab]
  }
  
  if(is.null(splits)){
    X2 = designAll(k, TRUE)
    X=X2[[1]]
  }
  else X = as.matrix(splits2design(splits))
  
  y = dm[lower.tri(dm)]
  if(is.null(splits))ind = c(2^(0:(k-2)),2^(k-1)-1)
  else ind = which(sapply(splits, length)==1)
  #   y2 = lm(y~X[,ind]-1)$res
  n = dim(X)[2]
  
  ridge <- lambda * diag(n) 
  ridge[ind,ind] <- 0
  if(!is.null(weight)) Dmat <- crossprod(X * sqrt(weight)) + ridge
  else Dmat <- crossprod(X) + ridge
  if(!is.null(weight)) dvec <- crossprod(X * sqrt(weight),y * sqrt(weight))
  else dvec <- crossprod(X, y)
  
  #    Dmat <- as.matrix(Dmat)
  #    dvec <- as.vector(dvec) 
  
  ind1       <- rep(1,n)
  ind1[ind]  <- 0 
  
  Amat       <- cbind(ind1,diag(n)) 
  bvec       <- c(gamma, rep(0,n))
  
  solution <- quadprog::solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=1)$sol   
  
  ind2 <- which(solution>1e-8)
  n2 <- length(ind2)
  
  ind3 = which(duplicated(c(ind2, ind), fromLast = TRUE)[1:n2])
  ridge2 <- lambda * diag(n2) 
  ridge2[ind3,ind3] <- 0
  
  if(!is.null(weight)) Dmat <- crossprod(X[, ind2] * sqrt(weight)) + ridge2
  else Dmat <- crossprod(X[, ind2]) + ridge2
  if(!is.null(weight)) dvec <- crossprod(X[, ind2] * sqrt(weight),y * sqrt(weight))
  else dvec <- crossprod(X[, ind2], y)
  
  Amat2 <- diag(n2)
  bvec2 <- rep(0, n2)
  solution2  <- quadprog::solve.QP(Dmat, dvec, Amat2)$sol
  
  RSS1 = sum((y-X[,ind2]%*%solution[ind2])^2)
  RSS2 = sum((y-X[,ind2]%*%solution2)^2)
  
  if(is.null(splits)){
    splits = vector("list", length(ind2))
    for(i in 1:length(ind2))splits[[i]] = which(X2[[2]][ind2[i],]==1)
  } 
  else splits = splits[ind2]
  attr(splits, "weights") = solution[ind2]
  attr(splits, "unrestricted") = solution2
  attr(splits, "stats") = c(df=n2, RSS_p = RSS1, RSS_u=RSS2)
  attr(splits,"labels") =dimnames(dm)[[1]]
  class(splits)='splits'
  return(splits)           
}


allSplits = function(k, labels=NULL){
  result <- lapply(1:(2^(k-1)-1),dec2Bin)
  if(is.null(labels)) labels=(as.character(1:k))
  attr(result, 'labels') =labels
  class(result)='splits'
  result
}   


allCircularSplits <- function(k, labels=NULL){
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
        for(i in 2:l){
            res[(ind+1):(ind+k)] <- lapply(0L:(k-1L), fun, i)
            ind <- ind+k
        }
        if((k%%2L)==0){
            m <- k%/%2
            res[(ind+1):(ind+m)] <- lapply(0L:(m-1L), fun, m)
        }
        
    }   
    if(is.null(labels)) labels=(as.character(1:k))
    attr(res, 'labels') =labels
    class(res)="splits"
    res   
}


getIndex = function(left, right, n){
  if(n<max(left) | n<max(right)) stop("Error")  
  left = as.integer(left)
  right = as.integer(right)
  ll = length(left)
  lr = length(right)
  .C("giveIndex", left, right, ll, lr, as.integer(n), integer(ll*lr))[[6]]+1
}


splits2design <- function(obj, weight=NULL){
  labels= attr(obj,'labels')
  m = length(labels)
  n=length(obj)
  l = 1:m 
  sl = sapply(obj, length)
  p0 = sl * (m-sl)
  p = c(0,cumsum(p0))
  i = numeric(max(p))
  for(k in 1:n){
    sp = obj[[k]]
    if(p0[k]!=0) i[(p[k]+1):p[k+1]] = getIndex(sp, l[-sp], m) 
  }
  dims=c(m*(m-1)/2,n)
  sparseMatrix(i=i, p=p, x=1.0, dims=dims) 
}


hC <- function(g, set){
    intersec = NULL
    allEdges = NULL
    fromTo <- set
    l = length(set)
    sptmp = shortest_paths(g, fromTo[l], fromTo[1], output=c("epath"))$epath[[1]]
    sptmp = as.vector(sptmp)
    allEdges = sptmp
    for(i in 2:length(set)){
        sptmp = shortest_paths(g, fromTo[i-1], fromTo[i], output=c("epath"))$epath[[1]]
        sptmp = as.vector(sptmp)
        intersec = c(intersec, intersect(allEdges, sptmp) )
        allEdges = c(allEdges, sptmp)
    }
    #    allEdges = unique(allEdges)
    list(allEdges, unique(allEdges), intersec)
}   


addEdge <- function(network, desc, spl){   
    edge <- network$edge
    parent <- edge[,1]
    child <- edge[,2]
    nTips <- length(network$tip.label)

    desc2 <- SHORTwise(desc, nTips)    
    split <- desc2[spl]
        
    index <- network$splitIndex
    ind <- which(compatible2(split, desc2[index]) == 1)
    if(is.null(ind) | (length(ind)==0)) return(network)
    add <- TRUE
  
    X <- as.matrix(desc2)
    rsX <- rowSums(X)
    z <- X %*% X[spl,]
    v <- which((rsX == z)[index] == TRUE) 

    
# intersection of shortest pathes of both partitions
# best with similar to circNetwork with shortest_paths 
    
    while(add){
        tmp = ind
        for(i in ind){          
            tmp2 = which(compatible2(desc2[index][i], desc2[index]) == 1)
            tmp = union(tmp, tmp2)
        }
        if(identical(ind, tmp)){
            ind=tmp           
            add=FALSE
        }
        ind=tmp
    }    
   

    g = graph(t(network$edge[ind,]), directed=FALSE)
    dec = decompose(g, min.vertices = 2)

    #    fromTo <- sort(match(split[[1]], attr(desc, "cycle")))
    #    sptmp = shortest_paths(g, fromTo[i-1], fromTo[i], 
    #                           output=c("epath"))$epath[[1]]
    #    sp2 = c(sp2, sptmp[-c(1, length(sptmp))])
    #    sp0 = c(sp0, sptmp)
    
    oldNodes = unique(as.vector(edge[ind,]))
    mNodes = max(network$edge)
    newNodes = (mNodes+1L) : (mNodes+length(oldNodes))

# duplicated splits
    dSpl = edge[ind,]
    edge2 = edge[v,] 
    for(i in 1:length(oldNodes)){
        edge2[edge2 == oldNodes[i]] = newNodes[i]
    } 
    edge[v,] = edge2    

  #alle Splits verdoppeln
    for(i in 1:length(oldNodes)) dSpl[dSpl==oldNodes[i]] = newNodes[i]
    edge = rbind(edge, dSpl, deparse.level = 0) # experimental: no labels
    index = c(index, index[ind])
  #neu zu alt verbinden   
    edge = rbind(edge, cbind(oldNodes, newNodes), deparse.level = 0) 
    index = c(index, rep(spl, length(oldNodes)) )
    network$splitIndex = index
    network$edge = edge
    network$Nnode = max(edge) - nTips
    network   
}

## as.splits.phylo
circNetwork <- function(x, ord=NULL){
    if(is.null(ord))ord = attr(x, "cycle")
    
    weight <- attr(x, "weights")
    if(is.null(weight)) weight = rep(1, length(x))
    nTips = length(ord)
    tmp = which(ord == 1)
    if(tmp!=1) ord = c(ord[tmp:nTips], ord[1:(tmp-1)])
    res = stree(nTips, tip.label = attr(x, "labels"))
    res$edge[, 2] = ord
    res$edge.length=NULL
    x <- SHORTwise(x, nTips)    
    spRes <- as.splits(res)[res$edge[,2]]
    index = match(spRes, x)
    
    if(any(is.na(index))){
        l.na = sum(is.na(index))
        x <- c(x, spRes[is.na(index)])    
        weight = c(weight, rep(0, l.na))
        index = match(spRes, x)
    }
    
    l = sapply(oneWise(x, nTips), length)
    l2 = sapply(x, length)
    #    dm <- as.matrix(compatible2(x))
    
    tmp <- countCycles(x, ord=ord)
    ind = which(tmp == 2 & l2>1) # & l<nTips changed with ordering
    
    ind = ind[order(l[ind])]
    
    dm2 <- as.matrix(compatible2(x, x[ind]))
    
    X = as.matrix(x)[,ord]
    Y = X    
    rsY = rowSums(Y)
    X = X[ind, ]
    
    for(k in 1: length(ind)){
        Vstart = ord[1]
        Vstop = ord[nTips]    
        ordStart = 1
        ordStop = nTips
        for(j in 2:nTips){
            
            if(X[k,j-1] < X[k,j]){ 
                Vstart = ord[j]
                ordStart = j                   
            }                       
            if(X[k,j-1] > X[k,j]){ 
                Vstop = ord[j-1]
                ordStop = j-1   
            }    
        } 
        
        fromTo <- ordStart:ordStop
        if(ordStart>ordStop) fromTo <- c(ordStart:nTips, 1:ordStop)
        fromTo = ord[fromTo] 
#        print(fromTo) 
        g = graph(t(res$edge), directed=FALSE)
        
        isChild = (rsY == (Y %*% X[k,]))[index]
        sp2 = NULL
        sp0 = NULL
        
        for(i in 2:length(fromTo)){
#            sptmp = get.shortest.paths(g, fromTo[i-1], fromTo[i], 
#                                       output=c("epath"))$epath[[1]]            
            sptmp = shortest_paths(g, fromTo[i-1], fromTo[i], 
                                       output=c("epath"))$epath[[1]]
            sp2 = c(sp2, sptmp[-c(1, length(sptmp))])
            sp0 = c(sp0, sptmp)
        }
        sp0 = unique(sp0)
        
        if(length(sp2)>0){
            #            blub = which(dm[index[sp2], ind[k]]>0)
            TMP = rowSums(dm2[index[sp2], 1:k, drop=FALSE])
            blub = which(TMP>0)
            sp2 = sp2[blub]
        }
        if(length(sp2)==0){
            isChild = (rsY == (Y %*% X[k,]))[index]  
            sp0 = which(isChild == TRUE)
            edge1 = unique(as.vector(res$edge[sp0,]))
            edge2 = as.vector(res$edge[-sp0,])
            asdf = edge1 %in% edge2
            sp = edge1[asdf]
        }
        if(length(sp2)>0)   sp = unique(as.vector(t(res$edge[sp2,])))     
        parent = res$edge[,1]
        child = res$edge[,2]    
        
        j = ord[which(X[k,]==1)]
        anc = unique(parent[match(j, child)])
        
        maxVert = max(parent)
        l = length(sp)
        
        newVert = (maxVert+1) : (maxVert+l)      
        sp01 = setdiff(sp0, sp2)
        for(i in 1:l) res$edge[sp01,][res$edge[sp01,]==sp[i]] = newVert[i] 
        
        newindex = rep(ind[k], l)        
        if(length(sp)>1)newindex = c(index[sp2], newindex)
        index = c(index, newindex)        
        # connect new and old vertices
        newEdge = matrix(cbind(sp, newVert), ncol=2) 
        if(length(sp)>1){
            # copy edges
            qwer = match(as.vector(res$edge[sp2,]), sp)
            newEdge = rbind(matrix(newVert[qwer], ncol=2), newEdge)
        }
        
        res$edge = rbind(res$edge, newEdge)      
        res$Nnode =  max(res$edge) - nTips
        
        res$splitIndex = index
        res$edge.length <- rep(1, nrow(res$edge))
        class(res) = c("networx", "phylo")
        attr(res, "order") = NULL
        #browser() 
    }
    res$Nnode =  max(res$edge) - nTips
    res$edge.length = weight[index]  # ausserhalb
    res$splitIndex = index 
#    attr(res, "splits") = x
    res$splits = x
    class(res) = c("networx", "phylo")
    attr(res, "order") = NULL
    res    
}


as.networx <- function (x, ...) 
{
    if (inherits(x, "networx")) 
        return(x)
    UseMethod("as.networx")
}


getOrdering <- function(x){
    tree = as.phylo(x)
    nTips = length(tree$tip)
    ord = reorder(tree)$edge[,2]
    ord = ord[ord<=nTips]
    ind = which(ord == 1L)
    if(ind>1) ord = c(ord[ind:nTips], ord[c(1:(ind-1L))])
    ord  
}

## as.splits.phylo
addTrivialSplits <- function(obj){
    label <- attr(obj, "label")
    nTips <- length(label)
    weight <- attr(obj, "weights")
    if(is.null(weight)) weight = rep(1, length(obj))
    STree = stree(nTips, tip.label = attr(obj, "labels"))
    STree$edge.length=NULL 
    spRes <- as.splits(STree)[STree$edge[,2]]
    tmpIndex = match(spRes, SHORTwise(obj, nTips))
    if(any(is.na(tmpIndex))){
        l.na = sum(is.na(tmpIndex))
        obj <- c(obj, spRes[is.na(tmpIndex)]) 
        weight = c(weight, rep(0, l.na))
        attr(obj, "weights") <- weight
    }
    obj
}


as.networx.splits <- function(x, planar=FALSE, coord = c("none", "2D", "3D"), ...){
  label <- attr(x, "label")
  
  x = addTrivialSplits(x)
  
  nTips <- length(label)
  weight <- attr(x, "weights")
  if(is.null(weight)) weight = rep(1, length(x))
  attr(x, "weights") <- weight
  
  x <- oneWise(x, nTips) 
  l <- sapply(x, length)
  if(any(l==nTips))x <- x[l!=nTips] # get rid of trivial splits
  ext <- sum(l==1 | l==(nTips-1))
  if(!is.null(attr(x, "cycle"))){  
      c.ord <- attr(x, "cycle") 
  }
  else c.ord <- getOrdering(x)
  attr(x, "cycle") = c.ord
  
  dm <- as.matrix(compatible2(x)) 
# which splits are in circular ordering  
    circSplits = which(countCycles(x, ord=c.ord)==2) 
    if(length(circSplits) == length(x)) planar=TRUE
    tmp = circNetwork(x, c.ord)  
    attr(tmp, "order") = NULL
    if(planar){
        return(reorder(tmp))
    }

    ll <- sapply(x, length)
    ind <- tmp$splitIndex     # match(sp, x)
    ind2 = union(ind, which(ll==0)) # which(duplicated(x))
    ind2 = union(ind2, which(ll==nTips))
    ord <- order(colSums(dm))
    ord <- setdiff(ord, ind2)
    if(length(ord)>0){    
        for(i in 1:length(ord)){ 
            tmp = addEdge(tmp, x, ord[i])
            tmp$edge.length = weight[tmp$splitIndex]
            tmp$Nnode = max(tmp$edge) - nTips
            class(tmp) = c("networx", "phylo")
        } 
    }
    tmp$Nnode = max(tmp$edge) - nTips
    tmp$edge.length = weight[tmp$splitIndex]
    attr(x, "cycle") <- c.ord
#    attr(tmp, "splits") = x 
    tmp$splits <- x
    class(tmp) = c("networx", "phylo")
    tmp <- reorder(tmp)
    coord <- match.arg(coord)
    vert <- switch(coord,
           "none" = NULL,
           "2D" = coords(tmp, dim="2D"),
           "3D" = coords(tmp, dim="3D"))
#    attr(tmp, "coords") <- coordinates
    tmp$plot <- list(vertices=vert)
    tmp
}


#as.networx.phylo <- function(x, ...){
#    spl <- as.splits(x)
#    as.networx(x, ...)
#}

## as.splits.phylo
as.networx.phylo <- function(x, ...){
    spl <- as.splits(x)
    spl <- spl[x$tree[,2]]
    x$splitIndex <- 1:nrow(x$edge)
    attr(x, "splits") = spl
    class(x) <- c("networx", "phylo")
    x
}


#as.igraph.networx <- function(x, directed=FALSE){
#    graph(t(x$edge), directed=directed)
#}


consensusNet <- function (obj, prob = 0.3, ...) 
{
    l = length(obj)
    spl = as.splits(obj)
    w = attr(spl, "weights")
    ind = (w/l) > prob
    spl = spl[ind]
    attr(spl, "confidences") = (w/l)[ind]
#    attr(spl, "weights") = w[ind]
    res = as.networx(spl)  
    res$edge.labels = as.character(res$edge.length / l * 100)
    res$edge.labels[res$edge[,2]<=length(res$tip.label)] = ""
    reorder(res)
}


createLabel <- function(x, y, label_y, type="edge", nomatch=NA){
    spl_x <- as.splits(x)
    if(inherits(x, "phylo", TRUE)==1) spl_x <- spl_x[x$edge[,2]]
    spl_y <- as.splits(y)
    if(inherits(y, "phylo", TRUE)==1) spl_y <- spl_y[y$edge[,2]]
    
    tiplabel <- attr(spl_x, "label")
    nTips <- length(tiplabel)
    
    spl_y <- changeOrder(spl_y, tiplabel)
    spl_y <- SHORTwise(spl_y, nTips)
    
    ind <- match(SHORTwise(spl_x, nTips), spl_y)
    pos <-  which(!is.na(ind))

    res <- rep(nomatch, length(spl_x))
    
    if(length(label_y)==1L) label_y <- rep(label_y, length(spl_y))
    res[pos] <- label_y[ind[pos]]
    if(type=="edge" && inherits(x, "networx")){
        return(res[x$splitIndex])
    }
    res  
}


addConfidences <- function (x, y, ...) UseMethod("addConfidences")


# y now more general 
addConfidences.splits <- function(x, y, ...){
    tiplabel <- attr(x, "label")
    nTips = length(tiplabel)
#    x = addTrivialSplits(x) 
    if(inherits(y,"phylo")){
        ind <- match(tiplabel, y$tip.label)
        if (any(is.na(ind)) | length(tiplabel) != length(y$tip.label)) 
            stop("trees have different labels")
        y$tip.label <- y$tip.label[ind]
        ind2 <- match(1:length(ind), y$edge[, 2])
        y$edge[ind2, 2] <- order(ind)
    }
    spl <- as.splits(y)
    spl <- changeOrder(spl, tiplabel)
    spl <- SHORTwise(spl, nTips)
#    ind <- match(SHORTwise(x, nTips), spl)
#    ind
    ind <- match(SHORTwise(x, nTips), spl)
    #    pos <-  which(ind > nTips)
    pos <-  which(!is.na(ind))
    confidences <- numeric(length(x))  #character 
    confidences[pos] <- attr(spl, "confidences")[ind[pos]]
    #        y$node.label[ind[pos] - nTips]
    attr(x, "confidences") <- confidences
    x  
}


addConfidences.networx <- function(x, y, ...){
    spl <- x$splits
    spl <- addConfidences(spl, y, ...)
    x$splits <- spl
    x    
}

## as.splits.phylo
addConfidences.phylo <- function(x, y, ...){
#    call <- x$call
    conf = attr(addConfidences(as.splits(x), y), "confidences")
    if(is.character(conf)) conf <- as.numeric(conf)
    nTips = length(x$tip.label)
    x$node.label = conf[-c(1:nTips)]  * 100
    x      
} 


reorder.networx <- function (x, order =  "cladewise", index.only = FALSE, ...) 
{
    order <- match.arg(order, c("cladewise", "postorder"))
    if (!is.null(attr(x, "order"))) 
        if (attr(x, "order") == order) 
            return(x)    
    g <- graph(t(x$edge))
#    if(order == "cladewise") neword <- topological.sort(g, "out")
#    else neword <- topological.sort(g, "in") 
    if(order == "cladewise") neword <- topo_sort(g, "out")
    else neword <- topo_sort(g, "in") 
    neworder <- order(match(x$edge[,1], neword))
    if(index.only) return(neworder)
    x$edge <- x$edge[neworder, ]
    if (!is.null(x$edge.length)) 
        x$edge.length <- x$edge.length[neworder]
    if (!is.null(x$edge.labels)) 
        x$edge.labels <- x$edge.labels[neworder]  
    if (!is.null(x$splitIndex))x$splitIndex <- x$splitIndex[neworder]
    attr(x, "order") <- order
    x
}


coords <- function(obj, dim="3D"){
#    if(is.null(attr(obj,"order")) || (attr(obj, "order")=="postorder") ) 
#        obj = reorder.networx(obj)

    l = length(obj$edge.length)
    ind1 = which(!duplicated(obj$splitIndex))

    n = max(obj$edge)
    adj = spMatrix(n, n, i = obj$edge[,2], j = obj$edge[,1], x = rep(1, length(obj$edge.length))) # Matrix::
    g = graph_from_adjacency_matrix(adj, "undirected")
#    g = graph.adjacency(adj, "undirected")
##########
#    add this 
#    g2 <- graph(t(obj$edge), directed=FALSE)
#    g2 <- set.edge.attribute(g, "weight", value=rep(1, nrow(obj$edge))
    if(dim=="3D"){
        coord <- layout_with_kk(g, dim=3)
#        coord <- layout.kamada.kawai(g, dim=3)
        k = matrix(0, max(obj$splitIndex), 3)
        for(i in ind1){
            tmp = coord[obj$edge[i, 2],] - coord[obj$edge[i, 1],]
            k[obj$splitIndex[i], ] = kart2kugel(tmp[1], tmp[2], tmp[3])
        }
        k[obj$splitIndex[ind1],1] = obj$edge.length[ind1] 

        res = matrix(0, vcount(g), 3)
        for(i in 1:l){# unique(obj$splitIndex)
            j = obj$edge[i,1]
            m = obj$edge[i,2]
            p = obj$splitIndex[i]
            res[m,] = res[j,] + kugel2kart(k[p,1], k[p,2], k[p,3])     
        }            
    }
    else{
        coord <- layout_with_kk(g, dim=2)
#        coord <- layout.kamada.kawai(g, dim=2)
        k = matrix(0, max(obj$splitIndex), 2)
        for(i in ind1){
            tmp = coord[obj$edge[i, 2],] - coord[obj$edge[i, 1],]
            k[obj$splitIndex[i], ] = kart2kreis(tmp[1], tmp[2])
        }
        k[obj$splitIndex[ind1],1] = obj$edge.length[ind1] 
        res = matrix(0, vcount(g), 2)
        for(i in 1:l){# unique(obj$splitIndex)
            j = obj$edge[i,1]
            m = obj$edge[i,2]
            p = obj$splitIndex[i]
            res[m,] = res[j,] + kreis2kart(k[p,1], k[p,2])     
        }
    }  
    res  
}


kart2kugel <- function(x,y,z){
    r = sqrt(x*x+y*y+z*z)
    alpha = atan(sqrt(x*x+y*y) / z)
    if(z<0) alpha = alpha+pi
    beta = atan(y/x)
    if(x<0) beta = beta+pi 
    c(r,alpha,beta)
}

	
kart2kreis <- function(x,y){
    r = sqrt(x*x+y*y)
    alpha = atan(y/x) 
    if(x<0) alpha = alpha+pi
    c(r,alpha)
}	
	

kreis2kart <- function(r,alpha){
	c(r*cos(alpha), r*sin(alpha))
}


kugel2kart <- function(r,alpha,beta){
    x = r * sin(alpha) * cos(beta) 
    y = r * sin(alpha) * sin(beta) 
    z = r * cos(alpha)
    c(x,y,z)
}


edgeLabels <- function(xx,yy,zz=NULL, edge){
        XX <- (xx[edge[, 1]] + xx[edge[, 2]])/2
        YY <- (yy[edge[, 1]] + yy[edge[, 2]])/2
        if(!is.null(zz)){
	        ZZ <- (zz[edge[, 1]] + zz[edge[, 2]])/2
	        return(cbind(XX, YY, ZZ))
        }  
        cbind(XX, YY)  
}


plot.networx = function(x, type="3D", use.edge.length = TRUE, show.tip.label=TRUE,
    show.edge.label=FALSE, edge.label=NULL, show.node.label = FALSE, node.label=NULL,
    show.nodes=FALSE, tip.color = "black", 
    edge.color="black", edge.width = 3, edge.lty = 1,
    split.color=NULL, split.width = NULL, split.lty = NULL,
    font = 3, cex = par("cex"), 
    cex.node.label=cex, cex.edge.label=cex,
    col.node.label = tip.color, col.edge.label = tip.color, 
    font.node.label = font, font.edge.label = font,
    ...){
    type = match.arg(type, c("3D", "2D")) 
    if(use.edge.length==FALSE) x$edge.length[] = 1
# test    
#    x = reorder(x)
    nTips = length(x$tip.label)
    conf = attr(x$splits,"confidences") * 100
    index = x$splitIndex
    if(is.null(edge.label) & !is.null(conf))edge.label = conf[index]
    if(is.null(node.label))node.label = as.character(1:max(x$edge))
    if(show.tip.label)node.label[1:nTips] = ""
    
    lspl <- max(x$splitIndex)
    if(!is.null(split.color)){
        if(length(split.color)!=lspl) stop("split.color must be same length as splits")
        else edge.color <- split.color[x$splitIndex]
    } 
    if(!is.null(split.width)){
        if(length(split.width)!=lspl) stop("split.color must be same length as splits")
        else edge.width <- split.width[x$splitIndex]
    } 
    if(!is.null(split.lty)){
        if(length(split.lty)!=lspl) stop("split.color must be same length as splits")
        else edge.lty <- split.lty[x$splitIndex]
    } 
    
    chk <- FALSE
    
    if(type=="3D") chk <- requireNamespace("rgl", quietly = TRUE) #.check.pkg("rgl")
    if(!chk && type=="3D"){
        warning("type=\"3D\" requires the package \"rgl\"\n, plotting =\"2D\" instead!\n")
        type="2D"
    }
    # use precomputed vertices when available
    coord <- NULL
    if(!is.null(x$.plot)) coord <- x$.plot$vertices
    
    if(type=="3D") {
        if(is.null(coord) || ncol(coord)!=3)        
            coord <- coords(x, dim="3D")
        plotRGL(coord, x, show.tip.label=show.tip.label, show.edge.label=show.edge.label, 
             edge.label = edge.label, show.node.label = show.node.label, node.label=node.label, 
             show.nodes=show.nodes, tip.color = tip.color, edge.color=edge.color, 
             edge.width = edge.width, font = font, cex = cex, 
             cex.node.label=cex.node.label, cex.edge.label=cex.edge.label,
             col.node.label = col.node.label, col.edge.label = col.edge.label,
             font.node.label = font.node.label, font.edge.label = font.edge.label)
    }
    else{
        if(is.null(coord) || ncol(coord)!=2)
     	    coord <- coords(x, dim="2D")
	    plot2D(coord, x, show.tip.label=show.tip.label, show.edge.label=show.edge.label, 
	        edge.label = edge.label, show.node.label = show.node.label, node.label=node.label,
	        show.nodes=show.nodes, tip.color = tip.color, edge.color=edge.color,
	        edge.width = edge.width, edge.lty=edge.lty,font = font, cex = cex, 
	        cex.node.label=cex.node.label, cex.edge.label=cex.edge.label,
	        col.node.label = col.node.label, col.edge.label = col.edge.label,
	        font.node.label = font.node.label, font.edge.label = font.edge.label,
	        add=FALSE)
    }   
    x$.plot <- list(vertices = coord, edge.color=edge.color, edge.width=edge.width, edge.lty = edge.lty)
    invisible(x)
}

    
plotRGL <- function(coords, net, show.tip.label=TRUE, 
        show.edge.label=FALSE, edge.label=NULL, show.node.label=FALSE, node.label=NULL,
        show.nodes=FALSE, tip.color = "blue", edge.color="grey", 
        edge.width = 3, font = 3, cex = par("cex"), 
        cex.node.label=cex,  cex.edge.label=cex,
        col.node.label=tip.color, col.edge.label=tip.color,
        font.node.label=font, font.edge.label=font,        
        ...){
    
#    chk <- .check.pkg("rgl")
#    if(!chk) open3d <- segments3d <- spheres3d <- rgl.texts <- function(...)NULL

    open3d <- rgl::open3d 
    segments3d  <- rgl::segments3d
    spheres3d <- rgl::spheres3d 
    rgl.texts <- rgl::rgl.texts
        
    edge = net$edge
  
    x = coords[,1]
    y = coords[,2]
    z = coords[,3]
     
    nTips = length(net$tip.label)
  
    segments3d(x[t(edge)],y[t(edge)],z[t(edge)], col=rep(edge.color, each=2), lwd=edge.width) 
    radius=0
    if(show.nodes){
        radius = sqrt((max(x)-min(x))^2 + (max(y)-min(y))^2 + (max(z)-min(z))^2) / 200    
        spheres3d(x[1:nTips], y[1:nTips],z[1:nTips], radius=2*radius, color="cyan")
        spheres3d(x[-c(1:nTips)], y[-c(1:nTips)],z[-c(1:nTips)], radius=radius, color="magenta")
    }
    if(show.tip.label){
      rgl.texts(x[1:nTips]+2.05*radius,y[1:nTips],z[1:nTips],net$tip.label, color=tip.color, cex=cex, font=font)
    }
    if(show.edge.label){
	    ec = edgeLabels(x, y, z, edge)
      if(is.null(edge.label)) edge.label = net$splitIndex
        #else edge.label = net$splitIndex    
	    rgl.texts(ec[,1], ec[,2], ec[,3], edge.label, color=col.edge.label, cex=cex.edge.label, font=font.edge.label)     
    } 
    if(show.node.label){
        rgl.texts(x, y, z, node.label, color=col.node.label, cex=cex.node.label, font=font.node.label) 
    }
}


plot2D <- function(coords, net, show.tip.label=TRUE,  
       show.edge.label=FALSE, edge.label=NULL, show.node.label=FALSE, node.label=NULL,
       tip.color = "blue", edge.color="grey",                   
       edge.width = 3, edge.lty=1, 
       font = 3, cex = par("cex"), 
       cex.node.label=cex,  cex.edge.label=cex,
       col.node.label=tip.color, col.edge.label=tip.color,
       font.node.label=font, font.edge.label=font,
       add=FALSE, ...){
   edge = net$edge
   label = net$tip.label
   xx = coords[,1]
   yy = coords[,2]
   nTips = length(label)

#   cex=1
   
   xlim <- range(xx)
   ylim <- range(yy)
     
   if(show.tip.label){
       offset <- max(nchar(label)) * 0.018 * cex * diff(xlim)
       xlim = c(xlim[1]-offset, xlim[2]+offset)
       ylim = c(ylim[1]-0.03 * cex * diff(ylim), ylim[2]+0.03 * cex * diff(ylim))
   }
   if(!add){ 
       plot.new() 
       plot.window(xlim, ylim, asp=1)
   }
   cladogram.plot(edge, xx, yy, edge.color, edge.width, edge.lty)
   if(show.tip.label){
        ind=match(1:nTips, edge[,2])
        pos = rep(4, nTips)
        XX <- xx[edge[ind, 1]] - xx[edge[ind, 2]]
        pos[XX>0] = 2
        YY <- yy[edge[ind, 1]] - yy[edge[ind, 2]]
        pos2 <- rep(3, nTips)
        pos2[YY>0] = 1
# needed if tiplabels are not at internal nodes        
        XX[is.na(XX)] = 0
        YY[is.na(YY)] = 0
        pos[abs(YY)>abs(XX)] <- pos2[abs(YY)>abs(XX)] 	
        text(xx[1:nTips], yy[1:nTips], labels=label, pos=pos, col=tip.color, cex=cex, font=font)
    }
    if(show.edge.label){
	    ec = edgeLabels(xx,yy, edge=edge)
	    if(is.null(edge.label))edge.label = net$splitIndex
	    
	    # show only one edge label
	    em <- apply(ec, 1, function(x)max(abs(x)))
	    si <- net$splitIndex
	    for(i in unique(si)){
	        tmp <- si==i
	        if(sum(tmp)>1){
	            w <- which(tmp)
	            wm <- which.max(em[w])
	            edge.label[w[-wm]] <- ""
	        }
	    }
	    
	    text(ec[,1], ec[,2], labels=edge.label, col=col.edge.label, cex=cex.edge.label, font=font.edge.label)     
	} 
    if(show.node.label){
         text(xx, yy, labels=node.label, col=col.node.label, cex=cex.node.label, font=font.node.label)    
    }   
}   
   
## as.splits.phylo    
lento <- function (obj, xlim = NULL, ylim = NULL, main = "Lento plot", 
    sub = NULL, xlab = NULL, ylab = NULL, bipart=TRUE, trivial=FALSE, col = rgb(0,0,0,.5), ...) 
{
    if (inherits(obj,"phylo")){ 
        if(inherits(obj,"phylo",TRUE)==1)  obj <- as.splits(obj)[obj$edge[,2]]
        obj <- as.splits(obj)
    }
    if (inherits(obj,"multiPhylo")) 
        obj = as.splits(obj)    
    labels = attr(obj, "labels") 
    l = length(labels)
    if(!trivial){
        triv = sapply(obj, length)
        ind = logical(length(obj)) 
        ind[(triv >1) & (triv < (l-1))] = TRUE
        if(length(col)==length(obj)) col=col[ind] 
        obj = obj[ind]
        }
    CM = compatible(obj)
    support = attr(obj, "weights")
    if (is.null(support)) 
        support = rep(1, length(obj))
    conflict = -as.matrix(CM) %*% support
    n = length(support)
    if (is.null(ylim)) {
        eps = (max(support) - min(conflict)) * 0.05
        ylim = c(min(conflict) - eps, max(support) + eps)
    }
    if (is.null(xlim)) {
        xlim = c(0, n + 1)
    }

    ord = order(support, decreasing = TRUE)
    support = support[ord]
    conflict = conflict[ord]
    if(length(col)==length(obj)) col=col[ord]
    plot.new()
    plot.window(xlim, ylim)
    title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
    segments(0:(n - 1), support, y1 = conflict, ...)
    segments(1:n, support, y1 = conflict, ...)
    segments(0:(n - 1), support, x1 = 1:n, ...)
    segments(0:(n - 1), conflict, x1 = 1:n, ...)
    abline(h = 0)
    axis(2, ...)
    aty = diff(ylim)/(l+1)
    at = min(ylim) + (1:l) * aty
    if(bipart){
        Y = rep(at, n)
        X = rep((1:n)-.5, each=l)
        Circles = matrix(1, l, n)
        for(i in 1:n) Circles[obj[[ord[i]]],i] = 19   
#    axis(4, labels=labels, at=at)
        col = rep(col, each=l)
        text(x=n+.1,y=at, labels, pos=4, ...) 
        points(X,Y,pch = as.numeric(Circles), col = col, ...)
        }
    invisible(list(support = cbind(support, conflict), splits=obj[ord]))
    }

    
write.splits = function (x, file = "", zero.print = ".", one.print = "|", print.labels = TRUE, ...) 
{
    labels = attr(x, "labels")
    x.orig <- x
    cx <- as.matrix(x, zero.print = zero.print, one.print = one.print)
    w = FALSE
    if (!is.null(attr(x, "names"))) {
        nam = TRUE
        vnames = format(attr(x, "names"))
    }
    nam = FALSE
    if (!is.null(attr(x, "weights"))) {
        w = TRUE
        weight = format(attr(x, "weights"))
    }
    d = FALSE
    if (!is.null(attr(x, "data"))) {
        d = TRUE
        data = attr(x, "data")
    }
    if(print.labels){for(i in 1:length(labels)) cat(labels[i], "\n", file = file, append = TRUE)}
    if (w) 
        cat("weight", "\t", file = file, append = TRUE)
    if (d) 
        cat(paste(colnames(data), "\t"), file = file, append = TRUE)
    cat("\n", file = file, append = TRUE) #"Matrix", 
    for (i in 1:length(x)) {
        if (nam) 
            cat(vnames[i], "\t", file = file, append = TRUE)
        if (d) 
            cat(paste(data[i, ], "\t"), file = file, append = TRUE)
        if (w) 
            cat(weight[i], "\t", file = file)
        cat("\n", paste(cx[i, ], collapse = ""),"\n",  file = file, append = TRUE)
    }
}
 

write.nexus.splits <- function (obj, file = "", weights=NULL, taxa=TRUE, append=FALSE) 
{
    taxa.labels <- attr(obj, "labels")
    ntaxa <- length(taxa.labels)
    obj <- oneWise(obj, ntaxa)
    ind <- which(sapply(obj, length)==ntaxa)
    if(length(ind))obj <- obj[-ind] 
    nsplits <- length(obj)    
    if(is.null(weights))weight <- attr(obj, "weights") 
    if (is.null(weight)) 
        weight = numeric(nsplits) + 100
    if(!append){cat("#NEXUS\n\n", file = file)
    cat("[Splits block for Spectronet or Splitstree]\n", file = file, append = TRUE)
    cat("[generated by phangorn:\n", file = file, append = TRUE)
    cat(format(citation("phangorn"), "text"), "]\n\n",
       file = file, append = TRUE)}
# TAXON BLOCK    
    if(taxa){
    cat(paste("BEGIN TAXA;\n\tDIMENSIONS ntax=", ntaxa, ";\n", 
        sep = ""), file = file, append = TRUE)
    cat("\tTAXLABELS", paste(taxa.labels, sep = " "), ";\nEND;\n\n", 
        file = file, append = TRUE)
    }
# SPLITS BLOCK      
    cat(paste("BEGIN SPLITS;\n\tDIMENSIONS ntax=", ntaxa, " nsplits=", nsplits,
        ";\n", sep = ""), file = file, append = TRUE)     
    format = "\tFORMAT labels=left weights=yes"
    fcon = fint = flab = FALSE
    if(!is.null(attr(obj, "confidences"))){ 
        format = paste(format, "confidences=yes")
        fcon=TRUE
        conf = attr(obj, "confidences")
        if(storage.mode(conf) == "character"){ 
            conf[conf==""] = "0"
            attr(obj, "confidences") = conf
        }                                       
    }
    else format = paste(format, "confidences=no") 
    if(!is.null(attr(obj, "intervals"))){ 
        format = paste(format, "intervals=yes")
        fint=TRUE
    }
    else format = paste(format, "intervals=no") 
    if(!is.null(attr(obj, "splitlabels"))) flab=TRUE
    format = paste(format, ";\n",  sep = "")
    cat(format, file = file, append = TRUE)
    if(!is.null(attr(obj, "cycle"))){
        cycle <- paste(attr(obj, "cycle"), collapse = " ")
        cat("\tCYCLE\t", cycle, ";\n", sep="", file = file, append = TRUE)
    }
    cat("\tMATRIX\n", file = file, append = TRUE)    
    
    for (i in 1:nsplits){
        slab <- ifelse(flab, attr(obj, "splitlabels")[i], i)
        scon <- ifelse(fcon, paste(attr(obj, "confidences")[i], "\t"), "")
        sint <- ifelse(fint, paste(attr(obj, "intervals")[i], "\t"), "")
        cat("\t\t", slab, "\t", weight[i], "\t", scon, sint, paste(obj[[i]], collapse=" "), 
            ",\n", file = file, append = TRUE, sep = "")  
    }
    cat("\t;\nEND;\n", file = file, append = TRUE)
}


write.nexus.networx <- function(obj, file = "", taxa=TRUE, splits=TRUE, append=FALSE){
    if(!append){
        cat("#NEXUS\n\n", file = file)
        cat("[Splits block for Spectronet or Splitstree]\n", file = file, append = TRUE)
        cat("[generated by phangorn:\n", file = file, append = TRUE)
        cat(format(citation("phangorn"), "text"), "]\n\n",
            file = file, append = TRUE)
        }
    ntaxa <- length(obj$tip.label)
    # TAXON BLOCK    
    if(taxa){
        cat(paste("BEGIN TAXA;\n\tDIMENSIONS NTAX=", ntaxa, ";\n", 
            sep = ""), file = file, append = TRUE)
        if(splits)taxalabel <- attr(obj$splits, "labels")
        else taxalabel <- obj$tip.label
        cat("\tTAXLABELS", paste(taxalabel, sep = " "), ";\nEND;\n\n", 
            file = file, append = TRUE)
    }
# SPLITS BLOCK    
    spl <- obj$splits
    if(splits){
#        spl <- changeOrder(spl, obj$tip.label) # orderSplitLabel
        write.nexus.splits(spl, file = file, weights=NULL, append = TRUE, taxa=FALSE) 
    }
    nvertices <- max(obj$edge)
    
    #    if(is.null(attr(obj, "coords")))   
    if(is.null(obj$.plot$vertices)) vertices <- coords(obj, "2D")
    else vertices <- obj$.plot$vertices
    
    # y-axis differs between R and SplitsTree
    vertices[,2] <- -vertices[,2]
    
    if(is.null(obj$.plot)) edge.col <- obj$.plot$edge.color 
    else edge.col=NULL 
    nedges <- nrow(obj$edge)
# NETWORK BLOCK
    cat(paste("BEGIN NETWORK;\nDIMENSIONS ntax=", ntaxa,
              "\tnvertices=", nvertices, "\tnedges=", nedges,";\n", sep = ""), file = file, append = TRUE)  
    cat("DRAW to_scale;\n", file = file, append = TRUE)
    cat("TRANSLATE\n", file = file, append = TRUE)
    if(is.null(obj$translate)){
        for(i in 1:ntaxa){
            cat(i, " ", obj$tip.label[i], ",\n", sep="", file = file, append = TRUE)
        }
    }
    else {
        translate <- obj$translate
        for(i in 1:nrow(translate)){
            cat(translate[i,1], " ", translate[i,2], ",\n", sep="", file = file, append = TRUE)
        }        
    }
    cat(";\nVERTICES\n", file = file, append = TRUE)
    for(i in 1:nvertices){
        cat(i, "\t", vertices[i,1], "\t", vertices[i,2], ",\n", sep="", file = file, append = TRUE)
    }
    if(!is.null(obj$tip.label)){
    cat(";\nVLABELS\n", file = file, append = TRUE)
    if(is.null(obj$translate)){    
        for(i in 1:ntaxa){
            cat(i, "\t", obj$tip.label[i], ",\n", sep="", file = file, append = TRUE)
        }
    }
    else{
        for(i in 1:nrow(translate)){
            cat(translate[i,1], " ", translate[i,2], ",\n", sep="", file = file, append = TRUE)
        }           
    }     
    }    
# cnet$splitIndex if splits = TRUE    
    cat(";\nEDGES\n", file = file, append = TRUE)
    
    if(is.null(obj$.plot$edge.color)) edge.col="black"
    else edge.col <- obj$.plot$edge.color
    if(length(edge.col)<nedges) edge.col <- rep(edge.col, length=nedges) 
    
    splI <- TRUE
    if(is.null(obj$splitIndex))splI <- FALSE
    for(i in 1:nedges){
        ecoli = edge.col[i]
        spInd <- ifelse(splI, paste("\ts=", obj$splitIndex[i], sep=""), "")
        edgeCol <- ifelse(ecoli=="black", "", paste("\tfg=", paste(col2rgb(ecoli), collapse=" "), sep=""))
        cat(i, "\t", obj$edge[i,1], "\t", obj$edge[i,2], spInd, edgeCol, ",\n", sep="", file = file, append = TRUE)
        }
    cat(";\n", file = file, append = TRUE)
    cat("END;\n", file = file, append = TRUE)
# force SplitsTree to accept the file    
    cat("\nBEGIN st_Assumptions;\n    uptodate;\nEND; [st_Assumptions]\n", file = file, append = TRUE)
}


read.nexus.networx <- function(file, splits=TRUE){
    spl <- NULL
    if(splits)spl <-read.nexus.splits(file)
    
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    semico <- grep(";", X)
    X=gsub("\\[(.*?)\\]", "", X) # get rid of comments
    
    netStart <- grep("BEGIN NETWORK;", X, ignore.case = TRUE)
    netEnd <- grep("END;", X, ignore.case = TRUE)
    netEnd <- netEnd[netEnd>netStart][1]
    dims <- grep("DIMENSION", X, ignore.case = TRUE)
    dims <- dims[(dims>netStart) & (dims<netEnd)]
    
    ntaxa = 0
    nvertices = 0 
    nedges = 0
    
    if(length(dims)>0){
        tmp = X[dims]    
        tmp = gsub("\\s+", "", tmp)
        
        ntaxa <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)", 
            "\\3", tmp, perl = TRUE, ignore.case = TRUE))
        nvertices  <- as.numeric(sub("(.+?)(nvertices\\s*\\=\\s*)(\\d+)(.+)", 
            "\\3", tmp, perl = TRUE, ignore.case = TRUE))
        nedges <- as.numeric(sub("(.+?)(nedges\\s*\\=\\s*)(\\d+)(.+)", 
            "\\3", tmp, perl = TRUE, ignore.case = TRUE))
    }
    transl <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- if (length(transl) == 1 && transl > netStart) TRUE
    else FALSE
    if (translation) {
        end <- semico[semico > transl][1]
        x <- X[(transl + 1):end]
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
#        n <- dim(TRANS)[1]
    }
    
    
    vert <- grep("VERTICES", X, ignore.case = TRUE)
    start <- vert[vert>max(dims, netStart)][1] + 1
    end <- semico[semico>start][1] -1
    VERT <- matrix(0, nvertices, 3, dimnames = list(NULL, c("id", "x", "y")))
    j=1
    for(i in start:end){
        tmp <- X[i]
#        tmp <- sub("\\s+", "", tmp) 
        tmp <- gsub("\\,", "", tmp)  
        tmp <- strsplit(tmp, "[[:space:]]")[[1]]
        VERT[j,1] <- as.numeric(tmp[1]) 
        VERT[j,2] <- as.numeric(tmp[2])
        VERT[j,3] <- as.numeric(tmp[3])
        j=j+1
    }
    
    edges <- grep("EDGES", X, ignore.case = TRUE)
    start <- edges[edges>max(dims, netStart)][1] + 1
    end <- semico[semico>start][1] -1
    EDGE <- NULL
    if(splits) EDGE <- matrix(0, nedges, 4, dimnames = list(NULL, c("id", "vert_id_2", "vert_id_2", "splits_id")))
    else EDGE <- matrix(0, nedges, 3, dimnames = list(NULL, c("id", "vert_id_2", "vert_id_2")))
    j=1
    for(i in start:end){
        tmp <- X[i]
        tmp <- gsub("\\,", "", tmp)
        #        tmp <- sub("\\s+", "", tmp) 
        tmp <- strsplit(tmp, "[[:space:]]")[[1]]
        EDGE[j,1] <- as.numeric(tmp[1]) 
        EDGE[j,2] <- as.numeric(tmp[2])
        EDGE[j,3] <- as.numeric(tmp[3])
        if(splits){
            EDGE[j,4] <- as.numeric(sub("s=", "", tmp[4], ignore.case = TRUE))
        }    
        j=j+1
    }
    
    swapEdge <- function(x, old, new) {
        x[x==new] <- -1L
        x[x==old] <- new
        x[x==-1L] <- old
        x     
    }
    swapRow <- function(x, old, new) {
        tmp <- x[old,]
        x[old,] <- x[new,]
        x[new,] <- tmp
        x     
    }
    splitIndex <- if(ncol(EDGE)==4) EDGE[,4]
    else NULL
# quick and dirty   
    el = sqrt(rowSums((VERT[EDGE[,2],c(2:3)] - VERT[EDGE[,3],c(2:3)])^2))
    edge <- EDGE[,c(2:3)]
    vert <- VERT[,c(2:3)]
    oldLabel <- as.integer(as.numeric(TRANS[,1]))
    for(i in 1:nrow(TRANS)){
        edge <- swapEdge(edge, oldLabel[i], i) 
        vert <- swapRow(vert, oldLabel[i], i)
    }
    # y-axis differs between in R and SplitsTree
    vert[,2] <- -vert[,2]  
    translate=data.frame(as.numeric(TRANS[,1]), TRANS[,2], stringsAsFactors=FALSE)
    plot <- list(vertices=vert)        
    obj <- list(edge=edge, tip.label=TRANS[,2], Nnode=max(edge)-ntaxa,
        edge.length=el, splitIndex=splitIndex, splits=spl ) 
    obj$.plot <- list(vertices = vert, edge.color="black", edge.width=3, edge.lty = 1)
    class(obj) <- c("networx", "phylo")
#    list(ntaxa=ntaxa, nvertices=nvertices, nedges=nedges, translate=TRANS, vertices=VERT, edges=EDGE, splits=spl)
    reorder(obj)
    obj
}


read.nexus.splits <- function(file)
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    semico <- grep(";", X)
    X=gsub("\\[(.*?)\\]", "", X) # get rid of comments
    i1 <- grep("TAXLABELS", X, ignore.case = TRUE)    
    taxlab <- ifelse(length(i1)>0, TRUE, FALSE) 
    if (taxlab) {
        end <- semico[semico >= i1][1]
        x <- X[(i1):end] # assumes not a 'new line' after "TRANSLATE"
        x <- gsub("TAXLABELS", "", x, ignore.case = TRUE)
        x <- unlist(strsplit(x, "[,; \t]"))   
        x <- x[nzchar(x)]
        x <- gsub("['\"]", "", x)
        xntaxa <- length(x)
    }
    sp <- grep("SPLITS;", X, ignore.case = TRUE)
    spEnd <- grep("END;", X, ignore.case = TRUE)
    spEnd <- spEnd[spEnd>sp][1]
    dims <- grep("DIMENSION", X, ignore.case = TRUE)
    cyc <- grep("CYCLE", X, ignore.case = TRUE)
    matr <- grep("MATRIX", X, ignore.case = TRUE)
    format <- grep("FORMAT", X, ignore.case = TRUE)
    start <- matr[matr>sp][1] + 1
    end <- semico[semico>start][1] -1
    format <- format[(format>sp) & (format<spEnd)]
    
    res <- vector("list", end - start + 1)
    weights = numeric(end - start + 1)
    j=1
    
    flab = fwei = fcon = fint = FALSE
    
    if(length(format)>0){
        tmp = X[format]    
        tmp = gsub("\\;", "", tmp)
        tmp = gsub("\\s+", "", tmp)
        flab = grepl("labels=left", tmp, ignore.case = TRUE) 
        fwei = grepl("weights=yes", tmp, ignore.case = TRUE) 
        fcon = grepl("confidences=yes", tmp, ignore.case = TRUE) 
        fint = grepl("intervals=yes", tmp, ignore.case = TRUE) 
        # = as.numeric(na.omit(as.numeric(strsplit(tmp, " ")[[1]])))        
        ind = cumsum(c(flab, fwei, fcon, fint))
        mformat = sum(c(flab, fwei, fcon, fint))
     }
    
    if(fint)intervals = numeric(end - start + 1)
    if(fcon)confidences = numeric(end - start + 1)
    if(flab)labels = vector("character", end - start + 1)
   
    for(i in start:end){
        tmp = X[i]
        tmp = sub("\\s+", "", tmp) 
        tmp = strsplit(tmp, "\t")[[1]]
        if(flab)labels[j] = as.numeric(tmp[ind[1]])        
        if(fwei)weights[j] = as.numeric(tmp[ind[2]])
        if(fcon)confidences[j] = as.numeric(tmp[ind[3]])
        if(fint)intervals[j] = as.numeric(tmp[ind[4]])
        tmp = tmp[length(tmp)]
        tmp = gsub("\\,", "", tmp)
        res[[j]] = as.integer(na.omit(as.numeric(strsplit(tmp, " ")[[1]])))
        j=j+1
    }
    if(length(cyc)>0){
        tmp = X[cyc]    
        tmp = gsub("\\;", "", tmp)
        tmp = gsub("CYCLE", "", tmp, ignore.case = TRUE)
        tmp = sub("\\s+", "", tmp)
        cyc = as.integer(na.omit(as.numeric(strsplit(tmp, " ")[[1]])))
    }
    attr(res, "labels") = x
    attr(res, "weights") = weights
    if(fint)attr(res, "intervals") = intervals
    if(fcon)attr(res, "confidences") = confidences
    if(flab)attr(res, "splitlabels") = labels
    attr(res, "cycle") = cyc 
    class(res) = c("splits", "prop.part")
    res
}



################################################################################
# delta.score
################################################################################
# Calculated from mathematical description given in Gray et al. (2010) Phil.
# Trans. Roy. Soc. B. 
# delta.score reference: Holland et al. (2002) Mol. Biol. Evol.
################################################################################ 


# Calculating Delta and Q-residual scores 
# internal
delta.quartet <-
    function(quartet,dist.dna) {
        m1 <- dist.dna[quartet[1],quartet[2]] + dist.dna[quartet[3],quartet[4]]
        m2 <- dist.dna[quartet[1],quartet[3]] + dist.dna[quartet[2],quartet[4]]
        m3 <- dist.dna[quartet[1],quartet[4]] + dist.dna[quartet[2],quartet[3]]
        m <- sort(c(m1,m2,m3),decreasing=T)
        if((m[1]-m[3])!=0) {
            ret <- (m[1]-m[2])/(m[1]-m[3])
        } else {
            ret <- 0
        }
        return(ret)
    }


delta.score <- function(x, arg="mean", ...) {
        # dist.dna <- as.matrix(dist.dna(dna,"raw"))   
        # dist.dna(dna,"raw") is equivalent to dist.hamming(as.phyDat(dna), exclude="all") 
        dist.dna <- as.matrix(dist.hamming(x, ...))
        # Number of quartets
        # choose(length(names(x)),4)
        # Create all quartets
        all.quartets <- t(combn(names(x),4))
        delta.values <- apply(all.quartets[,],1,delta.quartet,dist.dna)
        if (!arg%in%c("all", "mean","sd")) stop("return options are: all, mean, or sd")
        if (arg=='all') return(delta.values)
        if (arg=='mean') return(mean(delta.values))
        if (arg=='sd') return(sd(delta.values))
    }

