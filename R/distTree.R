#
# UPGMA, NJ, UNJ, nnls
#
"upgma" <- function(D,method="average",...){
    DD=as.dist(D)
    hc = hclust(DD,method=method,...)
    result = as.phylo(hc)
    result = reorder(result, "postorder")
    result
}


"wpgma" <- function(D,method="mcquitty",...){
    DD=as.dist(D)
    hc = hclust(DD,method=method,...)
    result = as.phylo(hc)
    result = reorder(result, "postorder")
    result
}


NJ_old <- function(x) 
{
    x = as.matrix(x)
    labels <- attr(x, "Labels")[[1]]
    edge.length = NULL
    edge = NULL
    d = as.matrix(x)
    if (is.null(labels)) 
        labels = colnames(d)
    l = dim(d)[1]
    m = l - 2
    nam = 1L:l
    k = 2L * l - 2L
    while (l > 2) {
        r = rowSums(d)/(l - 2)
        i = 0
        j = 0
        tmp <- .C("out", as.double(d), as.double(r), as.integer(l), 
                  as.integer(i), as.integer(j))
        e2 = tmp[[5]]
        e1 = tmp[[4]]
        l1 = d[e1, e2]/2 + (r[e1] - r[e2])/(2)
        l2 = d[e1, e2] - l1
        edge.length = c(l1, l2, edge.length)
        edge = rbind(c(k, nam[e2]), edge)
        edge = rbind(c(k, nam[e1]), edge)
        nam = c(nam[c(-e1, -e2)], k)
        dnew = (d[e1, ] + d[e2, ] - d[e1, e2])/2
        d = cbind(d, dnew)
        d = rbind(d, c(dnew, 0))
        d = d[-c(e1, e2), -c(e1, e2)]
        k = k - 1L
        l = l - 1L
    }
    edge.length = c(d[2, 1], edge.length)
    attr(edge.length,"names") = NULL
    result = list(edge = rbind(c(nam[2], nam[1]), edge), edge.length = edge.length,
                  tip.label = labels, Nnode = m)
    class(result) <- "phylo" 
    reorder(result, "postorder")
}


NJ <- function(x) reorder(nj(x), "postorder")


UNJ <- function(x) 
{
    x = as.matrix(x)
    labels <- attr(x, "Labels")[[1]]
    edge.length = NULL
    edge = NULL
    d = as.matrix(x)
    if (is.null(labels)) 
        labels = colnames(d)
    l = dim(d)[1]
    n = l
    nam = as.character(1:l)
    m=l-2
    nam = 1:l
    k = 2*l-2       
    w = rep(1,l)
    while (l > 2) {
        r = rowSums(d)/(l - 2)
        i = 0
        j = 0
        tmp <- .C("out", as.double(d), as.double(r), as.integer(l), as.integer(i), as.integer(j))
        e2 = tmp[[5]]
        e1 = tmp[[4]]
        l1 = d[e1, e2]/2 + sum((d[e1,-c(e1,e2)] - d[e2,-c(e1,e2)])*w[-c(e1,e2)])/(2*(n-w[e1]-w[e2]))
        l2 = d[e1, e2]/2 + sum((d[e2,-c(e1,e2)] - d[e1,-c(e1,e2)])*w[-c(e1,e2)])/(2*(n-w[e1]-w[e2]))
        edge.length = c(l1, l2, edge.length)
        edge = rbind(c(k, nam[e2]), edge)
        edge = rbind(c(k, nam[e1]), edge)
        nam = c(nam[c(-e1, -e2)], k)
        
        dnew = (w[e1]*d[e1, ] + w[e2]*d[e2, ] - w[e1]*l1 - w[e2]*l2)/(w[e1] + w[e2])
        d = cbind(d, dnew)
        d = rbind(d, c(dnew, 0))
        d = d[-c(e1, e2), -c(e1, e2)]
        w = c(w, w[e1] + w[e2])
        w = w[-c(e1, e2)]
        k = k - 1
        l = l - 1
    }
    edge.length=c(d[2,1],edge.length)
    result = list(edge = rbind(c(nam[2], nam[1]), edge), 
                  edge.length=edge.length, tip.label = labels, Nnode=m)
    class(result) <- "phylo"
    reorder(result)  
}


PNJ <- function (data) 
{
    q <- l <- r <- length(data)
    weight <- attr(data,"weight")
    
    height = NULL    
    parentNodes <- NULL
    childNodes <- NULL    
    nam <- names(data)
    tip.label <- nam
    edge = 1:q
    
    z = 0
    D = matrix(0, q, q)
    
    for (i in 1:(l - 1)) {
        for (j in (i + 1):l) {
            w = (data[[i]] * data[[j]]) %*% c(1, 1, 1, 1)
            D[i, j] = sum(weight[w==0])
        }
    }
    
    while (l > 1) {
        l = l - 1
        z = z + 1
        d = D + t(D)
        if(l>1) r = rowSums(d)/(l-1)
        if(l==1) r = rowSums(d)
        M = d - outer(r,r,"+")
        diag(M) = Inf
        
        e=which.min(M)
        e0=e%%length(r)
        e1 = ifelse(e0==0, length(r), e0)
        e2= ifelse(e0==0, e%/%length(r), e%/%length(r) + 1)
        
        ind = c(e1,e2)       
        len = d[e]/2
        nam = c(nam[-ind], as.character(-l))
        
        parentNodes = c(parentNodes,-l,-l)            
        childNodes = c(childNodes,edge[e1],edge[e2])        
        
        height = c(height, len, len)
        edge = c(edge[-ind], -l)
        w = (data[[e1]] * data[[e2]]) %*% c(1, 1, 1, 1)
        w = which(w == 0)
        newDat = data[[e1]] * data[[e2]]
        newDat[w, ] = data[[e1]][w, ] + data[[e2]][w, ]   
        data = data[-c(e1,e2)]
        data[[l]] = newDat 
        if (l > 1) {
            D = as.matrix(D[, -ind])
            D = D[-ind, ]
            dv = numeric(l - 1)
            for (i in 1:(l - 1)) {
                w = (data[[i]] * data[[l]]) %*% c(1, 1, 1, 1)
                dv[i] = sum(weight[w==0])
            }
            D = cbind(D, dv)
            D = rbind(D, 0)
        }
    }
    tree <- list(edge = cbind(as.character(parentNodes),as.character(childNodes)),tip.label=tip.label) 
    class(tree) <- "phylo"
    tree <- old2new.phylo(tree)   
    reorder(tree)    
}


#
# Distance Matrix methods
#

# as.Matrix, sparse = TRUE, 
designTree <- function(tree, method="unrooted", sparse=FALSE, ...){
    if (!is.na(pmatch(method, "all"))) 
        method <- "unrooted"
    METHOD <- c("unrooted", "rooted")
    method <- pmatch(method, METHOD)
    if (is.na(method)) stop("invalid method")
    if (method == -1) stop("ambiguous method")
    if(!is.rooted(tree) & method==2) stop("tree has to be rooted")  
    if(method==1){ X <- designUnrooted(tree,...)
                   if(sparse) X = Matrix(X)  
    }
    if(method==2) X <- designUltra(tree, sparse=sparse,...)
    X
}


# splits now work
designUnrooted = function (tree, order = NULL) 
{
    if(inherits(tree, "phylo")){ 
        if (is.rooted(tree)) 
            tree = unroot(tree)
        p = bipartition(tree)
    }
    if(inherits(tree, "splits")) p <- as.matrix(tree)
    if (!is.null(order)) 
        p = p[, order]
    
    m = dim(p)[2]
    ind = rowSums(p)
    p=p[ind!=m,]
    n = dim(p)[1]
    res = matrix(0, (m - 1) * m/2, n)
    k = 1
    for (i in 1:(m - 1)) {
        for (j in (i + 1):m) {
            res[k, ] = p[, i] != p[, j]
            k = k + 1
        }
    }
    if(inherits(tree, "phylo"))colnames(res) = paste(tree$edge[, 1], tree$edge[, 2], sep = "<->")
    res
}


designUltra <- function (tree, sparse=TRUE) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree = reorder(tree, "postorder")
    leri = allChildren(tree)
    bp = bip(tree)
    n = length(tree$tip)
    l = tree$Nnode   
    nodes = integer(l)
    k = 1L
    u=numeric( n * (n - 1)/2)
    v=numeric( n * (n - 1)/2)
    m = 1L
    for (i in 1:length(leri)) {
        if (!is.null(leri[[i]])) {
            if(length(leri[[i]])==2)ind = getIndex(bp[[leri[[i]][1] ]], bp[[leri[[i]][2] ]], n) 
            else {
                ind=NULL
                le=leri[[i]]
                nl = length(le)
                for(j in 1:(nl-1)) ind =c(ind, getIndex(bp[[le[j] ]], unlist(bp[ le[(j+1):nl] ]), n))
            }
            li = length(ind)
            v[m: (m+li-1)]=k
            u[m: (m+li-1)]=ind   
            nodes[k]=i
            m = m+li
            k = k + 1L
        }
    }
    if(sparse) X = sparseMatrix(i=u,j=v, x=2L)
    else{
        X = matrix(0L, n * (n - 1)/2, l)              
        X[cbind(u,v)]=2L
    }
    colnames(X) = nodes
    attr(X, "nodes") = nodes
    X
}


designUnrooted2 <- function (tree, sparse=TRUE) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree = reorder(tree, "postorder")
    leri = allChildren(tree)
    bp = bip(tree)
    n = length(tree$tip)
    l = tree$Nnode   
    nodes = integer(l)
    nTips = as.integer(length(tree$tip))
    k = nTips
    u=numeric( n * (n - 1)/2)
    v=numeric( n * (n - 1)/2)
    z=numeric( n * (n - 1)/2)
    y=numeric( n * (n - 1)/2)    
    p=1L
    m = 1L
    for (i in 1:length(leri)) {
        if (!is.null(leri[[i]])) {
            
            if(length(leri[[i]])==2){
                ind =  getIndex(bp[[leri[[i]][1] ]], bp[[leri[[i]][2] ]], n) 
                ytmp = rep(bp[[leri[[i]][1] ]], each = length(bp[[leri[[i]][2] ]]))
                ztmp = rep(bp[[leri[[i]][2] ]], length(bp[[leri[[i]][1] ]]))  
            }    
            else {
                #        browser()
                ind=NULL
                le=leri[[i]]
                nl = length(le)
                ytmp=NULL
                ztmp=NULL
                for(j in 1:(nl-1)){ 
                    bp1 = bp[[le[j] ]]
                    bp2 = unlist(bp[le[(j+1):nl] ])
                    ind =c(ind,  getIndex(bp1, unlist(bp2), n))
                    ytmp = c(ytmp, rep(bp1, each = length(bp2)))
                    ztmp = c(ztmp, rep(bp2, length(bp1)))                     
                }    
            }
            
            li = length(ind)
            v[m: (m+li-1)]=k
            u[m: (m+li-1)]=ind
            y[m: (m+li-1)]=ytmp
            z[m: (m+li-1)]=ztmp 
            
            nodes[p]=i
            m = m+li
            k = k + 1L
            p = p + 1L
        }
    }
    jj = c(y,z) #[ind],v)
    ii = c(u,u) #[ind],u)
    ind =  (jj < nTips)
    jj = c(jj[ind], v)
    ii = c(ii[ind], u)
    l1 = length(u)    
    l2 = sum(ind)
    x= rep(c(-1L,2L), c(l2, l1))

    X = sparseMatrix(i=ii,j=jj, x=x)
    if(!sparse){
        X = as.matrix(X)
    }
    nodes = c(1:(nTips-1L), nodes)
    colnames(X) = nodes
    attr(X, "nodes") = nodes
    X
}


nnls.tree <- function(dm, tree, rooted=FALSE, trace=1){
    if(is.rooted(tree) & rooted==FALSE){
        tree = unroot(tree)
        warning("tree was rooted, I unrooted the tree!")
    }
    tree = reorder(tree, "postorder")
    dm = as.matrix(dm)
    k = dim(dm)[1]
    labels = tree$tip
    dm = dm[labels,labels]
    y = dm[lower.tri(dm)]
    #computing the design matrix from the tree   
    if(rooted) X = designUltra(tree) 
    else X = designUnrooted2(tree)
    
    lab = attr(X, "nodes")
    
    # na.action
    if(any(is.na(y))){
        ind = which(is.na(y))
        X = X[-ind,,drop=FALSE]
        y= y[-ind]
    }
    # LS solution 
    Dmat <- crossprod(X) # cross-product computations
    dvec <- crossprod(X, y)
    betahat <- as.vector(solve(Dmat, dvec))
    betahattmp = betahat
    bhat = numeric(max(tree$edge))
    bhat[as.integer(lab)] = betahat
    betahat = bhat[tree$edge[,1]] - bhat[tree$edge[,2]]
    
    if(!any(betahat<0)){
        if(!rooted){
            RSS = sum((y-(X%*%betahattmp))^2)    
            if(trace)print(paste("RSS:", RSS))
            attr(tree, "RSS") = RSS
        }
        tree$edge.length = betahat 
        return(tree)
    }
    
    # non-negative LS
    n = dim(X)[2]
    l = nrow(tree$edge)
    Amat = matrix(0, n, l)
    
    lab = attr(X, "nodes")
    # vielleicht solve.QP.compact
    ind1 = match(tree$edge[,1], lab)
    Amat[cbind(ind1, 1:l)] = 1
    ind2 = match(tree$edge[,2], lab)
    Amat[cbind(ind2, 1:l)] = -1  
    
    betahat <- quadprog::solve.QP(as.matrix(Dmat),as.vector(dvec),Amat)$sol 
    
    # quadratic programing solving
    if(!rooted){
        RSS = sum((y-(X%*%betahat))^2) 
        if(trace)print(paste("RSS:", RSS))
        attr(tree, "RSS") = RSS
    }
    bhat = numeric(max(tree$edge))
    bhat[as.integer(lab)] = betahat
    betahat = bhat[tree$edge[,1]] - bhat[tree$edge[,2]]
    tree$edge.length = betahat
    tree
}


nnls.phylo <- function(x, dm, rooted=FALSE, trace=0){
    nnls.tree(dm, x, rooted, trace=trace)
}


nnls.splits <- function(x, dm, trace=0){
    labels=attr(x, "labels")
    dm = as.matrix(dm)
    k = dim(dm)[1]
    dm = dm[labels,labels]
    y = dm[lower.tri(dm)]
    
    x = SHORTwise(x, k)
    l <- sapply(x, length)
    if(any(l==0)) x = x[-which(l==0)]
    
    X = splits2design(x)
    
    if(any(is.na(y))){
        ind = which(is.na(y))
        X = X[-ind,,drop=FALSE]
        y= y[-ind]
    }
    X = as.matrix(X)
    n = dim(X)[2]
    int = sapply(x, length)
    Amat = diag(n) # (int)
    betahat <- nnls(X, y)  
    ind = (betahat$x > 1e-8) | int==1  
    x = x[ind]
    RSS <- betahat$deviance
    attr(x, "weights") = betahat$x[ind]
    if(trace)print(paste("RSS:", RSS))
    attr(x, "RSS") = RSS
    x
}    


nnls.splitsOld <- function(x, dm, trace=0){
    labels=attr(x, "labels")
    dm = as.matrix(dm)
    k = dim(dm)[1]
    dm = dm[labels,labels]
    y = dm[lower.tri(dm)]
    
    x = SHORTwise(x, k)
    l <- sapply(x, length)
    if(any(l==0)) x = x[-which(l==0)]
    
    X = splits2design(x)
    
    if(any(is.na(y))){
        ind = which(is.na(y))
        X = X[-ind,,drop=FALSE]
        y= y[-ind]
    }
    
    Dmat <- crossprod(X) # cross-product computations
    dvec <- crossprod(X, y)
    betahat <- as.vector(solve(Dmat, dvec))
    
    if(!any(betahat<0)){
        RSS = sum((y-(X%*%betahat))^2)    
        if(trace)print(paste("RSS:", RSS))
        attr(x, "RSS") = RSS
        attr(x, "weights") = betahat 
        return(x)
    }
    n = dim(X)[2]
    
    int = sapply(x, length)
    #    int = as.numeric(int==1)# (int>1)
    Amat = diag(n) # (int)
    betahat <- quadprog::solve.QP(as.matrix(Dmat),as.vector(dvec),Amat)$sol # quadratic programing solving
    RSS = sum((y-(X%*%betahat))^2)
    ind = (betahat > 1e-8) | int==1  
    x = x[ind]
    attr(x, "weights") = betahat[ind]
    if(trace)print(paste("RSS:", RSS))
    attr(x, "RSS") = RSS
    x
}  

nnls.networx <- function(x, dm){
    spl <- attr(x, "splits")
    spl2 <- nnls.splits(spl, dm)
    weight <- attr(spl, "weight")
    weight[] <- 0
    weight[match(spl2, spl)] = attr(spl2, "weight")
    attr(attr(x, "splits"), "weight") <- weight
    x$edge.length = weight[x$splitIndex]
    x
}


designSplits <- function (x, splits = "all", ...) 
{
    if (!is.na(pmatch(splits, "all"))) 
        splits <- "all"
    if(inherits(x, "splits")) return(designUnrooted(x))
    SPLITS <- c("all", "star") #,"caterpillar")
    splits <- pmatch(splits, SPLITS)
    if (is.na(splits)) stop("invalid splits method")
    if (splits == -1) stop("ambiguous splits method")  
    if(splits==1) X <-  designAll(x)
    if(splits==2) X <-  designStar(x)
    return(X)
}

# add return splits=FALSE
designAll <- function(n, add.split=FALSE){
    Y = matrix(0L, n*(n-1)/2, n)
    k = 1
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            Y[k,c(i,j)]=1L
            k=k+1L
        }
    }
    m <- n-1L
    X <- matrix(0L, m+1, 2^m)
    for(i in 1:m)
        X[i, ] <- rep(rep(c(0L,1L), each=2^(i-1)),2^(m-i))
    X <- X[,-1]
    if(!add.split) return((Y%*%X)%%2)
    list(X=(Y%*%X)%%2,Splits=t(X))
}


designStar = function(n){
    res=NULL
    for(i in 1:(n-1)) res = rbind(res,cbind(matrix(0,(n-i),i-1),1,diag(n-i)))
    res
}


