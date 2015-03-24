#
# Maximum Parsimony 
#
rowMin = function(X){
    d=dim(X)
    .Call("C_rowMin", X, as.integer(d[1]), as.integer(d[2]), PACKAGE = "phangorn") 
}


sankoff.quartet <- function (dat, cost, p, l, weight) 
{
    erg <- .Call("sankoffQuartet", sdat = dat, sn = p, scost = cost, 
        sk = l, PACKAGE = "phangorn")
    sum(weight * erg)
}


parsimony <- function(tree, data, method='fitch', ...){
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
    if(method=='sankoff') result <- sankoff(tree, data, ...)
    if(method=='fitch') result <- fitch(tree, data, ...)
    result 
}


ancestral.pars <- function (tree, data, type = c("MPR", "ACCTRAN"), cost=NULL) 
{
    call <- match.call()
    type <- match.arg(type)
    if (type == "ACCTRAN") 
        res = ptree(tree, data, retData = TRUE)[[2]]
#    if (type == "MPR") 
#        res = mpr(tree, data)
    if (type == "MPR"){ 
        res <- mpr(tree, data, cost=cost)
        attr(res, "call") = call
        return(res)
    }
    l = length(tree$tip)
    
    x = attributes(data)
    m = dim(res)[2]
    label = as.character(1:m)
    nam = tree$tip.label
    label[1:length(nam)] = nam
    x[["names"]] = label

    nc = attr(data, "nc")
    result = vector("list", m) 
    Z = unique(as.vector(res))
    tmp = t(sapply(Z, function(x)dec2bin(x, nc)))
    tmp = tmp / rowSums(tmp)
#    rownames(tmp) = Z
    dimnames(tmp) = list(Z, attr(data, "levels"))
    for(i in 1:m){ 
#        tmp = t(sapply(res[,i], function(x, k=4)dec2bin(x, nc)))
#        result[[i]] = tmp / rowSums(tmp) no indices
#         test = match(res[,i], Z) sollte stimmen wegen fitch
         result[[i]] = tmp[as.character(res[,i]),,drop=FALSE]
         rownames(result[[i]]) = NULL
        }

    attributes(result) = x
    attr(result, "call") <- call
    result
}


pace <- ancestral.pars


mpr.help = function (tree, data, cost=NULL) 
{   
    tree<- reorder(tree, "postorder")     
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")    
    levels <- attr(data, "levels")
    l = length(levels)
    if (is.null(cost)) {
        cost <- matrix(1, l, l)
        cost <- cost - diag(l)
    }   
    weight = attr(data, "weight")
    p = attr(data, "nr")
    kl = TRUE
    i = 1
    dat <- prepareDataSankoff(data)
    for (i in 1:length(dat)) storage.mode(dat[[i]]) = "double"    
    tmp = fit.sankoff(tree, dat, cost, returnData='data')
    p0 = tmp[[1]]    
    datf = tmp[[2]]
    datp = pnodes(tree, datf, cost) 

    nr = attr(data, "nr")
    nc = attr(data, "nc")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]

    node = as.integer(node - 1)      
    edge = as.integer(edge - 1) 

    res <- .Call("sankoffMPR", datf, datp, as.numeric(cost), as.integer(nr),as.integer(nc),
                 node, edge, PACKAGE="phangorn")    
    root = getRoot(tree)
    res[[root]] <- datf[[root]]
    res
}


mpr <- function(tree, data, cost=NULL){
    data = subset(data, tree$tip.label)
    att = attributes(data)
    nr = att$nr
    nc = att$nc
    res <- mpr.help(tree,data,cost)
    l = length(tree$tip)
    m = length(res)
    label = as.character(1:m)
    nam = tree$tip.label
    label[1:length(nam)] = nam
    att[["names"]] = label
    ntips = length(tree$tip)
    contrast = att$contrast
    eps=5e-6
    rm = apply(res[[ntips+1]], 1, min)
    RM = matrix(rm,nr, nc) + eps
    for(i in 1:ntips) res[[i]] = contrast[data[[i]],,drop=FALSE]
    for(i in (ntips+1):m) res[[i]][] = as.numeric(res[[i]] < RM)
    fun = function(X){
        rs = apply(X, 1, sum)
        X / rs
    }
    res <- lapply(res, fun)
    attributes(res) = att
    res
}



plotAnc <- function (tree, data, i = 1, col=NULL, cex.pie=par("cex"), pos="bottomright", ...)
{
   y = subset(data, , i)
#   args <- list(...)
#   CEX <- if ("cex" %in% names(args))
#       args$cex 
#   else par("cex")
   CEX = cex.pie
   xrad <- CEX * diff(par("usr")[1:2])/50
   levels = attr(data, "levels")
   nc = attr(data, "nc")
   y = matrix(unlist(y[]), ncol = nc, byrow = TRUE)
   l = dim(y)[1]
   dat = matrix(0, l, nc)
   for (i in 1:l) dat[i, ] = y[[i]]
   plot(tree, label.offset = 1.1 * xrad, plot = FALSE, ...)
   lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
   XX <- lastPP$xx
   YY <- lastPP$yy
   xrad <- CEX * diff(lastPP$x.lim * 1.1)/50
   par(new = TRUE)
   plot(tree, label.offset = 1.1 * xrad, plot = TRUE, ...)
   if(is.null(col)) col = rainbow(nc)
   if(length(col)!=nc) warning("Length of color vector differs from number of levels!")
   BOTHlabels(pie = y, XX = XX, YY = YY, adj = c(0.5, 0.5),
       frame = "rect", pch = NULL, sel = 1:length(XX), thermo = NULL,
       piecol = col, col = "black", bg = "lightblue", horiz = FALSE,
       width = NULL, height = NULL, cex=cex.pie)
   legend(pos, levels, text.col = col)
}


prepareDataFitch <- function (data) 
{
    lev <- attr(data, "levels")
    l <- length(lev)
    nr <- attr(data, "nr")
    nc <- length(data)
    contrast <- attr(data, "contrast")
    tmp = contrast %*% 2L^c(0L:(l - 1L))
    tmp = as.integer(tmp)
    attrData <- attributes(data)
    nam <- attrData$names
    attrData$names <- NULL
    data = unlist(data, FALSE, FALSE)
    X = tmp[data]  
    attributes(X) <- attrData
    attr(X, "dim") <- c(nr, nc)
    dimnames(X) <- list(NULL, nam)
    X
}


compressSites <- function(data){
    attrData <- attributes(data)
    lev <- attr(data, "levels")  
    LEV <- attr(data,"allLevels")
    l <- length(lev)
    nr <- attr(data, "nr")
    nc <- length(data)

    data = unlist(data, FALSE, FALSE)

    attr(data, "dim") <- c(nr, nc)
    uni <- match(lev, LEV)    
    fun = function(x, uni) {
        u = unique.default(x)
        res=  if(any(is.na(match(u, uni)))) return(x)
        match(x, u)
    }
    data = t(apply(data, 1, fun, uni))         
    ddd = fast.table(data)
    data = ddd$data
    class(data) = "list"
    attrData$weight = tapply(attrData$weight,ddd$index, sum)
    attrData$index=NULL
    attrData$nr <- length(attrData$weight)
    attrData$compressed <- TRUE
    attributes(data) <- attrData
    data
}



is.rooted2 = function(tree){
    length(tree$edge[, 1][!match(tree$edge[, 1], tree$edge[, 2], 0)]) < 3
}


#
# Branch and bound 
#

parsinfo <- function (x) 
{
    low = lowerBound(x)
    up = upperBound(x)
    ind = which(low==up)
    cbind(ind, low[ind])
}


lowerBound <- function(x, cost=NULL){
    tip <- names(x)
    att = attributes(x)
    nc = attr(x, "nc")
    nr = attr(x, "nr")
    contrast = attr(x, "contrast")
    rownames(contrast) = attr(x, "allLevels")
    colnames(contrast) = attr(x, "levels")
    attr(x, "weight") = rep(1, nr)
    attr(x, "index") = NULL
 
    y <- as.character(x)
    states <- apply(y, 2, unique.default)
# duplicated function(x)x[duplicated(x)]="?" avoids looping
    if(nr==1) nst <- length(states)   
    else nst <- sapply(states, length)

    res = numeric(nr)
    ust = sort(unique(nst))

    if(is.null(cost))cost <- 1 - diag(nc)



    if(any(ust>1)){ 
        ust = ust[ust>1]
        m <- max(ust)    
        tips = paste("t", 1:m, sep="") 
#         
        for(i in ust){
            dat = matrix(unlist(states[nst==i]), nrow=i, dimnames=list(tips[1:i], NULL))
            dat = phyDat(dat, type="USER", contrast=contrast)      
            tree = stree(i)
            res[nst==i] = sankoffNew(tree, dat, cost=cost, site="site")[attr(dat, "index")]
        }
    }
    res
}


upperBound <- function(x, cost=NULL){
    tree = stree(length(x), tip.label=names(x))
    if(is.null(cost))cost <- 1 - diag(attr(x, "nc")) 
    sankoffNew(tree, x, cost=cost, site="site")
}


CI <- function (tree, data, cost=NULL){
    pscore = sankoff(tree, data, cost=cost)
    weight = attr(data, "weight")
    data = subset(data, tree$tip.label) 
    m = lowerBound(data, cost=cost)    
    sum(m * weight)/pscore
}


RI <- function (tree, data, cost=NULL)
{
    pscore = sankoffNew(tree, data, cost=cost)
    data = subset(data, tree$tip.label)
    weight = attr(data, "weight")
    m = lowerBound(data, cost=cost)
    m = sum(m * weight)
    g = upperBound(data, cost=cost)
    g = sum(g * weight)
    (g - pscore)/(g - m)
}

# not used
add.one <- function (tree, tip.name, i){
    if (class(tree) != "phylo") 
        stop("tree should be an object of class 'phylo.'")
    nTips = length(tree$tip)
    tmpedge = tree$edge
    m = max(tmpedge)
    l = nrow(tmpedge)
    trees <- vector("list", l)
    tmp = tree
    tmp$tip.label = c(tree$tip.label, tip.name)
    tmpedge[tmpedge > nTips] <- tmpedge[tmpedge > nTips] + 1L
    tmp$Nnode = tmp$Nnode + 1L
    tmp$edge.length <- NULL
    tmpedge = rbind(tmpedge, matrix(c(m + 2L, m + 2L, 0L, nTips + 1L), 2, 2))
    edge = tmpedge
    edge[l + 1L, 2] <- edge[i, 2]
    edge[i, 2] <- m + 2L
    neworder = .C("C_reorder", edge[, 1], edge[, 2], as.integer(l + 
           2L), as.integer(m + 2L), integer(l + 2L), as.integer(nTips + 
           1L), PACKAGE = "phangorn")[[5]] 
    tmp$edge <- edge[neworder, ]
    tmp
}


mmsNew0 <- function (x, Y) 
{
    w <- attr(x, "weight")
    names(w) = NULL
    m = length(x)
    data <- matrix(unlist(x[1:m]), ncol = m)
    l = nrow(data)
    v = Y[,1] + 1L
#    v = numeric(l)
#    for (i in 1:l) v[i] = length(.Internal(unique(data[i, ], 
#        FALSE, FALSE)))
    result = matrix(NA, sum(w), 6)
    Res = matrix(NA, sum(w), m)
    Res2 = matrix(NA, sum(w), m)
    j = 1
    res = 0
    bin = as.integer(2L^c(0L:30L))
    for (i in 1:(l - 1L)) {
        if (w[i] > 0) {
            v2 = v[i] + v[(i + 1L):l] - 2L
            v3 = integer(l - i)
            ind = which(w[(i + 1):l] > 0)
            V3 = matrix(NA, m, l - i)
            k = length(ind)
            V3[, ind] <- t(data[ind + i, , drop = FALSE]) + 100L * 
                data[i, ]
            v3[ind] <- apply(V3[, ind, drop = FALSE], 2, function(x) {
                 length(unique.default(x, FALSE, FALSE)) - 1L })
#                length(.Internal(unique(x, FALSE, FALSE))) - 1L })
            r = v3 - v2
            while (any(r > 0) && w[i] > 0) {
                a = which.max(r)
                w0 = min(w[i], w[i + a])
                if (w0 == 0) {
                  r[a] = 0
                }
                else {
                  res = res + w0 * v3[a]
                  w[i] = w[i] - w0
                  w[i + a] = w[i + a] - w0
                  result[j, ] = c(i, a + i, w0, r[a], v3[a], 
                    v2[a])
                  abc = V3[, a]
                  Res[j, ] = bin[match(abc, unique(abc))]
                  Res2[j, ] = match(abc, unique(abc))
                  r[a] = 0
                  j = j + 1
                }
            }
        }
    }
    result = na.omit(result)
    mm = max(result[, 5])
    Res = na.omit(Res)
    Res2 = na.omit(Res2)
    maxr = max(Res2)
    resm = apply(Res2, 1, function(x) {
           length(unique.default(x, FALSE, FALSE)) - 1L })
#            length(.Internal(unique(x, FALSE, FALSE))) - 1L })
    Res2 = t(Res2)
    Res2 = phyDat(Res2, type="USER", levels=1:maxr)
    names(Res2) = as.character(1:m)
    resm = lowerBound(Res2)
    ind = which(w > 0)
#    data = data[ind, ]

    tmp = matrix(0, attr(Res2, "nr"), m)
    for (i in 4:m) {
        tmp[, i] = resm - upperBound(subset(Res2, 1:i))
    }
    tmp = tmp[attr(Res2, "index"), , drop=FALSE]
    tmp2 = Y[result[,1],] + Y[result[,2],]
    tmp3 = pmax(tmp, tmp2) 
#    Res = rbind(Res, data[ind, ])
    tmp = rbind(tmp3, Y[ind, ])
    weight = c(result[, 3], w[ind])
    res = t(tmp) %*% weight
    #res[m] - res
    res
}


#
# Sankoff 
#

#old2new.phyDat <- function(data){}
# works only for nucleotides
old2new.phyDat <- function(obj){
    att <- attributes(obj)
    l = length(obj)
    contrast <- attr(obj, "contrast")
    nr <- attr(obj, "nr")
    X = matrix(rep(rowSums(contrast), each=nr),nrow=nr)    
    res <- vector("list", l)
    for(i in 1:l){
        browser()
        tmp = X - tcrossprod(obj[[i]], contrast)
        res[[i]] = unlist(apply(tmp, 1, function(x)which(x<1e-6)[1]))
    }
    attributes(res) <- att
    res
}

old2new.phyDat <- function(obj){
    att <- attributes(obj)
    l = length(obj)
    contrast <- attr(obj, "contrast")
    nr <- attr(obj, "nr")
    X = matrix(rep(rowSums(contrast), each=nr),nrow=nr)   
    for(i in 1:l)obj[[i]][obj[[i]]>0] = 1
    res <- vector("list", l)
    contrast[contrast==0]=1e6   
    for(i in 1:l){
        tmp =  tcrossprod(obj[[i]], contrast) - X
        res[[i]] = unlist(apply(tmp, 1, function(x)which(x<1e-6)[1]))
    }
    attributes(res) <- att
    res
}



new2old.phyDat <- function(data){
    contrast = attr(data, "contrast")
    for(i in 1:length(data)) data[[i]] = contrast[data[[i]],,drop=FALSE]
    data
    }


prepareDataSankoff <- function(data){
    contrast = attr(data, "contrast")
    contrast[contrast == 0] = 1e+06
    contrast[contrast == 1] <- 0
    for (i in 1:length(data)) data[[i]] = contrast[data[[i]], , drop = FALSE]
    data
}



sankoff <- function (tree, data, cost = NULL, site = 'pscore') 
{
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")
    data <- prepareDataSankoff(data)
    levels <- attr(data, "levels")
    l = length(levels)  

    if (is.null(cost)) {
        cost <- matrix(1, l, l)
        cost <- cost - diag(l)
    }   
    for (i in 1:length(data)) storage.mode(data[[i]]) = "double"
    if(class(tree)=="phylo") return(fit.sankoff(tree, data, cost, returnData =site))
    if(class(tree)=="multiPhylo"){
	    if(is.null(tree$TipLabel))tree = unclass(tree)
	    return(sapply(tree, fit.sankoff, data, cost, site))
    }    
}


fit.sankoff <- function (tree, data, cost, returnData = c("pscore", "site", "data")) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    returnData <- match.arg(returnData) 
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = attr(data, "weight")
    nr = p = attr(data, "nr")
    q = length(tree$tip.label)
    nc = l = attr(data, "nc")
    m = length(edge) + 1
    dat = vector(mode = "list", length = m)
    dat[1:q] = data[tree$tip.label]
    node = as.integer(node - 1)
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    tips = as.integer((1:length(tree$tip))-1)
    res <- .Call("sankoff3", dat, as.numeric(cost), as.integer(nr),as.integer(nc),
         node, edge, mNodes, tips, PACKAGE="phangorn")  
    root <- getRoot(tree) 
    erg <- .Call("C_rowMin", res[[root]], as.integer(nr), as.integer(nc), PACKAGE = "phangorn")
    if (returnData=='site') return(erg)
    pscore <- sum(weight * erg)
    result = pscore
    if (returnData=="data"){ 
        result <- list(pscore = pscore, dat = res)
        }
    result
}


pnodes <- function (tree, data, cost) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    nr = nrow(data[[1]])
    nc = ncol(data[[1]])
    node = as.integer(node - 1)
    edge = as.integer(edge - 1)  
    .Call("pNodes", data, as.numeric(cost), as.integer(nr),as.integer(nc),
         node, edge, PACKAGE="phangorn")
}

           
indexNNI <- function(tree){
    parent = tree$edge[, 1]
    child = tree$edge[, 2]
 
    ind = which(child %in% parent)
    Nnode = tree$Nnode
    edgeMatrix = matrix(0,(Nnode-1),5)

    pvector <- numeric(max(parent))
    pvector[child] <- parent
    tips  <- !logical(max(parent))
    tips[parent] <-  FALSE
#    cvector <- allCildren(tree)  
    cvector <- vector("list",max(parent))   
    for(i in 1:length(parent))  cvector[[parent[i]]] <- c(cvector[[parent[i]]], child[i]) 
    k=0
    for(i in ind){        
            p1 = parent[i]
            p2 = child[i]
            e34 = cvector[[p2]]
            ind1 = cvector[[p1]]
            e12 = ind1[ind1 != p2]
            if(pvector[p1])e12=c(p1,e12)
            edgeMatrix[k+1, ] = c(e12,e34,p2)
            k=k+1
    } 
# vielleicht raus
    attr(edgeMatrix, 'root') <-cvector[[min(parent)]]  
    edgeMatrix
}
                   
        
sankoff.nni = function (tree, data, cost, ...) 
{   
    if(is.rooted(tree))tree<- reorder(unroot(tree), "postorder")     
    INDEX <-  indexNNI(tree)
    rootEdges <- attr(INDEX,"root")
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")
    levels <- attr(data, "levels")
    l = length(levels)
    weight = attr(data, "weight")
    p = attr(data, "nr")
    kl = TRUE
    i = 1
    tmp = fit.sankoff(tree, data, cost, returnData='data')
    p0 = tmp[[1]]
    datf = tmp[[2]]
    datp = pnodes(tree, datf, cost) 
    
    parent = tree$edge[,1]
    child = tree$edge[,2]
    m <- dim(INDEX)[1]
    k = min(parent)
    pscore = numeric(2*m)

    for(i in 1:m){
        ei = INDEX[i,]
        datn <- datf[ei[1:4]]
        if (!(ei[5] %in% rootEdges)) datn[1] = datp[ei[1]]
        pscore[(2*i)-1] <- sankoff.quartet(datn[ c(1, 3, 2, 4)], 
            cost, p, l, weight)
        pscore[(2*i)] <- sankoff.quartet(datn[ c(1, 4, 3, 2)], 
            cost, p, l, weight)    
    }
    swap <- 0
    candidates <- pscore < p0
    while(any(candidates)){
    
        ind = which.min(pscore)
        pscore[ind]=Inf
        if( ind %% 2 ) swap.edge = c(2,3)
        else swap.edge = c(2,4)

        tree2 <- changeEdge(tree, INDEX[(ind+1)%/%2,swap.edge])
        test <- fit.sankoff(tree2, data, cost, 'pscore')

        if(test >= p0) candidates[ind] = FALSE
        if(test < p0) {
            p0 <- test
            swap=swap+1
            tree <- tree2
            candidates[ind] = FALSE
            indi <- which(rep(colSums(apply(INDEX,1,match,INDEX[(ind+1)%/%2,],nomatch=0))>0,each=2))
            candidates[indi] <- FALSE
            pscore[indi] <- Inf
        }
    }
    list(tree = tree, pscore = p0, swap = swap)
}


optim.parsimony <- function(tree,data, method='fitch', cost=NULL, trace=1, rearrangements="SPR", ...){
    if(method=='fitch') result <- optim.fitch(tree=tree, data=data, trace=trace, rearrangements=rearrangements, ...) 
    if(method=='sankoff') result <- optim.sankoff(tree=tree, data=data, cost=cost, trace=trace, ...)
    result 
}


pratchet <- function (data, start=NULL, method="fitch", maxit=100, k=5, trace=1, all=FALSE, rearrangements="SPR", ...) 
{
    eps = 1e-08
#    if(method=="fitch" && (is.null(attr(data, "compressed")) || attr(data, "compressed") == FALSE)) 
#       data <- compressSites(data)
    trace = trace - 1
    uniquetree <- function(trees) {
        k = 1
        res = trees[[1]]
        result = list()
        result[[1]]=res
        k=2
        trees = trees[-1]
        while (length(trees) > 0) {
# test and replace            
# change RF to do this faster RF.dist(res, trees) class(tree) = "multiPhylo"
#            rf2 = RF.dist(res, trees, FALSE)
            rf = sapply(trees, RF.dist, res, FALSE) 
            if(any(rf==0))trees = trees[-which(rf == 0)]
            if (length(trees) > 0) {
                res = trees[[1]]
                result[[k]] = res
                k=k+1 
                trees = trees[-1]
            }
        }
        result
    }
    if (is.null(start)) 
        start = optim.parsimony(nj(dist.hamming(data)), data, trace = trace, method=method, rearrangements=rearrangements, ...)
    tree = start
    data = subset(data, tree$tip.label) 
    attr(tree, "pscore") = parsimony(tree, data, method=method, ...)
    mp <- attr(tree, "pscore")
    if (trace >= 0) 
        print(paste("Best pscore so far:",attr(tree, "pscore")))

    FUN = function(data, tree, method, rearrangements, ...) 
         optim.parsimony(tree, data = data, method=method, rearrangements=rearrangements, ...)
    result = list()
    result[[1]] = tree
    kmax = 1
    for (i in 1:maxit) {
        bstrees <- bootstrap.phyDat(data, FUN, tree = tree, bs = 1, trace = trace, method=method, rearrangements=rearrangements, ...)
        trees <- lapply(bstrees, optim.parsimony, data, trace = trace, method=method, rearrangements=rearrangements, ...)
        if(class(result)=="phylo")m=1
        else m = length(result)
        if(m>0) trees[2 : (1+m)] = result[1:m]
        pscores <- sapply(trees, function(data) attr(data, "pscore"))
        mp1 = min(pscores)
        if((mp1+eps) < mp) kmax=1
        else kmax=kmax+1
        mp=mp1

        if (trace >= 0) 
            print(paste("Best pscore so far:",mp))
        ind = which(pscores < mp + eps)
        if (length(ind) == 1) {
            result = trees[ind]
            tree = result[[1]]
        }
        else {
            result = uniquetree(trees[ind])
            l = length(result)
            tree = result[[sample(l, 1)]]
        }
        if(kmax == k) break()
    }# for
    if(!all) return(tree)
    if(length(result)==1) return(result[[1]])
    class(result) = "multiPhylo"
    result
}  # pratchet



optim.sankoff <- function(tree, data, cost=NULL, trace=1, ...) {
    if(class(tree)!="phylo") stop("tree must be of class phylo") 
    if(is.rooted(tree))tree <- unroot(tree)
    if(is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") tree <- reorder(tree, "postorder")
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
    
    rt = FALSE
    dat <- prepareDataSankoff(data)
    l <- attr(dat, "nc")
    if (is.null(cost)) {
        cost <- matrix(1, l, l)
        cost <- cost - diag(l)
        #       rt = TRUE
    }
    tree$edge.length=NULL
    swap = 0
    iter = TRUE
    pscore <- fit.sankoff(tree,dat,cost,'pscore')
    while (iter) {
        res <- sankoff.nni(tree,dat,cost,...)
        tree <- res$tree
        if(trace>1)cat("optimize topology: ", pscore , "-->", res$pscore, 
            "\n")
        pscore = res$pscore
        swap = swap + res$swap
        if (res$swap == 0) iter = FALSE
        }
    if(trace>0)cat("Final p-score",pscore,"after ",swap, "nni operations \n") 
    if(rt)tree <- ptree(tree, data)  
    attr(tree, "pscore") = pscore
    tree
}


#
# ACCTRAN
#
ptree <- function (tree, data, type = "ACCTRAN", retData = FALSE) 
{
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "pruningwise") 
 #   if (!is.binary.tree(tree)) 
 #       stop("Tree must be binary!")
    tmp = fitch(tree, data, site = "data")
    nr = attr(data, "nr")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = attr(data, "weight")
    m = length(edge) + 1
    q = length(tree$tip)
    l = as.integer(length(edge))
    nTips = length(tree$tip)
    dat = tmp[[2]]
    if (!is.rooted2(tree)) {
        root = getRoot(tree)
        ind = edge[node == root]
        rSeq = .C("fitchTriplet", integer(nr), dat[, ind[1]], 
            dat[, ind[2]], dat[, ind[3]], as.integer(nr))
        dat[, root] = rSeq[[1]]
    }
    result <- .C("ACCTRAN2", dat, as.integer(nr), numeric(nr), 
        as.integer(node[l:1L]), as.integer(edge[l:1L]), l, as.double(weight), 
        numeric(l), as.integer(nTips))
    el = result[[8]][l:1L]
    if (!is.rooted2(tree)) {
        ind2 = which(node[] == root)
        dat = matrix(result[[1]], nr, max(node))
        result <- .C("ACCTRAN3", result[[1]], as.integer(nr), 
            numeric(nr), as.integer(node[(l - 3L):1L]), as.integer(edge[(l - 
                3L):1L]), l - 3L, as.double(weight), numeric(l), 
            as.integer(nTips))
        el = result[[8]][(l - 3L):1L]
        pars = .C("fitchTripletACC4", dat[, root], dat[, ind[1]], 
            dat[, ind[2]], dat[, ind[3]], as.integer(nr), numeric(1), 
            numeric(1), numeric(1), as.double(weight), numeric(nr), 
            integer(nr))
        el[ind2[1]] = pars[[6]]
        el[ind2[2]] = pars[[7]]
        el[ind2[3]] = pars[[8]]
    }
    else {
        result <- .C("ACCTRAN3", result[[1]], as.integer(nr), 
            numeric(nr), as.integer(node[l:1L]), as.integer(edge[l:1L]), 
            l, as.double(weight), numeric(l), as.integer(nTips))
        el = result[[8]][l:1L]
    }
    tree$edge.length = el
    if (retData) 
        return(list(tree, matrix(result[[1]], nr, max(node))))
    tree
}


acctran <- function(tree, data) ptree(tree, data, type="ACCTRAN", retData=FALSE)


parsimony.plot <- function(tree, ...){
   x = numeric(max(tree$edge))
   x[tree$edge[,2]] = tree$edge.length 
   plot(tree, ...)
   ind <- get("last_plot.phylo", envir = .PlotPhyloEnv)$edge[, 2]
   edgelabels(prettyNum(x[ind]), frame = "none")
}

