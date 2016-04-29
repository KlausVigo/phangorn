#
# tree distance functions
#

allKids <- function(phy){
    nTips = as.integer(length(phy$tip))
    lp=nrow(phy$edge)
    nNode = phy$Nnode
    .C("AllKids", as.integer(phy$edge[,2]), as.integer(phy$edge[,1]), as.integer(nTips), 
       as.integer(nNode), as.integer(lp), integer(lp), integer(nNode+1L),integer(nNode))
} 


coph <- function(x, path=FALSE){ 
    if (is.null(attr(x, "order")) || attr(x, "order") == "cladewise") 
        x <- reorder(x, "postorder")
    nTips = as.integer(length(x$tip.label))   
    parents = as.integer(x$edge[,1]) 
    kids = as.integer(x$edge[,2])
    lp= as.integer(length(parents))
    nNode = as.integer(x$Nnode)
    m = as.integer(max(x$edge))
    el = double(m)
    if(path) el <- rep(1.0, m)
    else el[kids] = x$edge.length
    dm <- .C("C_cophenetic", kids, parents, as.double(el), lp, m, nTips, nNode, double(nTips*(nTips-1L)/2L))[[8]]
    attr(dm, "Size") <- nTips
    attr(dm, "Labels") <- x$tip.label
    attr(dm, "Diag") <- FALSE
    attr(dm, "Upper") <- FALSE
    class(dm) <- "dist"
    dm
} 


cophenetic.splits <- function(x){
    labels <- attr(x, "labels")
    X <- splits2design(x)
    dm <- as.vector(X%*%attr(x, "weight")) 
    attr(dm, "Size") <- length(labels)
    attr(dm, "Labels") <- labels
    attr(dm, "Diag") <- FALSE
    attr(dm, "Upper") <- FALSE
    class(dm) <- "dist"
    dm
}


cophenetic.networx <- function(x){
#    spl <- attr(x, "splits")
    spl <- x$splits
    cophenetic.splits(spl)
}


SHORTwise <- function (x, nTips, delete=FALSE) 
{
    v <- 1:nTips
    l <- sapply(x, length)
    lv = floor(nTips/2)  
    for (i in 1:length(x)) { 
        if(l[i]>lv){
            y <- x[[i]]
            x[[i]] <- v[-y]
        }        
        if(l[i]==nTips/2){ 
            y <- x[[i]]
            if (y[1] != 1) 
                x[[i]] <- v[-y]
        }
    }
    if(any(l==nTips) && delete){
        x=x[l!=nTips]
    }
    x
}


oneWise <- function (x, nTips=NULL) 
{
    if(is.null(nTips))nTips <- length(x[[1L]])
    v <- 1:nTips
    for (i in 2:length(x)) {
        y <- x[[i]]
        if (y[1] != 1) 
            x[[i]] <- v[-y]
    }
    x
}


treedist <- function (tree1, tree2, check.labels=TRUE) 
{
    tree1 = unroot(tree1)
    tree2 = unroot(tree2)
    
    if (check.labels) {
        ind <- match(tree1$tip.label, tree2$tip.label)
        if (any(is.na(ind)) | length(tree1$tip.label) !=
                length(tree2$tip.label))
            stop("trees have different labels")
        tree2$tip.label <- tree2$tip.label[ind]
        ind2 <- match(1:length(ind), tree2$edge[, 2])
        tree2$edge[ind2, 2] <- order(ind)
    }
    
    tree1 = reorder(tree1, "postorder")
    tree2 = reorder(tree2, "postorder")
    
    symmetric.difference = NULL
    branch.score.difference = NULL
    path.difference = NULL
    quadratic.path.difference = NULL
    if(!is.binary.tree(tree1) | !is.binary.tree(tree2))message("Trees are not binary!")
    
    bp1 = bip(tree1)
    bp2 = bip(tree2)
    bp1 <- SHORTwise(bp1, length(tree1$tip))
    bp2 <- SHORTwise(bp2, length(tree2$tip))
    bp1 <- sapply(bp1, paste, collapse = "_")
    bp2 <- sapply(bp2, paste, collapse = "_")
    
    l = length(tree1$tip.label)

    if (!is.null(tree1$edge.length) & !is.null(tree2$edge.length)) {      
        dv1 = coph(tree1)
        dv2 = coph(tree2)
        quadratic.path.difference = sqrt(sum((dv1 - dv2)^2))
        
    }
   
    RF = sum(match(bp1, bp2, nomatch=0L)==0L) + sum(match(bp2, bp1, nomatch=0L)==0L)

    symmetric.difference = RF #2 * (p - sum(r1))
    if (!is.null(tree1$edge.length) & !is.null(tree2$edge.length)) {
        w1 = numeric(max(tree1$edge))
        w2 = numeric(max(tree2$edge))
        w1[tree1$edge[,2]] = tree1$edge.length
        w2[tree2$edge[,2]] = tree2$edge.length
        
        v1 = tree1$edge.length
        v2 = tree2$edge.length
     
        ind3 = match(bp1, bp2, nomatch=0L)
        ind4 = ind3[ind3>0]
        ind3 = which(ind3>0)

        s1 = sum((w1[ind3] - w2[ind4])^2)

        s2 = sum(w1[-ind3]^2)
        s3 = sum(w2[-ind4]^2)
        branch.score.difference = sqrt(s1 + s2 + s3)
    }
    
    tree1$edge.length = rep(1, nrow(tree1$edge))
    tree2$edge.length = rep(1, nrow(tree2$edge))

    dt1 = coph(tree1)
    dt2 = coph(tree2)  
    path.difference = sqrt(sum((dt1 - dt2)^2))
    
    result = c(symmetric.difference = symmetric.difference, 
               branch.score.difference = branch.score.difference, 
               path.difference = path.difference, 
               quadratic.path.difference = quadratic.path.difference)
    result              
}


mRF2 <- function(tree, trees, check.labels = TRUE){
    if (!inherits(trees,"multiPhylo")) 
        stop("trees should be an object of class \"multiPhylo\"")
    if (!inherits(tree,"phylo")) 
        stop("trees should be an object of class \"phylo\"")
    trees <- .compressTipLabel(trees)
    tipLabel <- attr(trees, "TipLabel")
    if (check.labels) {
        ind <- match(tipLabel, tree$tip.label)
        if (any(is.na(ind)) | length(tipLabel) != length(tree$tip.label))
            stop("trees have different labels")
        tree$tip.label <- tree$tip.label[ind]
        ind2 <- match(1:length(ind), tree$edge[, 2])
        tree$edge[ind2, 2] <- order(ind)
    }
    nTips <- length(tipLabel)
    l <- length(trees)
    RF <- numeric(l)
    trees <- .uncompressTipLabel(trees)
    #    n <- length(attr(trees, "TipLabel"))
    trees <- unclass(trees)
    if (any(sapply(trees, is.rooted))) {
        message("Some trees are rooted. Unrooting all trees.\n")
        trees <- lapply(trees, unroot)
    }
    if (any(sapply(trees, function(x) !is.binary.tree(x)))) {
        message("Some trees are not binary. Result may not what you expect!")
    }
    tree <- reorder(tree, "postorder")
    trees <- lapply(trees, reorder, "postorder")
    xx <- lapply(trees, bipart)  
    xx <- lapply(xx, SHORTwise, nTips)
    xx <- lapply(xx,function(x)sapply(x, paste, collapse="_"))
    yy <- bipart(tree)  
    yy <- SHORTwise(yy, nTips)
    yy <- sapply(yy, paste, collapse="_")
    
    NnodeT <- Nnode(tree)
    Nnodes <- sapply(trees, Nnode)
    
    for (i in 1:l){   
        RF[i] <- Nnodes[i] + NnodeT - 2 * sum(fmatch(xx[[i]], yy, nomatch=0L)>0L)
#        RF[i] <- sum(match(xx[[i]], yy, nomatch=0L)==0L) + sum(match(yy, xx[[i]], nomatch=0L)==0L)
    }
    if(!is.null(names(trees)))names(RF) <- names(trees)
    return(RF)
}


mRF<-function(trees){
    if (!inherits(trees,"multiPhylo")) 
        stop("trees should be an object of class \"multiPhylo\"")
    trees <- .compressTipLabel(trees)
    tipLabel <- attr(trees, "TipLabel")
    nTips <- length(tipLabel)
    l <- length(trees)
    RF <- numeric((l * (l - 1))/2)
    trees <- .uncompressTipLabel(trees)
    #    n <- length(attr(trees, "TipLabel"))
    trees <- unclass(trees)
    if (any(sapply(trees, is.rooted))) {
        message("Some trees are rooted. Unrooting all trees.\n")
        trees <- lapply(trees, unroot)
    }
    if (any(sapply(trees, function(x) !is.binary.tree(x)))) {
        message("Some trees are not binary. Result may not what you expect!")
    }
    trees <- lapply(trees, reorder, "postorder")
    xx <- lapply(trees, bipart)  
    xx <- lapply(xx, SHORTwise, nTips)
    xx <- lapply(xx,function(x)sapply(x, paste, collapse="_")) 
    # returns list of character vectors
    Nnodes <- sapply(trees, Nnode)
    k=1
    for (i in 1:(l - 1)){
        tmp = xx[[i]]        
        for (j in (i + 1):l){
            RF[k] <- Nnodes[i] + Nnodes[j] - 2 * sum(fmatch(xx[[j]], tmp, nomatch=0L)>0L)
#            RF[k] <- sum(match(xx[[j]], tmp, nomatch=0L)==0L) + sum(match(tmp, xx[[j]], nomatch=0L)==0L)
            k=k+1
        }   
    }
    attr(RF, "Size") <- l
    if(!is.null(names(trees)))attr(RF, "Labels") <- names(trees)
    attr(RF, "Diag") <- FALSE
    attr(RF, "Upper") <- FALSE
    class(RF) <- "dist"
    return(RF)
}


RF.dist <- function(tree1, tree2=NULL, check.labels = TRUE, rooted=FALSE)
{
    if(class(tree1)=="multiPhylo" && is.null(tree2))return(mRF(tree1)) 
    if(class(tree1)=="phylo" && class(tree2)=="multiPhylo")return(mRF2(tree1, tree2, check.labels))
    if(class(tree2)=="phylo" && class(tree1)=="multiPhylo")return(mRF2(tree2, tree1, check.labels))
    r1 = is.rooted(tree1)
    r2 = is.rooted(tree2)
    if(r1 != r2){
        message("one tree is unrooted, unrooted both")
    }  
    if(!rooted){
        if(r1) tree1<-unroot(tree1)
        if(r2) tree2<-unroot(tree2)
    }
    
    if (check.labels) {
        ind <- match(tree1$tip.label, tree2$tip.label)
        if (any(is.na(ind)) | length(tree1$tip.label) !=
                length(tree2$tip.label))
            stop("trees have different labels")
        tree2$tip.label <- tree2$tip.label[ind]
        #       tree2$edge[match(ind, tree2$edge[, 2]), 2] <- 1:length(ind)
        ind2 <- match(1:length(ind), tree2$edge[, 2])
        tree2$edge[ind2, 2] <- order(ind)
    }
    
    if(!r1 | !r2){
        if(r1) tree1 = unroot(tree1)
        if(r2) tree2 = unroot(tree2)
    }
    if(!is.binary.tree(tree1) | !is.binary.tree(tree2))message("Trees are not binary!")
    bp1 = bipart(tree1)
    bp2 = bipart(tree2)
    if(!rooted){
        bp1 <- SHORTwise(bp1, length(tree1$tip))
        bp2 <- SHORTwise(bp2, length(tree2$tip))    
    }
    RF = sum(match(bp1, bp2, nomatch=0L)==0L) + sum(match(bp2, bp1, nomatch=0L)==0L)
    RF
}


kf0 <- function(tree1, tree2, check.labels = TRUE){    
    if(check.labels)tree2 <- checkLabels(tree2, tree1$tip.label)
    if(is.rooted(tree1)) tree1 <- unroot(tree1)
    if(is.rooted(tree2)) tree2 <- unroot(tree2)
    bp1 = bip(tree1)
    bp2 = bip(tree2)
    
    bp1 <- SHORTwise(bp1, length(tree1$tip))
    bp2 <- SHORTwise(bp2, length(tree2$tip))
    bp1 <- sapply(bp1, paste, collapse = "_")
    bp2 <- sapply(bp2, paste, collapse = "_")
    
    w1 <- numeric(max(tree1$edge))
    w2 <- numeric(max(tree2$edge))
    w1[tree1$edge[,2]] <- tree1$edge.length
    w2[tree2$edge[,2]] <- tree2$edge.length
        
    ind3 = match(bp1, bp2, nomatch=0L)
    ind4 = ind3[ind3>0]
    ind3 = which(ind3>0)
        
    s1 = sum((w1[ind3] - w2[ind4])^2)
    s2 = sum(w1[-ind3]^2)
    s3 = sum(w2[-ind4]^2)
    branch.score.difference = sqrt(s1 + s2 + s3)
    branch.score.difference
}


kf1 <- function(tree, trees, check.labels = TRUE){    
    if(check.labels){
        trees <- .compressTipLabel(trees)
        tree <- checkLabels(tree, attr(trees, "TipLabel"))
    }    
    trees <- .uncompressTipLabel(trees)
    unclass(trees) 
    
    if (any(sapply(trees, is.rooted))) {
        trees <- lapply(trees, unroot)
    }
    
    nTips <- length(tree$tip.label)
    
    fun1 <- function(x){
        w <- numeric(max(x$edge))
        w[x$edge[,2]] <- x$edge.length
        w
    }
    W <- lapply(trees, fun1)
    
    fun2 <- function(x, nTips){
        bp <- bip(x)
        bp <- SHORTwise(bp, nTips)
        bp <- sapply(bp, paste, collapse = "_")
        bp
    }
    BP <- lapply(trees, fun2, nTips)
    
    
    if(is.rooted(tree)) tree <- unroot(tree)
    bp = bip(tree)
    
    bp <- SHORTwise(bp, nTips)
    bp <- sapply(bp, paste, collapse = "_")

    w <- numeric(max(tree$edge))
    w[tree$edge[,2]] <- tree$edge.length

    l <- length(trees)
    branch.score.difference <- numeric(l)
    
    for(i in 1:l){
        ind3 = fmatch(BP[[i]], bp, nomatch=0L)
        ind4 = ind3[ind3>0]
        ind3 = which(ind3>0)
    
        s1 = sum((W[[i]][ind3] - w[ind4])^2)
        s2 = sum(W[[i]][-ind3]^2)
        s3 = sum(w[-ind4]^2)
        branch.score.difference[i] = sqrt(s1 + s2 + s3)
    }
    branch.score.difference
}


kf2 <- function(trees, check.labels = TRUE){
    if(check.labels) trees <- .compressTipLabel(trees)
    trees <- .uncompressTipLabel(trees)
    unclass(trees)     
    
    nTips <- length(trees[[1]]$tip.label)
    if (any(sapply(trees, is.rooted))) {
        trees <- lapply(trees, unroot)
    }
    fun1 <- function(x){
        w <- numeric(max(x$edge))
        w[x$edge[,2]] <- x$edge.length
        w
    }
    W <- lapply(trees, fun1)
    fun2 <- function(x, nTips){
        bp <- bip(x)
        bp <- SHORTwise(bp, nTips)
        bp <- sapply(bp, paste, collapse = "_")
        bp
    }
    BP <- lapply(trees, fun2, nTips)
    k <- 1
    l <- length(trees)
    KF <- numeric((l * (l - 1))/2)
    for (i in 1:(l - 1)){
        bp <- BP[[i]]
        w <- W[[i]]
        for (j in (i + 1):l){
            ind3 = fmatch(BP[[j]], bp, nomatch=0L)
            ind4 = ind3[ind3>0]
            ind3 = which(ind3>0)
            s1 = sum((W[[j]][ind3] - w[ind4])^2)
            s2 = sum(W[[j]][-ind3]^2)
            s3 = sum(w[-ind4]^2)
            KF[k] <- sqrt(s1 + s2 + s3)
            k=k+1
        }
    }
    attr(KF, "Size") <- l
    if(!is.null(names(trees)))attr(KF, "Labels") <- names(trees)
    attr(KF, "Diag") <- FALSE
    attr(KF, "Upper") <- FALSE
    class(KF) <- "dist"
    return(KF)
}


# TODO distance matrices
KF.dist <- function(tree1, tree2=NULL, check.labels = TRUE, rooted=FALSE){
    if(inherits(tree1, "multiPhylo") && is.null(tree2))return(kf2(tree1)) 
    if(inherits(tree1, "phylo") && inherits(tree2, "phylo"))return(kf0(tree1, tree2, check.labels))
    if(inherits(tree1, "phylo") && inherits(tree2, "multiPhylo"))return(kf1(tree1, tree2, check.labels))
    if(inherits(tree2, "phylo") && inherits(tree1, "multiPhylo"))return(kf1(tree2, tree1,check.labels))
    return(NULL)
}


path.dist <- function(tree1, tree2=NULL, check.labels = TRUE, use.weight=FALSE){
    if(inherits(tree1, "phylo") && inherits(tree2, "phylo"))
         return(pd0(tree1, tree2, check.labels, !use.weight))
    if(inherits(tree1, "phylo") && inherits(tree2, "multiPhylo"))
         return(pd1(tree1, tree2, check.labels, !use.weight))
    if(inherits(tree2, "phylo") && inherits(tree1, "multiPhylo"))
        return(pd1(tree2, tree1, check.labels, !use.weight))
    if(inherits(tree1, "multiPhylo") && is.null(tree2))
        return(pd2(tree1, check.labels, !use.weight))
    else return(NULL)
}


#wPath.dist <- function(tree1, tree2=NULL, check.labels = TRUE){
#    if(inherits(tree1, "phylo") && inherits(tree2, "phylo"))
#        return(pd0(tree1, tree2, check.labels, FALSE))
#    if(inherits(tree1, "phylo") && inherits(tree2, "multiPhylo"))
#        return(pd1(tree1, tree2, check.labels, FALSE))
#    if(inherits(tree2, "phylo") && inherits(tree1, "multiPhylo"))
#        return(pd1(tree2, tree1, check.labels, FALSE))
#    if(inherits(tree1, "multiPhylo") && is.null(tree2))
#        return(pd2(tree1, check.labels, FALSE))
#    else return(NULL)
#}


pd0 <- function(tree1, tree2, check.labels=TRUE, path=TRUE){
    if(check.labels)tree2 <- checkLabels(tree2, tree1$tip.label)
    if(path){
        tree1 <- unroot(tree1)
        tree2 <- unroot(tree2)
    }    
#    n <- nrow(tree1$edge)
#    tree1$edge.length <- tree2$edge.length <- rep(1, n)
    dt1 = coph(tree1, path)
    dt2 = coph(tree2, path)  
    sqrt(sum((dt1 - dt2)^2))
}


pd1 <- function(tree, trees, check.labels=TRUE, path=TRUE){
    if(check.labels){
        trees <- .compressTipLabel(trees)
        tree <- checkLabels(tree, attr(trees, "TipLabel"))
    }    
    trees <- .uncompressTipLabel(trees)
    unclass(trees)
    if(path)trees <- lapply(trees, unroot)
    trees <- lapply(trees, reorder, "postorder")
    l <- length(trees)
    dt <- coph(tree, path)
    res <- numeric(l)
    for(i in 1:l){
        dt2 <- coph(trees[[i]], path)  
        res[i] <- sqrt(sum((dt - dt2)^2))        
    }
    res
} 

pd2 <- function(trees, check.labels=TRUE, path=TRUE){
    if(check.labels) trees <- .compressTipLabel(trees)
    trees <- .uncompressTipLabel(trees)
    unclass(trees) 
    if(path)trees <- lapply(trees, unroot)
    trees <- lapply(trees, reorder, "postorder")
    l <- length(trees)
    CM <- lapply(trees, coph, path)
    res <- numeric(l)
    k <- 1 
    PD <- numeric((l * (l - 1))/2)
    for (i in 1:(l - 1)){
        for (j in (i + 1):l){
            PD[k] <- sqrt(sum((CM[[i]] - CM[[j]])^2))    
            k=k+1
        }
    }
    attr(PD, "Size") <- l
    if(!is.null(names(trees)))attr(PD, "Labels") <- names(trees)
    attr(PD, "Diag") <- FALSE
    attr(PD, "Upper") <- FALSE
    class(PD) <- "dist"
    return(PD)
}

