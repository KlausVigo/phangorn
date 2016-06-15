bootstrap.pml = function (x, bs = 100, trees = TRUE, multicore=FALSE, mc.cores = NULL, ...) 
{
    if(multicore && is.null(mc.cores)){
        mc.cores <- detectCores()
    }
    
    data = x$data
    weight = attr(data, "weight")
    v = rep(1:length(weight), weight)
    BS = vector("list", bs)
    for (i in 1:bs) BS[[i]] = tabulate(sample(v, replace = TRUE), 
                                       length(weight))
    pmlPar <- function(weights, fit, trees = TRUE, ...) {
        data = fit$data
        ind <- which(weights > 0)
        data <- getRows(data, ind)
        attr(data, "weight") <- weights[ind]
        fit = update(fit, data = data)
        fit = optim.pml(fit, ...)
        if (trees) {
            tree = fit$tree
            return(tree)
        }
        attr(fit, "data") = NULL
        fit
    }
    eval.success <- FALSE
    if (!eval.success & multicore) {
        res <- mclapply(BS, pmlPar, x, trees = trees, ..., mc.cores = mc.cores)
        eval.success <- TRUE
    }
    if (!eval.success) res <- lapply(BS, pmlPar, x, trees = trees, ...)
    if (trees) {
        class(res) = "multiPhylo"
        res = .compressTipLabel(res) # save memory
    }
    res
}


bootstrap.phyDat <- function (x, FUN, bs = 100, multicore=FALSE, mc.cores = NULL, jumble=TRUE, ...) 
{
    if(multicore && is.null(mc.cores)){
        mc.cores <- detectCores()
    }
    weight = attr(x, "weight")
    v = rep(1:length(weight), weight)
    BS = vector("list", bs)
    for(i in 1:bs)BS[[i]]=tabulate(sample(v, replace=TRUE),length(weight)) 
    if(jumble){
        J = vector("list", bs)
        l = length(x)
        for(i in 1:bs) J[[i]] = list(BS[[i]], sample(l))
    }
    fitPar <- function(weights, data, ...){     
        ind <- which(weights > 0)
        data <- getRows(data, ind)
        attr(data, "weight") <- weights[ind]
        FUN(data,...)        
    }
    fitParJumble <- function(J, data, ...){     
        ind <- which(J[[1]] > 0)
        data <- getRows(data, ind)
        attr(data, "weight") <- J[[1]][ind]
        data <- subset(data, J[[2]])
        FUN(data,...)        
    }
    if(multicore){
        if(jumble) res <- mclapply(J, fitPar, x, ..., mc.cores = mc.cores) 
        else res <- mclapply(BS, fitPar, x, ..., mc.cores = mc.cores) 
    }        
    else{
        if(jumble) res <- lapply(J, fitParJumble, x, ...) 
        else res <- lapply(BS, fitPar, x, ...)
    } 
    if(class(res[[1]]) == "phylo"){
        class(res) <- "multiPhylo"   
        res = .compressTipLabel(res) # save memory
    }
    res 
}


matchEdges = function(tree1, tree2){
    bp1 = bip(tree1)
    bp2 = bip(tree2)
    l = length(tree1$tip.label)
    fn = function(x, y){
        if(x[1]==1)return(x)
        else return(y[-x])
    } 
    bp1[] = lapply(bp1, fn, 1:l)
    bp2[] = lapply(bp2, fn, 1:l)
    match(bp1, bp2)
}


checkLabels <- function(tree, tip){
    ind <- match(tip, tree$tip.label)
    if (any(is.na(ind)) | length(tree$tip.label) != length(tip))
        stop("tree has different labels")
    tree$tip.label <- tree$tip.label[ind]
    ind2 <- match(1:length(ind), tree$edge[, 2])
    tree$edge[ind2, 2] <- order(ind)
    tree
}


plotBS <- function (tree, BStrees, type = "unrooted", bs.col = "black", 
                    bs.adj = NULL, p=50, frame="none",...) 
{
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                              "unrooted", "radial"))
    if (type == "phylogram" | type == "cladogram") {
        if (!is.rooted(tree) & !is.null(tree$edge.length)) 
            tree2 = midpoint(tree)
        else tree2=tree
        plot(tree2, type = type, ...)
    }
    else plot(tree, type = type, ...)
    
    if(hasArg(BStrees)){
        BStrees <- .uncompressTipLabel(BStrees)
        if(any(unlist(lapply(BStrees, is.rooted)))){
            BStrees <- lapply(BStrees, unroot)   
        }
        x = prop.clades(tree, BStrees)
        x = round((x/length(BStrees)) * 100)
        tree$node.label = x
    }
    else{
        if(is.null(tree$node.label))stop("You need to supply BStrees or tree needs 
        needs BS-values as node.label")
        x <- tree$node.label
    }
    
    label = c(rep(0, length(tree$tip.label)), x)
    ind <- get("last_plot.phylo", envir = .PlotPhyloEnv)$edge[, 
                                                              2]
    if (type == "phylogram" | type == "cladogram") {
        root = getRoot(tree)
        label = c(rep(0, length(tree$tip.label)), x)
        label[root] = 0
        ind2 = matchEdges(tree2, tree)
        label = label[ind2]
        ind = which(label > p)
        #        browser()
        if (is.null(bs.adj)) 
            bs.adj = c(1, 1)
        if(length(ind)>0)nodelabels(text = label[ind], node = ind, frame = frame, 
                                    col = bs.col, adj = bs.adj, ...)
    }
    else {
        if (is.null(bs.adj)) 
            bs.adj = c(0.5, 0.5)
        ind2 = which(label[ind]>p)
        if(length(ind2>0))edgelabels(label[ind][ind2],ind2, frame = frame, col = bs.col, 
                                     adj = bs.adj, ...)
    }
    invisible(tree)
}


maxCladeCred <- function(x, tree=TRUE, part=NULL, rooted=TRUE){
    if(inherits(x, "phylo")) x <- c(x)
    if(!rooted){
        x <- lapply(x, unroot) 
        class(x) <- "multiPhylo"
        x <- .compressTipLabel(x)
    }    
    if(is.null(part))pp <- prop.part(x)
    else pp <- part
    pplabel <- attr(pp, "labels")
    if(!rooted)pp <- oneWise(pp)
    x <- .uncompressTipLabel(x)
    class(x) <- NULL
    m <- max(attr(pp, "number"))
    nb <- log( attr(pp, "number") / m )
    l <-  length(x)
    res <- numeric(l)
    for(i in 1:l){
        tmp <- checkLabels(x[[i]], pplabel)
        ppi <- prop.part(tmp)  # trees[[i]]
        if(!rooted)ppi <- oneWise(ppi)
        indi <- fmatch(ppi, pp)
        if(any(is.na(indi))) res[i] <- -Inf
        else res[i] <- sum(nb[indi])
    }
    if(tree) {
        k <- which.max(res)
        tr <- x[[k]]
        attr(tr, "clade.credibility") <- res[k]
        return(tr)
    }    
    res
}


mcc <- maxCladeCred


cladeMatrix <- function(x, rooted=FALSE){
    if(!rooted){
        x <- .uncompressTipLabel(x)
        x <- lapply(x, unroot) 
        class(x) <- "multiPhylo"
        x <- .compressTipLabel(x)
    }    
    pp <- prop.part(x)
    pplabel <- attr(pp, "labels")
    if(!rooted)pp <- oneWise(pp)
    x <- .uncompressTipLabel(x)
    class(x) <- NULL
    nnodes <- sapply(x, Nnode)
    l <-  length(x)
    from <- cumsum(c(1, nnodes[-l]))
    to <- cumsum(nnodes)
    
    ivec <- integer(to[l])
    pvec <- c(0,to)
    
    res <- vector("list", l)
    k=1
    for(i in 1:l){
        ppi <- prop.part(x[[i]])  
        if(!rooted)ppi <- oneWise(ppi)
        indi <- sort(fmatch(ppi, pp))
        ivec[from[i]:to[i]] = indi
    }
    X <- sparseMatrix(i=ivec, p=pvec, dims=c(length(pp),l))
    list(X=X, prop.part=pp)
}


moving_average <- function(x, window=50){
     cx <- c(0, cumsum(x))
     (cx[(window+1):length(cx)] - cx[1:(length(cx)-window)])/(window)
}
