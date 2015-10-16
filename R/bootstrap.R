bootstrap.pml = function (x, bs = 100, trees = TRUE, multicore=FALSE, mc.cores = NULL, ...) 
{
#    multicore=FALSE,
#    multicore <- mc.cores > 1L
    
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
        #  !require(parallel) ||      
#        if (.Platform$GUI!="X11") {
#            warning("package 'parallel' not found or GUI is used, 
#                    bootstrapping is performed in serial")
#        } 
#        else {       
            res <- mclapply(BS, pmlPar, x, trees = trees, ..., mc.cores = mc.cores)
            eval.success <- TRUE
#        } 
    }
    if (!eval.success) res <- lapply(BS, pmlPar, x, trees = trees, ...)
    if (trees) {
        class(res) = "multiPhylo"
        res = .compressTipLabel(res) # save memory
    }
    res
}


bootstrap.phyDat <- function (x, FUN, bs = 100, multicore=FALSE, mc.cores = NULL, ...) 
{
    if(multicore && is.null(mc.cores)){
        mc.cores <- detectCores()
    }
    weight = attr(x, "weight")
    v = rep(1:length(weight), weight)
    BS = vector("list", bs)
    for(i in 1:bs)BS[[i]]=tabulate(sample(v, replace=TRUE),length(weight)) 
    fitPar <- function(weights, data, ...){     
        ind <- which(weights > 0)
        data <- getRows(data, ind)
        attr(data, "weight") <- weights[ind]
        FUN(data,...)        
    }
    if(multicore) res <- mclapply(BS, fitPar, x, ..., mc.cores = mc.cores) 
    else res <- lapply(BS, fitPar, x, ...)
    if(class(res[[1]]) == "phylo"){
        class(res) <- "multiPhylo"   
        res = .compressTipLabel(res) # save memory
    }
    res 
}


matchEdges = function(tree1, tree2){
    bp1 = bip(tree1)
    bp2 = bip(tree2)
    l = length(tree1$tip)
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
    tree$tip.label <- tree$tip.label[ind]
    ind2 <- match(1:length(ind), tree$edge[, 2])
    tree$edge[ind2, 2] <- order(ind)
    tree
}


plotBS <- function (tree, BStrees, type = "unrooted", bs.col = "black", 
                    bs.adj = NULL, p=80, ...) 
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
    BStrees <- .uncompressTipLabel(BStrees)
    if(any(unlist(lapply(BStrees, is.rooted)))){
        BStrees <- lapply(BStrees, unroot)   
    }
    x = prop.clades(tree, BStrees)
    x = round((x/length(BStrees)) * 100)
    tree$node.label = x
    label = c(rep(0, length(tree$tip)), x)
    ind <- get("last_plot.phylo", envir = .PlotPhyloEnv)$edge[, 
                                                              2]
    if (type == "phylogram" | type == "cladogram") {
        root = getRoot(tree)
        label = c(rep(0, length(tree$tip)), x)
        label[root] = 0
        ind2 = matchEdges(tree2, tree)
        label = label[ind2]
        ind = which(label > p)
        #        browser()
        if (is.null(bs.adj)) 
            bs.adj = c(1, 1)
        if(length(ind)>0)nodelabels(text = label[ind], node = ind, frame = "none", 
                                    col = bs.col, adj = bs.adj, ...)
    }
    else {
        if (is.null(bs.adj)) 
            bs.adj = c(0.5, 0.5)
        ind2 = which(label[ind]>p)
        if(length(ind2>0))edgelabels(label[ind][ind2],ind2, frame = "none", col = bs.col, 
                                     adj = bs.adj, ...)
    }
    invisible(tree)
}


plotBS.Old <- function (tree, BStrees, type = "unrooted", bs.col = "black", 
                        bs.adj = NULL, p=80, ...) 
{
    # prop.clades raus??
    prop.clades <- function(phy, ..., part = NULL, rooted = FALSE) {
        if (is.null(part)) {
            obj <- list(...)
            if (length(obj) == 1 && class(obj[[1]]) != "phylo") 
                obj <- unlist(obj, recursive = FALSE)
            if (!identical(phy$tip, obj[[1]]$tip)) 
                obj[[1]] = checkLabels(obj[[1]], phy$tip)
            part <- prop.part(obj, check.labels = TRUE)
        }
        bp <- prop.part(phy)
        if (!rooted) {
            bp <- postprocess.prop.part(bp)
            part <- postprocess.prop.part(part)
        }
        n <- numeric(phy$Nnode)
        for (i in seq_along(bp)) {
            for (j in seq_along(part)) {
                if (identical(bp[[i]], part[[j]])) {
                    n[i] <- attr(part, "number")[j]
                    done <- TRUE
                    break
                }
            }
        }
        n
    }
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                              "unrooted", "radial"))
    if (type == "phylogram" | type == "cladogram") {
        if (!is.rooted(tree) & !is.null(tree$edge.length)) 
            tree2 = midpoint(tree)
        else tree2=tree
        plot(tree2, type = type, ...)
    }
    else plot(tree, type = type, ...)
    BStrees <- .uncompressTipLabel(BStrees)
    x = prop.clades(tree, BStrees)
    x = round((x/length(BStrees)) * 100)
    tree$node.label = x
    label = c(rep(0, length(tree$tip)), x)
    ind <- get("last_plot.phylo", envir = .PlotPhyloEnv)$edge[, 
                                                              2]
    if (type == "phylogram" | type == "cladogram") {
        root = getRoot(tree)
        label = c(rep(0, length(tree$tip)), x)
        label[root] = 0
        ind2 = matchEdges(tree2, tree)
        label = label[ind2]
        ind = which(label > p)
        #        browser()
        if (is.null(bs.adj)) 
            bs.adj = c(1, 1)
        if(length(ind)>0)nodelabels(text = label[ind], node = ind, frame = "none", 
                                    col = bs.col, adj = bs.adj, ...)
    }
    else {
        if (is.null(bs.adj)) 
            bs.adj = c(0.5, 0.5)
        ind2 = which(label[ind]>p)
        if(length(ind2>0))edgelabels(label[ind][ind2],ind2, frame = "none", col = bs.col, 
                                     adj = bs.adj, ...)
    }
    invisible(tree)
}
