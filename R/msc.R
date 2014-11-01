mapping <- function(x, header=TRUE, ...){
    if(is.character(x)) x <- read.table(x, header=header, stringsAsFactors=FALSE, ...)
    lab = unique(x[,2])
    allNames = c(x[,1], lab)
    Y=matrix(c(x[,2], lab), ncol=1, dimnames=list(allNames, NULL))
    phyDat(Y, "USER", levels=lab)
}


coalBranch <- function(coaltimes, tau1, tau2, k, theta=1){
    if(k==1)return(0)  
    if(!is.null(coaltimes))tt = diff.default(c(tau1, sort(coaltimes), tau2))
    else return(0)
    res = 0 
    l = length(tt) - 1L  
    n=k
    if(l>0){
        for(i in 1:l){
            res = res + log(2) - log(theta) + (-n*(n-1)*tt[i] / theta)
            n = n - 1L
        }
    }
    if(n>1L) res = res + ((-n*(n-1) / theta) * tt[l+1])   
    res
}


MSC2 <- function(sTree, gTree, X, theta=100, edge=NULL, desc=NULL, sst=NULL, snh=NULL, gnh=NULL, gst=NULL){
    if(class(gTree)=="multiPhylo") gTree = unclass(gTree)
    if(class(gTree)=="phylo") gTree = list(gTree)
    
    nsp <- length(sTree$tip) 
    if(is.null(snh)) snh <- nodeHeight(sTree) 
    if(is.null(sst)) sst <- ancstat(sTree, X)
    
    edge <- sTree$edge
    if(is.null(desc))desc <- Descendants(sTree, 1:max(sTree$edge), type="all")
    root <- getRoot(sTree)
    
    # start a loop over the trees   
    # should be a list
    # assume gTree is a list of trees no multiphylo and not compressed
    
    if(is.null(gnh))gnh <- lapply(gTree, nodeHeight)  # scaling on these
    if(is.null(gst))gst <- lapply(gTree, ancstat, X)  # usually fixed    
    ntips = sapply(gTree, function(x)length(x$tip.label))
    res = numeric(length(gTree))
    for(i in 1:length(gTree)){
        res[i] = MSC.fit(sst, snh, gst[[i]], gnh[[i]], edge, ntips[i], desc, root, nsp, theta)       
    }
    res
}    



MSC3 <- function(sTree, gTree, X, theta=100, edge=NULL, desc=NULL, sst=NULL, snh=NULL, gnh=NULL, gst=NULL){
    if(class(gTree)=="multiPhylo") gTree = unclass(gTree)
    if(class(gTree)=="phylo") gTree = list(gTree)
    
    nsp <- length(sTree$tip) 
    if(is.null(snh)) snh <- nodeHeight(sTree) 
    if(is.null(sst)) sst <- ancstat(sTree, X)
    
    edge <- sTree$edge
    if(is.null(desc))desc <- Descendants(sTree, 1:max(sTree$edge), type="all")
    root <- getRoot(sTree)
    
    # start a loop over the trees   
    # should be a list
    # assume gTree is a list of trees no multiphylo and not compressed
    
    if(is.null(gnh))gnh <- lapply(gTree, nodeHeight)  # scaling on these
    if(is.null(gst))gst <- lapply(gTree, ancstat, X)  # usually fixed
    tmpfun = function(x){
        m <- max(x$edge)
        res <- logical(m)
        res[ (length(x$tip.label)+1L): m] = TRUE
        res
    }
    ind1 <- lapply(gTree, tmpfun)
    
    for(i in 1:length(gTree)){
        ind = order(gnh [[i]])
        gnh[[i]] = gnh[[i]][ind]
        gst[[i]] = gst[[i]][ind,] 
        ind1[[i]] = ind1[[i]][ind,]
    }
    
    ntips = sapply(gTree, function(x)length(x$tip.label))
    res = numeric(length(gTree))
    for(i in 1:length(gTree)){
        res[i] = MSC.fit(sst, snh, gst[[i]], gnh[[i]], edge, ntips[i], desc, root, nsp, theta)       
    }
    res
}    



MSC.fit <- function(sst, snh, gst, gnh, edge, ntips, desc, root, nsp, theta){
    br = comp(gst, sst)
    ind1 = logical(length(br))
    ind1[ (ntips+1): length(br) ] = TRUE
    # loop for trees, triplets (edge length optimisation) & quartets (NNI moves)  
    l = nrow(edge) 
    ct = vector("list", l+1L)  # ??? max(l)
    # Ziel ct muss sortiert sein
    for(i in 1:l){ # for(i in blub)
        ind2= br %in% desc[[edge[i,2]]] # vielleicht ausserhalb der Schleife
        ind3=(gnh>=snh[edge[i,2]] & gnh<snh[edge[i,1]])
        ind = which(ind1 & ind2 & ind3)
#        print(ind)
        ct[[edge[i,2]]] = gnh[ind]    
    }
    
    ct[[root]] = gnh[gnh>max(snh)]   
    nrc = sapply(ct, length)
    k = integer(l+1)
    for(i in 1:nsp)k[i] = sum(br[!ind1]==i)
    
    for(i in 1:l){
        ei = edge[i,1]
        k[ei] = k[ei]+k[edge[i,2]] - nrc[edge[i,2]]  
    }
    res = numeric(l+1L)
    nhm = rep(100, length(snh))
    for(i in 1:l) nhm[edge[i,2]] = snh[edge[i,1]]
    #   for(i in 1:(l+1L))res[i] = coalBranchNeu(ct[[i]], snh[i], nhm[i], k[i], theta)
    for(i in 1:(l+1L))res[i] = coalBranch(ct[[i]], snh[i], nhm[i], k[i], theta)
    sum(res)  
}



indEdges <- function(tree, x){
    h <- nodeHeight(tree)
    ind <- which( h[tree$edge[,1]]>x & h[tree$edge[,2]]<x )
    tree$edge[ind,2]
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
    res[1:nTips, ] = contrast[as.numeric(x)[match(phy$tip, names(x))],, drop=FALSE]
    for(i in 1:length(pa)){
        res[pa[i],] = res[pa[i],] | res[ch[i],]    
    }
    res
}


plotA = function (tree, data, i = 1, col = NULL, ...) 
{
    if(class(data)=="phyDAt"){
        y = subset(data, , i)
        nc = attr(data, "nc")
        y = matrix(unlist(y[]), ncol = nc, byrow = TRUE)
        l = dim(y)[1]
        dat = matrix(0, l, nc)
        for (i in 1:l) dat[i, ] = y[[i]]
        levels = attr(data, "levels")
    }
    else{
        dat = data
        y = data
        nc = ncol(y)
        levels=colnames(data) 
    } 
    
    args <- list(...)
    CEX <- if ("cex" %in% names(args)) 
        args$cex
    else par("cex")
    xrad <- CEX * diff(par("usr")[1:2])/50
    
    plot(tree, label.offset = 1.1 * xrad, plot = FALSE, ...)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    XX <- lastPP$xx
    YY <- lastPP$yy
    xrad <- CEX * diff(lastPP$x.lim * 1.1)/50
    par(new = TRUE)
    plot(tree, label.offset = 1.1 * xrad, plot = TRUE, ...)
    if (is.null(col)) 
        col = rainbow(nc)
    if (length(col) != nc) 
        warning("Length of color vector differs from number of levels!")
    BOTHlabels(pie = y, XX = XX, YY = YY, adj = c(0.5, 0.5), 
               frame = "rect", pch = NULL, sel = 1:length(XX), thermo = NULL, 
               piecol = col, col = "black", bg = "lightblue", horiz = FALSE, 
               width = NULL, height = NULL)
    legend("bottomright", levels, text.col = col)
}


msc = function(sTree, gTrees, X, optEdge=TRUE, rearrangement="NNI", optCoal="singleRate", optGeneRate=FALSE, theta=NULL){
    if(is.null(theta))theta = rep(1, length(sTree$Nnode))
#    gTrees = lapply(pm, function(x)x$tree)
#    treeLL = sapply(pm, logLik)
    coalLL = MSC2(sTree, gTrees, X, theta)
#    start = sum(coalLL + treeLL)

    start = sum(coalLL)
    logLik= start


# nur coalBranch von theta abhaengig von 5.4 to 2.5 ??
    optSingle = function(sc, st, gTrees, X){
#        theta = sc*theta  , theta
        sum(MSC2(st, gTrees, X, theta=sc))
    }
    
    # scale Edge 
    # scale Rates   
# optimise theta (global)
    if(optCoal=="singleRate"){
        res = optimize(optSingle, c(0.1, 100000), st=sTree, gTrees=gTrees, X=X, 
             maximum=TRUE)  #  theta=theta,
        if(res[[2]] > logLik){
            theta <- res[[1]] #theta
            logLik <- res[[2]]
        }
    }    

    if(optGeneRate){
        res = numeric(length(gTrees))
        tmpfun = function(sc, st, gTrees, X, theta){
            gTrees$edge.length = gTrees$edge.length*sc
            MSC2(st, gTrees, X, theta)
        } 
        for(i in 2:length(gTrees))res[i]=optimize(tmpfun, c(0.01, 100), st=sTree, gTrees=gTrees[[i]], X=X, theta=theta, maximum=TRUE)[[1]]
        for(i in 1:length(gTrees))gTrees[[i]]$edge.length = gTrees[[i]]$edge.length * res[i]
    }
    # optimize population size
    attr(sTree, "logLik") = logLik
    attr(sTree, "theta") = theta
    sTree
}


optimCoalEdge <- function(sTree, gTrees, Y, theta, trace=0){
    
    nTips = as.integer(length(sTree$tip.label))   
    
    if (is.null(attr(stree, "order")) || attr(sTree, "order") == "cladewise") 
        sTree <- reorder.phylo(sTree, "postorder") 
    if(!is.rooted(sTree))stop("species tree must be rooted!")
    
    getEL = function(t, nh){
        el = numeric(3)
        tnh = max(nh[1:2])
        l = nh[3] - tnh
        el[3] = (1-t) * l
        el[1:2] = t*l
        if(nh[1] > nh[2]) el[2] = el[2] + nh[1] - nh[2] 
        else el[1] = el[1] + nh[2] - nh[1] 
        el
    }    
    
    
    optEdge = function(t, sTree, gTrees, Y, theta, nh, kids){
        el = getEL(t, nh)  
        sTree <- changeEdgeLength(sTree, kids, el) 
        ll <- sum(MSC2(sTree, gTrees, Y, theta))
        ll
    }
    
    optEdge2 = function(t, sTree, gTrees, Y, theta, kids, indk, ch, i){
        el = sTree$edge.length  
        sTree <- changeEdgeLength(sTree, kids, el[indk]+t) 
        sTree <- changeEdgeLength(sTree, ch, el[i]-t)        
        ll <- sum(MSC2(sTree, gTrees, Y, theta))
        ll
    }
  
    
#    sTree2 <- changeEdgeLength(sTree, kid, el[indk]+upper) 
#    sTree2 <- changeEdgeLength(sTree2, ch, el[i]-upper)        
#    sum(MSC2(sTree2, gTrees, Y, theta))
    
    # mit reference or 
#    trace=0,
    scaleEdges = function(t=1, sTree, gTrees, Y, theta,...){
        fn = function(t, sTree, gTrees, Y, theta,...){
            sTree$edge.length = sTree$edge.length*t
            sum(MSC2(sTree, gTrees, Y, theta))
        }
        optimize(f=fn, interval=c(0.25,4), sTree=sTree, gTrees=gTrees, Y=Y, 
                 theta=theta, maximum = TRUE, tol = .00001)
    }
    
    
    child = sTree$edge[, 2]
    parent = sTree$edge[, 1]
    ll <- sum(MSC2(sTree, gTrees, Y, theta))
    llstart <- ll
    loglik <- llstart 
    eps=.00001
    iter = 1  
    
    EL = numeric(max(sTree$edge)) 
    EL[child] = sTree$edge.length  #child2
    
    change = numeric(length(parent)) + 1
    
    rootNode = getRoot(sTree)    
    anc = Ancestors(sTree, 1:max(sTree$edge), "parent")  
    cvector = allChildren(sTree)
    
    
    loli <- rootNode                
    pa <-rootNode
    nchanges = 0
    
    CH = rev(child[child>nTips])
    eps=1

    while(eps>.001){        
        t = scaleEdges(t=1, sTree=sTree, gTrees=gTrees, Y=Y, theta=theta) 
        if(t[[2]] > ll){
            ll <- t[[2]]
            sTree$edge.length = sTree$edge.length * t[[1]]
# checking            
            sum(MSC2(sTree, gTrees, Y, theta))
        }
        
        for(i in 1:length(child)){
            ch <-  child[i] #child2[i]
            if(ch>nTips){
#                browser()
                dad <- anc[ch]#parent[i] #parent2[i]
                #               cat(ch, dad, "\n")
                nh=nodeHeight(sTree)
                el <- sTree$edge.length
                kid = cvector[[ch]]
                indk = match(kid, child)
 #     i
                lower = min(el[indk])    
                upper = max(el[i])
            
                nhi = nh[c(cvector[[ch]], dad)]
                
                kids = c(cvector[[ch]], ch)        
 
                res= optimize(f=optEdge, interval=c(0,1), sTree, gTrees, Y, theta, nhi, kids, maximum=TRUE)
                res2 <- optimize(f=optEdge2, interval=c(-lower, upper), sTree, gTrees, Y, theta, kid, indk, ch, i, maximum=TRUE)
 
                if(res[[2]] > (ll + 1e-6)){
                    browser()
                    t = res[[1]]
                #            print(t)
                    el = getEL(t, nhi)  
                    sTree <- changeEdgeLength(sTree, kids, el)
                    sum(MSC2(sTree, gTrees, Y, theta))
                    ll <- res[[2]]
                
                } 
            }
            
        }
        test= sum(MSC2(sTree, gTrees, Y, theta))
        eps=abs(test-loglik)
        if(trace>0) cat("eps:", eps, "\n")
        if(trace>0) cat("logLik:", loglik, "==>", test, "\n")
        loglik <- test
        ll=test
    }
    list(sTree=sTree, end=test, start=llstart)
}




optimCoalNni <- function(tree, gTree, Y, theta){
    
    nTips = as.integer(length(tree$tip.label))   
    
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree <- reorder.phylo(tree, "postorder") 
    if(!is.rooted(tree))stop("tree must be rooted")
    
    getEL = function(t, nh){
        el = numeric(3)
        tnh = max(nh[1:2])
        l = nh[3] - tnh
        el[3] = (1-t) * l
        el[1:2] = t*l
        if(nh[1] > nh[2]) el[2] = el[2] + nh[1] - nh[2] 
        else el[1] = el[1] + nh[2] - nh[1] 
        el
    }    
    
    
    optEdge = function(t, tree, gTree, Y, theta, nh, kids){
        el = getEL(t, nh)  
        tree <- changeEdgeLength(tree, kids, el) 
        ll <- sum(MSC2(gTree, tree, Y, theta))
        ll
    }
    
    scaleEdges = function(t=1, trace=0, gTree, tree, Y, theta,...){
        fn = function(t, gTree, tree, Y, theta,...){
            tree$edge.length = tree$edge.length*t
            sum(MSC2(gTree, tree, Y, theta))
        }
        optimize(f=fn, interval=c(0.25,4), tree=tree, data=data, ..., maximum = TRUE,
                 tol = .00001)
    }
    
    
    child = tree$edge[, 2]
    parent = tree$edge[, 1]
    ll <- sum(MSC2(gTree, tree, Y, theta))
    llstart <- ll
    eps=.00001
    iter = 1  
    
    EL = numeric(max(tree$edge)) 
    EL[child] = tree$edge.length  #child2
    
    change = numeric(length(parent)) + 1
    
    rootNode = getRoot(tree)    
    anc = Ancestors(tree, 1:max(tree$edge), "parent")  
    cvector = allChildren(tree)
    
    
    loli <- rootNode                
    pa <-rootNode
    nchanges = 0
    
    CH = rev(child[child>nTips])
    eps=1
    while(eps>.001){
        t= scaleEdges(t=1, trace=0, gTree, tree, Y, theta)    
        for(i in 1:length(child)){
            ch <-  child[i] #child2[i]
            if(ch>nTips){
                
                dad <- anc[ch]#parent[i] #parent2[i]
                #               cat(ch, dad, "\n")
                nh=nodeHeight(tree)
                nhi = nh[c(cvector[[ch]], dad)]
                
                kids = c(cvector[[ch]], ch)
                res= optimize(f=optEdge, interval=c(0,1), tree, gTree, Y, theta, nhi, kids, maximum=TRUE)
                t = res[[1]]
                #            print(t)
                el = getEL(t, nhi)  
                tree <- changeEdgeLength(tree, kids, el)
                
            }
            
        }
        test= sum(MSC2(gTree, tree, Y, theta))
        eps=abs(test-ll)
        ll=test
    }
    list(tree=tree, end=test, start=llstart)
}

