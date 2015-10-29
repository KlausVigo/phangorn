# function to  extract quartets
# order cladewise
indexNNI3 <- function(tree){
    parent = tree$edge[, 1]
    child = tree$edge[, 2]
    ind = reorder(tree)$edge[,2]
    nTips <- length(tree$tip.label)
    ind = ind[ind>nTips]
    #    ind = which(child %in% parent)
    Nnode = tree$Nnode
    #     a         d
    #      \       /
    #       e-----f       c is closest to root, f is root from subtree
    #      /       \  
    #     b         c     c(a,b,c,d,e,f)     
    edgeMatrix = matrix(0,(Nnode-1), 6)
    
    pvector <- numeric(max(parent))
    pvector[child] <- parent
    tips  <- !logical(max(parent))
    tips[parent] <-  FALSE
    #    cvector <- allCildren(tree)  
    cvector <- vector("list",max(parent))   
    for(i in 1:length(parent))  cvector[[parent[i]]] <- c(cvector[[parent[i]]], child[i]) 
    k=1L
    for(i in ind){   
        f = pvector[i] #f
        ab = cvector[[i]] #a,b
        ind1 = cvector[[f]] #c,d
        cd = ind1[ind1 != i]
        if(pvector[f])cd=c(pvector[f], cd) # cd  
        edgeMatrix[k, 1:6] = c(ab,cd,i,f)
        k=k+1L
    } 
    edgeMatrix
}

# EL ausserhalb
index2tree <- function(x, tree, root=length(tree$tip.label)+1L){
    EL = numeric(max(tree$edge))
    EL[tree$edge[,2]] = tree$edge.length
    pa = c(5L,5L,6L,6L,6L)
    ch = c(1L,2L,5L,4L,3L)
 elind = c(1L,2L,5L,4L,6L)
    if(x[6L]==root) el = EL[x[ch]]
    else   el = EL[x[elind]]  
    structure(list(edge = structure(c(x[pa], x[ch]), .Dim = c(5L, 2L)), 
                   tip.label = tree$tip.label, edge.length = el, Nnode = 2L), 
              .Names = c("edge", "tip.label", "edge.length", "Nnode"), class = "phylo", order = "postorder")
}


index2tree2 <- function(x, tree, root=length(tree$tip.label)+1L){
    EL = numeric(max(tree$edge))
    EL[tree$edge[,2]] = tree$edge.length
    pa = c(6L,6L,5L,5L,5L)
    ch = c(3L,4L,6L,1L,2L)
    elr = c(3L,4L,5L,1L,2L)
    eln = c(6L,4L,5L,1L,2L)
    if(x[6L]==root) el = EL[x[elr]]
    else   el = EL[x[eln]]  
    structure(list(edge = structure(c(x[pa], x[ch]), .Dim = c(5L, 2L)), 
                   tip.label = tree$tip.label, edge.length = el, Nnode = 2L), 
              .Names = c("edge", "tip.label", "edge.length", "Nnode"), class = "phylo", order = "postorder")
}



reorderQuartet <- function(tree){
    ind <- c(4,5,3,1,2)
    tree$edge <-  tree$edge[ind,]    
    tree$edge.length <-  tree$edge.length[ind]
    attr(tree, "order") <- "cladewise"
    tree
}


optimQuartet <- function (tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, anc, 
          control = pml.control(epsilon = 1e-08, maxit = 10, trace=0), ...) 
{
    nTips <- length(tree$tip)
    el <- tree$edge.length
    tree$edge.length[el < 1e-08] <- 1e-08
    oldtree = tree
    k = length(w)    

    loglik = pml.fit2(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k)
    print(paste("first loglik", loglik) )    
    
    start.ll <- old.ll <- loglik 
    contrast <- attr(data, "contrast")
    contrast2 <- contrast %*% eig[[2]] 
    evi = (t(eig[[3]]) * bf)
    weight <- attr(data, "weight")
    eps = 1
    iter = 0
    
    treeP = tree
    tree = reorderQuartet(tree)

    child = tree$edge[, 2]
    parent = tree$edge[, 1]
    m <- max(tree$edge)
    pvec <- integer(m)
    pvec[child] <- parent
    
    EL = numeric(m)
    EL[child] = tree$edge.length
    
#    n = length(tree$edge.length)  
    
    nr = as.integer(length(weight))
    nc = as.integer(length(bf))
    nco = as.integer(nrow(contrast))
    eve = eig[[2]]
    lg = k
    rootNode = getRoot(tree)         
    ScaleEPS = 1.0/4294967296.0
#    anc = Ancestors(tree, 1:m, "parent")  
    anc0 = as.integer(c(0L, anc))
    ancQ = anc
    ancQ[child] <- parent
    ancQ <- as.integer(c(0L, ancQ))
    
    blub3 <- .Call("extractScale", as.integer(rootNode), w, g, as.integer(nr), as.integer(nc), as.integer(nTips))
    rowM = apply(blub3, 1, min)       
    blub3 = (blub3-rowM) 
    blub3 = ScaleEPS ^ (blub3) 
    #        blub3[]=1.0   
    
    while (eps > control$eps && iter < control$maxit) {
        
        print(paste("root", rootNode))
        
        
        EL <- .Call("optE", as.integer(parent), as.integer(child), 
                    as.integer(ancQ), eig, evi, EL, w, g, as.integer(nr), as.integer(nc), 
                    as.integer(nTips), as.double(contrast), 
                    as.double(contrast2), nco, blub3, data, as.double(weight), as.double(ll.0))       
        iter = iter + 1
        
#        print("edge length new")
#        print(EL[treeP$edge[,2]])
#        print("edge length old")
#        print(treeP$edge.length)
        
#browser()

        treeP$edge.length = EL[treeP$edge[,2]]
        newll <- pml.fit2(treeP, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k)
        print(paste("nach iteration:",iter,"ll in optimEdge", newll))
        eps = ( old.ll - newll ) / newll
        if( eps <0 ) return(list(oldtree, old.ll))
#        print(eps)
        oldtree = treeP
        if(control$trace>1) cat(old.ll, " -> ", newll, "\n") 
        old.ll = newll
    }
    if(control$trace>0) cat(start.ll, " -> ", newll, "\n")
    list(tree=treeP, logLik=newll, c(eps, iter))
}



pml.nni2 <- function (tree, data, w, g, eig, bf, ll.0, ll, ...) 
{        
    k = length(w)
    INDEX <-  indexNNI3(tree)
    
    
    tmpl <- pml.fit4(tree, data, bf=bf, g=g, w=w, eig=eig, INV=INV, ll.0=ll.0, k=k)
    
    nr <- attr(data, "nr")
    nc <- attr(data, "nc")
    parent = tree$edge[,1]
    child = tree$edge[,2]
    weight = attr(data, "weight")
    contrast <- attr(data, "contrast")
    contrast2 <- contrast %*% eig[[2]] 
    evi = (t(eig[[3]]) * bf)
    
    nTips = length(tree$tip.label)
    EL <- numeric(max(parent)) 
    EL[child] <- tree$edge.length
    m <- dim(INDEX)[1]
    loglik = numeric(2*m)
    edgeMatrix <- matrix(0, 2*m, 5)
    
    anc <- Ancestors(tree, 1:max(tree$edge), "parent")  
    loli <- getRoot(tree)
    
#    l = length(datp[, 1])
    for(i in 1:m){
        ei = INDEX[i,]
#        el0 = evector[INDEX[i,]]
        tree0 = index2tree(INDEX[i, ], tree, nTips+1L)
        
        ch = ei[5]
        pa = ei[6]

        while(pa != loli){
            tmpr = match(loli, INDEX[,5])
            treetmp = index2tree(INDEX[tmpr, ], tree, nTips+1L)
            tmpl <- pml.fit4(treetmp, data, bf=bf, g=g, w=w, eig=eig, INV=INV, ll.0=ll.0, k=k)
            print(paste("tmp:", tmpl))
            loli = anc[loli]
        }
        
        tree0$edge.length = tree0$edge.length+.01
        tmp <- pml.fit4(tree0, data, bf=bf, g=g, w=w, eig=eig, INV=INV, ll.0=ll.0, k=k)
           print(paste("vor optimQuartet", tmp))
        browser()
           fit1 <- optimQuartet(tree0, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, anc=anc)

           ll1 <- fit1$logLik  
           print(paste("nach optimierung", fit1$logLik))
           tmp <- pml.fit4(fit1$tree, data, bf=bf, g=g, w=w, eig=eig, INV=INV, ll.0=ll.0, k=k)
#       }
        
#        if(ch>nTips){
            tree2 <- index2tree2(INDEX[i, ], tree, nTips+1L)     
            tmp2 <- pml.fit4(tree2, data, bf=bf, g=g, w=w, eig=eig, INV=INV, ll.0=ll.0, k=k)
            loli <- getRoot(tree2)
#        }    
#        else{}
        
        print(paste("tmp", tmp))
        print(paste("tmp2", tmp2))
        print(paste("loli", loli))
        
        
        
#        ext = ei[1:4] < nTips+1L
#        if (!(ei[5] %in% rootEdges)) dat1 = datp[, ei[1], drop = FALSE]
#        else{ if(ext[1]) dat1 = data[[ ei[1] ]]
#        else dat1 = .dat[, ei[1], drop=FALSE]
#        } 
#        if(ext[2]) dat2 = data[[ ei[2] ]]
#        else dat2 = .dat[, ei[2], drop=FALSE] 
#        if(ext[3]) dat3 = data[[ ei[3] ]]
#        else dat3 = .dat[, ei[3], drop=FALSE]
#        if(ext[4]) dat4 = data[[ ei[4] ]]
#        else dat4 = .dat[, ei[4], drop=FALSE]
        
#        new1 <- optim.quartet2(el0[c(1, 3, 2, 4, 5)], eig, bf, 
#                               dat1, dat3, dat2, dat4, g, w, weight, ll.0, llcomp=ll, evi=evi, contrast=contrast, contrast2=contrast2, ext=ext[c(1, 3, 2, 4)])
#        new2 <- optim.quartet2(el0[c(1, 4, 3, 2, 5)], eig, bf,  
#                               dat1, dat4, dat3, dat2, g, w, weight, ll.0, llcomp=ll, evi=evi, contrast=contrast, contrast2=contrast2, ext=ext[c(1, 4, 3, 2)])
        
        
#        loglik[(2*i)-1]=new1[[2]]
#        loglik[(2*i)]=new2[[2]] 
#        edgeMatrix[(2*i)-1,]=new1[[1]]
#        edgeMatrix[(2*i),]=new2[[1]]           
    }
#    swap <- 0
#    eps0 <- 1e-6
#    candidates <- loglik > ll + eps0
    
#    nr <- as.integer(attr(data, "nr")) 
#    nc <- as.integer(attr(data, "nc"))
#    nTips <- as.integer(length(tree$tip.label))
    
#    while(any(candidates)){     
#        ind = which.max(loglik)
#        loglik[ind]=-Inf
#        if( ind %% 2 ) swap.edge = c(2,3)
#        else swap.edge = c(2,4)
#        tree2 <- changeEdge(tree, INDEX[(ind+1)%/%2,swap.edge], INDEX[(ind+1)%/%2,], edgeMatrix[ind,])
        
#        test <- pml.fit(tree2, data, bf = bf, k=k, g=g, w=w, eig=eig, ll.0=ll.0, ...) 
#        if(test <= ll + eps0) candidates[ind] = FALSE
#        if(test > ll + eps0) {
#            ll = test 
#            swap=swap+1
#            tree <- tree2
#            indi <- which(rep(colSums(apply(INDEX,1,match,INDEX[(ind+1)%/%2,],nomatch=0))>0,each=2))
#            candidates[indi] <- FALSE
#            loglik[indi] <- -Inf
#        }
#    } 
    list(tree=tree, ll=ll, swap=0)     
}




optim.pml2 <- function (object, optNni = FALSE, optBf = FALSE, optQ = FALSE, 
                       optInv = FALSE, optGamma = FALSE, optEdge = TRUE, optRate = FALSE, optRooted=FALSE, 
                       optRatchet = FALSE, 
                       control = pml.control(epsilon = 1e-8, maxit = 10, trace = 1L), 
                       model = NULL, subs = NULL, ratchet.par = list(prop = 1/3, iter = 10L, maxit = 100L), ...) 
{
    extras <- match.call(expand.dots = FALSE)$...
    pmla <- c("wMix", "llMix")
    wMix <- object$wMix
    llMix <- object$llMix
    if(is.null(llMix)) llMix=0
    if (!is.null(extras)) {
        names(extras) <- pmla[pmatch(names(extras), pmla)]
        existing <- match(pmla, names(extras))
        if (!is.na(existing[1])) 
            wMix <- eval(extras[[existing[1]]], parent.frame())
        if (!is.na(existing[2])) 
            llMix <- eval(extras[[existing[2]]], parent.frame())
    }
    tree = object$tree
    call = object$call
    ratchet=FALSE
    if(optRatchet == TRUE){
        if(optRooted==FALSE){
            optNni=TRUE
            optEdge=TRUE
            ratchet=TRUE
        }
    }
    
    data <- object$data
    addTaxa <- FALSE
    
    if(optNni) {
        dup <-  duplicated(data)
        if(any(dup)){ # && optNNI
            orig.data <- data
            addTaxa <- TRUE
            labels <- names(data)
            ldup <- labels[dup]
            ind2 <- match(subset(data, dup), data)
            tree2 <- drop.tip(tree, ldup)
            tree <- reorder(tree2, "postorder")
            mapping <- cbind(labels[dup], labels[ind2])
        }
        if(!is.binary.tree(tree)) 
            tree = multi2di(tree)
#        optEdge = TRUE     
    }
    if(is.rooted(tree)) {
        if(optRooted==FALSE && optEdge==TRUE){
            tree = unroot(tree)
            attr(tree, "order") <- NULL
            tree = reorder(tree, "postorder")
            warning("I unrooted the tree", call. = FALSE)
        }    
    }
    if(is.null(attr(tree, "order")) || attr(tree, "order") == 
       "cladewise") 
        tree <- reorder(tree, "postorder")
    if(any(tree$edge.length < 1e-08)) {
        tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
        # save to change to new update.pml       
        object <- update.pml(object, tree = tree)
    }
    if(optEdge & optRate) {
        warning("You can't optimise edges and rates at the same time, only edges are optimised!", call. = FALSE)
        optRate = FALSE
    }

    trace <- control$trace
    
    data = subset(data, tree$tip.label) 
    
    type <- attr(data, "type")
    if (type == "AA" & !is.null(model)){
        object = update(object, model=model)  
    }     
    if (type == "CODON") {
        dnds <- object$dnds 
        tstv <- object$tstv
        if(!is.null(model)){
            if(model == "codon0") optQ = FALSE
            else  optQ = TRUE
        }
    }       
    Q = object$Q
    if(is.null(subs)) subs = c(1:(length(Q) - 1), 0)
    bf = object$bf
    eig = object$eig
    inv = object$inv
    k = object$k
    if(k==1 & optGamma){
        optGamma = FALSE
        message('only one rate class, ignored optGamma')
    }
    shape = object$shape
    w = object$w
    g = object$g
    if (type == "DNA" & !is.null(model)) {
        tmp = subsChoice(model)
        optQ = tmp$optQ
        if (!optQ) 
            Q = rep(1, 6)
        optBf = tmp$optBf
        if (!optBf) 
            bf = c(0.25, 0.25, 0.25, 0.25)
        subs = tmp$subs
    }   
    ll0 <- object$logLik
    INV <- object$INV
    ll.0 <- object$ll.0
    rate <- object$rate
    ll = ll0
    ll1 = ll0
    opti = TRUE
    
    nr <- as.integer(attr(data, "nr")) 
    nc <- as.integer(attr(data, "nc"))
    nTips <- as.integer(length(tree$tip.label))
    
    #    on.exit(.C("ll_free"))
    #    .C("ll_init", nr, nTips, nc, as.integer(k))
    .INV <- .iind <- NULL
    on.exit({
        if(type=="CODON"){
            object$dnds = dnds
            object$tstv = tstv
        }
        
        tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q, 
                       levels = attr(data, "levels"), inv = inv, rate = rate, 
                       g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix, 
                       wMix = wMix, site = TRUE)
        
        df <- ifelse(optRooted, tree$Nnode, length(tree$edge.length))
        # length(tree$edge.length)    
        if (type == "CODON") {
            df <- df + (k > 1) + (inv > 0) + 
                length(unique(bf)) - 1 + (dnds != 1) + (tstv != 1) 
        }
        else df = df + (k > 1) + (inv > 0) + 
            length(unique(bf)) - 1 + length(unique(Q)) - 1
        
        if(addTaxa){
            #            pml.free()
            tree <- addAllTips(tree, mapping)
            data <- orig.data
            #            pml.init(subset(data, tree$tip.label), k) 
        }
        
        object = list(logLik = tmp$loglik, inv = inv, k = k, shape = shape, 
                      Q = Q, bf = bf, rate = rate, siteLik = tmp$siteLik, weight = attr(data, "weight"), 
                      g = g, w = w, eig = eig, data = data, model = model, 
                      INV = INV, ll.0 = ll.0, tree = tree, lv = tmp$resll, 
                      call = call, df = df, wMix = wMix, llMix = llMix)
        if (type == "CODON") {
            object$dnds <- dnds
            object$tstv <- tstv
        }
        class(object) = "pml"
        
        extras = pairlist(bf = bf, Q = Q, inv = inv, shape = shape, rate = rate)[c(optBf, optQ, optInv, optGamma, optRate)]
        if (length(extras)) {
            existing <- !is.na(match(names(extras), names(call)))
            for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
            if (any(!existing)) {
                call <- c(as.list(call), extras[!existing])
                call <- as.call(call)
            }
        }
        object$call = call   
        
        pml.free()
        return(object)
        #        rm(.INV, .iind)
    })
    pml.init(data, k)    
    
    if (optEdge) {
        
        # check if non-negative least-squares is better for start of optimisation
        treetmp <- nnls.phylo(tree, dist.ml(data))
        treetmp$edge.length[treetmp$edge.length < 1e-8] <- 1e-8
        tmplogLik <- pml.fit(treetmp, data, bf, k = k, inv = inv, g = g, w = w, 
                             eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix, wMix = wMix)
        if(tmplogLik>ll) tree <- treetmp
        
        res <- optimEdge(tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, INV=INV,
                         control = pml.control(epsilon = 1e-07, maxit = 5, trace=trace - 1)) 
        if(trace > 0) 
            cat("optimize edge weights: ", ll, "-->", res[[2]], "\n")  
        if (res[[2]] > ll){  
            ll <- res[[2]]
            tree <- res[[1]]
        }
    }
    if(optRooted){
        res <- optimRooted(tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, INV=INV, control = pml.control(epsilon = 1e-07, maxit = 10, trace = trace-1))
        if(trace > 0) 
            cat("optimize edge weights: ", ll, "-->", res[[2]], "\n")
        if(res[[2]] > ll){  
            ll <- res[[2]]
            tree <- res[[1]]
        }     
    }
    rounds = 1
    while (opti) {


        if (optEdge) {  
            res <- optimEdge(tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0,
                             control = pml.control(epsilon = 1e-08, maxit = 5, trace=trace - 1)) 
            if (trace > 0) 
                cat("optimize edge weights: ", ll, "-->", res[[2]], "\n")
            if (res[[2]] > ll){  
                ll <- res[[2]]
                tree <- res[[1]]
            }
        }

        if(optNni) {
            swap = 0
            iter = 1
            while (iter < 4) {
                
                    tmp <- pml.nni2(tree, data, w, g, eig, bf, ll.0, ll, ...) 
                    swap = swap + tmp$swap
                    res <- optimEdge(tmp$tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, control = pml.control(epsilon = 1e-08, maxit = 3, trace=0)) 
                    ll2 = res[[2]] 
                    tree <- res[[1]]
                
                if (trace > 0) 
                    cat("optimize topology: ", ll, "-->", ll2, "\n")
                ll = ll2
                iter = iter + 1
                if (tmp$swap == 0) {
                    iter = 4
                }
            }
            if (trace > 0) 
                cat(swap, "\n")
            if (swap > 0) 
                rounds = 1
            if (swap == 0) 
                optNni = FALSE
        }
        epsR <- 1e-8

        if(rounds > control$maxit) opti <- FALSE
        if ((( ll1 - ll ) / ll  < control$eps ) && rounds > 2 ) #abs(ll1 - ll)
            opti <- FALSE
        rounds = rounds + 1
        ll1 = ll
    }  
}

