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
    ch = c(1L,2L,3L,4L,5L)
 elind = c(1L,2L,6L,4L,5L)
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
    ch = c(3L,4L,1L,2L,6L)
    elr = c(3L,4L,1L,2L,5L)
    eln = c(6L,4L,1L,2L,5L)
    if(x[6L]==root) el = EL[x[elr]]
    else   el = EL[x[eln]]  
    structure(list(edge = structure(c(x[pa], x[ch]), .Dim = c(5L, 2L)), 
                   tip.label = tree$tip.label, edge.length = el, Nnode = 2L), 
              .Names = c("edge", "tip.label", "edge.length", "Nnode"), class = "phylo", order = "postorder")
}




optimQuartet <- function (tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0,
          control = pml.control(epsilon = 1e-08, maxit = 10, trace=0), ...) 
{
    nTips <- length(tree$tip)
    el <- tree$edge.length
    tree$edge.length[el < 1e-08] <- 1e-08
    oldtree = tree
    k = length(w)    
    data = subset(data, tree$tip) 
    loglik = pml.fit2(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k)
    start.ll <- old.ll <- loglik 
    contrast <- attr(data, "contrast")
    contrast2 <- contrast %*% eig[[2]] 
    evi = (t(eig[[3]]) * bf)
    weight <- attr(data, "weight")
    eps = 1
    iter = 0
    
    treeP = tree
    tree = reorder(tree)
    
    child = tree$edge[, 2]
    parent = tree$edge[, 1]
    m <- max(tree$edge)
    pvec <- integer(m)
    pvec[child] <- parent
    
    EL = numeric(m)
    EL[child] = tree$edge.length
    
    n = length(tree$edge.length)  
    
    nr = as.integer(length(weight))
    nc = as.integer(length(bf))
    nco = as.integer(nrow(contrast))
    eve = eig[[2]]
    lg = k
    rootNode = getRoot(tree)         
    ScaleEPS = 1.0/4294967296.0
    anc = Ancestors(tree, 1:m, "parent")  
    anc0 = as.integer(c(0L, anc))
    
    while (eps > control$eps && iter < control$maxit) {
        blub3 <- .Call("extractScale", as.integer(rootNode), w, g, as.integer(nr), as.integer(nc), as.integer(nTips))
        rowM = apply(blub3, 1, min)       
        blub3 = (blub3-rowM) 
        blub3 = ScaleEPS ^ (blub3) 
        EL <- .Call("optE", as.integer(parent), as.integer(child), 
                    as.integer(anc0), eig, evi, EL, w, g, as.integer(nr), as.integer(nc), 
                    as.integer(nTips), as.double(contrast), 
                    as.double(contrast2), nco, blub3, data, as.double(weight), as.double(ll.0))       
        iter = iter + 1
        treeP$edge.length = EL[treeP$edge[,2]]
        newll <- pml.fit2(treeP, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k)
        
        eps = ( old.ll - newll ) / newll
        if( eps <0 ) return(list(oldtree, old.ll))
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
#    rootEdges <- attr(INDEX,"root")
#    data = getCols(data, tree$tip)
    
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
        tree0 = index2tree2(INDEX[i, ], tree, nTips+1L)
        
        ch = ei[5]
        pa = ei[6]
        print(ei)
#        browser()        

        while(pa != loli){
            tmpr = match(loli, INDEX[,5])
            treetmp = index2tree(INDEX[tmpr, ], tree, nTips+1L)
            tmpl <- pml.fit4(treetmp, data, bf=bf, g=g, w=w, eig=eig, INV=INV, ll.0=ll.0, k=k)
            print(paste("tmp:", tmpl))
            loli = anc[loli]
        }
        
        
#        if(anc[pa]==loli){
#                browser()
#                blub <- .Call("moveloli", as.integer(anc[pa]), as.integer(pa), eig, EL[pa], w, g, as.integer(nr), as.integer(nc), as.integer(nTips))
#                loli <- pa 
#      }    

            
#        browser()   
#           while(loli != pa){      
#                blub <- .Call("moveloli", as.integer(loli), as.integer(anc[loli]), eig, EL[loli], w, g, as.integer(nr), as.integer(nc), as.integer(nTips))
#                loli=anc[loli] 
#            } 
            
#        if (old.el < 1e-8) old.el <- 1e-8
#        X <- .Call("moveDad", data, as.integer(pa), as.integer(ch), eig, evi, old.el, w, g, as.integer(nr), as.integer(nc), as.integer(nTips), as.double(contrast), as.double(contrast2), nco)         
        
        tmp <- pml.fit4(tree0, data, bf=bf, g=g, w=w, eig=eig, INV=INV, ll.0=ll.0, k=k)
        print(tmp)
        loli = getRoot(tree0)
        
        
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

