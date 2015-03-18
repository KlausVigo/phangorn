optimEdge2 <- function (tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0,
                       control = pml.control(epsilon = 1e-08, maxit = 10, trace=0), ...) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree <- reorder(tree, "postorder") 
    nTips <- length(tree$tip)
    el <- tree$edge.length
    tree$edge.length[el < 0] <- 1e-08
    oldtree = tree
    k = length(w)    
    data = subset(data, tree$tip) 
    loglik = pml.fit(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k)
    start.ll <- old.ll <- loglik 
    contrast <- attr(data, "contrast")
    contrast2 <- contrast %*% eig[[2]] 
    evi = (t(eig[[3]]) * bf)
    weight <- attr(data, "weight")
    eps = 1
    iter = 0
    
    treeP = tree
    tree = reorder(tree)
    el <- tree$edge.length
    
    child = tree$edge[, 2]
    parent = tree$edge[, 1]
    loli = parent[1]
    
    pvec <- integer(max(tree$edge))
    pvec[child] <- parent
    
    EL = numeric(max(child))
    EL[child] = tree$edge.length
    
    nTips = min(parent) - 1
    n = length(tree$edge.length)  
    lt = length(tree$tip)
    
    nr = as.integer(length(weight))
    nc = as.integer(length(bf))
    nco = as.integer(nrow(contrast))
    eve = eig[[2]]
    lg = k
    rootNode = getRoot(tree)         

    anc = Ancestors(tree, 1:max(tree$edge), "parent")        

    while (eps > control$eps && iter < control$maxit) {
#        LL = dat[, rootNode] # sollte root sein
        for(j in 1:n) {       
            ch = child[j]
            pa = parent[j]
                       
#            blub <- .Call("extractI", as.integer(pa), w, g, as.integer(nr), as.integer(nc), as.integer(nTips))
#            if(ch>nTips)blub2 <- .Call("extractI", as.integer(ch), w, g, as.integer(nr), as.integer(nc), as.integer(nTips))
# SEXP extractScale(SEXP CH, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP NTIPS){            
#            blub3 <- .Call("extractScale", as.integer(pa), w, g, as.integer(nr), as.integer(nc), as.integer(nTips))
             
            while(loli != pa){
#                browser()
                blub <- .Call("moveloli", as.integer(loli), as.integer(anc[loli]), eig, EL[loli], w, g, as.integer(nr), as.integer(nc), as.integer(nTips))
                
                loli=anc[loli] 
            } 
            
            #            parent.dat = dat[, parent[j]]  # sollte LL sein 
            old.el = tree$edge.length[j] 

            if (old.el < 1e-8) old.el <- 1e-8
# SEXP moveDad(SEXP dlist, SEXP PA, SEXP CH, SEXP eig, SEXP EVI, SEXP EL, SEXP W, SEXP G, SEXP NR,  SEXP NC, SEXP NTIPS, SEXP CONTRAST, SEXP contrast2, SEXP NCO)
#browser()
            X <- .Call("moveDad", data, as.integer(pa), as.integer(ch), eig, evi, old.el, w, g, as.integer(nr), as.integer(nc), as.integer(nTips), as.double(contrast), as.double(contrast2), nco)                
            newEL <- .Call("FS5", eig, nc, as.double(old.el), as.double(w), as.double(g), X, as.integer(length(w)), as.integer(length(weight)), as.double(bf), as.double(weight), as.double(ll.0))
                

            el[j] = newEL[[1]]
            EL[ch] = newEL[[1]]     
            if (child[j] > nTips) {
#                LL = newEL[[3]]
                loli  = child[j]
            }
            else{ 
#                LL = newEL[[2]]
                loli  = parent[j]   
            } 
#SEXP updateLL(SEXP dlist, SEXP PA, SEXP CH, SEXP eig, SEXP EL, SEXP W, SEXP G, SEXP NR,
#              SEXP NC, SEXP NTIPS, SEXP CONTRAST, SEXP NCO)

            blub <- .Call("updateLL", data, as.integer(pa), as.integer(ch), eig, newEL[[1]], w, g,
                as.integer(nr), as.integer(nc), as.integer(nTips), as.double(contrast), nco)

            tree$edge.length = el
        }
        tree$edge.length = el
        iter = iter + 1
        
        treeP$edge.length = EL[treeP$edge[,2]]
        newll <- pml.fit(treeP, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k)
        
        eps = ( old.ll - newll ) / newll
        if( eps <0 ) return(list(oldtree, old.ll))
        oldtree = tree
        if(control$trace>1) cat(old.ll, " -> ", newll, "\n") 
        old.ll = newll
        loli = parent[1] 
        
    }
    if(control$trace>0) cat(start.ll, " -> ", newll, "\n")
    list(tree=treeP, logLik=newll, c(eps, iter))
}



