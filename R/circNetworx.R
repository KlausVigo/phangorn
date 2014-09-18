#TODO add edge.length and other parameters

circNetwork <- function(x, ord){
    nTips = length(ord)
    tmp = which(ord == 1)
    if(tmp!=1) ord = c(ord[tmp:nTips], ord[1:(tmp-1)])
    res = stree(nTips, tip.label = attr(x, "labels"))
    res$edge[, 2] = ord
    
    x <- phangorn:::SHORTwise(x, nTips)
    
    res$edge.length=NULL
    
    l = sapply(x, length)
    tmp <- countCycles(x, ord=ord)
    ind = which(tmp == 2 & l>1)

    
    
    X = as.matrix(x)[,ord]
    X = X[ind, ]
    
    fun = function(x, table) all(match(anc, table, nomatch = 0L) > 0L)
#    browser()
    for(i in 1: length(ind)){
#        y = X[i,]
        Vstart = ord[1]
        Vstop = ord[nTips]
        for(j in 2:nTips){
            if(X[i,j-1] < X[i,j]) Vstart = ord[j]
            if(X[i,j-1] > X[i,j]) Vstop = ord[j-1]
        } 
        g = graph(t(res$edge), directed=FALSE)
        sp = get.all.shortest.paths(g, Vstart, Vstop)$res #$vpath[[1]]
        parent = res$edge[,1]
        child = res$edge[,2]
        j = which(X[i,]==1)
        browser()
        anc = unique(parent[match(j, child)])
        sp = sp[[which(sapply(sp, fun, anc))]]
        sp = sp[-c(1, length(sp))]
        maxVert = max(parent)
        l = length(sp)
        newVert = (maxVert+1) : (maxVert+l)
        if(length(sp>1))newEdge = cbind(newVert[-1], newVert[-length(sp)])
        newEdge = rbind(newEdge, cbind(sp, newVert))
        res$edge = rbind(res$edge, newEdge)
        tmp = numeric(maxVert+l)
        tmp[sp] = newVert
        indT = match(j, child)
        res$edge[indT,1]  = tmp[res$edge[indT,1]]
        res$Nnode =  max(res$edge) - nTips
#        if(!is.null(res$edge.length)) res$edge.length = c(res$edge.length, rep)
    }
    # which splits are in circular ordering  
    class(res) = c("networx", "phylo")
    res    
}

#sp = allSplits(5, tr$tip)
# tr = circNetwork(sp, 1:5)
#g = graph(t(tr$edge), directed=FALSE)
#plot(g)
#addCircSplit <- function(net){}
    
    