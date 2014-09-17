circNetwork <- function(x, ord){
    nTips = length(ord)
    tmp = which(ord == 1)
    if(tmp!=1) ord = c(ord[tmp:nTips], ord[1:(tmp-1)])
    res = stree(nTips, tip.label = attr(x, "labels"))
    res$edge[, 2] = ord
    
    x <- phangorn:::SHORTwise(x, nTips)
    l = sapply(x, length)
    tmp <- countCycles(x, ord=ord)
    ind = which(tmp == 2 & l>1)

    X = as.matrix(x)[,ord]
    X = X[ind, ]
#    browser()
    for(i in 1: length(ind)){
        
        Vstart = ord[1]
        Vstop = ord[nTips]
        for(j in 2:nTips){
            if(X[i,j-1] < X[i,j]) Vstart = ord[j]
            if(X[i,j-1] > X[i,j]) Vstop = ord[j-1]
        } 
        g = graph(t(res$edge), directed=FALSE)
        shortestPath = get.shortest.paths(g, Vstart, Vstop)$vpath[[1]]
        
        browser()
    }
    # which splits are in circular ordering  
    
    res    
}


#addCircSplit <- function(net){}
    
    