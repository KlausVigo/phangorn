circNetwork <- function(x, ord=NULL){
    if(is.null(ord))ord = attr(x, "cycle")
    weight <- attr(x, "weights")
    
    nTips = length(ord)
    tmp = which(ord == 1)
    if(tmp!=1) ord = c(ord[tmp:nTips], ord[1:(tmp-1)])
    res = stree(nTips, tip.label = attr(x, "labels"))
    res$edge[, 2] = ord
    
    x <- phangorn:::SHORTwise(x, nTips)
    dm <- as.matrix(phangorn:::compatible2(x))
    
    res$edge.length=NULL
    
    l = sapply(x, length)
    
    tmp <- countCycles(x, ord=ord)
    ind = which(tmp == 2 & l>1)
    
    index = match(as.splits(res)[res$edge[,2]], x)
    
    ind = ind[order(l[ind])]

    X = as.matrix(x)[,ord]
    X = X[ind, ]
    
    
#    fun = function(x, table) all(match(table, x, nomatch = 0L) > 0L)

    for(k in 1: length(ind)){
        Vstart = ord[1]
        Vstop = ord[nTips]
    
        ordStart = 1
        ordStop = nTips
        for(j in 2:nTips){
   
            if(X[k,j-1] < X[k,j]){ 
                Vstart = ord[j]
                ordStart = j                   
            }                       
            if(X[k,j-1] > X[k,j]){ 
                Vstop = ord[j-1]
                ordStop = j-1   
            }    
        } 
        fromTo <- ordStart:ordStop
        if(ordStart>ordStop) fromTo <- c(ordStart:nTips, 1:ordStop)
        fromTo = ord[fromTo] 

        if(min(res$edge)==0) browser()      
        g = graph(t(res$edge), directed=FALSE)
        g2= as.igraph.phylo(res)
 #       print(k)
 #       plot(g)
        
        sp2 = NULL
        sp0 = NULL
        for(i in 2:length(fromTo)){
            sptmp = get.shortest.paths(g, fromTo[i-1], fromTo[i], 
                output=c("epath"))$epath[[1]]
            sp2 = c(sp2, sptmp[-c(1, length(sptmp))])
            sp0 = c(sp0, sptmp)
        }
        sp0 = unique(sp0)
        print(dm[index[sp2], ind[k]])
 #       browser() 
        blub = which(dm[index[sp2], ind[k]]>0)
        sp2 = sp2[blub]
        
        if(length(sp2)==0){
            
            edge1 = unique(as.vector(res$edge[sp0,]))
            edge2 = as.vector(res$edge[-sp0,])
            asdf = edge1 %in% edge2
            sp = edge1[asdf]
 
    }
#        sp3 = NULL 
#        if(length(sp2>0))
#        else 
if(length(sp2)>0) {   sp = unique(as.vector(t(res$edge[sp2,])))
#browser()
}
#        print(sp)
 
#        sp = get.all.shortest.paths(g, Vstart, Vstop)$res #$vpath[[1]]  
        parent = res$edge[,1]
        child = res$edge[,2]

        
        j = ord[which(X[k,]==1)]
        anc = unique(parent[match(j, child)])
#        indi = which(sapply(sp, fun, anc))
#        if(length(indi > 1)) indi = indi[1]
#        sp = sp[[indi]]
#        sp = sp[-c(1, length(sp))]
        maxVert = max(parent)
        l = length(sp)

        newVert = (maxVert+1) : (maxVert+l)
#        for(i in 1:l) res$edge[sp0,][res$edge[sp0,]==sp[i]] = newVert[i]                
#        print(sp3)
#        print(res$edge[sp2,])
        sp01 = setdiff(sp0, sp2)
        for(i in 1:l) res$edge[sp01,][res$edge[sp01,]==sp[i]] = newVert[i]
        
        newindex = rep(ind[k], l)        
        if(length(sp)>1)newindex = c(index[sp2], newindex)
        index = c(index, newindex)        
        # verbinde vertices
        newEdge = matrix(cbind(sp, newVert), ncol=2) 
        if(length(sp)>1){
            # copy edges
            qwer = match(as.vector(res$edge[sp2,]), sp)
            newEdge = rbind(matrix(newVert[qwer], ncol=2), newEdge)
        }

       
        res$edge = rbind(res$edge, newEdge)
if(length(index)!= nrow(res$edge))browser()
        tmp = numeric(maxVert+l)
        tmp[sp] = newVert

 #       indT = match(j, child)
 #       res$edge[indT,1]  = tmp[res$edge[indT,1]]

        res$Nnode =  max(res$edge) - nTips
#browser()
#        g = graph(t(res$edge), directed=FALSE)
#        plot(g)     
    }
    res$split = index 
#print(index)
    res$edge.length = weight[index] 
    class(res) = c("networx", "phylo")
    res    
}

#sp = allSplits(5, tr$tip)
# tr = circNetwork(sp, 1:5)
#g = graph(t(tr$edge), directed=FALSE)
#plot(g)
#addCircSplit <- function(net){}
    



    