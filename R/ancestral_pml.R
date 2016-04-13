#
# ancestral sequences ML
#
ancestral.pml <- function (object, type=c("ml", "bayes")) 
{
    call <- match.call()
    type <- match.arg(type)
    pt <- match.arg(type, c("ml", "bayes"))   
    tree = object$tree 
    
    INV <- object$INV
    inv <- object$inv
    
    data = getCols(object$data, tree$tip) 
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1  # max(edge)
    w = object$w
    g = object$g
    l = length(w)    
    nr <- attr(data, "nr")
    nc <- attr(data, "nc")
    dat = vector(mode = "list", length = m*l)
    result = vector(mode = "list", length = m)
    dim(dat) <- c(l,m)
    
    x = attributes(data)
    label = as.character(1:m)
    nam = tree$tip.label
    label[1:length(nam)] = nam
    x[["names"]] = label
    
    
    tmp = length(data)
    result = new2old.phyDat(data) 
    eig = object$eig
    
    bf = object$bf
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])
    for(i in 1:l)dat[i,(q + 1):m] <- .Call("LogLik2", data, P[i,], nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
    
    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    nTips = min(parent) - 1
    
    for(i in 1:l){     
        for (j in (m - 1):1) {
            if (child[j] > nTips){
                tmp2 = (dat[[i, parent[j]]]/(dat[[i,child[j]]] %*% P[[i,j]]))
                dat[[i, child[j]]] = (tmp2 %*% P[[i,j]]) * dat[[i, child[j]]]  
            }
        }
    }
    for (j in unique(parent)) {
        tmp <- matrix(0, nr, nc)
        if(inv>0) tmp = as.matrix(INV) * inv
        for(i in 1:l){  
            tmp = tmp + w[i] * dat[[i, j]]                                 
        }
        if (pt == "bayes") tmp = tmp * rep(bf, each=nr)
        tmp = tmp / rowSums(tmp)
        result[[j]] = tmp
    } 
    attributes(result) = x
    attr(result, "call") <- call
    result 
}


fast.tree  = function(tree, node){
    parent = c(node, Ancestors(tree, node))
    children = Descendants(tree, parent, 'children')
    l = sapply(children, length)
    edge = cbind(rep(parent, l), unlist(children))
    obj = list(edge=edge, Nnode=sum(l>0), tip.label=as.character(edge[is.na(match(edge[,2], edge[,1])),2]))
    class(obj) = 'phylo'
    obj
}

# schneller ???
fast.tree2  = function(tree, node){
    parent = c(node, Ancestors(tree, node))
    edge = tree$edge 
    ind = match(edge[,1], parent)
    edge=edge[which(!is.na(ind)),] 
    obj = list(edge=edge, Nnode=length(parent), tip.label=as.character(edge[is.na(match(edge[,2], edge[,1])),2]))
    class(obj) = 'phylo'
    obj
}
