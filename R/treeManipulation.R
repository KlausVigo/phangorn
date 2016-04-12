#
# tree manipulation
# 

# from coalescenMCMC
getIndexEdge <- function(tip, edge)
    ## 'integer(1)' mustn't be substituted by '0L' except if 'DUP = TRUE':
    .C("get_single_index_integer", as.integer(edge[, 2L]),
       as.integer(tip), integer(1L), PACKAGE = "phangorn",
       NAOK = TRUE)[[3L]]

getIndexEdge2 <- function(node, edge)
    .C("get_two_index_integer", as.integer(edge[, 1L]),
       as.integer(node), integer(2L), PACKAGE = "phangorn",
       NAOK = TRUE)[[3L]]

# no checks for postorder
getRoot <- function (tree) 
{
    if(!is.null(attr(tree, "order")) && attr(tree, "order") == 
           "postorder"){
        return(tree$edge[nrow(tree$edge), 1])
    }    
    res = unique(tree$edge[, 1][!match(tree$edge[, 1], tree$edge[, 2], 0)])
    if (length(res) == 1) 
        return(res)
    else stop("There are apparently two root edges in your tree")
}


# renames root node 
reroot <-  function (tree, node) 
{
    anc = Ancestors(tree, node, "all")
    l = length(anc)
    if(is.na(match(node,tree$edge[,1])))stop("node not in tree")
    if(l==0)return(tree)
    ind = match(c(node, anc[-l]), tree$edge[, 2])
    tree$edge[ind, c(1, 2)] = tree$edge[ind, c(2, 1)]
    root = anc[l]
    tree$edge[tree$edge == root] = 0L
    tree$edge[tree$edge == node] = root
    tree$edge[tree$edge == 0L] = node
# needed for unrooted trees    
    tree <- collapse.singles(tree)
    attr(tree, "order") <- NULL
    reorder(tree, "postorder")
}


reroot2 <- function(tree, node) {
    if(node==getRoot(tree)) return(tree)
    anc = Ancestors(tree, node, "all")
    l = length(anc)
    ind = match(c(node, anc[-l]), tree$edge[, 2])
    tree$edge[ind, c(1, 2)] = tree$edge[ind, c(2, 1)]
    reorderPruning(tree)   
}    


changeEdge = function (tree, swap, edge = NULL, edge.length = NULL) 
{
    attr(tree, "order") = NULL
    child <- tree$edge[, 2]
    tmp = numeric(max(child))
    tmp[child] = 1:length(child)
    tree$edge[tmp[swap[1]], 2] = swap[2]
    tree$edge[tmp[swap[2]], 2] = swap[1]
    if (!is.null(edge)) {
        tree$edge.length[tmp[edge]] = edge.length
    }
    reorder(tree, "postorder")
}


changeEdgeLength = function (tree, edge, edge.length) 
{
    tree$edge.length[match(edge, tree$edge[,2])] = edge.length
    tree
}


# O(n) statt O(n^2) Speicher und Geschwindigkeit
midpoint <- function(tree, node.labels = "support"){
# distance from node to root
node2root <- function(x){
    x = reorder(x, "postorder")
    el = numeric(max(x$edge))   
    parents <- x$edge[, 1]
    child <- x$edge[, 2]
    el[child] = x$edge.length  
    l = length(parents)
    res <- numeric(max(x$edge))
    for(i in l:1){            
          res[child[i]] = el[child[i]]  + res[parents[i]]
     } 
     res
}
    oldtree <- tree
    tree = unroot(tree)   
    nTips = length(tree$tip)
    maxD1 = node2root(tree)[1:nTips] 
    ind = which.max(maxD1)
    tmproot = Ancestors(tree, ind, "parent")
#    tree = reroot(tree, tmproot)
    nTips  = length(tree$tip.label)
    if(tmproot>nTips) tree = root(tree, node=tmproot)
    else  tree = root(tree, tmproot)   
    el = numeric(max(tree$edge))
    el[tree$edge[,2]]=tree$edge.length  
    maxdm = el[ind]
    tree$edge.length[tree$edge[,2]==ind] = 0 
    maxD1 = node2root(tree)[1:nTips]  
    tree$edge.length[tree$edge[,2]==ind] = maxdm 
    ind = c(ind, which.max(maxD1) ) 
    maxdm = maxdm + maxD1[ind[2]]    
    rn = max(tree$edge)+1
    edge = tree$edge
    el = tree$edge.length
    children = tree$edge[,2]
    left = match(ind[1], children)
    tmp = Ancestors(tree, ind[2], "all")
    tmp= c(ind[2], tmp[-length(tmp)]) 
    right = match(tmp, children)
    if(el[left]>= (maxdm/2)){
         edge = rbind(edge, c(rn, ind[1]))       
         edge[left,2] = rn 
         el[left] = el[left] - (maxdm/2)
         el = c(el, maxdm/2) 
    }
    else{
        sel = cumsum(el[right]) 
        i = which(sel>(maxdm/2))[1]
        edge = rbind(edge, c(rn, tmp[i]))       
        edge[right[i],2] = rn  
        eltmp =  sel[i] - (maxdm/2)
        el = c(el, el[right[i]] - eltmp)
        el[right[i]] = eltmp
    }
    tree$edge.length = el
    tree$edge=edge
    tree$Nnode  = tree$Nnode+1
    attr(tree, "order") <- NULL
    tree <- reroot(tree, rn) 
    if(!is.null(tree$node.label)){
        node.label = tree$node.label
        tmp = node.label[1]
        node.label[1] = node.label[rn-nTips]
        node.label[rn-nTips] = tmp
        node.label[is.na(node.label)] = ""
        tree$node.label = node.label
    }
    attr(tree, "order") <- NULL
    tree <- reorder(tree, "postorder")
#    tree <- reorder(reroot(tree, rn), "postorder")
#    if(!is.null(oldtree$node.label))tree <- addConfidences.phylo(tree, oldtree)
    if(!is.null(oldtree$node.label)){
        type <- match.arg(node.labels, c("support", "label", "delete"))
        if(type=="support") tree <- addConfidences.phylo(tree, oldtree)
        if(type=="delete") tree$node.label <- NULL
    }
    tree 
}


pruneTree = function(tree, ..., FUN = ">="){
     if(is.null(tree$node)) stop("no node labels")
     if(is.rooted(tree)) tree = unroot(tree)
     m = max(tree$edge)
     nTips = length(tree$tip)
     bs = rep(TRUE, m)
     bs[ (nTips+1) : m] = sapply(as.numeric(as.character(tree$node)), FUN,...)    
     tree$edge.length[!bs[tree$edge[,2]]] = 0
   
     reorder(di2multi(tree), "postorder")
}


# requires postorder
# works fine with fit.fitch  
# for internal use in fitch.spr  
# pos statt i      
dropTip <- function(x, i, check.binary=FALSE, check.root=TRUE){
    edge <- x$edge
#    edge1 <- edge[,1]
#    edge2 <- edge[,2]
    root <- getRoot(x)
    ch <- which(edge[,2] == i)
    pa <- edge[ch,1] 
    edge = edge[-ch,]
    ind <- which(edge[,1] == pa) 
    if(root == pa){
        if(length(ind)==1){
            edge = edge[-ind,]
            x$Nnode=x$Nnode-1L
        }
        if(length(ind)==2){
            n = dim(edge)[1]
            newroot = edge[n-2L,1]
            newedge = edge[ind,2] 
            if(newedge[1]==newroot)edge[n-1,] <- newedge
            else edge[n-1,] <- newedge[2:1]
            edge = edge[-n,]   
            x$Nnode=x$Nnode-1L
            edge[edge==newroot] = root
            pa <- newroot
        }
    # todo handle unrooted trees  
    }
    else{
        nind <- which(edge[,2] == pa)         
    # normal binary case
        if(length(ind)==1){
            edge[nind,2] = edge[ind,2]
            edge <- edge[-ind,]
            x$Nnode <- x$Nnode-1L           
        }  
    }
    #
    edge[edge>pa]  = edge[edge>pa] -1L 
    x$edge <- edge
    x
}

# kind of works well too
dropTip2 <- function(x, i, check.binary=FALSE, check.root=TRUE){
  edge <- x$edge
  root <- getRoot(x)
  ch <- which(edge[,2] == i)
  pa <- edge[ch,1] 
  edge = edge[-ch,]
  ind <- which(edge[,1] == pa) 
  if(root == pa){
    if(length(ind)==1){
      edge = edge[-ind,]
      x$Nnode=x$Nnode-1L
    }
    if(length(ind)==2){
      n = dim(edge)[1]
      newroot = edge[n-2L,1]
      newedge = edge[ind,2] 
      if(newedge[1]==newroot)edge[n-1,] <- newedge
      else edge[n-1,] <- newedge[2:1]
      edge = edge[-n,]   
      x$Nnode=x$Nnode-1L
      edge[edge==newroot] = root
      pa <- newroot
    }
    # todo handle unrooted trees  
  }
  else{
    nind <- which(edge[,2] == pa)         
    # normal binary case
    if(length(ind)==1){
      edge[nind,2] = edge[ind,2]
      edge <- edge[-ind,]
      x$Nnode <- x$Nnode-1L           
    }  
  }
  #
#  edge[edge>pa]  = edge[edge>pa] -1L 
  x$edge <- edge
  x
}


# like drop tip and returns two trees, 
# to be used in fitch.spr
#  ch = allKids(edge, nTips)
descAll = function (x, node, nTips, ch) 
{
    edge = x[,1]
    m = max(x)
    isInternal = logical(m)
    isInternal[(nTips+1):m] = TRUE
    desc = function(node, isInternal) {
        if (!isInternal[node]) return(node)
        res = NULL
        while (length(node) > 0) {
            tmp = unlist(ch[node])
            res = c(res, tmp)
            node = tmp[isInternal[tmp]]
        }
        res
    }
    desc(node, isInternal)
}  


dropNode <- function(x, i, check.binary=FALSE, check.root=TRUE, all.ch=NULL){
  edge <- x$edge
  root <- getRoot(x)
  ch <- which(edge[,2] == i)
#  getIndexEdge(tip, edge) 
  
  nTips <- length(x$tip.label)    
  pa <- edge[ch,1] 
  if(i>nTips){
#    kids <- Descendants(x, i, "all")
    if(is.null(all.ch)) all.ch=allChildren(x)  
    kids <- descAll(edge, i, nTips, all.ch)  
    ind <- match(kids,edge[,2])
    edge2 <- edge[sort(ind),]            
    edge <- edge[-c(ch, ind),]
  }    
  else edge = edge[-ch,]
  if(nrow(edge)<3)return(NULL)  
  ind <- which(edge[,1] == pa) 
  if(root == pa){
    if(length(ind)==1){
      edge = edge[-ind,]
      x$Nnode=x$Nnode-1L
    }
    if(length(ind)==2){
      n = dim(edge)[1]
      newroot = edge[n-2L,1]
      newedge = edge[ind,2] 
      if(newedge[1]==newroot)edge[n-1,] <- newedge
      else edge[n-1,] <- newedge[2:1]
      edge = edge[-n,]   
      x$Nnode=length(unique(edge[,1]))
      edge[edge==newroot] = root
      pa <- newroot
    }
    # todo handle unrooted trees  
  }
  else{
    nind <- which(edge[,2] == pa)         
    # normal binary case
    if(length(ind)==1){
      edge[nind,2] = edge[ind,2]
      edge <- edge[-ind,]
      x$Nnode <- length(unique(edge[,1]))          
    }  
  }
  #
#  edge[edge>pa]  = edge[edge>pa] -1L 
  x$edge <- edge
  y <- x
  y$edge <- edge2
  y$Nnode <- length(unique(edge2[,1]))
  list(x, y, pa)
}



dropNodeNew <- function(edge, i, nTips, check.binary=FALSE, check.root=TRUE){
    root <- edge[nrow(edge),2]
    ch <- which(edge[,2] == i)
    pa <- edge[ch,1]
    edge2=NULL
    
    # einfachere allChildren Variante 2*schneller
    allKids = function (edge, nTips) 
    {
        parent = edge[, 1]
        children = edge[, 2]
        .Call("AllChildren", as.integer(children), as.integer(parent), as.integer(max(edge)), PACKAGE = "phangorn")
    }
    
    descAll = function (edge, node, nTips) 
    {
        ch = allKids(edge, nTips)
        isInternal = logical(max(edge))
        isInternal[unique(edge[, 1])] = TRUE
        desc = function(node, isInternal) {
            if (!isInternal[node]) return(node)
            res = NULL
            while (length(node) > 0) {
                tmp = unlist(ch[node])
                res = c(res, tmp)
                node = tmp[isInternal[tmp]]
            }
            res
        }
        desc(node, isInternal)
    }    
    
    if(i>nTips){
        kids <- descAll(edge, i, nTips)
        ind <- match(kids,edge[,2])
        edge2 <- edge[sort(ind),]            
        edge <- edge[-c(ch, ind),]
    }    
    else edge = edge[-ch,]
    if(nrow(edge)<3)return(NULL)  
    ind <- which(edge[,1] == pa) 
    if(root == pa){
        if(length(ind)==1){
            edge = edge[-ind,]
        }
        if(length(ind)==2){
            n = dim(edge)[1]
            newroot = edge[n-2L,1]
            newedge = edge[ind,2] 
            if(newedge[1]==newroot)edge[n-1,] <- newedge
            else edge[n-1,] <- newedge[2:1]
            edge = edge[-n,]   
            edge[edge==newroot] = root
            pa <- newroot
        }
        # todo handle unrooted trees  
    }
    else{
        nind <- which(edge[,2] == pa)         
        # normal binary case
        if(length(ind)==1){
            edge[nind,2] = edge[ind,2]
            edge <- edge[-ind,]          
        }  
    }
    #
    #  edge[edge>pa]  = edge[edge>pa] -1L 
    list(edge, edge2, pa)
}


dropTipNew <- function(edge, i, nTips, check.binary=FALSE, check.root=TRUE){
    root <- edge[nrow(edge),2]
    ch <- which(edge[,2] == i)
    pa <- edge[ch,1] 
    edge = edge[-ch,]
    ind <- which(edge[,1] == pa) 
    if(root == pa){
        if(length(ind)==1){
            edge = edge[-ind,]
        }
        if(length(ind)==2){
            n = dim(edge)[1]
            newroot = edge[n-2L,1]
            newedge = edge[ind,2] 
            if(newedge[1]==newroot)edge[n-1,] <- newedge
            else edge[n-1,] <- newedge[2:1]
            edge = edge[-n,]   
            edge[edge==newroot] = root
            pa <- newroot
        }
        # todo handle unrooted trees  
    }
    else{
        nind <- which(edge[,2] == pa)         
        # normal binary case
        if(length(ind)==1){
            edge[nind,2] = edge[ind,2]
            edge <- edge[-ind,]       
        }  
    }
    #
    edge[edge>pa]  = edge[edge>pa] -1L 
    edge
}


# postorder remained tip in 1:nTips
addOne <- function (tree, tip, i){
    edge = tree$edge
    parent = edge[,1]
    l = dim(edge)[1]
    m = max(edge)+1L 
    p = edge[i,1]
    k = edge[i,2] 
    edge[i, 2] = m
    ind = match(p, parent)
    if(ind==1) edge = rbind(matrix(c(m,m,k,tip), 2, 2), edge)
    else edge = rbind(edge[1:(ind-1), ], matrix(c(m,m,k,tip), 2, 2), edge[ind:l, ])  
    tree$edge = edge 
    tree$Nnode = tree$Nnode+1
    tree
}         


addOneTree <- function (tree, subtree, i, node){
  edge = tree$edge
  parent = edge[,1]
  l = dim(edge)[1]
  m = node #max(edge)+1L 
  p = edge[i,1]
  k = edge[i,2] 
  edge[i, 2] = m
  edge2 = subtree$edge
  ind = match(p, parent)
  r2 = edge2[nrow(edge2),1]
  if(ind==1) edge = rbind(edge2, matrix(c(m,m,r2,k), 2, 2), edge)
  else edge = rbind(edge[1:(ind-1), ], edge2, matrix(c(m,m,r2,k), 2, 2), edge[ind:l, ])  
  tree$edge = edge 
  tree$Nnode = tree$Nnode + subtree$Nnode + 1L
  attr(tree, "order") = NULL
  tips1 = as.integer(length(tree$tip) + 1L)
  tmproot = getRoot(tree)
  if(tmproot!=tips1){
      tree$edge[tree$edge==tmproot] = 0L
      tree$edge[tree$edge==tips1] = tmproot
      tree$edge[tree$edge==0L] = tips1    
  }
  tree <- reorder(tree, "postorder")
  if(tmproot!=tips1) tree <- unroot(tree)
  tree
}         


reorderPruning <- function (x, ...) 
{
    parents <- as.integer(x$edge[, 1])
    child <- as.integer(x$edge[, 2])
    root <- as.integer(parents[!match(parents, child, 0)][1])  # unique out
    if (length(root) > 2) 
        stop("more than 1 root found")
    n = length(parents)    
    m = max(x$edge)  # edge  parents 
    neworder = .C("C_reorder", parents, child, as.integer(n), as.integer(m), integer(n), as.integer(root-1L), PACKAGE = "phangorn")[[5]]    
    x$edge = x$edge[neworder,]
    x$edge.length = x$edge.length[neworder]
    attr(x, "order") <- "pruningwise"
    x
}


add.tip <- function(phy, n, edgeLength=NULL, tip=""){ 
     ind <- which(phy$edge[,2] == n)
     phy <- new2old.phylo(phy) 
     edge <- matrix(as.numeric(phy$edge),ncol=2)
     k <- min(edge)
     l <- max(edge)
     phy$edge <- rbind(phy$edge, c(k-1,phy$edge[ind,2]))
     phy$edge <- rbind(phy$edge, c(k-1,l+1))
     phy$edge[ind,2] = k-1 
     phy$edge.length[ind] = edgeLength[1]
     phy$edge.length <- c(phy$edge.length, edgeLength[-1])
     phy$tip.label <- c(phy$tip.label, tip) 
     phy <- old2new.phylo(phy)
     phy <- reorder(phy, "postorder") 
     phy
}


nnin <- function (tree, n) 
{
    attr(tree, "order") = NULL 
    tree1 = tree
    tree2 = tree
    edge = matrix(tree$edge, ncol = 2)
    parent = edge[, 1]
    child = tree$edge[, 2]
    k = min(parent) - 1
    ind = which(child > k)[n]
    if(is.na(ind))return(NULL)
    p1 = parent[ind]
    p2 = child[ind]
    ind1 = which(parent == p1)
    ind1 = ind1[ind1 != ind][1]
    ind2 = which(parent == p2)
    e1 = child[ind1]
    e2 = child[ind2[1]]
    e3 = child[ind2[2]]
    tree1$edge[ind1, 2] = e2
    tree1$edge[ind2[1], 2] = e1
    tree2$edge[ind1, 2] = e3
    tree2$edge[ind2[2], 2] = e1
    if(!is.null(tree$edge.length)){
        tree1$edge.length[c(ind1, ind2[1])] = tree$edge.length[c(ind2[1] ,ind1)]
        tree2$edge.length[c(ind1, ind2[2])] = tree$edge.length[c(ind2[2] ,ind1)]
        }
    tree1 <- reorder(tree1, "postorder")  
    tree2 <- reorder(tree2, "postorder")  
#    tree1$tip.label <- tree2$tip.label <- NULL    
    result = list(tree1, tree2)
    result
} 


nni <- function (tree) 
{
    tip.label <- tree$tip.label
    attr(tree, "order") = NULL
    k = min(tree$edge[, 1]) - 1
    n = sum(tree$edge[, 2] > k)
    result = vector("list", 2*n)
    l=1
    for (i in 1:n) {
          tmp = nnin(tree, i)
          tmp[[1]]$tip.label <- tmp[[2]]$tip.label <- NULL
          result[c(l, l+1)] = tmp
          l = l + 2
          }
    attr(result, "TipLabel") <- tip.label
    class(result) <- "multiPhylo"
    result
}


allTrees <- function (n, rooted = FALSE, tip.label = NULL) 
{
	n <- as.integer(n)  
    nt <- as.integer(round(dfactorial(2 * (n + rooted) - 5))) 
    if ((n + rooted) > 10) {
        nt <- dfactorial(2 * (n + rooted) - 5)
        stop("That would generate ", round(nt), " trees, and take up more than ", 
            round(nt/1000), " MB of memory!")
    }
    if (n < 2) {
        stop("A tree must have at least two taxa.")
    }
    if (!rooted && n == 2) {
        stop("An unrooted tree must have at least three taxa.")
    }

    if (rooted) {
        edge <- matrix(NA, 2*n-2, 2)
        edge[1:2,] <- c(n+1L, n+1L, 1L, 2L)
    }
    else {
        edge <- matrix(NA,  2*n-3, 2)
        edge[1:3,] <- c(n+1L, n+1L, n+1L, 1L, 2L, 3L)
    }
    edges <- list()
    edges[[1]] <- edge

    m=1     
    nedge = 1
    trees <- vector("list", nt)
    if ((n + rooted) > 3) {
        i = 3L  + (!rooted)    
        pa = n + 2L
        nr = 2L + (!rooted)
        while(i < (n+1L)){
            nedge = nedge+2
            m2 = m * nedge 
            
            newedges <- vector("list", m2)
            for (j in 1:m) {
                edge <- edges[[j]]
                l <- nr  # nrow(edge)

                    edgeA <- edge
                    edgeB <- edge

                    for (k in 1L:l) {
                       edge = edgeA
                       node <- edge[k, 1]
                       edge[k, 1] <- pa             
                       edge[l + 1, ] <- c(pa, i)
                       edge[l + 2, ] <- c(node, pa)

                       newedges[[(j - 1) * (l + rooted) + k]] <- edge
                       }

                if(rooted) { 
                  edgeB[] <- as.integer(sub(n+1L, pa, edgeB))
                  edge = edgeB
                  edge[l + 1, ] <- c(n+1L, i)
                  edge[l + 2, ] <- c(n+1L, pa) 
                  newedges[[j * (l + 1)]] <- edge
                }
            } # end for 
            edges <- newedges
            m=m2
            i = i + 1L
            pa = pa + 1L  
            nr = nr + 2L 
        } # end for m
    } # end if
    for (x in 1:m) {
        tree <- list(edge = edges[[x]])
        tree$Nnode <- n - 2L + rooted
        class(tree) <- "phylo"       
        trees[[x]] <- reorderPruning(tree)
    }
    attr(trees, "TipLabel") <- if (is.null(tip.label)) 
        paste("t", 1:n, sep = "")
    else tip.label
    class(trees) <- "multiPhylo"
    trees
}



dn <- function (x){
    if (!is.binary.tree(x) ) 
        x <- multi2di(x, random = FALSE)  
    x = reroot2(x, 1)       
    n <- length(x$tip.label)
    n.node <- x$Nnode
    N <- n + n.node
    x <- reorderPruning(x)
    res <- matrix(NA, N, N)
    res[cbind(1:N, 1:N)] <- 0
    res[x$edge] <- res[x$edge[, 2:1]] <- 1
    for (i in seq(from = 1, by = 2, length.out = n.node)) {
        j <- i + 1
        anc <- x$edge[i, 1]
        des1 <- x$edge[i, 2]
        des2 <- x$edge[j, 2]
        if (des1 > n) 
            des1 <- which(!is.na(res[des1, ]))
        if (des2 > n) 
            des2 <- which(!is.na(res[des2, ]))
        for (y in des1) res[y, des2] <- res[des2, y] <- res[anc, 
            y] + res[anc, des2]
        if (anc != 1) {
            ind <- which(x$edge[, 2] == anc)
            nod <- x$edge[ind, 1]
            l <- length(ind)
            res[des2, nod] <- res[nod, des2] <- res[anc, des2] + 
                l
            res[des1, nod] <- res[nod, des1] <- res[anc, des1] + 
                l
        }
    }
    dimnames(res)[1:2] <- list(1:N)
    res
}


# replace with dist.nodes
dn <- function(x){
#  if (!is.binary.tree(x) ) x <- multi2di(x, random = FALSE)
  if(is.null(x$edge.length))x$edge.length=rep(1,nrow(x$edge))
  else x$edge.length[] = 1
  dist.nodes(x)
}


rSPR = function (tree, moves = 1, n = length(moves), k=NULL) 
{
    if (n == 1) {
        trees = tree
        for (i in 1:moves) trees = kSPR(trees, k=k)
    }
    else {
        trees = vector("list", n)
        if(length(moves)==1) moves = rep(moves, n)
        
        for (j in 1:n) {
            tmp = tree
            if(moves[j]>0){
               for (i in 1:moves[j]) tmp = kSPR(tmp, k=k)
            }
            tmp$tip.label = NULL
            trees[[j]] = tmp
        }
        attr(trees, "TipLabel") <- tree$tip.label
        class(trees) <- "multiPhylo"
    }
    trees
}


kSPR = function(tree, k=NULL){  
    l <- length(tree$tip.label)
    root= getRoot(tree)
    distN = dn(tree)[-c(1:l), -c(1:l)]
    distN[upper.tri(distN)]=Inf
    dN = distN[lower.tri(distN)]
    tab = table(dN) 
    tab[1] = tab[1] * 2 
    tab[-1] = tab[-1] * 8   
    if(is.null(k)) k = 1:length(tab)
    k = na.omit((1:length(tab))[k])
    if(length(k)>1)k = sample((1:length(tab))[k], 1, prob=tab[k] / sum(tab[k]) )
    if(k==1) return(rNNI(tree, 1, 1))
    index = which(distN==k, arr.ind=TRUE) + l
    m = dim(index)[1]
    if(m==0)stop("k is chosen too big")
    ind = index[sample(m, 1),]
    s1 = sample(c(1,2),1) 
    if(s1==1)res = (oneOf4(tree, ind[1], ind[2], sample(c(1,2),1), sample(c(1,2),1)))
    if(s1==2)res = (oneOf4(tree, ind[2], ind[1], sample(c(1,2),1), sample(c(1,2),1))) 
    res=reroot2(res, root)
    reorderPruning(res)    
}


oneOf4 = function(tree, ind1, ind2, from=1, to=1){
    if (!is.binary.tree(tree)) 
        stop("Sorry, trees must be binary!")        
    tree=reroot2(tree, ind2)
    trees = vector('list', 8)
    kids1 = Children(tree, ind1)
    anc = Ancestors(tree, ind1, "all")
    l = length(anc)
    kids2 = Children(tree, ind2)
    kids2 = kids2[kids2!=anc[l-1]]

    child = tree$edge[,2]
    tmp = numeric(max(tree$edge))
    tmp[child] = 1:length(child)

    edge = tree$edge
    edge[tmp[kids1[-from]],1] = Ancestors(tree, ind1, "parent")
    edge[tmp[kids2[to]],1] = ind1
    edge[tmp[ind1]] = ind2
    tree$edge=edge
    reorderPruning(tree)   
}


# faster than kSPR
rSPR_Old <- function(tree, moves=1, n=1){
    k=length(tree$edge[,1])
    if(n==1){
        trees = tree
        for(i in 1:moves) trees = sprMove(trees,sample(k,1))  
    }  
    else{
        trees = vector("list", n)
        for(j in 1:n){
            tmp = tree 
            for(i in 1:moves) tmp = sprMove(tmp,sample(k,1))
            tmp$tip.label=NULL
            trees[[j]] = tmp
        }
        attr(trees, "TipLabel") <- tree$tip.label
        class(trees) <- "multiPhylo"   
    }
    trees
}


sprMove <- function(tree, m){
    if (is.rooted(tree)) tree <- unroot(tree)
    #stop("Sorry trees must be unrooted")
    if(!is.binary.tree(tree))stop("Sorry trees must be binary!")

    reroot2 <- function(tree, node){
        anc = Ancestors(tree, node, "all")
        l = length(anc)
        ind = match(c(node, anc[-l]), tree$edge[,2])
        tree$edge[ind,c(1,2)] = tree$edge[ind,c(2,1)]
        tree    
    }    
    changeEdge <- function(tree, new, old){
        tree$edge[tree$edge==old] = 0L
        tree$edge[tree$edge==new] = old
        tree$edge[tree$edge==0L] = new
    # needed for unrooted trees
        tree <- collapse.singles(tree)
        tree          
    }

    edge = tree$edge    
    k = max(edge)
    nTips = length(tree$tip)
    nEdges = 2*nTips-3 
    if(m > nEdges) stop("m to big")

    parent = edge[,1]
    child = edge[,2]
    pv = integer(k)      
    pv[child] = parent
    cv = list()
    for(i in unique(parent)) cv[[i]] = child[parent==i]
    bp = bip(tree)
    root <- parent[!match(parent, child, 0)][1]    
       
    ch = child[m]
    pa = parent[m] 

    candidates = !logical(k)
    candidates[root] = FALSE     
    candidates[cv[[ch]]] = FALSE
    candidates[cv[[pa]]] = FALSE
    candidates[pv[pa]] = FALSE
    candidates[pa] = FALSE

    ind = which(candidates)
    l = sample(ind,1)

    cr=FALSE 

    if(!any(is.na(match(bp[[l]], bp[[ch]]))) ){
        
        newroot = cv[[ch]] #[ 1]
        newroot = newroot[newroot>nTips][1]
        tree <- reroot2(tree, newroot)
        edge = tree$edge
        parent = tree$edge[,1]
        child = tree$edge[,2]
        pv = integer(k)      
        pv[child] = parent
        cv = list()
        for(i in unique(parent)) cv[[i]] = child[parent==i]
        
        tmp = pa
        pa=ch
        ch=tmp
        cr = TRUE
    }

    if(pa==root){
        cp = cv[[pa]]
        newroot = cp[cp!=ch]
        
        newroot = newroot[newroot>nTips][1]
#        if(length(newroot)==0)browser()
        #!newroot = cp[cp>nTips][1]
        tree = reroot2(tree, newroot)
        edge = tree$edge
        parent = tree$edge[,1]
        child = tree$edge[,2]
        pv = integer(k)      
        pv[child] = parent
        cv = list()
        for(i in unique(parent)) cv[[i]] = child[parent==i]
        
        cr = TRUE 
    }

    el = tree$edge.length
    cp = cv[[pa]]
    sib = cp[cp!=ch]

    edge[child==l,1] = pa
    edge[child==pa,1] = pv[l]  
    edge[child==sib,1] = pv[pa]

    el[child==sib] = el[child==sib] + el[child==pa]
    el[child==l] = el[child==l] / 2
    el[child==pa] = el[child==l]   

    tree$edge=edge
    tree$edge.length = el
    if(cr) tree <- changeEdge(tree,root,newroot)    
    tree <- reorder(tree, "postorder") 
    tree    
}
 

rNNI <- function (tree, moves = 1, n = length(moves)) 
{
    k = length(na.omit(match(tree$edge[, 2], tree$edge[, 1])))
    if (n == 1) {
        trees = tree
        if(moves>0){
            for (i in 1:moves) trees = nnin(trees, sample(k, 1))[[sample(2,1)]]
        }
        trees$tip.label <- tree$tip.label
    }
    else {
        trees = vector("list", n)
        if(length(moves)==1) moves = rep(moves, n)
        for (j in 1:n) {
            tmp = tree
            if(moves[j]>0){
                for (i in 1:moves[j]) tmp = nnin(tmp, sample(k, 1))[[sample(2,1)]]
            }
            tmp$tip.label = NULL
            trees[[j]] = tmp
        }
        attr(trees, "TipLabel") <- tree$tip.label
        class(trees) <- "multiPhylo"
    }
    trees
}


#
# some generic tree functions
#
allAncestors <- function(x){
    x = reorder(x, "postorder")
    parents <- x$edge[, 1]
    child <- x$edge[, 2]
    l = length(parents)
    res <- vector("list",max(x$edge))
    for(i in l:1){
          pa = parents[i]  
          res[[child[i]]] = c(pa, res[[pa]])
     } 
     res
}


Ancestors <- function (x, node, type = c("all", "parent")) 
{
    parents <- x$edge[, 1]
    child <- x$edge[, 2]
    pvector <- numeric(max(x$edge)) # parents
    pvector[child] <- parents    
    type <- match.arg(type)
    if (type == "parent") 
        return(pvector[node])
    anc <- function(pvector, node){
        res <- numeric(0)
        repeat {
            anc <- pvector[node]
            if (anc == 0) break
            res <- c(res, anc)
            node <- anc
        }
    res
    }
    if(length(node)==1) return(anc(pvector, node))
    else allAncestors(x)[node]
}


allChildren <- function(x){
   l = length(x$tip) 
   if(l<20){
       parent = x$edge[,1]
       children = x$edge[,2]
       res = vector("list", max(x$edge))
       for(i in 1:length(parent)) res[[parent[i]]] = c(res[[parent[i]]], children[i])
       return(res)
   }
   else{
       if (is.null(attr(x, "order")) || attr(x, "order") == "cladewise") 
           x <- reorder(x, "postorder")
       parent = x$edge[,1]
       children = x$edge[,2]
       res <- .Call("AllChildren", as.integer(children), as.integer(parent), as.integer(max(x$edge))) # , PACKAGE="phangorn"
       return(res)
   }
}


Children <- function(x, node){
   if(length(node)==1)return(x$edge[x$edge[,1]==node,2])
   allChildren(x)[node]
}


Descendants = function(x, node, type=c("tips","children","all")){
  type <- match.arg(type)
  if(type=="children") return(Children(x, node))
  if(type=="tips") return(bip(x)[node])
  ch = allChildren(x) # out of the loop
  isInternal = logical(max(x$edge))
  isInternal[ unique(x$edge[,1]) ] =TRUE  
  desc = function(node, isInternal){     
    if(!isInternal[node])return(node)   
    res = NULL
    while(length(node)>0){
      tmp = unlist(ch[node])
      res = c(res, tmp)
      node = tmp[isInternal[tmp]]
    }
    res
  }
  if(length(node)>1) return(lapply(node, desc, isInternal))
  desc(node, isInternal)
}


Siblings = function (x, node, include.self = FALSE) 
{
    l = length(node)
    if(l==1){
        v <- Children(x, Ancestors(x, node, "parent"))
        if (!include.self) 
            v <- v[v != node]
        return(v)
    }
    else{    
        parents <- x$edge[, 1]
        child <- x$edge[, 2]
        pvector <- integer(max(x$edge)) # parents
        pvector[child] <- parents
        root <- as.integer(parents[!match(parents, child, 0)][1])
        res = vector("list", l)
        ch = allChildren(x)
        if(include.self) return(ch[ pvector[node] ])
        k = 1
        for(i in node){
            if(i != root){
                tmp <- ch[[ pvector[i] ]]
                res[[k]] = tmp[tmp != i]
            } 
            k=k+1    
        }     
    }
    res
}



mrca.phylo <- function(x, node=NULL){
    if(is.null(node)) return(mrca2(x))
    anc <- Ancestors(x, node, type = "all")
    res <- Reduce(intersect, anc)[1]
    res
}


# should be in ape
mrca2 <- function(phy, full=FALSE){
    Nnode <- phy$Nnode
    Ntips <- length(phy$tip.label)      
    phy <- reorder(phy)
    nodes <- unique(phy$edge[,1])
    if(!full){
        res <- Descendants(phy, nodes, "tips")
        M <- matrix(nodes[1], Ntips, Ntips)
        for(i in 2:Nnode)M[res[[i]], res[[i]]]<-nodes[i]
        diag(M) <- 1:Ntips
        dimnames(M) <-  list(phy$tip.label, phy$tip.label)
    }
    else{
        res <- Descendants(phy, nodes, "all")
        M <- matrix(nodes[1], Ntips+Nnode, Ntips+Nnode)
        for(i in 2:Nnode){
            tmp  = c(res[[i]], nodes[i])
            M[tmp, tmp]<-nodes[i]
        }    
        diag(M) <- 1:(Ntips+Nnode)
        dimnames(M) <-  list(1:(Ntips+Nnode),1:(Ntips+Nnode))
    }
    M
}


# mrca.phylo <- getMRCA


# 1090
rNNI_Old <- function(tree, moves=1, n=1){   
    k = length(na.omit(match(tree$edge[,2], tree$edge[,1])))   
    if(n==1){
        trees = tree
        for(i in 1:moves) trees = nnin(trees,sample(k,1))[[sample(2,1)]] 
        trees$tip.label <- tree$tip.label
    }  
    else{
        trees = vector("list", n)
        for(j in 1:n){
            tmp = tree 
            for(i in 1:moves) tmp = nnin(tmp, sample(k,1))[[sample(2,1)]]
            tmp$tip.label=NULL
            trees[[j]] = tmp
        }
        attr(trees, "TipLabel") <- tree$tip.label
        class(trees) <- "multiPhylo"   
    }
    trees
}


