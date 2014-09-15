#
# splits format, networx, Matrix, lento plot 
#
as.splits <- function (x, ...){
    if(inherits(x, "splits")) return(x)
    UseMethod("as.splits")
}


as.Matrix <- function (x, ...){
    if (class(x) == "Matrix") return(x)
    UseMethod("as.Matrix")
}


as.matrix.splits <- function(x, zero.print = 0L, one.print=1L, ...){
   m = length(x)
   labels = attr(x, "labels")
   n = length(labels)    
   res = matrix(zero.print, m, n)
   for(i in 1:m)res[i,x[[i]]]=one.print
   dimnames(res) = list(names(x), labels)
   res
}


as.Matrix.splits <- function(x, ...){
    labels = attr(x, "labels")
    l = length(x)
    j = unlist(x)
    i = rep(1:l, sapply(x, length))
    sparseMatrix(i,j, x = rep(1L, length(i)), dimnames = list(NULL, labels)) # included x und labels
}


print.splits <- function (x, maxp = getOption("max.print"), 
    zero.print = ".", one.print="|", ...)
{
    x.orig <- x
    cx <- as.matrix(x, zero.print = zero.print, one.print=one.print)
    print(cx, quote = FALSE, right = TRUE, max = maxp)
    invisible(x.orig)
}


"[.splits" = function(x, i){
   result = unclass(x)[i]
   if(!is.null(attr(x, "weights"))) attr(result, "weights") = attr(x, "weights")[i] 
   if(!is.null(attr(x, "data"))) attr(result, "data") = attr(x, "data")[i,, drop=FALSE] 
   attr(result, "labels") = attr(x, "labels")
   class(result) = c("splits", "prop.part")
   result
}


orderSplitLabel = function(x, order){
    label = attr(x, "labels")
    nTips = length(label)
    ord = match(label, order)
    for(i in 1:length(x))
        x[[i]] = sort(ord[x[[i]]])
    attr(x, "labels") = label[ord]
    SHORTwise(x, nTips)
}


optCycle <- function(splits, tree){
    tips = tree$tip.label
    tree = reorder(tree)
    nodes = sort(unique(tree$edge[,1]))
    
    M = as.matrix(splits)
    
    l = as.integer(nrow(M))
    m = as.integer(ncol(M))

    tmp = tree$edge[,2]
    tmp = tmp[tmp<=m]

    start <- .C("countCycle", M[, tmp], l, m, integer(1))[[4]]
#  print(start)  
    best = start
    eps = 1
    if(eps>0){
        for(i in 1:length(nodes)){
           tmptree = rotate(tree, nodes[i])
           tmp = tmptree$edge[,2]
           tmp = tmp[tmp<=m]
           tmpC <- .C("countCycle", M[, tmp], l, m, integer(1))[[4]]
#           print(tmpC)
           if(tmpC < best){
              best <- tmpC
              tree = tmptree
           }
        }
        eps = start - best
    }
#    print(best)
    tree # list(best, tree)
}


countCycles <- function(splits, tree=NULL, ord=NULL){
    M = as.matrix(splits)
    l = as.integer(nrow(M))
    m = as.integer(ncol(M))
    if(!is.null(tree)){
        tree = reorder(tree)
        nodes = sort(unique(tree$edge[,1]))
        tmp = tree$edge[,2]
        tmp = tmp[tmp<=m]
        ord = tmp
    }
    res <- .C("countCycle2", M[, ord], l, m, integer(l))[[4]]
#    which(res<=2) weakly compatible splits with ordering 
    res
}

  
c.splits <- function (..., recursive=FALSE) 
{
    x <- list(...)
    n <- length(x)
    match.names <- function(a, b) {
        if (any(!(a %in% b))) 
            stop("names do not match previous names")
    }
    if (n == 1) 
        return(x[[1]])

    labels <- attr(x[[1]], "labels")
    for (i in 2:n) {
        match.names(labels, attr(x[[i]], "labels"))
    }
    res = structure(NextMethod("c"), class=c("splits", "prop.part"))
    attr(res, "labels") = labels
    attr(res, "weight") = as.vector(sapply(x, attr, "weight"))
    res
}


# computes splits from phylo
#as.splits.phylo <- function(x, ...){
#    result = bip(x)[x$edge[,2]]
#    attr(result, "weights") = x$edge.length
#    attr(result, "labels") <- x$tip
#    class(result) = c('splits', 'prop.part')
#    result 
#}



# computes splits from phylo
as.splits.phylo <- function(x, ...){
    result = bip(x)
    if(!is.null(x$edge.length)){
        edge.weights = numeric(max(x$edge))
        edge.weights[x$edge[,2]] = x$edge.length
        attr(result, "weights") = edge.weights
    }
    attr(result, "labels") <- x$tip
    class(result) = c('splits', 'prop.part')
    result 
}


# computes splits from multiPhylo object (e.g. bootstrap, MCMC etc.)
as.splits.multiPhylo <- function(x, ...){
    if(class(x)=="multiPhylo")x = .uncompressTipLabel(x)
    lx = length(x)
    if(class(x)=="multiPhylo")class(x)='list'  # prop.part allows not yet multiPhylo
    firstTip = x[[1]]$tip[1]
    x = lapply(x, root, firstTip) # old trick  
    splits <- prop.part(x)
    class(splits)='list'
    weights = attr(splits, 'number')
    lab = attr(splits,'labels')
    attr(splits,'labels') <- attr(splits, 'number') <- NULL
    l = length(lab)
    splitTips = vector('list', l)
    for(i in 1:l) splitTips[[i]] = i
    result = c(splitTips,splits)
    attr(result, "weights") = c(rep(lx, l), weights)
    attr(result, "labels") <- lab
    class(result) = c('splits', 'prop.part')
    result  
}


as.splits.prop.part <- function(x, ...){
    if(is.null(attr(x, "number"))) attr(x, "weights") = rep(1, length(x)) 
	else attr(x, "weights") = attr(x, "number")
    class(x) = c('splits', 'prop.part')	
    x
}


as.splits.networx <- function(x, ...){
    if(!is.null(attr(x, "splits")))attr(x, "splits")
    else warning("No split object included!")    
}


# now also defined in ape
#as.prop.part <- function (x, ...){
#    if (class(x) == "prop.part") return(x)
#    UseMethod("as.prop.part")
#}


as.prop.part.splits <- function(x, ...){
    attr(x, "number") = attr(x, "weights")
    attr(x, "weights") = NULL
    class(x) = c('prop.part')	
    x
}


as.phylo.splits <- function (x, result = "phylo", ...) 
{
    result <- match.arg(result, c("phylo", "all"))
    labels = attr(x, "labels")
    nTips = length(labels)
    weights = attr(x, "weights")
    nTips = length(labels)
    x = SHORTwise(x, nTips)
    dm = as.matrix(compatible(x))
    rs = rowSums(dm)
    ind = which(rs == 0)
    if (any(rs > 0)) {
        tmp = which(rs > 0)
        candidates = tmp[order(rs[tmp])]
        for (i in candidates) {
            if (sum(dm[ind, i]) == 0) 
                ind = c(ind, i)
        }
    }
    splits = x[ind]
    weights = weights[ind]
    l = length(ind)
    res = matrix(0L, l, nTips)
    for (i in 1:l) res[i, splits[[i]]] = 1L
    dm2 = (crossprod(res * weights, 1 - res))
    dm2 = dm2 + t(dm2)
    dimnames(dm2) = list(labels, labels)
    tree <- di2multi(NJ(dm2), tol = 1e-08)
    attr(tree, "order") = NULL
    tree <- reorder(tree)    
    tree <- optCycle(x, tree)
    tree <- reorder(tree, "postorder")
    if (result == "phylo") 
        return(tree)  
#    tree = reroot(tree, Ancestors(tree, 1, "parent")) 
    spl = as.splits(tree)
    spl = SHORTwise(spl, nTips)
    spl <- spl[tree$edge[,2]]
    list(tree = tree, index = tree$edge[, 2], split = spl, rest = x[-ind])
}


# computes compatible splits
compatible <- function(obj){
    labels = attr(obj, "labels")
    if(!inherits(obj, "splits"))stop("obj needs to be of class splits")
    
    l = length(labels)
    n = length(obj)
    
    bp = matrix(0L, n, l)
    for(i in 1:n)bp[i,obj[[i]]] = 1L
    bp[bp[, 1] == 0L, ] = 1L - bp[bp[, 1] == 0L, ]
    k=1
    res = matrix(0L, n, n) 
            
    tmp1 = tcrossprod(bp) #sum(bp[i,]* bp[j,])
    tmp2 = tcrossprod(1L - bp) #sum((1L - bp[i,])*(1L - bp[j,]))
    tmp3 = tcrossprod(bp, 1L - bp) #sum(bp[i,]*(1L - bp[j,]))
    tmp4 = tcrossprod(1L - bp, bp) #sum((1L - bp[i,])*bp[j,]) 
    res[(tmp1 * tmp2 * tmp3 * tmp4)>0]=1L
    k = k+1
    
    res = res[lower.tri(res)]
    attr(res, "Size") <- n
    attr(res, "Diag") <- FALSE
    attr(res, "Upper") <- FALSE
    class(res) <- "dist"
    return(res)
    }

    
compatible2 <- function (obj1, obj2=NULL) 
{   
    if (!inherits(obj1, "splits")) 
        stop("obj needs to be of class splits")
    labels = attr(obj1, "labels")    
    l = length(labels)
    n = length(obj1)
    bp1 = as.matrix(obj1)
    bp1[bp1[, 1] == 0L, ] = 1L - bp1[bp1[, 1] == 0L, ] 
    if(!is.null(obj2)){
        m = length(obj2) 
        bp2 = as.matrix(obj2)
        labels2 = attr(obj2, "labels")
        bp2 = bp2[, match(labels2, labels), drop=FALSE]
        bp2[bp2[, 1] == 0L, ] = 1L - bp2[bp2[, 1] == 0L, ]
    }
    else bp2 = bp1

    if(is.null(obj2)) res = matrix(0L, n, n)
    else res = matrix(0L, n, m)

    tmp1 = tcrossprod(bp1, bp2)
    tmp2 = tcrossprod(1L - bp1, 1L - bp2)
    tmp3 = tcrossprod(bp1, 1L - bp2)
    tmp4 = tcrossprod(1L - bp1, bp2)
    res[(tmp1 * tmp2 * tmp3 * tmp4) > 0] = 1L
    if(is.null(obj2)){
        res = res[lower.tri(res)]
        attr(res, "Size") <- n
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        class(res) <- "dist"
    }
    return(res)
}


compatible3 <- function(x, y=NULL) 
{
    if (!inherits(x, "splits")) 
        stop("x needs to be of class splits")
    if(is.null(y)) y <- x
        
    if (!inherits(y, "splits")) 
        stop("y needs to be of class splits")
    xlabels = attr(x, "labels")
    ylabels = attr(y, "labels")
    if(identical(xlabels, ylabels)) labels = xlabels 
    else labels = intersect(xlabels, ylabels)
#    if(length(labels) maybe warning
    nx = length(x)
    ny = length(y)   
    bp1 = as.matrix(x)[,labels, drop=FALSE]
    bp2 = as.matrix(y)[,labels, drop=FALSE]
    rs1 = rowSums(bp1)
    rs2 = rowSums(bp2)
    res = matrix(0L, nx, ny)
    tmp1 = tcrossprod(bp1, bp2)
    res = matrix(0L, nx, ny)
    for(i in 1:nx){
        for(j in 1:ny){            
            if(tmp1[i, j]==rs1[i]) res[i,j] = 1
            if(tmp1[i, j]==rs2[j]) res[i,j] = 2
            if(tmp1[i, j]==rs1[i] & tmp1[i, j]==rs2[j])res[i,j] = 3
        }
    }      
    if(is.null(y)){
        res = res[lower.tri(res)]
        attr(res, "Size") <- length(x)
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        class(res) <- "dist"
    }
    return(res)
}
    

#
# networx
#
addEdge <- function(network, dm, nsplit, w0, nTips, x){ # desc, split, index, dm, l){
  edge = network$edge
  parent = edge[,1]
  child = edge[,2]
  
  index = network$split
  uindex = unique(index)
  #  blub = match(index, uindex)
  
  ind = which(dm[index, nsplit]==1)
  if(length(ind)==0) return(network)
  add = TRUE
  
  while(add){
    tmp = ind
    for(i in ind){
      tmp2 = which(dm[index[i], index] == 1) 
      tmp = union(tmp, tmp2)      
    }
    if(identical(ind, tmp)){add=FALSE}
    ind=tmp 
  }    
  
  oldNodes = unique(as.vector(edge[ind,]))
  newNodes = (max(parent)+1L) : (max(parent)+length(oldNodes))
  
  ind2 = index[-ind]
  edge2 = edge[-ind,, drop=FALSE] 
  
  for(i in 1:length(oldNodes)){
    ind3 = which(edge2[,1] == oldNodes[i])
    #        ind3 = which(edge2 == oldNodes[i],arr=TRUE)[,1]
    for(j in ind3) {
      if(any( x[[ ind2[j] ]] %in% x[[nsplit]])){
        edge2[j,edge2[j,]==oldNodes[i]] = newNodes[i]
      }
    } 
  } 
  edge[-ind,] = edge2
  
  #alle Splits verdoppeln
  dSpl = edge[ind,]
  for(i in 1:length(oldNodes)) dSpl[dSpl==oldNodes[i]] = newNodes[i]
  edge = rbind(edge, dSpl, deparse.level = 0) # experimental: no labels
  network$edge.length = c(network$edge.length, network$edge.length[ind])   
  index = c(index, index[ind])
  
  #neu zu alt verbinden   
  edge = rbind(edge, cbind(oldNodes, newNodes), deparse.level = 0) #  rbind(edge, cbind(oldNodes, newNodes), deparse.level = 0)# experimental: no labels
  network$edge.length = c(network$edge.length, rep(w0, length(oldNodes))) # w0 = attr(split,"weights")
  #    index = c(index, rep(max(index)+1, length(oldNodes)) )
  index = c(index, rep(nsplit, length(oldNodes)) )
  network$edge = edge
  network$Nnode = length(unique(edge[,1]))
  network$split = index
  network   
}



as.networx <- function (x, ...) 
{
    if (inherits(x, "networx")) 
        return(x)
    UseMethod("as.networx")
}


as.networx.splits <- function(x, include.splits=TRUE, ...){
  label <- attr(x, "label")
  weight <- attr(x, "weights")
  nTips <- length(label)
  if(!is.null(attr(x, "cycle"))){
#    tmp <- stree(length(label), tip.label=label)
#    tmp$edge[,2] <- as.integer(attr(x, "cycle"))  
      tmp <- as.phylo(x)
  }
  else tmp <- as.phylo(x) 
  
  tmp = reorder(tmp)

  dm <- as.matrix(compatible2(x))
  
  x <- SHORTwise(x, nTips)
  sp <- as.splits(tmp)[tmp$edge[,2]]
  sp <- SHORTwise(sp, nTips)
# which splits are in circular ordering  
  circSplits = which(countCycles(x, tmp)==2)  

  cycord = tmp$edge[,2]
  cycord = cycord[cycord <= nTips] 

  ll <- sapply(x, length)
  ind <- match(sp, x)
  ind2 = union(ind, which(ll==0)) # which(duplicated(x))

  tmp$split = ind
  ord <- order(colSums(dm))
  ord <- setdiff(ord, ind2)
  #  browser()
  if(length(ord)>0){    
      for(i in 1:length(ord)){ 
          tmp = addEdge(tmp, dm, ord[i], weight[ord[i]], nTips, x)
          class(tmp) = c("networx", "phylo")
      } 
  }
  if(include.splits)attr(tmp, "splits") = x 
  class(tmp) = c("networx", "phylo")
  tmp
}

#consensusNet <- function(obj, prob=.3, ...){
#    l = length(obj)
#    spl = as.splits(obj)
#    w = attr(spl, "weight")
#    ind = (w/l) > prob 
#    spl = spl[ind] 
#    as.networx(spl)
#}

consensusNet <- function (obj, prob = 0.3, ...) 
{
    l = length(obj)
    spl = as.splits(obj)
    w = attr(spl, "weight")
    ind = (w/l) > prob
    spl = spl[ind]
    edge.labels = as.character(round((w/l)[ind]*100))
    attr(spl, "confidences") = round((w/l)[ind]*100)
    edge.labels[1:length(attr(spl,"labels"))]=""
    spl = as.networx(spl)
    spl$edge.labels = as.character(spl$edge.length / l * 100)
    spl$edge.labels[spl$edge[,2]<=length(spl$tip.label)] = ""
    spl
}


# rename X to obj
addConfidences <- function(X, phy){
    tiplabel <- attr(X, "label")
    ind <- match(tiplabel, phy$tip.label)
    if (any(is.na(ind)) | length(tiplabel) != length(phy$tip.label)) 
        stop("trees have different labels")
    phy$tip.label <- phy$tip.label[ind]
    ind2 <- match(1:length(ind), phy$edge[, 2])
    phy$edge[ind2, 2] <- order(ind)
    
    spl <- as.splits(phy)
    
    nTips <- length(tiplabel)
    spl <- phangorn:::SHORTwise(spl, nTips)
    ind <- match(phangorn:::SHORTwise(X, nTips), spl)
    pos <-  which(ind > nTips)
    confidences <- numeric(length(X))
    confidences[pos] <- phy$node.label[ind[pos] - nTips]
    attr(X, "confidences") <- confidences
    X  
}



reorder.networx <- function (x, order = "cladewise", ...) 
{
    order <- match.arg(order, c("cladewise"))
    if (!is.null(attr(x, "order"))) 
        if (attr(x, "order") == order) 
            return(x)
    nb.node <- x$Nnode
    if (nb.node == 1) 
        return(x)
    nb.tip <- length(x$tip.label)
    nb.edge <- dim(x$edge)[1]
    #neworder <- if (order == "cladewise") 
    neworder <- .C("neworder_cladewise", as.integer(nb.tip), as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
               as.integer(nb.edge), integer(nb.edge), PACKAGE = "phangorn")[[5]]

    x$edge <- x$edge[neworder, ]
    if (!is.null(x$edge.length)) 
        x$edge.length <- x$edge.length[neworder]
    if (!is.null(x$edge.labels)) 
        x$edge.labels <- x$edge.labels[neworder]  
    if (!is.null(x$split))x$split <- x$split[neworder]
    attr(x, "order") <- order
    x
}


coords <- function(obj, dim="3D"){
    if(is.null(attr(obj,"order")) || attr(obj, "order")=="pruningwise") #obj <- reorder(obj)
        obj = reorder.networx(obj)

    l = length(obj$edge.length)
    ind1 = which(!duplicated(obj$split))

    n = max(obj$edge)
    adj = Matrix::spMatrix(n, n, i = obj$edge[,2], j = obj$edge[,1], x = rep(1, length(obj$edge.length)))
    g = graph.adjacency(adj, "undirected")

    g2 <- graph(t(obj$edge), directed=FALSE)
    g2 <- set.edge.attribute(g, "weight", value=obj$edge.length)
    if(dim=="3D"){
        coord <- layout.kamada.kawai(g, dim=3)
        k = matrix(0, max(obj$split), 3)
        for(i in ind1){
            tmp = coord[obj$edge[i, 2],] - coord[obj$edge[i, 1],]
            k[obj$split[i], ] = kart2kugel(tmp[1], tmp[2], tmp[3])
        }
        k[obj$split[ind1],1] = obj$edge.length[ind1] 

        res = matrix(0, vcount(g), 3)
        for(i in 1:l){# unique(obj$split)
            j = obj$edge[i,1]
            m = obj$edge[i,2]
            p = obj$split[i]
            res[m,] = res[j,] + kugel2kart(k[p,1], k[p,2], k[p,3])     
        }            
    }
    else{
        coord <- layout.kamada.kawai(g, dim=2)
        k = matrix(0, max(obj$split), 2)
        for(i in ind1){
            tmp = coord[obj$edge[i, 2],] - coord[obj$edge[i, 1],]
            k[obj$split[i], ] = kart2kreis(tmp[1], tmp[2])
        }
        k[obj$split[ind1],1] = obj$edge.length[ind1] 
        res = matrix(0, vcount(g), 2)
        for(i in 1:l){# unique(obj$split)
            j = obj$edge[i,1]
            m = obj$edge[i,2]
            p = obj$split[i]
            res[m,] = res[j,] + kreis2kart(k[p,1], k[p,2])     
        }
    }  
    res  
}


kart2kugel <- function(x,y,z){
    r = sqrt(x*x+y*y+z*z)
    alpha = atan(sqrt(x*x+y*y) / z)
    if(z<0) alpha = alpha+pi
    beta = atan(y/x)
    if(x<0) beta = beta+pi 
    c(r,alpha,beta)
}

	
kart2kreis <- function(x,y){
    r = sqrt(x*x+y*y)
    alpha = atan(y/x) 
    if(x<0) alpha = alpha+pi
    c(r,alpha)
}	
	

kreis2kart <- function(r,alpha){
	c(r*cos(alpha), r*sin(alpha))
}


kugel2kart <- function(r,alpha,beta){
    x = r * sin(alpha) * cos(beta) 
    y = r * sin(alpha) * sin(beta) 
    z = r * cos(alpha)
    c(x,y,z)
}


edgeLabels <- function(xx,yy,zz=NULL, edge){
        XX <- (xx[edge[, 1]] + xx[edge[, 2]])/2
        YY <- (yy[edge[, 1]] + yy[edge[, 2]])/2
        if(!is.null(zz)){
	        ZZ <- (zz[edge[, 1]] + zz[edge[, 2]])/2
	        return(cbind(XX, YY, ZZ))
        }  
        cbind(XX, YY)  
}


plot.networx = function(x, type="3D", show.tip.label=TRUE, show.edge.label=FALSE, show.nodes=FALSE, 
    tip.color = "blue", edge.color="grey", edge.width = 3, 
    font = 3, cex = 1, ...){
    type = match.arg(type, c("3D", "2D")) 
    n = max(x$edge)
    tip = rep(NA, n)
    tips = x$tip.label
    tip[1:length(tips)] = tips
    
    x = reorder(x)
    
    adj = spMatrix(n, n, i = x$edge[,2], j = x$edge[,1], x = rep(1, length(x$edge.length)))
    g = graph.adjacency(adj, "undirected")
    
#    plot.success <- FALSE
#    if (!plot.success & type=="3D") {
    if (type=="3D") {
  #   require(rgl) &         
  #     if (!require(rgl)) {
  #          warning("package 'rgl' not found, can only plot in 2D")
  #      } else {       
             coord <- coords(x, dim="3D")
             plotRGL(coord, x, show.tip.label=show.tip.label, show.edge.label=show.edge.label, show.nodes=show.nodes, tip.color = tip.color,
             edge.color=edge.color, edge.width = edge.width, font = font, cex = cex)
             plot.success <- TRUE
#        } 
    }
    #if (!plot.success){
   else{
	    coord <- coords(x, dim="2D")
	    plot2D(coord, x, show.tip.label=show.tip.label, show.edge.label=show.edge.label, tip.color = tip.color, edge.color=edge.color, 
	    edge.width = edge.width, font = font, cex = cex, add=FALSE)
	    }    
}

    
plotRGL <- function(coords, net, show.tip.label=TRUE, show.edge.label=FALSE, show.nodes=FALSE, tip.color = "blue", edge.color="grey", 
    edge.width = 3, font = 3, cex = par("cex"), ...){
    edge = net$edge
  
    x = coords[,1]
    y = coords[,2]
    z = coords[,3]
     
    nTips = length(net$tip.label)
    
    segments3d(x[t(edge)],y[t(edge)],z[t(edge)], col=edge.color, lwd=edge.width) 
    radius=0
    if(show.nodes){
        radius = sqrt((max(x)-min(x))^2 + (max(y)-min(y))^2 + (max(z)-min(z))^2) / 200    
        spheres3d(x[1:nTips], y[1:nTips],z[1:nTips], radius=2*radius, color="cyan")
        spheres3d(x[-c(1:nTips)], y[-c(1:nTips)],z[-c(1:nTips)], radius=radius, color="magenta")
    }
    if(show.tip.label){
      rgl.texts(x[1:nTips]+2.05*radius,y[1:nTips],z[1:nTips],net$tip.label, color=tip.color, cex=cex, font=font)
    }
    if(show.edge.label){
	    ec = edgeLabels(x, y, z, edge)
        if(!is.null(net$edge.labels)) edge.labels = net$edge.labels
        else edge.labels = net$split    
	    rgl.texts(ec[,1], ec[,2], ec[,3], edge.labels, color=tip.color, cex=cex, font=font)     
    } 
}


plot2D <- function(coords, net, show.tip.label=TRUE, show.edge.label=FALSE, tip.color = "blue", edge.color="grey", edge.width = 3, 
    font = 3, cex = par("cex"), add=FALSE, ...){
   edge = net$edge
   label = net$tip.label
   xx = coords[,1]
   yy = coords[,2]
   nTips = length(label)

   cex=1
   
   xlim <- range(xx)
   ylim <- range(yy)
     
   if(show.tip.label){
       offset <- max(nchar(label)) * 0.018 * cex * diff(xlim)
       xlim = c(xlim[1]-offset, xlim[2]+offset)
       ylim = c(ylim[1]-0.03 * cex * diff(ylim), ylim[2]+0.03 * cex * diff(ylim))
   }
   if(!add){ 
       plot.new() 
       plot.window(xlim, ylim, asp=1)
   }
   cladogram.plot(edge, xx, yy, edge.color, edge.width, 1)
   if(show.tip.label){
        ind=match(1:nTips, edge[,2])
        pos = rep(4, nTips)
        XX <- xx[edge[ind, 1]] - xx[edge[ind, 2]]
        pos[XX>0] = 2
        YY <- yy[edge[ind, 1]] - yy[edge[ind, 2]]
        pos2 <- rep(3, nTips)
        pos2[YY>0] = 1
        pos[abs(YY)>abs(XX)] <- pos2[abs(YY)>abs(XX)] 	
        text(xx[1:nTips], yy[1:nTips], labels=label, pos=pos, col=tip.color, cex=cex, font=font)
    }
    if(show.edge.label){
	    ec = edgeLabels(xx,yy, edge=edge)
	    if(is.null(net$edge.labels))net$edge.labels = net$split
	    text(ec[,1], ec[,2], labels=net$edge.labels, col=tip.color, cex=cex, font=font)     
	    } 
}   
   
    
lento <- function (obj, xlim = NULL, ylim = NULL, main = "Lento plot", 
    sub = NULL, xlab = NULL, ylab = NULL, bipart=TRUE, trivial=FALSE, ...) 
{
    if (class(obj) == "phylo") 
        obj = as.splits(obj)
    if (class(obj) == "multiPhylo") 
        obj = as.splits(obj)    
    labels = attr(obj, "labels") 
    l = length(labels)
    if(!trivial){
        triv = sapply(obj, length)
        ind = logical(length(obj)) 
        ind[(triv >1) & (triv < (l-1))] = TRUE
        obj = obj[ind]
        }
    CM = compatible(obj)
    support = attr(obj, "weights")
    if (is.null(support)) 
        support = rep(1, length(obj))
    conflict = -as.matrix(CM) %*% support
    n = length(support)
    if (is.null(ylim)) {
        eps = (max(support) - min(conflict)) * 0.05
        ylim = c(min(conflict) - eps, max(support) + eps)
    }
    if (is.null(xlim)) {
        xlim = c(0, n + 1)
    }

    ord = order(support, decreasing = TRUE)
    support = support[ord]
    conflict = conflict[ord]
    plot.new()
    plot.window(xlim, ylim)
    title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
    segments(0:(n - 1), support, y1 = conflict, ...)
    segments(1:n, support, y1 = conflict, ...)
    segments(0:(n - 1), support, x1 = 1:n, ...)
    segments(0:(n - 1), conflict, x1 = 1:n, ...)
    abline(h = 0)
    axis(2, ...)
    aty = diff(ylim)/(l+1)
    at = min(ylim) + (1:l) * aty
    if(bipart){
        Y = rep(at, n)
        X = rep((1:n)-.5, each=l)
        Circles = matrix(1, l, n)
        for(i in 1:n) Circles[obj[[ord[i]]],i] = 19   
#    axis(4, labels=labels, at=at)
        text(x=n+.1,y=at, labels, pos=4, ...) 
        points(X,Y,pch = as.numeric(Circles), col = rgb(0,0,0,.5), ...)
        }
    invisible(cbind(support, conflict))
    }

    
write.splits = function (x, file = "", zero.print = ".", one.print = "|", print.labels = TRUE, ...) 
{
    labels = attr(x, "labels")
    x.orig <- x
    cx <- as.matrix(x, zero.print = zero.print, one.print = one.print)
    w = FALSE
    if (!is.null(attr(x, "names"))) {
        nam = TRUE
        vnames = format(attr(x, "names"))
    }
    nam = FALSE
    if (!is.null(attr(x, "weight"))) {
        w = TRUE
        weight = format(attr(x, "weight"))
    }
    d = FALSE
    if (!is.null(attr(x, "data"))) {
        d = TRUE
        data = attr(x, "data")
    }
    if(print.labels){for(i in 1:length(labels)) cat(labels[i], "\n", file = file, append = TRUE)}
    if (w) 
        cat("weight", "\t", file = file, append = TRUE)
    if (d) 
        cat(paste(colnames(data), "\t"), file = file, append = TRUE)
    cat("\n", file = file, append = TRUE) #"Matrix", 
    for (i in 1:length(x)) {
        if (nam) 
            cat(vnames[i], "\t", file = file, append = TRUE)
        if (d) 
            cat(paste(data[i, ], "\t"), file = file, append = TRUE)
        if (w) 
            cat(weight[i], "\t", file = file)
        cat("\n", paste(cx[i, ], collapse = ""),"\n",  file = file, append = TRUE)
    }
}
 

# tries to be compatible with splitstree and spectronet
# FORMAT labels=no weights=yes confidences=no intervals=no;
write.nexus.splits <- function (obj, file = "", weights=NULL) 
{
    if(is.null(weights))weight <- attr(obj, "weights")
    taxa.labels <- attr(obj, "labels")
    ntaxa = length(taxa.labels)
    nsplits = length(obj)
    
    if (is.null(weight)) 
        weight = numeric(nsplits) + 100
    cat("#NEXUS\n\n", file = file)
    cat("[Splits block for Spectronet or Splitstree]\n", file = file, append = TRUE)
    cat("[generated by phangorn:\n", file = file, append = TRUE)
    cat(format(citation("phangorn"), "text"), "]\n\n",
       file = file, append = TRUE)
    cat(paste("BEGIN TAXA;\n\tDIMENSIONS NTAX=", ntaxa, ";\n", 
        sep = ""), file = file, append = TRUE)
    cat("\tTAXLABELS", paste(taxa.labels, sep = " "), ";\nEND;\n\n", 
        file = file, append = TRUE)
    cat(paste("BEGIN ST_SPLITS;\n\tDIMENSIONS NSPLITS=", nsplits, 
        ";\n", sep = ""), file = file, append = TRUE)
# labels=YES/NO, WEIGHTS, CONFIDENCES, INTERVALS     
    format = "\tFORMAT labels=yes weights=yes"
    fcon = fint = flab = FALSE
    if(!is.null(attr(obj, "confidences"))){ 
        format = paste(format, "confidences=yes")
        fcon=TRUE
    }
    else format = paste(format, "confidences=no") 
    if(!is.null(attr(obj, "intervals"))){ 
        format = paste(format, "intervals=yes")
        fint=TRUE
    }
    else format = paste(format, "intervals=no") 
    if(!is.null(attr(obj, "splitlabels"))) flab=TRUE
    
    format = paste(format, ";\n",  sep = "")
#    cat("\tFORMAT LABELS WEIGHTS;\n", file = file, append = TRUE)
    cat(format, file = file, append = TRUE)
    cat("\tMATRIX\n", file = file, append = TRUE)    

    for (i in 1:nsplits){
        slab <- ifelse(flab, attr(obj, "splitlabels")[i], i)
        scon <- ifelse(fcon, paste(attr(obj, "confidences")[i], "\t"), "")
        sint <- ifelse(fint, paste(attr(obj, "intervals")[i], "\t"), "")
        cat("\t\t", slab, "\t", weight[i], "\t", scon, sint, paste(obj[[i]]), 
            ",\n", file = file, append = TRUE, sep = "")  
    }
    cat("\t;\nEND;\n", file = file, append = TRUE)
}


read.nexus.splits <- function(file)
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    semico <- grep(";", X)
    X=gsub("\\[(.*?)\\]", "", X) # get rid of comments
    i1 <- grep("TAXLABELS", X, ignore.case = TRUE)    
    taxlab <- TRUE 
    if (taxlab) {
        end <- semico[semico > i1][1]
        x <- X[(i1 + 1):end] # assumes there's a 'new line' after "TRANSLATE"
        ## x <- gsub("TRANSLATE", "", x, ignore.case = TRUE)
        x <- unlist(strsplit(x, "[,; \t]"))   
        x <- x[nzchar(x)]
        x <- gsub("['\"]", "", x)
        xntaxa <- length(x)
    }
    sp <- grep("SPLITS;", X, ignore.case = TRUE)
    dims <- grep("DIMENSION", X, ignore.case = TRUE)
    cyc <- grep("CYCLE", X, ignore.case = TRUE)
    matr <- grep("MATRIX", X, ignore.case = TRUE)
    format <- grep("FORMAT", X, ignore.case = TRUE)
    start <- matr[matr>sp][1] + 1
    end <- semico[semico>start][1] -1 
    res <- vector("list", end - start + 1)
    weights = numeric(end - start + 1)
    j=1
    
    flab = fwei = fcon = fint = FALSE
    
    if(length(format)>0){
        tmp = X[format]    
        tmp = gsub("\\;", "", tmp)
        tmp = gsub("\\s+", "", tmp)
        flab = grepl("labels=yes", tmp, ignore.case = TRUE) 
        fwei = grepl("weights=yes", tmp, ignore.case = TRUE) 
        fcon = grepl("confidences=yes", tmp, ignore.case = TRUE) 
        fint = grepl("intervals=yes", tmp, ignore.case = TRUE) 
        # = as.numeric(na.omit(as.numeric(strsplit(tmp, " ")[[1]])))        
        ind = cumsum(c(flab, fwei, fcon, fint))
        mformat = sum(c(flab, fwei, fcon, fint))
     }
    
    if(fint)intervals = numeric(end - start + 1)
    if(fcon)confidences = numeric(end - start + 1)
    if(flab)labels = vector("character", end - start + 1)
    
    for(i in start:end){
        tmp = X[i]
        tmp = sub("\\s+", "", tmp) 
        tmp = strsplit(tmp, "\t")[[1]]
        if(length(tmp)!=(mformat+1)) warning("blub")
        if(flab)labels[j] = as.numeric(tmp[ind[1]])        
        if(fwei)weights[j] = as.numeric(tmp[ind[2]])
        if(fcon)confidences[j] = as.numeric(tmp[ind[3]])
        if(fint)intervals[j] = as.numeric(tmp[ind[4]])
        tmp = tmp[length(tmp)]
        tmp = gsub("\\,", "", tmp)
        res[[j]] = as.numeric(na.omit(as.numeric(strsplit(tmp, " ")[[1]])))
        j=j+1
    }
# read in cycle   
    if(length(cyc)>0){
        tmp = X[cyc]    
        tmp = gsub("\\;", "", tmp)
        tmp = gsub("CYCLE", "", tmp, ignore.case = TRUE)
        tmp = sub("\\s+", "", tmp)
        cyc = as.numeric(na.omit(as.numeric(strsplit(tmp, " ")[[1]])))
    }

    attr(res, "labels") = x
    attr(res, "weights") = weights
    if(fint)attr(res, "intervals") = intervals
    if(fcon)attr(res, "confidences") = confidences
    if(flab)attr(res, "splitlabels") = labels
    attr(res, "cycle") = cyc 
    class(res) = c("splits", "prop.part")
    res
}


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


