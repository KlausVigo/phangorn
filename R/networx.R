#' Phylogenetic Network
#' 
#' \code{splitsNetwork} estimates weights for a splits graph from a distance
#' matrix.
#' 
#' \code{splitsNetwork} fits non-negative least-squares phylogenetic networks
#' using L1 (LASSO), L2(ridge regression) constraints.  The function minimizes
#' the penalized least squares 
#' \deqn{\beta = min \sum(dm - X\beta)^2 + \lambda \|\beta \|^2_2 }{ beta = sum(dm - X*beta)^2 + lambda |beta|^2_2 } 
#' with respect to \deqn{\|\beta \|_1 <= \gamma, \beta >= 0}{ |beta|_1 = gamma, beta >= 0} 
#' where \eqn{X} is a design matrix constructed with \code{designSplits}.
#' External edges are fitted without L1 or L2 constraints.
#' 
#' @param dm A distance matrix.
#' @param splits a splits object, containing all splits to consider, otherwise
#' all possible splits are used
#' @param gamma penalty value for the L1 constraint.
#' @param lambda penalty value for the L2 constraint.
#' @param weight a vector of weights.
#' @return \code{splitsNetwork} returns a splits object with a matrix added.
#' The first column contains the indices of the splits, the second column an
#' unconstrained fit without penalty terms and the third column the constrained
#' fit.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[phangorn]{distanceHadamard}},
#' \code{\link[phangorn]{designTree}} \code{\link[phangorn]{consensusNet}},
#' \code{\link[phangorn]{plot.networx}}
#' @references Efron, Hastie, Johnstone and Tibshirani (2004) Least Angle
#' Regression (with discussion) \emph{Annals of Statistics} \bold{32(2)}, 407--499
#' 
#' K. P. Schliep (2009). Some Applications of statistical phylogenetics (PhD
#' Thesis)
#' @keywords cluster
#' @examples
#' 
#' data(yeast)
#' dm = dist.ml(yeast)
#' fit = splitsNetwork(dm)
#' net = as.networx(fit)
#' plot(net, "2D")
#' write.nexus.splits(fit)
#' 
#' @export splitsNetwork
splitsNetwork <- function(dm, splits=NULL, gamma=.1, lambda=1e-6, weight=NULL){
  dm = as.matrix(dm)
  k = dim(dm)[1]
  
  if(!is.null(splits)){
#    tmp = which(sapply(splits, length)==k)
    tmp = which(lengths(splits)==k)
    splits = splits[-tmp]
    lab = attr(splits, "labels")
    dm = dm[lab, lab]
  }
  
  if(is.null(splits)){
    X2 = designAll(k, TRUE)
    X=X2[[1]]
  }
  else X = as.matrix(splits2design(splits))
  
  y = dm[lower.tri(dm)]
  if(is.null(splits))ind = c(2^(0:(k-2)),2^(k-1)-1)
  else ind = which(lengths(splits)==1)
#  else ind = which(sapply(splits, length)==1)  
  #   y2 = lm(y~X[,ind]-1)$res
  n = dim(X)[2]
  
  ridge <- lambda * diag(n) 
  ridge[ind,ind] <- 0
  if(!is.null(weight)) Dmat <- crossprod(X * sqrt(weight)) + ridge
  else Dmat <- crossprod(X) + ridge
  if(!is.null(weight)) dvec <- crossprod(X * sqrt(weight),y * sqrt(weight))
  else dvec <- crossprod(X, y)
  
  #    Dmat <- as.matrix(Dmat)
  #    dvec <- as.vector(dvec) 
  
  ind1       <- rep(1,n)
  ind1[ind]  <- 0 
  
  Amat       <- cbind(ind1,diag(n)) 
  bvec       <- c(gamma, rep(0,n))
  
# needs quadprog::solve.QP.compact  
  solution <- quadprog::solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=1)$sol   
  
  ind2 <- which(solution>1e-8)
  n2 <- length(ind2)
  
  ind3 = which(duplicated(c(ind2, ind), fromLast = TRUE)[1:n2])
  ridge2 <- lambda * diag(n2) 
  ridge2[ind3,ind3] <- 0
  
  if(!is.null(weight)) Dmat <- crossprod(X[, ind2] * sqrt(weight)) + ridge2
  else Dmat <- crossprod(X[, ind2]) + ridge2
  if(!is.null(weight)) dvec <- crossprod(X[, ind2] * sqrt(weight),y * sqrt(weight))
  else dvec <- crossprod(X[, ind2], y)
  
  Amat2 <- diag(n2)
  bvec2 <- rep(0, n2)
 # needs quadprog::solve.QP.compact 
 # bvec2 not used
  solution2  <- quadprog::solve.QP(Dmat, dvec, Amat2)$sol
  
  RSS1 = sum((y-X[,ind2]%*%solution[ind2])^2)
  RSS2 = sum((y-X[,ind2]%*%solution2)^2)
  
  if(is.null(splits)){
    splits = vector("list", length(ind2))
    for(i in 1:length(ind2))splits[[i]] = which(X2[[2]][ind2[i],]==1)
  } 
  else splits = splits[ind2]
  attr(splits, "weights") = solution[ind2]
  attr(splits, "unrestricted") = solution2
  attr(splits, "stats") = c(df=n2, RSS_p = RSS1, RSS_u=RSS2)
  attr(splits,"labels") =dimnames(dm)[[1]]
  class(splits)='splits'
  return(splits)           
}


#' @rdname as.splits
#' @export
allSplits = function(k, labels=NULL){
  result <- lapply(1:(2^(k-1)-1),dec2Bin)
  if(is.null(labels)) labels=(as.character(1:k))
  attr(result, 'labels') =labels
  class(result)='splits'
  result
}   


#' @rdname as.splits
#' @export
allCircularSplits <- function(k, labels=NULL){
    k = as.integer(k)
    l = (k-1L) %/% 2L
    res <- vector("list", k*(k-1L)/2)
    
    res[1:k] = 1L:k
    ind = k
    if(k>3){
        fun = function(x,y){
            tmp = (1L:y)+x
            tmp %% (k+1L) + tmp %/% (k+1L)
        }
        if(k>4L){
            for(i in 2:l){
                res[(ind+1):(ind+k)] <- lapply(0L:(k-1L), fun, i)
                ind <- ind+k
            }
        }
        if((k%%2L)==0){
            m <- k%/%2
            res[(ind+1):(ind+m)] <- lapply(0L:(m-1L), fun, m)
        }
        
    }   
    if(is.null(labels)) labels <- as.character(1:k)
    attr(res, 'labels') <- labels
    attr(res, "cycle") <- 1:k
    class(res)="splits"
    res   
}


getIndex = function(left, right, n){
  if(n<max(left) | n<max(right)) stop("Error")  
  left = as.integer(left)
  right = as.integer(right)
  ll = length(left)
  lr = length(right)
  .C("giveIndex", left, right, ll, lr, as.integer(n), integer(ll*lr))[[6]]+1
}


splits2design <- function(obj, weight=NULL){
  labels= attr(obj,'labels')
  m = length(labels)
  n=length(obj)
  l = 1:m 
  sl = lengths(obj)
#  sl = sapply(obj, length)
  p0 = sl * (m-sl)
  p = c(0,cumsum(p0))
  i = numeric(max(p))
  for(k in 1:n){
    sp = obj[[k]]
    if(p0[k]!=0) i[(p[k]+1):p[k+1]] = getIndex(sp, l[-sp], m) 
  }
  dims=c(m*(m-1)/2,n)
  sparseMatrix(i=i, p=p, x=1.0, dims=dims) 
}


hC <- function(g, set){
    intersec = NULL
    allEdges = NULL
    fromTo <- set
    l = length(set)
    sptmp = shortest_paths(g, fromTo[l], fromTo[1], output=c("epath"))$epath[[1]]
    sptmp = as.vector(sptmp)
    allEdges = sptmp
    for(i in 2:length(set)){
        sptmp = shortest_paths(g, fromTo[i-1], fromTo[i], output=c("epath"))$epath[[1]]
        sptmp = as.vector(sptmp)
        intersec = c(intersec, intersect(allEdges, sptmp) )
        allEdges = c(allEdges, sptmp)
    }
    #    allEdges = unique(allEdges)
    list(allEdges, unique(allEdges), intersec)
}   


addEdge <- function(network, desc, spl){   
    edge <- network$edge
    parent <- edge[,1]
    child <- edge[,2]
    nTips <- length(network$tip.label)

    desc2 <- SHORTwise(desc, nTips)    
    split <- desc2[spl]
        
    index <- network$splitIndex
    ind <- which(compatible2(split, desc2[index]) == 1)
    if(is.null(ind) | (length(ind)==0)) return(network)
    add <- TRUE
  
    X <- as.matrix(desc2)
    rsX <- rowSums(X)
    z <- X %*% X[spl,]
    v <- which((rsX == z)[index] == TRUE) 

    
# intersection of shortest pathes of both partitions
# best with similar to circNetwork with shortest_paths 
    
    while(add){
        tmp = ind
        for(i in ind){          
            tmp2 = which(compatible2(desc2[index][i], desc2[index]) == 1)
            tmp = union(tmp, tmp2)
        }
        if(identical(ind, tmp)){
            ind=tmp           
            add=FALSE
        }
        ind=tmp
    }    
   

    g = graph(t(network$edge[ind,]), directed=FALSE)
    dec = decompose(g, min.vertices = 2)

    #    fromTo <- sort(match(split[[1]], attr(desc, "cycle")))
    #    sptmp = shortest_paths(g, fromTo[i-1], fromTo[i], 
    #                           output=c("epath"))$epath[[1]]
    #    sp2 = c(sp2, sptmp[-c(1, length(sptmp))])
    #    sp0 = c(sp0, sptmp)
    
    oldNodes = unique(as.vector(edge[ind,]))
    mNodes = max(network$edge)
    newNodes = (mNodes+1L) : (mNodes+length(oldNodes))

# duplicated splits
    dSpl = edge[ind,]
    edge2 = edge[v,] 
    for(i in 1:length(oldNodes)){
        edge2[edge2 == oldNodes[i]] = newNodes[i]
    } 
    edge[v,] = edge2    

  #alle Splits verdoppeln
    for(i in 1:length(oldNodes)) dSpl[dSpl==oldNodes[i]] = newNodes[i]
    edge = rbind(edge, dSpl, deparse.level = 0) # experimental: no labels
    index = c(index, index[ind])
  #neu zu alt verbinden   
    edge = rbind(edge, cbind(oldNodes, newNodes), deparse.level = 0) 
    index = c(index, rep(spl, length(oldNodes)) )
    network$edge = edge
    network$Nnode = max(edge) - nTips
    network$splitIndex = index
    network   
}

## as.splits.phylo
circNetwork <- function(x, ord=NULL){
    if(is.null(ord))ord = attr(x, "cycle")
    
    weight <- attr(x, "weights")
    if(is.null(weight)) weight = rep(1, length(x))
    nTips = length(ord)
    tmp = which(ord == 1)
    if(tmp!=1) ord = c(ord[tmp:nTips], ord[1:(tmp-1)])
    res = stree(nTips, tip.label = attr(x, "labels"))
    res$edge[, 2] = ord
    res$edge.length=NULL
    x <- SHORTwise(x, nTips)    
    spRes <- as.splits(res)[res$edge[,2]]
    index = match(spRes, x)
    
    if(any(is.na(index))){
        l.na = sum(is.na(index))
        x <- c(x, spRes[is.na(index)])    
        weight = c(weight, rep(0, l.na))
        index = match(spRes, x)
    }
    
    l = lengths(oneWise(x, nTips))
    l2 = lengths(x)
#    l = sapply(oneWise(x, nTips), length)
#    l2 = sapply(x, length)

    #    dm <- as.matrix(compatible2(x))
    
    tmp <- countCycles(x, ord=ord)
    ind = which(tmp == 2 & l2>1) # & l<nTips changed with ordering
    
#    ind = ind[order(l[ind])]
    ind = ind[order(l2[ind], decreasing = TRUE)]
    
    dm2 <- as.matrix(compatible2(x, x[ind]))
    
    X = as.matrix(x)[,ord]
    Y = X    
    rsY = rowSums(Y)
    X = X[ind, ]
    
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
#        print(fromTo) 
        g = graph(t(res$edge), directed=FALSE)
        
        isChild = (rsY == (Y %*% X[k,]))[index]
        sp2 = NULL
        sp0 = NULL
        
        for(i in 2:length(fromTo)){
#            sptmp = get.shortest.paths(g, fromTo[i-1], fromTo[i], 
#                                       output=c("epath"))$epath[[1]]            
            sptmp = shortest_paths(g, fromTo[i-1], fromTo[i], 
                                       output=c("epath"))$epath[[1]]
            sp2 = c(sp2, sptmp[-c(1, length(sptmp))])
            sp0 = c(sp0, sptmp)
        }
        sp0 = unique(sp0)
        
        if(length(sp2)>0){
            #            blub = which(dm[index[sp2], ind[k]]>0)
            TMP = rowSums(dm2[index[sp2], 1:k, drop=FALSE])
            blub = which(TMP>0)
            sp2 = sp2[blub]
        }
        if(length(sp2)==0){
            isChild = (rsY == (Y %*% X[k,]))[index]  
            sp0 = which(isChild == TRUE)
            edge1 = unique(as.vector(res$edge[sp0,]))
            edge2 = as.vector(res$edge[-sp0,])
            asdf = edge1 %in% edge2
            sp = edge1[asdf]
        }
        if(length(sp2)>0)   sp = unique(as.vector(t(res$edge[sp2,])))     
        parent = res$edge[,1]
        child = res$edge[,2]    
        
        j = ord[which(X[k,]==1)]
        anc = unique(parent[match(j, child)])
        
        maxVert = max(parent)
        l = length(sp)
        
        newVert = (maxVert+1) : (maxVert+l)      
        sp01 = setdiff(sp0, sp2)
        for(i in 1:l) res$edge[sp01,][res$edge[sp01,]==sp[i]] = newVert[i] 
        
        newindex = rep(ind[k], l)        
        if(length(sp)>1)newindex = c(index[sp2], newindex)
        index = c(index, newindex)        
        # connect new and old vertices
        newEdge = matrix(cbind(sp, newVert), ncol=2) 
        if(length(sp)>1){
            # copy edges
            qwer = match(as.vector(res$edge[sp2,]), sp)
            newEdge = rbind(matrix(newVert[qwer], ncol=2), newEdge)
        }
        
        res$edge = rbind(res$edge, newEdge)      
        res$Nnode =  max(res$edge) - nTips
        
        res$splitIndex = index
        res$edge.length <- rep(1, nrow(res$edge))
        class(res) = c("networx", "phylo")
        attr(res, "order") = NULL
    }
    res$edge.length = weight[index]  # ausserhalb
    res$Nnode =  max(res$edge) - nTips
    res$splitIndex = index 
#    attr(res, "splits") = x
    res$splits = x
    class(res) = c("networx", "phylo")
    attr(res, "order") = NULL
    res    
}


#' Phylogenetic networks
#' 
#' \code{as.networx} convert \code{splits} objects into a \code{networx}
#' object. And most important there exists a generic \code{plot} function to plot 
#' phylogenetic network or split graphs.
#' 
#' @details A \code{networx} object hold the information for a phylogenetic network and
#' extends the \code{phylo} object. Therefore some generic function for
#' \code{phylo} objects will also work for \code{networx} objects.  The
#' argument \code{planar = TRUE} will create a planar split graph based on a cyclic
#' ordering. These objects can be nicely plotted in \code{"2D"}. So far not all
#' parameters behave the same on the the \code{rgl} \code{"3D"} and basic graphic \code{"2D"}
#' device.
#' 
#' Often it is easier and safer to supply vectors of graphical parameters for
#' splits (e.g. splits.color) than for edges. These overwrite values edge.color.
#' 
#' @aliases networx 
#' @param x an object of class \code{"splits"} (as.networx) or \code{"networx"}
#' (plot)
#' @param planar logical whether to produce a planar graph from only cyclic
#' splits (may excludes splits).
#' @param coord add coordinates of the nodes, allows to reproduce the plot.
#' @param type "3D" to plot using rgl or "2D" in the normal device.
#' @param use.edge.length a logical indicating whether to use the edge weights
#' of the network to draw the branches (the default) or not.
#' @param show.tip.label a logical indicating whether to show the tip labels on
#' the graph (defaults to \code{TRUE}, i.e. the labels are shown).
#' @param show.edge.label a logical indicating whether to show the tip labels
#' on the graph.
#' @param edge.label an additional vector of edge labels (normally not needed).
#' @param show.node.label a logical indicating whether to show the node labels
#' (see example).
#' @param node.label an additional vector of node labels (normally not needed).
#' @param show.nodes a logical indicating whether to show the nodes (see
#' example).
#' @param tip.color the colors used for the tip labels.
#' @param edge.color the colors used to draw edges.
#' @param edge.width the width used to draw edges.
#' @param edge.lty a vector of line types.
#' @param split.color the colors used to draw edges.
#' @param split.width the width used to draw edges.
#' @param split.lty a vector of line types.
#' @param font an integer specifying the type of font for the labels: 1 (plain
#' text), 2 (bold), 3 (italic, the default), or 4 (bold italic).
#' @param cex a numeric value giving the factor scaling of the labels.
#' @param cex.node.label a numeric value giving the factor scaling of the node
#' labels.
#' @param cex.edge.label a numeric value giving the factor scaling of the edge
#' labels.
#' @param col.node.label the colors used for the node labels.
#' @param col.edge.label the colors used for the edge labels.
#' @param font.node.label the font used for the node labels.
#' @param font.edge.label the font used for the edge labels.
#' @param \dots Further arguments passed to or from other methods.
#' @note The internal representation is likely to change.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{consensusNet}}, \code{\link{neighborNet}},
#' \code{\link{splitsNetwork}}, \code{\link{hadamard}},
#' \code{\link{distanceHadamard}}, \code{\link{layout_with_kk}},
#' \code{\link[ape]{evonet}}, \code{\link[ape]{as.igraph}},
#' \code{\link{densiTree}}
#' @references Dress, A.W.M. and Huson, D.H. (2004) Constructing Splits Graphs
#' \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics
#' (TCBB)}, \bold{1(3)}, 109--115
#' @keywords plot
#' @examples
#' 
#' set.seed(1)
#' tree1 = rtree(20, rooted=FALSE)
#' sp = as.splits(rNNI(tree1, n=10))
#' net = as.networx(sp)
#' plot(net, "2D")
#' \dontrun{
#' # also see example in consensusNet 
#' example(consensusNet)
#' }
#' 
#' @rdname as.networx
#' @export as.networx
as.networx <- function (x, ...) 
{
    if (inherits(x, "networx")) 
        return(x)
    UseMethod("as.networx")
}


getOrdering <- function(x){
    tree = as.phylo(x)
    nTips = length(tree$tip.label)
    ord = reorder(tree)$edge[,2]
    ord = ord[ord<=nTips]
    ind = which(ord == 1L)
    if(ind>1) ord = c(ord[ind:nTips], ord[c(1:(ind-1L))])
    ord  
}

## as.splits.phylo
addTrivialSplits <- function(obj){
    label <- attr(obj, "label")
    nTips <- length(label)
    weight <- attr(obj, "weights")
    if(is.null(weight)) weight = rep(1, length(obj))
    STree = stree(nTips, tip.label = attr(obj, "labels"))
    STree$edge.length=NULL 
    spRes <- as.splits(STree)[STree$edge[,2]]
    tmpIndex = match(spRes, SHORTwise(obj, nTips))
    if(any(is.na(tmpIndex))){
        l.na = sum(is.na(tmpIndex))
        obj <- c(obj, spRes[is.na(tmpIndex)]) 
        weight = c(weight, rep(0, l.na))
        attr(obj, "weights") <- weight
    }
    obj
}


removeTrivialSplits <- function(obj){
    nTips <- length(attr(obj, "label"))
    l <- lengths(obj)
    ind <- which((l == 0L) | (l == 1L) | (l == nTips) | (l == (nTips-1L)))
    obj[-ind]
}


#' @rdname as.networx
#' @method as.networx splits
#' @export
as.networx.splits <- function(x, planar=FALSE, coord = c("none", "2D", "3D"), ...){
  label <- attr(x, "label")
  
  x = addTrivialSplits(x)
  
  nTips <- length(label)
  weight <- attr(x, "weights")
  if(is.null(weight)) weight = rep(1, length(x))
  attr(x, "weights") <- weight
  
  x <- oneWise(x, nTips) 
  l <- lengths(x)
#  l <- sapply(x, length)
  if(any(l==nTips))x <- x[l!=nTips] # get rid of trivial splits
  ext <- sum(l==1 | l==(nTips-1))
  if(!is.null(attr(x, "cycle"))){  
      c.ord <- attr(x, "cycle") 
  }
  else c.ord <- getOrdering(x)
  attr(x, "cycle") = c.ord
  
  dm <- as.matrix(compatible2(x)) 
# which splits are in circular ordering  
    circSplits = which(countCycles(x, ord=c.ord)==2) 
    if(length(circSplits) == length(x)) planar=TRUE
    tmp = circNetwork(x, c.ord)  
    attr(tmp, "order") = NULL
    if(planar){
        return(reorder(tmp))
    }

    ll <- lengths(x)
#    ll <- sapply(x, length)    
    ind <- tmp$splitIndex     # match(sp, x)
    ind2 = union(ind, which(ll==0)) # which(duplicated(x))
    ind2 = union(ind2, which(ll==nTips))
    ord <- order(colSums(dm))
    ord <- setdiff(ord, ind2)
    if(length(ord)>0){    
        for(i in 1:length(ord)){ 
            tmp = addEdge(tmp, x, ord[i])
            tmp$edge.length = weight[tmp$splitIndex]
            tmp$Nnode = max(tmp$edge) - nTips
            class(tmp) = c("networx", "phylo")
        } 
    }
    tmp$edge.length = weight[tmp$splitIndex]
    tmp$Nnode = max(tmp$edge) - nTips
    attr(x, "cycle") <- c.ord
#    attr(tmp, "splits") = x 
    tmp$splits <- x
    class(tmp) = c("networx", "phylo")
    tmp <- reorder(tmp)
    coord <- match.arg(coord)
    vert <- switch(coord,
           "none" = NULL,
           "2D" = coords(tmp, dim="2D"),
           "3D" = coords(tmp, dim="3D"))
#    attr(tmp, "coords") <- coordinates
    tmp$plot <- list(vertices=vert)
    tmp
}


#as.networx.phylo <- function(x, ...){
#    spl <- as.splits(x)
#    as.networx(x, ...)
#}


#' @rdname as.networx
#' @method as.networx phylo
#' @export
as.networx.phylo <- function(x, ...){
    spl <- as.splits(x)
    spl <- spl[x$tree[,2]]
    x$splitIndex <- 1:nrow(x$edge)
#    attr(x, "splits") = spl
    x$splits <- spl
    class(x) <- c("networx", "phylo")
    x
}


#as.igraph.networx <- function(x, directed=FALSE){
#    graph(t(x$edge), directed=directed)
#}




#' Computes a consensusNetwork from a list of trees Computes a \code{networx}
#' object from a collection of splits.
#' 
#' Computes a consensusNetwork, i.e. an object of class \code{networx} from a
#' list of trees, i.e. an class of class \code{multiPhylo}. Computes a
#' \code{networx} object from a collection of splits.
#' 
#' 
#' @param obj An object of class multiPhylo.
#' @param prob the proportion a split has to be present in all trees to be
#' represented in the network.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{consensusNet} returns an object of class networx.  This is
#' just an intermediate to plot phylogenetic networks with igraph.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{splitsNetwork}}, \code{\link{neighborNet}},
#' \code{\link{lento}}, \code{\link{distanceHadamard}},
#' \code{\link{plot.networx}}, \code{\link{maxCladeCred}}
#' @references Holland B.R., Huber K.T., Moulton V., Lockhart P.J. (2004) Using
#' consensus networks to visualize contradictory evidence for species
#' phylogeny. \emph{Molecular Biology and Evolution}, \bold{21}, 1459--61
#' @keywords hplot
#' @examples
#' 
#' data(Laurasiatherian)
#' set.seed(1)
#' bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)nj(dist.hamming(x)), 
#'     bs=50)
#' class(bs) <- 'multiPhylo'
#' cnet = consensusNet(bs, .3)
#' plot(cnet, "2D")
#' \dontrun{
#' library(rgl)
#' open3d()
#' plot(cnet, show.tip.label=FALSE, show.nodes=TRUE)
#' plot(cnet, type = "2D", show.edge.label=TRUE)
#' 
#' tmpfile <- normalizePath(system.file("extdata/trees/RAxML_bootstrap.woodmouse", package="phangorn"))
#' trees <- read.tree(tmpfile)
#' cnet_woodmouse = consensusNet(trees, .3)
#' plot(cnet_woodmouse, type = "2D", show.edge.label=TRUE)
#' }
#' 
#' @export consensusNet
consensusNet <- function (obj, prob = 0.3, ...) 
{
    l = length(obj)
    spl = as.splits(obj)
    w = attr(spl, "weights")
    ind = (w/l) > prob
    spl = spl[ind]
    attr(spl, "confidences") = (w/l)[ind]
#    attr(spl, "weights") = w[ind]
    res = as.networx(spl)  
    res$edge.labels = as.character(res$edge.length / l * 100)
    res$edge.labels[res$edge[,2]<=length(res$tip.label)] = ""
    reorder(res)
}


#' @rdname addConfidences
#' @export
createLabel <- function(x, y, label_y, type="edge", nomatch=NA){
    spl_x <- as.splits(x)
    if(inherits(x, "phylo", TRUE)==1) spl_x <- spl_x[x$edge[,2]]
    spl_y <- as.splits(y)
    if(inherits(y, "phylo", TRUE)==1) spl_y <- spl_y[y$edge[,2]]
    
    tiplabel <- attr(spl_x, "label")
    nTips <- length(tiplabel)
    
    spl_y <- changeOrder(spl_y, tiplabel)
    spl_y <- SHORTwise(spl_y, nTips)
    
    ind <- match(SHORTwise(spl_x, nTips), spl_y)
    pos <-  which(!is.na(ind))

    res <- rep(nomatch, length(spl_x))
    
    if(length(label_y)==1L) label_y <- rep(label_y, length(spl_y))
    res[pos] <- label_y[ind[pos]]
    if(type=="edge" && inherits(x, "networx")){
        return(res[x$splitIndex])
    }
    res  
}




#' Compare splits and add support values to an object
#' 
#' Add support values to a \code{splits}, \code{phylo} or \code{networx}
#' object.
#' 
#' @param x an object of class \code{splits}, \code{phylo} or \code{networx}
#' @param y an object of class \code{splits}, \code{phylo}, \code{multiPhylo}
#' or \code{networx}
#' @param ...  Further arguments passed to or from other methods.
#' @param label_y label of y matched on x. Will be usually of
#' length(as.splits(x)).
#' @param type should labels returned for edges (in \code{networx}) or splits.
#' @param nomatch default value if no match between x and y is found.
#' @return The object \code{x} with added bootstrap / MCMC support values.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{as.splits}}, \code{\link{as.networx}},
#' \code{\link{RF.dist}}, \code{\link{plot.phylo}}
#' @keywords cluster
#' @examples
#' 
#' data(woodmouse)
#' woodmouse <- phyDat(woodmouse)
#' tmpfile <- normalizePath(system.file("extdata/trees/RAxML_bootstrap.woodmouse", package="phangorn"))
#' boot_trees <- read.tree(tmpfile)
#' 
#' dm <- dist.ml(woodmouse)
#' tree <- upgma(dm)
#' nnet <- neighborNet(dm)
#' 
#' tree <- addConfidences(tree, boot_trees)
#' nnet <- addConfidences(nnet, boot_trees)
#' 
#' plot(tree, show.node.label=TRUE)
#' plot(nnet, "2D", show.edge.label=TRUE)
#' 
#' @rdname addConfidences
#' @export 
addConfidences <- function (x, y, ...) UseMethod("addConfidences")


# some function to add confidences on splits if trees have different taxa
#
addConfidencesMultiPhylo <- function(spl, trees){
    
    fun <- function(spl, intersect_labels){
        spl2 <- spl
        
        index <- match(attr(spl, "labels"), intersect_labels)
        attr(spl2, "labels") <- intersect_labels 
        for(i in 1:length(spl2)){
            spl2[[i]] <- sort(na.omit(index[spl[[i]]]))
        }
        l_spl <- lengths(spl2)
        l <- length(intersect_labels)
        ind <- which( (l_spl > 1) & (l_spl < (l-1L)) ) 
        if(length(ind)==0)return(NULL)
        list(spl=spl2[ind], index=ind)
    }
    
    spl_labels <- attr(spl, "labels")
    zaehler <- numeric(length(spl))
    nenner <- numeric(length(spl))
    
    for(i in 1:length(trees)){
#        print(i)
        intersect_labels <- intersect(trees[[i]]$tip.label, spl_labels)
        if(length(intersect_labels) > 3){    
            tmp <- fun(spl, intersect_labels)
            if(!is.null(tmp)){
                tree_spl <- as.splits(trees[[i]])
                if(!identical(intersect_labels, trees[[i]]$tip.label))
                    tree_spl <- fun(tree_spl, intersect_labels)[[1]]
                comp <- compatible_2(as.bitsplits(tmp[[1]]), as.bitsplits(tree_spl))
                #            print(comp)
                ind <- tmp$index
                zaehler[ind] <- zaehler[ind] + comp 
                nenner[ind] <- nenner[ind] + 1L
            }
        }    
        #        print(zaehler)
    }
    confidences <- zaehler / nenner
    attr(spl, "confidences") <- confidences
    spl
}            



# y now more general 
addConfidences.splits <- function(x, y, scaler=1, ...){
    if (hasArg(add)) 
        add <- list(...)$add
    else add <- FALSE
    
    tiplabel <- attr(x, "label")
    nTips = length(tiplabel)
#    x = addTrivialSplits(x) 
    if(inherits(y,"phylo")){
        ind <- match(tiplabel, y$tip.label)
        if (any(is.na(ind)) | length(tiplabel) != length(y$tip.label)) 
            stop("trees have different labels")
        y$tip.label <- y$tip.label[ind]
        ind2 <- match(1:length(ind), y$edge[, 2])
        y$edge[ind2, 2] <- order(ind)
    }
    if(inherits(y, "multiPhylo")){
        if(inherits(try(.compressTipLabel(y), TRUE), "try-error")){
            res <- addConfidencesMultiPhylo(x, y) 
            return(res)
        }
    }
#    browser()
    spl <- as.splits(y)
    spl <- changeOrder(spl, tiplabel)
    spl <- oneWise(spl, nTips)
#    ind <- match(SHORTwise(x, nTips), spl)
#    ind
    ind <- match(oneWise(x, nTips), spl)
    #    pos <-  which(ind > nTips)
    pos <-  which(!is.na(ind))
    confidences <- rep(NA_real_, length(x)) #numeric(length(x))  #character 
    confidences[pos] <- attr(spl, "confidences")[ind[pos]] * scaler
    if(add==TRUE) confidences <- paste(prettyNum(attr(x, "confidences")) , prettyNum(confidences * scaler), sep="/")
    #        y$node.label[ind[pos] - nTips]
    attr(x, "confidences") <- confidences
    x  
}



addConfidences.networx <- function(x, y, scaler=1, ...){
    spl <- x$splits
    spl <- addConfidences(spl, y, scaler=scaler, ...)
    x$splits <- spl
    x    
}

## as.splits.phylo
addConfidences.phylo <- function(x, y, ...){
#    call <- x$call
    if (hasArg(as.is)) 
        as.is <- list(...)$as.is
    else as.is <- TRUE
    nTips <- length(x$tip.label)
    spl <- as.splits(x) %>% oneWise(nTips=nTips)
    conf = attr(addConfidences(spl, y), "confidences")
    l <- lengths(spl)
    if(is.character(conf)) conf <- as.numeric(conf)
    ind <- (l==1L) | (l==(nTips-1L)) | (l==nTips)
    conf[ind==TRUE] <- NA_real_
    nTips = length(x$tip.label)
    if(!as.is) conf <- conf * 100
    x$node.label = conf[-c(1:nTips)]
    x      
} 


# returns order of x$edge
presenceAbsenceOld <- function(x, y){
    X <- as.splits(x)
    Y <- as.splits(y)
    labels <- attr(X, "labels") 
    #    if(inherits(x,"phylo")) X <- X[x$edge[,2]]
    #    if(inherits(y,"phylo")) Y <- Y[y$edge[,2]]
    Y <- changeOrder(Y, labels) # orderSplitLabel
    nTips <- length(labels)
    X <- oneWise(X, nTips)
    Y <- oneWise(Y, nTips)
    res <- match(X, Y)    
    res <- !is.na(res)
    if(inherits(x, "networx")){
        res <- res[x$splitIndex]    
    }    
    if(class(x)[1]=="phylo"){
        # res <- res[x$edge[,2]]
        x$node.label = res[-c(1:length(labels))]
        return(x)
    }
    res            
}


#' @rdname addConfidences
#' @export
presenceAbsence <- function(x,y){
    spl <- as.splits(y)
    l <- length(spl)
    attr(spl, "confidences") <- rep(1, l)
    addConfidences(x, y)
}


reorder.networx <- function (x, order =  "cladewise", index.only = FALSE, ...) 
{
    order <- match.arg(order, c("cladewise", "postorder"))
    if (!is.null(attr(x, "order"))) 
        if (attr(x, "order") == order) 
            return(x)    
    g <- graph(t(x$edge))
#    if(order == "cladewise") neword <- topological.sort(g, "out")
#    else neword <- topological.sort(g, "in") 
    if(order == "cladewise") neword <- topo_sort(g, "out")
    else neword <- topo_sort(g, "in") 
    neworder <- order(match(x$edge[,1], neword))
    if(index.only) return(neworder)
    x$edge <- x$edge[neworder, ]
    if (!is.null(x$edge.length)) 
        x$edge.length <- x$edge.length[neworder]
    if (!is.null(x$edge.labels)) 
        x$edge.labels <- x$edge.labels[neworder]  
    if (!is.null(x$splitIndex))x$splitIndex <- x$splitIndex[neworder]
    attr(x, "order") <- order
    x
}


coords <- function(obj, dim="3D"){
#    if(is.null(attr(obj,"order")) || (attr(obj, "order")=="postorder") ) 
#        obj = reorder.networx(obj)

    l = length(obj$edge.length)
    ind1 = which(!duplicated(obj$splitIndex))

    n = max(obj$edge)
    adj = spMatrix(n, n, i = obj$edge[,2], j = obj$edge[,1], x = rep(1, length(obj$edge.length))) # Matrix::
    g = graph_from_adjacency_matrix(adj, "undirected")
#    g = graph.adjacency(adj, "undirected")
##########
#    add this 
#    g2 <- graph(t(obj$edge), directed=FALSE)
#    g2 <- set.edge.attribute(g, "weight", value=rep(1, nrow(obj$edge))
    if(dim=="3D"){
        coord <- layout_nicely(g, dim=3)
#        coord <- layout_with_kk(g, dim=3)
#        coord <- layout.kamada.kawai(g, dim=3)
        k = matrix(0, max(obj$splitIndex), 3)
        for(i in ind1){
            tmp = coord[obj$edge[i, 2],] - coord[obj$edge[i, 1],]
            k[obj$splitIndex[i], ] = kart2kugel(tmp[1], tmp[2], tmp[3])
        }
        k[obj$splitIndex[ind1],1] = obj$edge.length[ind1] 

        res = matrix(0, vcount(g), 3)
        for(i in 1:l){# unique(obj$splitIndex)
            j = obj$edge[i,1]
            m = obj$edge[i,2]
            p = obj$splitIndex[i]
            res[m,] = res[j,] + kugel2kart(k[p,1], k[p,2], k[p,3])     
        }            
    }
    else{
        coord <- layout_nicely(g, dim=2)
#        coord <- layout_with_kk(g, dim=2)
#        coord <- layout.kamada.kawai(g, dim=2)
        k = matrix(0, max(obj$splitIndex), 2)
        for(i in ind1){
            tmp = coord[obj$edge[i, 2],] - coord[obj$edge[i, 1],]
            k[obj$splitIndex[i], ] = kart2kreis(tmp[1], tmp[2])
        }
        k[obj$splitIndex[ind1],1] = obj$edge.length[ind1] 
        res = matrix(0, vcount(g), 2)
        for(i in 1:l){# unique(obj$splitIndex)
            j = obj$edge[i,1]
            m = obj$edge[i,2]
            p = obj$splitIndex[i]
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



#' @rdname as.networx
#' @method plot networx
#' @export
plot.networx <- function(x, type="3D", use.edge.length = TRUE, show.tip.label=TRUE,
    show.edge.label=FALSE, edge.label=NULL, show.node.label = FALSE, node.label=NULL,
    show.nodes=FALSE, tip.color = "black", 
    edge.color="black", edge.width = 3, edge.lty = 1,
    split.color=NULL, split.width = NULL, split.lty = NULL,
    font = 3, cex = par("cex"), 
    cex.node.label=cex, cex.edge.label=cex,
    col.node.label = tip.color, col.edge.label = tip.color, 
    font.node.label = font, font.edge.label = font,
    ...){
    type = match.arg(type, c("3D", "2D")) 
    if(use.edge.length==FALSE) x$edge.length[] = 1
# test    
#    x = reorder(x)
    nTips = length(x$tip.label)
    conf = attr(x$splits,"confidences") 
    index = x$splitIndex
    if(is.null(edge.label) & !is.null(conf)){
        conf = conf[index]
        if(!is.null(x$translate)) conf[match(x$translate$node, x$edge[,2])] = ""
        else conf[x$edge[,2] <= nTips] = ""    
        edge.label = conf
    }
    if(is.null(node.label))node.label = as.character(1:max(x$edge))
    if(show.tip.label)node.label[1:nTips] = ""
    
    lspl <- max(x$splitIndex)
    if(!is.null(split.color)){
        if(length(split.color)!=lspl) stop("split.color must be same length as splits")
        else edge.color <- split.color[x$splitIndex]
    } 
    if(!is.null(split.width)){
        if(length(split.width)!=lspl) stop("split.color must be same length as splits")
        else edge.width <- split.width[x$splitIndex]
    } 
    if(!is.null(split.lty)){
        if(length(split.lty)!=lspl) stop("split.color must be same length as splits")
        else edge.lty <- split.lty[x$splitIndex]
    } 
    
    chk <- FALSE
    
    if(type=="3D") chk <- requireNamespace("rgl", quietly = TRUE) #.check.pkg("rgl")
    if(!chk && type=="3D"){
        warning("type=\"3D\" requires the package \"rgl\"\n, plotting =\"2D\" instead!\n")
        type="2D"
    }
    # use precomputed vertices when available
    coord <- NULL
    if(!is.null(x$.plot)) coord <- x$.plot$vertices
    
    if(type=="3D") {
        if(is.null(coord) || ncol(coord)!=3)        
            coord <- coords(x, dim="3D")
        plotRGL(coord, x, show.tip.label=show.tip.label, show.edge.label=show.edge.label, 
             edge.label = edge.label, show.node.label = show.node.label, node.label=node.label, 
             show.nodes=show.nodes, tip.color = tip.color, edge.color=edge.color, 
             edge.width = edge.width, font = font, cex = cex, 
             cex.node.label=cex.node.label, cex.edge.label=cex.edge.label,
             col.node.label = col.node.label, col.edge.label = col.edge.label,
             font.node.label = font.node.label, font.edge.label = font.edge.label)
    }
    else{
        if(is.null(coord) || ncol(coord)!=2)
     	    coord <- coords(x, dim="2D")
	    plot2D(coord, x, show.tip.label=show.tip.label, show.edge.label=show.edge.label, 
	        edge.label = edge.label, show.node.label = show.node.label, node.label=node.label,
	        show.nodes=show.nodes, tip.color = tip.color, edge.color=edge.color,
	        edge.width = edge.width, edge.lty=edge.lty,font = font, cex = cex, 
	        cex.node.label=cex.node.label, cex.edge.label=cex.edge.label,
	        col.node.label = col.node.label, col.edge.label = col.edge.label,
	        font.node.label = font.node.label, font.edge.label = font.edge.label,
	        add=FALSE)
    }   
    x$.plot <- list(vertices = coord, edge.color=edge.color, edge.width=edge.width, edge.lty = edge.lty)
    invisible(x)
}

    
plotRGL <- function(coords, net, show.tip.label=TRUE, 
        show.edge.label=FALSE, edge.label=NULL, show.node.label=FALSE, node.label=NULL,
        show.nodes=FALSE, tip.color = "blue", edge.color="grey", 
        edge.width = 3, font = 3, cex = par("cex"), 
        cex.node.label=cex,  cex.edge.label=cex,
        col.node.label=tip.color, col.edge.label=tip.color,
        font.node.label=font, font.edge.label=font,        
        ...){
    
#    chk <- .check.pkg("rgl")
#    if(!chk) open3d <- segments3d <- spheres3d <- rgl.texts <- function(...)NULL

    open3d <- rgl::open3d 
    segments3d  <- rgl::segments3d
    spheres3d <- rgl::spheres3d 
    rgl.texts <- rgl::rgl.texts
        
    edge = net$edge
  
    x = coords[,1]
    y = coords[,2]
    z = coords[,3]
     
    nTips = length(net$tip.label)
  
    segments3d(x[t(edge)],y[t(edge)],z[t(edge)], col=rep(edge.color, each=2), lwd=edge.width) 
    radius=0
    if(show.nodes){
        radius = sqrt((max(x)-min(x))^2 + (max(y)-min(y))^2 + (max(z)-min(z))^2) / 200    
        spheres3d(x[1:nTips], y[1:nTips],z[1:nTips], radius=2*radius, color="cyan")
        spheres3d(x[-c(1:nTips)], y[-c(1:nTips)],z[-c(1:nTips)], radius=radius, color="magenta")
    }
    if(show.tip.label){
        if(is.null(net$translate))   
        rgl.texts(x[1:nTips]+2.05*radius,y[1:nTips],z[1:nTips],net$tip.label, color=tip.color, cex=cex, font=font)
        else
        rgl.texts(x[net$translate$node]+2.05*radius,y[net$translate$node],z[net$translate$node],net$tip.label, color=tip.color, cex=cex, font=font)    
    }
    if(show.edge.label){
	    ec = edgeLabels(x, y, z, edge)
      if(is.null(edge.label)) edge.label = net$splitIndex
        #else edge.label = net$splitIndex    
	    rgl.texts(ec[,1], ec[,2], ec[,3], edge.label, color=col.edge.label, cex=cex.edge.label, font=font.edge.label)     
    } 
    if(show.node.label){
        rgl.texts(x, y, z, node.label, color=col.node.label, cex=cex.node.label, font=font.node.label) 
    }
}


plot2D <- function(coords, net, show.tip.label=TRUE,  
       show.edge.label=FALSE, edge.label=NULL, show.node.label=FALSE, node.label=NULL,
       tip.color = "blue", edge.color="grey",                   
       edge.width = 3, edge.lty=1, 
       font = 3, cex = par("cex"), 
       cex.node.label=cex,  cex.edge.label=cex,
       col.node.label=tip.color, col.edge.label=tip.color,
       font.node.label=font, font.edge.label=font,
       add=FALSE, ...){
   edge = net$edge
   label = net$tip.label
   xx = coords[,1]
   yy = coords[,2]
   nTips = length(label)

#   cex=1
   
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
   cladogram.plot(edge, xx, yy, edge.color, edge.width, edge.lty)
   if(show.tip.label){
        if(is.null(net$translate)) ind=match(1:nTips, edge[,2])
        else ind=match(net$translate$node, edge[,2])
        pos = rep(4, nTips)
        XX <- xx[edge[ind, 1]] - xx[edge[ind, 2]]
        pos[XX>0] = 2
        YY <- yy[edge[ind, 1]] - yy[edge[ind, 2]]
        pos2 <- rep(3, nTips)
        pos2[YY>0] = 1
# needed if tiplabels are not at internal nodes        
        XX[is.na(XX)] = 0
        YY[is.na(YY)] = 0
        pos[abs(YY)>abs(XX)] <- pos2[abs(YY)>abs(XX)] 	
        if(is.null(net$translate)) text(xx[1:nTips], yy[1:nTips], labels=label, pos=pos, col=tip.color, cex=cex, font=font)
        else text(xx[net$translate$node], yy[net$translate$node], labels=label, pos=pos, col=tip.color, cex=cex, font=font)
    }
    if(show.edge.label){
	    ec = edgeLabels(xx,yy, edge=edge)
	    if(is.null(edge.label))edge.label = net$splitIndex
	    
	    # show only one edge label
	    em <- apply(ec, 1, function(x)max(abs(x)))
	    si <- net$splitIndex
	    for(i in unique(si)){
	        tmp <- si==i
	        if(sum(tmp)>1){
	            w <- which(tmp)
	            wm <- which.max(em[w])
	            edge.label[w[-wm]] <- ""
	        }
	    }
	    
	    text(ec[,1], ec[,2], labels=edge.label, col=col.edge.label, cex=cex.edge.label, font=font.edge.label)     
	} 
    if(show.node.label){
         text(xx, yy, labels=node.label, col=col.node.label, cex=cex.node.label, font=font.node.label)    
    }   
}   
   



