#
# ancestral sequences ML
#


#' Ancestral character reconstruction.
#' 
#' Marginal reconstruction of the ancestral character states.
#' 
#' The argument "type" defines the criterion to assign the internal nodes. For
#' \code{ancestral.pml} so far "ml" and (empirical) "bayes" and for
#' \code{ancestral.pars} "MPR" and "ACCTRAN" are possible.
#' 
#' With parsimony reconstruction one has to keep in mind that there will be
#' often no unique solution.
#' 
#' For further details see vignette("Ancestral").
#' 
#' @param object an object of class pml
#' @param tree a tree, i.e. an object of class pml
#' @param data an object of class phyDat
#' @param type method used to assign characters to internal nodes, see details.
#' @param i plots the i-th site pattern of the \code{data}.
#' @param col a vector containing the colors for all possible states.
#' @param cex.pie a numeric defining the size of the pie graphs
#' @param pos a character string defining the position of the legend
#' @param cost A cost matrix for the transitions between two states.
#' @param return return a \code{phyDat} object or matrix of probabilities. 
#' @param \dots Further arguments passed to or from other methods.
#' @return %A matrix containing the the estimates character states. An object
#' of class "phyDat", containing the ancestral states of all nodes.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}}, \code{\link{parsimony}}, \code{\link[ape]{ace}},
#' \code{\link[ape]{root}}
#' @references Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer
#' Associates, Sunderland.
#' 
#' Swofford, D.L., Maddison, W.P. (1987) Reconstructing ancestral character
#' states under Wagner parsimony. \emph{Math. Biosci.} \bold{87}: 199--229
#' 
#' Yang, Z. (2006). \emph{Computational Molecular evolution}. Oxford University
#' Press, Oxford.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' example(NJ)
#' fit = pml(tree, Laurasiatherian)
#' anc.ml = ancestral.pml(fit, type = "ml")
#' anc.p = ancestral.pars(tree, Laurasiatherian)
#' \dontrun{
#' require(seqLogo)
#' seqLogo( t(subset(anc.ml, 48, 1:20)[[1]]), ic.scale=FALSE)
#' seqLogo( t(subset(anc.p, 48, 1:20)[[1]]), ic.scale=FALSE)
#' }
#' # plot the first site pattern
#' plotAnc(tree, anc.ml, 1)
#' # plot the third character 
#' plotAnc(tree, anc.ml, attr(anc.ml, "index")[3])
#' 
#' @rdname ancestral.pml
#' @export 
ancestral.pml <- function (object, type=c("marginal", "ml", "bayes"), return="prob") 
{
    call <- match.call()
    pt <- match.arg(type, c("marginal", "joint", "ml", "bayes"))   
    tree = object$tree 
    
    INV <- object$INV
    inv <- object$inv
    
    data = getCols(object$data, tree$tip.label) 
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    nTips = length(tree$tip.label)
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
    
    if(return!="phyDat")result = new2old.phyDat(data) 
    else result[1:nTips] = data
    eig = object$eig
    
    bf = object$bf
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip.label))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(data, "contrast")
# proper format    
    eps <- 1.0e-5
    ind1 <- which( apply(contrast, 1, function(x)sum(x > eps)) == 1L)
    ind2 <- which( contrast[ind1, ] > eps, arr.ind = TRUE)
    pos <- ind2[match(as.integer(1L:ncol(contrast)),  ind2[,2]),1]
    
    nco = as.integer(dim(contrast)[1])
    for(i in 1:l)dat[i,(nTips + 1):m] <- .Call("LogLik2", data, P[i,], nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
    
    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    nTips = min(parent) - 1
# in C with scaling    
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
# scaling!!!            
            tmp = tmp + w[i] * dat[[i, j]]                                 
        }
        if ( (pt == "bayes") || (pt == "marginal")) tmp = tmp * rep(bf, each=nr)
        tmp = tmp / rowSums(tmp)
        
        if(return=="phyDat") tmp <- pos[apply(tmp, 1, which.max)]
        
        result[[j]] = tmp
    } 
    attributes(result) = x
    attr(result, "call") <- call
    result 
}


#joint_reconstruction <- function(object){
#
#}
    

# in ancestral.pml
ancestral2phyDat <- function(x) {
    eps <- 1.0e-5
    contr <- attr(x, "contrast")
    # a bit too complicated    
    ind1 <- which( apply(contr, 1, function(x)sum(x > eps)) == 1L)
    ind2 <- which( contr[ind1, ] > eps, arr.ind = TRUE)
    pos <- ind2[match(as.integer(1L:ncol(contr)),  ind2[,2]),1]
    # only first hit    
    res <- lapply(x, function(x, pos) pos[apply(x, 1, which.max)], pos)
    attributes(res) <- attributes(x)
    return(res)
}


# in ancestral.pml
# variante fuer parsimony und ambiguous DNA 


# raus ??
fast.tree  = function(tree, node){
    parent = c(node, Ancestors(tree, node))
    children = Descendants(tree, parent, 'children')
    l = lengths(children)
#    l = sapply(children, length)
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
