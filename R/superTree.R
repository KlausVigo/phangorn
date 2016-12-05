# now more memoryefficient
# from phytools code by Liam Revell with a few changes
my.supertree<-function(trees, trace=0, ...){
    # some minor error checking
    if(!inherits(trees,"multiPhylo")) stop("trees must be object of class 'multiPhylo.'")
    
    labels <- lapply(trees, function(x)sort(x$tip.label))
    ulabels <- unique(labels)
    lul <- length(ulabels)
    # compute matrix representation phylogenies
    X<-vector("list", lul) # list of bipartitions
    characters<-0 # number of characters
    weights <- NULL
    species<-trees[[1]]$tip.label
    for(i in 1:lul){
        pos <- match(labels, ulabels[i])
        ind <- which(!is.na(pos))  
        temp<-prop.part(trees[ind]) # find all bipartitions
        # create matrix representation of trees[[i]] in X[[i]]
        X[[i]]<-matrix(0,nrow=length(trees[[ind[1]]]$tip.label),ncol=length(temp)-1)
        for(j in 1:ncol(X[[i]])) X[[i]][c(temp[[j+1]]),j]<-1
        rownames(X[[i]])<-attr(temp,"labels") # label rows
        #    if(i==1) species<-trees[[ind[1]]]$tip.label
        #    else 
        species<-union(species,trees[[ind[1]]]$tip.label) # accumulate labels
        characters<-characters+ncol(X[[i]]) # count characters
        weights <- c(weights, attr(temp, "number")[-1])
    }
    XX<-matrix(data="?",nrow=length(species),ncol=characters,dimnames=list(species))
    j<-1
    for(i in 1:length(X)){
        # copy each of X into supermatrix XX
        XX[rownames(X[[i]]),c(j:((j-1)+ncol(X[[i]])))]<-X[[i]][1:nrow(X[[i]]),1:ncol(X[[i]])]
        j<-j+ncol(X[[i]])
    }
    # compute contrast matrix for phangorn
    contrast<-matrix(data=c(1,0,0,1,1,1),3,2,dimnames=list(c("0","1","?"),c("0","1")),byrow=TRUE)
    # convert XX to phyDat object
    XX<-phyDat(XX,type="USER",contrast=contrast, compress=FALSE) 
    attr(XX, "weight") <- weights 
    # estimate supertree
    supertree<-pratchet(XX,all=TRUE, trace=trace, ...)
    return(supertree)
}


# Robinson-Foulds supertree
fun.rf <- function(x, tree) sum(RF.dist(x, tree))
fun.spr <- function(x, tree) sum(SPR.dist(x, tree))

dist.superTree <- function(tree, trace=0, fun, start=NULL, multicore=FALSE, mc.cores = NULL){
    if(multicore && is.null(mc.cores)){
        mc.cores <- detectCores()
    }
    if(is.null(start))start <- superTree(tree, rooted=FALSE)
    if(inherits(start, "multiPhylo")) start <- start[[1]]
    best_tree <- unroot(start)
    best <- fun(best_tree, tree)
    if(trace>0)cat("best score so far:", best, "\n")
    eps <- TRUE
    while(eps){
        nni_trees <- nni(best_tree) 
        if (multicore) {
            tmp <- mclapply(nni_trees, fun, tree, mc.cores = mc.cores)
            tmp <-unlist(tmp)
        }
        else tmp <- sapply(nni_trees, fun, tree)
        if(min(tmp) < best){
            ind <- which.min(tmp)
            best_tree <- nni_trees[[ind]]
            best <- tmp[ind]
            if(trace>0)cat("best score so far:", best, "\n")
        } 
        else eps <- FALSE
    }
    attr(best_tree, "score") <- best
    best_tree
}


# we want a rooted supertree


#' Super Tree methods
#' 
#' These function \code{superTree} allows the estimation of a supertree from a
#' set of trees using either Matrix representation parsimony, Robinson-Foulds
#' or SPR as criterion.
#' 
#' The function \code{superTree} extends the function mrp.supertree from Liam
#' Revells, with artificial adding an outgroup on the root of the trees.  This
#' allows to root the supertree afterwards. The functions is internally used in
#' DensiTree. The implementation for the RF- and SPR-supertree are very basic
#' so far and assume that all trees share the same set of taxa.
#' 
#' @param tree an object of class \code{multiPhylo}
#' @param method An argument defining which algorithm is used to optimize the
#' tree.  Possible are "MRP", "NNI", and "SPR".
#' @param rooted should the resulting supertrees be rooted.
#' @param trace defines how much information is printed during optimization.
#' @param \dots further arguments passed to or from other methods.
#' @return The function returns an object of class \code{phylo}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com} Liam Revell
#' @seealso \code{mrp.supertree}, \code{\link{densiTree}},
#' \code{\link{RF.dist}}, \code{\link{SPR.dist}}
#' @references Baum, B. R., (1992) Combining trees as a way of combining data
#' sets for phylogenetic inference, and the desirability of combining gene
#' trees. \emph{Taxon}, \bold{41}, 3-10.
#' 
#' Ragan, M. A. (1992) Phylogenetic inference based on matrix representation of
#' trees. \emph{Molecular Phylogenetics and Evolution}, \bold{1}, 53-58.
#' @keywords cluster
#' @examples
#' 
#' data(Laurasiatherian)
#' set.seed(1)
#' bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)upgma(dist.hamming(x)), bs=50)
#' class(bs) <- 'multiPhylo'
#' 
#' mrp_st <- superTree(bs, rooted=TRUE)
#' plot(superTree(mrp_st))
#' \dontrun{
#' rf_st <- superTree(bs, method = "RF")
#' spr_st <- superTree(bs, method = "SPR")
#' }
#' 
#' @export superTree
superTree = function(tree, method="MRP", rooted=FALSE, trace=0, ...){
    fun = function(x){
        x=reorder(x, "postorder")
        nTips = length(x$tip.label)
        x$edge[x$edge>nTips] = x$edge[x$edge>nTips] + 2L
        l=nrow(x$edge)
        oldroot = x$edge[l,1L]
        x$edge=rbind(x$edge,matrix(c(rep(nTips+2,2),oldroot,nTips+1),2L,2L))
        x$edge.length=c(x$edge.length, 100, 100)
        x$tip.label=c(x$tip.label, "ZZZ")
        x$Nnode=x$Nnode+1L
        x
    }
    if(!is.null(attr(tree, "TipLabel")))tree = .uncompressTipLabel(tree)
    tree = unclass(tree)
    if(rooted) tree = lapply(tree, fun)    
    class(tree)="multiPhylo"
    res = my.supertree(tree, trace=trace, ...)
    if(rooted){
        if(inherits(res,"multiPhylo")){
            res = lapply(res, root, "ZZZ")
            res = lapply(res, drop.tip, "ZZZ")  
            class(res) = "multiPhylo"
        }
        else{
            res = root(res, "ZZZ")
            res = drop.tip(res, "ZZZ")  
        }
    }
    if(inherits(res,"multiPhylo")){
        fun = function(x){
            x$edge.length <- rep(.1, nrow(x$edge)) 
            x
        }
        res <- lapply(res, fun)
        res <- lapply(res, reorder, "postorder")
        class(res) = "multiPhylo"
    }       
    else{ 
        res$edge.length = rep(.1, nrow(res$edge))
        res <- reorder(res, "postorder")
    }
    if(method=="MRP") return(res)
    
    tree <- unroot(tree)
    tree <- reorder(tree, "postorder")
#    tree <- lapply(tree, unroot)
#    tree <- lapply(tree, reorder, "postorder")
#    class(tree) = "multiPhylo"
    
    if(method=="RF") res <- dist.superTree(tree, trace=trace, fun.rf, start=res)
    if(method=="SPR") res <- dist.superTree(tree, trace=trace, fun.spr, start=res)   
    res
}
