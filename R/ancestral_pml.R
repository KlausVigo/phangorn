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
#' fit <- pml(tree, Laurasiatherian)
#' anc.ml <- ancestral.pml(fit, type = "ml")
#' anc.p <- ancestral.pars(tree, Laurasiatherian)
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
ancestral.pml <- function(object, type="marginal", return="prob") 
{
    call <- match.call()
    pt <- match.arg(type, c("marginal", "joint", "ml", "bayes"))   
    tree <- object$tree 
    
    INV <- object$INV
    inv <- object$inv
    
    data <- getCols(object$data, tree$tip.label) 
    data_type <- attr(data, "type")    
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    nTips <- length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m <- length(edge) + 1  # max(edge)
    w <- object$w
    g <- object$g
    l <- length(w)    
    nr <- attr(data, "nr")
    nc <- attr(data, "nc")
    dat <- vector(mode = "list", length = m*l)
    result <- vector(mode = "list", length = m)
    dim(dat) <- c(l,m)
    
    x <- attributes(data)
    label <- as.character(1:m)
    nam <- tree$tip.label
    label[seq_along(nam)] <- nam
    x[["names"]] <- label
    

    tmp <- length(data)
    
    if(return!="phyDat")result <- new2old.phyDat(data) 
    else result[1:nTips] <- data
    eig <- object$eig
    
    bf <- object$bf
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
    node <- as.integer(node - min(node))
    edge <- as.integer(edge - 1)
    nTips <- as.integer(length(tree$tip.label))
    mNodes <- as.integer(max(node) + 1)
    contrast <- attr(data, "contrast")
# proper format    
    eps <- 1.0e-5
    ind1 <- which( apply(contrast, 1, function(x)sum(x > eps)) == 1L)
    ind2 <- which( contrast[ind1, ] > eps, arr.ind = TRUE)
    pos <- ind2[match(as.integer(1L:ncol(contrast)),  ind2[,2]),1]
    
    nco <- as.integer(dim(contrast)[1])
    for(i in 1:l)dat[i,(nTips + 1):m] <- .Call("LogLik2", data, P[i,], nr, nc, 
                node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
    
    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    nTips <- min(parent) - 1
# in C with scaling    
    for(i in 1:l){     
        for (j in (m - 1):1) {
            if (child[j] > nTips){
                tmp2 <- (dat[[i, parent[j]]]/(dat[[i,child[j]]] %*% P[[i,j]]))
                dat[[i, child[j]]] <- (tmp2 %*% P[[i,j]]) * dat[[i, child[j]]]  
            }
        }
    }
    for (j in unique(parent)) {
        tmp <- matrix(0, nr, nc)
        if(inv>0) tmp <- as.matrix(INV) * inv
        for(i in 1:l){  
# scaling!!!            
            tmp <- tmp + w[i] * dat[[i, j]]                                 
        }
        if ((pt == "bayes") || (pt == "marginal")) tmp <- tmp * rep(bf, each=nr)
        tmp <- tmp / rowSums(tmp)
        
        if(return=="phyDat"){
            if(data_type=="DNA"){
                tmp <- p2dna(tmp) # prob2fitchCoding(tmp)
                tmp <- fitchCoding2ambiguous(tmp)
            }
            else tmp <- pos[max.col(tmp)]   # [apply(tmp, 1, which.max)]
        }
        result[[j]] <- tmp
    } 
    attributes(result) <- x
    attr(result, "call") <- call
    result 
}


#joint_reconstruction <- function(object){
#
#}
    

# in ancestral.pml and ancestral.pars
ancestral2phyDat <- function(x) {
    eps <- 1.0e-5
    contr <- attr(x, "contrast")
    # a bit too complicated    
    ind1 <- which( apply(contr, 1, function(x)sum(x > eps)) == 1L)
    ind2 <- which( contr[ind1, ] > eps, arr.ind = TRUE)
    pos <- ind2[match(as.integer(1L:ncol(contr)),  ind2[,2]),1]
    # only first hit    
#    res <- lapply(x, function(x, pos) pos[apply(x, 1, which.max)], pos)
    res <- lapply(x, function(x, pos) pos[max.col(x)], pos)
    attributes(res) <- attributes(x)
    return(res)
}


# in ancestral.pml
# variante fuer parsimony und ambiguous DNA, ersetzt durch p2dna
prob2fitchCoding <- function(x, eps=0.999){
    row_max <- apply(x, 1, max)
    x <- x / row_max
    as.vector((x>eps) %*% c(1L, 2L, 4L, 8L))
}


fitchCoding2ambiguous <- function(x, type="DNA"){
    y <- c(1L, 2L, 4L, 8L, 8L, 3L, 5L, 9L, 6L, 10L, 12L, 7L, 11L, 13L, 
           14L, 15L, 15L, 15L)
    fmatch(x, y)
}    


fitchCoding2ambiguous2 <- function(x, type="DNA"){
    y <- c(1L, 2L, 4L, 8L, 8L, 3L, 5L, 9L, 6L, 10L, 12L, 7L, 11L, 13L, 14L, 15L)
    dna <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y", "k", "v", "h", 
             "d", "b", "n")    
    rna <- c("a", "c", "g", "u", "u", "m", "r", "w", "s", "y", "k", "v", "h", 
             "d", "b", "n") 
    res <- switch(type, 
                  "DNA" = dna[match(x,y)],
                  "DNA" = dna[match(x,y)])
    res
}


#' @rdname ancestral.pml
#' @export
ancestral.pars <- function (tree, data, type = c("MPR", "ACCTRAN"), cost=NULL,
                            return="prob") 
{
    call <- match.call()
    type <- match.arg(type)
    if (type == "ACCTRAN") 
        res <- ptree(tree, data, retData = TRUE)[[2]]
    #    if (type == "MPR") 
    #        res = mpr(tree, data)
    if (type == "MPR"){ 
        res <- mpr(tree, data, cost=cost, return=return)
        attr(res, "call") <- call
        return(res)
    }
    l <- length(tree$tip.label)
    data_type <- attr(data, "type")
    x <- attributes(data)
    m <- dim(res)[2]
    label <- as.character(1:m)
    nam <- tree$tip.label
    label[seq_along(nam)] <- nam
    x[["names"]] <- label
    
    nc <- attr(data, "nc")
    result <- vector("list", m) 
    
    if(data_type=="DNA"){
        result[1:l] <- subset(data, nam)
        for(i in (l+1) : m){
            result[[i]] <- fitchCoding2ambiguous(res[,i])
        }
        # phyDat needs new2old <-
        #if(return=="prob") tmp <- contrast[tmp, , drop=FALSE]
    }
    else {
        Z <- unique(as.vector(res))
        tmp <- t(sapply(Z, function(x)dec2bin(x, nc)))
        tmp <- tmp / rowSums(tmp)
        dimnames(tmp) <- list(Z, attr(data, "levels"))
        for(i in seq_len(m)){ 
        result[[i]] <- tmp[as.character(res[,i]),,drop=FALSE]
        rownames(result[[i]]) <- NULL
        }
    }
    attributes(result) <- x
    attr(result, "call") <- call
    if(data_type!="DNA" & return=="phyDat"){
        contrast <- attr(data, "contrast")
        # proper format    
        eps <- 1.0e-5
        ind1 <- which( apply(contrast, 1, function(x)sum(x > eps)) == 1L)
        ind2 <- which( contrast[ind1, ] > eps, arr.ind = TRUE)
        pos <- ind2[match(as.integer(1L:ncol(contrast)),  ind2[,2]),1]
        result[1:l] <- subset(data, nam)
        for(i in (l+1) : length(result)) {
            tmp <- result[[i]]
            result[[i]] <- pos[max.col(tmp)] 
        }    
    }
    if(data_type=="DNA" & return=="prob"){
        result <- new2old.phyDat(result)
        for(i in seq_len(length(result))) {
            tmp <- result[[i]]
            tmp <- tmp / rowSums(tmp)
            result[[i]] <- tmp
        }
    }    
    result
}


#' @rdname ancestral.pml
#' @export
pace <- ancestral.pars


mpr.help <- function (tree, data, cost=NULL) 
{   
    tree<- reorder(tree, "postorder")     
    if (!inherits(data,"phyDat")) 
        stop("data must be of class phyDat")    
    levels <- attr(data, "levels")
    l <- length(levels)
    if (is.null(cost)) {
        cost <- matrix(1, l, l)
        cost <- cost - diag(l)
    }   
    weight <- attr(data, "weight")
    p <- attr(data, "nr")
    kl <- TRUE
    i <- 1
    dat <- prepareDataSankoff(data)
    for (i in seq_along(dat)) storage.mode(dat[[i]]) <- "double"    
    tmp <- fit.sankoff(tree, dat, cost, returnData='data')
    p0 <- tmp[[1]]    
    datf <- tmp[[2]]
    datp <- pnodes(tree, datf, cost) 
    
    nr <- attr(data, "nr")
    nc <- attr(data, "nc")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    
    node <- as.integer(node - 1)      
    edge <- as.integer(edge - 1) 
    
    res <- .Call("sankoffMPR", datf, datp, as.numeric(cost), as.integer(nr),
                 as.integer(nc), node, edge, PACKAGE="phangorn")    
    root <- getRoot(tree)
    res[[root]] <- datf[[root]]
    res
}


mpr <- function(tree, data, cost=NULL, return="prob"){
    data <- subset(data, tree$tip.label)
    att <- attributes(data)
    type <- att$type
    nr <- att$nr
    nc <- att$nc
    res <- mpr.help(tree,data,cost)
    l <- length(tree$tip.label)
    m <- length(res)
    label <- as.character(1:m)
    nam <- tree$tip.label
    label[seq_along(nam)] <- nam
    att[["names"]] <- label
    ntips <- length(tree$tip.label)
    contrast <- att$contrast
    eps <- 5e-6
    rm <- apply(res[[ntips+1]], 1, min)
    RM <- matrix(rm,nr, nc) + eps

#    if(return!="prob" & type=="DNA"){
#        for(i in (ntips+1):m){
#            tmp <- as.numeric(res[[i]] < RM)
#            tmp <- prob2fitchCoding(tmp)
#            tmp <- fitchCoding2ambiguous(tmp)
#            res[[i]] <- tmp
#        }
#        browser()
#        res[1:ntips] <- data
#        attributes(res) <- att
#        return(res)
#    }
    fun <- function(X){
        rs <- rowSums(X) #apply(X, 1, sum)
        X / rs
    }
    for(i in 1:ntips) res[[i]] <- contrast[data[[i]],,drop=FALSE]
    for(i in (ntips+1):m) res[[i]][] <- as.numeric(res[[i]] < RM)
    if(return=="prob"){
#        for(i in 1:ntips) res[[i]] <- contrast[data[[i]],,drop=FALSE]
        if(return=="prob") res <- lapply(res, fun)
    }
#    else res[1:ntips] <- data[1:ntips]
    attributes(res) <- att
    fun2 <- function(x){
        x <- p2dna(x) # prob2fitchCoding(x)
        fitchCoding2ambiguous(x)
    }
    if(return!="prob"){
        if(type=="DNA"){
            res <- lapply(res, fun2)
            attributes(res) <- att
        }    
        else res <- ancestral2phyDat(res)
        res[1:ntips] <- data
    }    
    res
}


#' @rdname ancestral.pml
#' @param site.pattern logical, plot i-th site pattern or i-th site
#' @export
plotAnc <- function (tree, data, i = 1, site.pattern = TRUE, col=NULL, 
                     cex.pie=par("cex"), pos="bottomright", ...)
{
    y <- subset(data, , i, site.pattern=site.pattern)
    #   args <- list(...)
    #   CEX <- if ("cex" %in% names(args))
    #       args$cex 
    #   else par("cex")
    CEX <- cex.pie
    xrad <- CEX * diff(par("usr")[1:2])/50
    levels <- attr(data, "levels")
    nc <- attr(data, "nc")
    y <- matrix(unlist(y[]), ncol = nc, byrow = TRUE)
    l <- dim(y)[1]
    dat <- matrix(0, l, nc)
    for (i in 1:l) dat[i, ] <- y[[i]]
    plot(tree, label.offset = 1.1 * xrad, plot = FALSE, ...)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    XX <- lastPP$xx
    YY <- lastPP$yy
    xrad <- CEX * diff(lastPP$x.lim * 1.1)/50
    par(new = TRUE)
    plot(tree, label.offset = 1.1 * xrad, plot = TRUE, ...)
    if(is.null(col)) col <- rainbow(nc)
    if(length(col)!=nc) 
        warning("Length of color vector differs from number of levels!")
    BOTHlabels(pie = y, XX = XX, YY = YY, adj = c(0.5, 0.5),
               frame = "rect", pch = NULL, sel = seq_along(XX), thermo = NULL,
               piecol = col, col = "black", bg = "lightblue", horiz = FALSE,
               width = NULL, height = NULL, cex=cex.pie)
    if(!is.null(pos))legend(pos, levels, text.col = col)
}

