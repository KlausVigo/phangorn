fitch <- function (tree, data, site="pscore") 
{ 
    if (class(data) != "phyDat") 
        stop("data must be of class phyDat")
    levels <- attr(data, "levels")
    if(class(tree)=="multiPhylo"){ 
        TL = attr(tree,"TipLabel")
        if (!is.null(TL)){
            data <- subset(data, TL)
            nTips <- length(TL) 
            weight <- attr(data, "weight")   
            nr <- attr(data, "nr")
            m <- nr*(2L*nTips - 1L)
        } 
    }
    data <- prepareDataFitch(data) 
    d = attributes(data)
    data <- as.integer(data)
    attributes(data) <- d
    if(class(tree)=="phylo") return(fit.fitch(tree, data, site))
    {
        if(is.null(attr(tree,"TipLabel"))){
            tree = unclass(tree)
            return(sapply(tree, fit.fitch, data, site))
        }    
        else{
            tree = unclass(tree)
#            tree = lapply(tree, reorder, "postorder")
            site = ifelse(site == "pscore", 1L, 0L) 
            on.exit(.C(fitch_free))
            .C(fitch_init, as.integer(data), as.integer(nTips*nr), as.integer(m), as.double(weight), as.integer(nr)) 
            return(sapply(tree, fast.fitch, nr, site)) 
        }       
    }
}


fit.fitch <- function (tree, data, returnData = c("pscore", "site", "data")) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    returnData <- match.arg(returnData)
    nr = attr(data, "nr")
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    weight = attr(data, "weight")
    m = max(tree$edge) 
    q = length(tree$tip)
    result <- .Call(FITCH, data[, tree$tip.label], as.integer(nr), as.integer(node), as.integer(edge), as.integer(length(edge)), as.double(weight), as.integer(m), as.integer(q))
    if (returnData == "site") return(result[[2]])
    pscore <- result[[1]]
    res = pscore
    if (returnData == "data") 
        res <- list(pscore = pscore, dat = result[[3]], site = result[[2]])
    res
}   


# NNI
fnodesNew2 <- function (EDGE, nTips, nr) 
{
    node <- EDGE[, 1]
    edge <- EDGE[, 2]
    n = length(node)
    m= as.integer(max(EDGE)+1L)
    m2 = 2L*n
    root0 <- as.integer(node[n]) 
    .Call(FNALL_NNI, as.integer(nr), node, edge, as.integer(n), as.integer(m), as.integer(m2), as.integer(root0))
}   


# SPR und bab kompakter
fnodesNew5 <- function (EDGE, nTips, nr) 
{
    node <- EDGE[, 1]
    edge <- EDGE[, 2]
    n = length(node)
    m= as.integer(max(EDGE)+1L)
    m2 = 2L*n
    root0 <- as.integer(node[n]) 
    .Call(FNALL5, as.integer(nr), node, edge, as.integer(n), as.integer(m), as.integer(m2), as.integer(root0))
}   


random.addition <- function(data, method="fitch") 
{
    label <- names(data)
    nTips <- as.integer(length(label))
    remaining <- as.integer(sample(nTips))  
    tree <- structure(list(edge = structure(c(rep(nTips+1L, 3), remaining[1:3]), .Dim = c(3L, 2L)), 
    tip.label = label, Nnode = 1L), .Names = c("edge", "tip.label", "Nnode"), class = "phylo", order = "postorder")
    remaining <- remaining[-c(1:3)]
    
    if(nTips==3L) return(tree)
 
    nr <- attr(data, "nr")
    storage.mode(nr) <- "integer"
    n <- length(data) #- 1L
     
    data <- subset(data,,order(attr(data, "weight"), decreasing=TRUE))   
    data <- prepareDataFitch(data) 
    weight <- attr(data, "weight")

    m = nr*(2L*nTips - 2L)

    on.exit(.C(fitch_free))
    .C(fitch_init, as.integer(data), as.integer(nTips*nr), as.integer(m), as.double(weight), as.integer(nr))
    
    storage.mode(weight) <- "double"

    for (i in remaining) {               
        edge = tree$edge[,2]   
        score = fnodesNew5(tree$edge, nTips, nr)[edge]      
        score <- .Call(FITCHTRIP3, as.integer(i), as.integer(nr), as.integer(edge), as.double(score), as.double(Inf))    
        res = min(score) 
        nt = which.min(score)  
        tree = addOne(tree, i, nt) 
        }
    attr(tree, "pscore") = res
    tree 
}

 
fast.fitch <- function (tree,  nr, ps = TRUE) 
{
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = max(tree$edge) 
    .Call(FITCH345, as.integer(nr), as.integer(node), as.integer(edge), as.integer(length(edge)), as.integer(m), as.integer(ps))
}


fitch.spr <- function(tree, data){
  nTips = as.integer(length(tree$tip))
  nr = attr(data, "nr")
  minp = fast.fitch(tree, nr, TRUE)
  
  for(i in 1:nTips){
    treetmp = dropTip(tree, i)   
    edge = treetmp$edge[,2] 
    score = fnodesNew5(treetmp$edge, nTips, nr)[edge]   
    score <- .Call(FITCHTRIP3, as.integer(i), as.integer(nr), as.integer(edge),  as.double(score), as.double(minp))  
    
    if(min(score)<minp){
      nt = which.min(score)
      tree = addOne(treetmp, i, nt) 
      minp=min(score)
      #   print(paste("new",minp))
    }
  }
  m=max(tree$edge)
  root <- getRoot(tree) 
  ch = allChildren(tree)
  for(i in (nTips+1L):m){
      if(i!=root){
          tmp = dropNode(tree, i, all.ch=ch)
          if(!is.null(tmp)){
          edge = tmp[[1]]$edge[,2]                          
          blub = fast.fitch(tmp[[2]], nr, TRUE)
          score = fnodesNew5(tmp[[1]]$edge, nTips, nr)[edge] + blub
          score <- .Call(FITCHTRIP3, as.integer(i), as.integer(nr), as.integer(edge), as.double(score), as.double(minp))    
          if(min(score)<minp){
              nt = which.min(score)
              tree = addOneTree(tmp[[1]], tmp[[2]], nt, tmp[[3]])
              minp <- min(score)
              ch = allChildren(tree)
            }
        }
      }
  }
  tree
}

# raus 
fitch.spr2 <- function(tree, data){
    nTips = as.integer(length(tree$tip))
    nr = attr(data, "nr")
    minp = fast.fitch(tree, nr, TRUE)
    
    changeIndex <- function(x, i, j){
        x$edge[x$edge == i] = 0L
        x$edge[x$edge == j] = i
        x$edge[x$edge == 0L] = j
        x
    }
    
    
    for(i in 1:nTips){
        treetmp = dropTip(tree, i)   
        edge = treetmp$edge[,2] 
        score = fnodesNew5(treetmp$edge, nTips, nr)[edge]   
        score <- .Call(FITCHTRIP3, as.integer(i), as.integer(nr), as.integer(edge),  as.double(score), as.double(minp))  
        
        if(min(score)<minp){
            nt = which.min(score)
            tree = addOne(treetmp, i, nt) 
            minp=min(score)
            #            print(paste("new",minp))
        }
    }
    m=max(tree$edge)
    
    root <- getRoot(tree) 
    for(i in (nTips+1L):m){
        if(i!=root){
            tmp = dropNode(tree, i)
            print(i)
            if(!is.null(tmp)){
            edge = tmp[[1]]$edge[,2]                          
            blub = fast.fitch(tmp[[2]], nr, TRUE)
            score = fnodesNew5(tmp[[1]]$edge, nTips, nr)[edge] + blub
            score <- .Call(FITCHTRIP3, as.integer(i), as.integer(nr), as.integer(edge), as.double(score), as.double(minp))    
            if(min(score)<minp){
                nt = which.min(score)
                tree = addOneTree(tmp[[1]], tmp[[2]], nt, tmp[[3]])
                minp <- min(score)
            }
            }
#            browser()
#            j = Ancestors(tree, i, "parent")
#            tree2 = reroot(tree, node=i)
#            tree2 = unroot(tree2)
#            tree2 = reorder(tree2, "postorder")
#         if(j == (nTips+1L)) tree2 = changeIndex(tree2, as.integer(nTips+1L), as.integer(i))
#            tmp = dropNode(tree2, j)
#            if(!is.null(tmp)){
#            edge = tmp[[1]]$edge[,2]                          
#            blub = fast.fitch(tmp[[2]], nr, TRUE)
#            score = fnodesNew5(tmp[[1]]$edge, nTips, nr)[edge] + blub
#            score <- .Call("FITCHTRIP3", as.integer(j), as.integer(nr), as.integer(edge), as.double(score), as.double(minp), PACKAGE="phangorn")    
#            if(min(score)<minp){
#                nt = which.min(score)
#                tree = addOneTree(tmp[[1]], tmp[[2]], nt, tmp[[3]])
#                minp <- min(score)
#            }
#            }
        }
    }
    tree
}
           
indexNNI2 <- function(tree){
    parent = tree$edge[, 1]
    child = tree$edge[, 2]
 
    ind = which(child %in% parent)
    Nnode = tree$Nnode
    edgeMatrix = matrix(0L, 6, length(ind))

    pvector <- numeric(max(parent))
    pvector[child] <- parent
    cvector <- allChildren(tree)  

    k=0
    for(i in ind){        
            p1 = parent[i]          
            p2 = child[i]
            e34 = cvector[[p2]]
            ind1 = cvector[[p1]]
            e12 = ind1[ind1 != p2]
            if(pvector[p1]) edgeMatrix[, k+1] = c(p1,e12, e34, p2, 1L) 
            else edgeMatrix[, k+1] = c(e12, e34, p2, 0L)
            k=k+1
    } 
    cbind(edgeMatrix[c(1,3,2,4,5,6),], edgeMatrix[c(1,4,2,3,5,6),])
}
       
# nr statt data uebergeben, fitchQuartet ohne weight
# weniger Speicher 2 Zeilen weinger 
fitch.nni <- function (tree, data, ...) 
{
    nTips = as.integer(length(tree$tip)) # auskommentieren?
    INDEX <- indexNNI2(tree)    
    nr = attr(data, "nr")
    weight <- attr(data, "weight")
    p0 <- fast.fitch(tree, nr)
    m <- dim(INDEX)[2]    
    tmp = fnodesNew2(tree$edge, nTips, nr)
    pscore <- .C(fitchQuartet, as.integer(INDEX), as.integer(m), as.integer(nr), as.double(tmp[[1]]), as.double(tmp[[2]]), as.double(weight), double(m))[[7]]    
    swap <- 0
    candidates <- pscore < p0
    while (any(candidates)) {
        ind = which.min(pscore)
        pscore[ind] = Inf
        tree2 <- changeEdge(tree, INDEX[c(2,3), ind])        
        test <- fast.fitch(tree2, nr)
        if (test >= p0) 
            candidates[ind] = FALSE
        if (test < p0) {
            p0 <- test
            swap = swap + 1
            tree <- tree2
            indi <- which(INDEX[5,] %in% INDEX[1:5, ind])
            candidates[indi] <- FALSE
            pscore[indi] <- Inf
        }
    }
    list(tree = tree, pscore = p0, swap = swap)
}


optim.fitch <- function(tree, data, trace=1, rearrangements = "SPR", ...) {
    if(class(tree)!="phylo") stop("tree must be of class phylo") 
    if(!is.binary.tree(tree)){
        tree <- multi2di(tree)
        attr(tree, "order") <- NULL  
    }
    if(is.rooted(tree)){
        tree <- unroot(tree)
        attr(tree, "order") <- NULL
    }
    if(is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") tree <- reorder(tree, "postorder")  
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")

    rt = FALSE
    nTips = as.integer(length(tree$tip))

    nr = attr(data, "nr")    
    pis <- parsinfo(data)
    p0 <- sum(attr(data, "weight")[pis[, 1]] * pis[, 2])
    if (length(pis) > 0) 
        data <- getRows(data, c(1:nr)[-pis[, 1]], TRUE)    
    
    nr = attr(data, "nr")
   
    data <- subset(data,tree$tip,order(attr(data, "weight"), decreasing=TRUE))   
    dat <- prepareDataFitch(data) 
    weight <- attr(data, "weight")

    m = nr*(2L*nTips - 2L)
    on.exit(.C(fitch_free))
    .C(fitch_init, as.integer(dat), as.integer(nTips*nr), as.integer(m), as.double(weight), as.integer(nr))

    tree$edge.length=NULL
    swap = 0
    iter = TRUE
    pscore <- fast.fitch(tree, nr)  
    while (iter) {
        res <- fitch.nni(tree, dat, ...)
        tree <- res$tree
        if(trace>1)cat("optimize topology: ", pscore + p0, "-->", res$pscore + p0, 
            "\n")
        pscore = res$pscore
        swap = swap + res$swap
        if (res$swap == 0){
            if(rearrangements=="SPR"){
                tree <- fitch.spr(tree, dat)             
                psc <- fast.fitch(tree, nr)
                if(trace>1)cat("optimize topology (SPR): ", pscore + p0 , "-->", psc + p0, "\n")
                if(pscore < psc+1e-6) iter=FALSE
                pscore <- psc
            } 
            else iter = FALSE
        }
    }
    if(trace>0)cat("Final p-score",pscore + p0,"after ",swap, "nni operations \n") 
    if(rt)tree <- ptree(tree, data)  
    attr(tree, "pscore") = pscore + p0
    tree
}

# branch and bound
getOrder <- function (x) 
{
    label = names(x)
    dm = as.matrix(dist.hamming(x, FALSE))
    ind = as.vector(which(dm == max(dm), arr.ind = TRUE)[1, ])
    nTips = as.integer(length(label))
    added = ind
    remaining <- c(1:nTips)[-ind]

    tree <- structure(list(edge = structure(c(rep(nTips+1L, 3), c(ind, 0L)), .Dim = c(3L, 2L)), tip.label = label, Nnode = 1L), .Names = c("edge", "tip.label", "Nnode"), class = "phylo", order = "postorder")      

    l = length(remaining)
    res = numeric(l)

    nr <- attr(x, "nr")
    storage.mode(nr) <- "integer"
    n <- length(x) #- 1L
      
    data <- prepareDataFitch(x) 
    weight <- attr(data, "weight")
    storage.mode(weight) <- "double"

    m = nr*(2L*nTips - 2L)

    on.exit(.C(fitch_free))
    .C(fitch_init, as.integer(data), as.integer(nTips*nr), as.integer(m), as.double(weight), as.integer(nr))

    for(i in 1:length(remaining)){
        tree$edge[3,2]= remaining[i]     
        res[i] = fast.fitch(tree, nr) 
    }
    tmp = which.max(res)
    added = c(added, remaining[tmp])
    remaining <- remaining[-tmp]
    tree$edge[,2]= added

    for (i in 4:(nTips - 1L)) {
        edge = tree$edge[,2]                 
        score0 = fnodesNew5(tree$edge, nTips, nr)[edge]        
        
        l = length(remaining)
        res = numeric(l)
        nt = numeric(l)
        k = length(added)+1L
        for(j in 1:l){
            score <- .Call(FITCHTRIP3, as.integer(remaining[j]), as.integer(nr), as.integer(edge), as.double(score0), as.double(Inf))   
            
#            score = score0[edge] + psc
            res[j] = min(score) 
            nt[j] = which.min(score)
        }
        tmp = which.max(res)
        added = c(added, remaining[tmp])        
        tree = addOne(tree, remaining[tmp], nt[tmp])
        remaining <- remaining[-tmp]  
    }
    added = c(added, remaining) 
    added 
}


bab <- function (data, tree = NULL, trace = 1, ...) 
{
    o = order(attr(data, "weight"), decreasing = TRUE)
    data = subset(data, , o)
    nr <- attr(data, "nr")
    pis <- parsinfo(data)
    p0 <- sum(attr(data, "weight")[pis[, 1]] * pis[, 2])
    if (length(pis) > 0) 
        data <- getRows(data, c(1:nr)[-pis[, 1]], TRUE)
    tree <- pratchet(data, start = tree, trace = trace - 1, ...)
    data <- subset(data, tree$tip.label) 
    nr <- as.integer(attr(data, "nr"))
    inord <- getOrder(data)
    lb <- lowerBound(data)
    nTips <- m <- length(data)
    
    nr <- as.integer(attr(data, "nr"))
    TMP <- matrix(0, m, nr)
    for (i in 4:m) {
        TMP[i, ] = lowerBound(subset(data, inord[1:i]))
    }

    weight <- as.double(attr(data, "weight"))
    data <- prepareDataFitch(data)
    m = nr*(2L*nTips - 2L)
    on.exit(.C(fitch_free))
    .C(fitch_init, as.integer(data), as.integer(nTips*nr), as.integer(m), as.double(weight), as.integer(nr))
    mmsAmb = 0
    mmsAmb = TMP %*% weight  
    mmsAmb = mmsAmb[nTips] - mmsAmb
    mms0 = 0 
    mms0 = mms0 + mmsAmb

    minPars = mms0[1]
    kPars = 0

    if (trace) 
        print(paste("lower bound:", p0 + mms0[1]))
    bound <- fast.fitch(tree, nr)
    if (trace) 
        print(paste("upper bound:", bound + p0))

    startTree <- structure(list(edge = structure(c(rep(nTips+1L, 3), as.integer(inord)[1:3]), .Dim = c(3L, 2L)), 
        tip.label = tree$tip.label, Nnode = 1L), .Names = c("edge", "tip.label", "Nnode"), class = "phylo", order = "postorder")

    trees <- vector("list", nTips)
    trees[[3]] <- list(startTree$edge)
    for(i in 4:nTips) trees[[i]] <- vector("list", (2L*i) - 5L) # new

# index M[i] is neues node fuer edge i+1
# index L[i] is length(node) tree mit i+1 
    L = as.integer( 2L*(1L:nTips) -3L ) 
    M = as.integer( 1L:nTips + nTips - 1L )    

    PSC <- matrix(c(3,1,0), 1, 3)
    PSC[1,3] <- fast.fitch(startTree, nr)

    k = 4L
    Nnode = 1L
    npsc = 1

    result <- list() 
    while (npsc > 0) {
        a = PSC[npsc,1]
        b = PSC[npsc,2]
        PSC = PSC[-npsc,, drop=FALSE]  

        tmpTree <- trees[[a]][[b]]
        edge = tmpTree[,2]  
        score = fnodesNew5(tmpTree, nTips, nr)[edge] + mms0[a+1L] 
        score <- .Call(FITCHTRIP3, as.integer(inord[a+1L]), as.integer(nr), as.integer(edge), as.double(score), as.double(bound))    
                   
        ms = min(score)
        if(ms<=bound){
            if((a+1L)<nTips){
                ind = (1:L[a])[score<=bound]
                for(i in 1:length(ind))trees[[a+1]][[i]] <- .Call(AddOne, tmpTree, as.integer(inord[a+1L]), as.integer(ind[i]), as.integer(L[a]), as.integer(M[a])) 
                l = length(ind)
                os = order(score[ind], decreasing=TRUE)                 
                PSC = rbind(PSC, cbind(rep(a+1, l), os, score[ind][os] ))
            }
            else{
                ind = which(score==ms) 
                tmp <- vector("list", length(ind)) 
                for(i in 1:length(ind))tmp[[i]] <- .Call(AddOne, tmpTree, as.integer(inord[a+1L]), as.integer(ind[i]), as.integer(L[a]), as.integer(M[a]))

                if(ms < bound){
                     bound = ms
                     if(trace)cat("upper bound:", bound, "\n") 
                     result = tmp    
                     PSC = PSC[PSC[,3]<(bound+1e-8),]  

                }
                else result = c(result, tmp)  
            }
        }    
        npsc = nrow(PSC)
    }
    for(i in 1:length(result)){
        result[[i]] = structure(list(edge = result[[i]], Nnode = nTips-2L), .Names = c("edge", "Nnode"), class = "phylo", order = "postorder")
    }
    attr(result, "TipLabel") = tree$tip.label
    class(result) <- "multiPhylo"
    return(result)
}

