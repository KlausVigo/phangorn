###########################################################################################


getClans = function (tree) 
{
	if (is.rooted(tree)) 
        tree = unroot(tree)
    bp = bip(tree)
    nTips = length(tree$tip)
    root = nTips + 1
    bp[root] = NULL
    X = matrix(0, length(bp) - nTips, nTips)
    k = 1
    nl = NULL
    if (!is.null(tree$node.label)) {
        nl = c(rep("-1", nTips), rep("-1", nTips), tree$node.label[-1], 
            tree$node.label[-1])
    }
    if(root<=length(bp)){
        for (i in root:length(bp)) {
           X[k, bp[[i]]] = 1
           k = k + 1
        }
    }
    res <- rbind(diag(nTips), 1 - diag(nTips), X, 1 - X)
    colnames(res) <- tree$tip
    if (!is.null(nl)) 
        rownames(res) = nl
    res
}


getSlices <- function(tree){
    nTips = length(tree$tip)
    clans = getClans(tree)
    m = dim(clans)[1]
    X = tcrossprod(clans)
    z = rowSums(clans)
    Z1 = matrix(z,m,m)
    Z2 = t(Z1)
    Z = matrix(0,m,m)
    Z[Z1<=Z2] = Z1[Z1<=Z2]
    Z[Z2<Z1] = Z2[Z2<Z1]

    diag(X)=0
    X[upper.tri(X)] = 0
    X[X==1] = 0
    X[X==Z] = 0 
    index = which(X>0,arr.ind=TRUE)
    l = dim(index)[1]
    nSlices = 2 * nTips^2 -10 * nTips + 12
    result = matrix(0, nSlices, nTips)
    strClan = do.call("paste", c(as.data.frame(clans), sep = ""))
    k=1
    for(i in 1:l){
        tmp1 = as.numeric((clans[index[i,1],] + clans[index[i,2],])==2)
        tmp = paste(tmp1,sep="",collapse="")
        if(is.na(match(tmp,strClan))){ 
            result[k,] = tmp1
            k=k+1
        }
    }
    if(k<nSlices) result = result[1:(k-1),] 
    colnames(result) <- tree$tip
    result   
}


getClips = function (tree, all = TRUE) 
{
    if (any(is.na(tree$edge.length))) 
        return(NULL)
    dm = cophenetic(tree)
    tips = tree$tip
    nTips = length(tips)
    res = numeric(nTips)
    result = NULL
    for (i in 1:nTips) {
        ord = order(dm[i, ])
        for (j in 2:(nTips - 1)) {
            ind = ord[1:j]
            
            if(i>min(ind) ) break()
            within = max(dm[ind, ind])
            between = min(dm[ind, -ind])
            if (within < between) {
                res = numeric(nTips)
                res[ind] = 1L
                result = rbind(result, res)
            }
        }
    }
    dimnames(result) = list(NULL, tips)  
    if (all) 
        return(result)
    ind = which.max(rowSums(result))
    result[ind, ]
}


shannon <- function (x, norm=TRUE) 
{
    p = x/sum(x)
    ShD = -sum(p * log10(p))
    if(norm){
        if (sum(x) == 1) return(0)
        ShD = ShD/log10(sum(x))
    }
    ShD
}


shannon2 <- function (x, norm=TRUE) 
{
    p = x/sum(x)
    ShD = -sum(p * log(p))
    if(norm){
        if (sum(x) == 1) return(0)
        ShD = ShD/log(sum(x))
    }
    ShD
}


getE = function (tree, x, clans = NULL, norm = TRUE) 
{
    if (is.rooted(tree)) 
        tree = unroot(tree)
    if (is.null(clans)) 
        clans = getClans(tree)
    labels = tree$tip.label
    x = x[labels]
    result = rep(NA, 12)
    names(result) = c("E* tree", "# natives", "# intruder", "# unknown", 
        "E* clan", "# intruder", "# unknown", "E* slice", "# intruder", 
        "# unknown", "bs 1", "bs 2")
    result[2] = sum(x == 1)
    result[3] = sum(x == 2)
    result[4] = sum(x == 3)
    if (result[2] == 0 || result[3] == 0) {
        if (result[2] > 1) 
            return(list(result, labels))
        else return(list(result, integer(0)))
    }
    LHG = E_Intruder(clans, x)
    d = dim(LHG)[1]
    if (d == 1) {
        result[1] = 0
        if (!is.null(tree$node.label)) 
            result[11] = as.numeric(rownames(LHG))
        return(list(result, labels[LHG == 0]))
    }
    intr = drop(LHG %*% as.numeric(x == 2))
    result[1] = shannon2(intr, norm = norm)
    o <- order(intr, decreasing = TRUE)
    if (!is.null(tree$node.label)) 
        result[11:12] = as.numeric(rownames(LHG)[o[c(1, 2)]])
    ind = which(LHG[o[1], ] == 1)
    result[6] = sum(x[-ind] == 2)
    result[7] = sum(x[-ind] == 3)


    if (length(x[-ind]) < 4) 
        return(list(result, NULL))
    result[5] = shannon2(intr[-o[1]], norm = norm)
    ind2 = c(which(LHG[o[1], ] == 1), which(LHG[o[2], ] == 1))
  
    spl = structure(list(which(colSums(LHG)==0)), labels=labels, weights=1)
    class(spl)="splits"

    if (d == 2) {
         return(list(result, spl))
    } 
    result[9] = sum(x[-ind2] == 2)
    result[10] = sum(x[-ind2] == 3)
    if (length(x[-ind2]) < 4){ 
         return(list(result, spl))
    } 
    result[8] = shannon2(intr[-c(o[1], o[2])], norm = norm)
    return(list(result, spl))
}


E_Intruder <- function (clans, x) 
{
    cp = drop(clans %*% as.numeric(x == 1))
    ci = drop(clans %*% as.numeric(x == 2))
    homo = which(cp == 0 & ci > 0)
    l = length(homo)
    if (l > 0) {
        HG = clans[homo, , drop = FALSE]
        lhg = rep(TRUE, l)
        rsh = rowSums(HG)
        Z = tcrossprod(HG)>0
        Z = Z * rsh
        zmax = apply(Z,2,max)
        lhg = !(zmax > rsh)  
        LHG = HG[lhg, , drop = FALSE]
        return(LHG)
    }
    return(NULL)
}


E_Intruder_2 <- function (clans, x, native=NULL) 
{     
    contr = attr(x, "contr")
    d = dim(contr)
    if(d[1]>d[2])contr[(d[2]+1):d[1],]=0
    cp = clans %*% contr[as.numeric(x),]
    
	homo = which(rowSums(cp > 0) == 1)
		
    l = length(homo)
    if (l > 0) {
        HG = clans[homo, , drop = FALSE]
        lhg = rep(TRUE, l)
        rsh = rowSums(HG)
        Z = tcrossprod(HG)>0
        Z = Z * rsh
        zmax = apply(Z,2,max)
        lhg = !(zmax > rsh)  
        LHG = HG[lhg, , drop = FALSE]
        return(LHG)
    }
    return(NULL)
}


getDiv <- function(tree, x, native=NULL){
    clans = getClans(tree)
    labels = tree$tip.label
    x = subset(x, labels)
    LHG = E_Intruder_2(clans, subset(x,,1))
    if(!is.null(native)){
	    ll = match(native, attr(x, "allLevels"))
	    ind = (as.numeric(x) %in% ll)
	    }    	    
	if(!is.null(native)){    
	    rs = rowSums(clans)
	    intr = clans %*% ind    
	    clans = clans[intr==0,]
	    d = which.max(rs[intr==0])
	    tree2 = drop.tip(tree, tip=labels[which(clans[d, ]==1)])
    } 
    else tree2=NULL
	list(c(shannon(rowSums(LHG)),      
    summary(factor(attr(x, "allLevels"))[as.numeric(subset(x,,1))]), parsimony(tree, x)), tree2 )     
}


getDiversity <- function (tree, x, norm = TRUE, var.names = NULL, labels="new") 
{
    k = 1
    if(inherits(tree,"multiPhylo")) 
        k = length(tree)
    l = attr(x, "nr")
    tmp = matrix(0, k * l, 12)

    tnam = 1
    if (inherits(tree,"multiPhylo")) {
        tnam = names(tree)
        if (is.null(tnam)) 
            tnam = 1:length(tree)
    }
    if(is.null(var.names)) var.names = 1:l
    PM = data.frame("t1", "a", stringsAsFactors = FALSE)
    colnames(PM) = c("Tree", "Var")
    PM = PM[FALSE,] 
    PM[1 :(k*l), ] = NA 
    perfect = names(x)
    L = vector("list",k*l)
    m = 1
    o = 1
    ok= 0
    for (i in 1:k) {
        if (inherits(tree,"multiPhylo")) 
            tmptree = tree[[i]]
        else tmptree = tree
        if (is.rooted(tmptree)) 
            tmptree = unroot(tmptree)
        clans = getClans(tmptree) 
        for (j in 1:l) {          
            TMP = getE(tmptree, getRows(x, j), clans, norm = norm)
            tmp[m, ] = TMP[[1]]
            L[[m]] = TMP[[2]] # if class =splits else NULL
            PM[m, 1] = tnam[i]
            PM[m, 2] = var.names[j]
            m = m + 1   
        }
    }

    tnam = rep(tnam, each = l)
    dnam = var.names
    dnam = rep(dnam, k)
    pscore = as.numeric(sankoff(tree, x, site = "site"))
    res = data.frame(tnam, dnam, tmp, pscore)
    if(labels=="old")names(res) = c("tree", "variable", "E tree", "# natives", 
        "# intruder", "# unknown", "E clan", "# intruder", "# unknown", 
        "E slice", "# intruder", "# unknown", "bs 1", "bs 2", "p-score")
    else{
        names(res) = c("tree", "variable", "E clan", "# natives", 
            "# intruder", "# unknown", "E slice", "# intruder", "# unknown", 
            "E melange", "# intruder", "# unknown", "bs 1", "bs 2", "p-score")    
        warning("The variable names have changed")       
    }    
    attr(res, "Perfect") = L
    class(res) = c("clanistics", "data.frame")
    res
}


summary.clanistics <- function(object, ...){
    res <- matrix(FALSE, nrow(object), 5)
    res[,1] = object[,4]>0 & object[,"p-score"]==0 # "natives"
    res[,2] = object[,5]>0 & object[,"p-score"]==0 # "intruder"
    res[,3] = object[,"p-score"]==1
    res[,4] = ( (object[,"p-score"]==2) & (object[,7]==0) & (!is.na(object[,7])) ) | 
              ( (object[,"p-score"]==2) & (object[,4]==2) & (is.na(object[,7])) )  
    res[,5] = object[,"p-score"]>=2 & (object[,7]>0) & (!is.na(object[,7]))
    res[] = as.numeric(res)
    tmp = data.frame(factor(object[,"variable"]), res)	
    colnames(tmp) = c("Variable", "Natives_only", "Intruder_only", "Clan", "Slice", "Melange")
#        colnames(res) = c("Natives only", "Intruder only", "Clan", "Melange")
    class(tmp) <- c("summary.clanistics", "data.frame")
    tmp
    }
	

print.summary.clanistics <- function(x, ...){
    print(aggregate(x[,-1], list(Variable=x[,1]), sum), ...)
}


compareSplits <- function(res, nam1, nam2){
    wide <- reshape(res[, c("tree", "E tree", "variable")], v.names="E tree", idvar="tree", timevar="variable", direction="wide")
    wideI <- reshape(res[, c("tree", "# natives", "variable")], v.names="# natives", idvar="tree", timevar="variable", direction="wide")
    for(i in 2:dim(wide)[2])colnames(wide)[i] = strsplit(colnames(wide)[i],"E tree.")[[1]][2]
    for(i in 2:dim(wide)[2])colnames(wideI)[i] = strsplit(colnames(wideI)[i],"# natives.")[[1]][2]
	ntrees = wide[,1]
	splits = attr(res, "Perfect")
	dat = attr(attr(res, "Perfect"), "data")
    res = matrix(NA, length(ntrees), length(nam1)*length(nam2))
    for(m in 1:ntrees){
        k=1
        trnam=ntrees[m]
        for(i in nam1){
            E1 = wide[m, i]
            for(j in nam2){
                E2 = wide[m, j] 
                if(!is.na(E1) & !is.na(E2)){
                    if(E1 == E2){ # if(E1 == 0 & E2 == 0){
	                if( (wideI[m, i] >0) & (wideI[m, j]) >0){
                            ind1 = which(dat[,1]==trnam & dat[,2]==i)
                            sp1 = splits[[ind1]]                     
                            ind2 = which(dat[,1]==trnam & dat[,2]==j) 
                            sp2 = splits[[ind2]]
                            if(length(ind1)>0 & length(ind2)>0 )res[m, k] = drop(compatible3(sp1, sp2))
                        }
                    }
                }
            k=k+1 
            }
        }
    }    
    res
}


diversity <- function(tree, X){  
# from kknn
    contr.dummy <- function (n, contrasts = TRUE) 
    {
        if (length(n) <= 1) {
            if (is.numeric(n) && length(n) == 1 && n > 1) 
                levels <- 1:n
            else stop("contrasts are not defined for 0 degrees of freedom")
        }
        else levels <- n
        lenglev <- length(levels)
       cont <- array(0, c(lenglev, lenglev), list(levels, levels))
       cont[col(cont) == row(cont)] <- 1
       cont
    }


    l = dim(X)[2]
    m <- ifelse(inherits(tree,"multiPhylo"), length(tree), 1)
    
    contr = as.list(rep("contr.dummy", l))
    names(contr) = names(X)
    tmp = model.matrix(~.-1, X, contrast=contr)
    tmp1 <- phyDat.default(tmp, levels=c(1,0), compress = FALSE)
    attr(tmp1, "varnames")  = colnames(tmp)
    fd = sankoff(tree,tmp1,site = "site") 
    fd = matrix(fd, ncol=m) 

    if(m>1){
         if(is.null(names(tree))) tnames <- paste("tree", 1:m, sep=".")
         else tnames <- names(tree)  
    }
    else tnames = "tree"
    dimnames(fd) = list(colnames(tmp), tnames)
    res = stack(data.frame(fd))

    if(m>1)nt = rep(sapply(tree, function(x)length(x$tip)), each=dim(fd)[1])    
    else nt = rep(length(tree$tip), each=dim(fd)[1]) 
    if(m>1)res2 = as.vector(sapply(tree, function(x,y)colSums(y[x$tip,,drop=FALSE]) , y=tmp))
    else res2 = colSums(tmp[tree$tip,,drop=FALSE])
    result <- data.frame(tree = res[,2], variable=rep(colnames(tmp),m), pscore=res[,1], ntips=nt, natives=res2)
    result
}
