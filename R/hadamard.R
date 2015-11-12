dec2Bin = function (x) 
{
    res = NULL
    i = 1L
    while (x > 0) {
        if (x%%2L) 
            res = c(res, i)
        x = x%/%2L
        i = i + 1L
    }
    res
}


# returns binary (0, 1) vector of length k
dec2bin <- function (x, k=ceiling(log2(x))) 
{
    i = 1L
    res = integer(k)
    while (x > 0) {
        if (x%%2L) 
            res[i] = 1L
        x = x%/%2L
        i = i + 1L
    }
    res
}

# double factorial: log version
"ldfactorial" <- function(x){
    x = (x+1)/2
    res = lgamma(2*x)-(lgamma(x)+(x-1)*log(2))
    res
}

# double factorial
"dfactorial" <- function(x){exp(ldfactorial(x))}


#
# Hadamard Conjugation
#

hadamard <- function(x){
    res=1
    while(x>0){
        res=rbind(cbind(res,res),cbind(res,-res))
        x=x-1
    }
    res
}


fhm <- function(v){
    n = length(v)
    n = log2(n)
    res = .C("C_fhm", v = as.double(v), n = as.integer(n))$v # 
    res
}


seq2split = function(s){
    n=length(s)
    res= fhm(log(fhm(s)))/n
    res
}


split2seq = function(q){
    n=length(q)
    res= fhm(exp(fhm(q)))/n
    res
}


distanceHadamard <- function (dm, eps = 0.001) 
{
    if (inherits(dm,"dist")) {
        n <- attr(dm, "Size")
        Labels = attr(dm, "Labels")
    }
    if (inherits(dm,"matrix")) {
        n <- dim(dm)[1]
        Labels <- colnames(dm)
        dm <- dm[lower.tri(dm)]
    }
    ns <- 2^(n - 1)
    if (n > 23) 
        stop("Hadamard conjugation works only efficient for n < 24")
    result <- .Call("dist2spectra", dm, as.integer(n), as.integer(ns), 
                    PACKAGE = "phangorn")
    weights = -fhm(result)/2^(n - 2)    
    
    if(eps>0){
        weights = weights[-1]
        ind2 = which(weights>eps)
        n2 = length(ind2)
        splits = vector("list", n2)
        for(i in 1:n2)splits[[i]] = dec2Bin(ind2[i])
        attr(splits, "weights") = weights[ind2]
        attr(splits, "labels") = Labels
        attr(splits, 'dm') = dm
        class(splits)='splits'
        return(splits)      
    }  
    res <- data.frame(distance = result, edges = weights, index = 0:(ns - 1))
    attr(res, "Labels") <- Labels
    res
}


h4st = function(obj, levels=c('a','c','g','t')){
    if (is.matrix(obj)) 
        obj = as.data.frame(t(obj))
    if (inherits(obj,"phyDat")) 
        obj = as.data.frame(t(as.character(obj)))    
    #    if(is.matrix(obj)) obj = as.data.frame(t(obj))
    #    DNA = as.data.frame(obj)
    #    DNA = t(as.character(obj))
    
    n = dim(obj)[1]
    p = dim(obj)[2]
    
    if(p>11) stop("4-state Hadamard conjugation works only efficient for n < 12")
    
    DNAX = matrix(0,n,p)
    DNAY = matrix(0,n,p)
    
    DNAX[obj==levels[1]]=0
    DNAX[obj==levels[2]]=1
    DNAX[obj==levels[3]]=1
    DNAX[obj==levels[4]]=0
    
    DNAY[obj==levels[1]]=0
    DNAY[obj==levels[2]]=1
    DNAY[obj==levels[3]]=0
    DNAY[obj==levels[4]]=1
    
    DNAY = DNAY - DNAY[,p]
    DNAX = DNAX - DNAX[,p]
    
    DNAY = abs(DNAY[,-p])
    DNAX = abs(DNAX[,-p])
    dy = DNAY %*% (2^(0:(p-2))) 
    dx = DNAX %*% (2^(0:(p-2))) 
    
    INDEX =  dx + 2^(p-1) * dy
    blub = table(INDEX)
    index = as.numeric(rownames(blub)) + 1
    sv = numeric(4^(p-1))
    sv[index] = blub
    qv = matrix(seq2split(sv),2^(p-1),2^(p-1))
    sv = matrix(sv,2^(p-1),2^(p-1))
    #    q = cbind(transversion = qv[-1,1], transition.1 = diag(qv)[-1], transition.2 = qv[1,-1])
    transversion <- transition.1 <- transition.2 <- allSplits(p, colnames(obj)) 
    attr(transversion,"weights") = qv[-1,1]
    attr(transition.1,"weights") = diag(qv)[-1]
    attr(transition.2,"weights") = qv[1,-1]
    #    result = list(q = q, qv = qv, sv=sv, n=sum(sv), names=names(obj))
    result = list(transversion = transversion, transition.1=transition.1, transition.2 = transition.2, 
                  qv = qv, sv=sv, n=sum(sv), names=names(obj))
    result
}


h2st <- function (obj, eps=0.001) 
{
    if (!inherits(obj,"phyDat")) stop("Error") 
    if (attr(obj,"nc") != 2)stop("Error")
    nr = attr(obj, "nr") #n
    p = length(obj) #p
    weight = attr(obj, "weight")
    if (p > 23) 
        stop("Hadamard conjugation works only efficient for n < 24")
    DNAX = matrix(0, nr, p-1)
    for(i in 1:(p-1)) DNAX[,i] = obj[[i]]-1
    DNAX[obj[[p]]==2,] = 1 - DNAX[obj[[p]]==2,]
    
    index = DNAX %*% (2^(0:(p - 2))) + 1
    sv = numeric(2^(p - 1))
    for(i in 1:nr)sv[index[i]] = sv[index[i]]+ weight[i]
    qv = seq2split(sv)
    
    if(eps>0){
        qv = qv[-1]
        ind2 = which(qv>eps)
        indT= c(2L^(0:(p-2)), 2L^(p-1)-1) 
        ind2 = union(ind2, indT)
        n2 = length(ind2)
        splits = vector("list", n2)
        for(i in 1:n2)splits[[i]] = dec2Bin(ind2[i])
        attr(splits, "weights") = qv[ind2]
        attr(splits, "labels") = names(obj)
        class(splits)='splits'
        return(splits)    
    }
    result = data.frame(edges = qv, splits = sv, index = 0:(2^(p - 
                                                                   1) - 1))
    attr(result, "Labels") = names(obj)
    result
}
