#
# dist
#
dist.hamming <- function (x, ratio = TRUE, exclude = "none") 
{
    if (!inherits(x,"phyDat")) 
        stop("x has to be element of class phyDat")
    l = length(x)

    contrast <- attr(x, "contrast")
    nc <- as.integer(attr(x, "nc"))
    con = rowSums(contrast > 0) < 2
    if (exclude == "all") {
        index = con[x[[1]]]
        for (i in 2:l) index = index & con[x[[i]]]
        index = which(index)
        x = subset(x, , index)
    }
    weight <- attr(x, "weight")  
    d = numeric((l * (l - 1))/2)

    if(exclude == "pairwise"){
        k=1
        W <- numeric(l*(l-1)/2)
        for (i in 1:(l - 1)) {
            tmp = con[x[[i]]] 
            for (j in (i + 1):l) {
                W[k] = sum(weight[tmp & con[ x[[j]] ] ])
                k = k + 1
            }
        }  
             
    } 

    if(nc > 31){
#        contrast <- attr(x, "contrast")
        k = 1
        for (i in 1:(l - 1)) {
            X = contrast[x[[i]], , drop = FALSE]
            for (j in (i + 1):l) {
                d[k] = sum(weight * (rowSums(X * contrast[x[[j]], , drop = FALSE]) == 0))
                k = k + 1
            }
        }
    } # end if  
    else{
        nr <- attr(x, "nr")
        if(exclude == "pairwise")ind <- which(con[unlist(x)]==FALSE)  
        x <- prepareDataFitch(x) 
        if(exclude == "pairwise")x[ind] <- as.integer(2L^nc -1L) 
        res <- .C("distHamming", as.integer(x), as.double(weight), as.integer(nr), as.integer(l), as.double(d), PACKAGE = "phangorn")
        d <- res[[5]]
    }     

    if (ratio){
        if(exclude == "pairwise") d = d/W
        else d = d/sum(weight)
    }
    attr(d, "Size") <- l
    if (is.list(x)) 
        attr(d, "Labels") <- names(x)
    else attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- "hamming"
    class(d) <- "dist"
    return(d)
}



dist.ml <- function (x, model = "JC69", exclude = "none", bf = NULL, Q = NULL, k=1L, shape=1, ...) 
{
    if (!inherits(x,"phyDat")) 
        stop("x has to be element of class phyDat")
    l = length(x)
    d = numeric((l * (l - 1))/2)
    v = numeric((l * (l - 1))/2)
    contrast <- attr(x, "contrast")
    con = rowSums(contrast > 0) < 2
    if (exclude == "all") {
        index = con[x[[1]]]
        for (i in 2:l) index = index & con[x[[i]]]
        index = which(index)
        x = subset(x, , index)
    }
    nc <- as.integer(attr(x, "nc"))
    nr <- as.integer(attr(x, "nr"))
    #    model <- match.arg(model, c("JC69", "WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
    model <- match.arg(model, c("JC69", "F81", .aamodels))
    #    if (!is.na(match(model, c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24")))) 
    if (!is.na(match(model, .aamodels))) 
        getModelAA(model, bf = is.null(bf), Q = is.null(Q))
    if(is.null(bf) && model=="F81") bf <- baseFreq(x)
    if (is.null(bf)) 
        bf <- rep(1/nc, nc)
    if (is.null(Q)) 
        Q <- rep(1, (nc - 1) * nc/2L)
    
    bf = as.double(bf)
    eig <- edQt(Q = Q, bf = bf)
    pos = 1
    k = as.integer(k)
    w = as.double(w <- rep(1/k, k))
    g = as.double(discrete.gamma(shape,k))
    fun <- function(s) -(nc - 1)/nc * log(1 - nc/(nc - 1) * s)
    eps <- (nc - 1)/nc
    n = as.integer(dim(contrast)[1])
    ind1 = rep(1:n, n:1)
    ind2 = unlist(lapply(n:1, function(x) seq_len(x) + n - x))
    li <- as.integer(length(ind1))
    weight = as.double(attr(x, "weight"))
    ll.0 = as.double(weight * 0)
    if (exclude == "pairwise") {
        index = con[ind1] & con[ind2]
        index = which(!index)
    }
    tmp = (contrast %*% eig[[2]])[ind1, ] * (contrast %*% (t(eig[[3]]) * bf))[ind2, ]
    tmp2 = vector("list", k)
    wdiag = .Call("PWI", as.integer(1:n), as.integer(1:n), as.integer(n), 
                  as.integer(n), rep(1, n), as.integer(li), PACKAGE = "phangorn")
    wdiag = which(wdiag > 0)
    
    tmp2 = vector("list", k)
    for (i in 1:(l - 1)) {
        for (j in (i + 1):l) {
            w0 = .Call("PWI", as.integer(x[[i]]), as.integer(x[[j]]), 
                       nr, n, weight, li, PACKAGE = "phangorn")
            if (exclude == "pairwise") 
                w0[index] = 0.0
            ind = w0 > 0
            
            old.el <- 1 - (sum(w0[wdiag])/sum(w0))
            if (old.el > eps) 
                old.el <- 10
            else old.el <- fun(old.el)
            #        sind = sum(ind)
            #           browser()
            
            for(lk in 1:k) tmp2[[lk]] <- tmp[ind, , drop = FALSE]
            # FS0 verwenden!!!        
            res <- .Call("FS5", eig, nc, as.double(old.el), w, g, tmp2, as.integer(k), as.integer(sum(ind)), 
                         bf, w0[ind], ll.0[ind], PACKAGE = "phangorn")
            d[pos] <- res[1] # res[[1]]
            v[pos] <- res[2] # res[[2]]
            pos = pos + 1
        }
    }
    attr(d, "Size") <- l
    if (is.list(x)) 
        attr(d, "Labels") <- names(x)
    else attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "variance") <- v
    class(d) <- "dist"
    return(d)
}


dist.mlOld <- function (x, model = "JC69", exclude = "none", bf = NULL, Q = NULL, ...) 
{
    if (!inherits(x,"phyDat")) 
        stop("x has to be element of class phyDat")
    l = length(x)
    d = numeric((l * (l - 1))/2)
    v = numeric((l * (l - 1))/2)
    contrast <- attr(x, "contrast")
    nc <- as.integer(attr(x, "nc"))
    nr <- as.integer(attr(x, "nr"))
    con = rowSums(contrast > 0) < 2
    if (exclude == "all") {
        index = con[x[[1]]]
        for (i in 2:l) index = index & con[x[[i]]]
        index = which(index)
        x = subset(x, , index)
    }
#    model <- match.arg(model, c("JC69", "WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
    model <- match.arg(model, c("JC69", .aamodels))
#    if (!is.na(match(model, c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24")))) 
    if (!is.na(match(model, .aamodels))) 
        getModelAA(model, bf = is.null(bf), Q = is.null(Q))
    if (is.null(bf)) 
        bf <- rep(1/nc, nc)
    if (is.null(Q)) 
        Q <- rep(1, (nc - 1) * nc/2L)

    bf = as.double(bf)
    eig <- edQt(Q = Q, bf = bf)
    k = 1
    w = as.double(1)
    g = as.double(1)
    fun <- function(s) -(nc - 1)/nc * log(1 - nc/(nc - 1) * s)
    eps <- (nc - 1)/nc
    n = as.integer(dim(contrast)[1])
    ind1 = rep(1:n, n:1)
    ind2 = unlist(lapply(n:1, function(x) seq_len(x) + n - x))
    li <- as.integer(length(ind1))
    weight = as.double(attr(x, "weight"))
    ll.0 = as.double(weight * 0)
    if (exclude == "pairwise") {
        index = con[ind1] & con[ind2]
        index = which(!index)
    }
    tmp = (contrast %*% eig[[2]])[ind1, ] * (contrast %*% (t(eig[[3]]) * bf))[ind2, ]
    tmp2 = vector("list", k)
    wdiag = .Call("PWI", as.integer(1:n), as.integer(1:n), as.integer(n), 
        as.integer(n), rep(1, n), as.integer(li), PACKAGE = "phangorn")
    wdiag = which(wdiag > 0)
    for (i in 1:(l - 1)) {
        for (j in (i + 1):l) {
            w0 = .Call("PWI", as.integer(x[[i]]), as.integer(x[[j]]), 
                nr, n, weight, li, PACKAGE = "phangorn")
            if (exclude == "pairwise") 
                w0[index] = 0.0
            ind = w0 > 0
            
            old.el <- 1 - (sum(w0[wdiag])/sum(w0))
            if (old.el > eps) 
                old.el <- 10
            else old.el <- fun(old.el)
    #        sind = sum(ind)
    #        tmp2 = vector("list", k)
            tmp2[[1]] <- tmp[ind, , drop = FALSE]
    # FS0 verwenden!!!        
            res <- .Call("FS5", eig, nc, as.double(old.el), w, g, tmp2, 1L, as.integer(sum(ind)), 
                bf, w0[ind], ll.0, PACKAGE = "phangorn")
            d[k] <- res[1] # res[[1]]
            v[k] <- res[2] # res[[2]]
            k = k + 1
        }
    }
    attr(d, "Size") <- l
    if (is.list(x)) 
        attr(d, "Labels") <- names(x)
    else attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "variance") <- v
    class(d) <- "dist"
    return(d)
} 

   
dist.logDet = function (x) 
{
    if (!inherits(x,"phyDat")) 
        stop("x has to be element of class phyDat")
    weight <- attr(x, "weight")
    contrast <- attr(x, 'contrast')
    r <- attr(x, "nc")
    l = length(x)
    d = numeric((l * (l - 1))/2)
    k = 1
    for (i in 1:(l - 1)) {
        Xi = weight * contrast[x[[i]], , drop=FALSE]
        for (j in (i + 1):l) {
            tmp = crossprod(Xi, contrast[x[[j]], , drop=FALSE])
            class(tmp) = "matrix"
            z = determinant.matrix(tmp, logarithm=TRUE)  
            res = z$sign*z$modulus
            if (is.nan(res)) {
                d[k] = 10
            }
            else d[k] = (-res + sum(log(rowSums(tmp) * colSums(tmp)))/2)/r
            k = k + 1
        }
    }
    attr(d, "Size") <- l
    if (is.list(x)) 
        attr(d, "Labels") <- names(x)
    else attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- "logDet"
    class(d) <- "dist"
    return(d)
}


readDist <- function(file){ #, format="phylip"
    tmp <- read.table(file, skip=1, stringsAsFactors = FALSE)
    labels = tmp[,1]
    dm <- as.matrix(tmp[,-1]) 
    dimnames(dm)=list(labels, labels)    
    as.dist(dm)
}
    
    
writeDist <- function(dm, file=""){ # , format="phylip"
    dm <- as.matrix(dm)
    cat(ncol(dm), "\n", file=file)
    write.table(dm, file, append=TRUE, quote=FALSE, col.names=FALSE)
}    
    
    


