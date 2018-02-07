#
# dist
#


#' Pairwise Distances from Sequences
#' 
#' \code{dist.hamming}, \code{dist.ml} and \code{dist.logDet} compute pairwise
#' distances for an object of class \code{phyDat}.  \code{dist.ml} uses DNA /
#' AA sequences to compute distances under different substitution models.
#' 
#' So far 17 amino acid models are supported ("WAG", "JTT", "LG", "Dayhoff",
#' "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb",
#' "FLU", "Blossum62", "Dayhoff_DCMut" and "JTT_DCMut") and additional rate
#' matrices and frequencies can be supplied.
#' 
#' The "F81" model uses empirical base frequencies, the "JC69" equal base
#' frequencies. This is even the case if the data are not nucleotides.
#' 
#' @param x An object of class \code{phyDat}
#' @param ratio Compute uncorrected ('p') distance or character difference.
#' @param model One of "JC69", "F81" or one of 17 amino acid models see
#' details.
#' @param exclude One of "none", "all", "pairwise" indicating whether to delete
#' the sites with missing data (or ambiguous states). The default is handle
#' missing data as in pml.
#' @param bf A vector of base frequencies.
#' @param Q A vector containing the lower triangular part of the rate matrix.
#' @param k Number of intervals of the discrete gamma distribution.
#' @param shape Shape parameter of the gamma distribution.
#' @param \dots Further arguments passed to or from other methods.
#' @return an object of class \code{dist}
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso For more distance methods for nucleotide data see
#' \code{\link[ape]{dist.dna}} and \code{\link{dist.p}} for pairwise
#' polymorphism p-distances. \code{\link{writeDist}} for export and import 
#' distances.
#' @references Lockhart, P. J., Steel, M. A., Hendy, M. D. and Penny, D. (1994)
#' Recovering evolutionary trees under a more realistic model of sequence
#' evolution. \emph{Molecular Biology and Evolution}, \bold{11}, 605--602.
#' 
#' Jukes TH and Cantor CR (1969). \emph{Evolution of Protein Molecules}. New
#' York: Academic Press. 21--132.
#' @keywords cluster
#' @examples
#' 
#' data(Laurasiatherian)
#' dm1 <- dist.hamming(Laurasiatherian)
#' tree1 <- NJ(dm1)
#' dm2 <- dist.logDet(Laurasiatherian)
#' tree2 <- NJ(dm2)
#' treedist(tree1,tree2)
#' # JC model
#' dm3 <- dist.ml(Laurasiatherian)
#' tree3 <- NJ(dm3)
#' treedist(tree1,tree3)
#' # F81 + Gamma
#' dm4 <- dist.ml(Laurasiatherian, model="F81", k=4, shape=.4)
#' tree4 <- NJ(dm4)
#' treedist(tree1,tree4)
#' treedist(tree3,tree4)
#' 
#' @rdname dist.hamming 
#' @export 
dist.hamming <- function (x, ratio = TRUE, exclude = "none") 
{
    if (!inherits(x,"phyDat")) 
        stop("x must be of class phyDat")
    l <- length(x)

    contrast <- attr(x, "contrast")
    nc <- as.integer(attr(x, "nc"))
    con <- rowSums(contrast > 0) < 2
    if (exclude == "all") {
        index <- con[x[[1]]]
        for (i in 2:l) index <- index & con[x[[i]]]
        index <- which(index)
        x <- subset(x, , index)
    }
    weight <- attr(x, "weight")  
    d <- numeric((l * (l - 1))/2)

    if(exclude == "pairwise"){
        k <- 1
        W <- numeric(l*(l-1)/2)
        for (i in 1:(l - 1)) {
            tmp <- con[x[[i]]] 
            for (j in (i + 1):l) {
                W[k] <- sum(weight[tmp & con[ x[[j]] ] ])
                k <- k + 1
            }
        }  
             
    } 

    if(nc > 31){
#        contrast <- attr(x, "contrast")
        k <- 1
        for (i in 1:(l - 1)) {
            X <- contrast[x[[i]], , drop = FALSE]
            for (j in (i + 1):l) {
                d[k] <- sum(weight * (rowSums(X * contrast[x[[j]], , drop = FALSE]) == 0))
                k <- k + 1
            }
        }
    } # end if  
    else{
        nr <- attr(x, "nr")
        if(exclude == "pairwise")ind <- which(con[unlist(x)]==FALSE)  
        x <- prepareDataFitch(x) 
        if(exclude == "pairwise")x[ind] <- as.integer(2L^nc -1L) 
        res <- .C("distHamming", as.integer(x), as.double(weight), 
            as.integer(nr), as.integer(l), as.double(d), PACKAGE = "phangorn")
        d <- res[[5]]
    }     

    if (ratio){
        if(exclude == "pairwise") d <- d/W
        else d <- d/sum(weight)
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


#' @rdname dist.hamming
#' @export
dist.ml <- function (x, model = "JC69", exclude = "none", bf = NULL, Q = NULL, 
                     k=1L, shape=1, ...) 
{
    if (!inherits(x,"phyDat")) 
        stop("x must be of class phyDat")
    l <- length(x)
    d <- numeric((l * (l - 1))/2)
    v <- numeric((l * (l - 1))/2)
    contrast <- attr(x, "contrast")
    con <- rowSums(contrast > 0) < 2
    if (exclude == "all") {
        index <- con[x[[1]]]
        for (i in 2:l) index <- index & con[x[[i]]]
        index <- which(index)
        x <- subset(x, , index)
    }
    nc <- as.integer(attr(x, "nc"))
    nr <- as.integer(attr(x, "nr"))
    model <- match.arg(model, c("JC69", "F81", .aamodels))
    if (!is.na(match(model, .aamodels))) 
        getModelAA(model, bf = is.null(bf), Q = is.null(Q))
    if(is.null(bf) && model=="F81") bf <- baseFreq(x)
    if (is.null(bf)) 
        bf <- rep(1/nc, nc)
    if (is.null(Q)) 
        Q <- rep(1, (nc - 1) * nc/2L)
    
    bf <- as.double(bf)
    eig <- edQt(Q = Q, bf = bf)
    pos <- 1
    k <- as.integer(k)
    w <- as.double(w <- rep(1/k, k))
    g <- as.double(discrete.gamma(shape,k))
    fun <- function(s) -(nc - 1)/nc * log(1 - nc/(nc - 1) * s)
    eps <- (nc - 1)/nc
    n <- as.integer(dim(contrast)[1])
    ind1 <- rep(1:n, n:1)
    ind2 <- unlist(lapply(n:1, function(x) seq_len(x) + n - x))
    li <- as.integer(length(ind1))
    weight <- as.double(attr(x, "weight"))
    ll.0 <- as.double(weight * 0)
    if (exclude == "pairwise") {
        index <- con[ind1] & con[ind2]
        index <- which(!index)
    }
    tmp <- (contrast %*% eig[[2]])[ind1, ] * (contrast %*% (t(eig[[3]]) * bf))[ind2, ]
    tmp2 <- vector("list", k)
#    wdiag = .Call("PWI", as.integer(1:n), as.integer(1:n), as.integer(n), 
#                  as.integer(n), rep(1, n), as.integer(li), PACKAGE = "phangorn")
#    wdiag = which(wdiag > 0)
    
    wshared <- which(rowSums(contrast[ind1, ] * contrast[ind2, ]) > 0)
    
    tmp2 <- vector("list", k)
    for (i in 1:(l - 1)) {
        for (j in (i + 1):l) {
            w0 <- .Call("PWI", as.integer(x[[i]]), as.integer(x[[j]]), 
                       nr, n, weight, li, PACKAGE = "phangorn")
            if (exclude == "pairwise") 
                w0[index] <- 0.0
            ind <- w0 > 0
            
            old.el <- 1 - (sum(w0[wshared])/sum(w0))
            if (old.el > eps) 
                old.el <- 10
            else old.el <- fun(old.el)
            #        sind = sum(ind)
            #           browser()
            
            for(lk in 1:k) tmp2[[lk]] <- tmp[ind, , drop = FALSE]
            # FS0 verwenden!!!        
            res <- .Call("FS5", eig, nc, as.double(old.el), w, g, tmp2, 
                    as.integer(k), as.integer(sum(ind)), 
                    w0[ind], ll.0, PACKAGE = "phangorn") 
            d[pos] <- res[1] # res[[1]]
            v[pos] <- res[2] # res[[2]]
            pos <- pos + 1
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


#' @rdname dist.hamming
#' @export   
dist.logDet <- function (x) 
{
    if (!inherits(x,"phyDat")) 
        stop("x must be of class phyDat")
    weight <- attr(x, "weight")
    contrast <- attr(x, 'contrast')
    r <- attr(x, "nc")
    l <- length(x)
    d <- numeric((l * (l - 1))/2)
    k <- 1
    for (i in 1:(l - 1)) {
        Xi <- weight * contrast[x[[i]], , drop=FALSE]
        for (j in (i + 1):l) {
            tmp <- crossprod(Xi, contrast[x[[j]], , drop=FALSE])
            class(tmp) <- "matrix"
            z <- determinant.matrix(tmp, logarithm=TRUE)  
            res <- z$sign*z$modulus
            if (is.nan(res)) {
                d[k] <- 10
            }
            else d[k] <- (-res + sum(log(rowSums(tmp) * colSums(tmp)))/2)/r
            k <- k + 1
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




#' Writing and reading distances in phylip and nexus format
#' 
#' \code{readDist}, \code{writeDist} and \code{write.nexus.dist} are useful to 
#' exchange distance matrices with other phylogenetic programs.  
#' 
#' 
#' @param x A \code{dist} object.
#' @param file A file name.
#' @param format file format, default is "phylip", only other option so far is 
#' "nexus". 
#' @param \dots Further arguments passed to or from other methods. 
#' @param upper	logical value indicating whether the upper triangle of the 
#' distance matrix should be printed.
#' @param diag	logical value indicating whether the diagonal of the distance 
#' matrix should be printed.
#' @param digits passed to format inside of \code{write.nexus.dist}.
#' @param taxa logical. If TRUE a taxa block is added.
#' @param append logical. If TRUE the nexus blocks will be added to a file.
#' @return an object of class \code{dist}
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso To compute distance matrices see \code{\link{dist.ml}}
#' \code{\link[ape]{dist.dna}} and \code{\link{dist.p}} for pairwise
#' polymorphism p-distances
#' @references Maddison, D. R., Swofford, D. L. and Maddison, W. P. (1997) 
#' NEXUS: an extensible file format for systematic information. 
#' \emph{Systematic Biology}, \bold{46}, 590--621.
#' 
#' @keywords cluster
#' @examples
#' 
#' data(yeast)
#' dm <- dist.ml(yeast)
#' writeDist(dm)
#' write.nexus.dist(dm)
#' 
#' @rdname writeDist 
#' @export writeDist
writeDist <- function(x, file="", format="phylip", ...){ 
    format <- match.arg(format, c("phylip", "nexus"))
    if(format=="phylip"){
        x <- as.matrix(x)
#        x <- format(x, digits = digits, justify = "none")
        cat(ncol(x), "\n", file=file)
        write.table(x, file, append=TRUE, quote=FALSE, col.names=FALSE)
    }
    else write.nexus.dist(x, file=file, ...)
}    
    

#' @rdname writeDist
#' @export 
write.nexus.dist <- function(x, file="", append=FALSE, upper=FALSE, diag=TRUE, 
                             digits = getOption("digits"), taxa=!append){
    
    taxa.labels <- attr(x, "Labels")
    ntaxa <- length(taxa.labels)
    
    m <- as.matrix(x)
    cf <- format(m, digits = digits, justify = "none")    
    l <-  length(attr(x, "Labels"))
    if(upper) diag <- TRUE
    if (!upper) cf[row(cf) < col(cf)] <- ""
    if (!diag) cf[row(cf) == col(cf)] <- ""    
    
    cf2 <- apply(cf, 1, "paste0", collapse=" ")
    cf2 <- paste(attr(x, "Labels"), cf2)
    cf2 <- trimws(cf2, "right")
    
    if(!append)cat("#NEXUS\n\n", file = file)  
    
    if(taxa){
        cat(paste("BEGIN TAXA;\n\tDIMENSIONS ntax=", ntaxa, ";\n", 
                  sep = ""), file = file, append = TRUE)
        cat("\tTAXLABELS", paste(taxa.labels, sep = " "), ";\nEND;\n\n", 
            file = file, append = TRUE)
    }
    
    cat("BEGIN DISTANCES; \n", file = file, append = TRUE)
    if(upper) cat("\tFORMAT TRIANGLE = BOTH;\n", file = file, append = TRUE)
    else cat("\tFORMAT TRIANGLE = LOWER;\n", file = file, append = TRUE)
    if(!diag) cat("\tFORMAT NODIAGONAL;\n", file = file, append = TRUE)
    cat("\tMatrix \n", file = file, append = TRUE)
    for(i in 1:l) cat("\t", cf2[i], "\n", sep="", file = file, append = TRUE)
    cat("\t;\nEND; \n", file = file, append = TRUE)
}



RSS <- function(x, dm, trace=0){
    labels <- attr(x, "labels")
    dm <- as.matrix(dm)
    k <- dim(dm)[1]
    dm <- dm[labels,labels]
    y <- dm[lower.tri(dm)]
    betahat <- attr(x, "weights")
    X <- splits2design(x)
    RSS <- sum((y-(X%*%betahat))^2)
    RSS
}    



#' @rdname writeDist
#' @export  
readDist <- function(file){ #, format="phylip"
    tmp <- read.table(file, skip=1, stringsAsFactors = FALSE)
    labels <- tmp[,1]
    dm <- as.matrix(tmp[,-1]) 
    dimnames(dm) <- list(labels, labels)    
    as.dist(dm)
}


#' @rdname writeDist
#' @param incomparables Not used so far.
#' @export 
unique.dist <-  function(x, incomparables, ...){
    y <- as.matrix(x) 
    l <- nrow(y)
    z <- character(l)
    for(i in seq_len(l)) z[i] <- paste( round(y[i, ] ,8), collapse="_")
    if(any(duplicated(z))){
        ind <- !duplicated(z)
        y <- y[ind, ind]
        if(is.matrix(x))return(y)
        if(inherits(x, "dist")){
            #            attrib <- attributes(x)
            res <- as.dist(y, diag=attr(x, "Diag"), upper=attr(x, "Upper"))  
            return(res)
        }
    }
    x
}
