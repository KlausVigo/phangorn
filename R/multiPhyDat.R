
read.multiPhyDat <- function(files, ...){
    gene.names <- gsub(".fasta","",sapply(strsplit(files, "/"), tail, 1))
    l <- length(files)
    res <- vector("list", l)
    for(i in 1:l) res[[i]] <- read.phyDat(files[i], ...)
    class(res) <- "multiPhyDat"
    res
}


read.multiphydatFASTA <- function(files){
    gene.names <- gsub(".fasta","",sapply(strsplit(files, "/"), tail, 1))
    dna <- lapply(files, read.phyDatFASTA)
    names(dna) <- gene.names
    out <- new("multidna", dna=dna)
    return(out)
}


read.FASTA.AA <- function (file) 
{
    if (length(grep("^(ht|f)tp:", file))) {
        url <- file
        file <- tempfile()
        download.file(url, file)
    }
    sz <- file.info(file)$size
    x <- readBin(file, "raw", sz)
    icr <- which(x == as.raw(13))
    if (length(icr)) 
        x <- x[-icr]
    res <- .Call("rawStream2phyDat", x)
 
    aa <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", 
        "l", "k", "m", "f", "p", "s", "t", "w", "y", "v")
    aa2 <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", 
         "l", "k", "m", "f", "p", "s", "t", "w", "y", "v", "b", 
         "z", "x", "-", "?")
    AA <- diag(20)
    AA <- rbind(AA, matrix(0, 5, 20))
    AA[21, 3] <- AA[21, 4] <- 1 # Aspartate or Asparagine
    AA[22, 6] <- AA[22, 7] <- 1 #
    AA[23:25, ] = 1
    dimnames(AA) <- list(aa2, aa)

    ddd = fast.table(res)

    data = ddd$data
    names(data) <- sub("^ +", "", names(data))
    row.names(data) = NULL

    attr(data, "row.names") = NULL
    attr(data, "weight") = ddd$weight
    attr(data, "nr") = length(ddd$weight)
    attr(data, "nc") = 20
    attr(data, "index") = as.integer(ddd$index)
    attr(data, "levels") = aa
    attr(data, "allLevels") = aa2
    attr(data, "type") = "AA"
    attr(data, "contrast") = AA    
    class(data) = "phyDat"
    data
}


cbindPD <- function(..., gaps="-"){
    object <- as.list(substitute(list(...)))[-1]    
    x <- list(...)
    n <- length(x) 
    if (n == 1) 
        return(x[[1]])
    type <- attr(x[[1]], "type")
    nr = numeric(n)
    
    ATTR <- attributes(x[[1]])
    
    nr[1] <- sum(attr(x[[1]], "weight"))
    levels <- attr(x[[1]], "levels")
    allLevels <- attr(x[[1]], "allLevels")
    gapsInd <- match(gaps, allLevels)
    snames <- vector("list", n)  # names(x[[1]])
    vec = numeric(n+1)
    wvec = numeric(n+1)
    objNames<-as.character(object)
    if(any(duplicated(objNames))) objNames <- paste(objNames,1:n,sep="")
#    tmp <- as.character(x[[1]])
     
    for(i in 1:n){
        snames[[i]] = names(x[[i]]) 
        nr[i] <- attr(x[[i]], "nr") 
        vec[i+1] = attr(x[[i]], "nr")
        wvec[i+1] = sum(attr(x[[i]], "weight"))
    }
    vec = cumsum(vec)
    wvec = cumsum(wvec)
    snames = unique(unlist(snames))
    weight <- numeric(vec[n+1])

    index <- numeric(wvec[n+1]) 
    ATTR$names <- snames
    ATTR$nr <- vec[n+1]

    tmp = matrix(gapsInd, vec[n+1], length(snames), dimnames = list(NULL, snames))
    tmp <- as.data.frame(tmp)

    for(i in 1:n){
        nam = names(x[[i]])
        tmp[(vec[i]+1):vec[i+1], nam] <- x[[i]][nam]
        weight[(vec[i]+1):vec[i+1]] <- attr(x[[i]], "weight")
        index[(wvec[i]+1):wvec[i+1]] <- attr(x[[i]], "index")
    }
    ATTR$index <- index
    ATTR$weight <- weight
    attributes(tmp) <- ATTR
    tmp
}




as.phyDat.DNAbin = function (data) 
{
    if (is.matrix(data)) 
        nam = row.names(data)
    else nam = names(data)
    if (class(data) == "DNAbin") 
        data = as.character(data)
    if (is.matrix(data)) 
        data = as.data.frame(t(data), stringsAsFactors = FALSE)
    else data = as.data.frame(data, stringsAsFactors = FALSE)
    
    data = data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)
    
    ac = c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y", 
           "k", "v", "h", "d", "b", "n", "?", "-")
    AC = matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1), 
                  c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1), 
                  c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1), 
                  c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1)), 
                  18, 4, dimnames = list(NULL, c("a", "c", "g", "t")))
    ddd = fast.table(data)
    data = ddd$data
    index = ddd$index
    q = length(data)
    p = length(data[[1]])
    d = dim(data)
    att = attributes(data) 
    data = match(unlist(data), ac)
    attr(data, "dim") = d
    data = as.data.frame(data, stringsAsFactors=FALSE)
    attributes(data) = att
    
    row.names(data) = as.character(1:p)
    data = na.omit(data)
    rn = as.numeric(rownames(data))
    
    aaa = match(index, attr(data, "na.action"))
    index = index[is.na(aaa)] 
    index = match(index, unique(index))
    rn = as.numeric(rownames(data))
    attr(data, "na.action") = NULL
    
    weight = ddd$weight[rn]
    p = dim(data)[1]
    names(data) = nam
    attr(data, "row.names") = NULL 
    attr(data, "weight") = weight
    attr(data, "nr") = p
    attr(data, "nc") = 4
    if (return.index) 
        attr(data, "index") = index
    attr(data, "levels") = c("a", "c", "g", "t")
    attr(data, "allLevels") = ac
    attr(data, "type") = "DNA"
    attr(data, "contrast") = AC
    class(data) = "phyDat"
    data
}



# for Emmanuel: here is some faster code 

# this is not faster
as.character.DNAbin2 <- function (x, ...) 
{
    f <- function(xx) {
        ans <- character(length(xx))
        bs <- as.raw(ape:::._bs_)
        ans <- ape:::._cs_[match(xx, bs)]
#        for (i in 1:15) ans[which(xx == ._bs_[i])] <- ._cs_[i]
#        ans[which(xx == 4)] <- "-"
#        ans[which(xx == 2)] <- "?"
        if (is.matrix(xx)) {
            dim(ans) <- dim(xx)
            dimnames(ans) <- dimnames(xx)
        }
        ans
    }
    if (is.list(x)) lapply(x, f)
    else f(x)
}


as.DNAbin.character2 <- function (x, ...) 
{
    n <- length(x)
    ans <- raw(n)
#    for (i in 1:15) ans[which(x == ._cs_[i])] <- as.raw(._bs_[i])
#    ans[which(x == "-")] <- as.raw(4)
#    ans[which(x == "?")] <- as.raw(2)
    ans <- as.raw(ape:::._bs_)[match(x, ape:::._cs_)]
# maybe add some NA handling
# ans[is.na(ans)] <- as.raw(0) 
    if (is.matrix(x)) {
        dim(ans) <- dim(x)
        dimnames(ans) <- dimnames(x)
    }
    class(ans) <- "DNAbin"
    ans
}


#  cy <- as.character(yeast)
#  system.time(blub1 <- as.DNAbin(cy))
#  system.time(blub2 <- as.DNAbin.character2(cy))
#  all.equal(blub2, blub1)

#Browse[2]> names(attributes(L1))
#[1] "names"     "class"     "weight"    "nr"        "nc"        "index"    
#[7]   "levels"    "allLevels" "type"      "contrast" 

