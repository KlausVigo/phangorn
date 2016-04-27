#
# Data structures for ML and MP
# 
fast.table <- function (data)                                                            
{                                                                                 
    if(!is.data.frame(data)) 
        data = as.data.frame(data, stringsAsFactors = FALSE)                    
    da = do.call("paste", c(data, sep = "\r"))                                             
    ind = !duplicated(da)                                                                  
    levels = da[ind]                                                                       
    cat <- factor(da,levels = levels)                                                      
    nl <- length(levels(cat))                                                        
    bin <- (as.integer(cat) - 1)                                                           
    pd <- nl                                                                               
    bin <- bin[!is.na(bin)]                                                                
    if (length(bin)) bin <- bin + 1                                                        
    y <- tabulate(bin, pd)                                                                 
    result=list(index = bin, weights = y, data = data[ind,])                                                                                  
    result                                                                                 
}                                                                                        


phyDat.default <- function (data, levels = NULL, return.index = TRUE, contrast = NULL, 
    ambiguity = "?", compress=TRUE, ...) 
{
    if (is.matrix(data)) 
        nam = row.names(data)
    else nam = names(data)
    if(is.null(nam))stop("data object must contain taxa names")
    if (inherits(data,"DNAbin")) 
        data = as.character(data)
    if (is.matrix(data)) 
        data = as.data.frame(t(data), stringsAsFactors = FALSE)
# new 4.4.2016 bug fix (reported by Eli Levy Karin)     
    if (is.vector(data) && !is.list(data))data = as.data.frame(t(data), stringsAsFactors = FALSE)
    else data = as.data.frame(data, stringsAsFactors = FALSE)
#    data = data.frame(as.matrix(data), stringsAsFactors = FALSE)    
    
    
    if(length(data[[1]])==1) compress=FALSE 
    if(compress){
        ddd = fast.table(data)
        data = ddd$data
        weight = ddd$weight
        index = ddd$index
    }
    else{
        p = length(data[[1]])
        weight = rep(1, p)
        index = 1:p
    }
    q = length(data)
    p = length(data[[1]])
    tmp <- vector("list", q)
    if (!is.null(contrast)) {
        levels = colnames(contrast)
        all.levels = rownames(contrast)
        rownames(contrast) = NULL
    }
    else {
        if (is.null(levels)) 
            stop("Either argument levels or contrast has to be supplied")
        l = length(levels)
        contrast = diag(l)
        all.levels = levels
        if (!is.null(ambiguity)) {
            all.levels = c(all.levels, ambiguity)
            k = length(ambiguity)
            if (k > 0) 
                contrast = rbind(contrast, matrix(1, k, l))
        }
    }
    att = attributes(data) 
    data = lapply(data, match, all.levels) # avoid unlist   
    attributes(data) = att
    
    row.names(data) = as.character(1:p)
    data = na.omit(data)
    
    aaa = match(index, attr(data, "na.action"))
    index = index[is.na(aaa)] 
    index = match(index, unique(index))
    rn = as.numeric(rownames(data))
    attr(data, "na.action") = NULL        
    weight = weight[rn] 
    p = dim(data)[1]
    names(data) = nam
    attr(data, "row.names") = NULL
    attr(data, "weight") = weight
    attr(data, "nr") = p
    attr(data, "nc") = length(levels)
    if (return.index) 
        attr(data, "index") = index
    attr(data, "levels") = levels
    attr(data, "allLevels") = all.levels
    attr(data, "type") = "USER"
    attr(data, "contrast") = contrast
    class(data) = "phyDat"
    data
}


phyDat.DNA = function (data, return.index = TRUE) 
{
    if (is.matrix(data)) 
        nam = row.names(data)
    else nam = names(data)
    if (inherits(data,"DNAbin")) 
        data = as.character(data)
    if (is.matrix(data)) 
        data = as.data.frame(t(data), stringsAsFactors = FALSE)
    else data = as.data.frame(data, stringsAsFactors = FALSE)

    data = data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)
 
    ac = c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y", 
        "k", "v", "h", "d", "b", "n", "?", "-")
    AC = matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 
        0, 1, 1, 1), c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 
        0, 1, 1, 1, 1), c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 1, 1), c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 1, 1, 1)), 18, 4, dimnames = list(NULL, c("a", 
        "c", "g", "t")))
    
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


phyDat.AA <- function (data, return.index = TRUE) 
{
    if(is.matrix(data)) nam = row.names(data)
    else nam = names(data)  
    if (inherits(data,"DNAbin")) 
        data = as.character(data)
    if (is.matrix(data)) 
        data = as.data.frame(t(data), stringsAsFactors = FALSE)
    else data = as.data.frame(data, stringsAsFactors = FALSE)
  
    data = data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)

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
    
    ddd = fast.table(data)
    data = ddd$data
    index = ddd$index
    q = length(data)
    p = length(data[[1]])
    tmp <- vector("list", q)

    d = dim(data)
    att = attributes(data) 
    data = match(unlist(data), aa2)
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
    attr(data, "nc") = 20
    if (return.index) 
        attr(data, "index") = index
    attr(data, "levels") = aa
    attr(data, "allLevels") = aa2
    attr(data, "type") = "AA"
    attr(data, "contrast") = AA    
    class(data) = "phyDat"
    data
}



phyDat.codon <- function (data, return.index = TRUE) 
{
    if(is.matrix(data)) nam = row.names(data)
    else nam = names(data)  
    if (inherits(data,"DNAbin")) 
        data = as.character(data)

    if (is.matrix(data)) 
        data = as.data.frame(t(data), stringsAsFactors = FALSE)
    else data = as.data.frame(data, stringsAsFactors = FALSE)
    
    data = data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)

    data[data=="u"] = "t" 

    splseq = function (seq, frame = 0) 
    {
        starts <- seq(from = frame + 1, to = length(seq), by = 3L)
        sapply(starts, function(x) paste(seq[x:(x + 2L)], collapse=""))
    } 
 
    data = sapply(data, splseq)
    
    ddd = fast.table(data)
    codon = c("aaa", "aac", "aag", "aat", "aca", "acc", "acg", "act", 
      "aga", "agc", "agg", "agt", "ata", "atc", "atg", "att", 
      "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct", "cga", 
      "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt", "gaa", "gac", 
      "gag", "gat", "gca", "gcc", "gcg", "gct", "gga", "ggc", "ggg", 
      "ggt", "gta", "gtc", "gtg", "gtt", "tac", "tat", 
      "tca", "tcc", "tcg", "tct", "tgc", "tgg", "tgt", "tta", 
      "ttc", "ttg", "ttt")
# ohne Stopcodons "taa", "tag", "tga",     

    CODON <- diag(61)
    dimnames(CODON) <- list(codon, codon)

    data = ddd$data
    index = ddd$index
    q = length(data)
    p = length(data[[1]])
    tmp <- vector("list", q)

    d = dim(data)
    att = attributes(data) 
    data = match(unlist(data), codon)
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
    attr(data, "nc") = 61
    if (return.index) 
        attr(data, "index") = index
    attr(data, "levels") = codon
    attr(data, "allLevels") = codon
    attr(data, "type") = "CODON"
    attr(data, "contrast") = CODON    
    class(data) = "phyDat"
    data
}


as.phyDat <- function (x, ...){
    if (inherits(x,"phyDat")) return(x)
    UseMethod("as.phyDat")
}


as.phyDat.DNAbin <- function(x,...) phyDat.DNA(x,...)




as.phyDat.alignment <- function (x, type="DNA",...) 
{
    x$seq <- tolower(x$seq)
    data <- sapply(x$seq, strsplit, "")
    names(data) <- x$nam
    if(type=="DNA") dat <- phyDat.DNA(data,...)
    if(type=="AA") dat <- phyDat.AA(data, ...)
    if(type=="CODON") dat <- phyDat.codon(data, ...)
    if(type=="USER") dat <- phyDat.default(data, ...)
    dat
}


#as.alignment.phyDat <- function(x, ...) as.alignment(as.character(x))

phyDat2alignment <-  function(x){
    z = as.character(x)
    nam = rownames(z)
    type = attr(x, "type")
    seq <- switch(type, 
                  DNA = tolower(apply(z, 1, paste, collapse="")), 
                  AA = toupper(apply(z, 1, paste, collapse="")))                     
    names(seq) <- NULL
    res <- list(nb=length(seq), nam=nam, seq=seq, com=NA)
    class(res) = "alignment"
    res
}


as.phyDat.MultipleAlignment <- function(x, ...){
    if (requireNamespace('Biostrings')){
    if(inherits(x, "DNAMultipleAlignment"))
        res <- phyDat.DNA(Biostrings::as.matrix(x))
    if(inherits(x, "RNAMultipleAlignment"))
        res <- phyDat.DNA(Biostrings::as.matrix(x))
    if(inherits(x, "AAMultipleAlignment"))
        res <- phyDat.AA(Biostrings::as.matrix(x))
    return(res)
    }
    return(NULL)
}


as.MultipleAlignment <- function (x, ...){
    if (inherits(x,"MultipleAlignment")) return(x)
    UseMethod("as.MultipleAlignment")
}


as.MultipleAlignment.phyDat <- function(x){
    if (requireNamespace('Biostrings')){
    z = as.character(x)
    nam = rownames(z)
    type = attr(x, "type")
    seq <- switch(type, 
                  DNA = tolower(apply(z, 1, paste, collapse="")), 
                  AA = toupper(apply(z, 1, paste, collapse="")))
    if(type=="DNA") return(Biostrings::DNAMultipleAlignment(seq))
    if(type=="AA") return(Biostrings::AAMultipleAlignment(seq))
    }
    return(NULL)
}


phyDat2MultipleAlignment <- as.MultipleAlignment.phyDat


as.phyDat.matrix <- function (x, ...) phyDat(data=x, ...)


as.phyDat.data.frame <- function (x, ...) phyDat(data=x, ...)
 

acgt2ry <- function(obj){
   ac = c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y", 
        "k", "v", "h", "d", "b", "n", "?", "-")
   AC = matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 
        0, 1, 1, 1), c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 
        0, 1, 1, 1, 1), c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 1, 1), c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 1, 1, 1)), 18, 4, dimnames = list(NULL, c("a", 
        "c", "g", "t")))
   ry = AC[c(7,10),]
   RY = AC %*% t(ry)
   RY[RY==2] = 1
   dimnames(RY) = list(NULL, c("r", "y"))
   attr(obj, "levels") = c("r", "y")
   attr(obj, "nc") = 2
   attr(obj, "type") = "USER"
   attr(obj, "contrast") = RY
   obj=phyDat.default(as.character(obj, allLevels=FALSE), levels = c("r", "y"), ambiguity = NULL)
   obj  
}


as.character.phyDat <- function (x, allLevels=TRUE, ...) 
{
    nr <- attr(x, "nr")
    nc <- attr(x, "nc")
    type <- attr(x, "type")
    if (type == "DNA") {
        labels <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", 
            "y", "k", "v", "h", "d", "b", "n", "?", "-")
    }
    if (type == "AA") {
        labels <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", 
            "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", 
            "v", "b", "z", "x", "-", "?")
    }
    if (type == "USER") {
        #levels
        if(allLevels)labels = attr(x, "allLevels")
        else{
            tmp = attr(x, "levels")
            contrast = attr(x, "contrast") # contrast=AC
            contrast[contrast>0] = 1
            ind = which(rowSums(contrast)==1)
            contrast[rowSums(contrast)>1,] = 0 
            labels = rep(NA, length(attr(x, "allLevels")))
            labels[ind] = tmp[contrast%*%c(1:length(tmp))]
            }
    }
    result = matrix(NA, nrow = length(x), ncol = nr)
    for (i in 1:length(x)) result[i, ] <- labels[x[[i]]]
    if (is.null(attr(x, "index"))) 
        index = rep(1:nr, attr(x, "weight"))
    else {
        index = attr(x, "index")
        if (is.data.frame(index)) 
            index <- index[, 1]
    }
    result = result[, index, drop = FALSE]
    rownames(result) = names(x)
    result
}


# replace as.character.phyDat 20 Zeilen weniger
as.character.phyDat2 <- function (x, ...) 
{
    nr <- attr(x, "nr")
    nc <- attr(x, "nc")
    type <- attr(x, "type")
    labels = attr(x, "allLevels")
    result = matrix(NA, nrow = length(x), ncol = nr)
    for (i in 1:length(x)) result[i, ] <- labels[x[[i]]]
    if (is.null(attr(x, "index"))) 
        index = rep(1:nr, attr(x, "weight"))
    else {
        index = attr(x, "index")
        if (is.data.frame(index)) 
            index <- index[, 1]
    }
    result = result[, index, drop = FALSE]
    rownames(result) = names(x)
    result
}


#as.data.frame.phyDat <- function(x, ...){
#   data.frame(t(as.character(x, ...)), stringsAsFactors=FALSE)
#}

# much faster
# TODO as stringsAsFactors=FALSE
# result[[i]] <- x[[i]] + factor levels setzen
# 
as.data.frame.phyDatOld <- function(x, ...){
    nr <- attr(x, "nr")
    nc <- attr(x, "nc")
    labels <- attr(x, "allLevels")
    result <- vector("list", length(x))
    for (i in 1:length(x)) result[[i]] <- labels[x[[i]]]
    attr(result, "names") <- names(x)
    attr(result, "row.names") <- 1:nr
    attr(result, "class") <- "data.frame"
    result
}


as.data.frame.phyDat <- function(x, ...){
  nr <- attr(x, "nr")
  nc <- attr(x, "nc")
  labels <- attr(x, "allLevels")
  result <- vector("list", length(x))
  if (is.null(attr(x, "index"))) 
    index = rep(1:nr, attr(x, "weight"))
  else {
    index = attr(x, "index")
    if (is.data.frame(index)) 
      index <- index[, 1]
  }
  for (i in 1:length(x)) result[[i]] <- labels[x[[i]][index]]
  attr(result, "names") <- names(x)
  attr(result, "row.names") <- 1:length(index)
  attr(result, "class") <- "data.frame"
  result
}


#as.DNAbin.phyDat <- function(x,...) {
#   if(attr(x, "type")=="DNA") return(as.DNAbin(as.character(x, ...)))
#   else stop("x must be a nucleotide sequence")
#}

# quite abit faster
as.DNAbin.phyDat <- function (x, ...) 
{
    if(attr(x, "type")=="DNA"){

    nr <- attr(x, "nr")
    ac = attr(x, "allLevels")
    result = matrix(as.raw(0), nrow = length(x), ncol = nr)
    # from ape ._cs_
    cs <- c("a", "g", "c", "t", "r", "m", "w", "s", "k", "y", "v", "h", 
      "d", "b", "n", "-", "?")
    # from ape ._bs_
    bs <- as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48, 224, 176, 208, 
                   112, 240, 4, 2))
    ord <- match(ac, cs)
    ord[5] <- 4

    for (i in 1:length(x)){
        ind <- ord[x[[i]]]
        result[i,] <- bs[ind]    
    }    
    if (is.null(attr(x, "index"))) 
        index = rep(1:nr, attr(x, "weight"))
    else {
        index = attr(x, "index")
        if (is.data.frame(index)) 
            index <- index[, 1]
    }
    result = result[, index, drop = FALSE]
    rownames(result) = names(x)
    class(result) <- "DNAbin"
    return(result)
    }
    else stop("x must be a nucleotide sequence")
}


if (getRversion() >= "2.15.1") utils::globalVariables("as.AAbin")
as.AAbin.phyDat <- function(x,...) {
   if(attr(x, "type")=="AA") return(as.AAbin(as.character(x, ...)))
   else stop("x must be a amino acid sequence")
}

 
phyDat <- function (data, type="DNA", levels=NULL, return.index = TRUE,...) 
{
    if (inherits(data,"DNAbin")) type <- "DNA"
    pt <- match.arg(type, c("DNA", "AA", "CODON", "USER"))  
    if(pt=="DNA") dat <- phyDat.DNA(data, return.index=return.index,...)
    if(pt=="AA") dat <- phyDat.AA(data, return.index=return.index, ...)
    if(pt=="CODON") dat <- phyDat.codon(data, return.index=return.index, ...)
    if(pt=="USER") dat <- phyDat.default(data, levels = levels, return.index=return.index, ...)
    dat
}


print.phyDat = function (x, ...) 
{
    cat(length(x), "sequences with",sum(attr(x,"weight")), "character and",attr(x,"nr"),"different site patterns.\n")
    cat("The states are",attr(x,"levels"), "\n")
}


c.phyDat <- function(...){
    object <- as.list(substitute(list(...)))[-1]    
    x <- list(...)
    n <- length(x) 
    match.names <- function(a,b){
        if(any(!(a %in% b)))stop("names do not match previous names") 
        }
    if (n == 1) 
        return(x[[1]])
    type <- attr(x[[1]], "type")
    nr = numeric(n)
    nr[1] <- sum(attr(x[[1]], "weight"))
    levels <- attr(x[[1]], "levels")
    snames <- names(x[[1]])
    objNames<-as.character(object)
    if(any(duplicated(objNames))) objNames <- paste(objNames,1:n,sep="")
    tmp <- as.character(x[[1]])
    for(i in 2:n){
        match.names(snames,names(x[[i]]))
        x[[i]] <- getCols(x[[i]],snames)
        nr[i] <- sum(attr(x[[i]], "weight"))
        tmp <- cbind(tmp, as.character(x[[i]]))
    }
    if (type == "DNA") 
        dat <- phyDat.DNA(tmp, return.index = TRUE)
    if (type == "AA") 
        dat <- phyDat.AA(tmp, return.index = TRUE)
    if (type == "USER") 
        dat <- phyDat.default(tmp, levels = levels, return.index = TRUE)
     if (type == "CODON") 
        dat <- phyDat.codon(tmp, return.index = TRUE)       
    attr(dat,"index") <- data.frame(index=attr(dat,"index"), genes=rep(objNames, nr))   
    dat
}


# new cbind.phyDat 
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



cbind.phyDat <- function(..., gaps="-"){
    object <- as.list(substitute(list(...)))[-1]    
    x <- list(...)
    n <- length(x) 
    if (n == 1) 
        return(x[[1]])
    type <- attr(x[[1]], "type")
    nr = numeric(n)
    nr[1] <- sum(attr(x[[1]], "weight"))
    levels <- attr(x[[1]], "levels")
    snames <- vector("list", n)  # names(x[[1]])
    vec = numeric(n+1)
    objNames<-as.character(object)
    if(any(duplicated(objNames))) objNames <- paste(objNames,1:n,sep="")
    tmp <- as.character(x[[1]])
    for(i in 1:n){
        snames[[i]] = names(x[[i]]) #match.names(snames,names(x[[i]]))
        nr[i] <- sum(attr(x[[i]], "weight")) 
        vec[i+1] = sum(attr(x[[i]], "weight"))
    }
    vec = cumsum(vec)
    snames = unique(unlist(snames))

    tmp = matrix(gaps, length(snames), vec[n+1], dimnames = list(snames, NULL))

    for(i in 1:n){
        nam = names(x[[i]])
        tmp[nam,(vec[i]+1):vec[i+1] ] <- as.character(x[[i]])
    }
    if (type == "DNA") 
        dat <- phyDat.DNA(tmp, return.index = TRUE)
    if (type == "AA") 
        dat <- phyDat.AA(tmp, return.index = TRUE)
    if (type == "USER") 
        dat <- phyDat.default(tmp, levels = levels, 
            return.index = TRUE)
    if (type == "CODON") 
        dat <- phyDat.codon(tmp, return.index = TRUE)            
    attr(dat,"index") <- data.frame(index=attr(dat,"index"), genes=rep(objNames, nr))   
    dat
}


write.phyDat <- function(x, file, format="phylip",...){
    if(format=="fasta") write.dna(as.character(x), file, format="fasta", colsep = "", nbcol=-1, ...)
    if(format=="phylip") write.dna(as.character(x), file, format="sequential", ...)    
    if(format=="nexus"){   
         type = attr(x, "type")
         if(type=="DNA") write.nexus.data(as.list(as.data.frame(x)), file, format = "dna",...)
         else write.nexus.data(as.list(as.data.frame(x)), file, format = "protein", ...)
         }
    }


read.phyDat <- function(file, format="phylip", type="DNA", ...){
    
    formats <- c("phylip", "nexus", "interleaved", "sequential", "fasta", "clustal")
    format <- match.arg(tolower(format), formats)
    
    if(format=="nexus") data=read.nexus.data(file, ...)
    else {
        if(format=="phylip")format="interleaved"  #"sequential"
        if (type == "DNA" || type == "CODON"){ 
            data = read.dna(file, format, as.character = TRUE, ...)
        }
        if (type == "AA") data = read.aa(file, format=format, ...)
        # raus
    }
    phyDat(data, type, return.index = TRUE)
}


baseFreq <- function(obj, freq=FALSE, all=FALSE, drop.unused.levels = FALSE){
    if (!inherits(obj,"phyDat")) 
        stop("data must be of class phyDat")
    labels <- attr(obj, "allLevels")
    weight <- attr(obj,"weight")
    n <- length(obj)    
    res <- numeric(length(labels))  
    D = diag(length(labels))   
    for(i in 1:n)res <- res + colSums(D[obj[[i]],, drop=FALSE]*weight)
    names(res) <- labels
    if(!all) res <- res[attr(obj, "levels")]
    if(!freq)res <- res/sum(res)
    if(drop.unused.levels) return(res[res>0])    
    res    
}


phylo <- function(edge, tip, edge.length=NULL){
    res <- list(edge=edge, tip.label=tip, edge.length=edge.length)
    class(res)="phylo"
    res
    }


getCols <- function (data, cols) 
{
    attrib = attributes(data)
    attr(data, "class") <- "list"
    data = data[cols]
    if (is.character(cols)) 
        attrib$names = cols
    else attrib$names = attrib$names[cols]
    attributes(data) = attrib
    attr(data, "class") <- "phyDat" 
    data
}


# allows negative indexing subset(dat,,-c(3:5))
getRows <- function (data, rows, site.pattern = TRUE) 
{   
  index <- attr(data, "index")
  if(is.data.frame(index))index = index[,1]
  if(!site.pattern){ # & all(rows>0)
    weight = tabulate(index[rows])
    ind = which(weight>0)
    rows = ind   # rows[ind]
    weight = weight[ind]
  } 
  for (i in 1:length(data)){ 
    if(is.matrix(data[[i]]))data[[i]] = data[[i]][rows,]
    else data[[i]] = data[[i]][rows]
  }  
  attr(data, "weight") = attr(data, "weight")[rows]
  if(!site.pattern) attr(data, "weight") = weight    
  attr(data, "nr") = length(attr(data, "weight"))
  attr(data, "index") = NULL
  data
}


subset.phyDat <- function (x, subset, select, site.pattern = TRUE,...) 
{  
     
    if (!missing(subset)) x <- getCols(x, subset)
    if (!missing(select)){
#         if(!site.pattern){
#             if(is.data.frame(attr(x, "index"))) select <- attr(x, "index")[select,1]
#             else select <- attr(x, "index")[select]
#         }     
         if(any(is.na(select))) return(NULL) 
         x <- getRows(x, select, site.pattern=site.pattern)
    }    
    x 
}


duplicated_phyDat <- function(x, ...){
    dm <- as.matrix(dist.hamming(x))
    diag(dm) <- 1
    res <- logical(nrow(dm))
    if(all(dm>0)) return(res)
    tmp <- which(dm==0, arr.ind = TRUE, useNames = FALSE)
    tmp <- tmp[tmp[,1] < tmp[,2],2]
    res[unique(tmp)]=TRUE
    res
}


duplicated_map2 <- function(x, ...){
    dm <- dist.hamming(x)
    if(all(dm>0)) return(NULL)
    dm <- as.matrix(dm)
    tmp <- which(dm==0, arr.ind = TRUE, useNames = FALSE)
    tmp <- tmp[tmp[,1] < tmp[,2],] 
    tmp[] = names(x)[tmp]
    g <- graph_from_edgelist(tmp, directed = FALSE)
    clu <- components(g)
    gr <- groups(clu)
    mapping <- matrix(NA, 1,2)
    for(i in 1:length(gr)){
        if(length(gr[[i]])==2)mapping <- rbind(mapping, gr[[i]])
        else{
            tmplab = gr[[i]]
            sg <- sum(dm[tmplab, tmplab])
            if(sg==0)mapping <- rbind(mapping, cbind(tmplab[1], tmplab[-1]))
        }
    }
    mapping[-1,c(2,1)]
}


duplicated_map <- function(x, ...){
# exact matches    
   dup <-  duplicated(x)
   mapping <- NULL
   if(any(dup)){ # && optNNI
       labels <- names(x)
       ind <- match(subset(x, dup), x)
       mapping <- cbind(labels[dup], labels[ind])
   }
   mapping
}   


unique.phyDat <- function(x, incomparables=FALSE, identical=TRUE, ...){
    if(identical) return(getCols(x, !duplicated(x)))
    getCols(x, !duplicated_phyDat(x))
} 


removeUndeterminedSites <- function(x, use.contrast=TRUE, undetermined=c("?", "n", "-"), ...){
    nc <- attr(x, "nc")
    nr <- attr(x, "nr")
    contrast <- attr(x, "contrast")
    if(use.contrast) ind <- which( (contrast %*% rep(1, nc)) == nc )
    else ind <- sort(match(undetermined, attr(x, "allLevels")))
    tmp <- x[[1]] %in% ind
    for(i in 2:length(x)) tmp = tmp & (x[[i]] %in% ind)
    if(any(tmp)) x <- getRows(x, (1:nr)[!tmp]) #getRows(x, -which(tmp))
    x
}


removeParsUninfoSites <- function(data){
    nr <- attr(data, "nr")
    pis <- parsinfo(data)
    if (length(pis) > 0) 
        data <- getRows(data, c(1:nr)[-pis[, 1]], TRUE)
}



allSitePattern <- function(n,levels=c("a","c","g","t"), names=NULL){
    l=length(levels)
    X=vector("list", n)
    if(is.null(names))names(X) = paste("t",1:n, sep="") 
    else names(X)=names
    for(i in 1:n)
        X[[i]] = rep(rep(levels, each=l^(i-1)),l^(n-i)) 
    X = as.data.frame(X)
    phyDat.default(X, levels, compress=FALSE) 
} 


constSitePattern <- function(n,levels=c("a","c","g","t"), names=NULL){
    l=length(levels)
    X=matrix(0, l,n)
    X = matrix(rep(levels, each=n), n, l)
    if(is.null(names))rownames(X) = paste("t",1:n, sep="")
    else rownames(X)=names
    phyDat.default(X, levels)
} 


write.phylip <- function(data, weight, file=""){
        n = sum(weight)
        m = dim(data)[2]
        cat(m,n,"\n",file = file)
        for(i in 1:m)
        cat(colnames(data)[i],"   ",toupper(rep(data[,i],weight)),"\n", sep="", file=file, append=TRUE)
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


# throw out
read.aa <- function (file, format = "interleaved", skip = 0, nlines = 0, 
    comment.char = "#", seq.names = NULL) 
{
    getTaxaNames <- function(x) {
        x <- sub("^ +", "", x)
        x <- sub(" +$", "", x)
        x <- sub("^['\"]", "", x)
        x <- sub("['\"]$", "", x)
        x
    }
    format <- match.arg(format, c("interleaved", "sequential", "fasta"))
    phylip <- if (format %in% c("interleaved", "sequential")) 
        TRUE
    else FALSE
    
    
    if (format == "fasta") {
        obj <- read.FASTA.AA(file)
        return(obj)
    }
    X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE, 
        skip = skip, nlines = nlines, comment.char = comment.char)      
           
    if (phylip) {
        fl <- X[1]
        oop <- options(warn = -1)
        fl.num <- as.numeric(unlist(strsplit(gsub("^ +", "", fl), " +")))
        options(oop)
        if (all(is.na(fl.num))) 
            stop("the first line of the file must contain the dimensions of the data")
        if (length(fl.num) != 2) 
            stop("the first line of the file must contain TWO numbers")
        else {
            n <- fl.num[1]
            s <- fl.num[2]
        }
        X <- X[-1]
        obj <- vector("character", n * s)
        dim(obj) <- c(n, s)
    }
    if (format == "interleaved") {
        fl <- X[1]
        fl <- unlist(strsplit(fl, NULL))
        bases <- grep("[-AaRrNnDdCcQqEeGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzXx?]", fl)        
        z <- diff(bases)
        for (i in 1:length(z)) if (all(z[i:(i + 8)] == 1)) 
            break
        start.seq <- bases[i]
        if (is.null(seq.names)) 
            seq.names <- getTaxaNames(substr(X[1:n], 1, start.seq - 1))
        X[1:n] <- substr(X[1:n], start.seq, nchar(X[1:n]))
        X <- gsub(" ", "", X)
        nl <- length(X)
        for (i in 1:n) obj[i, ] <- unlist(strsplit(X[seq(i, nl, n)], NULL))
    }
    if (format == "sequential") {
        fl <- X[1]
        taxa <- character(n)
        j <- 1
        for (i in 1:n) {
            bases <- grep("[-AaRrNnDdCcQqEeGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzXx?]", 
                unlist(strsplit(X[j], NULL)))
            z <- diff(bases)
            for (k in 1:length(z)) if (all(z[k:(k + 8)] == 1)) 
                break
            start.seq <- bases[k]
            taxa[i] <- substr(X[j], 1, start.seq - 1)
            sequ <- substr(X[j], start.seq, nchar(X[j]))
            sequ <- gsub(" ", "", sequ)
            j <- j + 1
            while (nchar(sequ) < s) {
                sequ <- paste(sequ, gsub(" ", "", X[j]), sep = "")
                j <- j + 1
            }
            obj[i, ] <- unlist(strsplit(sequ, NULL))
        }
        if (is.null(seq.names)) 
            seq.names <- getTaxaNames(taxa)
    }
    if (format == "fasta") return(read.FASTA.AA(file))
#        start <- grep("^ {0,}>", X)
#        taxa <- X[start]
#        n <- length(taxa)
#        obj <- vector("list", n)
#        if (is.null(seq.names)) {
#            taxa <- sub("^ {0,}>", "", taxa)
#            seq.names <- getTaxaNames(taxa)
#        }
#        start <- c(start, length(X) + 1)
#        for (i in 1:n) obj[[i]] <- unlist(strsplit(gsub(" ", 
#            "", X[(start[i] + 1):(start[i + 1] - 1)]), NULL))
#    }
    if (phylip) {
        rownames(obj) <- seq.names
        obj <- tolower(obj)
    }
    else {
        names(obj) <- seq.names
        obj <- lapply(obj, tolower)
    }
    obj   
}


genlight2phyDat <- function(x, ambiguity=NA){
    tmp <- as.matrix(x)
    lev <- na.omit(unique(as.vector(tmp)))
    phyDat(tmp, "USER", levels=lev, ambiguity=ambiguity)
}


if (getRversion() >= "2.15.1") utils::globalVariables("image.AAbin")
image.phyDat <- function(x, ...){
    if(attr(x, "type")=="AA")image(as.AAbin(x), ...)
    if(attr(x, "type")=="DNA")image(as.DNAbin(x), ...)
    else return(NULL)
}

