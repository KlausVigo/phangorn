#' Generic functions for class phyDat
#'
#' These functions help to manipulate alignments of class phyDat.
#'
#'
#' \code{allSitePattern} generates all possible site patterns and can be useful
#' in simulation studies. For further details see the vignette
#' AdvancedFeatures.
#'
#' The generic function \code{c} can be used to to combine sequences and
#' \code{unique} to get all unique sequences or unique haplotypes.
#'
#' \code{phyDat} stores identical columns of an alignment only once and keeps an
#' index of the original positions. This saves memory and especially
#' computations as these are usually need to be done only once for each site
#' pattern.
#' In the example below the matrix x in the example has 8 columns, but column 1
#' and 2 and also 3 and 5 are identical. The \code{phyDat} object y has only 6
#' site pattern. If argument \code{site.pattern=FALSE} the indexing behaves like
#' on the original matrix x. \code{site.pattern=TRUE} can be useful inside
#' functions.
#'
#' @aliases
#' cbind.phyDat c.phyDat removeUndeterminedSites
#' @param x An object containing sequences.
#' @param levels Level attributes.
#' @param n Number of sequences.
#' @param names Names of sequences.
#' @param subset a subset of taxa.
#' @param select a subset of characters.
#' @param site.pattern select site pattern or sites (see details).
#' @param incomparables for compatibility with unique.
#' @param identical if TRUE (default) sequences have to be identical, if FALSE
#' sequences are considered duplicates if distance between sequences is zero
#' (happens frequently with ambiguous sites).
#' @param type Type of sequences ("DNA", "AA" or "USER").
#' @param ... further arguments passed to or from other methods.
#' @return The functions return an object of class \code{phyDat}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{DNAbin}}, \code{\link[ape]{as.DNAbin}},
#' \code{\link{baseFreq}}, \code{\link{glance.phyDat}}, \code{\link{dna2codon}},
#' \code{\link[ape]{read.dna}}, \code{\link[ape]{read.nexus.data}}
#' and the chapter 1 in the \code{vignette("AdvancedFeatures",
#' package="phangorn")} and the example of \code{\link{pmlMix}} for the use of
#' \code{\link{allSitePattern}}.
#' @keywords cluster
#' @examples
#'
#' data(Laurasiatherian)
#' class(Laurasiatherian)
#' Laurasiatherian
#' # base frequencies
#' baseFreq(Laurasiatherian)
#' # subsetting phyDat objects
#' # the first 5 sequences
#' subset(Laurasiatherian, subset=1:5)
#' # the first 5 characters
#' subset(Laurasiatherian, select=1:5, site.pattern = FALSE)
#' # subsetting with []
#' Laurasiatherian[1:5, 1:20]
#' # short for
#' subset(Laurasiatherian, subset=1:5, select=1:20, site.pattern = FALSE)
#' # the first 5 site patterns (often more than 5 characters)
#' subset(Laurasiatherian, select=1:5, site.pattern = TRUE)
#'
#' x <- matrix(c("a", "a", "c", "g", "c", "t", "a", "g",
#'               "a", "a", "c", "g", "c", "t", "a", "g",
#'               "a", "a", "c", "c", "c", "t", "t", "g"), nrow=3, byrow = TRUE,
#'             dimnames = list(c("t1", "t2", "t3"), 1:8))
#' (y <- phyDat(x))
#'
#' subset(y, 1:2)
#' subset(y, 1:2, compress=TRUE)
#'
#' subset(y, select=1:3, site.pattern = FALSE) |> as.character()
#' subset(y, select=1:3, site.pattern = TRUE) |> as.character()
#' y[,1:3] # same as subset(y, select=1:3, site.pattern = FALSE)
#'
#' # Compute all possible site patterns
#' # for nucleotides there $4 ^ (number of tips)$ patterns
#' allSitePattern(5)
#'
#' @rdname phyDat
#' @export
print.phyDat <- function (x, ...){
  cat(length(x), "sequences with",sum(attr(x,"weight")), "character and",
      attr(x,"nr"),"different site patterns.\n")
  cat("The states are",attr(x,"levels"), "\n")
}


#' @export cbind.phyDat
#' @export
cbind.phyDat <- function(..., gaps="-", compress=TRUE){
  object <- as.list(substitute(list(...)))[-1]
  x <- list(...)
  n <- length(x)
  if (n == 1) return(x[[1]])

  types <- sapply(x, function(x)attr(x, "type"))
  if(any(types!=types[1]))stop("Alignments must have same type!")
  nr <- numeric(n)
  ATTR <- attributes(x[[1]])
  nr[1] <- sum(attr(x[[1]], "weight"))
  allLevels <- attr(x[[1]], "allLevels")
  gapsInd <- match(gaps, allLevels)
  snames <- vector("list", n)
  vec <- numeric(n+1)
  wvec <- numeric(n+1)
  objNames <- as.character(object)
  if(any(duplicated(objNames))) objNames <- paste0(objNames, 1:n)

  for(i in 1:n){
    snames[[i]] <- names(x[[i]])
    nr[i] <- sum(attr(x[[i]], "weight"))
    vec[i+1] <- attr(x[[i]], "nr")
    wvec[i+1] <- sum(attr(x[[i]], "weight"))
  }
  vec <- cumsum(vec)
  wvec <- cumsum(wvec)
  snames <- unique(unlist(snames))
  weight <- numeric(vec[n+1])
  index <- numeric(wvec[n+1])
  ATTR$names <- snames
  ATTR$nr <- vec[n+1]
  tmp <- matrix(gapsInd, vec[n+1], length(snames),
                dimnames = list(NULL, snames))
  tmp <- as.data.frame(tmp)
  add.index <- TRUE
  for(i in 1:n){
    nam <- names(x[[i]])
    tmp[(vec[i]+1):vec[i+1], nam] <- x[[i]][nam]
    weight[(vec[i]+1):vec[i+1]] <- attr(x[[i]], "weight")
  }
  if(compress){
    dindex <- grp_duplicated(as.matrix(tmp), MARGIN=1)
    attr(dindex, "nlevels") <- NULL
    weight <- aggregate(weight, by=list(dindex), FUN=sum)$x
    pos <- which(!duplicated(dindex))
    tmp <- tmp[pos,]
  }
  if(any(sapply(x, function(x)is.null(attr(x, "index"))))) add.index <- FALSE
  if(add.index  & compress){
    for(i in 1:n){
      tmp2 <- attr(x[[i]], "index")
      if(!is.null(tmp2)){
        if(is.data.frame(tmp2))index[(wvec[i]+1):wvec[i+1]] <-
            dindex[(vec[i]+1):vec[i+1]][tmp2[,1]]
        else index[(wvec[i]+1):wvec[i+1]] <-
            dindex[(vec[i]+1):vec[i+1]][tmp2]
      }
      else add.index <- FALSE
    }
  }
  if(add.index) ATTR$index <- data.frame(index = index, genes=rep(objNames, nr))
  ATTR$weight <- weight
  ATTR$nr <- length(weight)
  attributes(tmp) <- ATTR
  tmp
}


# @rdname phyDat
#' @export c.phyDat
#' @export
rbind.phyDat <- function(...){
  x <- list(...)
  types <- sapply(x, function(x)attr(x, "type"))
  l <- sapply(x, function(x)sum(attr(x, "weight")))
  if(any(l!=l[1]))stop("Alignments have different # of characters!")
  if(any(types!=types[1]))stop("Alignments must have same type!")
  nam <- lapply(x, names) |> unlist()
  if(any(duplicated(nam)))stop("Duplicated names!")
  m <- lengths(x)
  mcs <- c(0, cumsum(m))
  res <- matrix(NA_character_, sum(m), l[1], dimnames=list(nam, NULL))
  for(i in seq_along(x)){
    res[(mcs[i]+1):mcs[i+1], ] <- as.character(x[[i]])
  }
  if(types[1]=="USER"){
    contrast <- attr(x[[1]], "contrast")
    dimnames(contrast) <- list(attr(x[[1]], "allLevels"),
                               attr(x[[1]], "levels"))
    return(phyDat(res, type="USER", contrast=contrast))
  }
  phyDat(res, type=types[1])
}


# @rdname phyDat
#' @export c.phyDat
#' @export
c.phyDat <- cbind.phyDat


compress.phyDat <- function(data){
  attrib <- attributes(data)
  attr(data, "class") <- "list"
  index <- grp_duplicated( matrix(unlist(data, use.names = FALSE), attrib$nr,
                                  length(data)))
  attrib$nr <- attr(index, "nlevels")
  attr(index, "nlevels") <- NULL
  pos <- which(!duplicated(index))
  weight <- tapply(attrib$weight, index, sum)
  names(weight) <- NULL
  attrib$weight <- as.vector(weight)
  if(is.null(attrib$index)) attrib$index <-index
  else attrib$index <- index[attrib$index]
  for(i in seq_len(length(data))) data[[i]] <- data[[i]][pos]
  attributes(data) <- attrib
  attr(data, "class") <- "phyDat"
  data
}



uncompress.phyDat <- function(data){
  attrib <- attributes(data)
  stopifnot(all.equal(sum(attrib$weight), length(attrib$index)))
  for(i in seq_along(data))
     data[[i]] <- data[[i]][attrib$index]
  attrib$nr <- length(attrib$index)
  attrib$index <- seq_along(attrib$index)
  attrib$weight <- rep(1, attrib$nr)
  attributes(data) <- attrib
  data
}


getCols <- function (data, cols, compress=FALSE){
  attrib <- attributes(data)
  if(inherits(attr(data, "index"), "data.frame")) compress <- FALSE
  attr(data, "class") <- "list"
  data <- data[cols]
  if (is.character(cols))
    attrib$names <- cols
  else attrib$names <- attrib$names[cols]
  attributes(data) <- attrib
  attr(data, "class") <- "phyDat"
  if(compress) return(compress.phyDat(data))
  data
}


getRows <- function (data, rows, site.pattern = TRUE){
  index <- attr(data, "index")
  if(is.data.frame(index))index <- index[,1]
  if(!site.pattern){
    if(is.null(index)) index <- seq_len(length(data[[1]]))
    weight <- tabulate(index[rows])
    ind <- which(weight>0)
# update index
    new_index <- integer(length(weight))
    new_index[ind] <- seq_along(ind)
    attr(data, "index") <- new_index[index[rows]]
    rows <- ind   # rows[ind]
    weight <- weight[ind]
  }
  for (i in seq_along(data)){
    if(is.matrix(data[[i]]))data[[i]] <- data[[i]][rows, ]
    else data[[i]] <- data[[i]][rows]
  }
  attr(data, "weight") <- attr(data, "weight")[rows]
  if(!site.pattern) attr(data, "weight") <- weight
  attr(data, "nr") <- length(attr(data, "weight"))
  if(site.pattern)attr(data, "index") <- NULL
  data
}


#' @rdname phyDat
#' @method subset phyDat
#' @export
subset.phyDat <- function (x, subset, select, site.pattern = TRUE, ...){
  if (!missing(subset)){
    if(is.numeric(subset) & any(subset>length(x)))
        stop("subscript out of bounds")
    x <- getCols(x, subset, ...)
  }
  if (!missing(select)){
    w <- attr(x, "weight")
    if(site.pattern){
      if(any(select > length(w))) stop("subscript out of bounds")
    }
    else if(any(select > sum(w))) stop("subscript out of bounds")
    if(any(is.na(select))) return(NULL)
    x <- getRows(x, select, site.pattern=site.pattern)
  }
  x
}


#' @param i,j	indices of the rows and/or columns to select or to drop. They
#' may be numeric, logical, or character (in the same way than for standard R
#' objects).
#' @param drop for compatibility with the generic (unused).
#' @rdname phyDat
#' @export
"[.phyDat" <- function(x, i, j, ..., drop=FALSE){
   subset(x, subset = i, select = j, site.pattern=FALSE, compress=FALSE)
}


#' @rdname phangorn-internal
#' @export
map_duplicates <-  function(x, dist=length(x)<500, ...){
  labels <- names(x)
  nr <- attr(x, "nr")

  if(dist){
    y <- as.matrix(dist.hamming(x, FALSE))
    l <- nrow(y)
    ind <- duplicated(y)
  }
  else ind <- duplicated(x)
  res <- NULL
  if(any(ind)){
    ind2 <- grp_duplicated( matrix(unlist(x, recursive = FALSE,
                            use.names = FALSE), nr, length(labels)), MARGIN=2)
    if(dist) ind2 <- grp_duplicated(y)
    ind2 <- ind2[ind]
    res <- data.frame(duplicates=labels[ind], where=labels[!ind][ind2],
                      stringsAsFactors = FALSE)
  }
  res
}


#' @rdname phyDat
#' @method unique phyDat
#' @export
unique.phyDat <- function(x, incomparables=FALSE, identical=TRUE, ...){
  tmp <- map_duplicates(x, !identical)
  if (!is.null(tmp)) {
    x <- getCols(x, setdiff(names(x), tmp[, 1]))
    attr(x, "duplicated") <- list(tmp)
  }
  x
}


addTaxa <- function(tree, dup_list) {
  fun <- function(x, dup_list){
    if(Ntip(x)==1L){
      dup <- dup_list[[1]][,1]
      x <- add.tips(x, dup, rep(2, length(dup)))
      dup_list[[1]] <- NULL
    }
    for (i in seq_along(dup_list)) {
      dup <- dup_list[[i]]
      x <- add.tips(x, dup[, 1], dup[, 2])
    }
    x
  }
  if(inherits(tree, "phylo")) return(fun(tree, dup_list))
  if(inherits(tree, "multiPhylo")){
    trees <- .uncompressTipLabel(tree)
    trees <- unclass(trees)
    trees <- lapply(trees, fun, dup_list)
    class(trees) <- "multiPhylo"
    trees <- .compressTipLabel(trees)
    return(trees)
  }
  NULL
}


#' @rdname phyDat
#' @export
removeUndeterminedSites <- function(x, ...){
# , use.contrast=TRUE, undetermined=c("?", "n", "-")
  nc <- attr(x, "nc")
  nr <- attr(x, "nr")
  contrast <- attr(x, "contrast")
#  if(use.contrast)
  ind <- which( (contrast %*% rep(1, nc)) == nc )
#  else ind <- sort(match(undetermined, attr(x, "allLevels")))
  tmp <- x[[1]] %in% ind
  for(i in 2:length(x)) tmp <- tmp & (x[[i]] %in% ind)
  if(any(tmp)) x <- getRows(x, (1:nr)[!tmp])
  x
}


hasAmbiguousSites <- function(x){
  contrast <- attr(x, "contrast")
  nc <- as.integer(attr(x, "nc"))
  con <- rowSums(contrast > 0) > 1
  for (i in seq_along(x)) {
    tmp <- con[x[[i]]]
    if(any(tmp)) return(TRUE)
  }
  FALSE
}


#' @rdname phyDat
#' @export
removeAmbiguousSites <- function(x){
  contrast <- attr(x, "contrast")
  nc <- as.integer(attr(x, "nc"))
  con <- rowSums(contrast > 0) < 2
  index <- con[x[[1]]]
  for (i in 2:length(x)) index <- index & con[x[[i]]]
  index <- which(index)
  if(length(index)==0) stop('each site contains at least one ambiguous state!')
  subset(x, select = index, site.pattern = TRUE)
}


removeParsUninfoSites <- function(data, exact=TRUE){
  nr <- attr(data, "nr")
  pis <- parsinfo(data, exact)
  if (length(pis) > 0){
    p0 <- sum(attr(data, "weight")[pis[, 1]] * pis[, 2])
    data <- getRows(data, c(1:nr)[-pis[, 1]], TRUE)
    if(is.null(attr(data, "informative")))
      attr(data, "informative") <- seq_len(nr)[-pis[, 1]]
    else attr(data, "informative") <- attr(data, "informative")[-pis[, 1]]
  }
  else p0 <- 0
  if(length(attr(data, "p0"))) p0 <- p0 + attr(data, "p0")
  attr(data, "p0") <- p0
  data
}


removeParsimonyUninfomativeSites <- function(data, recursive=FALSE, exact=TRUE){
  dup_list <- NULL
  tmp <- TRUE
  star_tree <- FALSE
  while (tmp) {
    nam <- names(data)
    data <- removeParsUninfoSites(data, exact)
    if(!recursive) return(data)
    if (attr(data, "nr") == 0) {
      star_tree <- TRUE
      break()
      tmp <- FALSE
    }
    # unique sequences
    dup <- map_duplicates(data)
    if (!is.null(dup)) {
      dup_list <- c(list(dup), dup_list)
      data <- subset(data, setdiff(names(data), dup[, 1]), site.pattern=TRUE)
    }
    else break() # tmp <- FALSE
  }
  attr(data, "duplicated") <- dup_list
  data
}


#' @rdname phyDat
#' @param code The ncbi genetic code number for translation.
#' By default the standard genetic code is used.
#' @export
allSitePattern <- function(n, levels=NULL, names=NULL, type="DNA", code=1){
  type <- match.arg(type, c("DNA", "AA", "CODON", "USER"))
  if(type=="DNA") levels <- c("a", "c", "g", "t")
  if(type=="AA") levels <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  if(type=="CODON"){
    tmp <- .CODON[, as.character(code)]
    levels <- rownames(.CODON)
    levels <- levels[tmp != "*"]
  }
  l <- length(levels)
  X <- matrix(NA_integer_, n, l^n)
  if(is.null(names))rownames(X) <- paste0("t", 1:n)
  else rownames(X) <- names
  for(i in 1:n)
    X[i, ] <- rep(rep(levels, each=l^(i-1)), l^(n-i))
  if(type=="DNA") return(phyDat.DNA(X, compress=FALSE, return.index=FALSE))
  if(type=="AA") return(phyDat.AA(X, return.index=FALSE))
  res <- phyDat.default(X, levels, compress=FALSE, return.index=FALSE)
  if(type=="CODON"){
    attr(res, "type") <- "CODON"
    attr(res, "code") <- code
  }
  res
}


constSitePattern <- function(n, names=NULL, type="DNA", levels=NULL){
  if(type=="DNA"){
    levels <- c("a", "c", "g", "t")
    l <- 4L
  } else if(type=="AA"){
    levels <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
                "F", "P", "S", "T", "W", "Y", "V")
    l <- 20L
  }
  else l <- length(levels)
  X <- matrix(0, l, n)
  X <- matrix(rep(levels, each=n), n, l)
  if(is.null(names)) rownames(X) <- paste0("t", seq_len(n))
  else rownames(X) <- names
  phyDat(X, type=type, levels=levels)
}



#' Read Amino Acid Sequences in a File
#'
#' This function reads amino acid sequences in a file, and returns a matrix
#' list of DNA sequences with the names of the taxa read in the file as row
#' names.
#'
#'
#' @param file a file name specified by either a variable of mode character, or
#' a double-quoted string.
#' @param format a character string specifying the format of the DNA sequences.
#' Three choices are possible: \code{"interleaved"}, \code{"sequential"}, or
#' \code{"fasta"}, or any unambiguous abbreviation of these.
#' @param skip the number of lines of the input file to skip before beginning
#' to read data.
#' @param nlines the number of lines to be read (by default the file is read
#' until its end).
#' @param comment.char a single character, the remaining of the line after this
#' character is ignored.
#' @param seq.names the names to give to each sequence; by default the names
#' read in the file are used.
#' @return a matrix of amino acid sequences.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{read.dna}}, \code{\link[ape]{read.GenBank}},
#' \code{\link[phangorn]{phyDat}}, \code{\link[seqinr]{read.alignment}}
#' @references % Anonymous. FASTA format description. %
#' \url{https://en.wikipedia.org/wiki/FASTA_format}
#'
#' Felsenstein, J. (1993) Phylip (Phylogeny Inference Package) version 3.5c.
#' Department of Genetics, University of Washington.
#' \url{https://phylipweb.github.io/phylip/}
#' @keywords IO
#' @noRd
read.aa <- function (file, format = "interleaved", skip = 0, nlines = 0,
                     comment.char = "#", seq.names = NULL){
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
    obj <- read.FASTA(file, type = "AA")
    return(obj)
  }
  X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE,
    skip = skip, nlines = nlines, comment.char = comment.char)

  if (phylip) {
    fl <- X[1]
    oop <- options(warn = -1)
    fl.num <- as.numeric(unlist(strsplit(gsub("^ +", "", fl), " +")))
    options(oop)
    if (all(is.na(fl.num))  || length(fl.num) != 2)
      stop("the first line of the file must contain the dimensions of the data")
#  if (length(fl.num) != 2)
#    stop("the first line of the file must contain the dimensions of the data")
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
    for (i in seq_along(z)) if (all(z[i:(i + 8)] == 1))
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
      for (k in seq_along(z)) if (all(z[k:(k + 8)] == 1))
        break
      start.seq <- bases[k]
      taxa[i] <- substr(X[j], 1, start.seq - 1)
      sequ <- substr(X[j], start.seq, nchar(X[j]))
      sequ <- gsub(" ", "", sequ)
      j <- j + 1
      while (nchar(sequ) < s) {
        sequ <- paste0(sequ, gsub(" ", "", X[j]))
        j <- j + 1
      }
      obj[i, ] <- unlist(strsplit(sequ, NULL))
    }
    if (is.null(seq.names))
      seq.names <- getTaxaNames(taxa)
  }
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

