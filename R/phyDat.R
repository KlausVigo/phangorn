#
# Data structures for ML and MP
#


#' Conversion among Sequence Formats
#'
#' These functions transform several DNA formats into the \code{phyDat} format.
#' \code{allSitePattern} generates an alignment of all possible site patterns.
#'
#' If \code{type} "USER" a vector has to be give to \code{levels}. For example
#' c("a", "c", "g", "t", "-") would create a data object that can be used in
#' phylogenetic analysis with gaps as fifth state.  There is a more detailed
#' example for specifying "USER" defined data formats in the vignette
#' "phangorn-specials".
#'
#' \code{allSitePattern} returns all possible site patterns and can be useful
#' in simulation studies. For further details see the vignette
#' phangorn-specials.
#'
#' The generic function \code{c} can be used to to combine sequences and
#' \code{unique} to get all unique sequences or unique haplotypes.
#'
#' \code{acgt2ry} converts a \code{phyDat} object of nucleotides into an binary
#' ry-coded dataset.
#'
#' @aliases
#' as.phyDat.character as.phyDat.data.frame as.phyDat.matrix
#' as.MultipleAlignment as.MultipleAlignment.phyDat cbind.phyDat c.phyDat
#' acgt2ry removeUndeterminedSites phyDat2MultipleAlignment
#' @param data An object containing sequences.
#' @param x An object containing sequences.
#' @param type Type of sequences ("DNA", "AA", "CODON" or "USER").
#' @param levels Level attributes.
#' @param return.index If TRUE returns a index of the site patterns.
#' @param n Number of sequences.
#' @param names Names of sequences.
#' @param subset a subset of taxa.
#' @param select a subset of characters.
#' @param site.pattern select site pattern or sites.
#' @param allLevels return original data.
#' @param obj as object of class phyDat
#' @param freq logical, if 'TRUE', frequencies or counts are returned otherwise
#' proportions
#' @param all all a logical; if all = TRUE, all counts of bases, ambiguous
#' codes, missing data, and alignment gaps are returned as defined in the
#' contrast.
#' @param drop.unused.levels logical, drop unused levels
#' @param incomparables for compatibility with unique.
#' @param identical if TRUE (default) sequences have to be identical, if FALSE
#' sequences are considered duplicates if distance between sequences is zero
#' (happens frequently with ambiguous sites).
#' @param ambiguity character for ambiguous character and no contrast is
#' provided.
#' @param ... further arguments passed to or from other methods.
#' @return The functions return an object of class \code{phyDat}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{DNAbin}}, \code{\link{as.DNAbin}},
#' \code{\link{read.dna}}, \code{\link{read.aa}}, \code{\link{read.nexus.data}}
#' and the chapter 1 in the \code{vignette("phangorn-specials",
#' package="phangorn")} and the example of \code{\link{pmlMix}} for the use of
#' \code{allSitePattern}
#' @keywords cluster
#' @examples
#'
#' data(Laurasiatherian)
#' class(Laurasiatherian)
#' Laurasiatherian
#' # base frequencies
#' baseFreq(Laurasiatherian)
#' baseFreq(Laurasiatherian, all=TRUE)
#' baseFreq(Laurasiatherian, freq=TRUE)
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
#' # transform into old ape format
#' LauraChar <- as.character(Laurasiatherian)
#' # and back
#' Laura <- phyDat(LauraChar)
#' all.equal(Laurasiatherian, Laura)
#' # Compute all possible site patterns
#' # for nucleotides there $4 ^ (number of tips)$ patterns
#' allSitePattern(5)
#'
#' @rdname phyDat
#' @export
phyDat <- function (data, type="DNA", levels=NULL, return.index = TRUE, ...){
  if (inherits(data, "DNAbin")) type <- "DNA"
  pt <- match.arg(type, c("DNA", "AA", "CODON", "USER"))
  if(pt=="DNA") dat <- phyDat.DNA(data, return.index=return.index, ...)
  if(pt=="AA") dat <- phyDat.AA(data, return.index=return.index, ...)
  if(pt=="CODON") dat <- phyDat.codon(data, return.index=return.index, ...)
  if(pt=="USER") dat <- phyDat.default(data, levels = levels,
                                       return.index=return.index, ...)
  dat
}


#' @rdname phyDat
#' @export
as.phyDat <- function (x, ...){
  if (inherits(x, "phyDat")) return(x)
  UseMethod("as.phyDat")
}


#' @rdname phyDat
#' @method as.phyDat factor
#' @export
as.phyDat.factor <- function(x, ...){
  nam <- names(x)
  lev <- levels(x)
  x <- as.character(x)
  names(x) <- nam
  phyDat(x, type="USER", levels = lev, ...)
}


#' @rdname phyDat
#' @method as.phyDat DNAbin
#' @export
as.phyDat.DNAbin <- function(x, ...) phyDat.DNA(x, ...)


#' @rdname phyDat
#' @method as.phyDat alignment
#' @export
as.phyDat.alignment <- function (x, type="DNA", ...){
  x$seq <- tolower(x$seq)
  data <- sapply(x$seq, strsplit, "")
  names(data) <- x$nam
  if(type=="DNA") dat <- phyDat.DNA(data, ...)
  if(type=="AA") dat <- phyDat.AA(data, ...)
  if(type=="CODON") dat <- phyDat.codon(data, ...)
  if(type=="USER") dat <- phyDat.default(data, ...)
  dat
}


#as.alignment.phyDat <- function(x, ...) as.alignment(as.character(x))
#' @rdname phyDat
#' @export
phyDat2alignment <-  function(x){
  z <- as.character(x)
  nam <- rownames(z)
  type <- attr(x, "type")
  seq <- switch(type,
                DNA = tolower(apply(z, 1, paste, collapse="")),
                AA = toupper(apply(z, 1, paste, collapse="")))
  names(seq) <- NULL
  res <- list(nb=length(seq), nam=nam, seq=seq, com=NA)
  class(res) <- "alignment"
  res
}


#' @rdname phyDat
#' @method as.phyDat MultipleAlignment
#' @export
as.phyDat.MultipleAlignment <- function(x, ...){
  if (requireNamespace("Biostrings")){
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


# @rdname phyDat
#' @export
as.MultipleAlignment <- function (x, ...){
  if (inherits(x, "MultipleAlignment")) return(x)
  UseMethod("as.MultipleAlignment")
}


#' @rdname phyDat
#' @export
as.MultipleAlignment.phyDat <- function(x, ...){
  if (requireNamespace("Biostrings")){
  z <- as.character(x)
#  nam <- rownames(z)
  type <- attr(x, "type")
  seq <- switch(type,
                DNA = tolower(apply(z, 1, paste, collapse="")),
                AA = toupper(apply(z, 1, paste, collapse="")))
  if(type=="DNA") return(Biostrings::DNAMultipleAlignment(seq))
  if(type=="AA") return(Biostrings::AAMultipleAlignment(seq))
  }
  return(NULL)
}


#' @export
phyDat2MultipleAlignment <- as.MultipleAlignment.phyDat


#' @export
as.phyDat.matrix <- function (x, ...) phyDat(data=x, ...)


#' @export
as.phyDat.character <- function (x, ...) phyDat(data=x, ...)

# @rdname phyDat
#' @export
as.phyDat.data.frame <- function (x, ...) phyDat(data=x, ...)


#' @rdname phyDat
#' @export
acgt2ry <- function(obj){
   ac <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y",
        "k", "v", "h", "d", "b", "n", "?", "-")
   AC <- matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
        0, 1, 1, 1), c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1,
        0, 1, 1, 1, 1), c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1,
        0, 1, 1, 1, 1, 1), c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1,
        0, 1, 1, 1, 1, 1, 1)), 18, 4, dimnames = list(NULL, c("a",
        "c", "g", "t")))
   ry <- AC[c(7, 10), ]
   RY <- AC %*% t(ry)
   RY[RY==2] <- 1
   dimnames(RY) <- list(NULL, c("r", "y"))
   attr(obj, "levels") <- c("r", "y")
   attr(obj, "nc") <- 2
   attr(obj, "type") <- "USER"
   attr(obj, "contrast") <- RY
   obj <- phyDat.default(as.character(obj, allLevels=FALSE),
                         levels = c("r", "y"), ambiguity = NULL)
   obj
}


#' @rdname phyDat
#' @export
as.character.phyDat <- function (x, allLevels=TRUE, ...){
  nr <- attr(x, "nr")
#  nc <- attr(x, "nc")
  type <- attr(x, "type")
  labels <- attr(x, "allLevels")

  if (!is.null(attr(x, "index"))) {
    index <- attr(x, "index")
    if (is.data.frame(index))
      index <- index[, 1]
  }
  else index <- rep(1:nr, attr(x, "weight"))
  if (type == "USER") {
    #levels in acgt2ry
    if(!allLevels){
      tmp <- attr(x, "levels")
      contrast <- attr(x, "contrast")
      contrast[contrast>0] <- 1
      ind <- which(rowSums(contrast)==1)
      contrast[rowSums(contrast)>1, ] <- 0
      labels <- rep(NA, length(attr(x, "allLevels")))
      labels[ind] <- tmp[contrast%*%c(seq_along(tmp))]
    }
  }
  if(type == "AA") labels <- toupper(labels)
  if(type == "CODON"){
    nr <- length(index)
    result <- matrix(NA, nrow = length(x), ncol = 3L*nr)
    labels <- strsplit(labels, "")
    for (i in seq_along(x)) result[i, ] <- unlist(labels[ x[[i]][index] ])
  }
  else {
    result <- matrix(NA, nrow = length(x), ncol = nr)
    for (i in seq_along(x)) result[i, ] <- labels[x[[i]]]
    result <- result[, index, drop = FALSE]
  }
  rownames(result) <- names(x)
  result
}


#' @rdname phyDat
#' @export
as.data.frame.phyDat <- function(x, ...){
  nr <- attr(x, "nr")
#  nc <- attr(x, "nc")
  labels <- attr(x, "allLevels")
  if(attr(x, "type") == "AA") labels <- toupper(labels)
  result <- vector("list", length(x))
  if (is.null(attr(x, "index")))
    index <- rep(1:nr, attr(x, "weight"))
  else {
    index <- attr(x, "index")
    if (is.data.frame(index))
      index <- index[, 1]
  }
  for (i in seq_along(x)) result[[i]] <- labels[x[[i]][index]]
  attr(result, "names") <- names(x)
  attr(result, "row.names") <- seq_along(index)
  attr(result, "class") <- "data.frame"
  result
}



#' @rdname phyDat
#' @export
as.DNAbin.phyDat <- function (x, ...){
  if(attr(x, "type")=="DNA"){
  nr <- attr(x, "nr")
  ac <- attr(x, "allLevels")
  result <- matrix(as.raw(0), nrow = length(x), ncol = nr)
  # from ape ._cs_
  cs <- c("a", "g", "c", "t", "r", "m", "w", "s", "k", "y", "v", "h",
    "d", "b", "n", "-", "?")
  # from ape ._bs_
  bs <- as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48, 224, 176, 208,
                 112, 240, 4, 2))
  ord <- match(ac, cs)
  ord[5] <- 4

  for (i in seq_along(x)){
    ind <- ord[x[[i]]]
    result[i,] <- bs[ind]
  }
  if (is.null(attr(x, "index")))
    index <- rep(1:nr, attr(x, "weight"))
  else {
    index <- attr(x, "index")
    if (is.data.frame(index))
      index <- index[, 1]
  }
  result <- result[, index, drop = FALSE]
  rownames(result) <- names(x)
  class(result) <- "DNAbin"
  return(result)
  }
  else stop("x must be a nucleotide sequence")
}


#' @rdname phyDat
#' @export
as.AAbin.phyDat <- function(x,...) {
   if(attr(x, "type")=="AA") return(as.AAbin(as.character(x, ...)))
   else stop("x must be a amino acid sequence")
}


#' @export
print.phyDat <- function (x, ...){
  cat(length(x), "sequences with",sum(attr(x,"weight")), "character and",
    attr(x,"nr"),"different site patterns.\n")
  cat("The states are",attr(x,"levels"), "\n")
}


## @export
summary.phyDat <- function (x, ...){
  nc  <- attr(x, "nc")
  nseq <- length(x)
  nchar <- sum(attr(x,"weight"))
  unique_sites <- attr(x, "nr")
  tmp <- logical(unique_sites)
  for(i in 2:nseq) tmp <- tmp | (x[[1]] != x[[i]])
  const_sites <- sum(attr(x, "weight")[tmp==0])
  list(nseq=nseq, nchar=nchar, unique_sites=unique_sites, const_sites=const_sites)
}


#' @export cbind.phyDat
#' @export
cbind.phyDat <- function(..., gaps="-", compress=TRUE){
  object <- as.list(substitute(list(...)))[-1]
  x <- list(...)
  n <- length(x)
  if (n == 1) return(x[[1]])
  #  type <- attr(x[[1]], "type")
  nr <- numeric(n)
  ATTR <- attributes(x[[1]])
  nr[1] <- sum(attr(x[[1]], "weight"))
  #  levels <- attr(x[[1]], "levels")
  allLevels <- attr(x[[1]], "allLevels")
  gapsInd <- match(gaps, allLevels)
  snames <- vector("list", n)
  vec <- numeric(n+1)
  wvec <- numeric(n+1)
  objNames <- as.character(object)
  if(any(duplicated(objNames))) objNames <- paste0(objNames, 1:n)
  #  tmp <- as.character(x[[1]])

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
c.phyDat <- cbind.phyDat


#' @rdname phyDat
#' @export
baseFreq <- function(obj, freq=FALSE, all=FALSE, drop.unused.levels = FALSE){
  if (!inherits(obj,"phyDat"))
    stop("data must be of class phyDat")
  labels <- attr(obj, "allLevels")
  weight <- attr(obj,"weight")
  n <- length(obj)
  res <- numeric(length(labels))
  D <- diag(length(labels))
  for(i in 1:n)res <- res + colSums(D[obj[[i]],, drop=FALSE]*weight)
  names(res) <- labels
  if(!all) res <- res[attr(obj, "levels")]
  if(!freq)res <- res/sum(res)
  if(drop.unused.levels) return(res[res>0])
  res
}


phylo <- function(edge, tip, edge.length=NULL){
  res <- list(edge=edge, tip.label=tip, edge.length=edge.length)
  class(res) <- "phylo"
  res
  }

compress.phyDat <- function(data){
  attrib <- attributes(data)
  attr(data, "class") <- "list"
  index <- grp_duplicated( matrix(unlist(data, use.names = FALSE), attrib$nr, length(data)))
  attrib$nr <- attr(index, "nlevels")
  attr(index, "nlevels") <- NULL
  pos <- which(!duplicated(index))
  attrib$weight <- tapply(attrib$weight, index, sum)
  if(is.null(attrib$index)) attrib$index <-index
  else attrib$index <- index[attrib$index]
  for(i in seq_len(length(data))) data[[i]] <- data[[i]][pos]
  attributes(data) <- attrib
  attr(data, "class") <- "phyDat"
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
  if(compress){
    index <- grp_duplicated( matrix(unlist(data, use.names = FALSE), attrib$nr, length(data)))
    attrib$nr <- attr(index, "nlevels")
    attr(index, "nlevels") <- NULL
    pos <- which(!duplicated(index))
    attrib$weight <- tapply(attrib$weight, index, sum)
    if(is.null(attrib$index)) attrib$index <-index
    else attrib$index <- index[attrib$index]
    for(i in seq_len(length(data))) data[[i]] <- data[[i]][pos]
  }
  attributes(data) <- attrib
  attr(data, "class") <- "phyDat"
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
    if(is.numeric(subset) & any(subset>length(x))) stop("subscript out of bounds")
    x <- getCols(x, subset)
  }
  if (!missing(select)){
    w <- attr(x, "weight")
    if(site.pattern) if(any(select > length(w))) stop("subscript out of bounds")
    else if(any(select > sum(w))) stop("subscript out of bounds")
    if(any(is.na(select))) return(NULL)
    x <- getRows(x, select, site.pattern=site.pattern)
  }
  x
}


## Needs testing that it is not used e.g. prepareDataFitch returns no class
#' @param i,j	indices of the rows and/or columns to select or to drop. They
#' may be numeric, logical, or character (in the same way than for standard R
#' objects).
#' @param drop for compatibility with the generic (unused).
#' @rdname phyDat
#' @export
"[.phyDat" <- function(x, i, j, ..., drop=FALSE){
   subset(x, subset = i, select = j, site.pattern=FALSE)
}


#' @rdname phangorn-internal
#' @importFrom fastmatch fmatch
#' @export
map_duplicates <-  function(x, dist=length(x)<500, ...){
  labels <- names(x)
  nr <- attr(x, "nr")

  if(dist){
    y <- as.matrix(dist.hamming(x, FALSE))
    l <- nrow(y)
#    z <- character(l)
#    for(i in seq_len(l)) z[i] <- paste( round(y[i, ], 8), collapse="_")
    ind <- duplicated(y)
  }
  else ind <- duplicated(x)
  res <- NULL
  if(any(ind)){
    ind2 <- grp_duplicated( matrix(unlist(x, recursive = FALSE, use.names = FALSE), nr, length(labels)), MARGIN=2)
    if(dist) ind2 <- grp_duplicated(y)
    ind2 <- ind2[ind]
    res <- data.frame(duplicates=labels[ind], where=labels[!ind][ind2], stringsAsFactors = FALSE)
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
  addTaxa <- FALSE
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
      addTaxa <- TRUE
      data <- subset(data, setdiff(names(data), dup[, 1]))
    }
    else break() # tmp <- FALSE
  }
  attr(data, "duplicated") <- dup_list
  data
}


#' @rdname phyDat
#' @export
allSitePattern <- function(n, levels=c("a", "c", "g", "t"), names=NULL){
  l <- length(levels)
  X <- matrix(NA_integer_, n, l^n)
  if(is.null(names))rownames(X) <- paste0("t", 1:n)
  else rownames(X) <- names
  for(i in 1:n)
    X[i, ] <- rep(rep(levels, each=l^(i-1)), l^(n-i))
  phyDat.default(X, levels, compress=FALSE, return.index=FALSE)
}


constSitePattern <- function(n, levels=c("a", "c", "g", "t"), names=NULL){
  l <- length(levels)
  X <- matrix(0, l, n)
  X <- matrix(rep(levels, each=n), n, l)
  if(is.null(names))rownames(X) <- paste0("t", 1:n)
  else rownames(X) <- names
  phyDat.default(X, levels)
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
#' \url{https://evolution.genetics.washington.edu/phylip/phylip.html}
#' @keywords IO
#' @export read.aa
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


#' @rdname phyDat
#' @export
genlight2phyDat <- function(x, ambiguity=NA){
  tmp <- as.matrix(x)
  lev <- na.omit(unique(as.vector(tmp)))
  phyDat(tmp, "USER", levels=lev, ambiguity=ambiguity)
}


#' @rdname phyDat
#' @method image phyDat
#' @export
image.phyDat <- function(x, ...){
  if(attr(x, "type") == "AA") image(as.AAbin(x), ...)
  if(attr(x, "type") == "DNA") image(as.DNAbin(x), ...)
  else return(NULL)
}
