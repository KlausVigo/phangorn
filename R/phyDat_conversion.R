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
#' \code{acgt2ry} converts a \code{phyDat} object of nucleotides into an binary
#' ry-coded dataset.
#'
#' @aliases
#' as.phyDat.character as.phyDat.data.frame as.phyDat.matrix
#' as.MultipleAlignment as.MultipleAlignment.phyDat
#' acgt2ry phyDat2MultipleAlignment
#' @param data An object containing sequences.
#' @param x An object containing sequences.
#' @param obj as object of class phyDat
#' @param type Type of sequences ("DNA", "AA", "CODON" or "USER").
#' @param levels Level attributes.
#' @param return.index If TRUE returns a index of the site patterns.
#' @param allLevels return original data.
#' @param ambiguity character for ambiguous character and no contrast is
#' provided.
#' @param ... further arguments passed to or from other methods.
#' @return The functions return an object of class \code{phyDat}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso [DNAbin()], [as.DNAbin()],
#' \code{\link{baseFreq}}, \code{\link{glance.phyDat}},
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
#' # transform as characters
#' LauraChar <- as.character(Laurasiatherian)
#' # and back
#' Laura <- phyDat(LauraChar)
#' all.equal(Laurasiatherian, Laura)
#' LauraDNAbin <- as.DNAbin(Laurasiatherian)
#' all.equal(Laurasiatherian, as.phyDat(LauraDNAbin))
#'
#' @rdname as.phyDat
#' @export
phyDat <- function (data, type="DNA", levels=NULL, return.index = TRUE, ...){
  #  if (inherits(data, "DNAbin")) type <- "DNA"
  pt <- match.arg(type, c("DNA", "AA", "CODON", "USER"))
  if(pt=="DNA") dat <- phyDat.DNA(data, return.index=return.index, ...)
  if(pt=="AA") dat <- phyDat.AA(data, return.index=return.index, ...)
  if(pt=="CODON") dat <- phyDat.codon(data, return.index=return.index, ...)
  if(pt=="USER") dat <- phyDat.default(data, levels = levels,
                                       return.index=return.index, ...)
  dat
}


#' @rdname as.phyDat
#' @export
as.phyDat <- function (x, ...){
  if (inherits(x, "phyDat")) return(x)
  UseMethod("as.phyDat")
}


#' @rdname as.phyDat
#' @method as.phyDat factor
#' @export
as.phyDat.factor <- function(x, ...){
  nam <- names(x)
  lev <- levels(x)
  x <- as.character(x)
  names(x) <- nam
  phyDat(x, type="USER", levels = lev, ...)
}


#' @rdname as.phyDat
#' @method as.phyDat DNAbin
#' @export
as.phyDat.DNAbin <- function(x, ...) phyDat.DNA(x, ...)


#' @rdname as.phyDat
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
#' @rdname as.phyDat
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


#' @rdname as.phyDat
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


#' @rdname as.phyDat
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


#' @rdname as.phyDat
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


#' @rdname as.phyDat
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



#' @rdname as.phyDat
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


#' @rdname as.phyDat
#' @export
as.AAbin.phyDat <- function(x,...) {
  if(attr(x, "type")=="AA") return(as.AAbin(as.character(x, ...)))
  else stop("x must be a amino acid sequence")
}


#' @rdname as.phyDat
#' @export
genlight2phyDat <- function(x, ambiguity=NA){
  tmp <- as.matrix(x)
  lev <- na.omit(unique(as.vector(tmp)))
  phyDat(tmp, "USER", levels=lev, ambiguity=ambiguity)
}



#' @rdname as.phyDat
#' @export
acgt2ry <- function(obj){
  ac <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y",
          "k", "v", "h", "d", "b", "n", "?", "-")
  AC <- matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1),
                 c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1),
                 c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1)),
               18, 4, dimnames = list(NULL, c("a", "c", "g", "t")))
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
