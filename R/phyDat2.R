# replaces fast.match & fast.match2
grp_duplicated <- function(x, MARGIN = 1, factor=FALSE, fromLast = FALSE, ...)
{
  ans <- .Call(grpDupAtomMat, x, as.integer(MARGIN), as.logical(fromLast))
  if(fromLast) ans[] <- (attr(ans, 'nlevels'):1L)[ans] # ensure the group ids agree with row/col index of result from "unique"
  if(factor) {
    attr(ans, 'levels') <- as.character(seq_len(attr(ans, 'nlevels')))
    class(ans) <- 'factor'
  }
  ans
}


phyDat.DNA <- function (data, return.index = TRUE, compress = TRUE){
  if (is.matrix(data)) nam <- row.names(data)
  else nam <- names(data)
  if(inherits(data, "list")) data <- as.data.frame(data)
  if(inherits(data, "data.frame")) data <- t(as.matrix(data))
  if(inherits(data, "character")){
    data <- as.matrix(data)
  }
  if (inherits(data,"DNAbin")){
    if(is.list(data)) data <- as.matrix(data)
  }
  if(mode(data)=="character") data <- tolower(data)
  ac <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y",
          "k", "v", "h", "d", "b", "n", "?", "-")
  AC <- matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1),
                 c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1),
                 c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1)),
               18, 4, dimnames = list(NULL, c("a", "c", "g", "t")))
  if(ncol(data)==1) compress <- FALSE
  if(compress){
    index <- grp_duplicated(data, MARGIN=2)
    attr(index, "nlevels") <- NULL
    weight <- tabulate(index, max(index))
    ind <- which(!duplicated(index))
  }
  else{
    p <- ncol(data)
    weight <- rep(1, p)
    index <- ind <- 1:p
  }
  n <- nrow(data)
  res <- vector("list", nrow(data))
  ind_na <- logical(length(weight))
  if(inherits(data,"DNAbin")){
    cs <- rep(NA_integer_, 256)
    # from ape ._bs_ and ._cs_
    .cs <- c("a", "g", "c", "t", "r", "m", "w", "s", "k", "y", "v", "h",
            "d", "b", "n", "-", "?")
    .bs <- c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48, 224, 176, 208,
                   112, 240, 4, 2)
    cs[.bs] <- match(.cs, ac)
    for(i in seq_len(n)){
      res[[i]] <- cs[ as.integer(data[i, ind]) ]
      #      maybe check for NAs
      #      ind_na[is.na(res[[i]])] <- TRUE
    }
  } else {
    for(i in seq_len(n)){
      res[[i]] <- match(data[i, ind], ac)
      ind_na[is.na(res[[i]])] <- TRUE
    }
  }
  # this needs testing
  if(any(ind_na)){
    res <- lapply(res, function(x, ind_na)x[!ind_na], ind_na)
    weight <- weight[!ind_na]
    index <- index[which(ind_na == FALSE)]
    index <- match(index, unique(index))
  }
  names(res) <- nam
  attr(res, "row.names") <- NULL
  attr(res, "weight") <- weight
  attr(res, "nr") <- length(weight)
  attr(res, "nc") <- 4
  if (return.index) attr(res, "index") <- index
  attr(res, "levels") <- c("a", "c", "g", "t")
  attr(res, "allLevels") <- ac
  attr(res, "type") <- "DNA"
  attr(res, "contrast") <- AC
  class(res) <- "phyDat"
  res
}


phyDat.default <- function (data, levels = NULL, return.index = TRUE,
                            contrast = NULL, ambiguity = "?", compress=TRUE, ...){
  if (is.matrix(data))
    nam <- row.names(data)
  else nam <- names(data)
  if(is.null(nam))stop("data object must contain taxa names")
  if(inherits(data, "list")) data <- as.data.frame(data)
  if(inherits(data, "data.frame")) data <- t(as.matrix(data))
  if(inherits(data, "character") | inherits(data, "numeric"))
    data <- as.matrix(data)
  if (inherits(data, "DNAbin") | inherits(data, "AAbin") |
      inherits(data, "phyDat")) data <- as.character(data)
  if(ncol(data)==1) compress <- FALSE
  if(compress){
    index <- grp_duplicated(data, MARGIN=2)
    attr(index, "nlevels") <- NULL
    weight <- tabulate(index, max(index))
    ind <- which(!duplicated(index))
  }
  else{
    p <- ncol(data)
    weight <- rep(1, p)
    index <- ind <- 1:p
  }

  if (!is.null(contrast)) {
    levels <- colnames(contrast)
    all.levels <- rownames(contrast)
    rownames(contrast) <- NULL
  }
  else {
    if (is.null(levels))
      stop("Either argument levels or contrast has to be supplied")
    l <- length(levels)
    contrast <- diag(l)
    all.levels <- levels
    if (!is.null(ambiguity)) {
      all.levels <- c(all.levels, ambiguity)
      k <- length(ambiguity)
      if (k > 0)
        contrast <- rbind(contrast, matrix(1, k, l))
    }
  }
  n <- nrow(data)
  res <- vector("list", nrow(data))
  ind_na <- logical(length(weight))
  for(i in seq_len(n)){
    res[[i]] <- match(data[i, ind], all.levels)
    ind_na[is.na(res[[i]])] <- TRUE
  }
  if(any(ind_na)){
    warning("Found unknown characters (not supplied in levels). Deleted sites with with unknown states.")
    res <- lapply(res, function(x, ind_na)x[!ind_na], ind_na)
    weight <- weight[!ind_na]
    index <- index[which(ind_na == FALSE)]
    index <- match(index, unique(index))
  }
  p <- length(weight)
  names(res) <- nam
  attr(res, "weight") <- weight
  attr(res, "nr") <- p
  attr(res, "nc") <- length(levels)
  if (return.index) attr(res, "index") <- index
  attr(res, "levels") <- levels
  attr(res, "allLevels") <- all.levels
  attr(res, "type") <- "USER"
  attr(res, "contrast") <- contrast
  class(res) <- "phyDat"
  res
}


phyDat.AA <- function (data, return.index = TRUE){
  if(is.matrix(data)) nam <- row.names(data)
  else nam <- names(data)
  if(inherits(data, "list")) data <- as.data.frame(data)
  if(inherits(data, "data.frame")) data <- t(as.matrix(data))
  # AAbin
  if (inherits(data,"AAbin")){
    if(is.list(data)) data <- as.matrix(data)
    data <- as.character(data)
  }
  if(inherits(data, "character")) data <- as.matrix(data)
  data <- toupper(data)
#  aa <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i",
#          "l", "k", "m", "f", "p", "s", "t", "w", "y", "v")
#  aa2 <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i",
#           "l", "k", "m", "f", "p", "s", "t", "w", "y", "v", "b",
#           "z", "x", "-", "?")
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
          "P", "S", "T", "W", "Y", "V")
  aa2 <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
           "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "-", "?")
  AA <- diag(20)
  AA <- rbind(AA, matrix(0, 5, 20))
  AA[21, 3] <- AA[21, 4] <- 1 # Aspartate or Asparagine
  AA[22, 6] <- AA[22, 7] <- 1 #
  AA[23:25, ] <- 1
  dimnames(AA) <- list(aa2, aa)
  compress <- TRUE
  if(ncol(data)==1) compress <- FALSE
  if(compress){
    index <- grp_duplicated(data, MARGIN=2)
    attr(index, "nlevels") <- NULL
    weight <- tabulate(index, max(index))
    ind <- which(!duplicated(index))
  }
  else{
    p <- ncol(data)
    weight <- rep(1, p)
    index <- ind <- 1:p
  }

  n <- nrow(data)
  res <- vector("list", n)
  ind_na <- logical(length(weight))
  for(i in seq_len(n)){
    res[[i]] <- match(data[i, ind], aa2)
    ind_na[is.na(res[[i]])] <- TRUE
  }
  if(any(ind_na)){
    warning("Found unknown characters (not supplied in levels). Deleted sites with with unknown states.")
    res <- lapply(res, function(x, ind_na)x[!ind_na], ind_na)
    weight <- weight[!ind_na]
    index <- index[which(ind_na == FALSE)]
    index <- match(index, unique(index))
  }
  names(res) <- nam
  attr(res, "weight") <- weight
  attr(res, "nr") <- length(weight)
  attr(res, "nc") <- 20
  if (return.index)
    attr(res, "index") <- index
  attr(res, "levels") <- aa
  attr(res, "allLevels") <- aa2
  attr(res, "type") <- "AA"
  attr(res, "contrast") <- AA
  class(res) <- "phyDat"
  res
}


phyDat.codon <- function (data, return.index = TRUE, ambiguity = "---",
                          NA_as_ambiguous=TRUE, code=1, stopcodon="exclude"){
  if(is.matrix(data)) nam <- row.names(data)
  else nam <- names(data)
  if(inherits(data, "list")) data <- as.data.frame(data)
  if(inherits(data, "data.frame")) data <- t(as.matrix(data))
  if(inherits(data, "character")){
    data <- as.matrix(data)
    data <- tolower(data)
  }
  if (inherits(data,"DNAbin") || inherits(data, "phyDat")) data <- as.character(data)

  data[data=="u"] <- "t"
  stopcodon <- match.arg(stopcodon, c("exclude", "include"))

  splseq <- function (seq, frame = 0){
    starts <- seq(from = frame + 1, to = length(seq), by = 3L)
    sapply(starts, function(x) paste(seq[x:(x + 2L)], collapse=""))
  }

  data <- apply(data, 1, splseq)
  tmp <- .CODON[, as.character(code)]
  codon <- rownames(.CODON)
  stop_codons <- codon[tmp == "*"]
  if(stopcodon!="include") codon <- codon[tmp != "*"] # no stop codons

  if(stopcodon=="exclude"){
    for(i in stopcodon){
      data[data==i] <- NA
    }
    rm_stop <- which(is.na(data), arr.ind = TRUE)[,1]
    if(length(rm_stop)>0) data[-unique(rm_stop), ]
  }

  compress <- TRUE
  if(nrow(data)==1) compress <- FALSE
  if(compress){
    index <- grp_duplicated(data, MARGIN=1)
    attr(index, "nlevels") <- NULL
    weight <- tabulate(index, max(index))
    ind <- which(!duplicated(index))
  }
  else{
    p <- nrow(data)
    weight <- rep(1, p)
    index <- ind <- 1:p
  }

  CODON <- diag(length(codon))
  if(NA_as_ambiguous){
    ambiguity <- unique(c("---", ambiguity))
  }
  if(ambiguity!=""){
    codon_amb <- c(codon, ambiguity)
    CODON <- rbind(CODON, matrix(1, length(ambiguity), length(codon)))
  }
  else codon_amb <- codon
  dimnames(CODON) <- list(codon_amb, codon)

  n <- ncol(data)
  res <- vector("list", n)
  ind_na <- logical(length(weight))
  for(i in seq_len(n)){
    res[[i]] <- match(data[ind, i], codon_amb)
    ind_na[is.na(res[[i]])] <- TRUE
  }
  if(NA_as_ambiguous){
    ind <- match("---", codon_amb)
    res <- lapply(res, function(x, ind){x[is.na(x)] <- ind; x}, ind)
  }
  else if(any(ind_na)){
    warning("Found unknown characters. Deleted sites with with unknown states.")
    res <- lapply(res, function(x, ind_na)x[!ind_na], ind_na)
    weight <- weight[!ind_na]
    index <- index[which(ind_na == FALSE)]
    index <- match(index, unique(index))
  }

  p <- length(weight)
  names(res) <- nam
  attr(res, "row.names") <- NULL
  attr(res, "weight") <- weight
  attr(res, "nr") <- p
  attr(res, "nc") <- length(codon)
  if (return.index) attr(res, "index") <- index
  attr(res, "levels") <- codon
  attr(res, "allLevels") <- codon_amb
  attr(res, "type") <- "CODON"
  attr(res, "code") <- code
  attr(res, "contrast") <- CODON
  class(res) <- "phyDat"
  res
}

