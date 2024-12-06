#' Summaries of alignments
#'
#' \code{baseFreq} computes the frequencies (absolute or relative) of the states
#' from a sample of sequences.
#' \code{glance} computes some useful information about the alignment.
#' \code{composition\_test} computes a \eqn{\chi^2}-test testing if the state
#' composition for a species differs.
#'
#' @param object,obj,x as object of class phyDat
#' @param freq logical, if 'TRUE', frequencies or counts are returned otherwise
#' proportions
#' @param all all a logical; if all = TRUE, all counts of bases, ambiguous
#' codes, missing data, and alignment gaps are returned as defined in the
#' contrast.
#' @param drop.unused.levels logical, drop unused levels
#' @param ... further arguments passed to or from other methods.
#' @param digits minimal number of significant digits.
#'
#' @return  \code{baseFreq} returns a named vector and \code{glance} a one row
#' \code{data.frame}.
#' @seealso \code{\link{phyDat}, \link[ape]{base.freq}, \link{glance}}
#' @author Klaus Schliep
#' @examples
#'
#' data(Laurasiatherian)
#' data(chloroplast)
#' # base frequencies
#' baseFreq(Laurasiatherian)
#' baseFreq(Laurasiatherian, all=TRUE)
#' baseFreq(Laurasiatherian, freq=TRUE)
#' baseFreq(chloroplast)
#' glance(Laurasiatherian)
#' glance(chloroplast)
#' composition_test(Laurasiatherian)[1:10,]
#' @rdname baseFreq
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
  if(!all) res <- res[as.character( attr(obj, "levels") )]
  if(!freq)res <- res/sum(res)
  if(drop.unused.levels) return(res[res>0])
  res
}


const_site <- function(x){
  tmp <- lli(x)
  ind <- which(rowSums(tmp)>1e-6)
  sw <- sum(attr(x, "weight"))
  weight <- attr(x, "weight")[ind]
  # list(index=ind, weight=attr(x, "weight")[ind], M=tmp[ind,])
  sum(weight)
}

#' @importFrom generics glance
#' @export
generics::glance


#' @rdname baseFreq
#' @export
glance.phyDat <- function (x, ...){
  nc  <- attr(x, "nc")
  nseq <- length(x)
  nchar <- sum(attr(x,"weight"))
  unique_sites <- attr(x, "nr")
  pis <- parsinfo(x)
  type <- attr(x, "type")

  gap <- "-" # needs check
  contr <- attr(x, "contrast")
  levels <- attr(x, "levels")
  allLevels <- attr(x, "allLevels")
  ambiguous_states <- allLevels[rowSums(contr>0) > 1]
  ambiguous_states <- ambiguous_states[-match(gap, ambiguous_states)]
  BF <- baseFreq(x, freq=TRUE, all=TRUE)
  amb <- sum(BF[ambiguous_states])
  gaps <- sum(BF[gap])

  parsimony_informative_sites <- sum(attr(x, "weight")[-pis[, 1]])
  duplicats <- sum(duplicated(x))
  data.frame(nseq=nseq, nchar=nchar, unique_site_pattern=unique_sites,
             parsimony_informative_sites=parsimony_informative_sites,
             const_sites=const_site(x), duplicated_seq=duplicats,
             gaps=gaps, ambiguous=amb, type=type)
}


#' @rdname baseFreq
#' @importFrom stats chisq.test
#' @export
composition_test <- function(obj){
  stopifnot(inherits(obj,"phyDat"))
  labels <- attr(obj, "allLevels")
  levs <- attr(obj, "levels")
  weight <- attr(obj,"weight")
  n <- length(obj)
  ALL <- baseFreq(obj, freq=TRUE)
  res <- matrix(0, n, 3, dimnames = list(names(obj),
                                    c("statistic", "df", "p-value")))
  for(i in seq_len(n)){
    tmp <- baseFreq(obj[i], freq=TRUE)
    res[i, ] <- unlist(chisq.test(rbind(ALL-tmp, tmp))[1:3])
  }
  res
}


#' @rdname baseFreq
#' @export
summary.phyDat <- function(object, ...){
  bf <- baseFreq(object)
  gl <- glance(object, ...)
  ct <- composition_test(object)
  ans <- list(glance=gl, composition=ct)
  class(ans) <- "summary.phyDat"
  ans
}

#' @rdname baseFreq
#' @export
print.summary.phyDat <- function(x, ...,
                                 digits = max(3L, getOption("digits") - 3L)){
  cat("Alignment statistics \n\n")
  cat("Type: ", x$glance$type, "\n")
  cat("Sequences: ", x$glance$nseq, "\n")
  cat("Columns: ", x$glance$nchar, "\n")
  cat("Site pattern: ", x$glance$unique_site_pattern, "\n")
  cat("Parsimony informative sites: ", x$glance$parsimony_informative_sites,
      "\n")
  cat("Constant sites: ", x$glance$const_sites, "\n")
  cat("Number of gaps: ", x$glance$gaps, "\n")
  cat("Number of ambiguous states: ", x$glance$amb, "\n")
  cat("Duplicated sequences: ", x$glance$duplicated_seq, "\n\n")
  cat("Composition Test (Chisq) \n")
  comp <- x$composition
  comp <- comp[order(comp[,3]),]
  if(any(comp[,3] < 0.05)) ind <- max(which(comp[,3] < 0.05))
  else ind <- min(6, x$glance$nseq)
  print(comp[seq_len(ind),], digits=digits, ...)
}
