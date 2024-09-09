#' Summaries of alignments
#'
#' \code{baseFreq} computes the frequencies (absolute or relative) of the states
#' from a sample of sequences.
#' \code{glance} computes some useful information about the alignment.
#' \code{composition\_test} computes a \eqn{\chi^2}-test testing if the state
#' composition for a species differs.
#'
#' @param obj,x as object of class phyDat
#' @param freq logical, if 'TRUE', frequencies or counts are returned otherwise
#' proportions
#' @param all all a logical; if all = TRUE, all counts of bases, ambiguous
#' codes, missing data, and alignment gaps are returned as defined in the
#' contrast.
#' @param drop.unused.levels logical, drop unused levels
#' @param ... further arguments passed to or from other methods.
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
  parsimony_informative_sites <- sum(attr(x, "weight")[-pis[, 1]])
  data.frame(nseq=nseq, nchar=nchar, unique_site_pattern=unique_sites,
             parsimony_informative_sites=parsimony_informative_sites,
             const_sites=const_site(x))
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
                                    c("statistic", "parameter df", "p-value")))
  for(i in seq_len(n)){
    tmp <- baseFreq(obj[i], freq=TRUE)
    res[i, ] <- unlist(chisq.test(rbind(ALL-tmp, tmp))[1:3])
  }
  res
}

