#' Treat gaps as a state
#'
#' The function \code{gap_as_state} changes the contrast of an phyDat object to
#' treat as its own state. Internally \code{phyDat} are stored similar to a
#' \code{factor} objects and only the contrast matrix and some attributes
#' change.
#'
#' @param obj An object of class phyDat.
#' @param gap a character which codes for the gaps (default is "-").
#' @param ambiguous a character which codes for the ambiguous state
#' @return The functions return an object of class \code{phyDat}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{phyDat}}, \code{\link{latag2n.phyDat}},
#' \code{\link[ape]{latag2n}}, \code{\link{ancestral.pml}},
#' \code{\link{gap_as_state}}
#' @keywords cluster
#' @examples
#' data(Laurasiatherian)
#' tmp <- gap_as_state(Laurasiatherian)
#' contr <- attr(tmp, "contrast")
#' rownames(contr) <- attr(tmp, "allLevels")
#' contr
#' @rdname gap_as_state
#' @export
gap_as_state <- function(obj, gap="-", ambiguous="?"){
  if(has_gap_state(obj)) return(obj)
#  if(!is.null(attr(obj, "gap_is_state")) & isTRUE(attr(obj, "gap_is_state")))
#    return(obj)
  contrast <- cbind(attr(obj, "contrast"), gap = 0)
  levels <- c(attr(obj, "levels"), gap)
  colnames(contrast) <- levels
  rownames(contrast) <- attr(obj, "allLevels")
  contrast[gap, ] <- 0
  contrast[gap, gap] <- 1
  # todo check for ambiguous
  contrast[ambiguous, "-"] <- 1
  rownames(contrast) <- NULL
  attr(obj, "levels") <- levels
  attr(obj, "nc") <- attr(obj, "nc") + 1L
  attr(obj, "contrast") <- contrast
#  attr(obj, "gap_is_state") <- TRUE
  obj
}


#' @rdname gap_as_state
#' @export
gap_as_ambiguous <- function(obj, gap="-"){
  if(!has_gap_state(obj)) return(obj)
#  if(is.null(attr(obj, "gap_is_state")) | !isTRUE(attr(obj, "gap_is_state")))
#    return(obj)
  contrast <- attr(obj, "contrast")
  levels <- attr(obj, "levels")
#  colnames(contrast) <- levels
  rownames(contrast) <- attr(obj, "allLevels")
  contrast[gap, ] <- 1
  rownames(contrast) <- NULL
  ind <- match(gap, levels)
  contrast <- contrast[,-ind]
  attr(obj, "levels") <- levels[-ind]
  attr(obj, "contrast") <- contrast
  attr(obj, "nc") <- attr(obj, "nc") - 1L
#  attr(obj, "gap_is_state") <- FALSE
  obj
}


#' @rdname gap_as_state
#' @export
has_gap_state <- function(obj){
  type <- attr(obj, "type")
  if(type=="DNA" && attr(obj, "nc")==5) return(TRUE)
  if(type=="AA" && attr(obj, "nc")==21) return(TRUE)
  FALSE
}


add_gap_Q_AA <- function(Q, rate_gap=0.1){
  res <- matrix(0, 20, 20)
  res[lower.tri(res)] <- Q
  res <- cbind(rbind(res,rate_gap), 0)
  res <- res[lower.tri(res)]
  res
}


add_gap_bf_AA <- function(bf, gap=.01){
  bf <- c(bf, gap)
  bf <- bf / sum(bf)
  bf
}


remove_similar <- function(x, k=3, index=FALSE){
  dm <- dist.hamming(x, FALSE)
  dm <- as.matrix(dm)
  ind_dist <- which(dm < (k+1), arr.ind = TRUE)
  ind_dist <- ind_dist[ind_dist[,1] < ind_dist[,2], ]
  dist_i <- unique(ind_dist[, 2])
  if(index) return(dist_i)
  x[-dist_i, ]
}
