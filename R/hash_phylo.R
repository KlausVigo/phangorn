#' @rdname phangorn-internal
#' @export
hash <- function (x, ...) UseMethod("hash")

#' @export
hash.phylo <- function(x, rooted=FALSE, ...){
  if(has.singles(x)) x <- collapse.singles(x)
  if(!rooted)x <- unroot(x)
  if(rooted && !is.rooted(x)) stop("x must be rooted")
  x <- reorder(x, "postorder")
  nTips <- as.integer(length(x$tip.label))
  if(rooted){
    res <- .Call('_phangorn_sorted_bipartCPP', x$edge, nTips)
    l <- lengths(res)
    res <- res[l<nTips]
  }
  else {
    res <- .Call('_phangorn_short_bipCPP', x$edge, nTips)
    l <- lengths(res)
    res <- res[l>1]
  }
  digest::digest(res, ...)
}

#' @export
hash.multiPhylo <- function(x, ...){
  x <- .compressTipLabel(x)
  sapply(x, hash, ...)
}
