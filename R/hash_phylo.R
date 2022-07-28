#' @rdname phangorn-internal
#' @export
hash_phylo <- function(x, ...){
  if(has.singles(x)) x <- collapse.singles(x)
  x <- unroot(x)
  x <- reorder(x, "postorder")
  nTips <- as.integer(length(x$tip.label))
  res <- .Call('_phangorn_short_bipCPP', x$edge, nTips)
  digest::digest(res, ...)
}
