#' @rdname phangorn-internal
#' @export
hash <- function (x, ...) UseMethod("hash")

#' @export
hash.phylo <- function(x, rooted=FALSE, ...){
  assert_phylo(x, is_rooted = rooted)
  x <- clean_phylo(x, unroot=!rooted, collapse.singles = TRUE, reorder = TRUE)
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


my_hash <- function(x, rooted=FALSE, ...){
  assert_phylo(x, is_rooted = rooted)
  x <- clean_phylo(x, unroot=!rooted, collapse.singles = TRUE, reorder = TRUE)
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

#' @rdname phangorn-internal
#' @export
hash_topo <- function(x, rooted=FALSE, ...){
  assert_treeish(x)
  if(inherits(x, "phylo")) return(my_hash(x, rooted=rooted, ...))
  if(inherits(x, "multiPhylo")){
    x <- .compressTipLabel(x)
    return(sapply(x, my_hash, rooted=rooted, ...))
  }
}
