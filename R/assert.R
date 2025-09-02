#' @rdname phangorn-internal
#' @export
assert_phylo <- function(x, has_edge_length=FALSE, is_rooted=FALSE,
                         is_ultrametric=FALSE, unique_tiplabel=TRUE){
  txt <-  deparse(substitute(x))
  if (!inherits(x, "phylo")) stop(gettextf("%s must be of class phylo", txt))
  if (unique_tiplabel && any(duplicated(x$tip.label)))
    stop(gettextf("%s must have unique labels!", txt))
  if (has_edge_length && is.null(x$edge.length))
    stop(gettextf("%s must have edge weights", txt))
  if(is_rooted && !is.rooted(x)) stop(gettextf("%s must be rooted!", txt))
  if(is_ultrametric && !is.ultrametric(x))
    stop(gettextf("%s must be ultrametric!", txt))
  invisible(x)
}


#' @rdname phangorn-internal
#' @export
assert_treeish <- function(x, null.ok=FALSE){
  txt <-  deparse(substitute(x))
  if (missing(x))
    stop(gettextf("argument \"%s\" is missing, with no default", txt))
  treeish <- inherits(x, "phylo") || inherits(x, "multiPhylo")
  if(null.ok) treeish <- treeish || is.null(x)
  if(!treeish) stop(gettextf("%s must be of class 'phylo' or 'multiPhylo'", txt))
  invisible(x)
}



#' @rdname phangorn-internal
#' @export
assert_multiPhylo <- function(x, is_ultrametric=FALSE, is_rooted=FALSE,
                              same_labels=FALSE,  has_edge_length=FALSE){
  txt <-  deparse(substitute(x))
  if (!inherits(x, "multiPhylo"))
    stop(gettextf("%s must be of class 'multiPhylo'", txt))
  invisible(x)
}


#' @rdname phangorn-internal
#' @export
assert_phyDat <- function(x, contains_label=!is.null(label), label=NULL){
  txt <-  deparse(substitute(x))
  if (!inherits(x, "phyDat")) stop(gettextf("%s must be of class 'phyDat'", txt))
  if (contains_label && anyNA(match(label, attr(x, "names"))))
    stop(gettextf("tip labels are not in %s", txt))
  invisible(x)
}


#' @rdname phangorn-internal
#' @export
assert_pml <- function(x){
  txt <-  deparse(substitute(x))
  if (!inherits(x, "pml")) stop(gettextf("%s must be of class 'pml'", txt))
  invisible(x)
}


#' @rdname phangorn-internal
#' @export
clean_phylo <- function(x, unroot=FALSE, multi2di=FALSE,
                        collapse.singles=FALSE, reorder=FALSE, compress=FALSE){
  assert_treeish(x)
  is_multiPhylo <- inherits(x, "multiPhylo")
  compressed <- !is.null(attr(x, "TipLabel"))
  if(collapse.singles && any(hs <- has.singles(x))){
    if(is_multiPhylo){
      x <- .uncompressTipLabel(x)
      for(i in which(hs)) x[[i]] <- collapse.singles(x[[i]])
    } else x <- collapse.singles(x)
  }
  if(compress && is_multiPhylo) x <- .compressTipLabel(x)
  if(multi2di && !all(is.binary(x))) x <- multi2di(x)
  if(unroot && any(is.rooted(x))) x <- unroot(x)
  if(reorder) x <- reorder(x, "postorder")
  x
}


