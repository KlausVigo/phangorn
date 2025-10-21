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
assert_phylo <- function(x, has_edge_length=FALSE, is_rooted=FALSE,
                         is_ultrametric=FALSE, is_binary=FALSE,
                         unique_tiplabel=TRUE){
  txt <-  deparse(substitute(x))
  if (!inherits(x, "phylo")) stop(gettextf("%s must be of class phylo", txt))
  if (has_edge_length && is.null(x$edge.length))
    stop(gettextf("%s must have edge weights", txt))
  if(is_rooted && !is.rooted(x)) stop(gettextf("%s must be rooted!", txt))
  if(is_ultrametric && !is.ultrametric(x))
    stop(gettextf("%s must be ultrametric!", txt))
  if(is_binary && !is.binary(x)) stop(gettextf("%s must be binary!", txt))
  if (unique_tiplabel && any(duplicated(x$tip.label)))
    stop(gettextf("%s must have unique labels!", txt))
  invisible(x)
}


#' @rdname phangorn-internal
#' @export
assert_multiPhylo <- function(x, has_edge_length=FALSE, is_rooted=FALSE,
                              is_ultrametric=FALSE, is_binary =FALSE){
  # same_labels=FALSE,
  txt <-  deparse(substitute(x))
  if (!inherits(x, "multiPhylo"))
    stop(gettextf("%s must be of class 'multiPhylo'", txt))
  if(has_edge_length && any(sapply(x, \(x)is.null(x$edge.length))) )
    stop(gettextf("All trees in %s must have edge weights", txt))
  if(is_rooted && !all(is.rooted(x)))
    stop(gettextf("All trees in %s must be rooted!", txt))
  if(is_ultrametric && !all(is.ultrametric(x)))
    stop(gettextf("All trees in %s must be ultrametric!", txt))
  if(is_binary && !all(is.binary(x)))
     stop(gettextf("All trees in %s must be binary!", txt))
  invisible(x)
}




has_edge_length <- function(x){
  assert_treeish(x)
  if(inherits(x, "phylo")) return(!is.null(x$edge.length))
  if(inherits(x, "multiPhylo"))
    return(!vapply(x, \(x)is.null (x$edge.length), TRUE))
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
clean_phylo <- function(x, ..., unroot=FALSE, collapse.singles=FALSE,
                        reorder=FALSE, order="postorder",
                        multi2di=FALSE, di2multi=FALSE, tol=1e-7,
                        compress=FALSE){
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
  if(multi2di && di2multi) stop("di2multi and multi2di cannot both be TRUE!")
  if(multi2di && !all(is.binary(x))) x <- multi2di(x)
  if(di2multi) x <- di2multi(x, tol=tol)
  if(unroot && any(is.rooted(x))) x <- unroot(x)
  if(reorder) x <- reorder(x, order=order)
  x
}

## TODO is.binary di2multi


