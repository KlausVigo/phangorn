#' @rdname phangorn-internal
#' @export
assert_phylo <- function(x, has_edge_length=FALSE, is_rooted=FALSE,
                         is_ultrametric=FALSE){
  txt <-  deparse(substitute(x))
  if (!inherits(x, "phylo")) stop(gettextf("%s must be of class phylo", txt))
  if (any(duplicated(x$tip.label)))
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
assert_pml <- function(x, contains_label=!is.null(label), label=NULL){
  txt <-  deparse(substitute(x))
  if (!inherits(x, "pml")) stop(gettextf("%s must be of class 'pml'", txt))
  invisible(x)
}

