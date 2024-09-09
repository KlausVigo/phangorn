#' @export
init_pml <- function(obj, k, ..., tree=NULL){
  attr(obj, "nSeq") <- length(obj) # add lengths
  if(is.null(attr(obj, "weight"))){
    attr(obj, "weight") <- rep(1, length(obj[[1]]))
  }
  if(!is.null(tree)) obj <- obj[tree$tip.label]
  weight <- attr(obj, "weight")
  f <- new(PML, obj, k)
  f
}
