#' Plot of a Sequence Alignment
#'
#' This function plots an image of an alignment of sequences.
#'
#' A wrapper for using \code{\link{image.DNAbin}} and \code{\link{image.AAbin}}.
#' @param x	 an object containing sequences, an object of class \code{phyDat}.
#' @param ... further arguments passed to or from other methods.
#' @seealso \code{\link{image.DNAbin}}, \code{\link{image.AAbin}}
#' @method image phyDat
#' @export
image.phyDat <- function(x, ...){
  if(attr(x, "type") == "AA") image(as.AAbin(x), ...)
  if(attr(x, "type") == "DNA") image(as.DNAbin(x), ...)
  else return(NULL)
}
