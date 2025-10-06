#' Plot of a Sequence Alignment
#'
#' This function plots an image of an alignment of sequences.
#'
#' A wrapper for using \code{\link[ape]{image.DNAbin}} and
#' \code{\link[ape]{image.AAbin}}.
#' Codons triplets are transformed and handled as nucleotide sequences. So far
#' it is not yet possible to plot data with \code{ type="USER"}.
#' @param x	 an object containing sequences, an object of class \code{phyDat}.
#' @param ... further arguments passed to or from other methods.
#' @returns \code{image.phyDat} returns silently x.
#' @seealso \code{\link[ape]{image.DNAbin}}, \code{\link[ape]{image.AAbin}}
#' @examples
#' data("chloroplast")
#' image(chloroplast[, 1:50], scheme="Clustal", show.aa = TRUE)
#' @rdname image.phyDat
#' @method image phyDat
#' @export
image.phyDat <- function(x, ...){
  x <- gap_as_ambiguous(x)
  if(attr(x, "type") == "CODON") x <- codon2dna(x)
  if(attr(x, "type") == "AA") image(as.AAbin(x), ...)
  if(attr(x, "type") == "DNA") image(as.DNAbin(x), ...)
  if(attr(x, "type") == "USER") message('No plot function available for type="USER"')
  invisible(x)
}


#' @rdname image.phyDat
#' @method image ancestral
#' @export
image.ancestral <- function(x, ...) image(as.phyDat(x), ...)


#' @srrstats {EA4.0} is implicitly done e.g. distinguish between DNA, AA or user type data
NULL

