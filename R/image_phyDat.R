#' Plot of a Sequence Alignment
#'
#' This function plots an image of an alignment of sequences.
#'
#' A wrapper for using \code{\link[ape]{image.DNAbin}} and \code{\link[ape]{image.AAbin}}.
#' Codons triplets are handled as nucleotide sequences.
#' @param x	 an object containing sequences, an object of class \code{phyDat}.
#' @param ... further arguments passed to or from other methods.
#' @returns Nothing. The function is called for plotting.
#' @seealso \code{\link[ape]{image.DNAbin}}, \code{\link[ape]{image.AAbin}}
#' @examples
#' data("chloroplast")
#' image(chloroplast[, 1:50], scheme="Clustal", show.aa = TRUE)
#' @method image phyDat
#' @export
image.phyDat <- function(x, ...){
  x <- gap_as_ambiguous(x)
  if(attr(x, "type") == "CODON") x <- codon2dna(x)
  if(attr(x, "type") == "AA") image(as.AAbin(x), ...)
  if(attr(x, "type") == "DNA") image(as.DNAbin(x), ...)
  if(attr(x, "type") == "USER") return(NULL)
}
