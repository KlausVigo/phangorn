#' Replace leading and trailing alignment gaps with an ambiguous state
#'
#' Substitutes leading and trailing alignment gaps in aligned sequences into N
#' (i.e., A, C, G, or T) or ?. The gaps in the middle of the sequences are left
#' unchanged.
#' @param x an object of class \code{phyDat}.
#' @param amb character of the ambiguous state t replace the gaps.
#' @param gap gap parameter to replace.
#' @param ... Further arguments passed to or from other methods.
#' @returns returns an object of class \code{phyDat}.
#' @seealso \code{\link[ape]{latag2n}}, \code{\link{ancestral.pml}},
#' \code{\link{gap_as_state}}
#' @examples
#' x <- phyDat(matrix(c("-", "A", "G", "-", "T", "C"), 2, 3))
#' y <- latag2n.phyDat(x)
#' image(x)
#' image(y)
#' @keywords cluster
#' @export
latag2n.phyDat <- function(x, amb=ifelse(attr(x,"type")=="DNA", "N", "?"),
                           gap="-", ...){
  X <- as.character(x)
  d <- dim(X)
  for(i in seq_len(d[1])){
    # from start
    j <- 1
    while(X[i,j]==gap){
      X[i,j] <- amb
      j <- j+1
      if(j>d[2])break()
    }
    j <- d[2]
    while(X[i,j]==gap){
      X[i,j] <- amb
      j <- j-1
      if(j<1)break()
    }
  }
  type <- attr(x, "type")
  if(type %in% c("DNA", "AA", "CODON")){
    y <- phyDat(X, type=type)
  } else y <- phyDat(X, type="USER", levels = attr(x, "levels"))
  y
}
# same as:
# dna2 <- as.DNAbin(dna) |> latag2n() |> as.phyDat()
# but works for AA
