#' Translate nucleic acid sequences into codons
#'
#' The function transforms \code{dna2codon} DNA sequences to codon sequences,
#' \code{codon2dna} transform the other way.
#'
#' The following genetic codes are described here. The number preceding each
#' corresponds to the code argument.
#' \tabular{rl}{
#'   1 \tab standard   \cr
#'   2 \tab vertebrate.mitochondrial \cr
#'   3 \tab yeast.mitochondrial \cr
#'   4 \tab protozoan.mitochondrial+mycoplasma \cr
#'   5 \tab invertebrate.mitochondrial \cr
#'   6 \tab ciliate+dasycladaceal \cr
#'   9 \tab echinoderm+flatworm.mitochondrial \cr
#'   10 \tab euplotid \cr
#'   11 \tab bacterial+plantplastid \cr
#'   12 \tab alternativeyeast \cr
#'   13 \tab ascidian.mitochondrial \cr
#'   14 \tab alternativeflatworm.mitochondrial \cr
#'   15 \tab blepharism \cr
#'   16 \tab chlorophycean.mitochondrial \cr
#'   21 \tab trematode.mitochondrial \cr
#'   22 \tab scenedesmus.mitochondrial \cr
#'   23 \tab thraustochytrium.mitochondria \cr
#'   24 \tab Pterobranchia.mitochondrial \cr
#'   25 \tab CandidateDivision.SR1+Gracilibacteria \cr
#'   26 \tab Pachysolen.tannophilus
#' }
#' Alignment gaps and ambiguities are currently ignored and sites containing
#' these are deleted.
#'
#' @aliases
#' dna2codon codon2dna
#' @param x An object containing sequences.
#' @param codonstart an integer giving where to start the translation. This
#' should be 1, 2, or 3, but larger values are accepted and have for effect to
#' start the translation further within the sequence.
#' @param code The ncbi genetic code number for translation (see details).
#' By default the standard genetic code is used.
#' @param ambiguity character for ambiguous character and no contrast is
#' provided.
#' @param ... further arguments passed to or from other methods.
#' @return The functions return an object of class \code{phyDat}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes}
#' @seealso \code{\link{trans}}, \code{\link{phyDat}} and the chapter 4 in the
#' \code{vignette("phangorn-specials", package="phangorn")}
#' @keywords cluster
#' @examples
#'
#' data(Laurasiatherian)
#' class(Laurasiatherian)
#' Laurasiatherian
#' dna2codon(Laurasiatherian)
#'
#' @rdname dna2codon
#' @export
dna2codon <- function(x, codonstart=1, code=1, ambiguity="---", ...){
  if(!inherits(x, "phyDat"))stop("x needs to be of class phyDat!")
  if(attr(x, "type")=="AA")stop("x needs to be a nucleotide sequence!")

  if(codonstart>1){
    del <- -seq_len(codonstart)
    x <- subset(x, select=del, site.pattern=FALSE)
  }
  n_sites <- sum(attr(x, "weight"))
  if( (n_sites %% 3) ){
    keep <- seq_len( (n_sites %/% 3) * 3 )
    x <- subset(x, select=keep, site.pattern=FALSE)
  }
  phyDat.codon(as.character(x), ambiguity=ambiguity, code=code, ...)
}


#' @rdname dna2codon
#' @export
codon2dna <- function(x){
  if(!inherits(x, "phyDat"))stop("x needs to be of class phyDat!")
  phyDat.DNA(as.character(x))
}


synonymous_subs <- function(code=1, stop.codon=FALSE){
  tmp <- .CODON[, as.character(code)]
  label <- rownames(.CODON)
  l <- length(tmp)
  res <- matrix(1, 64, 64, dimnames = list(label, label))
  for(i in seq_len(64)){
    for(j in seq_len(64)) {
      if(tmp[i] == tmp[j]) res[i, j] <- 0
    }
  }
  res[.ONE_TRANSITION == FALSE] <- -1
  if(!stop.codon){
    label <- label[tmp != "*"]
    res <- res[label, label]
  }
  res[lower.tri(res)]
}


tstv_subs <- function(code=1, stop.codon=FALSE){
  tmp <- .CODON[, as.character(code)]
  label <- rownames(.CODON)
  res <- .SUB
  if(!stop.codon){
    label <- label[tmp != "*"]
    res <- res[label, label]
  }
  res[lower.tri(res)]
}

