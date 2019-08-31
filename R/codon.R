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
#' Alignment gaps are simply and base ambiguities are currently ignored and
#' sites containing these are deleted.
#'
#' @aliases
#' dna2codon codon2dna
#' @param x An object containing sequences.
#' @param code The ncbi genetic code number for translation (see details).
#' By default the standard genetic code is used.
#' @param codonstart an integer giving where to start the translation. This
#' should be 1, 2, or 3, but larger values are accepted and have for effect to
#' start the translation further within the sequence.
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
dna2codon <- function(x, code=1, codonstart=1, ambiguity="---", ...){
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


phyDat.codon <- function (data, return.index = TRUE, ambiguity = "---",
                          NA_as_ambiguous=TRUE, code=1){
  if(is.matrix(data)) nam <- row.names(data)
  else nam <- names(data)
  if (inherits(data,"DNAbin"))
    data <- as.character(data)
  if(inherits(data, "character")) data <- as.matrix(data)
  if (is.matrix(data))
    data <- as.data.frame(t(data), stringsAsFactors = FALSE)
  else data <- as.data.frame(data, stringsAsFactors = FALSE)

  data <- data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)

  data[data=="u"] <- "t"

  splseq <- function (seq, frame = 0){
    starts <- seq(from = frame + 1, to = length(seq), by = 3L)
    sapply(starts, function(x) paste(seq[x:(x + 2L)], collapse=""))
  }

  data <- data.frame(lapply(data, splseq))
  compress <- TRUE
  if(nrow(data)==1) compress <- FALSE
  if(compress){
    ddd <- fast.table(data)
    data <- ddd$data
    weight <- ddd$weight
    index <- ddd$index
  }
  else{
    p <- length(data[[1]])
    weight <- rep(1, p)
    index <- 1:p
  }


  tmp <- .CODON[, as.character(code)]
  codon <- rownames(.CODON)
  codon <- codon[tmp != "*"] # no stop codons

  CODON <- diag(length(codon))

  if(NA_as_ambiguous){
    ambiguity <- unique(c("---", ambiguity))
  }
  if(ambiguity!=""){
    codon_amb <- c(codon, ambiguity)
    CODON <- rbind(CODON, matrix(1, length(ambiguity), length(codon)))
  }
  else codon_amb <- codon
  dimnames(CODON) <- list(codon_amb, codon)

  q <- length(data)
  p <- length(data[[1]])
  tmp <- vector("list", q)

  d <- dim(data)
  att <- attributes(data)
  data <- match(unlist(data), codon_amb)
  if(NA_as_ambiguous){
    ind <- match("---", codon_amb)
    data[is.na(data)] <- ind
  }
  attr(data, "dim") <- d
  data <- as.data.frame(data, stringsAsFactors=FALSE)
  attributes(data) <- att

  row.names(data) <- as.character(1:p)

  data <- na.omit(data)
  rn <- as.numeric(rownames(data))

  if(!is.null(attr(data, "na.action"))) warning("Found unknown characters. Deleted sites with with unknown states.")

  aaa <- match(index, attr(data, "na.action"))
  index <- index[is.na(aaa)]
  index <- match(index, unique(index))
  rn <- as.numeric(rownames(data))
  attr(data, "na.action") <- NULL

  weight <- weight[rn]
  p <- dim(data)[1]
  names(data) <- nam
  attr(data, "row.names") <- NULL
  attr(data, "weight") <- weight
  attr(data, "nr") <- p
  attr(data, "nc") <- length(codon)
  if (return.index)
    attr(data, "index") <- index
  attr(data, "levels") <- codon
  attr(data, "allLevels") <- codon_amb
  attr(data, "type") <- "CODON"
  attr(data, "code") <- code
  attr(data, "contrast") <- CODON
  class(data) <- "phyDat"
  data
}


