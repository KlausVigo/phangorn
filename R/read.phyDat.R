
read.fasta.user <- function (file, skip = 0, nlines = 0,
                             comment.char = "#", seq.names = NULL){
  getTaxaNames <- function(x) {
    x <- sub("^ +", "", x)
    x <- sub(" +$", "", x)
    x <- sub("^['\"]", "", x)
    x <- sub("['\"]$", "", x)
    x
  }

  X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE,
            skip = skip, nlines = nlines, comment.char = comment.char)

  start <- grep("^ {0,}>", X)
  taxa <- X[start]
  n <- length(taxa)
  obj <- vector("list", n)
  if (is.null(seq.names)) {
    taxa <- sub("^ {0,}>", "", taxa)
    seq.names <- getTaxaNames(taxa)
  }
  start <- c(start, length(X) + 1)
  for (i in 1:n) obj[[i]] <- unlist(strsplit(gsub(" ", "",
                                  X[(start[i] + 1):(start[i + 1] - 1)]), NULL))
  names(obj) <- seq.names
  obj <- lapply(obj, tolower)
  obj
}


#' Import and export sequence alignments
#'
#' These functions read and write sequence alignments.
#'
#' \code{write.phyDat} calls the function \code{\link[ape]{write.dna}} or
#' \code{\link[ape]{write.nexus.data}} and \code{read.phyDat} calls the function
#' \code{\link[ape]{read.dna}}, \code{read.aa} or \code{read.nexus.data}, so see
#' for more details over there.
#'
#' You may import data directly with \code{\link[ape]{read.dna}} or
#' \code{\link[ape]{read.nexus.data}} and convert the data to class phyDat.
#'
#' @param file a file name specified by either a variable of mode character, or
#' a double-quoted string.
#' @param format File format of the sequence alignment (see details).  Several
#' popular formats are supported: "phylip", "interleaved", "sequential",
#' "clustal", "fasta" or "nexus", or any unambiguous abbreviation of these.
#' @param type Type of sequences ("DNA", "AA", "CODON" or "USER").
#' @param ... further arguments passed to or from other methods.
#' @return \code{read.phyDat} returns an object of class phyDat,
#' \code{write.phyDat} write an alignment to a file.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{read.dna}}, \code{\link[ape]{read.GenBank}},
#' \code{\link[phangorn]{phyDat}}, \code{\link[seqinr]{read.alignment}}
#' @references % Anonymous. FASTA format description. %
#' \url{https://www.ncbi.nlm.nih.gov/blast/fasta.shtml} Felsenstein, J. (1993)
#' Phylip (Phylogeny Inference Package) version 3.5c. Department of Genetics,
#' University of Washington.
#' \url{https://evolution.genetics.washington.edu/phylip/phylip.html}
#' @examples
#' fdir <- system.file("extdata/trees", package = "phangorn")
#' primates <- read.phyDat(file.path(fdir, "primates.dna"),
#'                         format = "interleaved")
#' @keywords IO
#' @rdname read.phyDat
#' @export
read.phyDat <- function(file, format="phylip", type="DNA", ...){
  formats <- c("phylip", "nexus", "interleaved", "sequential", "fasta",
               "clustal")
  format <- match.arg(tolower(format), formats)

  if(format=="nexus") data <- read.nexus.data(file, ...)
  else {
    if(format=="phylip") format <- "sequential" #"interleaved"
    if (type == "DNA" || type == "CODON"){
      data <- read.dna(file, format, as.character = (format!="fasta"), ...)
    }
    if (type == "AA") data <- read.aa(file, format=format, ...)
    if (type == "USER"){
      if(format=="fasta") data <- read.fasta.user(file)
      else data <- read.dna(file, format, as.character = TRUE)
      extras <- match.call(expand.dots = FALSE)$...
      extras <- lapply(extras, eval)
      return(phyDat(data, type, levels=extras$levels,
                    ambiguity = extras$ambiguity, contrast = extras$contrast))
    }
  }
  if(is.list(data)){
    ll <- lengths(data)
    if(!all(ll == ll[[1]])) stop("sequences have different length")
  }
  phyDat(data, type, return.index = TRUE)
}


#' @param x An object of class \code{phyDat}.
#' @param colsep a character used to separate the columns (a single space by
#' default).
#' @param nbcol a numeric specifying the number of columns per row (-1 by
#' default); may be negative implying that the nucleotides are printed on a
#' single line.
#' @rdname read.phyDat
#' @export
write.phyDat <- function(x, file, format="phylip", colsep = "", nbcol=-1, ...){
  formats <- c("phylip", "nexus", "interleaved", "sequential", "fasta")
  format <- match.arg(tolower(format), formats)
  if(format=="nexus"){
    type <- attr(x, "type")
    if(type=="DNA") write.nexus.data(as.list(as.data.frame(x)), file,
                                     format = "dna",...)
    else write.nexus.data(as.list(as.data.frame(x)), file,
                          format = "protein", ...)
  }
  else{
    if(format=="phylip") format <- "interleaved"
    write.dna(as.character(x), file, format=format, colsep = colsep,
              nbcol=nbcol, ...)
  }
}


