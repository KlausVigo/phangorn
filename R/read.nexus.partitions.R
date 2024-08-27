read.nexus.charset <- function(file){
  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  X <- gsub("\\[(.*?)\\]", "", X)
  setStart <- grep("BEGIN SETS;", X, ignore.case = TRUE)
  if (length(setStart) == 0) return(NULL)
  setEnd <- grep("END;", X, ignore.case = TRUE)
  setEnd <- setEnd[setEnd > setStart][1]
  X <- gsub(";", "", X)
  tmp <- grep("Charset", X, ignore.case = TRUE)
  if(length(tmp) == 0) return(NULL)
  charset <- X[tmp]
  charset <- gsub("charset ", "", charset, ignore.case = TRUE)
  nam <- character(length(charset))
  cset <- character(length(charset))
  for(i in seq_along(charset)){
    tmp <- strsplit(charset[i], "=")[[1]]
    nam[i] <- trimws(tmp[1])
    cset[i] <- trimws(tmp[2])
  }
  res <- vector("list", length(nam))
  names(res) <- nam
  for(i in seq_along(cset)){
    tmp <- strsplit(cset[i], " ")[[1]]
    for (j in tmp){
      if(grepl("/", j)){
        y123 <- strsplit(j, "/")[[1]]
        y12 <-  as.numeric(strsplit(y123[1], "-")[[1]])
        res[[i]] <- c(res[[i]], seq(y12[1], y12[2], as.numeric(y123[2])))
      } else if(grepl("-", j)){
        y12 <-  as.numeric(strsplit(j, "-")[[1]])
        res[[i]] <- c(res[[i]], seq(y12[1], y12[2]))
      } else res[[i]] <- c(res[[i]], as.numeric(j))
    }
  }
  res
}


#' Function to import partitioned data from nexus files
#'
#' \code{read.nexus.partitions} reads in sequences in NEXUS format and splits
#' the data according to the charsets given in the SETS block.
#'
#' @param file a file name.
#' @param return either returns a list where each element is a 'phyDat' object
#' or an object of class 'multiphyDat'
#' @param \dots Further arguments passed to or from other methods.
#' @return a list where each element is a 'phyDat' object or an object of class
#' 'multiphyDat'.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{read.nexus.data}}, \code{\link{read.phyDat}}
#' @keywords cluster
#' @examples
#' tree <- rtree(10)
#' dat <- simSeq(tree, l=24)
#' fcat <- function(..., file = zz) cat(..., file=file, sep="", append=TRUE)
#' zz <- tempfile(pattern="file", tmpdir=tempdir(), fileext=".nex")
#' write.phyDat(dat, file=zz, format="nexus")
#' fcat("BEGIN SETS;\n")
#' fcat("  Charset codon1 = 1-12/3;\n")
#' fcat("  Charset codon2 = 2-12/3;\n")
#' fcat("  Charset codon3 = 3-12/3;\n")
#' fcat("  Charset range = 16-18;\n")
#' fcat("  Charset range2 = 13-15 19-21;\n")
#' fcat("  Charset singles = 22 23 24;\n")
#' fcat("END;\n")
#'
#' tmp <- read.nexus.partitions(zz)
#' tmp
#' unlink(zz)
#' @rdname read.nexus.partitions
#' @export
read.nexus.partitions <- function(file, return="list", ...){
  return <- match.arg(return, c("list", "multiphyDat"))
  dat <- read.phyDat(file, format="nexus", ...)
  genes <- read.nexus.charset(file)
  if(is.null(genes)) stop(paste(file, "does not contain Charset!"))
  seq <- lapply(genes, \(x, dat)dat[,x], dat)
  names(seq) <- names(genes)
  if(return=="multiphyDat" && requireNamespace("apex"))
    seq <- new("multiphyDat", seq = seq, add.gaps = FALSE)
  seq
}
