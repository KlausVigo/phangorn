#
# Data structures for ML and MP
#
fast.table <- function (data)
{
    if(!is.data.frame(data))
        data <- as.data.frame(data, stringsAsFactors = FALSE)
    da <- do.call("paste", c(data, sep = "\r"))
    ind <- !duplicated(da)
    levels <- da[ind]
    cat <- factor(da,levels = levels)
    nl <- length(levels(cat))
    bin <- (as.integer(cat) - 1)
    pd <- nl
    bin <- bin[!is.na(bin)]
    if (length(bin)) bin <- bin + 1
    y <- tabulate(bin, pd)
    result <- list(index = bin, weights = y, data = data[ind,])
    result
}


phyDat.default <- function (data, levels = NULL, return.index = TRUE,
                    contrast = NULL, ambiguity = "?", compress=TRUE, ...)
{
    if (is.matrix(data))
        nam <- row.names(data)
    else nam <- names(data)
    if(is.null(nam))stop("data object must contain taxa names")
    if(inherits(data, "character") | inherits(data, "numeric"))
        data <- as.matrix(data)
    if (inherits(data,"DNAbin"))
        data <- as.character(data)
    if (is.matrix(data))
        data <- as.data.frame(t(data), stringsAsFactors = FALSE)
    # new 4.4.2016 bug fix (reported by Eli Levy Karin)
    #    if (is.vector(data) && !is.list(data))data = as.data.frame(data,
    #    stringsAsFactors = FALSE)
    else data <- as.data.frame(data, stringsAsFactors = FALSE)
    #    data = data.frame(as.matrix(data), stringsAsFactors = FALSE)

    if(length(data[[1]])==1) compress <- FALSE
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
    q <- length(data)
    p <- length(data[[1]])
    tmp <- vector("list", q)
    if (!is.null(contrast)) {
        levels <- colnames(contrast)
        all.levels <- rownames(contrast)
        rownames(contrast) <- NULL
    }
    else {
        if (is.null(levels))
            stop("Either argument levels or contrast has to be supplied")
        l <- length(levels)
        contrast <- diag(l)
        all.levels <- levels
        if (!is.null(ambiguity)) {
            all.levels <- c(all.levels, ambiguity)
            k <- length(ambiguity)
            if (k > 0)
                contrast <- rbind(contrast, matrix(1, k, l))
        }
    }
#    row.names(data) = as.character(1:p)
#    data = na.omit(data)
#    rn = as.numeric(rownames(data))

    d <- dim(data)
    att <- attributes(data)
    data <- lapply(data, match, all.levels)
#    data <- match(unlist(data), all.levels)
#    attr(data, "dim") <- d
#    data <- as.data.frame(data, stringsAsFactors=FALSE)
    attributes(data) <- att

    row.names(data) <- as.character(1:p)
    data <- na.omit(data)
    aaa <- match(index, attr(data, "na.action"))

    if(!is.null(attr(data, "na.action"))) warning("Found unknown characters (not supplied in levels). Deleted sites with with unknown states.")

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
    attr(data, "nc") <- length(levels)
    if (return.index)
        attr(data, "index") <- index
    attr(data, "levels") <- levels
    attr(data, "allLevels") <- all.levels
    attr(data, "type") <- "USER"
    attr(data, "contrast") <- contrast
    class(data) <- "phyDat"
    data
}


phyDat.DNA <- function (data, return.index = TRUE)
{
  if (is.matrix(data))
    nam <- row.names(data)
  else nam <- names(data)
  if (inherits(data,"DNAbin"))
    data <- as.character(data)
  if(inherits(data, "character")) data <- as.matrix(data)
  if (is.matrix(data))
    data <- as.data.frame(t(data), stringsAsFactors = FALSE)
  else data <- as.data.frame(data, stringsAsFactors = FALSE)

  data <- data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)

  ac <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y",
          "k", "v", "h", "d", "b", "n", "?", "-")
  AC <- matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1),
                 c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1),
                 c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1)),
               18, 4, dimnames = list(NULL, c("a", "c", "g", "t")))

  compress <- TRUE
  if(length(data[[1]])==1) compress <- FALSE
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
  q <- length(data)
  p <- length(data[[1]])
  d <- dim(data)
  att <- attributes(data)

  data <- lapply(data, match, ac)
#  data <- match(unlist(data), ac)
#  attr(data, "dim") <- d
#  data <- as.data.frame(data, stringsAsFactors=FALSE)
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
  attr(data, "nc") <- 4
  if (return.index)
    attr(data, "index") <- index
  attr(data, "levels") <- c("a", "c", "g", "t")
  attr(data, "allLevels") <- ac
  attr(data, "type") <- "DNA"
  attr(data, "contrast") <- AC
  class(data) <- "phyDat"
  data
}


phyDat.AA <- function (data, return.index = TRUE)
{
    if(is.matrix(data)) nam <- row.names(data)
    else nam <- names(data)
    # AAbin
    if (inherits(data,"AAbin"))
        data <- as.character(data)
    if(inherits(data, "character")) data <- as.matrix(data)
    if (is.matrix(data))
        data <- as.data.frame(t(data), stringsAsFactors = FALSE)
    else data <- as.data.frame(data, stringsAsFactors = FALSE)

    data <- data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)

    aa <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i",
        "l", "k", "m", "f", "p", "s", "t", "w", "y", "v")
    aa2 <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i",
        "l", "k", "m", "f", "p", "s", "t", "w", "y", "v", "b",
        "z", "x", "-", "?")
    AA <- diag(20)
    AA <- rbind(AA, matrix(0, 5, 20))
    AA[21, 3] <- AA[21, 4] <- 1 # Aspartate or Asparagine
    AA[22, 6] <- AA[22, 7] <- 1 #
    AA[23:25, ] <- 1
    dimnames(AA) <- list(aa2, aa)
    compress <- TRUE
    if(length(data[[1]])==1) compress <- FALSE
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
    q <- length(data)
    p <- length(data[[1]])
    tmp <- vector("list", q)

    d <- dim(data)
    att <- attributes(data)

    data <- lapply(data, match, aa2)
    # data <- match(unlist(data), aa2)
    # attr(data, "dim") <- d
    # data <- as.data.frame(data, stringsAsFactors=FALSE)
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
    attr(data, "nc") <- 20
    if (return.index)
        attr(data, "index") <- index
    attr(data, "levels") <- aa
    attr(data, "allLevels") <- aa2
    attr(data, "type") <- "AA"
    attr(data, "contrast") <- AA
    class(data) <- "phyDat"
    data
}


phyDat.codon <- function (data, return.index = TRUE, ambiguity = "---",
                          NA_as_ambiguous=TRUE)
{
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

    splseq <- function (seq, frame = 0)
    {
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
    codon <- c("aaa", "aac", "aag", "aat", "aca", "acc", "acg", "act",
      "aga", "agc", "agg", "agt", "ata", "atc", "atg", "att",
      "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct", "cga",
      "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt", "gaa", "gac",
      "gag", "gat", "gca", "gcc", "gcg", "gct", "gga", "ggc", "ggg",
      "ggt", "gta", "gtc", "gtg", "gtt", "tac", "tat",
      "tca", "tcc", "tcg", "tct", "tgc", "tgg", "tgt", "tta",
      "ttc", "ttg", "ttt")
# ohne Stopcodons "taa", "tag", "tga",

    CODON <- diag(61)

    if(NA_as_ambiguous){
        ambiguity <- unique(c("---", ambiguity))
    }
    if(ambiguity!=""){
        codon_amb <- c(codon, ambiguity)
        CODON <- rbind(CODON, matrix(1, length(ambiguity), 61))
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
    attr(data, "nc") <- 61
    if (return.index)
        attr(data, "index") <- index
    attr(data, "levels") <- codon
    attr(data, "allLevels") <- codon_amb
    attr(data, "type") <- "CODON"
    attr(data, "contrast") <- CODON
    class(data) <- "phyDat"
    data
}


#' Conversion among Sequence Formats
#'
#' These functions transform several DNA formats into the \code{phyDat} format.
#' \code{allSitePattern} generates an alignment of all possible site patterns.
#'
#' If \code{type} "USER" a vector has to be give to \code{levels}. For example
#' c("a", "c", "g", "t", "-") would create a data object that can be used in
#' phylogenetic analysis with gaps as fifth state.  There is a more detailed
#' example for specifying "USER" defined data formats in the vignette
#' "phangorn-specials".
#'
#' \code{allSitePattern} returns all possible site patterns and can be useful
#' in simulation studies. For further details see the vignette
#' phangorn-specials.
#'
#' \code{write.phyDat} calls the function write.dna or write.nexus.data and
#' \code{read.phyDat} calls the function \code{read.dna}, \code{read.aa} or
#' \code{read.nexus.data} see for more details over there.
#'
#' You may import data directly with \code{\link[ape]{read.dna}} or
#' \code{\link[ape]{read.nexus.data}} and convert the data to class phyDat.
#'
#' The generic function \code{c} can be used to to combine sequences and
#' \code{unique} to get all unique sequences or unique haplotypes.
#'
#' \code{acgt2ry} converts a \code{phyDat} object of nucleotides into an binary
#' ry-coded dataset.
#'
#' @aliases
#' as.phyDat.character as.phyDat.data.frame as.phyDat.matrix
#' as.MultipleAlignment as.MultipleAlignment.phyDat cbind.phyDat c.phyDat
#' acgt2ry removeUndeterminedSites phyDat2MultipleAlignment
#' @param data An object containing sequences.
#' @param x An object containing sequences.
#' @param type Type of sequences ("DNA", "AA", "CODON" or "USER").
#' @param levels Level attributes.
#' @param return.index If TRUE returns a index of the site patterns.
#' @param file A file name.
#' @param format File format of the sequence alignment (see details).  Several
#' popular formats are supported: "phylip", "interleaved", "sequential",
#' "clustal", "fasta" or "nexus", or any unambiguous abbreviation of these.
#' @param colsep a character used to separate the columns (a single space by
#' default).
#' @param nbcol a numeric specifying the number of columns per row (-1 by
#' default); may be negative implying that the nucleotides are printed on a
#' single line.
#' @param n Number of sequences.
#' @param names Names of sequences.
#' @param subset a subset of taxa.
#' @param select a subset of characters.
#' @param site.pattern select site pattern or sites.
#' @param allLevels return original data.
#' @param obj as object of class phyDat
#' @param freq logical, if 'TRUE', frequencies or counts are returned otherwise
#' proportions
#' @param all all a logical; if all = TRUE, all counts of bases, ambiguous
#' codes, missing data, and alignment gaps are returned as defined in the
#' contrast.
#' @param drop.unused.levels logical, drop unused levels
#' @param incomparables for compatibility with unique.
#' @param identical if TRUE (default) sequences have to be identical, if FALSE
#' sequences are considered duplicates if distance between sequences is zero
#' (happens frequently with ambiguous sites).
#' @param ambiguity character for ambiguous character and no contrast is
#' provided.
#' @param codonstart an integer giving where to start the translation. This
#' should be 1, 2, or 3, but larger values are accepted and have for effect to
#' start the translation further within the sequence.
#' @param ... further arguments passed to or from other methods.
#' @return The functions return an object of class \code{phyDat}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{DNAbin}}, \code{\link{as.DNAbin}},
#' \code{\link{read.dna}}, \code{\link{read.aa}}, \code{\link{read.nexus.data}}
#' and the chapter 1 in the \code{vignette("phangorn-specials",
#' package="phangorn")} and the example of \code{\link{pmlMix}} for the use of
#' \code{allSitePattern}
#' @keywords cluster
#' @examples
#'
#' data(Laurasiatherian)
#' class(Laurasiatherian)
#' Laurasiatherian
#' # base frequencies
#' baseFreq(Laurasiatherian)
#' baseFreq(Laurasiatherian, all=TRUE)
#' baseFreq(Laurasiatherian, freq=TRUE)
#' # subsetting phyDat objects
#' # the first 5 sequences
#' subset(Laurasiatherian, subset=1:5)
#' # the first 5 characters
#' subset(Laurasiatherian, select=1:5, site.pattern = FALSE)
#' # the first 5 site patterns (often more than 5 characters)
#' subset(Laurasiatherian, select=1:5, site.pattern = TRUE)
#' # transform into old ape format
#' LauraChar <- as.character(Laurasiatherian)
#' # and back
#' Laura <- phyDat(LauraChar)
#' all.equal(Laurasiatherian, Laura)
#' # Compute all possible site patterns
#' # for nucleotides there $4 ^ (number of tips)$ patterns
#' allSitePattern(5)
#'
#' @rdname phyDat
#' @export
phyDat <- function (data, type="DNA", levels=NULL, return.index = TRUE,...)
{
    if (inherits(data,"DNAbin")) type <- "DNA"
    pt <- match.arg(type, c("DNA", "AA", "CODON", "USER"))
    if(pt=="DNA") dat <- phyDat.DNA(data, return.index=return.index,...)
    if(pt=="AA") dat <- phyDat.AA(data, return.index=return.index, ...)
    if(pt=="CODON") dat <- phyDat.codon(data, return.index=return.index, ...)
    if(pt=="USER") dat <- phyDat.default(data, levels = levels,
                                         return.index=return.index, ...)
    dat
}


#' @rdname phyDat
#' @export
dna2codon <- function(x, codonstart=1, ambiguity="---", ...){
    if(!inherits(x, "phyDat"))stop("x needs to be of class phyDat!")
    if(attr(x, "type")=="AA")stop("x needs to be a nucleotide sequence!")

    if(codonstart>1){
        del <- -seq_len(codonstart)
        x <- subset(x, select=del, site.pattern=FALSE)
    }
    n_sites <- sum(attr(x,"weight"))
    if( (n_sites %% 3) ){
      keep <- seq_len( (n_sites %/% 3) * 3 )
      x <- subset(x, select=keep, site.pattern=FALSE)
    }
    phyDat.codon(as.character(x), ambiguity=ambiguity, ...)
}


#' @rdname phyDat
#' @export
codon2dna <- function(x){
    if(!inherits(x, "phyDat"))stop("x needs to be of class phyDat!")
    phyDat.DNA(as.character(x))
}


#' @rdname phyDat
#' @export
as.phyDat <- function (x, ...){
    if (inherits(x,"phyDat")) return(x)
    UseMethod("as.phyDat")
}


#' @rdname phyDat
#' @method as.phyDat factor
#' @export
as.phyDat.factor <- function(x, ...){
    nam <- names(x)
    lev <- levels(x)
    x <- as.character(x)
    names(x) <- nam
    phyDat(x, type="USER", levels = lev, ...)
}


#' @rdname phyDat
#' @method as.phyDat DNAbin
#' @export
as.phyDat.DNAbin <- function(x,...) phyDat.DNA(x,...)


#' @rdname phyDat
#' @method as.phyDat alignment
#' @export
as.phyDat.alignment <- function (x, type="DNA",...)
{
    x$seq <- tolower(x$seq)
    data <- sapply(x$seq, strsplit, "")
    names(data) <- x$nam
    if(type=="DNA") dat <- phyDat.DNA(data,...)
    if(type=="AA") dat <- phyDat.AA(data, ...)
    if(type=="CODON") dat <- phyDat.codon(data, ...)
    if(type=="USER") dat <- phyDat.default(data, ...)
    dat
}


#as.alignment.phyDat <- function(x, ...) as.alignment(as.character(x))
#' @rdname phyDat
#' @export
phyDat2alignment <-  function(x){
    z <- as.character(x)
    nam <- rownames(z)
    type <- attr(x, "type")
    seq <- switch(type,
                  DNA = tolower(apply(z, 1, paste, collapse="")),
                  AA = toupper(apply(z, 1, paste, collapse="")))
    names(seq) <- NULL
    res <- list(nb=length(seq), nam=nam, seq=seq, com=NA)
    class(res) <- "alignment"
    res
}


#' @rdname phyDat
#' @method as.phyDat MultipleAlignment
#' @export
as.phyDat.MultipleAlignment <- function(x, ...){
    if (requireNamespace('Biostrings')){
    if(inherits(x, "DNAMultipleAlignment"))
        res <- phyDat.DNA(Biostrings::as.matrix(x))
    if(inherits(x, "RNAMultipleAlignment"))
        res <- phyDat.DNA(Biostrings::as.matrix(x))
    if(inherits(x, "AAMultipleAlignment"))
        res <- phyDat.AA(Biostrings::as.matrix(x))
    return(res)
    }
    return(NULL)
}


# @rdname phyDat
#' @export
as.MultipleAlignment <- function (x, ...){
    if (inherits(x,"MultipleAlignment")) return(x)
    UseMethod("as.MultipleAlignment")
}


#' @rdname phyDat
#' @export
as.MultipleAlignment.phyDat <- function(x, ...){
    if (requireNamespace('Biostrings')){
    z <- as.character(x)
    nam <- rownames(z)
    type <- attr(x, "type")
    seq <- switch(type,
                  DNA = tolower(apply(z, 1, paste, collapse="")),
                  AA = toupper(apply(z, 1, paste, collapse="")))
    if(type=="DNA") return(Biostrings::DNAMultipleAlignment(seq))
    if(type=="AA") return(Biostrings::AAMultipleAlignment(seq))
    }
    return(NULL)
}


#' @export
phyDat2MultipleAlignment <- as.MultipleAlignment.phyDat


#' @export
as.phyDat.matrix <- function (x, ...) phyDat(data=x, ...)


#' @export
as.phyDat.character <- function (x, ...) phyDat(data=x, ...)

# @rdname phyDat
#' @export
as.phyDat.data.frame <- function (x, ...) phyDat(data=x, ...)


#' @rdname phyDat
#' @export
acgt2ry <- function(obj){
   ac <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y",
        "k", "v", "h", "d", "b", "n", "?", "-")
   AC <- matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
        0, 1, 1, 1), c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1,
        0, 1, 1, 1, 1), c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1,
        0, 1, 1, 1, 1, 1), c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1,
        0, 1, 1, 1, 1, 1, 1)), 18, 4, dimnames = list(NULL, c("a",
        "c", "g", "t")))
   ry <- AC[c(7,10),]
   RY <- AC %*% t(ry)
   RY[RY==2] <- 1
   dimnames(RY) <- list(NULL, c("r", "y"))
   attr(obj, "levels") <- c("r", "y")
   attr(obj, "nc") <- 2
   attr(obj, "type") <- "USER"
   attr(obj, "contrast") <- RY
   obj <- phyDat.default(as.character(obj, allLevels=FALSE),
                         levels = c("r", "y"), ambiguity = NULL)
   obj
}


#' @rdname phyDat
#' @export
# replace as.character.phyDat weniger Zeilen, works also for codons
as.character.phyDat <- function (x, allLevels=TRUE, ...)
{
    nr <- attr(x, "nr")
    nc <- attr(x, "nc")
    type <- attr(x, "type")
    labels <- attr(x, "allLevels")

    if (!is.null(attr(x, "index"))) {
        index <- attr(x, "index")
        if (is.data.frame(index))
            index <- index[, 1]
    }
    else index <- rep(1:nr, attr(x, "weight"))
    if (type == "USER") {
        #levels in acgt2ry
        if(!allLevels){
            tmp <- attr(x, "levels")
            contrast <- attr(x, "contrast") # contrast=AC
            contrast[contrast>0] <- 1
            ind <- which(rowSums(contrast)==1)
            contrast[rowSums(contrast)>1,] <- 0
            labels <- rep(NA, length(attr(x, "allLevels")))
            labels[ind] <- tmp[contrast%*%c(seq_along(tmp))]
        }
    }
    if(type == "AA") labels <- toupper(labels)
    if(type == "CODON"){
        nr <- length(index)
        result <- matrix(NA, nrow = length(x), ncol = 3L*nr)
        labels <- strsplit(labels, "")
        for (i in seq_along(x)) result[i, ] <- unlist(labels[ x[[i]][index] ])
    }
    else {
        result <- matrix(NA, nrow = length(x), ncol = nr)
        for (i in seq_along(x)) result[i, ] <- labels[x[[i]]]
        result <- result[, index, drop = FALSE]
    }
    rownames(result) <- names(x)
    result
}


#' @rdname phyDat
#' @export
as.data.frame.phyDat <- function(x, ...){
  nr <- attr(x, "nr")
  nc <- attr(x, "nc")
  labels <- attr(x, "allLevels")
  if(attr(x, "type") == "AA") labels <- toupper(labels)
  result <- vector("list", length(x))
  if (is.null(attr(x, "index")))
    index <- rep(1:nr, attr(x, "weight"))
  else {
    index <- attr(x, "index")
    if (is.data.frame(index))
      index <- index[, 1]
  }
  for (i in seq_along(x)) result[[i]] <- labels[x[[i]][index]]
  attr(result, "names") <- names(x)
  attr(result, "row.names") <- seq_along(index)
  attr(result, "class") <- "data.frame"
  result
}


#as.DNAbin.phyDat <- function(x,...) {
#   if(attr(x, "type")=="DNA") return(as.DNAbin(as.character(x, ...)))
#   else stop("x must be a nucleotide sequence")
#}

# quite abit faster
#' @rdname phyDat
#' @export
as.DNAbin.phyDat <- function (x, ...)
{
    if(attr(x, "type")=="DNA"){

    nr <- attr(x, "nr")
    ac <- attr(x, "allLevels")
    result <- matrix(as.raw(0), nrow = length(x), ncol = nr)
    # from ape ._cs_
    cs <- c("a", "g", "c", "t", "r", "m", "w", "s", "k", "y", "v", "h",
      "d", "b", "n", "-", "?")
    # from ape ._bs_
    bs <- as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48, 224, 176, 208,
                   112, 240, 4, 2))
    ord <- match(ac, cs)
    ord[5] <- 4

    for (i in seq_along(x)){
        ind <- ord[x[[i]]]
        result[i,] <- bs[ind]
    }
    if (is.null(attr(x, "index")))
        index <- rep(1:nr, attr(x, "weight"))
    else {
        index <- attr(x, "index")
        if (is.data.frame(index))
            index <- index[, 1]
    }
    result <- result[, index, drop = FALSE]
    rownames(result) <- names(x)
    class(result) <- "DNAbin"
    return(result)
    }
    else stop("x must be a nucleotide sequence")
}


#' @rdname phyDat
#' @export
as.AAbin.phyDat <- function(x,...) {
   if(attr(x, "type")=="AA") return(as.AAbin(as.character(x, ...)))
   else stop("x must be a amino acid sequence")
}


#' @export
print.phyDat <- function (x, ...)
{
    cat(length(x), "sequences with",sum(attr(x,"weight")), "character and",
        attr(x,"nr"),"different site patterns.\n")
    cat("The states are",attr(x,"levels"), "\n")
}


# in C++ to replace aggregate or use of distinct from dplyr / data.table
aggr <- function(weight, ind){
    res <- numeric(max(ind))
    for(i in seq_along(weight))
        res[ind[i]] <- res[ind[i]] + weight[i]
    res
}


# data has to be a data.frame in cbind.phyDat
fast.table2 <- function (data)
{
    if(!is.data.frame(data))
        data <- as.data.frame(data, stringsAsFactors = FALSE)
    da <- do.call("paste", data)
    ind <- !duplicated(da)
    levels <- da[ind]
    cat <- factor(da,levels = levels)
    nl <- length(levels)
    bin <- (as.integer(cat) - 1L)
    bin <- bin[!is.na(bin)]
    if (length(bin)) bin <- bin + 1L
    result <- list(index = bin, pos = ind)
    result
}


# @rdname phyDat
#' @export cbind.phyDat
#' @export
cbind.phyDat <- function(..., gaps="-", compress=TRUE){
    object <- as.list(substitute(list(...)))[-1]
    x <- list(...)
    n <- length(x)
    if (n == 1)
        return(x[[1]])
    type <- attr(x[[1]], "type")
    nr <- numeric(n)

    ATTR <- attributes(x[[1]])

    nr[1] <- sum(attr(x[[1]], "weight"))
    levels <- attr(x[[1]], "levels")
    allLevels <- attr(x[[1]], "allLevels")
    gapsInd <- match(gaps, allLevels)
    snames <- vector("list", n)  # names(x[[1]])
    vec <- numeric(n+1)
    wvec <- numeric(n+1)
    objNames <- as.character(object)
    if(any(duplicated(objNames))) objNames <- paste0(objNames, 1:n)
    #    tmp <- as.character(x[[1]])

    for(i in 1:n){
        snames[[i]] <- names(x[[i]])
        nr[i] <- sum(attr(x[[i]], "weight"))
        vec[i+1] <- attr(x[[i]], "nr")
        wvec[i+1] <- sum(attr(x[[i]], "weight"))
    }
    vec <- cumsum(vec)
    wvec <- cumsum(wvec)
    snames <- unique(unlist(snames))
    weight <- numeric(vec[n+1])

    index <- numeric(wvec[n+1])

    ATTR$names <- snames
    ATTR$nr <- vec[n+1]

    tmp <- matrix(gapsInd, vec[n+1], length(snames),
                  dimnames = list(NULL, snames))
    tmp <- as.data.frame(tmp)
    add.index <- TRUE
    for(i in 1:n){
        nam <- names(x[[i]])
        tmp[(vec[i]+1):vec[i+1], nam] <- x[[i]][nam]
        weight[(vec[i]+1):vec[i+1]] <- attr(x[[i]], "weight")
    }
    if(compress){
        ddd <- fast.table2(tmp)
        tmp <- tmp[ddd$pos,]
        weight <- aggregate(weight, by=list(ddd$index), FUN=sum)$x
    }
    if(any(sapply(x, function(x)is.null(attr(x, "index"))))) add.index <- FALSE
    if(add.index  & compress){
        for(i in 1:n){
            tmp2 <- attr(x[[i]], "index")
            if(!is.null(tmp2)){
                if(is.data.frame(tmp2))index[(wvec[i]+1):wvec[i+1]] <-
                        ddd$index[(vec[i]+1):vec[i+1]][tmp2[,1]]
                else index[(wvec[i]+1):wvec[i+1]] <-
                        ddd$index[(vec[i]+1):vec[i+1]][tmp2]
            }
            else add.index <- FALSE
        }
    }
    if(add.index)ATTR$index <- data.frame(index = index, genes=rep(objNames, nr))
    ATTR$weight <- weight
    ATTR$nr <- length(weight)
    attributes(tmp) <- ATTR
    tmp
}

# @rdname phyDat
#' @export c.phyDat
#' @export
c.phyDat <- cbind.phyDat


#' @rdname phyDat
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



read.fasta.user <- function (file, skip = 0, nlines = 0,
                             comment.char = "#", seq.names = NULL)
{
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
    for (i in 1:n) obj[[i]] <- unlist(strsplit(gsub(" ",
                        "", X[(start[i] + 1):(start[i + 1] - 1)]), NULL))
    names(obj) <- seq.names
    obj <- lapply(obj, tolower)
    obj
}


#' @rdname phyDat
#' @export
read.phyDat <- function(file, format="phylip", type="DNA", ...){

    formats <- c("phylip", "nexus", "interleaved", "sequential", "fasta",
                 "clustal")
    format <- match.arg(tolower(format), formats)

    if(format=="nexus") data <- read.nexus.data(file, ...)
    else {
        if(format=="phylip") format <- "sequential" #"interleaved"
        if (type == "DNA" || type == "CODON"){
            data <- read.dna(file, format, as.character = TRUE, ...)
        }
        if (type == "AA") data <- read.aa(file, format=format, ...)
        if (type == "USER"){
            if(format=="fasta")
                data <- read.fasta.user(file)
            else data <- read.dna(file, format, as.character = TRUE)
            extras <- match.call(expand.dots = FALSE)$...
            extras <- lapply(extras, eval)
            return(phyDat(data, type, levels=extras$levels,
                    ambiguity = extras$ambiguity, contrast = extras$contrast))
        }
        # raus
    }
    if(is.list(data)){
        ll <- lengths(data)
        if(!all(ll == ll[[1]])) stop("sequences have different length")
    }
    phyDat(data, type, return.index = TRUE)
}


#' @rdname phyDat
#' @export
baseFreq <- function(obj, freq=FALSE, all=FALSE, drop.unused.levels = FALSE){
    if (!inherits(obj,"phyDat"))
        stop("data must be of class phyDat")
    labels <- attr(obj, "allLevels")
    weight <- attr(obj,"weight")
    n <- length(obj)
    res <- numeric(length(labels))
    D <- diag(length(labels))
    for(i in 1:n)res <- res + colSums(D[obj[[i]],, drop=FALSE]*weight)
    names(res) <- labels
    if(!all) res <- res[attr(obj, "levels")]
    if(!freq)res <- res/sum(res)
    if(drop.unused.levels) return(res[res>0])
    res
}


phylo <- function(edge, tip, edge.length=NULL){
    res <- list(edge=edge, tip.label=tip, edge.length=edge.length)
    class(res) <- "phylo"
    res
    }


getCols <- function (data, cols)
{
    attrib <- attributes(data)
    attr(data, "class") <- "list"
    data <- data[cols]
    if (is.character(cols))
        attrib$names <- cols
    else attrib$names <- attrib$names[cols]
    attributes(data) <- attrib
    attr(data, "class") <- "phyDat"
    data
}


# allows negative indexing subset(dat,,-c(3:5))
getRows <- function (data, rows, site.pattern = TRUE)
{
  index <- attr(data, "index")
  if(is.data.frame(index))index <- index[,1]
  if(!site.pattern){ # & all(rows>0)

    weight <- tabulate(index[rows])
    ind <- which(weight>0)
# update index
    new_index <- integer(length(weight))
    new_index[ind] <- seq_along(ind)
    attr(data, "index") <- new_index[index[rows]]

    rows <- ind   # rows[ind]
    weight <- weight[ind]
  }
  for (i in seq_along(data)){
    if(is.matrix(data[[i]]))data[[i]] <- data[[i]][rows,]
    else data[[i]] <- data[[i]][rows]
  }
  attr(data, "weight") <- attr(data, "weight")[rows]
  if(!site.pattern) attr(data, "weight") <- weight
  attr(data, "nr") <- length(attr(data, "weight"))
  if(site.pattern)attr(data, "index") <- NULL
  data
}


#' @rdname phyDat
#' @method subset phyDat
#' @export
subset.phyDat <- function (x, subset, select, site.pattern = TRUE,...)
{

    if (!missing(subset)) x <- getCols(x, subset)
    if (!missing(select)){
         if(any(is.na(select))) return(NULL)
         x <- getRows(x, select, site.pattern=site.pattern)
    }
    x
}

## Needs testing that it is not used e.g. prepareDataFitch returns no class
#' @param i,j	indices of the rows and/or columns to select or to drop. They
#' may be numeric, logical, or character (in the same way than for standard R
#' objects).
#' @param drop for compatibility with the generic (unused).
#' @rdname phyDat
#' @export
"[.phyDat" <- function(x, i, j, ..., drop=FALSE){
   subset(x, subset = i, select = j, site.pattern=FALSE)
}


#' @rdname phangorn-internal
#' @export
map_duplicates <-  function(x, dist=TRUE, ...){
    labels <- names(x)
    if(dist){
        y <- as.matrix(dist.hamming(x, FALSE))
        l <- nrow(y)
        z <- character(l)
        for(i in seq_len(l)) z[i] <- paste( round(y[i, ] ,8), collapse="_")
        ind <- duplicated(z)
    }
    else ind <- duplicated(x)
    res <- NULL
    if(any(ind)){
        if(dist) ind2 <- match(z[ind], z)
        else ind2 <- match(x[ind], x)
        res <- data.frame(duplicates=labels[ind], where=labels[ind2],
                          stringsAsFactors = FALSE)
    }
    res
}


#duplicated_phyDat <- function(x, ...){
#    tmp <- map_duplicates(x)[,1]
#    getCols(x, setdiff(names(x), tmp))
#}


#' @rdname phyDat
#' @method unique phyDat
#' @export
unique.phyDat <- function(x, incomparables=FALSE, identical=TRUE, ...){
    if(identical) return(getCols(x, !duplicated(x)))
    tmp <- map_duplicates(x)[,1]
    getCols(x, setdiff(names(x), tmp))
}


#' @rdname phyDat
#' @export
removeUndeterminedSites <- function(x, ...){
# , use.contrast=TRUE, undetermined=c("?", "n", "-")
    nc <- attr(x, "nc")
    nr <- attr(x, "nr")
    contrast <- attr(x, "contrast")
#    if(use.contrast)
    ind <- which( (contrast %*% rep(1, nc)) == nc )
#    else ind <- sort(match(undetermined, attr(x, "allLevels")))
    tmp <- x[[1]] %in% ind
    for(i in 2:length(x)) tmp <- tmp & (x[[i]] %in% ind)
    if(any(tmp)) x <- getRows(x, (1:nr)[!tmp]) #getRows(x, -which(tmp))
    x
}


removeParsUninfoSites <- function(data){
    nr <- attr(data, "nr")
    pis <- parsinfo(data)
    if (length(pis) > 0){
        p0 <- sum(attr(data, "weight")[pis[, 1]] * pis[, 2])
        data <- getRows(data, c(1:nr)[-pis[, 1]], TRUE)
    }
    else p0 <- 0
    if(length(attr(data, "p0"))) p0 <- p0 + attr(data, "p0")
    attr(data, "p0") <-  p0
    data
}


#' @rdname phyDat
#' @export
allSitePattern <- function(n,levels=c("a","c","g","t"), names=NULL){
    l <- length(levels)
    X <- vector("list", n)
    if(is.null(names))names(X) <- paste0("t", 1:n)
    else names(X) <- names
    for(i in 1:n)
        X[[i]] <- rep(rep(levels, each=l^(i-1)),l^(n-i))
    X <- as.data.frame(X)
    phyDat.default(X, levels, compress=FALSE, return.index=FALSE)
}


constSitePattern <- function(n,levels=c("a","c","g","t"), names=NULL){
    l <- length(levels)
    X <- matrix(0, l,n)
    X <- matrix(rep(levels, each=n), n, l)
    if(is.null(names))rownames(X) <- paste0("t", 1:n)
    else rownames(X) <- names
    phyDat.default(X, levels)
}



#' Read Amino Acid Sequences in a File
#'
#' This function reads amino acid sequences in a file, and returns a matrix
#' list of DNA sequences with the names of the taxa read in the file as row
#' names.
#'
#'
#' @param file a file name specified by either a variable of mode character, or
#' a double-quoted string.
#' @param format a character string specifying the format of the DNA sequences.
#' Three choices are possible: \code{"interleaved"}, \code{"sequential"}, or
#' \code{"fasta"}, or any unambiguous abbreviation of these.
#' @param skip the number of lines of the input file to skip before beginning
#' to read data.
#' @param nlines the number of lines to be read (by default the file is read
#' until its end).
#' @param comment.char a single character, the remaining of the line after this
#' character is ignored.
#' @param seq.names the names to give to each sequence; by default the names
#' read in the file are used.
#' @return a matrix of amino acid sequences.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{read.dna}}, \code{\link[ape]{read.GenBank}},
#' \code{\link[phangorn]{phyDat}}, \code{\link[seqinr]{read.alignment}}
#' @references % Anonymous. FASTA format description. %
#' \url{https://www.ncbi.nlm.nih.gov/blast/fasta.shtml} Felsenstein, J. (1993)
#' Phylip (Phylogeny Inference Package) version 3.5c. Department of Genetics,
#' University of Washington.
#' \url{http://evolution.genetics.washington.edu/phylip/phylip.html}
#' @keywords IO
#' @export read.aa
read.aa <- function (file, format = "interleaved", skip = 0, nlines = 0,
    comment.char = "#", seq.names = NULL)
{
    getTaxaNames <- function(x) {
        x <- sub("^ +", "", x)
        x <- sub(" +$", "", x)
        x <- sub("^['\"]", "", x)
        x <- sub("['\"]$", "", x)
        x
    }
    format <- match.arg(format, c("interleaved", "sequential", "fasta"))
    phylip <- if (format %in% c("interleaved", "sequential"))
        TRUE
    else FALSE


    if (format == "fasta") {
#        obj <- read.FASTA.AA(file)
        obj <- read.FASTA(file, type = "AA")
        return(obj)
    }
    X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE,
        skip = skip, nlines = nlines, comment.char = comment.char)

    if (phylip) {
        fl <- X[1]
        oop <- options(warn = -1)
        fl.num <- as.numeric(unlist(strsplit(gsub("^ +", "", fl), " +")))
        options(oop)
        if (all(is.na(fl.num)))
            stop("the first line of the file must contain the dimensions of the data")
        if (length(fl.num) != 2)
            stop("the first line of the file must contain TWO numbers")
        else {
            n <- fl.num[1]
            s <- fl.num[2]
        }
        X <- X[-1]
        obj <- vector("character", n * s)
        dim(obj) <- c(n, s)
    }
    if (format == "interleaved") {
        fl <- X[1]
        fl <- unlist(strsplit(fl, NULL))
        bases <- grep("[-AaRrNnDdCcQqEeGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzXx?]", fl)
        z <- diff(bases)
        for (i in seq_along(z)) if (all(z[i:(i + 8)] == 1))
            break
        start.seq <- bases[i]
        if (is.null(seq.names))
            seq.names <- getTaxaNames(substr(X[1:n], 1, start.seq - 1))
        X[1:n] <- substr(X[1:n], start.seq, nchar(X[1:n]))
        X <- gsub(" ", "", X)
        nl <- length(X)
        for (i in 1:n) obj[i, ] <- unlist(strsplit(X[seq(i, nl, n)], NULL))
    }
    if (format == "sequential") {
        fl <- X[1]
        taxa <- character(n)
        j <- 1
        for (i in 1:n) {
            bases <- grep("[-AaRrNnDdCcQqEeGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzXx?]",
                unlist(strsplit(X[j], NULL)))
            z <- diff(bases)
            for (k in seq_along(z)) if (all(z[k:(k + 8)] == 1))
                break
            start.seq <- bases[k]
            taxa[i] <- substr(X[j], 1, start.seq - 1)
            sequ <- substr(X[j], start.seq, nchar(X[j]))
            sequ <- gsub(" ", "", sequ)
            j <- j + 1
            while (nchar(sequ) < s) {
                sequ <- paste0(sequ, gsub(" ", "", X[j]))
                j <- j + 1
            }
            obj[i, ] <- unlist(strsplit(sequ, NULL))
        }
        if (is.null(seq.names))
            seq.names <- getTaxaNames(taxa)
    }
#    if (format == "fasta") return(read.FASTA.AA(file))
#        start <- grep("^ {0,}>", X)
#        taxa <- X[start]
#        n <- length(taxa)
#        obj <- vector("list", n)
#        if (is.null(seq.names)) {
#            taxa <- sub("^ {0,}>", "", taxa)
#            seq.names <- getTaxaNames(taxa)
#        }
#        start <- c(start, length(X) + 1)
#        for (i in 1:n) obj[[i]] <- unlist(strsplit(gsub(" ",
#            "", X[(start[i] + 1):(start[i + 1] - 1)]), NULL))
#    }
    if (phylip) {
        rownames(obj) <- seq.names
        obj <- tolower(obj)
    }
    else {
        names(obj) <- seq.names
        obj <- lapply(obj, tolower)
    }
    obj
}


#' @rdname phyDat
#' @export
genlight2phyDat <- function(x, ambiguity=NA){
    tmp <- as.matrix(x)
    lev <- na.omit(unique(as.vector(tmp)))
    phyDat(tmp, "USER", levels=lev, ambiguity=ambiguity)
}


#' @rdname phyDat
#' @method image phyDat
#' @export
image.phyDat <- function(x, ...){
    if(attr(x, "type")=="AA")image(as.AAbin(x), ...)
    if(attr(x, "type")=="DNA")image(as.DNAbin(x), ...)
    else return(NULL)
}

