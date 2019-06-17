#' Pairwise Polymorphism P-Distances from DNA Sequences
#'
#' This function computes a matrix of pairwise uncorrected polymorphism
#' p-distances. Polymorphism p-distances include intra-individual site
#' polymorphisms (2ISPs; e.g. "R") when calculating genetic distances.
#'
#' The polymorphism p-distances (Potts et al. 2014) have been developed to
#' analyse intra-individual variant polymorphism. For example, the widely used
#' ribosomal internal transcribed spacer (ITS) region (e.g. Alvarez and Wendel,
#' 2003) consists of 100's to 1000's of units within array across potentially
#' multiple nucleolus organising regions (Bailey et al., 2003; Goeker and
#' Grimm, 2008). This can give rise to intra-individual site polymorphisms
#' (2ISPs) that can be detected from direct-PCR sequencing or cloning . Clone
#' consensus sequences (see Goeker and Grimm, 2008) can be analysed with this
#' function.
#'
#' @param x a matrix containing DNA sequences; this must be of class "phyDat"
#' (use as.phyDat to convert from DNAbin objects).
#' @param cost A cost matrix or "polymorphism" for a predefined one.
#' @param ignore.indels a logical indicating whether gaps are treated as fifth
#' state or not. Warning, each gap site is treated as a characters, so an an
#' indel that spans a number of base positions would be treated as multiple
#' character states.
#' @return an object of class \code{dist}.
#' @author Klaus Schliep and Alastair Potts
#' @seealso \code{\link[ape]{dist.dna}}, \code{\link[phangorn]{dist.hamming}}
#' @references Alvarez, I., and J. F. Wendel. (2003) Ribosomal ITS sequences
#' and plant phylogenetic inference. \emph{ Molecular Phylogenetics and
#' Evolution}, \bold{29}, 417--434.
#'
#' Bailey, C. D., T. G. Carr, S. A. Harris, and C. E. Hughes. (2003)
#' Characterization of angiosperm nrDNA polymorphism, paralogy, and
#' pseudogenes. \emph{Molecular Phylogenetics and Evolution} \bold{29},
#' 435--455.
#'
#' Goeker, M., and G. Grimm. (2008) General functions to transform associate
#' data to host data, and their use in phylogenetic inference from sequences
#' with intra-individual variability. \emph{BMC Evolutionary Biology},
#' \bold{8}:86.
#'
#' Potts, A.J., T.A. Hedderson, and G.W. Grimm. (2014) Constructing phylogenies
#' in the presence of intra-individual site polymorphisms (2ISPs) with a focus
#' on the nuclear ribosomal cistron. \emph{Systematic Biology}, \bold{63},
#' 1--16
#' @keywords cluster
#' @examples
#'
#' data(Laurasiatherian)
#' laura <- as.DNAbin(Laurasiatherian)
#'
#' dm <- dist.p(Laurasiatherian, "polymorphism")
#'
#' ########################################################
#' # Dealing with indel 2ISPs
#' # These can be coded using an "x" in the alignment. Note
#' # that as.character usage in the read.dna() function.
#' #########################################################
#' cat("3 5",
#'     "No305     ATRA-",
#'     "No304     ATAYX",
#'     "No306     ATAGA",
#'     file = "exdna.txt", sep = "\n")
#' (ex.dna <- read.dna("exdna.txt", format = "sequential", as.character=TRUE))
#' dat <- phyDat(ex.dna, "USER", levels=unique(as.vector(ex.dna)))
#' dist.p(dat)
#'
#' unlink("exdna.txt")
#'
#' @export dist.p
dist.p <- function(x, cost = "polymorphism", ignore.indels = TRUE) {
  if(inherits(x, "DNAbin")) x <- as.phyDat(x)
  if (!inherits(x, "phyDat")) {
    stop("x must be of class phyDat")
  }

  l <- length(x)
  weight <- attr(x, "weight")
  n <- length(attr(x, "allLevels"))
  d <- numeric((l * (l - 1)) / 2)
  lev <- attr(x, "allLevels")
  if (is.null(cost)) {
    cost <- 1 - diag(n)
    dimnames(cost) <- list(lev, lev)
  }
  #    if(cost=="polymorphism" && attr(x, "type")=="DNA"){
  if (cost == "polymorphism") {
    costLev <- c("a", "c", "t", "u", "g", "x", "m", "r", "w", "s", "y", "k",
                 "v", "h", "d", "b", "-", "?", "n")

    cost <- matrix(c(
      # a,c,t,u,g,X,m,r,w,s,y,k,v,h,d,b,-,?,n,
      0, 2, 2, 2, 2, 1, 1, 1, 1, 3, 3, 3, 2, 2, 2, 4, 2, 0, 0, # a
      2, 0, 2, 2, 2, 1, 1, 3, 3, 1, 1, 3, 2, 2, 4, 2, 2, 0, 0, # c
      2, 2, 0, 0, 2, 1, 3, 3, 1, 3, 1, 1, 4, 2, 2, 2, 2, 0, 0, # t
      2, 2, 0, 0, 2, 1, 3, 3, 1, 3, 1, 1, 4, 2, 2, 2, 2, 0, 0, # u
      2, 2, 2, 2, 0, 1, 3, 1, 3, 1, 3, 1, 2, 4, 2, 2, 2, 0, 0, # g
      1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, # X
      1, 1, 3, 3, 3, 1, 0, 2, 2, 2, 2, 4, 1, 1, 3, 3, 3, 0, 0, # m
      1, 3, 3, 3, 1, 1, 2, 0, 2, 2, 4, 2, 1, 3, 1, 3, 3, 0, 0, # r
      1, 3, 1, 1, 3, 1, 2, 2, 0, 4, 2, 2, 3, 1, 1, 3, 3, 0, 0, # w
      3, 1, 3, 3, 1, 1, 2, 2, 4, 0, 2, 2, 1, 3, 3, 1, 3, 0, 0, # s
      3, 1, 1, 1, 3, 1, 2, 4, 2, 2, 0, 2, 3, 1, 3, 1, 3, 0, 0, # y
      3, 3, 1, 1, 1, 1, 4, 2, 2, 2, 2, 0, 3, 3, 1, 1, 3, 0, 0, # k
      2, 2, 4, 4, 2, 1, 1, 1, 3, 1, 3, 3, 0, 2, 2, 2, 4, 0, 0, # v
      2, 2, 2, 2, 4, 1, 1, 3, 1, 3, 1, 3, 2, 0, 2, 2, 4, 0, 0, # h
      2, 4, 2, 2, 2, 1, 3, 1, 1, 3, 3, 1, 2, 2, 0, 2, 4, 0, 0, # d
      4, 2, 2, 2, 2, 1, 3, 3, 3, 1, 1, 1, 2, 2, 2, 0, 4, 0, 0, # b
      2, 2, 2, 2, 2, 1, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 0, 0, 0, #-
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # ?
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ), # n
    ncol = 19, nrow = 19, dimnames = list(costLev, costLev)
    )
  }

  lev1 <- dimnames(cost)[[1]]

  if (any(is.na(match(lev, lev1))))
    stop("Levels of x are not in levels of cost matrix!")

  if (ignore.indels) {
    cost["-", ] <- 0
    cost[, "-"] <- 0
  }

  cost <- cost[lev, lev]
  k <- 1
  for (i in 1:(l - 1)) {
    for (j in (i + 1):l) {
      d[k] <- sum(weight * cost[cbind(x[[i]], x[[j]])])
      k <- k + 1
    }
  }
  attr(d, "Size") <- l
  if (is.list(x)) {
    attr(d, "Labels") <- names(x)
  } else {
    attr(d, "Labels") <- colnames(x)
  }
  attr(d, "Diag") <- FALSE
  attr(d, "Upper") <- FALSE
  attr(d, "call") <- match.call()
  attr(d, "method") <- "p"
  class(d) <- "dist"
  return(d)
}
