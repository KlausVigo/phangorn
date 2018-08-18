################################################################################
# delta.score
################################################################################
# Calculated from mathematical description given in Gray et al. (2010) Phil.
# Trans. Roy. Soc. B.
# delta.score reference: Holland et al. (2002) Mol. Biol. Evol.
################################################################################


# Calculating Delta and Q-residual scores
# internal
delta.quartet <-
  function(quartet, dist.dna) {
    m1 <- dist.dna[quartet[1], quartet[2]] + dist.dna[quartet[3], quartet[4]]
    m2 <- dist.dna[quartet[1], quartet[3]] + dist.dna[quartet[2], quartet[4]]
    m3 <- dist.dna[quartet[1], quartet[4]] + dist.dna[quartet[2], quartet[3]]
    m <- sort(c(m1, m2, m3), decreasing = TRUE)
    if ((m[1] - m[3]) != 0) {
      ret <- (m[1] - m[2]) / (m[1] - m[3])
    } else {
      ret <- 0
    }
    return(ret)
  }




#' Computes the \eqn{\delta} score
#'
#' Computes the treelikeness
#'
#'
#' @param x an object of class \code{phyDat}
#' @param arg Specifies the return value, one of "all", "mean" or "sd"
#' @param ...  further arguments passed through \code{dist.hamming}
#' @return A vector containing the \eqn{\delta} scores.
#' @author Alastair Potts and Klaus Schliep
#' @seealso \code{\link{dist.hamming}}
#' @references BR Holland, KT Huber, A Dress, V Moulton (2002) \eqn{\delta}
#' Plots: a tool for analyzing phylogenetic distance data Russell D. Gray,
#' David Bryant, Simon J. Greenhill (2010) On the shape and fabric of human
#' history \emph{Molecular Biology and Evolution}, \bold{19(12)} 2051--2059
#'
#' Russell D. Gray, David Bryant, Simon J. Greenhill (2010) On the shape and
#' fabric of human history \emph{Phil. Trans. R. Soc. B}, \bold{365}
#' 3923--3933; DOI: 10.1098/rstb.2010.0162
#' @keywords cluster
#' @examples
#'
#' data(yeast)
#' hist(delta.score(yeast, "all"))
#'
#' @export delta.score
delta.score <- function(x, arg = "mean", ...) {
  dist.dna <- as.matrix(dist.hamming(x, ...))
  # Create all quartets
  all.quartets <- t(combn(names(x), 4))
  delta.values <- apply(all.quartets[, ], 1, delta.quartet, dist.dna)
  if (!arg %in% c("all", "mean", "sd"))
    stop("return options are: all, mean, or sd")
  if (arg == "all") return(delta.values)
  if (arg == "mean") return(mean(delta.values))
  if (arg == "sd") return(sd(delta.values))
}
