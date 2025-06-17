#' UPGMA, WPGMA and sUPGMA
#'
#' UPGMA and WPGMA clustering. UPGMA (Sokal and Michener 1958) and WPGMA
#' (McQuitty 1966) are a wrapper function around \code{\link[stats]{hclust}}
#' returning a \code{phylo} object.
## By default UPGMA additionally performs nearest neighbor interchange (NNI)
## tree rearrangements to improve the phylogeny in respect to the munimum
## evolution criterium.
#' \code{supgma} perform serial sampled UPGMA similar to Drummond and Rodrigo
#' (2000).
#'
#' UPGMA and WPGMA return ultrametric trees, it is implicitly assumed that the
#' distances supplied are close to ultrametric, e.g. hold the molecular clock
#' assumption. Neighbor Joining (NJ) \code{\link[ape]{nj}} and fastME
#' \code{\link[ape]{fastme}} relax this assumption to additive distances.
#' sUPGMA assumes tip dated data.
#'
#' \code{tip.dates} should be a vector of sampling times, in any time unit, with
#' time increasing toward the present. For example, this may be in units of
#' "days since study start" or "years since 10.000 BCE", but not "millions of
#' years ago".
#'
#' @param D A distance matrix, i.e. an object of class \code{dist}. If an matrix
#' is supplied it is tried to covert it do a \code{dist} object.
#' @param method The agglomeration method to be used. This should be (an
#' unambiguous abbreviation of) one of "ward", "single", "complete", "average",
#' "mcquitty", "median" or "centroid". The default is "average" for UPGMA and
#' "mcquitty" for WPGMA.
## @param NNI logical whether make nearest neighbor rearrangements to improve
## the tree. Currently only available for UPGMA, i.e. \code{method="average"}.
#' @param trace	 Show output during optimization (see details).
#' @param tip.dates	 A named vector of sampling times associated to the tips.
#' @param \dots Further arguments passed to or from other methods.
#' @return A phylogenetic tree of class \code{phylo}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{hclust}}, \code{\link{dist.hamming}}, \code{\link{NJ}},
#' \code{\link[ape]{as.phylo}}, \code{\link[ape]{fastme}},
#' \code{\link{nnls.tree}}, \code{\link[ape]{rtt}}
#' @references Sneath, P. H., & Sokal, R. R. (1973). \emph{Numerical taxonomy.
#' The principles and practice of numerical classification.}
#'
#' Sokal, R. R., & Michener, C. D. (1958). A statistical method for evaluating
#' systematic relationships. \emph{University of Kansas Scientific Bulletin},
#' v. 38.
#'
#' Drummond, A., & Rodrigo, A. G. (2000). Reconstructing genealogies of serial
#' samples under the assumption of a molecular clock using serial-sample UPGMA.
#' \emph{Molecular Biology and Evolution}, \bold{17(12)}, 1807-1815.
#'
#' McQuitty, L.L. (1966). Similarity Analysis by Reciprocal Pairs for Discrete
#' and Continuous Data. \emph{Educational and Psychological Measurement},
#' \bold{26}, 825–831.
#' @keywords cluster
#' @examples
#'
#' data(Laurasiatherian)
#' dm <- dist.ml(Laurasiatherian)
#' tree <- upgma(dm)
#' plot(tree)
#'
#' @rdname upgma
#' @export
upgma <- function(D, method = "average", ...) {
  if(anyNA(D)) stop("missing values are not allowed in the distance matrix")
  if(any(is.infinite(D)))
    stop("infinite values are not allowed in the distance matrix")
  if (hasArg(NNI)) NNI <- list(...)$NNI
  else NNI <- FALSE
  method <- match.arg(method, c("average", "ward.D", "ward.D2", "single",
                      "complete", "mcquitty", "median", "centroid"))
  DD <- as.dist(D)
  hc <- hclust(DD, method = method)
  result <- as.phylo(hc)
  if(NNI && method == "average") result <- upgma_nni(DD, tree = result)
  result <- reorder(result, "postorder")
  result
}


#' @rdname upgma
#' @export
wpgma <- function(D, method = "mcquitty", ...) {
  if(anyNA(D)) stop("missing values are not allowed in the distance matrix")
  if(any(is.infinite(D)))
    stop("infinite values are not allowed in the distance matrix")
  method <- match.arg(method, c("average", "ward.D", "ward.D2", "single",
                      "complete", "mcquitty", "median", "centroid"))
  DD <- as.dist(D)
  hc <- hclust(DD, method = method, ...)
  result <- as.phylo(hc)
  result <- reorder(result, "postorder")
  result
}


#' @rdname upgma
#' @export
supgma <- function(D, tip.dates, trace=0, ...){
  if(anyNA(D)) stop("missing values are not allowed in the distance matrix")
  if(any(is.infinite(D)))
    stop("infinite values are not allowed in the distance matrix")
  D <- as.dist(D)
  tree <- fastme.ols(D)
  tree <- relabel(tree, attr(D, "Labels"))
  tree <- rtt(tree, tip.dates[tree$tip.label])
  tree <- nnls.tree(D, tree, method = "tipdated",
                  tip.dates=tip.dates[tree$tip.label])
  rate_0 <- attr(tree, "rate")
  dm_td <- designTipDated(tree, tip.dates[tree$tip.label])
  dm_td <- dm_td[, ncol(dm_td)]
  me_start <- sum(tree$edge.length* rate_0)
  if(trace) cat("ME: ", me_start, "rate:", rate_0, "\n")
  me_best <- me_start
  me <- 0
  iter <- TRUE
  swapi <- 1
  while(iter){
    D_tmp <- D - rate_0 * dm_td
    tree_tmp <- upgma(D_tmp, NNI=TRUE)
    swapi <- attr(tree_tmp, "swap")
    tree_tmp <- nnls.tree(D, tree_tmp, method = "tipdated",
                       tip.dates=tip.dates[tree$tip.label])
    rate_0 <- attr(tree_tmp, "rate")
    me <- sum(tree_tmp$edge.length * rate_0)
    if(trace) cat("ME: ", me, "rate:", rate_0, "\n")
    if(me < me_best){
      tree <- tree_tmp
      me_best <- me
    } else iter <- FALSE
  }
  tree$tip.dates <- tip.dates[tree$tip.label]
  tree
}


# SRR tags ---------------------------------------------------------------------

#' @srrstats {G1.1} UPGMA(NNI=TRUE) is the first implementation of a new
#' algorithm. UPGMA(NNI=FALSE), WPGMA() are wrapper around hclust. sUPGMA is the
#' first implementation in R (I am aware of).
#' @srrstats {G1.4} This file's functions are all documented with `{roxygen2}`.
#' @srrstats {G2.15, G2.16} Checks for missing, infinite values
#' @srrstats {UL1.0} explicitly document expected format (types or classes)
#'    for input data
#' @srrstats {UL1.1} test for distance matrix.
#' @srrstats {UL1.3, UL1.3a} tip labels correspond to labels from the distances.
#' @srrstats {UL1.4, UL1.4a, UL1.4b} Cite assumption about ultrametric and tip
#' dated distances.
#'
NULL

#' @srrstats {G1.0} in the lines folloing: 41
#' @srrstats {G2.3, G2.3a} in lines: 69, 86
