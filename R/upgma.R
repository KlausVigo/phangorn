#' UPGMA, WPGMA and sUPGMA
#'
#' UPGMA and WPGMA clustering. UPGMA and WPGMA are a wrapper function around
#' \code{\link[stats]{hclust}} returning a \code{phylo} object.
## UPGMA additionally performs nearest neighbor interchange (NNI) tree
## rearrangements to improve the phylogeny (Schliep et al. 2024).
#' \code{supgma} perform serial sampled UPGMA similar to Drummond and Rodrigo
#' (2000).
##  and also performing NNI rearrangements.
#'
#' @param D A distance matrix.
#' @param method The agglomeration method to be used. This should be (an
#' unambiguous abbreviation of) one of "ward", "single", "complete", "average",
#' "mcquitty", "median" or "centroid". The default is "average".
## @param NNI logical whether make nearest neighbor rearrangements to improve
## the tree. Currently only available for \code{method="average"}.
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
  method <- match.arg(method, c("average", "ward.D", "ward.D2", "single",
                      "complete", "average", "mcquitty", "median", "centroid"))
  DD <- as.dist(D)
  hc <- hclust(DD, method = method)
  result <- as.phylo(hc)
#  if(NNI && method=="average"){
#    result <- upgma_nni(DD, tree=result, ...)
#  }
  result <- reorder(result, "postorder")
  result
}


#' @rdname upgma
#' @export
wpgma <- function(D, method = "mcquitty", ...) {
  method <- match.arg(method, c("average", "ward.D", "ward.D2", "single",
                      "complete", "average", "mcquitty", "median", "centroid"))
  DD <- as.dist(D)
  hc <- hclust(DD, method = method, ...)
  result <- as.phylo(hc)
  result <- reorder(result, "postorder")
  result
}


#' @rdname upgma
#' @export
supgma <- function(D, tip.dates, trace=0){
  tree <- fastme.ols(D)
  tree <- checkLabels(tree, attr(D, "Labels"))
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
    tree_tmp <- upgma(D_tmp)
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
