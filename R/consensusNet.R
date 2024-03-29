#' Computes a consensusNetwork from a list of trees Computes a \code{networx}
#' object from a collection of splits.
#'
#' Computes a consensusNetwork, i.e. an object of class \code{networx} from a
#' list of trees, i.e. an class of class \code{multiPhylo}. Computes a
#' \code{networx} object from a collection of splits.
#'
#'
#' @param obj An object of class multiPhylo.
#' @param prob the proportion a split has to be present in all trees to be
#' represented in the network.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{consensusNet} returns an object of class networx.  This is
#' just an intermediate to plot phylogenetic networks with igraph.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{splitsNetwork}}, \code{\link{neighborNet}},
#' \code{\link{lento}}, \code{\link{distanceHadamard}},
#' \code{\link{plot.networx}}, \code{\link{maxCladeCred}}
#' @references Holland B.R., Huber K.T., Moulton V., Lockhart P.J. (2004) Using
#' consensus networks to visualize contradictory evidence for species
#' phylogeny. \emph{Molecular Biology and Evolution}, \bold{21}, 1459--61
#' @keywords hplot
#' @examples
#'
#' data(Laurasiatherian)
#' set.seed(1)
#' bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)nj(dist.hamming(x)),
#'     bs=50)
#' cnet <- consensusNet(bs, .3)
#' plot(cnet, angle=-60, direction="axial")
#' \dontrun{
#' library(rgl)
#' open3d()
#' plot(cnet, type = "3D", show.tip.label=FALSE, show.nodes=TRUE)
#' plot(cnet, type = "equal angle", show.edge.label=TRUE)
#'
#' tmpfile <- normalizePath(system.file(
#'               "extdata/trees/RAxML_bootstrap.woodmouse", package="phangorn"))
#' trees <- read.tree(tmpfile)
#' cnet_woodmouse <- consensusNet(trees, .3)
#' plot(cnet_woodmouse, type = "equal angle", show.edge.label=TRUE)
#' }
#'
#' @export consensusNet
consensusNet <- function(obj, prob = 0.3, ...) {
  l <- length(obj)
  spl <- as.splits(obj)
  w <- attr(spl, "weights")
  ind <- (w / l) > prob
  spl <- spl[ind]
  attr(spl, "confidences") <- (w / l)[ind]
  #    attr(spl, "weights") = w[ind]
  res <- as.networx(spl)
  res$edge.labels <- res$edge.length / l * 100
  #as.character(res$edge.length / l * 100)
  res$edge.labels[res$edge[, 2] <= length(res$tip.label)] <- NA_real_ #""
  reorder(res)
}
