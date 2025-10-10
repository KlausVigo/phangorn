#' Discrete Gamma and Beta distribution
#'
#' \code{terraces} visualizes in likelihood surface for the tree space (Sanderson ) Usually trees are from a bootstrap or MCMC sample. There the first two
#' axis are the principle components of distances between trees and the third
#' axis is the likelihood value.
## (it could also parsiomony score, minimum evolution criteria or least squares)
#'
#' @param x an object of class \code{pml}
#' @param bs_trees an object of class \code{multiPhylo}
#' @param dist_fun a function to compute distances between trees see e.g.
#' \code{\link{RF.dist}}
#' @param di2multi logical, should trees multichotomies get collapsed.
#' @param tol a numeric value giving the tolerance to consider a branch length
#' significantly greater than zero.
#' @param plot loggical if TRUE a 3D scatter is shown.
#' @param add whether to add the points to an existing plot.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{terraces} silently returns a matrix.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' #' @references Sanderson, M.J., McMahon, M.M. and Steel, M. (2011). Terraces in phylogenetic tree space. \emph{Science},
#' \bold{333}, 448--450.
#'
#' @importFrom stats cmdscale
#' @seealso \code{\link{pml_bb}, \link{optim.pml}, \link{pratchet},
#' \link{RF.dist}, \link[ape]{di2multi}, \link{cmdscale}}.
#' @examples
#' \dontrun{
#' data(woodmouse)
#' library(rgl)
#' fit <- pml_bb(woodmouse, model="JC")
#' library(rgl)
#' open3d()
#' terraces(fit)
#' }
#' @export
terraces <- function(x, bs_trees=NULL, dist_fun="RF.dist", di2multi=TRUE,
                     tol=2e-8, plot=TRUE, add=FALSE, ...) {
  assert_pml(x)
  tree <- x$tree
  if(is.null(bs_trees)) bs_trees <- x$bs
  trees <- c(tree, bs_trees)

  if(di2multi)trees <- multi2di(trees, tol)
  tmp <- hash(trees)
  trees <- trees[!duplicated(tmp)]

  if(length(trees) < 3) stop("less than 3 different trees found!")
  dm <- do.call(dist_fun, list(trees))

  xy <- cmdscale(dm)
  fun <- function(x, tree){
    tmp <- update(x, tree=tree)
    logLik(tmp)
  }
  z <- vapply(trees, fun, -Inf,  x=x)
  xyz <- cbind(xy, z)
  colnames(xyz) = c("prc_1", "prc_2", "log-likelihood")
  if (plot){
#    chk <- requireNamespace("rgl", quietly = TRUE)
#    if (!chk) {
#      warning("plot=TRUE requires the package 'rgl'!\n")
#    } else {
#      nr <- nrow(xyz)
#      col <- rep("black", nr)
#      col[1] <- "red"
#      if (!add) rgl::next3d()
#      rgl::plot3d(xyz, type="h", zlim = range(xyz[,3]), zlab="log(L)", lwd=2)
#      rgl::points3d(xyz, size=10, col=col)
#      invisible(xyz)
#    }
    plot_terraces(xyz, size=site, lwd=lwd, ...)
  } else return(xyz)
  #  scatterplot3d::scatterplot3d(xyz[,1], xyz[,2], xyz[,3], type = "h")
  #  plot3D::scatter3D(xyz[,1], xyz[,2], xyz[,3], type = "h")
# xyz
}


plot_terraces <- function(xyz, size=10, lwd=2, add=FALSE, ...){
  chk <- requireNamespace("rgl", quietly = TRUE)
  if (!chk) {
    warning("plot=TRUE requires the package 'rgl'!\n")
  } else {
    nr <- nrow(xyz)
    col <- rep("black", nr)
    col[1] <- "red"
    if (!add) rgl::next3d()
    rgl::plot3d(xyz, type="h", zlim = range(xyz[,3]), lwd=2, ...)
    rgl::points3d(xyz, size=10, col=col, ...)

  }
  invisible(xyz)
}


