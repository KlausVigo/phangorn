#' Explore likelihood parsimony surface
#'
#' \code{terraces} visualizes in likelihood surface for the tree space
#' (Sanderson et al. 2011). Usually trees are from a bootstrap or MCMC sample.
#' There the first two axis are the principle components of distances between
#' trees and the third axis is the likelihood value or parsimony score.
## (it could also parsiomony score, minimum evolution criteria or least squares)
#'
#' @param x an object of class \code{pml}
#' @param trees an object of class \code{multiPhylo}
#' @param dist_fun a function to compute distances between trees see e.g.
#' \code{\link{RF.dist}}
#' @param di2multi logical, should trees multichotomies get collapsed. Useful
#' for Robinson-Foulds distance. If edge length are used to compute the
#' distance, e.g. Kuhner-Felsenstein distance, this is not needed.
#' @param tol a numeric value giving the tolerance to consider a branch length
#' significantly greater than zero.
#' @param plot loggical if TRUE a 3D scatter is shown.
#' @param add whether to add the points to an existing plot.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{terraces} silently returns a matrix.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @references Sanderson, M.J., McMahon, M.M. and Steel, M. (2011). Terraces in phylogenetic tree space. \emph{Science},
#' \bold{333}, 448--450.
#'
#' @importFrom stats cmdscale
#' @seealso \code{\link{pml_bb}, \link{optim.pml}, \link{pratchet},
#' \link{RF.dist}, \link[ape]{di2multi}, \link{cmdscale}}.
#' @examples
#' \dontrun{
#' data(woodmouse)
#' fit <- pml_bb(woodmouse, model="JC")
#' library(rgl)
#' open3d()
#' terraces(fit)
#' }
#' @export
terraces <- function(x, trees=NULL, dist_fun="SPR.dist", di2multi=FALSE,
                     tol=2e-8, plot=TRUE, add=FALSE, ...) {
  assert_pml(x)
  tree <- x$tree
  if(is.null(trees)) trees <- x$bs
  else assert_multiPhylo(trees)
  trees <- c(tree, trees)
  trees <- .compressTipLabel(trees)

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
    plot_terraces(xyz, add=add, ...)
  } else return(xyz)
}


terraces_pars <- function(x, trees, dist_fun="RF.dist", di2multi=TRUE,
                          tol=2e-8){
  assert_phyDat(x)
  assert_multiPhylo(trees)
  clean_phylo(trees, compress = TRUE)
  if (di2multi) trees <- multi2di(trees, tol)
  if(length(trees) < 3) stop("less than 3 different trees found!")
  dm <- do.call(dist_fun, list(trees))
  xy <- cmdscale(dm)
  z <- parsimony(trees, x)
  xyz <- cbind(xy, z)
  colnames(xyz) = c("prc_1", "prc_2", "parsimony score")
}

plot_terraces <- function(xyz, size=10, lwd=2, pkg="rgl",
                          max=TRUE, add=FALSE, ...){
  match.arg <- match.arg(pkg, c("rgl", "plot3D"))
  nr <- nrow(xyz)
  col <- rep("black", nr)
  if(max) {
    ind <- which(xyz[,3] > (max(xyz[,3] - 1e-8)))
    col[ind] <- "red"
  } else{
    ind <- which(xyz[,3] < (min(xyz[,3] - 1e+8)))
    col[ind] <- "red"
  }
  if(pkg=="rgl"){
    chk <- requireNamespace("rgl", quietly = TRUE)
    if (!chk) {
      warning("package 'rgl' is required!\n")
    } else {
      if (!add) rgl::next3d()
      rgl::plot3d(xyz, type="h", zlim = range(xyz[,3]), lwd=2, ...)
      rgl::points3d(xyz, size=10, col=col, ...)
    }
  } else if(pkg=="plot3D"){
    chk <- requireNamespace("plot3D", quietly = TRUE)
    if (!chk) {
      warning("package 'plot3D' is required!\n")
    } else {
      attachNamespace("plot3D")
      col_names <- colnames(xyz)
      plot3D::scatter3D(xyz[,1], xyz[,2], xyz[,3], type = "h",
                        xlab = col_names[1], ylab = col_names[2],
                        zlab = col_names[3], ...)
    }
  }
  invisible(xyz)
}

#  scatterplot3d::scatterplot3d(xyz[,1], xyz[,2], xyz[,3], type = "h")

