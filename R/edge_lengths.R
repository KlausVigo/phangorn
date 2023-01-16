edge_length_matrix <- function(tree, trees, rooted=TRUE){
  if(!inherits(trees, "multiPhylo")) stop("Trees must be of class multiPhylo!")
  trees <- .uncompressTipLabel(trees) |> .compressTipLabel(ref=tree$tip.label)
  if(!rooted){
    trees <- di2multi(trees) |> unroot()
    tree <- di2multi(tree) |> unroot()
  }
  else{
    if(!is.rooted(tree) || any(!is.rooted(trees))) stop("All trees need to be rooted!")
  }
  fun <- function(x){
    el <- numeric(max(x$edge))
    el[x$edge[,2]] <- x$edge.length
  }
  bp <- bip(tree)
  if(!rooted) bp <- SHORTwise(bp)
  m <- length(bp)
  res <- matrix(NA_real_, length(trees), m)
  for(i in seq_along(trees)){
    tmp <- bip(trees[[i]])
    if(!rooted){
      l <- fun(trees[[i]])
      tmp <- SHORTwise(tmp)
    }
    else l <- nodeHeight(trees[[i]])
    ind <- match(tmp, bp)
    nna <- !is.na(ind)
    res[i, ind[nna]] <- l[nna]
  }
  res
}


##' @title Assign and compute edge lengths from a sample of trees
##' @description This command can infer some average edge lengths and assign
##' them from a (bootstrap/MCMC) sample.
##' @param tree tree where edge lengths are assigned to.
##' @param trees an object of class multiPhylo, where the average for the edges
##' is computed from.
##' @param fun a function to compute the average (default is median).
##' @param rooted the height of the boxes.
##' @param \dots arguments passed to [graphics::legend()]
##' @return NULL
##' @author Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep
##' @importFrom graphics legend rect
##' @examples
##' \dontrun{
##' data("Laurasiatherian")
##' dm <- dist.hamming(Laurasiatherian)
##' set.seed(123)
##' trees <- bootstrap.phyDat(Laurasiatherian,
##'                           FUN=function(x)upgma(dist.hamming(x)), bs=100)
##'                           tree <- plotBS(tree, trees, "phylogram")
##' tree_ultra <- allCompat(trees, rooted=TRUE) |> add_edge_length(trees)
##' tree_unrooted <- allCompat(trees, rooted=FALSE) |>
##'                  add_edge_length(trees, rooted=FALSE)
##' plot(tree_ultra)
##' plot(tree_unrooted)
##' get
##' }
##' @seealso [ape::node.depth.edgelength][ape::consensus]
##' @keywords aplot
add_edge_length <- function(tree, trees, fun=\(x)median(na.omit(x)),
                            rooted=TRUE){
  X <- edge_length_matrix(tree, trees, rooted)
  nh <- apply(X, 2, fun)
  if(rooted) tree$edge.length <- nh[tree$edge[,1]] - nh[tree$edge[,2]]
  else tree$edge.length <- nh[tree$edge[,2]]
  tree
}




##' @title Draw Confidences Intervals on Phylogenies
##' @description This is a low-level plotting command to draw the confidence
##' intervals on the node of a tree as rectangles with coloured backgrounds.
##' @param CI output from [chronosCI()] or a similar matrix.
##' @param col95 colour used for the 95% intervals; by default: transparent
##' red.
##' @param col50 colour used for the 50% intervals; by default: transparent
##' blue.
##' @param height the height of the boxes.
##' @param legend a logical value.
##' @param \dots arguments passed to [graphics::legend()]
##' @details The matrix \code{CI} must have four rows and as many columns as the
##' number of nodes of the tree. The first and fourth rows give the lower and
##' upper bounds of the 95% confidence intervals. The second and third rows
##' give the lower and upper bounds of the 50% confidence intervals. The Trees
##' should to be rooted, either ultrametric or tip dated.
##' @return NULL
##' @author Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep
##' @importFrom graphics legend rect
##' @importFrom stats median
##' @examples
##' \dontrun{
##' data("Laurasiatherian")
##' dm <- dist.hamming(Laurasiatherian)
##' tree <- upgma(dm)
##' set.seed(123)
##' trees <- bootstrap.phyDat(Laurasiatherian,
##'                           FUN=function(x)upgma(dist.hamming(x)), bs=100)
##'                           tree <- plotBS(tree, trees, "phylogram")
##' Y <- get_CI(tree, trees)
##' tree <- plotBS(tree, trees, "phylogram")
##' draw_ci(Y)
##' }
##' @seealso [plot.phylo, plotBS]
##' @keywords aplot
##' @export
draw_ci <- function(CI, col95 = "#FF00004D", col50 = "#0000FF4D",
                          height = 0.7, legend = TRUE, ...)
{
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  direction <- lastPP$direction
  left_right <- FALSE
  if(direction=="rightwards" || direction=="leftwards"){
    left_right <- TRUE
    if(direction=="rightwards") CI <- max(lastPP$xx) - CI
    if(direction=="leftwards") CI <- min(lastPP$xx) + CI
    Y <- lastPP$yy - height / 2
  }
  else {
    if(direction=="downwards") CI <- min(lastPP$yy) + CI
    if(direction=="upwards") CI <- max(lastPP$yy) - CI
    Y <- lastPP$xx - height / 2
  }
  L <- CI[1, ]
  R <- CI[4, ]
  B <- Y[ - seq_len(lastPP$Ntip)]
  T <- B + height
  if(left_right) graphics::rect(L, B, R, T, col = col95, border = NULL)
  else graphics::rect(B, L, T, R, col = col95, border = NULL)
  L <- CI[2, ]
  R <- CI[3, ]
  if(left_right) graphics::rect(L, B, R, T, col = col50, border = NULL)
  else graphics::rect(B, L, T, R, col = col50, border = NULL)
  if (!identical(legend, FALSE)) {
    loc <- if (is.logical(legend)) "topleft" else legend
    graphics::legend(loc, legend = c("95% CI", "50% CI"), pch = 22,
                     pt.bg = c(col95, col50), col = c(col95, col50), ...)
  }
}


#' @rdname draw_ci
#' @export
get_CI <- function(tree, trees, fun=\(x)median(na.omit(x)), fun=NULL){
  if(!is.rooted(tree) || !is.rooted(trees)) stop("Trees need to be rooted!")
  X <- edge_length_matrix(tree, trees, rooted)
  fun <- function(x) quantile(na.omit(x), probs=c(.025,.25,.75,.975))
  apply(X, 2, fun)
}
