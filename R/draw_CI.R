edge_length_matrix <- function(tree, trees, rooted=TRUE){
  if(!inherits(trees, "multiPhylo")) stop("Trees must be of class multiPhylo!")
  trees <- .uncompressTipLabel(trees) |> .compressTipLabel(ref=tree$tip.label)
  if(!rooted){
    trees <- unroot(trees)
    tree <- unroot(tree)
  }
  else{
    if(!is.rooted(tree) || any(!is.rooted(trees))) stop("All trees need to be rooted!")
  }
  fun <- function(x){
    el <- numeric(max(x$edge))
    el[x$edge[,2]] <- x$edge.length
    el
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
    ind <- fmatch(tmp, bp)
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
##' @param rooted rooted logical, if FALSE edge lengths is a function of the
##' observed splits, if TRUE edge lengths are estimated from height for the
##' observed clades.
##' @return NULL
##' @author Klaus Schliep
##' @importFrom graphics legend rect
##' @examples
##' data("Laurasiatherian")
##' set.seed(123)
##' bs <- bootstrap.phyDat(Laurasiatherian,
##'                 FUN=function(x)upgma(dist.ml(x)), bs=100)
##' tree_compat <- allCompat(bs, rooted=TRUE) |>
##'               add_edge_length(bs)
##' plot(tree_compat)
##' add_boxplot(tree_compat, bs)
##' @seealso \code{\link{node.depth.edgelength}}, \code{\link{consensus}},
##' \code{\link{maxCladeCred}}
##' @keywords aplot
##' @export
add_edge_length <- function(tree, trees, fun=\(x)median(na.omit(x)),
                            rooted=TRUE){
  if(!rooted) tree <- unroot(tree)
  X <- edge_length_matrix(tree, trees, rooted)
  nh <- apply(X, 2, fun)
  if(rooted) tree$edge.length <- nh[tree$edge[,1]] - nh[tree$edge[,2]]
  else tree$edge.length <- nh[tree$edge[,2]]
  tree
}


##' @title Draw Confidences Intervals on Phylogenies
##' @description These are low-level plotting commands to draw the confidence
##' intervals on the node of a tree as rectangles with coloured backgrounds or
##' add boxplots to ultrametric or tipdated trees.
##' @param tree a phylogenetic tree to which the confidences should be added.
##' @param trees phylogenetic trees, i.e. an object of class `multiPhylo`
##' @param col95 colour used for the 95% intervals; by default: transparent
##' red.
##' @param col50 colour used for the 50% intervals; by default: transparent
##' blue.
##' @param height the height of the boxes.
##' @param legend a logical value.
##' @param \dots arguments passed to other functions, \code{\link{legend}} or
##' \code{\link{bxp}}.
##' @details All trees should to be rooted, either ultrametric or tip dated.
##' @return NULL
##' @author Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep
##' @importFrom graphics legend rect bxp boxplot
##' @importFrom stats median
##' @examples
##' data("Laurasiatherian")
##' dm <- dist.hamming(Laurasiatherian)
##' tree <- upgma(dm)
##' set.seed(123)
##' trees <- bootstrap.phyDat(Laurasiatherian,
##'                           FUN=function(x)upgma(dist.hamming(x)), bs=100)
##'                           tree <- plotBS(tree, trees, "phylogram")
##' tree <- plotBS(tree, trees, "phylogram")
##' add_ci(tree, trees)
##' plot(tree, direction="downwards")
##' add_boxplot(tree, trees, boxwex=.7)
##' @seealso \code{\link{plot.phylo}}, \code{\link{plotBS}}
##' @keywords aplot
##' @rdname add_ci
##' @export
add_ci <- function(tree, trees, col95 = "#FF00004D", col50 = "#0000FF4D",
                    height = 0.7, legend = TRUE, ...)
{
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  direction <- lastPP$direction
  if(!is.rooted(tree) || !all(is.rooted(trees))) stop("Trees need to be rooted!")
  X <- edge_length_matrix(tree, trees, rooted=TRUE)[, -(seq_along(Ntip(tree)))]
  CI <- apply(X, 2, FUN=\(x)quantile(na.omit(x), probs=c(.025,.25,.75,.975)))
  horizontal <- FALSE
  if(direction=="rightwards" || direction=="leftwards"){
    horizontal <- TRUE
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
  if(horizontal) graphics::rect(L, B, R, T, col = col95, border = NULL)
  else graphics::rect(B, L, T, R, col = col95, border = NULL)
  L <- CI[2, ]
  R <- CI[3, ]
  if(horizontal) graphics::rect(L, B, R, T, col = col50, border = NULL)
  else graphics::rect(B, L, T, R, col = col50, border = NULL)
  if (!identical(legend, FALSE)) {
    loc <- if (is.logical(legend)) "topleft" else legend
    graphics::legend(loc, legend = c("95% CI", "50% CI"), pch = 22,
                     pt.bg = c(col95, col50), col = c(col95, col50), ...)
  }
}

##' @rdname add_ci
##' @export
add_boxplot <- function(tree, trees, ...)
{
  X <- edge_length_matrix(tree, trees, rooted=TRUE)
  X <- X[, -c(1:Ntip(tree))]
  tmp <- boxplot(X, plot=FALSE)
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  CI <- tmp$stats
  out <- tmp$out
  direction <- lastPP$direction
  horizontal <- FALSE
  if(direction=="rightwards" || direction=="leftwards"){
    horizontal <- TRUE
    if(direction=="rightwards"){
      CI <- max(lastPP$xx) - CI
      out <- max(lastPP$xx) - out
    }
    if(direction=="leftwards"){
      CI <- min(lastPP$xx) + CI
      out <- min(lastPP$xx) + out
    }
    Y <- lastPP$yy # - height / 2
  }
  else {
    if(direction=="downwards") {
      CI <- min(lastPP$yy) + CI
      out <- min(lastPP$yy) + out
    }
    if(direction=="upwards"){
      CI <- max(lastPP$yy) - CI
      out <- max(lastPP$yy) - out
    }
    Y <- lastPP$xx #- height / 2
  }
  tmp$stats <- CI
  tmp$out <- out
  Y <- Y[-c(1:Ntip(tree))]
  bxp(tmp, at=Y, horizontal=horizontal, add=TRUE, axes=FALSE, ...)
}
