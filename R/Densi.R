getAges <- function(x) {
  fun <- function(x) max(node.depth.edgelength(x))
  height <- NULL
  if (inherits(x, "phylo")) height <- fun(x)
  if (inherits(x, "multiPhylo")) {
    if (!is.null(attr(x, "TipLabel"))) x <- .uncompressTipLabel(x)
    x <- unclass(x)
    height <- vapply(x, fun, 0)
  }
  height
}


add_tiplabels <- function(xy, tip.label, direction, adj, font, srt = 0, cex = 1,
                          col = 1, label.offset = 0, maxBT=1) {
  direction <- match.arg(direction, c("rightwards", "leftwards",  "upwards",
    "downwards"))
  horizontal <- direction %in% c("rightwards", "leftwards")
  nTips <- length(tip.label)
  xx <- rep(maxBT, nrow(xy))
  yy <- xy[, 2 ]
  if (direction == "leftwards" | direction == "downwards") xx <- xx * 0
  if (!horizontal) {
#    tmp <- yy
    yy <- xx
    xx <- xy[, 1]
  }
  MAXSTRING <- max(strwidth(tip.label, cex = cex))
  loy <- 0
  if (direction == "rightwards") lox <- label.offset + MAXSTRING * 1.05 * adj
  if (direction == "leftwards")
    lox <- -label.offset - MAXSTRING * 1.05 * (1 - adj)
  if (!horizontal) {
    psr <- par("usr")
    MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3]) / (psr[2] - psr[1])
    loy <- label.offset + MAXSTRING * 1.05 * adj
    lox <- 0
    srt <- 90 + srt
    if (direction == "downwards") {
      loy <- -loy
      srt <- 180 + srt
    }
  }
  text(xx[1:nTips] + lox, yy[1:nTips] + loy, tip.label, adj = adj,
    font = font, srt = srt, cex = cex, col = col)
}



#' Plots a densiTree.
#'
#' An R function to plot trees similar to those produced by DensiTree.
#'
#' If no consensus tree is provided \code{densiTree} computes a consensus tree,
#' and if the input trees have different labels a mrp.supertree as a backbone.
#' This should avoid too many unnecessary crossings of edges.
#' Trees should be rooted, other wise the output may not be visually pleasing.
#' \code{jitter} shifts trees a bit so that they are not exactly on top of each
#' other.
#' If \code{amount == 0}, it is ignored. If \code{random=TRUE} the result of the
#' permutation is \code{runif(n, -amount, amount)}, otherwise
#' \code{seq(-amount, amount, length=n)}, where \code{n <- length(x)}.
#' @param x an object of class \code{multiPhylo}.
#' @param type a character string specifying the type of phylogeny, so far
#' "cladogram" (default) or "phylogram" are supported.
#' @param alpha parameter for semi-transparent colors.
#' @param consensus A tree or character vector which is used to define the order
#' of the tip labels.
#' @param direction a character string specifying the direction of the tree.
#' Four values are possible: "rightwards" (the default), "leftwards", "upwards",
#' and "downwards".
#' @param optim not yet used.
#' @param scaleX scale trees to have identical heights.
#' @param col a scalar or vector giving the colours used to draw the edges for
#' each plotted phylogeny. These are taken to be in the same order than input
#' trees x. If fewer colours are given than the number of trees, then the
#' colours are recycled.
#' @param width edge width.
#' @param lty line type.
#' @param cex a numeric value giving the factor scaling of the tip labels.
#' @param font an integer specifying the type of font for the labels:
#' 1 (plain text), 2 (bold), 3 (italic, the default), or 4 (bold italic).
#' @param tip.color color of the tip labels.
#' @param adj a numeric specifying the justification of the text strings of the
#' labels: 0 (left-justification), 0.5 (centering), or 1 (right-justification).
#' @param srt a numeric giving how much the labels are rotated in degrees.
#' @param underscore a logical specifying whether the underscores in tip labels
#' should be written as spaces (the default) or left as are (if TRUE).
#' @param label.offset a numeric giving the space between the nodes and the tips
#' of the phylogeny and their corresponding labels.
#' @param scale.bar a logical specifying whether add scale.bar to the plot.
#' @param jitter allows to shift trees. a list with two arguments: the amount of
#' jitter and random or equally spaced (see details below)
#' @param tip.dates A named vector of sampling times associated with the tips.
#' @param xlim the x limits of the plot.
#' @param ylim the y limits of the plot.
#' @param \dots further arguments to be passed to plot.
#' @returns \code{densiTree} returns silently x.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{plot.phylo}}, \code{\link{plot.networx}},
#' \code{\link{jitter}}, \code{\link[ape]{rtt}}
#' @references densiTree is inspired from the great
#' \href{https://www.cs.auckland.ac.nz/~remco/DensiTree/}{DensiTree} program of
#' Remco Bouckaert.
#'
#' Remco R. Bouckaert (2010) DensiTree: making sense of sets of phylogenetic
#' trees \emph{Bioinformatics}, \bold{26 (10)}, 1372-1373.
#' @keywords plot
#' @examples
#'
#' data(Laurasiatherian)
#' set.seed(1)
#' bs <- bootstrap.phyDat(Laurasiatherian, FUN =
#'    function(x) upgma(dist.hamming(x)), bs=25)
#' # cladogram nice to show topological differences
#' densiTree(bs, type="cladogram", col="blue")
#' densiTree(bs, type="phylogram", col="green", direction="downwards", width=2)
#' # plot five trees slightly shifted, no transparent color
#' densiTree(bs[1:5], type="phylogram", col=1:5, width=2, jitter=
#'     list(amount=.3, random=FALSE), alpha=1)
#' \dontrun{
#' # phylograms are nice to show different age estimates
#' require(PhyloOrchard)
#' data(BinindaEmondsEtAl2007)
#' BinindaEmondsEtAl2007 <- .compressTipLabel(BinindaEmondsEtAl2007)
#' densiTree(BinindaEmondsEtAl2007, type="phylogram", col="red")
#' }
#'
#' @export
densiTree <- function(x, type = "phylogram", ..., alpha = 1 / length(x),
                      consensus = NULL, direction = "rightwards", optim = FALSE,
                      scaleX = FALSE, col = 1, width = 1, lty = 1, cex = .8,
                      font = 3, tip.color = 1, adj = 0, srt = 0,
                      underscore = FALSE, label.offset = 0, scale.bar = TRUE,
                      jitter = list(amount = 0, random = TRUE), tip.dates=NULL,
                      xlim=NULL, ylim=NULL) {
  if (!inherits(x, "multiPhylo")) stop("x must be of class multiPhylo")
  if (is.character(consensus)) {
    consensus <- stree(length(consensus), tip.label = consensus)
    consensus$edge.length <- rep(1.0, nrow(consensus$edge))
  }
  if (is.null(consensus)) {
    consensus <- tryCatch(consensus(x, p = .5),
                          error = function(e) unroot(midpoint(superTree(x))))
  }
  if (inherits(consensus, "multiPhylo")) consensus <- consensus[[1]]

  sort_tips <- function(x) {
    x <- reorder(x)
    nTip <- as.integer(length(x$tip.label))
    e2 <- x$edge[, 2]
    x$tip.label <- x$tip.label[e2[e2 <= nTip]]
    x$edge[e2 <= nTip, 2] <- as.integer(1L:nTip)
    x
  }

  type <- match.arg(type, c("phylogram", "cladogram"))
  direction <- match.arg(direction, c("rightwards", "leftwards",  "upwards",
    "downwards"))
  horizontal <- direction %in% c("rightwards", "leftwards")

  consensus <- reorder(consensus)
  nTip <- as.integer(length(consensus$tip.label))
  consensus <- sort_tips(consensus)
  consensus <- reorder(consensus, "postorder")
  at <- NULL
  maxBT <- max(getAges(x))
  if(!is.null(tip.dates)){
    root_time <- max(tip.dates) - maxBT
    label <- pretty(c(root_time, max(tip.dates)), min.n = 3)
    label <- label[label < max(tip.dates)]
    maxBT <- max(maxBT, max(tip.dates) - min(label))
    at <- maxBT - (max(tip.dates) - label) #/ maxBT
    if(direction=="leftwards" || direction=="downwards") {
      at <- at + maxBT - max(at)
    }
    scaleX <- FALSE
  }
  else {
    if (scaleX) maxBT <- 1.0
    label <- rev(pretty(c(maxBT, 0)))
    maxBT <- max(label, maxBT)
    at <- seq(0, maxBT, length.out = length(label))
  }
  xy <- plotPhyloCoor(consensus, direction = direction, ...)
  yy <- xy[, 2]
  tl <- which.max(nchar(consensus$tip.label))
  if(horizontal) pin1 <- par("pin")[1]
  else pin1 <- par("pin")[2]
  sw <- strwidth(consensus$tip.label[tl], "inch",cex = cex) / pin1 * 1.1 * maxBT
  if(is.null(xlim)){
    xlim <- switch(direction,
                   rightwards = c(0, maxBT + sw),
                   leftwards = c(0 - sw, maxBT),
                   downwards = c(0, nTip + 1),
                   upwards = c(0, nTip + 1))
  }
  if(is.null(ylim)){
    ylim <- switch(direction,
                   rightwards = c(0, nTip + 1),
                   leftwards = c(0, nTip + 1),
                   downwards = c(0 - sw, maxBT),
                   upwards = c(0, maxBT + sw))
  }
  if (direction == "rightwards") {
    plot.default(0, type = "n", xlim = xlim, ylim = ylim,
                 xlab = "", ylab = "", axes = FALSE, ...)
    if (scale.bar) axis(side = 1, at = at, labels = label, cex.axis=cex)
  }
  if (direction == "leftwards") {
    plot.default(0, type = "n", xlim = xlim, ylim = ylim,
                 xlab = "", ylab = "", axes = FALSE, ...)
    if (scale.bar) axis(side = 1, at = at, labels = rev(label), cex.axis=cex)
  }
  if (direction == "downwards") {
    plot.default(0, type = "n", xlim = xlim, ylim = ylim,
                 xlab = "", ylab = "", axes = FALSE, ...)
    if (scale.bar) axis(side = 2, at = at, labels = rev(label), cex.axis=cex)
  }
  if (direction == "upwards") {
    plot.default(0, type = "n", xlim = xlim, ylim = ylim,
                 xlab = "", ylab = "", axes = FALSE, ...)
    if (scale.bar) axis(side = 2, at = at, labels = label, cex.axis=cex)
  }
  tip_labels <- consensus$tip.label
  if (is.expression(consensus$tip.label))
    underscore <- TRUE
  if (!underscore)
    tip_labels <- gsub("_", " ", tip_labels)
  add_tiplabels(xy, tip_labels, direction, adj = adj, font = font, srt = srt,
    cex = cex, col = tip.color, label.offset = label.offset, maxBT = maxBT)

  col <- rep(col, length.out = length(x))
  tiporder <-  1:nTip
  names(tiporder) <- consensus$tip.label

  if (jitter$amount > 0) {
    if (jitter$random) jit <- runif(length(x), -jitter$amount, jitter$amount)
    else jit <- seq(-jitter$amount, jitter$amount, length = length(x))
  }

  for (treeindex in seq_along(x)) {
    tmp <- reorder(x[[treeindex]])

    tmp <- sort_tips(tmp)
    #    if(!compressed) tiporder <- match(tmp$tip.label, consensus$tip.label)
    xy <- plotPhyloCoor(tmp, tip.height = tiporder, direction = direction, ...)
    xx <- xy[, 1]
    yy <- xy[, 2]

    if (horizontal) {
      if (scaleX) xx <- xx / max(xx)
      else xx <- xx #/ maxBT
      if (direction == "rightwards") xx <- xx + (maxBT - max(xx))
      if (jitter$amount > 0) yy <- yy + jit[treeindex]
    }
    else {
      if (scaleX) yy <- yy / max(yy)
      #else yy <- yy
      if (direction == "upwards") yy <- yy + (maxBT - max(yy))
      if (jitter$amount > 0) xx <- xx + jit[treeindex]
    }
    e1 <- tmp$edge[, 1]
    if (type == "cladogram") cladogram.plot(tmp$edge, xx, yy, edge.color =
        adjustcolor(col[treeindex], alpha.f = alpha), edge.width = width,
        edge.lty = lty)
    if (type == "phylogram") {
      Ntip <- min(e1) - 1L
      Nnode <- tmp$Nnode
      phylogram.plot(tmp$edge, Ntip, Nnode, xx, yy, horizontal, edge.color =
        adjustcolor(col[treeindex], alpha.f = alpha), edge.width = width,
        edge.lty = lty)
    }
  }
  L <- list(type = type, font = font, cex = cex,
            adj = adj, srt = srt, #no.margin = no.margin,
            label.offset = label.offset,
            x.lim = xlim, y.lim = ylim, direction = direction,
            tip.color = tip.color, Ntip = nTip #, Nnode = Nnode,
            #root.time = x$root.time, align.tip.label = align.tip.label
            )
  assign("last_plot.phylo", L, envir = .PlotPhyloEnv)
  invisible(x)
}
