#' Lento plot
#'
#' The lento plot represents support and conflict of splits/bipartitions.
#'
#'
#' @param obj an object of class phylo, multiPhylo or splits
#' @param xlim graphical parameter
#' @param ylim graphical parameter
#' @param main graphical parameter
#' @param sub graphical parameter
#' @param xlab graphical parameter
#' @param ylab graphical parameter
#' @param bipart plot bipartition information.
#' @param trivial logical, whether to present trivial splits (default is
#' FALSE).
#' @param col color for the splits / bipartition.
#' @param \dots Further arguments passed to or from other methods.
#' @return lento returns a plot.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{as.splits}, \link{hadamard}}
#' @references Lento, G.M., Hickson, R.E., Chambers G.K., and Penny, D. (1995)
#' Use of spectral analysis to test hypotheses on the origin of pinninpeds.
#' \emph{Molecular Biology and Evolution}, \bold{12}, 28-52.
#' @keywords cluster plot
#' @examples
#'
#' data(yeast)
#' yeast.ry <- acgt2ry(yeast)
#' splits.h <- h2st(yeast.ry)
#' lento(splits.h, trivial=TRUE)
#'
#' @export lento
lento <- function(obj, xlim = NULL, ylim = NULL, main = "Lento plot",
                  sub = NULL, xlab = NULL, ylab = NULL, bipart = TRUE,
                  trivial = FALSE, col = rgb(0, 0, 0, .5), ...) {
  if (inherits(obj, "phylo")) {
    if (inherits(obj, "phylo", TRUE) == 1) obj <- as.splits(obj)[obj$edge[, 2]]
    obj <- as.splits(obj)
  }
  if (inherits(obj, "multiPhylo"))
    obj <- as.splits(obj)
  labels <- attr(obj, "labels")
  l <- length(labels)
  if (!trivial) {
    triv <- lengths(obj)
    ind <- logical(length(obj))
    ind[(triv > 1) & (triv < (l - 1))] <- TRUE
    if (length(col) == length(obj)) col <- col[ind]
    obj <- obj[ind]
  }
  CM <- compatible(obj)
  support <- attr(obj, "weights")
  if (is.null(support))
    support <- rep(1, length(obj))
  conflict <- -as.matrix(CM) %*% support
  n <- length(support)
  if (is.null(ylim)) {
    eps <- (max(support) - min(conflict)) * 0.05
    ylim <- c(min(conflict) - eps, max(support) + eps)
  }
  if (is.null(xlim)) {
    xlim <- c(0, n + 1)
  }

  ord <- order(support, decreasing = TRUE)
  support <- support[ord]
  conflict <- conflict[ord]
  if (length(col) == length(obj)) col <- col[ord]
  plot.new()
  plot.window(xlim, ylim)
  title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
  segments(0:(n - 1), support, y1 = conflict, ...)
  segments(1:n, support, y1 = conflict, ...)
  segments(0:(n - 1), support, x1 = 1:n, ...)
  segments(0:(n - 1), conflict, x1 = 1:n, ...)
  abline(h = 0)
  axis(2, ...)
  aty <- diff(ylim) / (l + 1)
  at <- min(ylim) + (1:l) * aty
  if (bipart) {
    Y <- rep(at, n)
    X <- rep( (1:n) - .5, each = l)
    Circles <- matrix(1, l, n)
    for (i in 1:n) Circles[obj[[ord[i]]], i] <- 19
    col <- rep(col, each = l)
    text(x = n + .1, y = at, labels, pos = 4, ...)
    points(X, Y, pch = as.numeric(Circles), col = col, ...)
  }
  invisible(list(support = cbind(support, conflict), splits = obj[ord]))
}
