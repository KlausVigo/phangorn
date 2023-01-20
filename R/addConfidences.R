#' Compare splits and add support values to an object
#'
#' Add support values to a \code{splits}, \code{phylo} or \code{networx}
#' object.
#'
#' @param x an object of class \code{splits}, \code{phylo} or \code{networx}
#' @param y an object of class \code{splits}, \code{phylo}, \code{multiPhylo}
#' or \code{networx}
#' @param rooted logial, if FALSE bipartitions are considered, if TRUE clades.
#' @param ...  Further arguments passed to or from other methods.
#' @param label_y label of y matched on x. Will be usually of
#' length(as.splits(x)).
#' @param type should labels returned for edges (in \code{networx}) or splits.
#' @param nomatch default value if no match between x and y is found.
#' @return The object \code{x} with added bootstrap / MCMC support values.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{as.splits}}, \code{\link{as.networx}},
#' \code{\link{RF.dist}}, \code{\link{plot.phylo}}
#' @references Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W.
#' (2017), Intertwining phylogenetic trees and networks.
#' \emph{Methods Ecol Evol}.\bold{8}, 1212--1220. doi:10.1111/2041-210X.12760
#' @keywords cluster
#' @examples
#'
#' data(woodmouse)
#' woodmouse <- phyDat(woodmouse)
#' tmpfile <- normalizePath(system.file(
#'              "extdata/trees/RAxML_bootstrap.woodmouse", package="phangorn"))
#' boot_trees <- read.tree(tmpfile)
#'
#' dm <- dist.ml(woodmouse)
#' tree <- upgma(dm)
#' nnet <- neighborNet(dm)
#'
#' tree <- addConfidences(tree, boot_trees)
#' nnet <- addConfidences(nnet, boot_trees)
#'
#' plot(tree, show.node.label=TRUE)
#' plot(nnet, show.edge.label=TRUE)
#'
#' @rdname addConfidences
#' @export
addConfidences <- function(x, y, ...) UseMethod("addConfidences")


# some function to add confidences on splits if trees have different taxa
# used in addConfidences.splits
addConfidencesMultiPhylo <- function(spl, trees) {
  fun <- function(spl, intersect_labels) {
    spl2 <- spl
    index <- match(attr(spl, "labels"), intersect_labels)
    attr(spl2, "labels") <- intersect_labels
    for (i in seq_along(spl2)) {
      spl2[[i]] <- sort(na.omit(index[spl[[i]]]))
    }
    l_spl <- lengths(spl2)
    l <- length(intersect_labels)
    ind <- which((l_spl > 1) & (l_spl < (l - 1L)))
    if (length(ind) == 0) return(NULL)
    list(spl = spl2[ind], index = ind)
  }

  spl_labels <- attr(spl, "labels")
  zaehler <- numeric(length(spl))
  nenner <- numeric(length(spl))
  for (i in seq_along(trees)) {
    intersect_labels <- intersect(trees[[i]]$tip.label, spl_labels)
    if (length(intersect_labels) > 3) {
      tmp <- fun(spl, intersect_labels)
      if (!is.null(tmp)) {
        tree_spl <- as.splits(trees[[i]])
        if (!identical(intersect_labels, trees[[i]]$tip.label))
          tree_spl <- fun(tree_spl, intersect_labels)[[1]]
        comp <- compatible_2(as.bitsplits(tmp[[1]]), as.bitsplits(tree_spl))
        ind <- tmp$index
        zaehler[ind] <- zaehler[ind] + comp
        nenner[ind] <- nenner[ind] + 1L
      }
    }
  }
  confidences <- zaehler / nenner
  attr(spl, "confidences") <- confidences
  spl
}


#' @export
addConfidences.splits <- function(x, y, scaler = 1, rooted=FALSE, ...) {
  if (hasArg(add))
    add <- list(...)$add
  else add <- FALSE

  tiplabel <- attr(x, "labels")
  nTips <- length(tiplabel)
  #    x = addTrivialSplits(x)
  if (inherits(y, "phylo")) {
    ind <- match(tiplabel, y$tip.label)
    if (any(is.na(ind)) | length(tiplabel) != length(y$tip.label))
      stop("trees have different labels")
    y$tip.label <- y$tip.label[ind]
    ind2 <- match(seq_along(ind), y$edge[, 2])
    y$edge[ind2, 2] <- order(ind)
  }
  if (inherits(y, "multiPhylo")) {
    if (inherits(try(.compressTipLabel(y), TRUE), "try-error")) {
      res <- addConfidencesMultiPhylo(x, y)
      return(res)
    }
  }
  spl <- as.splits(y)
  spl <- changeOrder(spl, tiplabel)
  if(!rooted){
    spl <- SHORTwise(spl)
    ind <- match(SHORTwise(x), spl)
  }
  else ind <- match(x, spl)
  pos <-  which(!is.na(ind))
  confidences <- rep(NA_real_, length(x)) # numeric(length(x))
  if(is.numeric(attr(spl, "confidences")))
    confidences[pos] <- attr(spl, "confidences")[ind[pos]] * scaler
  else confidences[pos] <- attr(spl, "confidences")[ind[pos]]
  if (add == TRUE) confidences <- paste(prettyNum(attr(x, "confidences")),
                                        prettyNum(confidences * scaler), sep = "/")
  attr(x, "confidences") <- confidences
  x
}


#' @export
addConfidences.networx <- function(x, y, scaler = 1, ...) {
  spl <- x$splits
  spl <- addConfidences(spl, y, scaler = scaler, ...)
  x$splits <- spl
  x
}

#' @rdname addConfidences
#' @export
addConfidences.phylo <- function(x, y, rooted=FALSE, ...) {
  #    call <- x$call
  if (hasArg(as.is))
    as.is <- list(...)$as.is
  else as.is <- TRUE
  nTips <- length(x$tip.label)
  spl <- as.splits(x)
  if(!rooted) spl <- SHORTwise(spl)
  conf <- attr(addConfidences(spl, y, rooted=rooted, ...), "confidences")
  l <- lengths(spl)
  if (is.character(conf)) as.is <- TRUE
  ind <- (l == 1L) | (l == (nTips - 1L)) | (l == nTips)
  conf[ind == TRUE] <- NA_real_
  nTips <- length(x$tip.label)
  if (!as.is) conf <- conf * 100
  x$node.label <- conf[-c(1:nTips)]
  x
}


#' @export
addConfidences.multiPhylo <- function(x, y, ...) {
  x <- .uncompressTipLabel(x)
  x <- unclass(x)
  x <- lapply(x, addConfidences, y, ...)
  class(x) <- "multiPhylo"
  x
}

#' @rdname addConfidences
#' @export
presenceAbsence <- function(x, y) {
  spl <- as.splits(y)
  l <- length(spl)
  attr(spl, "confidences") <- rep(1, l)
  addConfidences(x, y)
}


#' @rdname addConfidences
#' @export
createLabel <- function(x, y, label_y, type = "edge", nomatch = NA) {
  spl_x <- as.splits(x)
  if (inherits(x, "phylo", TRUE) == 1) spl_x <- spl_x[x$edge[, 2]]
  spl_y <- as.splits(y)
  if (inherits(y, "phylo", TRUE) == 1) spl_y <- spl_y[y$edge[, 2]]

  tiplabel <- attr(spl_x, "labels")
  nTips <- length(tiplabel)

  spl_y <- changeOrder(spl_y, tiplabel)
  spl_y <- SHORTwise(spl_y)

  ind <- match(SHORTwise(spl_x), spl_y)
  pos <-  which(!is.na(ind))

  res <- rep(nomatch, length(spl_x))

  if (length(label_y) == 1L) label_y <- rep(label_y, length(spl_y))
  res[pos] <- label_y[ind[pos]]
  if (type == "edge" && inherits(x, "networx")) {
    return(res[x$splitIndex])
  }
  res
}
