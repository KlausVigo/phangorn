#' Maximum agreement subtree
#'
#' \code{mast} computes the maximum agreement subtree (MAST).
#'
#' The code is derived from the code example in Valiente (2009).
## for the original code see \url{http://www.cs.upc.edu/~valiente/comput-biol/}.
#' The version for the unrooted trees is much slower.
#'
#' @param x a tree, i.e. an object of class \code{phylo}.
#' @param y a tree, i.e. an object of class \code{phylo}.
#' @param tree a logicical, if TRUE returns a tree other wise the tip labels
#' of the the maximum agreement subtree.
#' @param rooted logical if TRUE treats trees as rooted otherwise unrooted.
#' @return \code{mast} returns a vector of the tip labels in the MAST.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com} based on code of Gabriel Valiente
#' @seealso \code{\link{SPR.dist}}
#' @references
#' G. Valiente (2009). \emph{Combinatorial Pattern Matching Algorithms in Computational Biology using Perl and R}. Taylor & Francis/CRC Press
#'
#'
#' @keywords cluster
#' @examples
#' tree1 <- rtree(100)
#' tree2 <- rSPR(tree1, 5)
#' tips <- mast(tree1, tree2)
#'
#' @rdname mast
#' @export
mast <- function(x, y, tree = TRUE, rooted = TRUE) {

  if (!is.rooted(x) | !is.rooted(y)) {
    rooted <- FALSE
  }

  shared_tips <- intersect(x$tip.label, y$tip.label)
  if (length(shared_tips) < length(x$tip.label))
    x <- drop.tip(x, setdiff(x$tip.label, shared_tips))
  if (length(shared_tips) < length(y$tip.label))
    y <- drop.tip(y, setdiff(y$tip.label, shared_tips))

  # make order of labels the same
  y <- .compressTipLabel(c(y), x$tip.label)[[1]]
  x <- reorder(x, "postorder")

  if (rooted) res <- mast.fit(x, y)
  else {
    x <- reorder(x, "postorder")
    bipart_x <- bipartCPP(x$edge, Ntip(x))[[2]]
    x <- root(x, bipart_x, resolve.root = TRUE)
    x <- reorder(x, "postorder")

    y <- unroot(y)
    y <- reorder(y, "postorder")

    res <- NULL

    bipart_y <- bipartCPP(y$edge, Ntip(y))
    for (i in 2:length(bipart_y)) {
      y <- root(y, bipart_y[[i]], resolve = TRUE)
      tmp <- mast.fit(x, y)
      if (length(tmp) > length(res)) res <- tmp
    }
  }
  if (tree) res <- drop.tip(x, setdiff(x$tip.label, res))
  res
}


mast.fit <- function(x, y) {
  y <- reorder(y, "postorder")

  po1 <- c(x$edge[, 2], x$edge[nrow(x$edge), 1])
  po2 <- c(y$edge[, 2], y$edge[nrow(y$edge), 1])

  nTips <- length(x$tip.label)
  # vielleicht ausserhalb
  p_vec_1 <- Ancestors(x, 1L:max(x$edge))  # nTips
  p_vec_2 <- Ancestors(y, 1L:max(y$edge))  # nTips

  # vielleicht ausserhalb
  CH1 <- allChildren(x)
  CH2 <- allChildren(y)

  m <- matrix(list(), nrow = length(po1), ncol = length(po2))

  for (i in seq_len(nTips)) {
    m[i, i] <- c(i)
    m[cbind(i, p_vec_2[[i]])] <- c(i)
    m[cbind(p_vec_1[[i]], i)] <- c(i)
  }

  for (i in po1) {
    for (j in po2) {
      if (i > nTips & j > nTips) {
        l1 <- CH1[[i]][1]
        r1 <- CH1[[i]][2]
        l2 <- CH2[[j]][1]
        r2 <- CH2[[j]][2]
        mm <- c(m[[l1, l2]], m[[r1, r2]])
        if (length(m[[l1, r2]]) + length(m[[r1, l2]]) > length(mm))
          mm <- c(m[[l1, r2]], m[[r1, l2]])
        if (length(m[[i, l2]]) > length(mm))
          mm <- m[[i, l2]]
        if (length(m[[i, r2]]) > length(mm))
          mm <- m[[i, r2]]
        if (length(m[[l1, j]]) > length(mm))
          mm <- m[[l1, j]]
        if (length(m[[r1, j]]) > length(mm))
          mm <- m[[r1, j]]
        if (!is.null(mm)) m[[i, j]] <- mm
      }
    }
  }
  x$tip.label[m[[i, j]]]
}
