#' Utility function to plot.phylo
#'
#' cladePar can help you coloring (choosing edge width/type) of clades.
#'
#'
#' @param tree an object of class phylo.
#' @param node the node which is the common ancestor of the clade.
#' @param edge.color see plot.phylo.
#' @param tip.color see plot.phylo.
#' @param edge.width see plot.phylo.
#' @param edge.lty see plot.phylo.
#' @param x the result of a previous call to cladeInfo.
#' @param plot logical, if TRUE the tree is plotted.
#' @param \dots Further arguments passed to or from other methods.
#' @return A list containing the information about the edges and tips.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{plot.phylo}}
#' @keywords plot
#' @examples
#'
#' tree <- rtree(10)
#' plot(tree)
#' nodelabels()
#' x <- cladePar(tree, 12)
#' cladePar(tree, 18, "blue", "blue", x=x, plot=TRUE)
#'
#' @export cladePar
cladePar <- function(tree, node, edge.color = "red", tip.color = edge.color,
              edge.width = 1, edge.lty = "solid", x = NULL, plot = FALSE, ...) {
  if (is.null(x)) {
    m <- max(tree$edge)
    x <- list(
      edge = data.frame(
        color = rep("black", m), width = rep(1, m),
        lty = rep("solid", m)),
      tip = rep("black", length(tree$tip.label))
    )
  }

  lty_str <- c("blank", "solid", "dashed", "dotted", "dotdash", "longdash",
               "twodash")
  if(is.numeric(edge.lty)) edge.lty <- lty_str[edge.lty + 1]

  ind <- Descendants(tree, node, "all")
  x$edge$color[ind] <- edge.color
  x$edge$width[ind] <- edge.width
  x$edge$lty[ind] <- edge.lty
  x[[2]][Descendants(tree, node, "tips")[[1]]] <- tip.color
  if (plot) {
    tree <- reorder(tree)
    plot(tree,
      edge.color = x$edge$color[tree$edge[, 2]],
      edge.width = x$edge$width[tree$edge[, 2]],
      edge.lty = x$edge$lty[tree$edge[, 2]], tip.color = x[[2]], ...
    )
  }
  else {
    return(x)
  }
}
