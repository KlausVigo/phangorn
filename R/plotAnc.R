getTransition <- function(scheme, levels){
  l <- length(scheme$properties)
  P <- matrix(0, length(levels), l,
              dimnames = list(levels, names(scheme$properties)))
  for(i in seq_along(scheme$properties)){
    ind <- match(scheme$properties[[i]], levels)
    P[ind,i] <- 1
  }
  P
}


#' Plot ancestral character on a tree
#'
#' \code{plotAnc} plots a phylogeny and adds character to the nodes. Either
#' takes output from  \code{ancestral.pars} or \code{ancestral.pml} or from an
#' alignment where there are node labels in the tree match the constructed
#' sequences in the alignment.
#'
#' For further details see vignette("Ancestral").
#'
#' @param tree a tree, i.e. an object of class pml
#' @param data an object of class \code{phyDat} or \code{ancestral}.
#' @param site.pattern logical, plot i-th site pattern or i-th site
#' @param i plots the i-th site of the \code{data}.
#' @param col a vector containing the colors for all possible states.
#' @param cex.pie a numeric defining the size of the pie graphs.
#' @param pos a character string defining the position of the legend.
#' @param scheme color scheme, for amino acids this can be "Ape_AA",
#' "Clustal_AA", "Hydrophobicity_AA" and "Zappo_AA" are available.
#' @param \dots Further arguments passed to or from other methods.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{ancestral.pml}}, \code{\link[ape]{plot.phylo}}
#' @keywords plot
#' @examples
#'
#' example(NJ)
#' # generate node labels to ensure plotting will work
#' tree <- makeNodeLabel(tree)
#' anc.p <- ancestral.pars(tree, Laurasiatherian)
#' # plot the third character
#' plotAnc(tree, anc.p, 3)
#'
#' data(chloroplast)
#' tree <- pratchet(chloroplast,  maxit=10, trace=0)
#' tree <- makeNodeLabel(tree)
#' anc.ch <- ancestral.pars(tree, chloroplast)
#' image(chloroplast[, 1:25])
#' plotAnc(tree, anc.ch, 21, scheme="Ape_AA")
#' plotAnc(tree, anc.ch, 21, scheme="Clustal_AA")
#' @importFrom grDevices hcl.colors
#' @export
plotAnc <- function(tree, data, i = 1, site.pattern = FALSE, col = NULL,
                    cex.pie = .5, pos = "bottomright", scheme=NULL,
                    ...) {
  stopifnot(inherits(data, "phyDat"))
  y <- subset(data, select = i, site.pattern = site.pattern)
  if(is.null(tree$node.label) || any(is.na(match(tree$node.label, names(y)))) ||
     is.numeric(tree$node.label))
    tree <- makeNodeLabel(tree)
  if(any(is.na(match(c(tree$tip.label, tree$node.label), names(y)))))
    stop("Tree needs nodelabel, which match the labels of the alignment!")
  y <- y[c(tree$tip.label, tree$node.label),]
  CEX <- cex.pie
  xrad <- CEX * diff(par("usr")[1:2]) / 50
  levels <- attr(data, "levels")
  nc <- attr(data, "nc")
  if(inherits(data, "ancestral")){
    y <- matrix(unlist(y[]), ncol = nc, byrow = TRUE)
  } else y <- attr(data, "contrast")[unlist(y),]
  if(!is.null(scheme)){
    sc <- get(scheme, environment(pml))
    P <- getTransition(sc, levels)
    y <- y %*% P
    levels <- colnames(P)
    col <- sc$col
    nc <- ncol(y)
  }
  #  l <- dim(y)[1]
  #  dat <- matrix(0, l, nc)
  #  for (i in 1:l) dat[i, ] <- y[[i]]
  plot(tree, label.offset = 1.1 * xrad, plot = FALSE, ...)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  XX <- lastPP$xx
  YY <- lastPP$yy
  xrad <- CEX * diff(lastPP$x.lim * 1.1) / 50
  par(new = TRUE)
  plot(tree, label.offset = 1.1 * xrad, plot = TRUE, ...)
  if (is.null(col)) col <-  hcl.colors(nc) #rainbow(nc)
  if (length(col) != nc) {
    warning("Length of color vector differs from number of levels!")
  }
  BOTHlabels(
    pie = y, XX = XX, YY = YY, adj = c(0.5, 0.5), frame = "rect", pch = NULL,
    sel = seq_along(XX), thermo = NULL, piecol = col, col = "black",
    bg = "lightblue", horiz = FALSE, width = NULL, height = NULL, cex = cex.pie
  )
  if (!is.null(pos)) legend(pos, legend=levels, pch=19, col = col)
}
