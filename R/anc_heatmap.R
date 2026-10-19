nh_anc <- function(x){
  x <- reorder(x)
  nt <- Ntip(x)
  nn <- Nnode(x)
  res <- integer(nn + nt)
  edge <- x$edge
  e2 <- edge[edge[,2] <= nt,2]
  #  last_anc <- 0
    res[e2[1]] <- 1L
  m <- 2L
  for(i in seq_len(nt - 1)){
    anc_i <- getMRCA(x, c(e2[i], e2[i+1]))
    if(res[anc_i] == 0L){ # != last_anc){
      res[anc_i] <- m
      m <- m + 1L
      }
    res[e2[i+1]] <- m
    m <- m + 1L
    last_anc <- anc_i
  }
  res
}



#' heatmapp
#'
#' Plots a heat map with ancestral states.
#'
#' \code{heatmapp} plots a heatmap, a phylogeny together with the
#' alignment. For ancestral reconstructions it additionally plots the
#' joint distribution or the most likely state of the internal nodes.
#' This eases comparison of ancestors with related Descendant sequences.
#'
#' \code{heatmapp} is based on the ape function s\code{plot.phylo} and
#' \code{image.AAbin} or \code{image.AAbin}.
#' @aliases heatmapp
#' @param x an object of class ancestral or phyDat.
#' @param y an object  of class phyDat.
#' @param use.edge.length a logical indicating whether to use the edge lengths
#' of the phylogeny to draw the branches (the default) or not (if FALSE). This
#' option has no effect if the object of class "phylo" has no ‘edge.length’
#' element.
#' @param align_label a logical value or an integer. If TRUE, the tips are
#' aligned and dotted lines are drawn between the nodes of the tree and the
#' labels.
#' @param clade a node number or label to extract the clade from the tree.
#' @param select a subset of characters, columns in the alignment.
#' @param cex.lab font size of the labels.
#' @param nodes Show ancestral sequences at internal nodes.
#' @param ... Further arguments passed to or from other methods.
#' @return \code{heatmapp} returns silently a \code{x}
#' @importFrom graphics layout mtext
#' @importFrom grDevices dev.flush dev.hold
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{anc_pml}}, \code{\link{plotAnc}},
#' \code{\link[ape]{image.DNAbin}}, \code{\link[ape]{image.AAbin}},
#' \code{\link[ape]{plot.phylo}}
#' @examples
#' data(woodmouse)
#' wm <- as.phyDat(woodmouse)
#' tree <- pratchet(wm)
#' tree <- makeNodeLabel(tree)
#' anc_mp <- anc_pars(tree, wm)
#' op <- par(mar=c(4,1,4,5))
#' heatmapp(anc_mp)
#' heatmapp(tree, woodmouse)
#' par(op)
#' @keywords hplot
#' @rdname heatmapp
#' @export
heatmapp <- function(x, ...) UseMethod("heatmapp")


#' @rdname heatmapp
#' @export
heatmapp.phylo <- function(x, y, use.edge.length = TRUE, align_label = TRUE,
                          clade=NULL, select=NULL, cex.lab=1, nodes = FALSE, ...){
  assert_phylo(x)
  y <- assert_phyDat(y, transform = TRUE)
  if(is.null(x$edge.length)) use.edge.length <- FALSE
  if(nodes) lab <- c(x$tip.label, x$node.label)
  else lab <- x$tip.label
  if(anyNA(match(lab, names(y))))
    stop("labels between tree and alignment mismatch")
  if(!is.null(clade)) x <- extract.clade(x, clade)
  if(!is.null(select)) y <- subset(y, select = select, site.pattern = FALSE)
  nt <- Ntip(x)
  nn <- Nnode(x)
  if(nodes) nl <- nn + nt
  else nl <- nt
  active <- seq_len(nl)
  if(nodes) yy <- nh_anc(x)
  else yy <- node.height(x)
  if(use.edge.length) xx <- node.depth.edgelength(x)
  else {
    xx <- node.depth(x, 2L)
    xx <- max(xx) - xx
  }
  dev.hold()
  on.exit(dev.flush())
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  nf <- layout(matrix(c(1,2),1,2,byrow = TRUE), c(1,3)) #layout.show(nf)
  ind <- rev(lab[order(yy[active])])
  align_label_lty <- 3
  if (is.numeric(align_label)) {
    align_label_lty <- align_label
    align_label <- TRUE
  }
  mar_phylo <- mar_image <- mar <- par("mar")
  mar_phylo[2] <- 1
  mar_phylo[4] <- 0
  par(mar = mar_phylo)
  mar_image[2] <- 0
  plot.default(0, type = "n", xlim = range(xx), ylim = range(yy),
               xlab = "", ylab = "", axes = FALSE)
  usr <- par("usr")
  if(nodes) {
    z <- (usr[4] - usr[3]) / (nt + nn)
    yy_scaled <- usr[3] - (z / 2) + yy * z
  }
  else{
    z <- (usr[4] - usr[3]) / (nt)
    yy_scaled <- usr[3] - (z / 2) + yy * z
  }
  phylogram.plot(x$edge, nt, nn, xx, yy_scaled, TRUE)
  #  if (!use.edge.length || is.ultrametric()) align_label <- FALSE
  if (align_label) {
    xx.tmp <- max(xx)
    segments(xx[active], yy_scaled[active], xx.tmp, yy_scaled[active], lty = align_label_lty)
  }
  par(mar = mar_image)
  y <- y[ind]
  image(y, show.label = FALSE, yaxt = "n", ...)
  n <- length(lab)
  mtext(ind, side = 4, line = 0.1, at = n:1, cex = cex.lab, adj = 0, las = 1)
  invisible(x)
}


#' @rdname heatmapp
#' @export
heatmapp.pml <- function(x, ...){
  assert_pml(x)
  phy <- x$tree
  y <- x$data
  heatmapp(phy, y, nodes=FALSE, ...)
  invisible(x)
}


#' @rdname heatmapp
#' @export
heatmapp.ancestral <- function(x, nodes=TRUE, ...){
  assert_ancestral(x)
  y <- as.phyDat(x)
  phy <- x$tree
  heatmapp(phy, y, nodes=nodes, ...)
  invisible(x)
}


#' @rdname phangorn-internal
#' @export
reference_position <- function(x, pos, ref=1, gap="-"){
  tmp <- as.character(x[ref,])
  ind <- which(tmp == gap)
  y <- rep(1L, length(tmp))
  y[ind] <- 0L
  y <- cumsum(y)
  match(pos, y)
}

