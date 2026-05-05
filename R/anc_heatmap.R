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


#' anc_heatmap
#'
#' Plots a heat map with ancestral states.
#'
#' \code{anc_heatmap} plot the joint distribution or the most likely state.

#' @aliases anc_heatmap
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
#' @param ... Further arguments passed to or from other methods.
#' @importFrom graphics layout
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{anc_pml}}, \code{\link{plotAnc}}
#' @examples
#' data(woodmouse)
#' wm <- as.phyDat(woodmouse)
#' tree <- pratchet(wm)
#' tree <- makeNodeLabel(tree)
#' anc_mp <- anc_pars(tree, wm)
#' anc_heatmap(anc_mp)
#' @keywords hplot
#' @export
anc_heatmap <- function(x, y=NULL, use.edge.length = FALSE,
                        align_label = TRUE, clade=NULL, select=NULL, ...){
  if(!is.null(y)){
    if(inherits(y, "AAbin") || inherits(y, "DNAbin")) y <- as.phyDat(y)
    }
  if(inherits(x, "ancestral")){
    if(is.null(y)) y <- as.phyDat(x)
    x <- x$tree
    }
  if(is.null(x$edge.length)) use.edge.length <- FALSE
  if(anyNA(match(c(x$tip.label, x$node.label), names(y))))
    stop("labels between treee and alignment mismatch")

  if(!is.null(clade)) x <- extract.clade(x, clade)
  if(!is.null(select)) y <-subset(y, select=select, site.pattern = FALSE)

  nt <- Ntip(x)
  nn <- Nnode(x)
  yy <- nh_anc(x)
  if(use.edge.length) xx <- node.depth.edgelength(x)
  else {
    xx <- node.depth(x, 2L)
    xx <- max(xx) - xx
  }
  op <- par()
  nf <- layout(matrix(c(1,2),1,2,byrow = TRUE), c(1,2))
  #layout.show(nf)
  #text(x=3, y=xx, labels=lab)
    lab <- c(x$tip.label, x$node.label)
  ind <- rev(lab[order(yy)])

  if (is.numeric(align_label)) {
    align_label_lty <- align_label
    align_label <- TRUE
    } else { # assumes is.logical(align.tip.labels) == TRUE
      if (align_label) align_label_lty <- 3
      }
  plot.default(0, type = "n", xlim = range(xx), ylim = range(yy),
                             xlab = "", ylab = "", axes = FALSE)
  usr <- par("usr")
  z <- (usr[4] - usr[3]) / (nt + nn)
  yy_scaled <- usr[3] - (z / 2) + yy * z

  mar_phylo <- mar <- par("mar")
  mar_phylo[2] <- 1
  mar_phylo[4] <- 0
  par(mar=mar_phylo)

  phylogram.plot(x$edge, nt, nn, xx, yy_scaled, TRUE)
  if (align_label) {
    if (!use.edge.length || is.ultrametric(x)) align_tip_label <- FALSE
    xx.tmp <- max(xx)
    segments(xx, yy_scaled, xx.tmp, yy_scaled, lty = align_label_lty)
  }
  par(mar=mar)

  y <- y[ind]
  image(y, ...)
  par("mfrow"=op$mfrow)
#  par(op)
  invisible(x)
}


#' @rdname phangorn-internal
#' @export
reference_position <- function(x, pos, ref=1, gap="-"){
  tmp <- as.character(x[ref,])
  ind <- which(tmp==gap)
  y <- rep(1L, length(tmp))
  y[ind] <- 0L
  y <- cumsum(y)
  match(pos, y)
}


#example(anc_pml)
