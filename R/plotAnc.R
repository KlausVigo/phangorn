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
## @param tree a tree, i.e. an object of class pml or an object of class
## ancestral.
#' @param x an object of class \code{ancestral}.
## @param site.pattern logical, plot i-th site pattern or i-th site
#' @param i plots the i-th site.
## ,site
## @param which either "pie" or "seqlogo"
#' @param node to plot for which the probabilities should be plotted.
#' @param type a character string specifying the type of phylogeny to be drawn;
#' it must be one of "phylogram" (the default), "cladogram", "fan", "unrooted",
#' "radial", "tidy", or any unambiguous abbreviation of these.
#' @param start start position to plot.
#' @param end end position to plot.
#' @param col a vector containing the colors for all possible states.
#' @param cex.pie a numeric defining the size of the pie graphs.
#' @param pos a character string defining the position of the legend.
#' @param scheme a predefined color scheme. For amino acid options are "Ape_AA",
#' "Zappo_AA", "Clustal", "Polarity" and "Transmembrane_tendency", for
#' nucleotides "Ape_NT" and"RY_NT". Names can be abbreviated.
#' @param \dots Further arguments passed to or from other methods.
#' @returns \code{plotAnc} returns silently x.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{ancestral.pml}}, \code{\link[ape]{plot.phylo}},
#' \code{\link[ape]{image.DNAbin}}, \code{\link[ape]{image.AAbin}}
#' \code{\link[ggseqlogo]{ggseqlogo}}, \code{\link[ape]{edgelabels}}
#' @keywords plot
#' @examples
#'
#' example(NJ)
#' # generate node labels to ensure plotting will work
#' tree <- makeNodeLabel(tree)
#' anc.p <- anc_pars(tree, Laurasiatherian)
#' # plot the third character
#' plotAnc(anc.p, 3, pos="topright")
#' plotSeqLogo(anc.p, node="Node10", 1, 25)
#'
#' data(chloroplast)
#' tree <- pratchet(chloroplast,  maxit=10, trace=0)
#' tree <- makeNodeLabel(tree)
#' anc.ch <- anc_pars(tree, chloroplast)
#' image(as.phyDat(anc.ch)[, 1:25])
#' plotAnc(anc.ch, 21, scheme="Ape_AA", pos="topleft")
#' plotAnc(anc.ch, 21, scheme="Clustal", pos="topleft")
#' plotSeqLogo(anc.ch, node="Node1", 1, 25, scheme="Clustal")
#'
#'
#' data(woodmouse)
#' tree <- pml_bb(woodmouse, "JC", rearrangement="NNI")$tree |> midpoint()
#' woodmouse_aa <- trans(woodmouse, 2) |> as.phyDat()
#' anc_aa <- anc_pars(tree, woodmouse_aa)
#' plot(tree, direction="downwards")
#' add_mutations(anc_aa)
#'
#' @importFrom grDevices hcl.colors
## @importFrom ggseqlogo make_col_scheme ggseqlogo
#' @rdname plot.ancestral
#' @export
plotAnc <- function(x, i = 1, type="phylogram", ..., col = NULL,
                    cex.pie = .5, pos = "bottomright", scheme=NULL) {
  stopifnot(inherits(x, "ancestral"))
  type <- match.arg(type, c("phylogram", "cladogram", "fan", "unrooted",
                            "radial", "tidy"))
  phylo_clado <- type %in% c("phylogram", "cladogram")
  df <- as.data.frame(x)
  data <- x$data
  tree <- x$tree
  subset <- df[,"Site"] == i
  Y <- df[subset & !is.na(subset),]
#  Y <- subset(df, Site==i)
  y <- as.matrix(Y[, -c(1:3)])
  #  y <- y[, -c(1:3)]
  colnames(y) <- gsub("p_", "", colnames(y))
  row.names(y) <- Y$Node
  y <- y[c(tree$tip.label, tree$node.label), ]
  if(is.null(tree$node.label) || any(is.na(match(tree$node.label, rownames(y)))) ||
     is.numeric(tree$node.label))
    tree <- makeNodeLabel(tree)

  if(any(is.na(match(c(tree$tip.label, tree$node.label), rownames(y)))))
    stop("Tree needs nodelabel, which match the labels of the alignment!")
  CEX <- cex.pie
  xrad <- CEX * diff(par("usr")[1:2]) / 50
  levels <- attr(data, "levels")
  nc <- attr(data, "nc")
  if(is.null(scheme) & attr(data, "type")=="AA") scheme <- "Ape_AA"
  if(is.null(scheme) & attr(data, "type")=="DNA") scheme <- "Ape_NT"
  if(!is.null(scheme)){
    scheme <- match.arg(scheme, c("Ape_AA", "Zappo_AA", "Clustal", "Polarity",
                                  "Transmembrane_tendency", "Ape_NT", "RY_NT"))
    sc <- get(scheme, environment(ace))
    if(has_gap_state(data) && attr(data, "type")=="AA"){
      sc$properties <- c(sc$properties, Gap="-")
      sc$color <- c(sc$color, "#FFFFFF")
    }
    if(attr(data, "type")=="DNA"){
      ind <- match("n", names(sc$properties))
      if(!is.na(ind)){
        sc$properties <- sc$properties[-ind]
        sc$color <- sc$color[-ind]
      }
    }
    P <- getTransition(sc, levels)
    y <- y %*% P
    levels <- colnames(P)
    col <- sc$col
    nc <- ncol(y)
  }
  plot(tree, label.offset = 1.1 * xrad, plot = FALSE, type=type, ...)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  XX <- lastPP$xx
  YY <- lastPP$yy
  xrad <- CEX * diff(lastPP$x.lim * 1.1) / 50
  par(new = TRUE)
  plot(tree, label.offset = 1.1 * xrad, plot = TRUE, type=type, ...)
  if (is.null(col)) col <- hcl.colors(nc) #rainbow(nc)
  if(!is.null(names(col))) col <- col[attr(data, "levels")]
  if (length(col) != nc) {
    warning("Length of color vector differs from number of levels!")
  }
  BOTHlabels(
    pie = y, XX = XX, YY = YY, adj = c(0.5, 0.5), frame = "rect", pch = NULL,
    sel = seq_along(XX), thermo = NULL, piecol = col, col = "black",
    bg = "lightblue", horiz = FALSE, width = NULL, height = NULL, cex = cex.pie
  )
#  if(legend) legend(pos, legend=levels, pch=21, pt.bg = col)
  if (!is.null(pos)) legend(pos, legend=levels, pch=21, pt.bg = col)
  invisible(x)
}


my_ggseqlogo <-function (data, facet = "wrap", scales = "free_x", ncol = NULL,
          nrow = NULL, start=NULL, end=NULL, ...)
{
  x <- ggseqlogo::geom_logo(data = data, ...)
  x[[2]] <- ggplot2::scale_x_continuous(limits = c(start-0.5, end+.5) ,
                               breaks=pretty(seq(start, end)))
  p <- ggplot2::ggplot() + x + ggseqlogo::theme_logo()
  if (!"list" %in% class(data)) return(p)
  facet <- match.arg(facet, c("grid", "wrap"))
  if (facet == "grid") {
    p <- p + ggplot2::facet_grid(~seq_group, scales = scales)
  }
  else if (facet == "wrap") {
    p <- p + ggplot2::facet_wrap(~seq_group, scales = scales, nrow = nrow,
                                 ncol = ncol)
  }
  return(p)
}


#' @rdname plot.ancestral
## @importFrom ggplot2 scale_x_continuous ggplot facet_grid facet_wrap
## @importFrom ggseqlogo geom_logo theme_logo
#' @returns \code{plotSeqLogo} returns a ggplot object.
#' @export
plotSeqLogo <- function(x, node=getRoot(x$tree), start=1, end=10,
                        scheme="Ape_NT", ...){
  stopifnot(inherits(x, "ancestral"))
  chk <- requireNamespace("ggseqlogo", quietly = TRUE)
  if (!chk) stop("package ggseqlogo needs to be not installed!\n")
  type <- attr(x$data, "type")
  tree <- x$tree
  df <- as.data.frame(x)
  nodes <- c(tree$tip.label, tree$node.label)
  if(is.numeric(node)) node <- nodes[node]
  subset <- df[,"Node"] == node
  X <- df[subset & !is.na(subset),]
#  X2 <- subset(df, subset=Node==node)
  end <- min(end, nrow(X))
  X <- X[seq_len(end), , drop=FALSE]  # creating the whole plot is slow
  X <- t(as.matrix(X[, -c(1:3)]))
  tmp <- gsub("p_", "", rownames(X))
  lev <- rownames(X) <- toupper(tmp)
  if(is.null(scheme) & type=="AA") scheme <- "Ape_AA"
  if(is.null(scheme) & type=="DNA") scheme <- "Ape_NT"
  if(!is.null(scheme)){
    scheme <- match.arg(scheme, c("Ape_AA", "Zappo_AA", "Clustal", "Polarity",
                                  "Transmembrane_tendency", "Ape_NT", "RY_NT"))
    sc <- get(scheme, environment(ace))
    if(has_gap_state(x$data) && type=="AA"){
      sc$properties <- c(sc$properties, Gap="-")
      sc$color <- c(sc$color, "#FFFFFF")
    }
    l <- lengths(sc$properties)
    SC <- ggseqlogo::make_col_scheme(chars = toupper(unlist(sc$properties)),
                          groups = rep(names(sc$properties), l),
                          cols=rep(sc$color, l))

  }
  else SC <- ggseqlogo::make_col_scheme(chars=lev, cols= hcl.colors(length(lev)))
  my_ggseqlogo(X, col_scheme=SC, method='p', start=start, end=end)
}


#' @rdname plot.ancestral
#' @param frame a character string specifying the kind of frame to be printed
#' around the text. See \code{\link[ape]{edgelabels}}.
#' @returns \code{add_mutations} adds the position and and changes of possible
#' mutations to a phylogeny.
#' @export
add_mutations <- function(x, frame="none", ...){
  y <- rbind(x$data, x$state)
  mt <- map_mutations(x$tree, y)
  ind <- which(lengths(mt)>0)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  edge <- lastPP$edge
  pos <- match(ind, edge[,2])
  edgelabels(lapply(mt[ind], paste, collapse=" "), pos, frame=frame, ...)
}


