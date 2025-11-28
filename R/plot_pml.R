#' Plot phylogeny of a pml object
#'
#' \code{plot.pml} is a wrapper around \code{plot.phylo} with different default
#' values for unrooted, ultrametric and tip dated phylogenies.
#'
#' @param x an object of class \code{pml}.
#' @param type a character string specifying the type of phylogeny to be drawn;
#' it must be one of "phylogram" (the default), "cladogram", "fan", "unrooted",
#' "radial", "tidy", or any unambiguous abbreviation of these.
#' @param direction a character string specifying the direction of the tree.
#' Four values are possible: "rightwards" (the default), "leftwards", "upwards",
#' and "downwards".
#' @param adj	one or two numeric values specifying the horizontal and vertical
#' justification of the text or symbols of the support values.
#' @param method either "FBP" the classical bootstrap (default), "TBE"
#' (transfer bootstrap) or "MCC" for assigning clade credibilities.
#' @param digits integer indicating the number of decimal places.
#' @param \dots further parameters to be passed to \code{plot.phylo}.
#' @return \code{plot.pml} returns the \code{pml} object x.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{plot.phylo}}, \code{\link[ape]{axisPhylo}},
#' \code{\link[ape]{add.scale.bar}}
#' @keywords hplot
#' @examples
#' fdir <- system.file("extdata/trees", package = "phangorn")
#' tmp <- read.csv(file.path(fdir,"H3N2_NA_20.csv"))
#' H3N2 <- read.phyDat(file.path(fdir,"H3N2_NA_20.fasta"), format="fasta")
#' dates <- setNames(tmp$numdate_given, tmp$name)
#'
#' fit_td <- pml_bb(H3N2, model="JC", method="tipdated", tip.dates=dates,
#'                  rearrangement="none", control = pml.control(trace = 0))
#' plot(fit_td, show.tip.label = FALSE)
#' # Same as:
#' # root_time <- max(dates) - max(node.depth.edgelength(fit_td$tree))
#' # plot(fit_td$tree, show.tip.label = FALSE)
#' # axisPhylo(root.time = root_time, backward = FALSE)
#' plot(fit_td, show.tip.label = FALSE, direction="up")
#'
#' fit_unrooted <- pml_bb(H3N2, model="JC", rearrangement="none",
#'                        control = pml.control(trace = 0))
#' plot(fit_unrooted, cex=.5)
#'
#' @export
plot.pml <- function(x, type="phylogram", direction = "rightwards",
                     ..., adj = NULL, digits=2, method="FBP"){
  type <- match.arg(type, c("phylogram","cladogram", "fan", "unrooted",
                            "radial", "tidy"))
  tree <- x$tree
  extras <- match.call(expand.dots = FALSE)$...
  cex <- ifelse(is.null(extras$cex), par("cex"), extras$cex)
  cex.axis <- ifelse(is.null(extras$cex.axis), cex, extras$cex.axis)
  if(!is.rooted(tree) && (type != "unrooted") ) tree <- midpoint(tree)
  tree <- ladderize(tree)
  L <- plot.phylo(tree, type=type, direction=direction, ...)
  if(is.rooted(tree) && (type %in% c("phylogram","cladogram"))){
    direction <- match.arg(direction, c("rightwards", "leftwards", "upwards",
                                        "downwards"))
    side <-   switch(direction,
                     "rightwards" = 1,
                     "leftwards" = 1,
                     "upwards" = 2,
                     "downwards" = 2)
    if(!is.null(x$tip.dates) && x$method=="tipdated"){
      root_time <- max(x$tip.dates) - max(node.depth.edgelength(x$tree))
      axisPhylo(side, root.time = root_time, backward = FALSE,
                cex.axis=cex.axis)
    }
    else if(!is.null(x$method) && x$method=="ultrametric")
      axisPhylo(side, cex.axis=cex.axis)
    else add.scale.bar(cex=cex)
    if(!is.null(x$bs)) {
      if(is.null(adj)){
        adj <- c(0.5, 0)
        if(side==2) adj <- c(1, 0.5)
      }
      add_support(tree, x$bs, cex=cex, adj=adj, method=method, digits=digits)
    }
  }
  else{
    add.scale.bar(cex=cex)
    if(!is.null(x$bs)) {
      if(is.null(adj)) adj <- c(0.5, 0.5)
      add_support(tree, x$bs, cex=cex, adj=adj, method=method, digits=digits)
    }
  }
  invisible(x)
}
#' @srrstats {G2.3, G2.3a} in lines: 46, 55
