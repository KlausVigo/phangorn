#' Plot phylogeny of a pml object
#'
#' \code{plot.pml} is a wrapper around \code{plot.phylo} with different default
#' values for unrooted, ultrametric and tip dated phylogenies.
#'
#' @param x an object of class \code{pml} or \code{phyDat}.
#' @param type a character string specifying the type of phylogeny to be drawn;
#' it must be one of "phylogram" (the default), "cladogram", "fan", "unrooted",
#' "radial", "tidy", or any unambiguous abbreviation of these.
#' @param direction a character string specifying the direction of the tree.
#' Four values are possible: "rightwards" (the default), "leftwards", "upwards",
#' and "downwards".
#' @param \dots further parameters to be passed to \code{plot.phylo}.
#' @return \code{plot.pml} returns invisibly a list with arguments dexcribing
#' the plot. For further details see the \code{plot.phylo}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{plot.phylo}}, \code{\link[ape]{axisPhylo}},
#' \code{\link[ape]{add.scale.bar}}
#' @keywords plot
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
plot.pml <- function(x, type="phylogram", direction = "rightwards", ...){
  type <- match.arg(type, c("phylogram","cladogram", "fan", "unrooted",
                            "radial", "tidy"))
  tree <- x$tree
  extras <- match.call(expand.dots = FALSE)$...
  cex <- ifelse(is.null(extras$cex), par("cex"), extras$cex)
  cex.axis <- ifelse(is.null(extras$cex.axis), cex, extras$cex.axis)
  if(!is.rooted(tree) && (type != "unrooted") ) tree <- midpoint(tree)
  L <- plot.phylo(tree, type=type, direction=direction, ...)
  if(is.rooted(tree) && (type %in% c("phylogram","cladogram"))){
    direction <- match.arg(direction, c("rightwards", "leftwards", "upwards",
                                        "downwards"))
    side <-   switch(direction,
                     rightwards = 1,
                     leftwards = 1,
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
  }
  else add.scale.bar(cex=cex)
  if(!is.null(x$bs)) add_support(tree, x$bs, cex=cex)
  invisible(L)
}
