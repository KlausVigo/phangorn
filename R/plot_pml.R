#' Plot phylogeny of a pml object
#'
#' \code{plot.pml} is a warapper around \code{plot.phylo} with different default
#' values for unrooted, ultrametric and tip dated phylogenies.
#'
#' @param x an object of class \code{pml} or \code{phyDat}.
#' @param type a character string specifying the type of phylogeny to be drawn;
#' it must be one of "phylogram" (the default), "cladogram", "fan", "unrooted",
#' "radial", "tidy", or any unambiguous abbreviation of these.
#' @param \dots further parameters to be passed to \code{plot.phylo}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{plot.phylo}}, \code{\link{axisPhylo}},
#' \code{\link{add.scale.bar}}
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
#'
#' fit_unrooted <- pml_bb(H3N2, model="JC", rearrangement="none",
#'                        control = pml.control(trace = 0))
#' plot(fit_unrooted, cex=.5)
#'
#' @export
plot.pml <- function(x, type="phylogram", ...){
  type <- match.arg(type, c("phylogram","cladogram", "fan", "unrooted",
                            "radial", "tidy"))
  plot.phylo(x$tree, type=type, ...)
  if(is.rooted(x$tree) && (type %in% c("phylogram","cladogram"))){
    if(!is.null(x$tip.dates) && x$method=="tipdated"){
      root_time <- max(x$tip.dates) - max(node.depth.edgelength(x$tree))
      axisPhylo(root.time = root_time, backward = FALSE)
    }
    else axisPhylo()
  }
  else add.scale.bar()
}
