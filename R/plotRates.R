#' Plot rate heterogeneity across sites
#'
#' \code{plotRates} is useful for exploratory data anlaysis (EDA) tool and
#' and can help judging rate heterogeneity across sites models (RHAS) and see
#' where it limitations are \bibcitep{Ferretti2026}.
#' \code{plotRates} plots the cumulative distribution function based on the
#' parsimony score for the data and the discrete rate model.
#' @param obj an object of class pml
#' @param append logical; if TRUE only add to an existing plot.
#' @param xlim the x limits of the plot.
#' @param \dots Further arguments passed to or from other methods.
#' @param main a main title for the plot.
#' @param col color of the cdf.
#' @param rug logical; if TRUE a \code{\link[graphics]{rug}} is added to the
#' plot.
#' @return \code{discrete.gamma} returns a matrix.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{discrete.gamma}, \link{pml}, \link{stepfun},
#' \link{modelTest}}
#' @references
#' \bibshow{*}
#' @examples
#' data(Laurasiatherian)
#' fit <- pml_bb(Laurasiatherian, "JC+G(4)", rearrangement = "none")
#' plotRates(fit)
#' \dontrun{
#' mt <- modelTest(Laurasiatherian, model=c("JC", "HKY"))
#' plotRates(mt)
#' plotRates(mt, substitution="JC")
#' plotRates(mt, model=c("HKY+G(4)+I", "HKY+G(4)", "JC+G(4)+I", "JC+G(4)"))
#' }
#' @keywords hplot
#' @rdname plotRates
#' @importFrom stats ecdf
#' @importFrom stats plot.ecdf
#' @importFrom graphics rug
#' @export
plotRates <- function(obj, ...) UseMethod("plotRates")

#' @rdname plotRates
#' @export
plotRates.pml <- function(obj, col="blue", main="", rug=FALSE, xlim=NULL,
                          append=FALSE, ...){
  # cdf.color
  pscores <- parsimony(obj$tree, obj$data, site="site")
  if(!is.null(attr(obj$data, "index"))){
    ecdf_pscores <- ecdf( pscores[attr(obj$data, "index")])
  } else{
    tmp <- aggregate(attr(obj$data, "weight"), by=list(pscores), FUN=sum)
    x <- tmp[, 1]
    y <- cumsum(c(0, tmp[,2]))
    y <- y / max(y)
    ecdf_pscores <- stepfun(x, y)
  }
  el <- obj$tree$edge.length * obj$rate
  el <- sum(el)
  w <- obj$w
  g <- obj$g
  inv <- obj$inv
  if(is.null(xlim)) xlim <- c(-0.25, 1.1 * max(c(pscores, g*el)))
  if(!append){
    plot.ecdf(ecdf_pscores, verticals = TRUE, do.points = FALSE, main = main,
              xlim=xlim, ...)
    if(rug) rug(jitter(pscores))
  }
  if (inv > 1e-8) {
    g <- c(0, g)
    w <- c(inv, w)
  }
  plot_gamma_plus_inv(w = w, g = g, append = TRUE, xlim = xlim,
                      edge.length = el, verticals = TRUE, col = col,
                      site.rate = obj$site.rate, ...)
  invisible(obj)
}

#' @rdname plotRates
#' @param substitution substitution model for which the cdf is plotted.
#' @param model vector containing the models to plot.
#' @param pos a character string defining the position of the legend.
#' @param lty the line type, see \code{\link{par}}.
#' @export
plotRates.modelTest <- function(obj, substitution = NULL, model = NULL,
                                append = FALSE, pos = "bottomright",
                                col = NULL, lty = 1, ...){
  if(!is.null(substitution)){
    model <- obj$Model[obj$Substitution %in% substitution]
  }
  if(is.null(model)){
    substitution <- obj$Substitution[which.min(obj$BIC)]
    model <- obj$Model[obj$Substitution == substitution]
  }
  l <- length(model)
  if(is.null(col)){
    if(l==1) col="blue"
    else col <- hcl.colors(l, "Set 2")
  }
  col <- rep(col, length.out = l)
  lty <- rep(lty, length.out = l)
  plotRates(as.pml(obj, model[1]), append=append, col = col[1],
            lty = lty[1], ...)
  if(length(model) > 1) {
    for(i in 2: length(model)){
      plotRates(as.pml(obj, model[i]), append = TRUE, col = col[i],
                lty = lty[i], ...)
    }
  }
  if (!is.null(pos)) legend(pos, legend = model, col = col, lty = lty)
  invisible(obj)
}

