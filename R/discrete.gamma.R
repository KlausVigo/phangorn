#' Discrete Gamma and Beta distribution
#'
#' \code{discrete.gamma} internally used for the likelihood computations in
#' \code{pml} or \code{optim.pml}. It is useful to understand how it works
#' for simulation studies or in cases where .
#'
#' These functions are exported to be used in different packages so far only in
#' the package coalescentMCMC, but are not intended for end user. Most of the
#' functions call C code and are far less forgiving if the import is not what
#' they expect than \code{pml}.
#'
#' @param shape Shape parameter of the gamma distribution.
#' @param alpha Shape parameter of the gamma distribution.
#' @param shape1,shape2 non-negative parameters of the Beta distribution.
#' @param k Number of intervals of the discrete gamma distribution.
#' @param inv Proportion of invariable sites.
#' @param g rates of discrete distribution.
#' @param w proportion of rate g.
#' @param site.rate Indicates what type of gamma distribution to use. Options
#' are "gamma" (Yang 1994) and "gamma_quadrature" using Laguerre quadrature
#' approach of Felsenstein (2001)
## or "free_rate" "gamma_phangorn".
#' @param edge.length Total edge length (sum of all edges in a tree).
#' @param discrete logical whether to plot discrete (default) or continuous pdf
#' or cdf.
#' @param cdf logical whether to plot the cumulative distribution function
#' or density / probability function.
#' @param append logical; if TRUE only add to an existing plot.
#' @param xlab a label for the x axis, defaults to a description of x.
#' @param ylab a label for the y axis, defaults to a description of y.
#' @param xlim the x limits of the plot.
#' @param verticals logical; if TRUE, draw vertical lines at steps.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{discrete.gamma} returns a matrix.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml.fit}, \link{stepfun}, \link{pgamma}, \link{pbeta}},
#' @examples
#' discrete.gamma(1, 4)
#'
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow = c(2,1))
#' plot_gamma_plus_inv(shape=2, discrete = FALSE, cdf=FALSE)
#' plot_gamma_plus_inv(shape=2, append = TRUE, cdf=FALSE)
#'
#' plot_gamma_plus_inv(shape=2, discrete = FALSE)
#' plot_gamma_plus_inv(shape=2, append = TRUE)
#' par(old.par)
#'
#' data(Laurasiatherian)
#' fit <- pml_bb(Laurasiatherian, "JC+G(4)", rearrangement = "none")
#' plotRates(fit)
#' @keywords distribution
#' @rdname discrete.gamma
#' @export
discrete.gamma <- function(alpha, k) {
  if (k == 1) return(1)
  quants <- qgamma( (1:(k - 1)) / k, shape = alpha, rate = alpha)
  diff(c(0, pgamma(quants * alpha, alpha + 1), 1)) * k
}


#' @rdname discrete.gamma
#' @importFrom stats pbeta qbeta
#' @export
discrete.beta <- function(shape1, shape2, k) {
  quants <- qbeta( (1:(k - 1)) / k, shape1, shape2)
  diff(c(0, pbeta(quants, shape1 + 1, shape2), 1)) * k * shape1 /
    (shape1 + shape2)
}


#' @rdname discrete.gamma
#' @importFrom stats dgamma qgamma stepfun
#' @importFrom graphics curve
#' @keywords hplot
#' @export
plot_gamma_plus_inv <- function(w=NULL, g=NULL, shape=1, inv=0, k=4, discrete=TRUE, cdf=TRUE,
                                append=FALSE, xlab = "x",
                                ylab=ifelse(cdf, "F(x)", "f(x)"), xlim=NULL,
                                verticals=FALSE, edge.length=NULL,
                                site.rate="gamma", ...){

  l <- max(lengths(list(shape, k, inv)))
  if(l>1){
    shape <- rep_len(shape, l)
    k <- rep_len(k, l)
    inv <- rep_len(inv, l)
    verticals <- rep_len(verticals, l)
  }

  gw <-  function(shape, k, inv, site.rate){
    gw <- rates_n_weights(shape, k, site.rate)
    g <- gw[, 1]
    w <- gw[, 2]
    if (inv > 0){
      w <- c(inv, (1 - inv) * w)
      g <- c(0, g/(1 - inv))
    }
    cbind(w=w, g=g)
  }

  cr <- TRUE
  if(!is.null(w) && !is.null(g)){
    tmp <- cbind(w=w, g=g)
    cr <- FALSE
  }

#  step_GpI <- function(alpha=1, k=4, inv=0, edge.length=NULL,
#                       site.rate="gamma"){
#    rw <- rates_n_weights(alpha, k, site.rate)
#    g <- rw[, 1]
#    w <- rw[, 2]
#    if (inv > 0) w <- c(inv, (1 - inv) * w)
#    cw <- cumsum(w)
#    if (inv > 0) g <- c(0, g/(1 - inv))
#    if(!is.null(edge.length)) g <- g * edge.length
#    stepfun(g, c(0, cw))
#  }

  step_gw <- function(g, w, edge.length=NULL){
    cw <- c(0, cumsum(w))
    if(!is.null(edge.length)) g <- g * edge.length
    stepfun(g, cw)
  }

  g <- mapply(function(shape, k, inv) max(gw(shape, k, inv, site.rate)[,"g"]),
              shape, k, inv) |> max()

  if(is.null(xlim)) xlim <- c(-0.25, 1.25 * g)

  # pgamma_invariant
  cdf_fun <- function(x, shape=1, inv=0){
    if(inv > 0) x <- x * (1-inv)
    y <- (1-inv) * pgamma(x, shape = shape, rate = shape) + inv
    y[x<0] <- 0
    y
  }

  # dgamma_invariant
  density_fun <- function(x, shape, inv){
    if(inv > 0) x <- x * (1-inv)
    (1-inv) * dgamma(x, shape = shape, rate = shape)
  }


#shape=1, inv=0, k=4,
  plot_cdf_discrete <- function(g, w, verticals=FALSE, append=FALSE, ylab=ylab,
                                xlim=xlim, edge.length=edge.length,
                                site.rate="gamma", ...){
#    sf <- step_GpI(shape, inv, k=k, edge.length=edge.length,
#                   site.rate = site.rate)
    sf <- step_gw(g=g, w=w, edge.length=edge.length)
    if(append) plot(sf, verticals=verticals, add=TRUE, xlim=xlim, ...)
    else plot(sf, verticals=verticals, ylab=ylab, xlim=xlim, ...)
  }

#shape=1, inv=0, k=4
  plot_density_discrete <- function(g, w, append=FALSE,  xlab=xlab, ylab=ylab,
                                    xlim=xlim, site.rate="gamma", ...){
#    g_w <- gw(shape, k, inv, site.rate)
#    g <- g_w[, "g"]
#    w <- g_w[, "w"]
    if(!append) plot(g, w, xlim = xlim, ylim=c(0, 1), type="n",
                     xlab=xlab, ylab=ylab, ...)
    segments(g, 0, g, w, ...)
    points(g, w, ...)
  }

  plot_cdf_continuos <- function(shape=1, inv=0, k=4, verticals=FALSE,
                                 append=FALSE, ylab=ylab, xlim=xlim,
                                 site.rate="gamma",...){
    if(inv==0){
      if(!append) plot(function(x)cdf_fun(x, shape, inv), xlim[1], xlim[2],
                       ylim=c(0, 1), xlab=xlab, ylab=ylab, ...)
      else plot(function(x)cdf_fun(x, shape, inv), add=TRUE, ...)
    }
    else{
      g_w <- gw(shape, k, inv, site.rate)
      g <- g_w[, "g"]
      w <- g_w[, "w"]
      if(!append) plot(g, w, xlim = xlim, #c(-.5, 1.25 * max(g)),
                       ylim=c(0, 1),
                       type="n", xlab=xlab, ylab=ylab, ...)
      plot(function(x)cdf_fun(x, shape, inv),  xlim[1], -0.001, add=TRUE, ...)
      plot(function(x)cdf_fun(x, shape, inv),  0, xlim[2], add=TRUE, ...)
      points(0, inv, ...)
      if(verticals) segments(0, 0, 0, inv, ...)
    }
  }

  plot_density_continuous <- function(shape=1, inv=0, k=4, append=FALSE,
                                      xlab=xlab, ylab=ylab, xlim=xlim, ...){
    if(!append) plot(function(x)density_fun(x, shape = shape, inv=inv),
                     xlim[1], xlim[2], xlab=xlab, ylab=ylab, ...)
    else{
      plot(function(x)density_fun(x, shape = shape, inv=inv),
           xlim[1], xlim[2], add=TRUE, ...)
    }
    if(inv>0){
      segments(0, 0, 0, inv, ...)
      points(0, inv, ...)
    }
  }

  i <- 1L
  while(l >= i ){

    if(cdf){
      if(discrete){
        if(cr)tmp <- gw(shape=shape[i], k=k[i], inv=inv[i], site.rate=site.rate)
        plot_cdf_discrete(g=tmp[,"g"], w=tmp[,"w"], verticals = verticals[i],
                          append = append, ylab = ylab, xlim=xlim,
                          edge.length=edge.length, site.rate=site.rate, ...)
      }
      else{
        plot_cdf_continuos(shape[i], inv[i], k[i], verticals = verticals[i],
                           append = append, ylab = ylab, xlim=xlim,
                           site.rate=site.rate,...)
      }
    }
    else{
      if(discrete){
        if(cr)tmp <- gw(shape=shape[i], k=k[i], inv=inv[i], site.rate=site.rate)
        plot_density_discrete(g=tmp[,"g"], w=tmp[,"w"], append=append,
                              xlab=xlab, ylab=ylab, xlim=xlim,
                              site.rate=site.rate, ...)
      }
      else{
        plot_density_continuous(shape[i], inv[i], k[i], append=append,
                                xlab=xlab, ylab=ylab, xlim=xlim, ...)
      }
    }
    append <- TRUE
    i <- i+1L
  }
}



#' @rdname discrete.gamma
#' @importFrom stats ecdf
#' @importFrom graphics rug
## @importFrom statmod gauss.quad.prob
#' @param obj an object of class pml
#' @param main a main title for the plot.
#' @param cdf.color color of the cdf.
#' @param rug logical; if TRUE a \code{\link[graphics]{rug}} is added to the
#' plot.
#' @export
plotRates <- function(obj, cdf.color="blue", main="cdf", rug=FALSE, xlim=NULL, ...){
  pscores <- parsimony(obj$tree, obj$data, site="site")[attr(obj$data, "index")]
  ecdf_pscores <- ecdf(pscores)
  if(is.null(xlim)) xlim <- c(-0.25, 1.1 * max(pscores))
  plot(ecdf_pscores, verticals = TRUE, do.points=FALSE, main=main, xlim=xlim, ...)
  if(rug) rug(jitter(pscores))
  el <- obj$tree$edge.length * obj$rate
  if(obj$site.rate == "free_rate"){
    plot_gamma_plus_inv(w=obj$w, g=obj$g, append=TRUE, xlim=xlim,
                        edge.length=sum(el), verticals=TRUE, col=cdf.color,
                        site.rate=obj$site.rate, ...)
  } else {
      plot_gamma_plus_inv(k=obj$k, shape=obj$shape, inv=obj$inv, append=TRUE,
                      xlim = xlim,
                      edge.length=sum(el), verticals=TRUE, col=cdf.color,
                      site.rate=obj$site.rate, ...)
  }
  invisible(obj)
}


# from selac
LaguerreQuad <- function(shape=1, ncats=4) {
  # Determine rates based on alpha and the number of bins
  # bins roots normalized to 1 of the General Laguerre Quadrature
  # first ncats elements are rates with mean 1
  # second ncats elements are probabilities with sum 1
  roots <- findRoots(shape - 1, ncats)
  weights <- numeric(ncats)
  f <- prod(1 + (shape - 1)/(1:ncats))

  for (i in 1:ncats) {
    weights[i] <- f * roots[i] / ((ncats + 1)^2 *
                                    Laguerre(roots[i], shape - 1, ncats + 1)^2)
  }
  roots <- roots/shape
  return(matrix(c(roots, weights), ncol=2L,
                dimnames = list(NULL, c("rate", "weight"))))
}


findRoots <- function(shape, ncats) {
  # Determine rates based on Gamma's alpha and the number of bins
  # bins roots normalized to 1 of the General Laguerre Polynomial (GLP)
  coeff  <- integer(ncats + 1)
  for (i in 0:ncats) {
    coeff[i + 1] <- (-1)^i * exp(lchoose(ncats + shape, ncats - i) -
                                   lfactorial(i))
  }
  return(sort(Re(polyroot(coeff))))
}


Laguerre <- function(x, shape, degree) {
  y <- 0
  for (i in 0:degree) {
    y <- y + (-1)^i * x^i *
         exp(lchoose(degree + shape, degree - i) - lfactorial(i))
  }
  return(y)
}

# w currently not used! Out?
rates_n_weights <- function(shape, k, site.rate = "gamma", w=NULL, inv=0){
  site.rate <- match.arg(site.rate, c("gamma", "gamma_phangorn",
                                      "gamma_quadrature", "free_rate"))
  if(site.rate == "gamma_quadrature")
    return(LaguerreQuad(shape=shape, k))
  if(k==1){
    g <- 1
    w <- 1
  }
    #rates.and.weights <- matrix(c(1,1), ncol=2L,
          #                        dimnames = list(NULL, c("rate", "weight")))
  else{
    if(site.rate == "gamma"){
      g <- discrete.gamma(shape, k=k)
      w <- rep(1 / k, k)
#      rates.and.weights <- matrix( c(g, w), ncol=2L,
#                          dimnames = list(NULL, c("rate", "weight")))
    }
#    if(site.rate == "gamma_phangorn"){
#      rates.and.weights <- discrete.gamma.2(alpha=shape, k=k)
#    }

    if(site.rate == "free_rate"){
      g <- discrete.gamma(1, k=k) # rep(1, k)
      w <- rep(1 / k, k)
#      rates.and.weights <- matrix( c(g, w), ncol=2L,
#                                   dimnames = list(NULL, c("rate", "weight")))
    }
  }
  if (inv > 0){
    w <- (1 - inv) * w
    g <- g / (1 - inv)
  }
  rates.and.weights <- matrix( c(g, w), ncol=2L,
                               dimnames = list(NULL, c("rate", "weight")))
  rates.and.weights
}


discrete.gamma.2 <- function(alpha, k){
  if (k == 1) return(list(w=1, g=1))
  bin <- c(rep(0, k), 1)
  quants <- rep(0, k+1)
  quants[k+1] <- 1
  for(i in 2:k){
    old_bin <- bin[i-1]
    fun <- function(x, k, alpha, old_bin){
      quants <- qgamma(c(old_bin, x), shape = alpha, rate = alpha)
      tmp <- diff(pgamma(quants * alpha, alpha + 1)) * (1 / (x-old_bin))
      abs( (x - old_bin) * tmp - 1/k )
    }
    res <- optimize(fun, k=k, alpha=alpha, old_bin=old_bin,
                    interval=c(old_bin, 1), tol = .Machine$double.eps^0.5)
    bin[i] <- res$minimum
  }
  w <- diff(bin)
  quants <- qgamma( bin[seq_len(k)], shape = alpha, rate = alpha)
  g <- diff(c(pgamma(quants * alpha, alpha + 1), 1)) * (1/w)
  matrix(c(g, w), ncol=2L, dimnames = list(NULL, c("rate", "weight")))
}


# free_rate <- function(k) nur optimisieren

#' @srrstats {G2.3, G2.3a} in lines: 285
