#' Discrete Gamma function
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
#' @param k Number of intervals of the discrete gamma distribution.
#' @param inv Proportion of invariable sites.
#' @param discrete logical wether to plot discrete (default) or continous pdf or
#' cdf.
#' @param cdf logical wether to plot the cummulative distribution function
#' or density / probability function.
#' @param append logical; if TRUE only add to an existing plot.
#' @param xlab a label for the x axis, defaults to a description of x.
#' @param ylab a label for the y axis, defaults to a description of y.
#' @param xlim the x limits of the plot.
#' @param verticals ogical; if TRUE, draw vertical lines at steps.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{discrete.gamma} returns a matrix.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml.fit}, \link{stepfun}}
#' @examples
#' discrete.gamma(1, 4)
#'
#' par(mfrow = c(2,1))
#'
#' plot_gamma_plus_inv(shape=2, discrete = FALSE, cdf=FALSE)
#' plot_gamma_plus_inv(shape=2, append = TRUE, cdf=FALSE)
#'
#' plot_gamma_plus_inv(shape=2, discrete = FALSE)
#' plot_gamma_plus_inv(shape=2, append = TRUE)
#'
#' par(mfrow = c(1,1))
#'
#' @keywords cluster
#'
#' @rdname discrete.gamma
#' @export
discrete.gamma <- function(alpha, k) {
  if (k == 1) return(1)
  quants <- qgamma( (1:(k - 1)) / k, shape = alpha, rate = alpha)
  diff(c(0, pgamma(quants * alpha, alpha + 1), 1)) * k
}


discrete.beta <- function(shape1, shape2, k) {
  qbeta( ( (0:(k - 1)) + .5) / k, shape1, shape2)
}

#' @rdname discrete.gamma
#' @importFrom stats dgamma qgamma stepfun
#' @importFrom graphics curve
#' @export
plot_gamma_plus_inv <- function(shape=1, inv=0, k=4, discrete=TRUE, cdf=TRUE,
                                append=FALSE, xlab = "x",
                                ylab=ifelse(cdf, "F(x)", "f(x)"), xlim=NULL,
                                verticals=FALSE, ...){

  l <- max(lengths(list(shape, k, inv)))
  if(l>1){
    shape <- rep_len(shape, l)
    k <- rep_len(k, l)
    inv <- rep_len(inv, l)
    verticals <- rep_len(verticals, l)
  }

  gw <-  function(shape, k, inv){
    w <- rep(1/k, k)
    g <- discrete.gamma(shape, k)
    if (inv > 0){
      w <- c(inv, (1 - inv) * w)
      g <- c(0, g/(1 - inv))
    }
    cbind(w=w, g=g)
  }

  g <- mapply(function(shape, k, inv) max(gw(shape, k, inv)[,"g"]),
              shape, k, inv) %>% max()

  if(is.null(xlim)) xlim <- c(-0.25, 1.25 * g)

  cdf_fun <- function(x, shape=1, inv=0){
    if(inv > 0) x <- x * (1-inv)
    y <- (1-inv) * pgamma(x, shape = shape, rate = shape) + inv
    y[x<0] <- 0
    y
  }

  density_fun <- function(x, shape, inv){
    if(inv > 0) x <- x * (1-inv)
    (1-inv) * dgamma(x, shape = shape, rate = shape)
  }

  step_GpI <- function(alpha=1, k=4, inv=0){
    w <- rep(1/k, k)
    if (inv > 0) w <- c(inv, (1 - inv) * w)
    cw <- cumsum(w)
    g <- discrete.gamma(alpha, k)
    if (inv > 0) g <- c(0, g/(1 - inv))
    stepfun(g, c(0, cw))
  }

  plot_cdf_discrete <- function(shape=1, inv=0, k=4, verticals=FALSE,
                                append=FALSE, ylab=ylab, xlim=xlim, ...){
    sf <- step_GpI(shape, inv, k=k)
    if(append) plot(sf, verticals=verticals, add=TRUE, xlim=xlim, ...)
    else plot(sf, verticals=verticals, ylab=ylab, xlim=xlim, ...)
  }

  plot_density_discrete <- function(shape=1, inv=0, k=4, append=FALSE,
                                    xlab=xlab, ylab=ylab, xlim=xlim,...){
    g_w <- gw(shape, k, inv)
    g <- g_w[, "g"]
    w <- g_w[, "w"]
    if(!append) plot(g, w, xlim = xlim, ylim=c(0, 1), type="n",
                     xlab=xlab, ylab=ylab, ...)
    segments(g, 0, g, w, ...)
    points(g, w, ...)
  }

  plot_cdf_continuos <- function(shape=1, inv=0, k=4, verticals=FALSE,
                                 append=FALSE, ylab=ylab, xlim=xlim, ...){
    if(inv==0){
      if(!append) plot(function(x)cdf_fun(x, shape, inv), xlim[1], xlim[2],
                       ylim=c(0, 1), xlab=xlab, ylab=ylab, ...)
      else plot(function(x)cdf_fun(x, shape, inv), add=TRUE, ...)
    }
    else{
      g_w <- gw(shape, k, inv)
      g <- g_w[, "g"]
      w <- g_w[, "w"]
      if(!append) plot(g, w, xlim = xlim, #c(-.5, 1.25 * max(g)),
                       ylim=c(0, 1),
                       type="n", xlab=xlab, ylab=ylab, ...)
      plot(function(x)cdf_fun(x, shape, inv),  xlim[1], -.001, add=TRUE, ...)
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
        plot_cdf_discrete(shape[i], inv[i], k[i], verticals = verticals[i],
                          append = append, ylab = ylab, xlim=xlim, ...)
      }
      else{
        plot_cdf_continuos(shape[i], inv[i], k[i], verticals = verticals[i],
                           append = append, ylab = ylab, xlim=xlim, ...)
      }
    }
    else{
      if(discrete){

        plot_density_discrete(shape[i], inv[i], k[i], append=append,
                              xlab=xlab, ylab=ylab, xlim=xlim, ...)
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
