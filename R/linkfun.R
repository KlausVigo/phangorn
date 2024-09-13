#' Family objects for evolutionary models
#'
#' Family objects provide a convenient way to specify the details of the models
#' used by functions such as \code{glm}. \code{binomial_mk} extends the
#' \code{binomial} family with the Mk model. The 2 state model is also known
#' Neyman and 4 state model as the Jukes Cantor model.
#' See the documentation for \code{\link{family}} for more details and
#' \code{\link{glm}} for the details on how such model fitting takes place.
#'
#' The link function for the Jukes Cantor 4 state model is
#' \deqn{g(\mu) = -0.75 \cdot log(1-(4/3)\cdot\mu)}{-0.75 * log(1 - (4/3) * mu)}
#'
#' @seealso \code{\link{family}}, \code{\link{binomial}}, \code{\link{glm}}
#' @param k number of states.
#' @importFrom stats binomial dbinom
#' @examples
#' plot(function(x) binomial_mk()$linkinv(x), 0, 1.5, , ylim=c(0, 1), asp=1,
#'      main = "inverse link function (JC69)", ylab="f(x)")
#' abline(0, 1)
#' abline(h=.75)
#' plot(function(x) binomial_mk()$linkfun(x), 0, .75,
#'      main = "link function (JC69)", ylab="f(x)")
#'
#' data(yeast)
#' dm <- dist.hamming(yeast, FALSE)
#' y <- cbind(dm, 127026 - dm)
#' tree <- nj(dm)
#' X <- designTree(tree)
#' glm(y ~ X -1, binomial_mk())
#' @noRd
binomial_mk <- function (k=4)
{
  stats <- make_link_mk(k)
  variance <- function(mu) mu * (1 - mu)

  tmp <- binomial()
  aic <- function(y, n, mu, wt, dev) {
    m <- if (any(n > 1))
      n
    else wt
    -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * y), round(m), mu,
                                               log = TRUE))
  }
  #  initialize fuer verschiedene Modelle veraendern !!!
  initialize <- expression({
    if (NCOL(y) == 1) {
      if (is.factor(y)) y <- y != levels(y)[1]
      n <- rep.int(1, nobs)
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
      mustart <- (weights * y + 0.5)/(2*(weights + 1))
      m <- weights * y
      if (any(abs(m - round(m)) > 0.001)) warning("non-integer #successes in a binomial glm!")
    } else if (NCOL(y) == 2) {
      if (any(abs(y - round(y)) > 0.001)) warning("non-integer counts in a binomial glm!")
      n <- y[, 1] + y[, 2]
      y <- ifelse(n == 0, 0, y[, 1]/n)
      weights <- weights * n
      mustart <- (n * y + 0.5)/(2*(n + 1))
    } else stop("for the binomial family, y must be a vector of 0 and 1's\n",
                "or a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
  })
  structure(list(family = "binomial", link = "Mk", linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = variance,
                 dev.resids = tmp$dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = tmp$initialize,
                 validmu = stats$validmu, valideta = stats$valideta,
                 simulate = tmp$simfun, dispersion = 1), class = "family")
}


make_link_mk <- function (n=4L){
  f <- substitute(ifelse(mu>0, -(k-1)/k * log(1-(k/(k-1))*mu),0), list(k=n))
  linkfun <- as.function(list(f))
  formals(linkfun) <- alist(mu=)
  f <- substitute(ifelse(eta>0, (k-1)/k*(1-exp(-(k/(k-1))*eta)),0), list(k=n))
  linkinv <-  as.function(list(f))
  formals(linkinv) <- alist(eta=)
  f <- substitute(exp(-(k/(k-1))*eta), list(k=n))
  mu.eta <-  as.function(list(f))
  formals(mu.eta) <- alist(eta=)
  #  validmu <- function(mu) all(mu > 0) && all(mu < .5)
  f <- substitute(all(mu > 0) && all(mu < (k-1)/k), list(k=n))
  validmu <- as.function(list(f))
  formals(validmu) <- alist(mu=)
  valideta <- function(eta) all(eta > 0)
  list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
       valideta = valideta)
}

