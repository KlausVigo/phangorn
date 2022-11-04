#' Auxiliary for Controlling Fitting
#'
#' Auxiliary functions for \code{\link{optim.pml}} fitting. Use it to construct
#' a \code{control} or \code{ratchet.par} argument.
#'
#' \code{pml.control} controls the fitting process. \code{epsilon} and
#' \code{maxit} are only defined for the most outer loop, this affects
#' \code{pmlCluster}, \code{pmlPart} and \code{pmlMix}.  \code{epsilon} is
#' defined as (logLik(k)-logLik(k+1))/logLik(k+1), this seems to be a good
#' heuristics which works reasonably for small and large trees or alignments.
#' If \code{trace} is set to zero than no out put is shown, if functions are
#' called internally than the trace is decreased by one, so a higher of trace
#' produces more feedback.
#'
### @param control A list of parameters for controlling the fitting process.
#' @param epsilon Stop criterion for optimization (see details).
#' @param maxit Maximum number of iterations (see details).
#' @param trace Show output during optimization (see details).
#' @param tau minimal edge length.
#' @param minit Minimum number of iterations.
#' @param iter Number of iterations to stop if there is no change.
#' @param prop Only used if \code{rearrangement=stochstic}. How many NNI moves
#' should be added to the tree in proportion of the number of taxa.´
#' @param rell logical, if TRUE approximate bootstraping similar Minh et al.
#' (2013) is performed.
#' @param bs number of approximate bootstrap samples.
#' @return A list with components named as the arguments for controlling the
#' fitting process.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{optim.pml}}
#' @references Minh, B. Q., Nguyen, M. A. T., & von Haeseler, A. (2013).
#' Ultrafast approximation for phylogenetic bootstrap. \emph{Molecular biology
#' and evolution}, \bold{30(5)}, 1188-1195.
#' @examples
#' pml.control()
#' pml.control(maxit=25)
#' @export
pml.control <- function(epsilon = 1e-08, maxit = 10, trace = 1, tau = 1e-8) {
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(tau) || tau <= 0)
    stop("tau must be > 0")
  list(epsilon = epsilon, maxit = maxit, trace = trace, tau = tau)
}

#' @rdname pml.control
#' @export
ratchet.control <- function(iter = 20L, maxit = 200L, minit = 50L, prop = 1/2,
                            rell = TRUE, bs=1000L){
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(minit) || minit <= 0)
    stop("minimum number of iterations must be > 0")
  if (!is.numeric(iter) || iter <= 0)
    stop("number of iterations must be > 0")
  if (!is.numeric(iter) || iter <= 0)
    stop("proportion of rearrangenemts must be > 0")
  list(iter=iter, maxit=maxit, minit = minit, prop = prop, rell = rell, bs=bs)
}
