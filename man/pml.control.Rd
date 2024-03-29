% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pml_control.R
\name{pml.control}
\alias{pml.control}
\alias{ratchet.control}
\title{Auxiliary for Controlling Fitting}
\usage{
pml.control(epsilon = 1e-08, maxit = 10, trace = 1, tau = 1e-08,
  statefreq = "empirical")

ratchet.control(iter = 20L, maxit = 200L, minit = 100L, prop = 1/2,
  rell = TRUE, bs = 1000L)
}
\arguments{
\item{epsilon}{Stop criterion for optimization (see details).}

\item{maxit}{Maximum number of iterations (see details).}

\item{trace}{Show output during optimization (see details).}

\item{tau}{minimal edge length.}

\item{statefreq}{take "empirical" or "estimate" state frequencies.}

\item{iter}{Number of iterations to stop if there is no change.}

\item{minit}{Minimum number of iterations.}

\item{prop}{Only used if \code{rearrangement=stochastic}. How many NNI moves
should be added to the tree in proportion of the number of taxa.´}

\item{rell}{logical, if TRUE approximate bootstraping similar Minh et al.
(2013) is performed.}

\item{bs}{number of approximate bootstrap samples.}
}
\value{
A list with components named as the arguments for controlling the
fitting process.
}
\description{
Auxiliary functions for providing \code{\link{optim.pml}, \link{pml_bb}}
fitting. Use it to construct a \code{control} or \code{ratchet.par} argument.
}
\details{
\code{pml.control} controls the fitting process. \code{epsilon} and
\code{maxit} are only defined for the most outer loop, this affects
\code{pmlCluster}, \code{pmlPart} and \code{pmlMix}.

\code{epsilon} is not an absolute difference between, but instead is
defined as (logLik(k)-logLik(k+1))/logLik(k+1). This seems to be a good
compromise and to work reasonably well for small and large trees or
alignments.

If \code{trace} is set to zero than no out put is shown, if functions are
called internally than the trace is decreased by one, so a higher of trace
produces more feedback. It can be useful to figure out how long an run will
take and for debugging.

\code{statefreq} controls if base/state frequencies are optimized or
empirical estimates are taken, when this applies. For some nucleotide models
(e.g. JC, SYM) equal base frequencies and for amino acid models precomputed
state frequencies are used, if not '+F' is specified.

\code{tau} might be exactly zero if duplicated sequences in the alignment are
observed. In this case the analysis is performed only on unique sequences and
duplicated taxa are added to the tree with zero edge length. This may lead to
multifurcations if there are three or more identical sequences. After
optimization it is good practice to prune away edges of length \code{tau}
using \code{di2multi}. See also Janzen et al. (2021).
}
\examples{
pml.control()
pml.control(maxit=25)
}
\references{
Minh, B. Q., Nguyen, M. A. T., & von Haeseler, A. (2013).
Ultrafast approximation for phylogenetic bootstrap. \emph{Molecular biology
and evolution}, \bold{30(5)}, 1188-1195.

Janzen, T., Bokma, F.,Etienne,  R. S. (2021) Nucleotide Substitutions during
Speciation may Explain Substitution Rate Variation,
\emph{Systematic Biology}, \bold{71(5)}, 1244–1254.
}
\seealso{
\code{\link{pml_bb}}, \code{\link{optim.pml}}
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
