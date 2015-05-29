\name{splitsNetwork}

\alias{splitsNetwork}
\title{Phylogenetic Network}
\description{
\code{splitsNetwork} estimates weights for a splits graph from a distance matrix.
}
\usage{
splitsNetwork(dm, splits=NULL, gamma=.1, lambda=1e-6, weight=NULL)
}
\arguments{
  \item{dm}{A distance matrix.}
  \item{splits}{a splits object, containing all splits to consider, 
   otherwise all possible splits are used}
  \item{gamma}{penalty value for the L1 constraint.}
  \item{lambda}{penalty value for the L2 constraint.}
  \item{weight}{a vector of weights.}
}

\details{
\code{splitsNetwork} fits non-negative least-squares phylogenetic networks using L1 (LASSO), L2(ridge regression) constraints.  
The function minimizes the penalized least squares
\deqn{\beta = min \sum(dm - X\beta)^2 + \lambda \|\beta \|^2_2 }{%
beta = sum(dm - X*beta)^2 + lambda |beta|^2_2 }
with respect to \deqn{\|\beta \|_1 <= \gamma, \beta >= 0}{%
|beta|_1 = gamma, beta >= 0}
where X is a design matrix constructed with \code{designSplits}. 
External edges are fitted without L1 or L2 constraints. 
}
\value{
\code{splitsNetwork} returns a splits object with a matrix added. 
The first column contains the indices of the splits, the second 
column an unconstrained fit without penalty terms and the third 
column the constrained fit.
}
\references{
Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression" (with discussion) Annals of Statistics

K. P. Schliep (2009). Some Applications of statistical phylogenetics (PhD Thesis)
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{\code{\link[phangorn]{distanceHadamard}}, \code{\link[phangorn]{designTree}}
\code{\link[phangorn]{consensusNet}}, \code{\link[phangorn]{plot.networx}}  
}

\examples{
data(yeast)
dm = dist.ml(yeast)
fit = splitsNetwork(dm)
net = as.networx(fit)
plot(net, "2D")
write.nexus.splits(fit)
}
\keyword{ cluster }% at least one, from doc/KEYWORDS


