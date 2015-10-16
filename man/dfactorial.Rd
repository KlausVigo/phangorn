\name{dfactorial}
\alias{dfactorial}
\alias{ldfactorial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Arithmetic Operators }
\description{
  double factorial function
}
\usage{
dfactorial(x)
ldfactorial(x)
}
\arguments{
  \item{x}{ a numeric scalar or vector }
}

\value{
\code{dfactorial(x)} returns the double factorial, that is \eqn{x\!\! = 1 *  3 * 5 *  \ldots * x } 
and \code{ldfactorial(x)} is the natural logarithm of it. 
}

\author{ Klaus Schliep \email{klaus.schliep@gmail.com}  }

\seealso{ \code{\link[base:Special]{factorial}}, \code{\link[ape]{howmanytrees}} }
\examples{
dfactorial(1:10)
}
\keyword{ classif }% at least one, from doc/KEYWORDS
