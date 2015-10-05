\name{cophenetic.networx}
\alias{cophenetic.networx}
\alias{cophenetic.splits}
\title{Pairwise Distances from a Phylogenetic Network}
\usage{
\method{cophenetic}{networx}(x)
\method{cophenetic}{splits}(x)
}
\arguments{
  \item{x}{an object of class \code{"networx"}.}
}
\description{
  \code{cophenetic.networx} computes the pairwise distances between the 
  pairs of tips from a phylogenetic network using its branch lengths.
}
\value{
  an object of class dist names set to the names of the
  tips (as given by the element \code{tip.label} of the argument \code{x}).
}
\author{Klaus Schliep}
\seealso{
  \code{\link[stats]{cophenetic}} for the generic function, 
  \code{neighborNet} to construct a network from a distance matrix
}
\keyword{manip}
