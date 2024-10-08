% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gap_as_state.R
\name{gap_as_state}
\alias{gap_as_state}
\alias{gap_as_ambiguous}
\alias{has_gap_state}
\title{Treat gaps as a state}
\usage{
gap_as_state(obj, gap = "-", ambiguous = "?")

gap_as_ambiguous(obj, gap = "-")

has_gap_state(obj)
}
\arguments{
\item{obj}{An object of class phyDat.}

\item{gap}{a character which codes for the gaps (default is "-").}

\item{ambiguous}{a character which codes for the ambiguous state}
}
\value{
The functions return an object of class \code{phyDat}.
}
\description{
The function \code{gap_as_state} changes the contrast of an phyDat object to
treat as its own state. Internally \code{phyDat} are stored similar to a
\code{factor} objects and only the contrast matrix and some attributes
change.
}
\examples{
data(Laurasiatherian)
tmp <- gap_as_state(Laurasiatherian)
contr <- attr(tmp, "contrast")
rownames(contr) <- attr(tmp, "allLevels")
contr
}
\seealso{
\code{\link{phyDat}}, \code{\link{latag2n.phyDat}},
\code{\link[ape]{latag2n}}, \code{\link{ancestral.pml}},
\code{\link{gap_as_state}}
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\keyword{cluster}
