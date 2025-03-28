% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pml.R
\name{plot.pml}
\alias{plot.pml}
\title{Plot phylogeny of a pml object}
\usage{
\method{plot}{pml}(x, type = "phylogram", direction = "rightwards", ...,
  adj = NULL, digits = 2, method = "FBP")
}
\arguments{
\item{x}{an object of class \code{pml}.}

\item{type}{a character string specifying the type of phylogeny to be drawn;
it must be one of "phylogram" (the default), "cladogram", "fan", "unrooted",
"radial", "tidy", or any unambiguous abbreviation of these.}

\item{direction}{a character string specifying the direction of the tree.
Four values are possible: "rightwards" (the default), "leftwards", "upwards",
and "downwards".}

\item{\dots}{further parameters to be passed to \code{plot.phylo}.}

\item{adj}{one or two numeric values specifying the horizontal and vertical
justification of the text or symbols of the support values.}

\item{digits}{integer indicating the number of decimal places.}

\item{method}{either "FBP" the classical bootstrap (default), "TBE"
(transfer bootstrap) or "MCC" for assigning clade credibilities.}
}
\value{
\code{plot.pml} returns the \code{pml} object x.
}
\description{
\code{plot.pml} is a wrapper around \code{plot.phylo} with different default
values for unrooted, ultrametric and tip dated phylogenies.
}
\examples{
fdir <- system.file("extdata/trees", package = "phangorn")
tmp <- read.csv(file.path(fdir,"H3N2_NA_20.csv"))
H3N2 <- read.phyDat(file.path(fdir,"H3N2_NA_20.fasta"), format="fasta")
dates <- setNames(tmp$numdate_given, tmp$name)

fit_td <- pml_bb(H3N2, model="JC", method="tipdated", tip.dates=dates,
                 rearrangement="none", control = pml.control(trace = 0))
plot(fit_td, show.tip.label = FALSE)
# Same as:
# root_time <- max(dates) - max(node.depth.edgelength(fit_td$tree))
# plot(fit_td$tree, show.tip.label = FALSE)
# axisPhylo(root.time = root_time, backward = FALSE)
plot(fit_td, show.tip.label = FALSE, direction="up")

fit_unrooted <- pml_bb(H3N2, model="JC", rearrangement="none",
                       control = pml.control(trace = 0))
plot(fit_unrooted, cex=.5)

}
\seealso{
\code{\link[ape]{plot.phylo}}, \code{\link[ape]{axisPhylo}},
\code{\link[ape]{add.scale.bar}}
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\keyword{hplot}
