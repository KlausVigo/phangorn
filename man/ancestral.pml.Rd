% Please edit documentation in R/ancestral_pml.R, R/parsimony.R
\name{ancestral.pml}
\alias{ancestral.pars}
\alias{ancestral.pml}
\alias{pace}
\alias{plotAnc}
\title{Ancestral character reconstruction.}
\usage{
ancestral.pml(object, type = c("marginal", "ml", "bayes"))

ancestral.pars(tree, data, type = c("MPR", "ACCTRAN"), cost = NULL)

plotAnc(tree, data, i = 1, col = NULL, cex.pie = par("cex"), pos = "bottomright", 
    ...)
}
\arguments{
\item{object}{an object of class pml}

\item{type}{method used to assign characters to internal nodes, see details.}

\item{tree}{a tree, i.e. an object of class pml}

\item{data}{an object of class phyDat}

\item{cost}{A cost matrix for the transitions between two states.}

\item{i}{plots the i-th character of the \code{data}.}

\item{col}{a vector containing the colors for all possible states.}

\item{cex.pie}{a numeric defining the size of the pie graphs}

\item{pos}{a character string defining the position of the legend}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
%A matrix containing the the estimates character states.  An object
of class "phyDat", containing the ancestral states of all nodes.
}
\description{
Marginal reconstruction of the ancestral character states.
}
\details{
The argument "type" defines the criterion to assign the internal nodes. For
\code{ancestral.pml} so far "ml" and (empirical) "bayes" and for
\code{ancestral.pars} "MPR" and "ACCTRAN" are possible.

With parsimony reconstruction one has to keep in mind that there will be
often no unique solution.

For further details see vignette("Ancestral").
}
\examples{

example(NJ)
fit = pml(tree, Laurasiatherian)
anc.ml = ancestral.pml(fit, type = "ml")
anc.p = ancestral.pars(tree, Laurasiatherian)
\dontrun{
require(seqLogo)
seqLogo(t(subset(anc.ml, 48, 1:20)[[1]]), ic.scale = FALSE)
seqLogo(t(subset(anc.p, 48, 1:20)[[1]]), ic.scale = FALSE)
}
plotAnc(tree, anc.ml, 1)

}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\references{
Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer
Associates, Sunderland.

Swofford, D.L., Maddison, W.P. (1987) Reconstructing ancestral character
states under Wagner parsimony. \emph{Math. Biosci.} \bold{87}: 199--229

Yang, Z. (2006). \emph{Computational Molecular evolution}. Oxford University
Press, Oxford.
}
\seealso{
\code{\link{pml}}, \code{\link{parsimony}}, \code{\link[ape]{ace}},
\code{\link[ape]{root}}
}
\keyword{~kwd1}
\keyword{~kwd2}
