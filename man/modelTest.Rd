\name{modelTest}
\alias{modelTest}
\title{
ModelTest
}
\description{
Comparison of different substition models 
}
\usage{
modelTest(object, tree=NULL, model = c("JC", "F81", "K80", "HKY", "SYM", "GTR"), 
    G = TRUE, I = TRUE, k = 4, control = pml.control(epsilon = 1e-08, maxit = 3, 
    trace = 1), multicore = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class phyDat or pml}
  \item{tree}{a phylogenetic tree.}
  \item{model}{a vector containing the substitution models to compare with each other}
  \item{G}{logical, TRUE (default) if (discrete) Gamma modelshould be tested}
  \item{I}{logical, TRUE (default) if invariant sites should be tested}
  \item{k}{number of rate classes}
  \item{control}{A list of parameters for controlling the fitting process.}
  \item{multicore}{logical, whether models should estimated in parallel.}
}
\details{
\code{modelTest} estimates all the specified models for a given tree and data. 
When the multicore package is available, the computations are done in parallel.
This is only possible without GUI interface and under linux.  
Only nucleotide models are tested so far.  
}
\value{
A data.frame containing the log-likelihood, AIC and BIC all tested models. 
The data.frame has an attributes "env" which is an environment which contains all the trees, the data and the calls to allow get the estimated models, e.g. as a starting point for further analysis (see example).
}
\references{
Posada, D. and Crandall, K.A. (1998) MODELTEST: testing the model of DNA substitution. \emph{Bioinformatics} \bold{14(9)}: 817-818

Posada, D. (2008) jModelTest: Phylogenetic Model Averaging. \emph{Molecular Biology and Evolution} \bold{25}: 1253-1256

Darriba D., Taboada G.L., Doallo R and Posada D. (2011) ProtTest 3: fast selection of best-fit models of protein evolution. . \emph{Bioinformatics} \bold{27}: 1164-1165
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pml}}, \code{\link{anova}}
}
\examples{
\dontrun{    
example(NJ)
(mT <- modelTest(Laurasiatherian, tree))

# some R magic
env = attr(mT, "env")
ls(env=env)
(F81 <- get("F81+G", env)) # a call  
eval(F81, env=env)

data(chloroplast)
(mTAA <- modelTest(chloroplast, model=c("JTT", "WAG", "LG")))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

