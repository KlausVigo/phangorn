\name{modelTest}
\alias{modelTest}
\alias{AICc}
\title{
ModelTest
}
\description{
Comparison of different nucleotide or amino acid substitution models 
}
\usage{
modelTest(object, tree=NULL, model = c("JC", "F81", "K80", "HKY", "SYM", "GTR"), 
    G = TRUE, I = TRUE, FREQ = FALSE, k = 4, control = pml.control(epsilon = 1e-08, 
    maxit = 10, trace = 1), multicore = FALSE, mc.cores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class phyDat or pml}
  \item{tree}{a phylogenetic tree.}
  \item{model}{a vector containing the substitution models to compare with each other
  or "all" to test all available models}
  \item{G}{logical, TRUE (default) if (discrete) Gamma model should be tested}
  \item{I}{logical, TRUE (default) if invariant sites should be tested}
  \item{FREQ}{logical, FALSE (default) if TRUE amino acid frequencies will be estimated.}
  \item{k}{number of rate classes}
  \item{control}{A list of parameters for controlling the fitting process.}
  \item{multicore}{logical, whether models should estimated in parallel.}
  \item{mc.cores}{The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one, and parallelization requires at least two cores.}
}
\details{
\code{modelTest} estimates all the specified models for a given tree and data. 
When the mclapply is available, the computations are done in parallel. \code{modelTest} 
runs each model in one thread. 
This is may not work within a GUI interface and will not work under Windows.  
}
\value{
A data.frame containing the log-likelihood, number of estimated parameters, AIC, AICc and BIC all tested models. 
The data.frame has an attributes "env" which is an environment which contains all the trees, the data and the calls to allow get the estimated models, e.g. as a starting point for further analysis (see example).
}
\references{
Burnham, K. P. and Anderson, D. R (2002)  \emph{Model selection and multimodel inference: a practical information-theoretic approach}. 2nd ed. Springer, New York

Posada, D. and Crandall, K.A. (1998) MODELTEST: testing the model of DNA substitution. \emph{Bioinformatics} \bold{14(9)}: 817-818

Posada, D. (2008) jModelTest: Phylogenetic Model Averaging. \emph{Molecular Biology and Evolution} \bold{25}: 1253-1256

Darriba D., Taboada G.L., Doallo R and Posada D. (2011) ProtTest 3: fast selection of best-fit models of protein evolution. . \emph{Bioinformatics} \bold{27}: 1164-1165
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pml}}, \code{\link{anova}}, \code{\link[stats]{AIC}}
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

# test all available amino acid models
(mTAA_all <- modelTest(chloroplast, model="all", multicore=TRUE, mc.cores=2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

