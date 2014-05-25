\name{pml}
\alias{pml}
\alias{optim.pml}
\alias{pml.control}
\title{Likelihood of a tree.}
\description{
\code{pml} computes the likelihood of a phylogenetic tree 
given a sequence alignment and a model. \code{optim.pml} optimizes the 
different model parameters.
}
\usage{
pml(tree, data, bf=NULL, Q=NULL, inv=0, k=1, shape=1, rate=1, model="", ...)     
optim.pml(object, optNni=FALSE, optBf=FALSE, optQ=FALSE, optInv=FALSE, optGamma=FALSE,
    optEdge=TRUE, optRate=FALSE, optRooted=FALSE, control = pml.control(epsilon=1e-08, 
    maxit=10, trace=1), model = NULL, subs = NULL, ...)  
pml.control(epsilon = 1e-08, maxit = 10, trace = 1)
}
\arguments{
  \item{tree}{A phylogenetic \code{tree}, object of class \code{phylo}. }
  \item{data}{An alignment, object of class \code{phyDat}.}
  \item{bf}{Base frequencies.}
  \item{Q}{A vector containing the lower triangular part of the rate matrix.}
  \item{inv}{Proportion of invariable sites.}
  \item{k}{Number of intervals of the discrete gamma distribution.}
  \item{shape}{Shape parameter of the gamma distribution.}
  \item{rate}{Rate.}
  \item{model}{allows to choose an amino acid models or nucleotide model, see details.}
  \item{object}{An object of class \code{pml}.}
  \item{optNni}{Logical value indicating whether toplogy gets optimized (NNI).}
  \item{optBf}{Logical value indicating whether base frequencies gets optimized.}
  \item{optQ}{Logical value indicating whether rate matrix gets optimized.}
  \item{optInv}{Logical value indicating whether proportion of variable size gets optimized.}
  \item{optGamma}{Logical value indicating whether gamma rate parameter gets optimized.}
  \item{optEdge}{Logical value indicating the edge lengths gets optimized.}
  \item{optRate}{Logical value indicating the overall rate gets optimized.}
  \item{optRooted}{Logical value indicating if the edge lengths of a rooted tree get optimized.}
  \item{control}{A list of parameters for controlling the fitting process.}
  \item{subs}{A (integer) vector same length as Q to specify the optimization of Q}
  \item{\dots}{Further arguments passed to or from other methods.}
  \item{epsilon}{Stop criterion for optimisation (see details).} 
  \item{maxit}{Maximum number of iterations (see details).}  
  \item{trace}{Show output during otimization (see details).}    
}
\details{
The topology search uses a nearest neighbor interchange (NNI) 
and the implementation is similar to phyML. 
The option model in pml is only used for amino acid models. 
The option model defines the nucleotide model which is getting optmised, 
all models which are included in modeltest can be chosen. Setting this option 
(e.g. "K81" or "GTR") overrules options optBf and optQ.  
Here is a overview how to estimate different phylogenetic models 
with \code{pml}:  
\tabular{lll}{
model \tab optBf \tab optQ \cr
Jukes-Cantor \tab FALSE \tab FALSE \cr
F81 \tab TRUE \tab FALSE \cr
symmetric \tab FALSE \tab TRUE \cr
GTR \tab TRUE \tab TRUE
}
Via model in optim.pml the following nucleotide models can be specified:  
JC, F81, K80, HKY, TrNe, TrN, TPM1, K81, TPM1u, TPM2, TPM2u, TPM3, TPM3u, 
TIM1e, TIM1, TIM2e, TIM2, TIM3e, TIM3, TVMe, TVM, SYM and GTR. 
These models are specified as in Posada (2008).

So far 17 amino acid models are supported ("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb", "FLU", "Blossum62", "Dayhoff_DCMut" and "JTT_DCMut") and additionally rate matrices and amino acid frequences can be supplied. 

If the option 'optRooted' is set to TRUE than the edge lengths of rooted tree are optimized.
The tree has to be rooted and by now ultrametric! Optimising rooted trees is generally much slower.
  
\code{pml.control} controls the fitting process. \code{epsilon} and \code{maxit} are only defined 
for the most outer loop, this affects \code{pmlCluster}, \code{pmlPart} and \code{pmlMix}. 
\code{epsilon} is defined as (logLik(k)-logLik(k+1))/logLik(k+1), this seems to be a good 
heuristics which works reasonalby for small and large trees or alignments. 
If \code{trace} is set to zero than no out put is shown, if functions are called internally 
than the trace is decreased by one, so a higher of trace produces more feedback. 
}
\value{
Returns a list of class \code{ll.phylo}
  \item{logLik}{Log likelihood of the tree.}
  \item{siteLik}{Site log likelihoods.}
  \item{root}{likelihood in the root node.}
  \item{weight}{Weight of the site patterns.}
}
\references{ 
Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maxumum
likelihood approach. \emph{Journal of Molecular Evolution}, \bold{17}, 368--376. 

Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer Associates, Sunderland.

Yang, Z. (2006). \emph{Computational Molecular evolution}. Oxford University Press, Oxford.

Adachi, J., P. J. Waddell, W. Martin, and M. Hasegawa (2000) 
Plastid genome phylogeny and a model of amino acid substitution for proteins
encoded by chloroplast DNA.  \emph{Journal of Molecular Evolution}, \bold{50}, 348--358                             

Rota-Stabelli, O., Z. Yang, and M. Telford. (2009) MtZoa: a general mitochondrial 
amino acid substitutions model for animal evolutionary studies. \emph{Mol. Phyl. Evol}, \bold{52(1)}, 268--72       

Whelan, S. and Goldman, N. (2001) A general empirical model of 
protein evolution derived from multiple protein families using 
a maximum-likelihood approach. \emph{Molecular Biology and Evolution},  \bold{18}, 691--699                       

Le, S.Q. and Gascuel, O. (2008) LG: An Improved, General 
Amino-Acid Replacement Matrix \emph{Molecular Biology and Evolution}, \bold{25(7)}, 1307--1320                     

Yang, Z., R. Nielsen, and M. Hasegawa (1998) Models of amino acid                                                
substitution and applications to Mitochondrial protein evolution. 
\emph{Molecular Biology and Evolution}, \bold{15}, 1600--1611

Abascal, F., D. Posada, and R. Zardoya (2007) MtArt: A new Model of amino acid 
replacement for Arthropoda. \emph{Molecular Biology and Evolution}, \bold{24}, 1--5                                          

Kosiol, C, and Goldman, N (2005) Different versions of the Dayhoff rate matrix -                               
\emph{Molecular Biology and Evolution}, \bold{22}, 193--199 
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{
\code{\link{bootstrap.pml}}, \code{\link{modelTest}}, \code{\link{pmlPart}}, \code{\link{pmlMix}}, \code{\link{plot.phylo}}, \code{\link{SH.test}}
}
% \note{For small trees the likelihood seems to be very similar to Paup* or PhyML.}
\examples{
  example(NJ)
# Jukes-Cantor (starting tree from NJ)  
  fitJC <- pml(tree, Laurasiatherian)  
# optimize edge length parameter     
  fitJC <- optim.pml(fitJC)
  fitJC 
  
\dontrun{    
# search for a better tree using NNI rearrangements     
  fitJC <- optim.pml(fitJC, optNni=TRUE)
  fitJC   
  plot(fitJC$tree)

# JC + Gamma + I - model
  fitJC_GI <- update(fitJC, k=4, inv=.2)
# optimize shape parameter + proportion of invariant sites     
  fitJC_GI <- optim.pml(fitJC_GI, optGamma=TRUE, optInv=TRUE)
# GTR + Gamma + I - model
  fitGTR <- optim.pml(fitJC_GI, optNni=TRUE, optGamma=TRUE, optInv=TRUE, model="GTR") 
}


# 2-state data (RY-coded)  
  dat <- acgt2ry(Laurasiatherian) 
  fit2ST <- pml(tree, dat) 
  fit2ST <- optim.pml(fit2ST,optNni=TRUE) 
  fit2ST
# show some of the methods available for class pml
  methods(class="pml")  
}
\keyword{ cluster }% at least one, from doc/KEYWORDS
