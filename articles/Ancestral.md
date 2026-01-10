# Ancestral Sequence Reconstruction

## Introduction

These notes describe the ancestral sequence reconstruction using the
*phangorn* package (Schliep 2011). *phangorn* provides several methods
to estimate ancestral character states with either Maximum Parsimony
(MP) or Maximum Likelihood (ML). For more background on all the methods
see e.g. (Felsenstein 2004) or (Yang 2006).

## Parsimony reconstructions

To reconstruct ancestral sequences we first load some data and
reconstruct a tree:

``` r
library(phangorn)
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")
tree <- pratchet(primates, trace=0) |> acctran(primates) |> makeNodeLabel()
parsimony(tree, primates)
```

    ## [1] 746

For parsimony analysis of the edge length represent the observed number
of changes. Reconstructing ancestral states therefore defines also the
edge lengths of a tree. However there can exist several equally
parsimonious reconstructions or states can be ambiguous and therefore
edge length can differ. In *phangorn* all the ancestral reconstructions
for parsimony used to be based on the fitch algorithm (below version
3.0) and needed bifurcating trees. However trees can get pruned
afterwards using the function `multi2di` from *ape*. Recently we
replaced the acctran routine with a method based on the sankoff
algorithm adopting the algorithm for joint reconstruction (Pupko et al.
2000) and breaking ties at random. This has the additional benefit that
it allows us to infer phylogenies with multifurcations.

“MPR” reconstructs the ancestral states for each (internal) node as if
the tree would be rooted in that node. However the nodes are not
independent of each other. If one chooses one state for a specific node,
this can restrict the choice of neighboring nodes (figures 2 and 3).
There is also an option “POSTORDER” which is only a one pass algorithm,
which is useful for teaching purposes. The function acctran (accelerated
transformation) assigns edge length and internal nodes to the tree
(Swofford and Maddison 1987).

``` r
anc.pars <- anc_pars(tree, primates)
```

The `plotSeqLogo` function is a wrapper around the from the *ggseqlogo*
function in the *ggseqlogo* package (Wagih 2024) and provides a simple
way to show proportions of a nucleotides of ancestral states (see figure
1).

``` r
plotSeqLogo(anc.pars, node=getRoot(tree), 1, 20)
```

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## ℹ The deprecated feature was likely used in the ggseqlogo package.
    ##   Please report the issue at <https://github.com/omarwagih/ggseqlogo/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![Fig 1. Ancestral reconstruction for a
node.](Ancestral_files/figure-html/seqLogo-1.png)

Fig 1. Ancestral reconstruction for a node.

``` r
plotAnc(anc.pars, 17)
title("MPR")
```

![Fig 2. Ancestral reconstruction using
MPR.](Ancestral_files/figure-html/MPR-1.png)

Fig 2. Ancestral reconstruction using MPR.

## Likelihood reconstructions

*phangorn* also offers the possibility to estimate ancestral states
using ML. The advantages of ML over parsimony is that the reconstruction
accounts for different edge lengths. Currently marginal construction
(see (Yang 2006; Koshi and Goldstein 1996)) and joint reconstruction
(Pupko et al. 2000) are implemented. Joint reconstructions is only
available for models without rate variation (e.g. gamma models) or
invariant sites.

``` r
fit <- pml(tree, primates)
fit <- optim.pml(fit, model="F81")
```

We can assign the ancestral states according to the highest likelihood
(“ml”):
$$P\left( x_{r} = A \right) = \frac{L\left( x_{r} = A \right)}{\sum\limits_{k \in \{ A,C,G,T\}}L\left( x_{r} = k \right)}$$
and the highest posterior probability (“bayes”) criterion:
$$P\left( x_{r} = A \right) = \frac{\pi_{A}L\left( x_{r} = A \right)}{\sum\limits_{k \in \{ A,C,G,T\}}\pi_{k}L\left( x_{r} = k \right)},$$
where $L\left( x_{r} \right)$ is the joint probability of states at the
tips and the state at the root $x_{r}$ and $\pi_{i}$ are the estimated
base frequencies of state $i$. Both methods agree if all states (base
frequencies) have equal probabilities.

``` r
anc.ml <- anc_pml(fit)
```

The differences of the two approaches for a specific site (17) are
represented in the following figures.

``` r
plotAnc(anc.ml, 17)
title("ML")
```

![Fig 4. Ancestral reconstruction the using the maximum
likelihood.](Ancestral_files/figure-html/plotML-1.png)

Fig 4. Ancestral reconstruction the using the maximum likelihood.

``` r
#plotAnc(anc.bayes, 17)
#title("Bayes")
```

## Fitting for discrete comparative data

Often have already a phylogeny and only want estimate the ancestral
reconstruction for this tree. This is a common problem in phylogentic
comparative methods and we can use the function *ace* in the ape
(Paradis and Schliep 2019), *fitDiscrete* in the geiger (Pennell et al.
2014) or *fitMK* in the phytools (Revell 2012) package. Here we want to
show how to fit these models using *optim.pml*.

First we load a tree and create some data.

``` r
data("bird.orders")
x <- c(rep(0, 5), rep(1, 18))
x[c(20,22,23)] <- 2
x <- factor(x)
names(x) <- bird.orders$tip.label
dat <- phyDat(x, "USER", levels=c(0,1,2))
```

We than set up the *pml* object and optimize the model. Instead of
optimizing the edge length we only optimize the rate.

``` r
fit <- pml(bird.orders, dat)
fit_ER <- optim.pml(fit, optEdge = FALSE, optRate=TRUE)
fit_ER
```

    ## model: Mk 
    ## loglikelihood: -16.47 
    ## unconstrained loglikelihood: 0 
    ## 
    ## Rate matrix:
    ##   0 1 2
    ## 0 0 1 1
    ## 1 1 0 1
    ## 2 1 1 0
    ## 
    ## Base frequencies:  
    ##      0      1      2 
    ## 0.3333 0.3333 0.3333 
    ## 
    ## Rate: 0.007846

We can also fit the symmetric (*model=“SYM”*) or ordered metristic model
(*model=“ORDERED”*).

``` r
fit_SYM <- optim.pml(fit, optEdge = FALSE, optRate=TRUE, model="SYM")
fit_SYM
```

    ## model: SYM 
    ## loglikelihood: -15.31 
    ## unconstrained loglikelihood: 0 
    ## 
    ## Rate matrix:
    ##           0      1         2
    ## 0 0.000e+00 0.2747 1.604e-06
    ## 1 2.747e-01 0.0000 1.000e+00
    ## 2 1.604e-06 1.0000 0.000e+00
    ## 
    ## Base frequencies:  
    ##      0      1      2 
    ## 0.3333 0.3333 0.3333 
    ## 
    ## Rate: 0.00678

We can compare the estimate with the one from *ace* from *ape*.

``` r
fit_ace <- ace(x, bird.orders, model="SYM", type = "d")
```

    ## Warning in sqrt(diag(solve(h))): NaNs produced

The log-likelihood values differ slightly as in phangorn the values get
multiplied by the state frequencies. Thus if we add `log(1/3)` as we
have three states to ace estimate the two estimates are almost
identical.

``` r
fit_SYM$logLik
```

    ## [1] -15.31

``` r
fit_ace$loglik+log(1/3)
```

    ## [1] -15.31

``` r
all.equal(fit_SYM$logLik, fit_ace$loglik+log(1/3))
```

    ## [1] "Mean relative difference: 1.229e-07"

``` r
anc_SYM <- anc_pml(fit_SYM)
plotAnc(anc_SYM)
```

![Marginal reconstruciton for comparative data
example.](Ancestral_files/figure-html/SYM_reconstruction-1.png)

Marginal reconstruciton for comparative data example.

More complicated models can be applied using defining the rate matrix as
shown in the vignette *Markov models and transition rate matrices*. The
“ARD” model is currently not available as *phangorn* only fits
reversible models.

## Session info

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] phangorn_2.12.1.3 ape_5.8-1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sass_0.4.10         future_1.68.0       generics_0.1.4     
    ##  [4] lattice_0.22-7      listenv_0.10.0      digest_0.6.39      
    ##  [7] magrittr_2.0.4      evaluate_1.0.5      grid_4.5.2         
    ## [10] RColorBrewer_1.1-3  fastmap_1.2.0       jsonlite_2.0.0     
    ## [13] Matrix_1.7-4        backports_1.5.0     scales_1.4.0       
    ## [16] codetools_0.2-20    textshaping_1.0.4   jquerylib_0.1.4    
    ## [19] cli_3.6.5           rlang_1.1.7         parallelly_1.46.1  
    ## [22] future.apply_1.20.1 withr_3.0.2         cachem_1.1.0       
    ## [25] yaml_2.3.12         otel_0.2.0          tools_4.5.2        
    ## [28] ggseqlogo_0.2.2     parallel_4.5.2      checkmate_2.3.3    
    ## [31] dplyr_1.1.4         ggplot2_4.0.1       fastmatch_1.1-6    
    ## [34] globals_0.18.0      vctrs_0.6.5         R6_2.6.1           
    ## [37] lifecycle_1.0.5     fs_1.6.6            htmlwidgets_1.6.4  
    ## [40] ragg_1.5.0          pkgconfig_2.0.3     desc_1.4.3         
    ## [43] pillar_1.11.1       pkgdown_2.2.0       bslib_0.9.0        
    ## [46] gtable_0.3.6        glue_1.8.0          Rcpp_1.1.0         
    ## [49] systemfonts_1.3.1   tidyselect_1.2.1    tibble_3.3.0       
    ## [52] xfun_0.55           knitr_1.51          farver_2.1.2       
    ## [55] htmltools_0.5.9     nlme_3.1-168        igraph_2.2.1       
    ## [58] labeling_0.4.3      rmarkdown_2.30      compiler_4.5.2     
    ## [61] S7_0.2.1

## References

Felsenstein, Joseph. 2004. *Inferring Phylogenies*. Sunderland: Sinauer
Associates.

Koshi, Jeffrey M., and Richard A. Goldstein. 1996. “Probabilistic
reconstruction of ancestral protein sequences.” *Journal of Molecular
Evolution* 42 (2): 313–20. <https://doi.org/10.1007/BF02198858>.

Paradis, Emmanuel, and Klaus Schliep. 2019. “Ape 5.0: An Environment for
Modern Phylogenetics and Evolutionary Analyses in r.” *Bioinformatics*
35 (3): 526–28. <https://doi.org/10.1093/bioinformatics/bty633>.

Pennell, M. W., J. M. Eastman, G. J. Slater, J. W. Brown, J. C. Uyeda,
R. G. Fitzjohn, M. E. Alfaro, and L. J. Harmon. 2014. “Geiger V2.0: An
Expanded Suite of Methods for Fitting Macroevolutionary Models to
Phylogenetic Trees.” *Bioinformatics* 30: 2216–18.
<https://doi.org/10.1093/bioinformatics/btu181>.

Pupko, Tal, Itsik Pe, Ron Shamir, and Dan Graur. 2000. “A Fast Algorithm
for Joint Reconstruction of Ancestral Amino Acid Sequences.” *Molecular
Biology and Evolution* 17 (6): 890–96.
<https://doi.org/10.1093/oxfordjournals.molbev.a026369>.

Revell, Liam J. 2012. “Phytools: An r Package for Phylogenetic
Comparative Biology (and Other Things).” *Methods in Ecology and
Evolution* 3: 217–23.

Schliep, Klaus Peter. 2011. “Phangorn: Phylogenetic Analysis in R.”
*Bioinformatics* 27 (4): 592–93.
<https://doi.org/10.1093/bioinformatics/btq706>.

Swofford, D. L., and W. P. Maddison. 1987. “Reconstructing Ancestral
Character States Under Wagner Parsimony.” *Math. Biosci.* 87: 199–229.

Wagih, Omar. 2024. *Ggseqlogo: A ’Ggplot2’ Extension for Drawing
Publication-Ready Sequence Logos*.
<https://CRAN.R-project.org/package=ggseqlogo>.

Yang, Ziheng. 2006. *Computational Molecular Evolution*. Oxford: Oxford
University Press.
