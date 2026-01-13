# Maximum likelihood by hand

## Maximum likelihood by hand

With the function `pml_bb` from *phangorn* (Schliep 2011) a lot of steps
have become easier and shorter. If you want to have more control over
all of the used parameters, it is also possible to use the older
functions, e.g. `optim_pml`. The data is the same as in the vignette
*Estimating phylogenetic trees with phangorn*:

``` r
library(ape)
library(phangorn)
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")
```

As a starting tree, we calculate a neighbor joining tree:

``` r
dm <- dist.ml(primates)
treeNJ  <- NJ(dm)
```

``` r
fit <- pml(treeNJ, data=primates)
fit
```

    ## model: JC 
    ## loglikelihood: -3075 
    ## unconstrained loglikelihood: -1230 
    ## 
    ## Rates:
    ## a <-> c : 1 
    ## a <-> g : 1 
    ## a <-> t : 1 
    ## c <-> g : 1 
    ## c <-> t : 1 
    ## g <-> t : 1 
    ## 
    ## Base frequencies:  
    ##    a    c    g    t 
    ## 0.25 0.25 0.25 0.25

The function `pml` returns an object of class `pml`. This object
contains the data, the tree and many different parameters of the model
like the likelihood. There are many generic functions for the class
`pml` available, which allow the handling of these objects.

``` r
methods(class="pml")
```

    ##  [1] AICc     anova    BIC      glance   logLik   plot     print    simSeq  
    ##  [9] terraces update   vcov    
    ## see '?methods' for accessing help and source code

The object fit just estimated the likelihood for the tree it got
supplied, but the branch length are not optimized for the Jukes-Cantor
(Jukes and Cantor 1969) model yet, which can be done with the function
`optim.pml`.

``` r
fitJC  <- optim.pml(fit, rearrangement="NNI")
logLik(fitJC)
```

    ## 'log Lik.' -3068 (df=25)

With the default values `pml` will estimate a Jukes-Cantor model. That
means equal base frequencies and all transition rates are equal. The
generic function `update` allows to change parameters manually. This is
not what we usually want to do. However we might want to supply a
different tree or change the number of rate categories.

``` r
fitF81 <- update(fitJC, k=4, inv=0.2, bf=baseFreq(primates))
fitF81
```

    ## model: F81+G(4)+I 
    ## loglikelihood: -3037 
    ## unconstrained loglikelihood: -1230 
    ## Proportion of invariant sites: 0.2 
    ## Model of rate heterogeneity: Discrete gamma model
    ## Number of rate categories: 4 
    ## Shape parameter: 1 
    ##     Rate Proportion
    ## 1 0.0000        0.2
    ## 2 0.1712        0.2
    ## 3 0.5959        0.2
    ## 4 1.2500        0.2
    ## 5 2.9829        0.2
    ## 
    ## Rates:
    ## a <-> c : 1 
    ## a <-> g : 1 
    ## a <-> t : 1 
    ## c <-> g : 1 
    ## c <-> t : 1 
    ## g <-> t : 1 
    ## 
    ## Base frequencies:  
    ##       a       c       g       t 
    ## 0.37481 0.40160 0.03911 0.18448

In the line above we changed the model to a (discrete) rate across site
model with 4 rate categories (using the default shape parameter of 1),
to 0.2 invariant sites and supply empirical base frequencies.

``` r
fitGTR <- optim.pml(fitF81, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "NNI")
fitGTR
```

    ## model: GTR+G(4)+I 
    ## loglikelihood: -2611 
    ## unconstrained loglikelihood: -1230 
    ## Proportion of invariant sites: 0.006978 
    ## Model of rate heterogeneity: Discrete gamma model
    ## Number of rate categories: 4 
    ## Shape parameter: 3.081 
    ##     Rate Proportion
    ## 1 0.0000   0.006978
    ## 2 0.3982   0.248256
    ## 3 0.7411   0.248256
    ## 4 1.0905   0.248256
    ## 5 1.7982   0.248256
    ## 
    ## Rates:
    ## a <-> c : 0.9472 
    ## a <-> g : 63.59 
    ## a <-> t : 0.807 
    ## c <-> g : 0.003986 
    ## c <-> t : 24.63 
    ## g <-> t : 1 
    ## 
    ## Base frequencies:  
    ##       a       c       g       t 
    ## 0.37481 0.40160 0.03911 0.18448

We will change the model to the GTR + $\Gamma(4)$ + I model and then
optimize all the parameters.

With the control parameters the thresholds for the fitting process can
be changed. Here we want just to suppress output during the fitting
process. For larger trees the NNI rearrangements often get stuck in a
local maximum. We added two stochastic algorithms to improve topology
search. The first (set `rearrangement="stochastic"`) performs stochastic
rearrangements similar as in (Nguyen et al. 2015), which makes random
NNI permutation to the tree, which than gets optimized to escape local
optima. The second option (`rearrangement="ratchet"`) perform the
likelihood ratchet (Vos 2003).

While these algorithms may find better trees they will also take more
time.

``` r
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "stochastic")
```

``` r
fitGTR
```

    ## model: GTR+G(4)+I 
    ## loglikelihood: -2608 
    ## unconstrained loglikelihood: -1230 
    ## Proportion of invariant sites: 0.00741 
    ## Model of rate heterogeneity: Discrete gamma model
    ## Number of rate categories: 4 
    ## Shape parameter: 2.994 
    ##     Rate Proportion
    ## 1 0.0000    0.00741
    ## 2 0.3917    0.24815
    ## 3 0.7366    0.24815
    ## 4 1.0906    0.24815
    ## 5 1.8110    0.24815
    ## 
    ## Rates:
    ## a <-> c : 0.7197 
    ## a <-> g : 73.86 
    ## a <-> t : 0.5978 
    ## c <-> g : 0.003286 
    ## c <-> t : 25.85 
    ## g <-> t : 1 
    ## 
    ## Base frequencies:  
    ##       a       c       g       t 
    ## 0.37481 0.40160 0.03911 0.18448

### Model comparison

We can compare nested models for the JC and GTR + $\Gamma(4)$ + I model
using likelihood ratio statistic

``` r
anova(fitJC, fitGTR)
```

    ## Likelihood Ratio Test Table
    ##   Log lik. Df Df change Diff log lik. Pr(>|Chi|)    
    ## 1    -3068 25                                       
    ## 2    -2608 35        10           921     <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

with the Shimodaira-Hasegawa test

``` r
SH.test(fitGTR, fitJC)
```

    ##      Trees  ln L Diff ln L p-value
    ## [1,]     1 -2608       0.0  0.5023
    ## [2,]     2 -3068     460.5  0.0000

or with the AIC

``` r
AIC(fitJC)
```

    ## [1] 6187

``` r
AIC(fitGTR)
```

    ## [1] 5286

``` r
AICc(fitGTR)
```

    ## [1] 5298

``` r
BIC(fitGTR)
```

    ## [1] 5406

### Bootstrap

At last we may want to apply standard bootstrap to test how well the
edges of the tree are supported. This has already been shown in the
vignette *Estimating phylogenetic trees with phangorn*.

``` r
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE)
```

Now we can plot the tree with the bootstrap support values on the edges
and also look at `consensusNet` to identify potential conflict.

``` r
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
```

![Tree with bootstrap support. Unrooted tree (midpoint rooted) with
bootstrap support values.](MLbyHand_files/figure-html/plotBS-1.png)

Tree with bootstrap support. Unrooted tree (midpoint rooted) with
bootstrap support values.

``` r
cnet <- consensusNet(bs, p=0.2)
plot(cnet, show.edge.label=TRUE)
```

![ConsensusNet from the bootstrap
sample.](MLbyHand_files/figure-html/ConsensusNet-1.png)

ConsensusNet from the bootstrap sample.

## Generating trees

*phangorn* has several functions to generate tree topologies, which may
are interesting for simulation studies. `allTrees` computes all possible
bifurcating tree topologies either rooted or unrooted for up to 10 taxa.
One has to keep in mind that the number of trees is growing
exponentially, use `howmanytrees` from *ape* as a reminder.

``` r
trees <- allTrees(5)
par(mfrow=c(3,5), mar=rep(0,4))
for(i in 1:15)plot(trees[[i]], cex=1, type="u")
```

![All 15 unrooted trees with five tip
labels.](MLbyHand_files/figure-html/allTrees-1.png)

All 15 unrooted trees with five tip labels.

`nni` returns a list of all trees which are one nearest neighbor
interchange away.

``` r
nni(trees[[1]])
```

    ## 4 phylogenetic trees

`rNNI` and `rSPR` generate trees which are a defined number of NNI
(nearest neighbor interchange) or SPR (subtree pruning and regrafting)
away.

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
    ## [1] future_1.68.0     phangorn_2.12.1.3 ape_5.8-1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Matrix_1.7-4        future.apply_1.20.1 jsonlite_2.0.0     
    ##  [4] compiler_4.5.2      Rcpp_1.1.1          parallel_4.5.2     
    ##  [7] jquerylib_0.1.4     globals_0.18.0      systemfonts_1.3.1  
    ## [10] textshaping_1.0.4   yaml_2.3.12         fastmap_1.2.0      
    ## [13] lattice_0.22-7      R6_2.6.1            generics_0.1.4     
    ## [16] igraph_2.2.1        knitr_1.51          htmlwidgets_1.6.4  
    ## [19] backports_1.5.0     checkmate_2.3.3     desc_1.4.3         
    ## [22] bslib_0.9.0         rlang_1.1.7         cachem_1.1.0       
    ## [25] xfun_0.55           fs_1.6.6            sass_0.4.10        
    ## [28] otel_0.2.0          cli_3.6.5           pkgdown_2.2.0      
    ## [31] magrittr_2.0.4      digest_0.6.39       grid_4.5.2         
    ## [34] lifecycle_1.0.5     nlme_3.1-168        evaluate_1.0.5     
    ## [37] listenv_0.10.0      codetools_0.2-20    ragg_1.5.0         
    ## [40] parallelly_1.46.1   rmarkdown_2.30      pkgconfig_2.0.3    
    ## [43] tools_4.5.2         htmltools_0.5.9

## References

Jukes, Thomas H., and Charles R. Cantor. 1969. “{CHAPTER} 24 - Evolution
of Protein Molecules.” In *Mammalian Protein Metabolism*, edited by H.
N. Munro, 21–132. Academic Press.

Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang
Minh. 2015. “IQ-TREE: A Fast and Effective Stochastic Algorithm for
Estimating Maximum-Likelihood Phylogenies.” *Molecular Biology and
Evolution* 32 (1): 268–74. <https://doi.org/10.1093/molbev/msu300>.

Schliep, Klaus Peter. 2011. “Phangorn: Phylogenetic Analysis in R.”
*Bioinformatics* 27 (4): 592–93.
<https://doi.org/10.1093/bioinformatics/btq706>.

Vos, R. A. 2003. “Accelerated Likelihood Surface Exploration: The
Likelihood Ratchet.” *Systematic Biology* 52 (3): 368–73.
