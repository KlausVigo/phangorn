# ModelTest

Comparison of different nucleotide or amino acid substitution models

## Usage

``` r
modelTest(object, tree = NULL, model = NULL, G = TRUE, I = TRUE,
  FREQ = FALSE, k = 4, control = pml.control(), RHAS = "gamma", ...,
  mt_control = mt.control(crit = "BIC", n_model = 5, n_rhas = 100))

write.modelTest(x, file = "modelTest", digits = 10)
```

## Arguments

- object:

  an object of class phyDat or pml.

- tree:

  a phylogenetic tree.

- model:

  a vector containing the substitution models to compare with each other
  or "all" to test all available models.

- G:

  logical, TRUE (default) if (discrete) Gamma model should be tested.

- I:

  logical, TRUE (default) if invariant sites should be tested.

- FREQ:

  logical, FALSE (default) if TRUE amino acid frequencies will be
  estimated.

- k:

  number of rate classes. Can be a list with a vector for each RHAS
  term.

- control:

  A list of parameters for controlling the fitting process.

- RHAS:

  a character vector specifying the rate heterogeneity among sites
  models. Option are "gamma" for discrete gamma model with equal
  weights, "variable_gamma" for discrete gamma model with estimated
  weights, "free_rate" and "gamma_quadrature".

- ...:

  Further arguments passed to or from other methods.

- mt_control:

  a list with some options.

- x:

  an object of class modelTest.

- file:

  a file name. File endings are added.

- digits:

  default is 10, i.e. edge length for the bootstrap trees are exported.
  For digits larger smaller than zero no edge length are exported.

## Value

A data.frame containing the log-likelihood, number of estimated
parameters, AIC, AICc and BIC all tested models. The data.frame has an
attributes "env" which is an environment which contains all the trees,
the data and the calls to allow get the estimated models, e.g. as a
starting point for further analysis (see example).

## Details

`modelTest` estimates all the specified models for a given tree and data
similar to ([Posada and Crandall 1998](#reference++Posada1998); [Posada
2008](#reference++Posada2008)) .

model="3Di" tests the three models for structural alignments mention in
[Garg and Hochberg (2025)](#reference++Garg2025) .

## References

Garg SG, Hochberg GKA (2025). “A General Substitution Matrix for
Structural Phylogenetics.” *Molecular Biology and Evolution*, **42**(6),
msaf124. ISSN 1537-1719.
[doi:10.1093/molbev/msaf124](https://doi.org/10.1093/molbev/msaf124) .

Posada D (2008). “jModelTest: Phylogenetic Model Averaging.” *Molecular
Biology and Evolution*, **25**(7), 1253–1256.
[doi:10.1093/molbev/msn083](https://doi.org/10.1093/molbev/msn083) .

Posada D, Crandall K (1998). “MODELTEST: testing the model of DNA
substitution.” *Bioinformatics*, **14**(9), 817–818.

Burnham, K. P. and Anderson, D. R (2002) *Model selection and multimodel
inference: a practical information-theoretic approach*. 2nd ed.
Springer, New York

Darriba D., Taboada G.L., Doallo R and Posada D. (2011) ProtTest 3: fast
selection of best-fit models of protein evolution. . *Bioinformatics*
**27**: 1164-1165

Puente-Lelievre, Caroline, et al. (2023) "Tertiary-interaction
characters enable fast, model-based structural phylogenetics beyond the
twilight zone." *bioRxiv*: 2023-12

## See also

[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md),
[`anova`](https://rdrr.io/r/stats/anova.html),
[`AIC`](https://rdrr.io/r/stats/AIC.html),
[`codonTest`](https://klausvigo.github.io/phangorn/reference/codonTest.md),
[`plan`](https://future.futureverse.org/reference/plan.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r

if (FALSE) { # \dontrun{
data(Laurasiatherian)
mT <- modelTest(Laurasiatherian, model = c("JC", "K80", "HKY", "GTR"))

# Some exploratory data analysis
plot(mT$TL, mT$logLik, xlim=c(3,6.5),
   xlab="total tree length", ylab="log likelihood")
text(mT$TL, mT$logLik, labels=mT$Model, pos=4)

fit_GTR_G <- as.pml(mt, "GTR+G(4)")
fit_GTR_VG <- as.pml(mt, "GTR+VG(4)")
fit_GTR_R <- as.pml(mt, "GTR+R(4)")
plotRates(fit_GTR_G)
plotRates(fit_GTR_VWG, append=TRUE, col="red")
plotRates(fit_GTR_R, append=TRUE, col="green")

# extract best model
(best_model <- as.pml(mT))


data(chloroplast)
(mTAA <- modelTest(chloroplast, model=c("JTT", "WAG", "LG")))

# test all available amino acid models
plan(multisession, workers = 2)
(mTAA_all <- modelTest(chloroplast, model="all"))
plan(sequential)
} # }
```
