# ModelTest

Comparison of different nucleotide or amino acid substitution models

## Usage

``` r
modelTest(object, tree = NULL, model = NULL, G = TRUE, I = TRUE,
  FREQ = FALSE, R = FALSE, k = 4, control = pml.control(),
  multicore = FALSE, mc.cores = NULL)
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

- R:

  logical, TRUE (default) if free rate model should be tested.

- k:

  number of rate classes.

- control:

  A list of parameters for controlling the fitting process.

- multicore:

  Currently not used.

- mc.cores:

  Currently not used.

## Value

A data.frame containing the log-likelihood, number of estimated
parameters, AIC, AICc and BIC all tested models. The data.frame has an
attributes "env" which is an environment which contains all the trees,
the data and the calls to allow get the estimated models, e.g. as a
starting point for further analysis (see example).

## Details

`modelTest` estimates all the specified models for a given tree and
data. When the mclapply is available, the computations are done in
parallel. `modelTest` runs each model in one thread. This is may not
work within a GUI interface and will not work under Windows.

model="3Di" tests the three models for structural alignments mention in
Garg(2025).

## References

Burnham, K. P. and Anderson, D. R (2002) *Model selection and multimodel
inference: a practical information-theoretic approach*. 2nd ed.
Springer, New York

Posada, D. and Crandall, K.A. (1998) MODELTEST: testing the model of DNA
substitution. *Bioinformatics* **14(9)**: 817-818

Posada, D. (2008) jModelTest: Phylogenetic Model Averaging. *Molecular
Biology and Evolution* **25**: 1253-1256

Darriba D., Taboada G.L., Doallo R and Posada D. (2011) ProtTest 3: fast
selection of best-fit models of protein evolution. . *Bioinformatics*
**27**: 1164-1165

Puente-Lelievre, Caroline, et al. (2023) "Tertiary-interaction
characters enable fast, model-based structural phylogenetics beyond the
twilight zone." *bioRxiv*: 2023-12

Garg, S.G., Hochberg, G.K.A. (2025) A General Substitution Matrix for
Structural Phylogenetics, *Molecular Biology and Evolution*, **42(6)**,
https://doi.org/10.1093/molbev/msaf124

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
mT <- modelTest(Laurasiatherian, model = c("JC", "K80", "HKY", "GTR"),
                R=TRUE)

# Some exploratory data analysis
plot(mT$TL, mT$logLik, xlim=c(3,6.5))
text(mT$TL, mT$logLik, labels=mT$Model, pos=4)

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
