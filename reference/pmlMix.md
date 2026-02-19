# Phylogenetic mixture model

Phylogenetic mixture model.

## Usage

``` r
pmlMix(formula, fit, m = 2, omega = rep(1/m, m),
  control = pml.control(epsilon = 1e-08, maxit = 20), ...)
```

## Arguments

- formula:

  a formula object (see details).

- fit:

  an object of class `pml`.

- m:

  number of mixtures.

- omega:

  mixing weights.

- control:

  A list of parameters for controlling the fitting process.

- ...:

  Further arguments passed to or from other methods.

## Value

`pmlMix` returns a list with elements

- logLik:

  log-likelihood of the fit

- omega:

  mixing weights.

- fits:

  fits for the final mixtures.

## Details

The `formula` object allows to specify which parameter get optimized.
The formula is generally of the form
`edge + bf + Q ~ rate + shape + ...{}`, on the left side are the
parameters which get optimized over all mixtures, on the right the
parameter which are optimized specific to each mixture. The parameters
available are `"nni", "bf", "Q", "inv", "shape", "edge", "rate"`. Each
parameters can be used only once in the formula. `"rate"` and `"nni"`
are only available for the right side of the formula. On the other hand
parameters for invariable sites are only allowed on the left-hand side.
The convergence of the algorithm is very slow and is likely that the
algorithm can get stuck in local optima.

## See also

[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md),[`pmlPart`](https://klausvigo.github.io/phangorn/reference/pmlPart.md),[`pmlCluster`](https://klausvigo.github.io/phangorn/reference/pmlCluster.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
fit <- pml(tree,X, k=4)
weights <- 1000*exp(fit$siteLik)
attr(X, "weight") <- weights
fit1 <- update(fit, data=X, k=1)
fit2 <- update(fit, data=X)

(fitMixture <- pmlMix(edge~rate, fit1 , m=4))
(fit2 <- optim.pml(fit2, optGamma=TRUE))


data(Laurasiatherian)
dm <- dist.logDet(Laurasiatherian)
tree <- NJ(dm)
fit <- pml(tree, Laurasiatherian)
fit <- optim.pml(fit)

fit2 <- update(fit, k=4, site.rate="free_rate")
fit2 <- optim.pml(fit2, optGamma=TRUE)

fitMix <- pmlMix(edge ~ rate, fit, m=4)
fitMix


#
# simulation of mixture models
#
X <- allSitePattern(5)
tree0 <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
tree1 <- read.tree(text = "((t1:0.1,t2:0.5):0.1,(t3:0.1,t4:0.5):0.1,t5:0.5);")
tree2 <- read.tree(text = "((t1:0.5,t2:0.1):0.1,(t3:0.5,t4:0.1):0.1,t5:0.5);")
tree1 <- unroot(tree1)
tree2 <- unroot(tree2)
fit1 <- pml(tree1,X)
fit2 <- pml(tree2,X)

weights <- 2000*exp(fit1$siteLik) + 1000*exp(fit2$siteLik)
attr(X, "weight") <- weights

fit0 <- pml(tree0, X)
fit1 <- pml(tree1, X)
fit2 <- optim.pml(fit1)
logLik(fit2)
AIC(fit2, k=log(3000))

fitMixEdge <- pmlMix( ~ edge, fit1, m=2)
logLik(fitMixEdge)
AIC(fitMixEdge, k=log(3000))

fit.p <- pmlPen(fitMixEdge, .25)
logLik(fit.p)
AIC(fit.p, k=log(3000))
} # }
```
