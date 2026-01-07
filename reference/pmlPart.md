# Partition model.

Model to estimate phylogenies for partitioned data.

## Usage

``` r
multiphyDat2pmlPart(x, method = "unrooted", tip.dates = NULL, ...)

pmlPart2multiPhylo(x)

pmlPart(formula, object, control = pml.control(epsilon = 1e-08, maxit = 10),
  model = NULL, method = "unrooted", ...)
```

## Arguments

- x:

  an object of class `pmlPart`

- method:

  One of "unrooted", "ultrametric" or "tiplabeled". Only unrooted is
  properly supported right now.

- tip.dates:

  A named vector of sampling times associated to the tips/sequences.
  Leave empty if not estimating tip dated phylogenies.

- ...:

  Further arguments passed to or from other methods.

- formula:

  a formula object (see details).

- object:

  an object of class `pml` or a list of objects of class `pml` .

- control:

  A list of parameters for controlling the fitting process.

- model:

  A vector containing the models containing a model for each partition.

## Value

`kcluster` returns a list with elements

- logLik:

  log-likelihood of the fit

- trees:

  a list of all trees during the optimization.

- object:

  an object of class `"pml"` or `"pmlPart"`

## Details

The `formula` object allows to specify which parameter get optimized.
The formula is generally of the form
`edge + bf + Q ~ rate + shape + ...{}`, on the left side are the
parameters which get optimized over all partitions, on the right the
parameter which are optimized specific to each partition. The parameters
available are `"nni", "bf", "Q", "inv", "shape", "edge", "rate"`. Each
parameters can be used only once in the formula. `"rate"` is only
available for the right side of the formula.

For partitions with different edge weights, but same topology, `pmlPen`
can try to find more parsimonious models (see example).

`pmlPart2multiPhylo` is a convenience function to extract the trees out
of a `pmlPart` object.

## See also

[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md),[`pmlCluster`](https://klausvigo.github.io/phangorn/reference/pmlCluster.md),[`pmlMix`](https://klausvigo.github.io/phangorn/reference/pmlMix.md),
[`SH.test`](https://klausvigo.github.io/phangorn/reference/SH.test.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(yeast)
dm <- dist.logDet(yeast)
tree <- NJ(dm)
fit <- pml(tree,yeast)
fits <- optim.pml(fit)
#> optimize edge weights:  -737063 --> -734615.7 
#> optimize edge weights:  -734615.7 --> -734615.7 
#> optimize edge weights:  -734615.7 --> -734615.7 

weight=xtabs(~ index+genes,attr(yeast, "index"))[,1:10]

sp <- pmlPart(edge ~ rate + inv, fits, weight=weight)
#> loglik: -61530.38 --> -59834.29 
#> loglik: -59834.29 --> -59833.25 
#> loglik: -59833.25 --> -59833.25 
#> loglik: -59833.25 --> -59833.25 
sp
#> 
#> loglikelihood: -59833.25 
#> 
#> loglikelihood of partitions:
#>   -9827.497 -8159.024 -8056.932 -5237.677 -3809.733 -5503.277 -2752.2 -7200.052 -4632.422 -4654.434 
#> AIC:  119730.5  BIC:  119963.5 
#> 
#> Proportion of invariant sites: 0.400014 0.3179422 0.4746719 0.44901 0.412222 0.2912371 0.2419212 0.3097902 0.4794213 0.3884183 
#> 
#> Rates:
#> 1.108054 0.9656349 0.8692647 0.8928479 0.8097379 1.266643 1.296063 1.212791 0.842597 0.8275076 
#> 
#> Base frequencies:  
#>      [,1] [,2] [,3] [,4]
#> [1,] 0.25 0.25 0.25 0.25
#> 
#> Rate matrix:
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    1    1    1    1    1

if (FALSE) { # \dontrun{
sp2 <- pmlPart(~ edge + inv, fits, weight=weight)
sp2
AIC(sp2)

sp3 <- pmlPen(sp2, lambda = 2)
AIC(sp3)
} # }
```
