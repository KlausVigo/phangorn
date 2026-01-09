# Explore likelihood parsimony surface

`terraces` visualizes in likelihood surface for the tree space
(Sanderson et al. 2011). Usually trees are from a bootstrap or MCMC
sample. There the first two axis are the principle components of
distances between trees and the third axis is the likelihood value or
parsimony score.

## Usage

``` r
terraces(x, ...)

# S3 method for class 'pml'
terraces(x, trees = x$bs, dist_fun = "RF.dist",
  di2multi = FALSE, tol = 2e-08, plot = TRUE, ...)

# S3 method for class 'phyDat'
terraces(x, trees, dist_fun = "RF.dist", di2multi = TRUE,
  tol = 2e-08, plot = TRUE, ...)
```

## Arguments

- x:

  an object of class `pml`

- ...:

  Further arguments passed to or from other methods.

- trees:

  an object of class `multiPhylo`

- dist_fun:

  a function to compute distances between trees see e.g.
  [`RF.dist`](https://klausvigo.github.io/phangorn/reference/treedist.md)

- di2multi:

  logical, should trees multichotomies get collapsed. Useful for
  Robinson-Foulds distance. If edge length are used to compute the
  distance, e.g. Kuhner-Felsenstein distance, this is not needed.

- tol:

  a numeric value giving the tolerance to consider a branch length
  significantly greater than zero.

- plot:

  loggical if TRUE a 3D scatter is shown.

## Value

`terraces` silently returns a matrix.

## References

Sanderson, M.J., McMahon, M.M. and Steel, M. (2011). Terraces in
phylogenetic tree space. *Science*, **333**, 448–450.

## See also

[`pml_bb`](https://klausvigo.github.io/phangorn/reference/pml_bb.md)`, `[`optim.pml`](https://klausvigo.github.io/phangorn/reference/pml.md)`, `[`pratchet`](https://klausvigo.github.io/phangorn/reference/parsimony.md)`, `[`RF.dist`](https://klausvigo.github.io/phangorn/reference/treedist.md)`, `[`di2multi`](https://rdrr.io/pkg/ape/man/multi2di.html)`, `[`cmdscale`](https://rdrr.io/r/stats/cmdscale.html).

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data(woodmouse)
fit <- pml_bb(woodmouse, model="JC")
terraces(fit)
} # }
```
