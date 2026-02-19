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
data(woodmouse)
pratchet(woodmouse)
#> Parsimony score of initial tree: 68 
#> Iteration: 10. Best parsimony score so far: 68Iteration: 20. Best parsimony score so far: 68Iteration: 30. Best parsimony score so far: 68Iteration: 40. Best parsimony score so far: 68Iteration: 50. Best parsimony score so far: 68Iteration: 60. Best parsimony score so far: 68Iteration: 70. Best parsimony score so far: 68Iteration: 80. Best parsimony score so far: 68Iteration: 90. Best parsimony score so far: 68Iteration: 100. Best parsimony score so far: 68
#> 
#> Phylogenetic tree with 15 tips and 13 internal nodes.
#> 
#> Tip labels:
#>   No305, No304, No306, No0906S, No0908S, No0909S, ...
#> Node labels:
#>   1, 0.73, 0.88, 0.8, 0.99, 0.7, ...
#> 
#> Unrooted; no branch length.
trs <- pratchet(woodmouse, all=TRUE)
#> Parsimony score of initial tree: 68 
#> Iteration: 10. Best parsimony score so far: 68Iteration: 20. Best parsimony score so far: 68Iteration: 30. Best parsimony score so far: 68Iteration: 40. Best parsimony score so far: 68Iteration: 50. Best parsimony score so far: 68Iteration: 60. Best parsimony score so far: 68Iteration: 70. Best parsimony score so far: 68Iteration: 80. Best parsimony score so far: 68Iteration: 90. Best parsimony score so far: 68Iteration: 100. Best parsimony score so far: 68
start_trs <- get("start_trees", envir = attr(trs, "env"))
terraces(woodmouse, c(trs, start_trs))
#> Error in UseMethod("terraces"): no applicable method for 'terraces' applied to an object of class "DNAbin"

if (FALSE) { # \dontrun{
fit <- pml_bb(woodmouse, model="JC")
terraces(fit, dist_fun="KF.dist")
terraces(fit, pkg="scatterplot3d")
terraces(fit, pkg="plot3D")
terraces(fit, pkg="rgl")
} # }
```
