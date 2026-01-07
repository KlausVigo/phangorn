# Shimodaira-Hasegawa Test

This function computes the Shimodaira–Hasegawa test for a set of trees.

## Usage

``` r
SH.test(..., B = 10000, data = NULL, weight = NULL)
```

## Arguments

- ...:

  either a series of objects of class `"pml"` separated by commas, a
  list containing such objects or an object of class `"pmlPart"` or a
  matrix containing the site-wise likelihoods in columns.

- B:

  the number of bootstrap replicates.

- data:

  an object of class `"phyDat"`.

- weight:

  if a matrix with site (log-)likelihoods is is supplied an optional
  vector containing the number of occurrences of each site pattern.

## Value

a numeric vector with the P-value associated with each tree given in
`...`.

## References

Shimodaira, H. and Hasegawa, M. (1999) Multiple comparisons of
log-likelihoods with applications to phylogenetic inference. *Molecular
Biology and Evolution*, **16**, 1114–1116.

## See also

[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md),
[`pmlPart`](https://klausvigo.github.io/phangorn/reference/pmlPart.md),
[`pmlCluster`](https://klausvigo.github.io/phangorn/reference/pmlCluster.md),
[`SOWH.test`](https://klausvigo.github.io/phangorn/reference/SOWH.test.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(Laurasiatherian)
dm <- dist.logDet(Laurasiatherian)
tree1 <- NJ(dm)
tree2 <- unroot(upgma(dm))
fit1 <- pml(tree1, Laurasiatherian)
fit2 <- pml(tree2, Laurasiatherian)
fit1 <- optim.pml(fit1) # optimize edge weights
#> optimize edge weights:  -54807.68 --> -54290.26 
#> optimize edge weights:  -54290.26 --> -54290.26 
#> optimize edge weights:  -54290.26 --> -54290.26 
fit2 <- optim.pml(fit2)
#> optimize edge weights:  -55623.41 --> -54911.33 
#> optimize edge weights:  -54911.33 --> -54911.33 
#> optimize edge weights:  -54911.33 --> -54911.33 
# with pml objects as input
SH.test(fit1, fit2, B=1000)
#>      Trees      ln L Diff ln L p-value
#> [1,]     1 -54290.26    0.0000   0.493
#> [2,]     2 -54911.33  621.0767   0.000
# in real analysis use larger B, e.g. 10000

# with matrix as input
X <- matrix(c(fit1$siteLik, fit2$siteLik), ncol=2)
SH.test(X, weight=attr(Laurasiatherian, "weight"), B=1000)
#>      Trees      ln L Diff ln L p-value
#> [1,]     1 -54290.26    0.0000   0.489
#> [2,]     2 -54911.33  621.0767   0.000
if (FALSE) { # \dontrun{
example(pmlPart)
SH.test(sp, B=1000)
} # }
```
