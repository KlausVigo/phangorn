# Swofford-Olsen-Waddell-Hillis Test

This function computes the Swofford–Olsen–Waddell–Hillis (SOWH) test, a
parametric bootstrap test. The function is computational very demanding
and likely to be very slow.

## Usage

``` r
SOWH.test(x, n = 100, restricted = list(optNni = FALSE), optNni = TRUE,
  trace = 1, ...)
```

## Arguments

- x:

  an object of class `"pml"`.

- n:

  the number of bootstrap replicates.

- restricted:

  list of restricted parameter settings.

- optNni:

  Logical value indicating whether topology gets optimized (NNI).

- trace:

  Show output during computations.

- ...:

  Further arguments passed to `"optim.pml"`.

## Value

an object of class SOWH. That is a list with three elements, one is a
matrix containing for each bootstrap replicate the (log-) likelihood of
the restricted and unrestricted estimate and two pml objects of the
restricted and unrestricted model.

## Details

`SOWH.test` performs a parametric bootstrap test to compare two trees.
It makes extensive use `simSeq` and `optim.pml` and can take quite long.

## References

Goldman, N., Anderson, J. P., and Rodrigo, A. G. (2000) Likelihood
-based tests of topologies in phylogenetics. *Systematic Biology* **49**
652-670.

Swofford, D.L., Olsen, G.J., Waddell, P.J. and Hillis, D.M. (1996)
Phylogenetic Inference in Hillis, D.M., Moritz, C. and Mable, B.K.
(Eds.) *Molecular Systematics* (2nd ed.) 407-514, Sunderland, MA:
Sinauer

## See also

[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md),
[`pmlPart`](https://klausvigo.github.io/phangorn/reference/pmlPart.md),
[`pmlCluster`](https://klausvigo.github.io/phangorn/reference/pmlCluster.md),
[`simSeq`](https://klausvigo.github.io/phangorn/reference/simSeq.md),
[`SH.test`](https://klausvigo.github.io/phangorn/reference/SH.test.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
# in real analysis use larger n, e.g. 500 preferably more
if (FALSE) { # \dontrun{
data(Laurasiatherian)
dm <- dist.logDet(Laurasiatherian)
tree <- NJ(dm)
fit <- pml(tree, Laurasiatherian)
fit <- optim.pml(fit, TRUE)
set.seed(6)
tree <- rNNI(fit$tree, 1)
fit <- update(fit, tree = tree)
(res <- SOWH.test(fit, n=100))
summary(res)
} # }
```
