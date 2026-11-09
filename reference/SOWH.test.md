# Swofford-Olsen-Waddell-Hillis Test

This function computes the Swofford–Olsen–Waddell–Hillis (SOWH) test
[Swofford, Olsen, Waddell, and Hillis (1996)](#reference++Swofford1996);
[Goldman, Anderson, and Rodrigo (2000)](#reference++Goldman2000) , a
parametric bootstrap test. The function is computational very demanding
and likely to be very slow.

## Usage

``` r
SOWH.test(x, n = 100, restricted = list(optNni = FALSE), optNni = TRUE,
  trace = !getOption("quiet"), ...)
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

Goldman N, Anderson JP, Rodrigo AG (2000). “Likelihood-based tests of
topologies in phylogenetics.” *Systematic biology*, **49**(4), 652–670.

Swofford D, Olsen G, Waddell P, Hillis D (1996). “Phylogenetic
Inference.” In Hillis D, Moritz C, Mable B (eds.), *Molecular
Systematics*, Second edition. Sinauer, Sunderland, MA.

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
