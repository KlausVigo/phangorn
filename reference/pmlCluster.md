# Stochastic Partitioning

Stochastic Partitioning of genes into p cluster.

## Usage

``` r
pmlCluster(formula, fit, weight, p = 1:5, part = NULL, nrep = 10,
  control = pml.control(epsilon = 1e-08, maxit = 10, trace = 1), ...)
```

## Arguments

- formula:

  a formula object (see details).

- fit:

  an object of class `pml`.

- weight:

  `weight` is matrix of frequency of site patterns for all genes.

- p:

  number of clusters.

- part:

  starting partition, otherwise a random partition is generated.

- nrep:

  number of replicates for each p.

- control:

  A list of parameters for controlling the fitting process.

- ...:

  Further arguments passed to or from other methods.

## Value

`pmlCluster` returns a list with elements

- logLik:

  log-likelihood of the fit

- trees:

  a list of all trees during the optimization.

- fits:

  fits for the final partitions

## Details

The `formula` object allows to specify which parameter get optimized.
The formula is generally of the form
`edge + bf + Q ~ rate + shape + ...{}`, on the left side are the
parameters which get optimized over all cluster, on the right the
parameter which are optimized specific to each cluster. The parameters
available are `"nni", "bf", "Q", "inv", "shape", "edge", "rate"`. Each
parameter can be used only once in the formula. There are also some
restriction on the combinations how parameters can get used. `"rate"` is
only available for the right side. When `"rate"` is specified on the
left hand side `"edge"` has to be specified (on either side), if
`"rate"` is specified on the right hand side it follows directly that
`edge` is too.

## References

K. P. Schliep (2009). Some Applications of statistical phylogenetics
(PhD Thesis)

Lanfear, R., Calcott, B., Ho, S.Y.W. and Guindon, S. (2012)
PartitionFinder: Combined Selection of Partitioning Schemes and
Substitution Models for Phylogenetic Analyses. *Molecular Biology and
Evolution*, **29(6)**, 1695-1701

## See also

[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md),[`pmlPart`](https://klausvigo.github.io/phangorn/reference/pmlPart.md),[`pmlMix`](https://klausvigo.github.io/phangorn/reference/pmlMix.md),
[`SH.test`](https://klausvigo.github.io/phangorn/reference/SH.test.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data(yeast)
dm <- dist.logDet(yeast)
tree <- NJ(dm)
fit <- pml(tree,yeast)
fit <- optim.pml(fit)

weight <- xtabs(~ index+genes,attr(yeast, "index"))
set.seed(1)

sp <- pmlCluster(edge~rate, fit, weight, p=1:4)
sp
SH.test(sp)
} # }
```
