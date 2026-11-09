# UPGMA, WPGMA and sUPGMA

UPGMA and WPGMA clustering. UPGMA ([Sokal and Michener
1958](#reference++Sokal1958)) and WPGMA ([McQuitty
1966](#reference++McQuitty1966)) are a wrapper function around
[`hclust`](https://rdrr.io/r/stats/hclust.html) returning a `phylo`
object. `supgma` perform serial sampled UPGMA similar to [Drummond and
Rodrigo (2000)](#reference++Drummond2000) .

## Usage

``` r
upgma(D, ...)

wpgma(D, ...)

supgma(D, tip.dates, trace = 0, ...)
```

## Arguments

- D:

  A distance matrix, i.e. an object of class `dist`. If an matrix is
  supplied it is tried to covert it do a `dist` object.

- ...:

  Further arguments passed to or from other methods.

- tip.dates:

  A named vector of sampling times associated to the tips.

- trace:

  Show output during optimization (see details).

## Value

A phylogenetic tree of class `phylo`.

## Details

UPGMA and WPGMA return ultrametric trees, it is implicitly assumed that
the distances supplied are close to ultrametric, e.g. hold the molecular
clock assumption. Neighbor Joining (NJ)
[`nj`](https://rdrr.io/pkg/ape/man/nj.html) and fastME
[`fastme`](https://rdrr.io/pkg/ape/man/fastme.html) relax this
assumption to additive distances. sUPGMA assumes tip dated data.

`tip.dates` should be a vector of sampling times, in any time unit, with
time increasing toward the present. For example, this may be in units of
"days since study start" or "years since 10.000 BCE", but not "millions
of years ago".

## References

Drummond A, Rodrigo AG (2000). “Reconstructing Genealogies of Serial
Samples Under the Assumption of a Molecular Clock Using Serial-Sample
UPGMA.” *Molecular Biology and Evolution*, **17**(12), 1807-1815. ISSN
0737-4038.
[doi:10.1093/oxfordjournals.molbev.a026281](https://doi.org/10.1093/oxfordjournals.molbev.a026281)
.

McQuitty LL (1966). “Similarity analysis by reciprocal pairs for
discrete and continuous data.” *Educational and Psychological
measurement*, **26**(4), 825–831.

Sokal R, Michener C (1958). “A statistical method for evaluating
systematic relationships.” *University of Kansas Science Bulletin*,
**38**, 1409–1438.

Sneath, P. H., & Sokal, R. R. (1973). *Numerical taxonomy. The
principles and practice of numerical classification.*

## See also

[`hclust`](https://rdrr.io/r/stats/hclust.html),
[`dist.hamming`](https://klausvigo.github.io/phangorn/reference/dist.hamming.md),
[`NJ`](https://klausvigo.github.io/phangorn/reference/NJ.md),
[`as.phylo`](https://rdrr.io/pkg/ape/man/as.phylo.html),
[`fastme`](https://rdrr.io/pkg/ape/man/fastme.html),
[`nnls.tree`](https://klausvigo.github.io/phangorn/reference/designTree.md),
[`rtt`](https://rdrr.io/pkg/ape/man/rtt.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r

data(Laurasiatherian)
dm <- dist.ml(Laurasiatherian)
tree <- upgma(dm)
plot(tree)

```
