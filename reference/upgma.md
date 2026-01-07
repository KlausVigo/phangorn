# UPGMA, WPGMA and sUPGMA

UPGMA and WPGMA clustering. UPGMA (Sokal and Michener 1958) and WPGMA
(McQuitty 1966) are a wrapper function around
[`hclust`](https://rdrr.io/r/stats/hclust.html) returning a `phylo`
object. `supgma` perform serial sampled UPGMA similar to Drummond and
Rodrigo (2000).

## Usage

``` r
upgma(D, method = "average", ...)

wpgma(D, method = "mcquitty", ...)

supgma(D, tip.dates, trace = 0, ...)
```

## Arguments

- D:

  A distance matrix, i.e. an object of class `dist`. If an matrix is
  supplied it is tried to covert it do a `dist` object.

- method:

  The agglomeration method to be used. This should be (an unambiguous
  abbreviation of) one of "ward", "single", "complete", "average",
  "mcquitty", "median" or "centroid". The default is "average" for UPGMA
  and "mcquitty" for WPGMA.

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

Sneath, P. H., & Sokal, R. R. (1973). *Numerical taxonomy. The
principles and practice of numerical classification.*

Sokal, R. R., & Michener, C. D. (1958). A statistical method for
evaluating systematic relationships. *University of Kansas Scientific
Bulletin*, v. 38.

Drummond, A., & Rodrigo, A. G. (2000). Reconstructing genealogies of
serial samples under the assumption of a molecular clock using
serial-sample UPGMA. *Molecular Biology and Evolution*, **17(12)**,
1807-1815.

McQuitty, L.L. (1966). Similarity Analysis by Reciprocal Pairs for
Discrete and Continuous Data. *Educational and Psychological
Measurement*, **26**, 825–831.

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
