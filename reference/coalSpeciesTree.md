# Species Tree

`coalSpeciesTree` estimates species trees and can handle multiple
individuals per species.

## Usage

``` r
coalSpeciesTree(tree, X = NULL, sTree = NULL)
```

## Arguments

- tree:

  an object of class `multiPhylo`

- X:

  A `phyDat` object to define which individual belongs to which species.

- sTree:

  A species tree which fixes the topology.

## Value

The function returns an object of class `phylo`.

## Details

`coalSpeciesTree` estimates a single linkage tree as suggested by Liu et
al. (2010) from the element wise minima of the cophenetic matrices of
the gene trees. It extends `speciesTree` in ape as it allows that have
several individuals per gene tree.

## References

Liu, L., Yu, L. and Pearl, D. K. (2010) Maximum tree: a consistent
estimator of the species tree. *Journal of Mathematical Biology*,
**60**, 95–106.

## See also

[`speciesTree`](https://rdrr.io/pkg/ape/man/speciesTree.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com> Emmanuel Paradies

## Examples

``` r
## example in Liu et al. (2010)
tr1 <- read.tree(text = "(((B:0.05,C:0.05):0.01,D:0.06):0.04,A:0.1);")
tr2 <- read.tree(text = "(((A:0.07,C:0.07):0.02,D:0.09):0.03,B:0.12);")
TR <- c(tr1, tr2)
sp_tree <- coalSpeciesTree(TR)
```
