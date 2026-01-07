# Compare splits and add support values to an object

Add support values to a `splits`, `phylo` or `networx` object.

## Usage

``` r
addConfidences(x, y, ...)

# S3 method for class 'phylo'
addConfidences(x, y, rooted = FALSE, ...)

presenceAbsence(x, y)

createLabel(x, y, label_y, type = "edge", nomatch = NA)
```

## Arguments

- x:

  an object of class `splits`, `phylo` or `networx`

- y:

  an object of class `splits`, `phylo`, `multiPhylo` or `networx`

- ...:

  Further arguments passed to or from other methods.

- rooted:

  logial, if FALSE bipartitions are considered, if TRUE clades.

- label_y:

  label of y matched on x. Will be usually of length(as.splits(x)).

- type:

  should labels returned for edges (in `networx`) or splits.

- nomatch:

  default value if no match between x and y is found.

## Value

The object `x` with added bootstrap / MCMC support values.

## References

Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W. (2017),
Intertwining phylogenetic trees and networks. *Methods Ecol Evol*.**8**,
1212–1220. doi:10.1111/2041-210X.12760

## See also

[`as.splits`](https://klausvigo.github.io/phangorn/reference/as.splits.md),
[`as.networx`](https://klausvigo.github.io/phangorn/reference/as.networx.md),
[`RF.dist`](https://klausvigo.github.io/phangorn/reference/treedist.md),
[`plot.phylo`](https://rdrr.io/pkg/ape/man/plot.phylo.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(woodmouse)
woodmouse <- phyDat(woodmouse)
tmpfile <- normalizePath(system.file(
             "extdata/trees/RAxML_bootstrap.woodmouse", package="phangorn"))
boot_trees <- read.tree(tmpfile)

dm <- dist.ml(woodmouse)
tree <- upgma(dm)
nnet <- neighborNet(dm)

tree <- addConfidences(tree, boot_trees)
nnet <- addConfidences(nnet, boot_trees)

plot(tree, show.node.label=TRUE)

plot(nnet, show.edge.label=TRUE)

```
