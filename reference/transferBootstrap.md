# Transfer Bootstrap

`transferBootstrap` assigns transfer bootstrap (Lemoine et al. 2018)
values to the (internal) edges.

## Usage

``` r
transferBootstrap(tree, trees, phylo = TRUE, scale = TRUE)
```

## Arguments

- tree:

  The tree on which edges the bootstrap values are plotted.

- trees:

  a list of trees (object of class "multiPhylo").

- phylo:

  Logical, return a phylogentic tree with support value or a vector of
  bootstrap values.

- scale:

  scale the values.

## Value

a phylogentic tree (a phylo object) with bootstrap values assigned to
the node labels.

## References

Lemoine, F., Entfellner, J. B. D., Wilkinson, E., Correia, D., Felipe,
M. D., De Oliveira, T., & Gascuel, O. (2018). Renewing Felsenstein’s
phylogenetic bootstrap in the era of big data. *Nature*, **556(7702)**,
452–456.

## See also

[`plotBS`](https://klausvigo.github.io/phangorn/reference/plotBS.md),
[`maxCladeCred`](https://klausvigo.github.io/phangorn/reference/maxCladeCred.md),
[`drawSupportOnEdges`](https://rdrr.io/pkg/ape/man/plot.phyloExtra.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
fdir <- system.file("extdata/trees", package = "phangorn")
# RAxML best-known tree with bipartition support (from previous analysis)
raxml.tree <- read.tree(file.path(fdir,"RAxML_bipartitions.woodmouse"))
# RAxML bootstrap trees (from previous analysis)
raxml.bootstrap <- read.tree(file.path(fdir,"RAxML_bootstrap.woodmouse"))

tree_tbe <- transferBootstrap(raxml.tree,  raxml.bootstrap)
par(mfrow=c(1,2))
plotBS(tree_tbe)
# same as
plotBS(raxml.tree,  raxml.bootstrap, "p", "TBE")
```
