# Utility function to plot.phylo

cladePar can help you coloring (choosing edge width/type) of clades.

## Usage

``` r
cladePar(tree, node, edge.color = "red", tip.color = edge.color,
  edge.width = 1, edge.lty = "solid", x = NULL, plot = FALSE, ...)
```

## Arguments

- tree:

  an object of class phylo.

- node:

  the node which is the common ancestor of the clade.

- edge.color:

  see plot.phylo.

- tip.color:

  see plot.phylo.

- edge.width:

  see plot.phylo.

- edge.lty:

  see plot.phylo.

- x:

  the result of a previous call to cladeInfo.

- plot:

  logical, if TRUE the tree is plotted.

- ...:

  Further arguments passed to or from other methods.

## Value

A list containing the information about the edges and tips.

## See also

[`plot.phylo`](https://rdrr.io/pkg/ape/man/plot.phylo.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
tree <- rtree(10)
plot(tree)
nodelabels()

x <- cladePar(tree, 12)
cladePar(tree, 18, "blue", "blue", x=x, plot=TRUE)

```
