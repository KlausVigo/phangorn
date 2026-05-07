# anc_heatmap

Plots a heat map with ancestral states.

## Usage

``` r
anc_heatmap(x, y = NULL, use.edge.length = FALSE, align_label = TRUE,
  clade = NULL, select = NULL, cex.lab = 1, ...)
```

## Arguments

- x:

  an object of class ancestral or phyDat.

- y:

  an object of class phyDat.

- use.edge.length:

  a logical indicating whether to use the edge lengths of the phylogeny
  to draw the branches (the default) or not (if FALSE). This option has
  no effect if the object of class "phylo" has no ‘edge.length’ element.

- align_label:

  a logical value or an integer. If TRUE, the tips are aligned and
  dotted lines are drawn between the nodes of the tree and the labels.

- clade:

  a node number or label to extract the clade from the tree.

- select:

  a subset of characters, columns in the alignment.

- ...:

  Further arguments passed to or from other methods.

## Details

`anc_heatmap` plot the joint distribution or the most likely state.

## See also

[`anc_pml`](https://klausvigo.github.io/phangorn/reference/ancestral.pml.md),
[`plotAnc`](https://klausvigo.github.io/phangorn/reference/plot.ancestral.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(woodmouse)
wm <- as.phyDat(woodmouse)
tree <- pratchet(wm)
tree <- makeNodeLabel(tree)
anc_mp <- anc_pars(tree, wm)
op <- par(mar=c(4,1,4,5))
anc_heatmap(anc_mp)

par(op)
```
