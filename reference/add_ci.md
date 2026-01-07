# Draw Confidences Intervals on Phylogenies

These are low-level plotting commands to draw the confidence intervals
on the node of a tree as rectangles with coloured backgrounds or add
boxplots to ultrametric or tipdated trees.

## Usage

``` r
add_ci(tree, trees, col95 = "#FF00004D", col50 = "#0000FF4D",
  height = 0.7, legend = TRUE, ...)

add_boxplot(tree, trees, ...)
```

## Arguments

- tree:

  a phylogenetic tree to which the confidences should be added.

- trees:

  phylogenetic trees, i.e. an object of class `multiPhylo`

- col95:

  colour used for the 95% intervals; by default: transparent red.

- col50:

  colour used for the 50% intervals; by default: transparent blue.

- height:

  the height of the boxes.

- legend:

  a logical value.

- ...:

  arguments passed to other functions,
  [`legend`](https://rdrr.io/r/graphics/legend.html) or
  [`bxp`](https://rdrr.io/r/graphics/bxp.html).

## Value

`add_ci` and `add_boxplot` return silently the tree object.

## Details

All trees should to be rooted, either ultrametric or tip dated.

## See also

[`plot.phylo`](https://rdrr.io/pkg/ape/man/plot.phylo.html),
[`plotBS`](https://klausvigo.github.io/phangorn/reference/plotBS.md),
[`add_edge_length`](https://klausvigo.github.io/phangorn/reference/add_edge_length.md),
[`maxCladeCred`](https://klausvigo.github.io/phangorn/reference/maxCladeCred.md)

## Author

Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep

## Examples

``` r
data("Laurasiatherian")
dm <- dist.hamming(Laurasiatherian)
tree <- upgma(dm)
set.seed(123)
trees <- bootstrap.phyDat(Laurasiatherian,
                          FUN=function(x)upgma(dist.hamming(x)), bs=100)
tree <- plotBS(tree, trees, "phylogram")
add_ci(tree, trees, bty="n")

plot(tree, direction="downwards")
add_boxplot(tree, trees, boxwex=.7)
```
