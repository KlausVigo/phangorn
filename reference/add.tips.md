# Add tips to a tree

This function binds tips to nodes of a phylogenetic trees.

## Usage

``` r
add.tips(tree, tips, where, edge.length = NULL)
```

## Arguments

- tree:

  an object of class "phylo".

- tips:

  a character vector containing the names of the tips.

- where:

  an integer or character vector of the same length as tips giving the
  number of the node or tip of the tree where to add the new tips.

- edge.length:

  optional numeric vector with edge length

## Value

an object of class phylo

## See also

[`bind.tree`](https://rdrr.io/pkg/ape/man/bind.tree.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
tree <- rcoal(10)
plot(tree)
nodelabels()
tiplabels()

tree1 <- add.tips(tree, c("A", "B", "C"), c(1,2,15))
plot(tree1)
```
