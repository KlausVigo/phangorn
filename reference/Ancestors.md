# tree utility function

Functions for describing relationships among phylogenetic nodes.

## Usage

``` r
Ancestors(x, node, type = c("all", "parent"))

allDescendants(x)

Children(x, node)

Descendants(x, node, type = c("tips", "children", "all"))

Siblings(x, node, include.self = FALSE)

mrca.phylo(x, node = NULL, full = FALSE)
```

## Arguments

- x:

  a tree (a phylo object).

- node:

  an integer or character vector (or scalar) corresponding to a node ID

- type:

  specify whether to return just direct children / parents or all

- include.self:

  whether to include self in list of siblings

- full:

  a logical indicating whether to return the MRCAs among all tips and
  nodes (if TRUE); the default is to return only the MRCAs among tips.

## Value

a vector or a list containing the indices of the nodes.

## Details

These functions are inspired by `treewalk` in phylobase package, but
work on the S3 `phylo` objects. The nodes are the indices as given in
edge matrix of an phylo object. From taxon labels these indices can be
easily derived matching against the `tip.label` argument of an phylo
object, see example below. All the functions allow `node` to be either a
scalar or vector. `mrca` is a faster version of the mrca in ape, in
phangorn only because of dependencies. If the argument node is missing
the function is evaluated for all nodes.

## Functions

- `allDescendants()`: list all the descendant nodes of each node

## See also

`treewalk`, [`as.phylo`](https://rdrr.io/pkg/ape/man/as.phylo.html),
[`nodelabels`](https://rdrr.io/pkg/ape/man/nodelabels.html)

## Examples

``` r
tree <- rtree(10)
plot(tree, show.tip.label = FALSE)
nodelabels()
tiplabels()

Ancestors(tree, 1:3, "all")
#> [[1]]
#> [1] 14 13 12 11
#> 
#> [[2]]
#> [1] 15 14 13 12 11
#> 
#> [[3]]
#> [1] 15 14 13 12 11
#> 
Children(tree, 11)
#> [1] 12 18
Descendants(tree, 11, "tips")
#> [[1]]
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
Siblings(tree, 3)
#> [1] 2
# Siblings of all nodes
Siblings(tree)
#> [[1]]
#> [1] 15
#> 
#> [[2]]
#> [1] 3
#> 
#> [[3]]
#> [1] 2
#> 
#> [[4]]
#> [1] 14
#> 
#> [[5]]
#> [1] 6
#> 
#> [[6]]
#> [1] 5
#> 
#> [[7]]
#> [1] 17
#> 
#> [[8]]
#> [1] 19
#> 
#> [[9]]
#> [1] 10
#> 
#> [[10]]
#> [1] 9
#> 
#> [[11]]
#> NULL
#> 
#> [[12]]
#> [1] 18
#> 
#> [[13]]
#> [1] 16
#> 
#> [[14]]
#> [1] 4
#> 
#> [[15]]
#> [1] 1
#> 
#> [[16]]
#> [1] 13
#> 
#> [[17]]
#> [1] 7
#> 
#> [[18]]
#> [1] 12
#> 
#> [[19]]
#> [1] 8
#> 
mrca.phylo(tree, 1:3)
#> [1] 14
mrca.phylo(tree, match(c("t1", "t2", "t3"), tree$tip))
#> [1] 11
mrca.phylo(tree)
#>     t4 t6 t1 t5 t7 t9 t2 t10 t8 t3
#> t4   1 14 14 13 12 12 12  11 11 11
#> t6  14  2 15 13 12 12 12  11 11 11
#> t1  14 15  3 13 12 12 12  11 11 11
#> t5  13 13 13  4 12 12 12  11 11 11
#> t7  12 12 12 12  5 17 16  11 11 11
#> t9  12 12 12 12 17  6 16  11 11 11
#> t2  12 12 12 12 16 16  7  11 11 11
#> t10 11 11 11 11 11 11 11   8 18 18
#> t8  11 11 11 11 11 11 11  18  9 19
#> t3  11 11 11 11 11 11 11  18 19 10
# same as mrca(tree), but faster for large trees
```
