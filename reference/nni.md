# Tree rearrangements.

`nni` returns a list of all trees which are one nearest neighbor
interchange away. `rNNI` and `rSPR` are two methods which simulate
random trees which are a specified number of rearrangement apart from
the input tree. Both methods assume that the input tree is bifurcating.
These methods may be useful in simulation studies.

## Usage

``` r
nni(tree)

rNNI(tree, moves = 1, n = length(moves))

rSPR(tree, moves = 1, n = length(moves), k = NULL)
```

## Arguments

- tree:

  A phylogenetic `tree`, object of class `phylo`.

- moves:

  Number of tree rearrangements to be transformed on a tree. Can be a
  vector

- n:

  Number of trees to be simulated.

- k:

  If defined just SPR of distance k are performed.

## Value

an object of class multiPhylo.

## See also

[`allTrees`](https://klausvigo.github.io/phangorn/reference/allTrees.md),
[`SPR.dist`](https://klausvigo.github.io/phangorn/reference/treedist.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
tree <- rtree(20, rooted = FALSE)
trees1 <- nni(tree)
trees2 <- rSPR(tree, 2, 10)
```
