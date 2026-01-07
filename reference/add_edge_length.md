# Assign and compute edge lengths from a sample of trees

This command can infer some average edge lengths and assign them from a
(bootstrap/MCMC) sample.

## Usage

``` r
add_edge_length(tree, trees, fun = function(x) median(na.omit(x)),
  rooted = all(is.rooted(trees)))
```

## Arguments

- tree:

  a phylogenetic tree or splitnetwork where edge lengths are assigned
  to.

- trees:

  an object of class multiPhylo, where the average for the edges is
  computed from.

- fun:

  a function to compute the average (default is median).

- rooted:

  rooted logical, if FALSE edge lengths is a function of the observed
  splits, if TRUE edge lengths are estimated from height for the
  observed clades.

## Value

The tree with newly assigned edge length.

## See also

[`node.depth`](https://rdrr.io/pkg/ape/man/node.depth.html),
[`consensus`](https://rdrr.io/pkg/ape/man/consensus.html),
[`maxCladeCred`](https://klausvigo.github.io/phangorn/reference/maxCladeCred.md),
[`add_boxplot`](https://klausvigo.github.io/phangorn/reference/add_ci.md)

## Author

Klaus Schliep

## Examples

``` r
data("Laurasiatherian")
set.seed(123)
bs <- bootstrap.phyDat(Laurasiatherian,
                FUN=function(x)upgma(dist.ml(x)), bs=100)
tree_compat <- allCompat(bs, rooted=TRUE) |>
              add_edge_length(bs)
plot(tree_compat)
add_boxplot(tree_compat, bs, boxwex=.7)
```
