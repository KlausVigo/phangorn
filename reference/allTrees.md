# Compute all trees topologies.

`allTrees` computes all bifurcating tree topologies for rooted or
unrooted trees with up to 10 tips. The number of trees grows fast.

## Usage

``` r
allTrees(n, rooted = FALSE, tip.label = NULL)
```

## Arguments

- n:

  Number of tips (\<=10).

- rooted:

  Rooted or unrooted trees (default: rooted).

- tip.label:

  Tip labels.

## Value

an object of class `multiPhylo`.

## See also

[`rtree`](https://rdrr.io/pkg/ape/man/rtree.html),
[`nni`](https://klausvigo.github.io/phangorn/reference/nni.md),
[`howmanytrees`](https://rdrr.io/pkg/ape/man/howmanytrees.html),
[`dfactorial`](https://klausvigo.github.io/phangorn/reference/dfactorial.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
trees <- allTrees(5)

old.par <- par(no.readonly = TRUE)
par(mfrow = c(3,5))
for(i in 1:15)plot(trees[[i]])

par(old.par)
```
