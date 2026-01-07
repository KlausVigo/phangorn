# Conversion among phylogenetic network objects

`as.networx` convert `splits` objects into a `networx` object. And most
important there exists a generic `plot` function to plot phylogenetic
network or split graphs.

## Usage

``` r
as.networx(x, ...)

# S3 method for class 'splits'
as.networx(x, planar = FALSE, coord = "none", ...)

# S3 method for class 'phylo'
as.networx(x, ...)
```

## Arguments

- x:

  an object of class `"splits"` or `"phylo"`

- ...:

  Further arguments passed to or from other methods.

- planar:

  logical whether to produce a planar graph from only cyclic splits (may
  excludes splits).

- coord:

  add coordinates of the nodes, allows to reproduce the plot.

## Value

an object of class `networx`.

## Details

A `networx` object hold the information for a phylogenetic network and
extends the `phylo` object. Therefore some generic function for `phylo`
objects will also work for `networx` objects. The argument
`planar = TRUE` will create a planar split graph based on a cyclic
ordering. These objects can be nicely plotted in `"2D"`. The argument
"coord" allows to create coordinates, with options are "none", "equal
angle", "3D", "2D" and "outline" (Bagci et al. 2021).

## Note

The internal representation is likely to change.

## References

Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W. (2017),
Intertwining phylogenetic trees and networks. *Methods Ecol Evol*.
**8**, 1212–1220. doi:10.1111/2041-210X.12760

Bagci, C., Bryant, D., Cetinkaya, B. and Huson, D.H. (2021), Microbial
Phylogenetic Context Using Phylogenetic Outlines. *Genome Biology and
Evolution*. Volume 13. Issue 9. evab213

## See also

[`consensusNet`](https://klausvigo.github.io/phangorn/reference/consensusNet.md),
[`neighborNet`](https://klausvigo.github.io/phangorn/reference/neighborNet.md),
[`splitsNetwork`](https://klausvigo.github.io/phangorn/reference/splitsNetwork.md),
[`hadamard`](https://klausvigo.github.io/phangorn/reference/hadamard.md),
[`distanceHadamard`](https://klausvigo.github.io/phangorn/reference/distanceHadamard.md),
[`plot.networx`](https://klausvigo.github.io/phangorn/reference/plot.networx.md),
[`evonet`](https://rdrr.io/pkg/ape/man/evonet.html),
[`as.phylo`](https://rdrr.io/pkg/ape/man/as.phylo.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
set.seed(1)
tree1 <- rtree(20, rooted=FALSE)
sp <- as.splits(rNNI(tree1, n=10))
net <- as.networx(sp)
plot(net)


spl <- allCircularSplits(5)
plot(as.networx(spl))

plot(as.networx(spl, coord="outline"))

if (FALSE) { # \dontrun{
# also see example in consensusNet
example(consensusNet)
} # }
```
