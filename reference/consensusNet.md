# Computes a consensusNetwork from a list of trees

`consensusNet` computes a split network, i.e. an object of class
`networx` from a list of trees, i.e. an class of class `multiPhylo`
([Holland, Huber, Moulton, and Lockhart 2004](#reference++Holland2004))
.

## Usage

``` r
consensusNet(obj, prob = 0.3, ...)
```

## Arguments

- obj:

  An object of class multiPhylo.

- prob:

  the proportion a split has to be present in all trees to be
  represented in the network.

- ...:

  Further arguments passed to or from other methods.

## Value

`consensusNet` returns an object of class networx. This is just an
intermediate to plot phylogenetic networks with igraph.

## References

Holland BR, Huber KT, Moulton V, Lockhart PJ (2004). “Using Consensus
Networks to Visualize Contradictory Evidence for Species Phylogeny.”
*Molecular Biology and Evolution*, **21**(7), 1459-1461.
[doi:10.1093/molbev/msh145](https://doi.org/10.1093/molbev/msh145) .

## See also

[`splitsNetwork`](https://klausvigo.github.io/phangorn/reference/splitsNetwork.md),
[`neighborNet`](https://klausvigo.github.io/phangorn/reference/neighborNet.md),
[`lento`](https://klausvigo.github.io/phangorn/reference/lento.md),
[`distanceHadamard`](https://klausvigo.github.io/phangorn/reference/distanceHadamard.md),
[`plot.networx`](https://klausvigo.github.io/phangorn/reference/plot.networx.md),
[`maxCladeCred`](https://klausvigo.github.io/phangorn/reference/maxCladeCred.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r

data(Laurasiatherian)
set.seed(1)
bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)nj(dist.hamming(x)),
    bs=50)
cnet <- consensusNet(bs, 0.3)
plot(cnet, angle=-60)

if (FALSE) { # \dontrun{
library(rgl)
open3d()
plot(cnet, type = "3D", show.tip.label=FALSE, show.nodes=TRUE)
plot(cnet, type = "equal angle", show.edge.label=TRUE)

tmpfile <- normalizePath(system.file(
              "extdata/trees/RAxML_bootstrap.woodmouse", package="phangorn"))
trees <- read.tree(tmpfile)
cnet_woodmouse <- consensusNet(trees, 0.3)
plot(cnet_woodmouse, type = "equal angle", show.edge.label=TRUE)
} # }
```
