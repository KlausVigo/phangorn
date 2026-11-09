# Computes a neighborNet from a distance matrix

Computes a neighborNet, i.e. an object of class `networx` from a
distance matrix.

## Usage

``` r
neighborNet(x, ord = NULL)
```

## Arguments

- x:

  a distance matrix.

- ord:

  a circular ordering.

## Value

`neighborNet` returns an object of class networx.

## Details

`neighborNet` is still experimental. The cyclic ordering sometimes
differ from the SplitsTree implementation, the *ord* argument can be
used to enforce a certain circular ordering.

## References

Bryant D, Huson DH (2023). “NeighborNet: improved algorithms and
implementation.” *Frontiers in bioinformatics*, **3**, 1178600.
[doi:10.3389/fbinf.2023.1178600](https://doi.org/10.3389/fbinf.2023.1178600)
.

Bryant D, Moulton V (2004). “Neighbor-Net: An Agglomerative Method for
the Construction of Phylogenetic Networks.” *Molecular Biology and
Evolution*, **21**(2), 255-265.
[doi:10.1093/molbev/msh018](https://doi.org/10.1093/molbev/msh018) .

## See also

[`splitsNetwork`](https://klausvigo.github.io/phangorn/reference/splitsNetwork.md),
[`consensusNet`](https://klausvigo.github.io/phangorn/reference/consensusNet.md),
[`plot.networx`](https://klausvigo.github.io/phangorn/reference/plot.networx.md),
[`lento`](https://klausvigo.github.io/phangorn/reference/lento.md),
[`cophenetic.networx`](https://klausvigo.github.io/phangorn/reference/cophenetic.networx.md),
[`distanceHadamard`](https://klausvigo.github.io/phangorn/reference/distanceHadamard.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r

data(yeast)
dm <- dist.ml(yeast)
nnet <- neighborNet(dm)
plot(nnet)

plot(nnet, type="outline")

```
