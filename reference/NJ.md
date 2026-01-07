# Neighbor-Joining

This function performs the neighbor-joining tree estimation of Saitou
and Nei (1987). UNJ is the unweighted version from Gascuel (1997).

## Usage

``` r
NJ(x)

UNJ(x)
```

## Arguments

- x:

  A distance matrix.

## Value

an object of class `"phylo"`.

## Details

NJ is a wrapper around nj from ape.

## References

Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new method
for reconstructing phylogenetic trees. *Molecular Biology and
Evolution*, **4**, 406–425.

Studier, J. A and Keppler, K. J. (1988) A Note on the Neighbor-Joining
Algorithm of Saitou and Nei. *Molecular Biology and Evolution*, **6**,
729–731.

Gascuel, O. (1997) Concerning the NJ algorithm and its unweighted
version, UNJ. in Birkin et. al. *Mathematical Hierarchies and Biology*,
149–170.

## See also

[`nj`](https://rdrr.io/pkg/ape/man/nj.html),
[`dist.dna`](https://rdrr.io/pkg/ape/man/dist.dna.html),
[`dist.hamming`](https://klausvigo.github.io/phangorn/reference/dist.hamming.md),
[`upgma`](https://klausvigo.github.io/phangorn/reference/upgma.md),
[`fastme`](https://rdrr.io/pkg/ape/man/fastme.html)

## Author

Klaus P. Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(Laurasiatherian)
dm <- dist.ml(Laurasiatherian)
tree <- NJ(dm)
plot(tree)

```
