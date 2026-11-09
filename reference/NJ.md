# Neighbor-Joining

This function performs the neighbor-joining tree estimation of [Saitou
and Nei (1987)](#reference++Saitou1987); [Studier and Keppler
(1988)](#reference++Studier1988) . UNJ is the unweighted version from
[Gascuel (1997)](#reference++Gascuel1997) .

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

Gascuel O (1997). “Concerning the NJ algorithm and its unweighted
version UNJ. 21 Pages 149-171 in: Mathematical hierarchies and biology.
B. Mirkin, FR McMorris, 22 F. Roberts and A. Rzhetsky, eds.” *American
Mathematical Society, Providence, RI*, **23**.

Saitou N, Nei M (1987). “The Neighbor-Joining Method - a New Method for
Reconstructing Phylogenetic Trees.” *Molecular Biology and Evolution*,
**4**(4), 406–425.

Studier JA, Keppler KJ (1988). “A Note on the Neighbor-Joining Algorithm
of Saitou and Nei.” *Molecular Biology and Evolution*, **5**(6),
729–731.

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
