# Replace leading and trailing alignment gaps with an ambiguous state

Substitutes leading and trailing alignment gaps in aligned sequences
into N (i.e., A, C, G, or T) or ?. The gaps in the middle of the
sequences are left unchanged.

## Usage

``` r
latag2n.phyDat(x, amb = ifelse(attr(x, "type") == "DNA", "N", "?"),
  gap = "-", ...)
```

## Arguments

- x:

  an object of class `phyDat`.

- amb:

  character of the ambiguous state t replace the gaps.

- gap:

  gap parameter to replace.

- ...:

  Further arguments passed to or from other methods.

## Value

returns an object of class `phyDat`.

## See also

[`latag2n`](https://rdrr.io/pkg/ape/man/latag2n.html),
[`ancestral.pml`](https://klausvigo.github.io/phangorn/reference/ancestral.pml.md),
[`gap_as_state`](https://klausvigo.github.io/phangorn/reference/gap_as_state.md)

## Examples

``` r
x <- phyDat(matrix(c("-", "A", "G", "-", "T", "C"), 2, 3))
y <- latag2n.phyDat(x)
image(x)

image(y)
```
