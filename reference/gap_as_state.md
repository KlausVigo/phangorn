# Treat gaps as a state

The function `gap_as_state` changes the contrast of an phyDat object to
treat as its own state. Internally `phyDat` are stored similar to a
`factor` objects and only the contrast matrix and some attributes
change.

## Usage

``` r
gap_as_state(obj, gap = "-", ambiguous = "?")

gap_as_ambiguous(obj, gap = "-")

has_gap_state(obj)
```

## Arguments

- obj:

  An object of class phyDat.

- gap:

  a character which codes for the gaps (default is "-").

- ambiguous:

  a character which codes for the ambiguous state

## Value

The functions return an object of class `phyDat`.

## See also

[`phyDat`](https://klausvigo.github.io/phangorn/reference/as.phyDat.md),
[`latag2n.phyDat`](https://klausvigo.github.io/phangorn/reference/latag2n.phyDat.md),
[`latag2n`](https://rdrr.io/pkg/ape/man/latag2n.html),
[`ancestral.pml`](https://klausvigo.github.io/phangorn/reference/ancestral.pml.md),
`gap_as_state`

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(Laurasiatherian)
tmp <- gap_as_state(Laurasiatherian)
contr <- attr(tmp, "contrast")
rownames(contr) <- attr(tmp, "allLevels")
contr
#>   a c g t -
#> a 1 0 0 0 0
#> c 0 1 0 0 0
#> g 0 0 1 0 0
#> t 0 0 0 1 0
#> u 0 0 0 1 0
#> m 1 1 0 0 0
#> r 1 0 1 0 0
#> w 1 0 0 1 0
#> s 0 1 1 0 0
#> y 0 1 0 1 0
#> k 0 0 1 1 0
#> v 1 1 1 0 0
#> h 1 1 0 1 0
#> d 1 0 1 1 0
#> b 0 1 1 1 0
#> n 1 1 1 1 0
#> ? 1 1 1 1 1
#> - 0 0 0 0 1
```
