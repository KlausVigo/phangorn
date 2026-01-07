# Computes the \\\delta\\ score

Computes the treelikeness

## Usage

``` r
delta.score(x, arg = "mean", ...)
```

## Arguments

- x:

  an object of class `phyDat`

- arg:

  Specifies the return value, one of "all", "mean" or "sd"

- ...:

  further arguments passed through `dist.hamming`

## Value

A vector containing the \\\delta\\ scores.

## References

BR Holland, KT Huber, A Dress, V Moulton (2002) \\\delta\\ Plots: a tool
for analyzing phylogenetic distance data Russell D. Gray, David Bryant,
Simon J. Greenhill (2010) On the shape and fabric of human history
*Molecular Biology and Evolution*, **19(12)** 2051–2059

Russell D. Gray, David Bryant, Simon J. Greenhill (2010) On the shape
and fabric of human history *Phil. Trans. R. Soc. B*, **365** 3923–3933;
DOI: 10.1098/rstb.2010.0162

## See also

[`dist.hamming`](https://klausvigo.github.io/phangorn/reference/dist.hamming.md)

## Author

Alastair Potts and Klaus Schliep

## Examples

``` r
data(yeast)
hist(delta.score(yeast, "all"))

```
