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

Gray RD, Bryant D, Greenhill SJ (2010). “On the shape and fabric of
human history.” *Philosophical Transactions of the Royal Society B:
Biological Sciences*, **365**(1559), 3923.

Holland BR, Huber KT, Dress A, Moulton V (2002). “\\\delta\\ Plots: A
Tool for Analyzing Phylogenetic Distance Data.” *Molecular Biology and
Evolution*, **19**(12), 2051-2059. ISSN 0737-4038.
[doi:10.1093/oxfordjournals.molbev.a004030](https://doi.org/10.1093/oxfordjournals.molbev.a004030)
.

## See also

[`dist.hamming`](https://klausvigo.github.io/phangorn/reference/dist.hamming.md)

## Author

Alastair Potts and Klaus Schliep

## Examples

``` r

data(yeast)
hist(delta.score(yeast, "all"))

```
