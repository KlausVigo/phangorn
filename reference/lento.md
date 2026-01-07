# Lento plot

The lento plot represents support and conflict of splits/bipartitions.

## Usage

``` r
lento(obj, xlim = NULL, ylim = NULL, main = "Lento plot", sub = NULL,
  xlab = NULL, ylab = NULL, bipart = TRUE, trivial = FALSE,
  col = rgb(0, 0, 0, 0.5), ...)
```

## Arguments

- obj:

  an object of class phylo, multiPhylo or splits

- xlim:

  graphical parameter

- ylim:

  graphical parameter

- main:

  graphical parameter

- sub:

  graphical parameter

- xlab:

  graphical parameter

- ylab:

  graphical parameter

- bipart:

  plot bipartition information.

- trivial:

  logical, whether to present trivial splits (default is FALSE).

- col:

  color for the splits / bipartition.

- ...:

  Further arguments passed to or from other methods.

## Value

lento returns a plot.

## References

Lento, G.M., Hickson, R.E., Chambers G.K., and Penny, D. (1995) Use of
spectral analysis to test hypotheses on the origin of pinninpeds.
*Molecular Biology and Evolution*, **12**, 28-52.

## See also

[`as.splits`](https://klausvigo.github.io/phangorn/reference/as.splits.md)`, `[`hadamard`](https://klausvigo.github.io/phangorn/reference/hadamard.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(yeast)
yeast.ry <- acgt2ry(yeast)
#> Warning: Found unknown characters (not supplied in levels). Deleted sites with unknown states.
splits.h <- h2st(yeast.ry)
lento(splits.h, trivial=TRUE)

```
