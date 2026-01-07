# Distance Hadamard

Distance Hadamard produces spectra of splits from a distance matrix.

## Usage

``` r
distanceHadamard(dm, eps = 0.001)
```

## Arguments

- dm:

  A distance matrix.

- eps:

  Threshold value for splits.

## Value

`distanceHadamard` returns a matrix. The first column contains the
distance spectra, the second one the edge-spectra. If eps is positive an
object of with all splits greater eps is returned.

## References

Hendy, M. D. and Penny, D. (1993). Spectral Analysis of Phylogenetic
Data. *Journal of Classification*, **10**, 5-24.

## See also

[`hadamard`](https://klausvigo.github.io/phangorn/reference/hadamard.md),
[`lento`](https://klausvigo.github.io/phangorn/reference/lento.md),
[`plot.networx`](https://klausvigo.github.io/phangorn/reference/plot.networx.md),
[`neighborNet`](https://klausvigo.github.io/phangorn/reference/neighborNet.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>, Tim White

## Examples

``` r
data(yeast)
dm <- dist.hamming(yeast)
dm <- as.matrix(dm)
fit <- distanceHadamard(dm)
lento(fit)

plot(as.networx(fit))

```
