# Internal maximum likelihood functions.

These functions are internally used for the likelihood computations in
`pml` or `optim.pml`.

## Usage

``` r
lli(data, tree = NULL, ...)

pml.fit(tree, data, bf = rep(1/length(levels), length(levels)), shape = 1,
  k = 1, Q = rep(1, length(levels) * (length(levels) - 1)/2),
  levels = attr(data, "levels"), inv = 0, rate = 1, g = NULL,
  w = NULL, eig = NULL, INV = NULL, ll.0 = NULL, llMix = NULL,
  wMix = 0, ..., site = FALSE, ASC = FALSE, site.rate = "gamma")
```

## Arguments

- data:

  An alignment, object of class `phyDat`.

- tree:

  A phylogenetic `tree`, object of class `phylo`.

- ...:

  Further arguments passed to or from other methods.

- bf:

  Base frequencies.

- shape:

  Shape parameter of the gamma distribution.

- k:

  Number of intervals of the discrete gamma distribution.

- Q:

  A vector containing the lower triangular part of the rate matrix.

- levels:

  The alphabet used e.g. c("a", "c", "g", "t") for DNA

- inv:

  Proportion of invariable sites.

- rate:

  Rate.

- g:

  vector of quantiles (default is NULL)

- w:

  vector of probabilities (default is NULL)

- eig:

  Eigenvalue decomposition of Q

- INV:

  Sparse representation of invariant sites

- ll.0:

  default is NULL

- llMix:

  default is NULL

- wMix:

  default is NULL

- site:

  return the log-likelihood or vector of sitewise likelihood values

- ASC:

  ascertainment bias correction (ASC), allows to estimate models like
  Lewis' Mkv.

- site.rate:

  Indicates what type of gamma distribution to use. Options are "gamma"
  approach of Yang 1994 (default), "gamma_quadrature" after the Laguerre
  quadrature approach of Felsenstein 2001 and "freerate".

## Value

`pml.fit` returns the log-likelihood.

## Details

These functions are exported to be used in different packages so far
only in the package coalescentMCMC, but are not intended for end user.
Most of the functions call C code and are far less forgiving if the
import is not what they expect than `pml`.

## References

Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maximum
likelihood approach. *Journal of Molecular Evolution*, **17**, 368–376.

## See also

[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md)`, `[`pml_bb`](https://klausvigo.github.io/phangorn/reference/pml_bb.md)`, `[`pmlPart`](https://klausvigo.github.io/phangorn/reference/pmlPart.md)`, `[`pmlMix`](https://klausvigo.github.io/phangorn/reference/pmlMix.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(Laurasiatherian)
tree <- NJ(dist.ml(Laurasiatherian))
bf <- rep(0.25, 4)
eig <- edQt()
pml.init(Laurasiatherian)
#> NULL
pml.fit(tree, Laurasiatherian, bf=bf, eig=eig)
#> [1] -54808.83
pml.free()
#> NULL
pml(tree, Laurasiatherian) |> logLik()
#> 'log Lik.' -54808.83 (df=91)
```
