# Discrete Gamma and Beta distribution

`discrete.gamma` internally used for the likelihood computations in
`pml` or `optim.pml`. The function `plot_gamma_plus_inv` and `plotRates`
are useful to understand how discrete rate models work and is used in
`plotRates`.

## Usage

``` r
discrete.gamma(alpha, k, w = NULL)

discrete.beta(shape1, shape2, k)

plot_gamma_plus_inv(w = NULL, g = NULL, shape = 1, inv = 0, k = 4,
  discrete = TRUE, cdf = TRUE, append = FALSE, xlab = "x",
  ylab = ifelse(cdf, "F(x)", "f(x)"), xlim = NULL, verticals = FALSE,
  edge.length = NULL, site.rate = "gamma", ...)
```

## Arguments

- alpha:

  Shape parameter of the gamma distribution.

- k:

  Number of intervals of the discrete gamma distribution.

- w:

  proportion of rate g.

- shape1, shape2:

  non-negative parameters of the Beta distribution.

- g:

  rates of discrete distribution.

- shape:

  Shape parameter of the gamma distribution.

- inv:

  Proportion of invariable sites.

- discrete:

  logical whether to plot discrete (default) or continuous pdf or cdf.

- cdf:

  logical whether to plot the cumulative distribution function or
  density / probability function.

- append:

  logical; if TRUE only add to an existing plot.

- xlab:

  a label for the x axis, defaults to a description of x.

- ylab:

  a label for the y axis, defaults to a description of y.

- xlim:

  the x limits of the plot.

- verticals:

  logical; if TRUE, draw vertical lines at steps.

- edge.length:

  Total edge length (sum of all edges in a tree).

- site.rate:

  Indicates what type of gamma distribution to use. Options are "gamma"
  (Yang 1994), "gamma_quadrature" using Laguerre quadrature approach of
  Felsenstein (2001) or "variable_gamma".

- ...:

  Further arguments passed to or from other methods.

## Value

`discrete.gamma` returns a matrix.

## See also

[`pml.fit`](https://klausvigo.github.io/phangorn/reference/pml.fit.md)`, `[`stepfun`](https://rdrr.io/r/stats/stepfun.html)`, `[`pgamma`](https://rdrr.io/r/stats/GammaDist.html)`, `[`pbeta`](https://rdrr.io/r/stats/Beta.html)`, `[`plotRates`](https://klausvigo.github.io/phangorn/reference/plotRates.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
discrete.gamma(1, 4)
#> [1] 0.1369538 0.4767519 1.0000000 2.3862944

old.par <- par(no.readonly = TRUE)
par(mfrow = c(2,1))
plot_gamma_plus_inv(shape=2, discrete = FALSE, cdf=FALSE)
plot_gamma_plus_inv(shape=2, append = TRUE, cdf=FALSE)

plot_gamma_plus_inv(shape=2, discrete = FALSE)
plot_gamma_plus_inv(shape=2, append = TRUE)

par(old.par)

data(Laurasiatherian)
fit <- pml_bb(Laurasiatherian, "JC+G(4)", rearrangement = "none")
plotRates(fit)
```
