# Auxiliary for Controlling Fitting

Auxiliary functions for providing
[`optim.pml`](https://klausvigo.github.io/phangorn/reference/pml.md)`, `[`pml_bb`](https://klausvigo.github.io/phangorn/reference/pml_bb.md)
fitting. Use it to construct a `control` or `ratchet.par` argument.

## Usage

``` r
pml.control(epsilon = 1e-08, maxit = 10, trace = !getOption("quiet"),
  tau = 1e-08, statefreq = "empirical")

ratchet.control(iter = 20L, maxit = 200L, minit = 100L, prop = 1/2,
  rell = TRUE, bs = 100L)
```

## Arguments

- epsilon:

  Stop criterion for optimization (see details).

- maxit:

  Maximum number of iterations (see details).

- trace:

  Show output during optimization (see details).

- tau:

  minimal edge length.

- statefreq:

  take "empirical" or "estimate" state frequencies.

- iter:

  Number of iterations to stop if there is no change.

- minit:

  Minimum number of iterations.

- prop:

  Only used if `rearrangement=stochastic`. How many NNI moves should be
  added to the tree in proportion of the number of taxa.´

- rell:

  logical, if TRUE approximate bootstraping similar Minh et al. (2013)
  is performed.

- bs:

  number of approximate bootstrap samples.

## Value

A list with components named as the arguments for controlling the
fitting process.

## Details

`pml.control` controls the fitting process. `epsilon` and `maxit` are
only defined for the most outer loop, this affects `pmlCluster`,
`pmlPart` and `pmlMix`.

`epsilon` is not an absolute difference between, but instead is defined
as (logLik(k)-logLik(k+1))/logLik(k+1). This seems to be a good
compromise and to work reasonably well for small and large trees or
alignments.

If `trace` is set to zero than no out put is shown, if functions are
called internally than the trace is decreased by one, so a higher of
trace produces more feedback. It can be useful to figure out how long an
run will take and for debugging.

`statefreq` controls if base/state frequencies are optimized or
empirical estimates are taken, when this applies. For some nucleotide
models (e.g. JC, SYM) equal base frequencies and for amino acid models
precomputed state frequencies are used, if not '+F' is specified.

`tau` might be exactly zero if duplicated sequences in the alignment are
observed. In this case the analysis is performed only on unique
sequences and duplicated taxa are added to the tree with zero edge
length. This may lead to multifurcations if there are three or more
identical sequences. After optimization it is good practice to prune
away edges of length `tau` using `di2multi`. See also Janzen et al.
(2021).

## References

Minh, B. Q., Nguyen, M. A. T., & von Haeseler, A. (2013). Ultrafast
approximation for phylogenetic bootstrap. *Molecular biology and
evolution*, **30(5)**, 1188-1195.

Janzen, T., Bokma, F.,Etienne, R. S. (2021) Nucleotide Substitutions
during Speciation may Explain Substitution Rate Variation, *Systematic
Biology*, **71(5)**, 1244–1254.

## See also

[`pml_bb`](https://klausvigo.github.io/phangorn/reference/pml_bb.md),
[`optim.pml`](https://klausvigo.github.io/phangorn/reference/pml.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
pml.control()
#> $epsilon
#> [1] 1e-08
#> 
#> $maxit
#> [1] 10
#> 
#> $trace
#> [1] TRUE
#> 
#> $tau
#> [1] 1e-08
#> 
#> $statefreq
#> [1] "empirical"
#> 
pml.control(maxit=25)
#> $epsilon
#> [1] 1e-08
#> 
#> $maxit
#> [1] 25
#> 
#> $trace
#> [1] TRUE
#> 
#> $tau
#> [1] 1e-08
#> 
#> $statefreq
#> [1] "empirical"
#> 
```
