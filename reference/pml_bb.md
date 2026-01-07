# Likelihood of a tree.

`pml_bb` for pml black box infers a phylogenetic tree infers a tree
using maximum likelihood (ML).

## Usage

``` r
pml_bb(x, model = NULL, rearrangement = ifelse(method == "unrooted",
  "stochastic", "ratchet"), method = "unrooted", start = NULL,
  tip.dates = NULL, ...)
```

## Arguments

- x:

  An alignment of class (either class `phyDat`, `DNAbin` or `AAbin`) or
  an object of class `modelTest`.

- model:

  A string providing model (e.g. "GTR+G(4)+I"). Not necessary if a
  modelTest object is supplied.

- rearrangement:

  Type of tree tree rearrangements to perform, one of "none", "NNI",
  "stochastic" or "ratchet"

- method:

  One of "unrooted", "ultrametric" or "tiplabeled".

- start:

  A starting tree can be supplied.

- tip.dates:

  A named vector of sampling times associated to the tips / sequences.

- ...:

  Further arguments passed to or from other methods.

## Value

`pml_bb` returns an object of class pml.

## Details

`pml_bb` is a convenience function combining `pml` and `optim.pml`. If
no tree is supplied, the function will generate a starting tree. If a
modelTest object is supplied the model will be chosen according to BIC.

`tip.dates` should be a named vector of sampling times, in any time
unit, with time increasing toward the present. For example, this may be
in units of “days since study start” or “years since 10,000 BCE”, but
not “millions of years ago”.

`model` takes a string and tries to extract the model. When an
`modelTest` object the best BIC model is chosen by default. The string
should contain a substitution model (e.g. JC, GTR, WAG) and can
additional have a term "+I" for invariant sites, "+G(4)" for a discrete
gamma model, "+R(4)" for a free rate model. In case of amino acid models
a term "+F" for estimating the amino acid frequencies. Whether
nucleotide frequencies are estimated is defined by
[`pml.control`](https://klausvigo.github.io/phangorn/reference/pml.control.md).

Currently very experimental and likely to change.

## See also

[`optim.pml`](https://klausvigo.github.io/phangorn/reference/pml.md),
[`modelTest`](https://klausvigo.github.io/phangorn/reference/modelTest.md),
[`rtt`](https://rdrr.io/pkg/ape/man/rtt.html),
[`pml.control`](https://klausvigo.github.io/phangorn/reference/pml.control.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(woodmouse)
tmp <- pml_bb(woodmouse, model="HKY+I", rearrangement="NNI")
#> optimize edge weights:  -1812.646 --> -1810.473 
#> optimize rate matrix:  -1810.473 --> -1758.757 
#> optimize invariant sites:  -1758.757 --> -1744.355 
#> optimize edge weights:  -1744.355 --> -1744.199 
#> optimize topology:  -1744.199 --> -1744.199  NNI moves:  0 
#> optimize rate matrix:  -1744.199 --> -1744.186 
#> optimize invariant sites:  -1744.186 --> -1744.186 
#> optimize edge weights:  -1744.186 --> -1744.186 
#> optimize rate matrix:  -1744.186 --> -1744.186 
#> optimize invariant sites:  -1744.186 --> -1744.186 
#> optimize edge weights:  -1744.186 --> -1744.186 

if (FALSE) { # \dontrun{
data(Laurasiatherian)
mt <- modelTest(Laurasiatherian)
fit <- pml_bb(mt)

# estimate free rate model with 2 rate categories
fit_HKY_R2 <- pml_bb(woodmouse, model="HKY+R(2)")
} # }
```
