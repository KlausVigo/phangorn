# Pairwise Distances from Sequences

`dist.hamming`, `dist.ml` and `dist.logDet` compute pairwise distances
for an object of class `phyDat`. `dist.ml` uses DNA / AA sequences to
compute distances under different substitution models.

## Usage

``` r
dist.hamming(x, ratio = TRUE, exclude = "none")

dist.ml(x, model = "JC69", exclude = "none", bf = NULL, Q = NULL,
  k = 1L, shape = 1, ...)

dist.logDet(x)
```

## Arguments

- x:

  An object of class `phyDat`

- ratio:

  Compute uncorrected ('p') distance or character difference.

- exclude:

  One of "none", "all", "pairwise" indicating whether to delete the
  sites with gaps, missing data (or ambiguous states). See details
  below.

- model:

  One of "JC69", "F81" or one of 17 amino acid models see details.

- bf:

  A vector of base frequencies.

- Q:

  A vector containing the lower triangular part of the rate matrix.

- k:

  Number of intervals of the discrete gamma distribution.

- shape:

  Shape parameter of the gamma distribution.

- ...:

  Further arguments passed to or from other methods.

## Value

an object of class `dist`

## Details

So far 17 amino acid models are supported ("WAG", "JTT", "LG",
"Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV",
"HIVw", "HIVb", "FLU", "Blosum62", "Dayhoff_DCMut" and "JTT_DCMut") and
additional rate matrices and frequencies can be supplied.

The "F81" model uses empirical base frequencies, the "JC69" equal base
frequencies. This is even the case if the data are not nucleotides.

The argument `exclude` decides how gaps / ambiguous data / missing data
are treated. Usually gaps are treated as ambiguous states, but you can
give gaps its on state
[`gap_as_state`](https://klausvigo.github.io/phangorn/reference/gap_as_state.md).
`exclude="none"` keeps all ambiguous data. The behavior of `dist.ml` is
in this case these same you would achieve using `optim.pml` to compute
pairwise distances, it might be a bit odd. `exclude="all"` removes all
sites with ambiguous states and all gaps if these are coded as ambiguous
states. This can lead to the situation that there only few sites if any
fo the alignment left. Safer is therefore to use `exclude="pairwise"`
which only removes sites which are ambiguous for each pair of sequences.

## References

Lockhart, P. J., Steel, M. A., Hendy, M. D. and Penny, D. (1994)
Recovering evolutionary trees under a more realistic model of sequence
evolution. *Molecular Biology and Evolution*, **11**, 605–602.

Jukes TH and Cantor CR (1969). *Evolution of Protein Molecules*. New
York: Academic Press. 21–132.

McGuire, G., Prentice, M. J. and Wright, F. (1999). Improved error
bounds for genetic distances from DNA sequences. *Biometrics*, **55**,
1064–1070.

## See also

For more distance methods for nucleotide data see
[`dist.dna`](https://rdrr.io/pkg/ape/man/dist.dna.html) and
[`dist.p`](https://klausvigo.github.io/phangorn/reference/dist.p.md) for
pairwise polymorphism p-distances.
[`writeDist`](https://klausvigo.github.io/phangorn/reference/writeDist.md)
for export and import distances.

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(Laurasiatherian)
dm1 <- dist.hamming(Laurasiatherian)
tree1 <- NJ(dm1)
dm2 <- dist.logDet(Laurasiatherian)
tree2 <- NJ(dm2)
treedist(tree1,tree2)
#>      symmetric.difference   branch.score.difference           path.difference 
#>                4.00000000                0.05705091               30.95157508 
#> quadratic.path.difference 
#>                0.80097967 
# JC model
dm3 <- dist.ml(Laurasiatherian)
tree3 <- NJ(dm3)
treedist(tree1,tree3)
#>      symmetric.difference   branch.score.difference           path.difference 
#>                 6.0000000                 0.0412520                30.3644529 
#> quadratic.path.difference 
#>                 0.6106899 
# F81 + Gamma
dm4 <- dist.ml(Laurasiatherian, model="F81", k=4, shape=.4)
tree4 <- NJ(dm4)
treedist(tree1,tree4)
#>      symmetric.difference   branch.score.difference           path.difference 
#>                12.0000000                 0.1356107                40.7676342 
#> quadratic.path.difference 
#>                 2.0709714 
treedist(tree3,tree4)
#>      symmetric.difference   branch.score.difference           path.difference 
#>                8.00000000                0.09494752               39.52214569 
#> quadratic.path.difference 
#>                1.46345381 
```
