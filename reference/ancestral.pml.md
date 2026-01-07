# Ancestral character reconstruction.

Marginal reconstruction of the ancestral character states.

## Usage

``` r
ancestral.pml(object, type = "marginal", return = "prob", ...)

anc_pml(object, type = "marginal", ...)

ancestral.pars(tree, data, type = c("MPR", "ACCTRAN", "POSTORDER"),
  cost = NULL, return = "prob", ...)

anc_pars(tree, data, type = c("MPR", "ACCTRAN", "POSTORDER"), cost = NULL,
  ...)

pace(tree, data, type = c("MPR", "ACCTRAN", "POSTORDER"), cost = NULL,
  return = "prob", ...)
```

## Arguments

- object:

  an object of class pml

- type:

  method used to assign characters to internal nodes, see details.

- return:

  return a `phyDat` object or matrix of probabilities.

- ...:

  Further arguments passed to or from other methods.

- tree:

  a tree, i.e. an object of class pml

- data:

  an object of class phyDat

- cost:

  A cost matrix for the transitions between two states.

## Value

An object of class ancestral. This is a list containing the tree with
node labels, the original alignment as an `phyDat` object, a
`data.frame` containing the probabilities belonging to a state for all
(internal nodes) and the most likely state.

## Details

The argument "type" defines the criterion to assign the internal nodes.
For `ancestral.pml` so far "ml and marginal (empirical) "bayes" and for
`ancestral.pars` "MPR" and "ACCTRAN" are possible.

The function return a list containing the tree with node labels, the
original alignment as an `phyDat` object, a data.frame containing the
probabilities belonging to a state for all (internal nodes) and the most
likely state. For parsimony and nucleotide data the most likely state
might be ambiguous. For ML this is very unlikely to be the case.

If the input tree does not contain unique node labels the function
`ape::MakeNodeLabel` is used to create them.

With parsimony reconstruction one has to keep in mind that there will be
often no unique solution.

The functions use the node labels of the provided tree (also if part of
the `pml` object) if these are unique. Otherwise the function
`ape::MakeNodeLabel` is used to create them.

For further details see vignette("Ancestral").

## References

Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates,
Sunderland.

Swofford, D.L., Maddison, W.P. (1987) Reconstructing ancestral character
states under Wagner parsimony. *Math. Biosci.* **87**: 199–229

Yang, Z. (2006). *Computational Molecular evolution*. Oxford University
Press, Oxford.

## See also

[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md),
[`parsimony`](https://klausvigo.github.io/phangorn/reference/parsimony.md),
[`ace`](https://rdrr.io/pkg/ape/man/ace.html),
[`plotAnc`](https://klausvigo.github.io/phangorn/reference/plot.ancestral.md),
[`latag2n.phyDat`](https://klausvigo.github.io/phangorn/reference/latag2n.phyDat.md),
[`latag2n`](https://rdrr.io/pkg/ape/man/latag2n.html),
[`gap_as_state`](https://klausvigo.github.io/phangorn/reference/gap_as_state.md),
[`root`](https://rdrr.io/pkg/ape/man/root.html),
[`makeNodeLabel`](https://rdrr.io/pkg/ape/man/makeNodeLabel.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
example(NJ)
#> 
#> NJ> data(Laurasiatherian)
#> 
#> NJ> dm <- dist.ml(Laurasiatherian)
#> 
#> NJ> tree <- NJ(dm)
#> 
#> NJ> plot(tree)

# generate node labels to ensure plotting will work
tree <- makeNodeLabel(tree)
fit <- pml(tree, Laurasiatherian)
anc.ml <- anc_pml(fit)
anc.p <- anc_pars(tree, Laurasiatherian)
# plot ancestral sequences at the root
plotSeqLogo( anc.ml, 48, 1, 20)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the ggseqlogo package.
#>   Please report the issue at <https://github.com/omarwagih/ggseqlogo/issues>.

plotSeqLogo( anc.p, 48, 1, 20)

# plot the first character
plotAnc(anc.ml)

# plot the third character
plotAnc(anc.ml, 3)

```
