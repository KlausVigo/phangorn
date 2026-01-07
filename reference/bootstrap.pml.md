# Bootstrap

`bootstrap.pml` performs (non-parametric) bootstrap analysis and
`bootstrap.phyDat` produces a list of bootstrapped data sets. `plotBS`
plots a phylogenetic tree with the bootstrap values assigned to the
(internal) edges.

## Usage

``` r
bootstrap.pml(x, bs = 100, trees = TRUE, multicore = FALSE,
  mc.cores = NULL, tip.dates = NULL, ...)

bootstrap.phyDat(x, FUN, bs = 100, multicore = FALSE, mc.cores = NULL,
  jumble = TRUE, ...)
```

## Arguments

- x:

  an object of class `pml` or `phyDat`.

- bs:

  number of bootstrap samples.

- trees:

  return trees only (default) or whole `pml` objects.

- multicore:

  logical, whether models should estimated in parallel.

- mc.cores:

  The number of cores to use during bootstrap. Only supported on
  UNIX-alike systems.

- tip.dates:

  A named vector of sampling times associated to the tips/sequences.
  Leave empty if not estimating tip dated phylogenies.

- ...:

  further parameters used by `optim.pml` or `plot.phylo`.

- FUN:

  the function to estimate the trees.

- jumble:

  logical, jumble the order of the sequences.

## Value

`bootstrap.pml` returns an object of class `multi.phylo` or a list where
each element is an object of class `pml`. `plotBS` returns silently a
tree, i.e. an object of class `phylo` with the bootstrap values as node
labels. The argument `BStrees` is optional and if not supplied the tree
with labels supplied in the `node.label` slot.

## Details

It is possible that the bootstrap is performed in parallel, with help of
the multicore package. Unfortunately the multicore package does not work
under windows or with GUI interfaces ("aqua" on a mac). However it will
speed up nicely from the command line ("X11").

## References

Felsenstein J. (1985) Confidence limits on phylogenies. An approach
using the bootstrap. *Evolution* **39**, 783–791

Lemoine, F., Entfellner, J. B. D., Wilkinson, E., Correia, D., Felipe,
M. D., De Oliveira, T., & Gascuel, O. (2018). Renewing Felsenstein’s
phylogenetic bootstrap in the era of big data. *Nature*, **556(7702)**,
452–456.

Penny D. and Hendy M.D. (1985) Testing methods evolutionary tree
construction. *Cladistics* **1**, 266–278

Penny D. and Hendy M.D. (1986) Estimating the reliability of
evolutionary trees. *Molecular Biology and Evolution* **3**, 403–417

## See also

[`optim.pml`](https://klausvigo.github.io/phangorn/reference/pml.md),
[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md),
[`plot.phylo`](https://rdrr.io/pkg/ape/man/plot.phylo.html),
[`maxCladeCred`](https://klausvigo.github.io/phangorn/reference/maxCladeCred.md)
[`nodelabels`](https://rdrr.io/pkg/ape/man/nodelabels.html),[`consensusNet`](https://klausvigo.github.io/phangorn/reference/consensusNet.md)
and
[`SOWH.test`](https://klausvigo.github.io/phangorn/reference/SOWH.test.md)
for parametric bootstrap

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data(Laurasiatherian)
dm <- dist.hamming(Laurasiatherian)
tree <- NJ(dm)
# NJ
set.seed(123)
NJtrees <- bootstrap.phyDat(Laurasiatherian,
     FUN=function(x)NJ(dist.hamming(x)), bs=100)
treeNJ <- plotBS(tree, NJtrees, "phylogram")

# Maximum likelihood
fit <- pml(tree, Laurasiatherian)
fit <- optim.pml(fit, rearrangement="NNI")
set.seed(123)
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
treeBS <- plotBS(fit$tree,bs)

# Maximum parsimony
treeMP <- pratchet(Laurasiatherian)
treeMP <- acctran(treeMP, Laurasiatherian)
set.seed(123)
BStrees <- bootstrap.phyDat(Laurasiatherian, pratchet, bs = 100)
treeMP <- plotBS(treeMP, BStrees, "phylogram")
add.scale.bar()

# export tree with bootstrap values as node labels
# write.tree(treeBS)
} # }
```
