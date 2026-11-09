# Distances between trees

`treedist` computes different tree distance methods and `RF.dist` the
Robinson-Foulds or symmetric distance. The Robinson-Foulds distance
([Robinson and Foulds 1979](#reference++Robinson1979); [Robinson and
Foulds 1981](#reference++Robinson1981)) only depends on the topology of
the trees. If edge weights should be considered `wRF.dist` calculates
the weighted RF distance ([Robinson and Foulds
1981](#reference++Robinson1981)) and `KF.dist` calculates the branch
score distance ([Kuhner and Felsenstein 1994](#reference++Kuhner1994)) .
`path.dist` computes the path difference metric as described in [Steel
and Penny (1993)](#reference++Steel1993) . `sprdist` computes the
approximate SPR distance ([de Oliveira Martins, Leal, and Kishino
2008](#reference++Leo2008); [De Oliveira Martins, Mallo, and Posada
2016](#reference++Leo2016)) .

## Usage

``` r
treedist(tree1, tree2, check.labels = TRUE)

sprdist(tree1, tree2)

SPR.dist(tree1, tree2 = NULL)

RF.dist(tree1, tree2 = NULL, normalize = FALSE, check.labels = TRUE,
  rooted = FALSE)

wRF.dist(tree1, tree2 = NULL, normalize = FALSE, check.labels = TRUE,
  rooted = FALSE)

KF.dist(tree1, tree2 = NULL, check.labels = TRUE, rooted = FALSE)

path.dist(tree1, tree2 = NULL, check.labels = TRUE, use.weight = FALSE)
```

## Arguments

- tree1:

  A phylogenetic tree (class `phylo`) or vector of trees (an object of
  class `multiPhylo`). See details

- tree2:

  A phylogenetic tree.

- check.labels:

  compares labels of the trees.

- normalize:

  compute normalized RF-distance, see details.

- rooted:

  take bipartitions for rooted trees into account, default is unrooting
  the trees.

- use.weight:

  use edge.length argument or just count number of edges on the path
  (default)

## Value

`treedist` returns a vector containing the following tree distance
methods

- symmetric.difference:

  symmetric.difference or Robinson-Foulds distance

- branch.score.difference:

  branch.score.difference

- path.difference:

  path.difference

- weighted.path.difference:

  weighted.path.difference

## Details

The Robinson-Foulds distance between two trees \\T_1\\ and \\T_2\\ with
\\n\\ tips is defined as (following the notation Steel and Penny 1993):
\$\$d(T_1, T_2) = i(T_1) + i(T_2) - 2v_s(T_1, T_2)\$\$ where \\i(T_1)\\
denotes the number of internal edges and \\v_s(T_1, T_2)\\ denotes the
number of internal splits shared by the two trees. The normalized
Robinson-Foulds distance is derived by dividing \\d(T_1, T_2)\\ by the
maximal possible distance \\i(T_1) + i(T_2)\\. If both trees are
unrooted and binary this value is \\2n-6\\.

Functions like `RF.dist` returns the Robinson-Foulds distance between
either 2 trees or computes a matrix of all pairwise distances if a
`multiPhylo` object is given.

For large number of trees the distance functions can use a lot of
memory!

## References

de Oliveira Martins L, Leal É, Kishino H (2008). “Phylogenetic Detection
of Recombination with a Bayesian Prior on the Distance between Trees.”
*PLoS ONE*, **3**(7), 1-13.
[doi:10.1371/journal.pone.0002651](https://doi.org/10.1371/journal.pone.0002651)
.

De Oliveira Martins L, Mallo D, Posada D (2016). “A Bayesian Supertree
Model for Genome-Wide Species Tree Reconstruction.” *Systematic
Biology*, **65**(3), 397-416.
[doi:10.1093/sysbio/syu082](https://doi.org/10.1093/sysbio/syu082) .

Kuhner MK, Felsenstein J (1994). “A simulation comparison of phylogeny
algorithms under equal and unequal evolutionary rates.” *Molecular
Biology and Evolution*, **11**(3), 459-468.
[doi:10.1093/oxfordjournals.molbev.a040126](https://doi.org/10.1093/oxfordjournals.molbev.a040126)
.

Robinson D, Foulds L (1981). “Comparison of phylogenetic trees.”
*Mathematical Biosciences*, **53**(1), 131 - 147. ISSN 0025-5564.
[doi:10.1016/0025-5564(81)90043-2](https://doi.org/10.1016/0025-5564%2881%2990043-2)
.

Robinson DF, Foulds LR (1979). “Comparison of weighted labelled trees.”
In Horadam AF, Wallis WD (eds.), *Combinatorial Mathematics VI:
Proceedings of the Sixth Australian Conference on Combinatorial
Mathematics, Armidale, Australia, August 1978*, 119–126. Springer Berlin
Heidelberg, Berlin, Heidelberg. ISBN 978-3-540-34857-3.
[doi:10.1007/BFb0102690](https://doi.org/10.1007/BFb0102690) .

Steel MA, Penny D (1993). “Distributions of Tree Comparison Metrics—Some
New Results.” *Systematic Biology*, **42**(2), 126-141. ISSN 1063-5157.
[doi:10.1093/sysbio/42.2.126](https://doi.org/10.1093/sysbio/42.2.126) .

## See also

[`dist.topo`](https://rdrr.io/pkg/ape/man/dist.topo.html),
[`nni`](https://klausvigo.github.io/phangorn/reference/nni.md),
[`superTree`](https://klausvigo.github.io/phangorn/reference/superTree.md),
[`mast`](https://klausvigo.github.io/phangorn/reference/mast.md)

## Author

Klaus P. Schliep <klaus.schliep@gmail.com>, Leonardo de Oliveira Martins

## Examples

``` r

tree1 <- rtree(100, rooted=FALSE)
tree2 <- rSPR(tree1, 3)
RF.dist(tree1, tree2)
#> [1] 52
treedist(tree1, tree2)
#>      symmetric.difference   branch.score.difference           path.difference 
#>                 52.000000                  4.000749                160.483644 
#> quadratic.path.difference 
#>                 70.116070 
sprdist(tree1, tree2)
#>       spr spr_extra        rf     hdist 
#>         3         0        52        34 
trees <- rSPR(tree1, 1:5)
SPR.dist(tree1, trees)
#> [1] 1 2 2 4 5
```
