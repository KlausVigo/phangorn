# Distances between trees

`treedist` computes different tree distance methods and `RF.dist` the
Robinson-Foulds or symmetric distance. The Robinson-Foulds distance only
depends on the topology of the trees. If edge weights should be
considered `wRF.dist` calculates the weighted RF distance (Robinson &
Foulds 1981). and `KF.dist` calculates the branch score distance (Kuhner
& Felsenstein 1994). `path.dist` computes the path difference metric as
described in Steel and Penny 1993). `sprdist` computes the approximate
SPR distance (Oliveira Martins et al. 2008, de Oliveira Martins 2016).

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

Functions like `RF.dist` returns the Robinson-Foulds distance (Robinson
and Foulds 1981) between either 2 trees or computes a matrix of all
pairwise distances if a `multiPhylo` object is given.

For large number of trees the distance functions can use a lot of
memory!

## References

de Oliveira Martins L., Leal E., Kishino H. (2008) *Phylogenetic
Detection of Recombination with a Bayesian Prior on the Distance between
Trees*. PLoS ONE **3(7)**. e2651. doi: 10.1371/journal.pone.0002651

de Oliveira Martins L., Mallo D., Posada D. (2016) *A Bayesian Supertree
Model for Genome-Wide Species Tree Reconstruction*. Syst. Biol.
**65(3)**: 397-416, doi:10.1093/sysbio/syu082

Steel M. A. and Penny P. (1993) *Distributions of tree comparison
metrics - some new results*, Syst. Biol., **42(2)**, 126–141

Kuhner, M. K. and Felsenstein, J. (1994) *A simulation comparison of
phylogeny algorithms under equal and unequal evolutionary rates*,
Molecular Biology and Evolution, **11(3)**, 459–468

D.F. Robinson and L.R. Foulds (1981) *Comparison of phylogenetic trees*,
Mathematical Biosciences, **53(1)**, 131–147

D.F. Robinson and L.R. Foulds (1979) Comparison of weighted labelled
trees. In Horadam, A. F. and Wallis, W. D. (Eds.), *Combinatorial
Mathematics VI: Proceedings of the Sixth Australian Conference on
Combinatorial Mathematics, Armidale, Australia*, 119–126

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
#> [1] 50
treedist(tree1, tree2)
#>      symmetric.difference   branch.score.difference           path.difference 
#>                 50.000000                  4.364967                161.291661 
#> quadratic.path.difference 
#>                 88.148561 
sprdist(tree1, tree2)
#>       spr spr_extra        rf     hdist 
#>         3         0        50        46 
trees <- rSPR(tree1, 1:5)
SPR.dist(tree1, trees)
#> [1] 1 2 3 4 5
```
