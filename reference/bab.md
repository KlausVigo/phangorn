# Branch and bound for finding all most parsimonious trees

`bab` finds all most parsimonious trees.

## Usage

``` r
bab(data, tree = NULL, trace = 1, ...)
```

## Arguments

- data:

  an object of class phyDat.

- tree:

  a phylogenetic tree an object of class phylo, otherwise a pratchet
  search is performed.

- trace:

  defines how much information is printed during optimization.

- ...:

  Further arguments passed to or from other methods

## Value

`bab` returns all most parsimonious trees in an object of class
`multiPhylo`.

## Details

This implementation is very slow and depending on the data may take very
long time. In the worst case all \\(2n-5)!! = 1 \times 3 \times 5 \times
\ldots \times (2n-5)\\ possible trees have to be examined, where n is
the number of species / tips. For ten species there are already 2027025
tip-labelled unrooted trees. It only uses some basic strategies to find
a lower and upper bounds similar to penny from phylip. `bab` uses a very
basic heuristic approach of MinMax Squeeze (Holland et al. 2005) to
improve the lower bound.  
`bab` might return multifurcating trees. These multifurcations could be
resolved in all ways.  
On the positive side `bab` is not like many other implementations
restricted to binary or nucleotide data.

## References

Hendy, M.D. and Penny D. (1982) Branch and bound algorithms to determine
minimal evolutionary trees. *Math. Biosc.* **59**, 277-290

Holland, B.R., Huber, K.T. Penny, D. and Moulton, V. (2005) The MinMax
Squeeze: Guaranteeing a Minimal Tree for Population Data, *Molecular
Biology and Evolution*, **22**, 235–242

White, W.T. and Holland, B.R. (2011) Faster exact maximum parsimony
search with XMP. *Bioinformatics*, **27(10)**,1359–1367

## See also

[`pratchet`](https://klausvigo.github.io/phangorn/reference/parsimony.md),
[`dfactorial`](https://klausvigo.github.io/phangorn/reference/dfactorial.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com> based on work on Liam Revell

## Examples

``` r
data(yeast)
dfactorial(11)
#> [1] 10395
# choose only the first two genes
gene12 <- yeast[, 1:3158]
trees <- bab(gene12)
#> Compute starting tree
#> lower bound: 2798 
#> upper bound: 3045 
#> Search Baumraum (tree space)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |======================================================================| 100%
```
