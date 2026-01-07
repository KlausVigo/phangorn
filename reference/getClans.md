# Clans, slices and clips

Functions for clanistics to compute clans, slices, clips for unrooted
trees and functions to quantify the fragmentation of trees.

## Usage

``` r
getClans(tree)

getSlices(tree)

getClips(tree, all = TRUE)

getDiversity(tree, x, norm = TRUE, var.names = NULL, labels = "new")

# S3 method for class 'clanistics'
summary(object, ...)

diversity(tree, X)
```

## Arguments

- tree:

  An object of class phylo or multiPhylo (getDiversity).

- all:

  A logical, return all or just the largest clip.

- x:

  An object of class phyDat.

- norm:

  A logical, return Equitability Index (default) or Shannon Diversity.

- var.names:

  A vector of variable names.

- labels:

  see details.

- object:

  an object for which a summary is desired.

- ...:

  Further arguments passed to or from other methods.

- X:

  a data.frame

## Value

getClans, getSlices and getClips return a matrix of partitions, a matrix
of ones and zeros where rows correspond to a clan, slice or clip and
columns to tips. A one indicates that a tip belongs to a certain
partition.  
getDiversity returns a list with tree object, the first is a data.frame
of the equitability index or Shannon divergence and parsimony scores
(p-score) for all trees and variables. The data.frame has two
attributes, the first is a splits object to identify the taxa of each
tree and the second is a splits object containing all partitions that
perfectly fit.

## Details

Every split in an unrooted tree defines two complementary clans. Thus
for an unrooted binary tree with \\n\\ leaves there are \\2n - 3\\
edges, and therefore \\4n - 6\\ clans (including \\n\\ trivial clans
containing only one leave).

Slices are defined by a pair of splits or tripartitions, which are not
clans. The number of distinguishable slices for a binary tree with \\n\\
tips is \\2n^2 - 10n + 12\\.

A clip is a different type of partition, defining groups of leaves that
are related in terms of evolutionary distances and not only topology.
Namely, clips are groups of leaves for which all pairwise path-length
distances are smaller than a given threshold value (Lapointe et al.
2010). There exists different numbers of clips for different thresholds,
the largest (and trivial) one being the whole tree. There is always a
clip containing only the two leaves with the smallest pairwise distance.

Clans, slices and clips can be used to characterize how well a vector of
categorial characters (natives/intruders) fit on a tree. We will follow
the definitions of Lapointe et al.(2010). A complete clan is a clan that
contains all leaves of a given state/color, but can also contain leaves
of another state/color. A clan is homogeneous if it only contains leaves
of one state/color.

`getDiversity` computes either the  
Shannon Diversity: \\H = -\sum\_{i=1}^{k}(N_i/N) log(N_i/N),
N=\sum\_{i=1}^{k} N_i\\  
or the  
Equitability Index: \\E = H / log(N)\\  
where \\N_i\\ are the sizes of the \\k\\ largest homogeneous clans of
intruders. If the categories of the data can be separated by an edge of
the tree then the E-value will be zero, and maximum equitability (E=1)
is reached if all intruders are in separate clans. getDiversity computes
these Intruder indices for the whole tree, complete clans and complete
slices. Additionally the parsimony scores (p-scores) are reported. The
p-score indicates if the leaves contain only one color (p-score=0), if
the the leaves can be separated by a single split (perfect clan,
p-score=1) or by a pair of splits (perfect slice, p-score=2).

So far only 2 states are supported (native, intruder), however it is
also possible to recode several states into the native or intruder state
using contrasts, for details see section 2 in
vignette("phangorn-specials"). Furthermore unknown character states are
coded as ambiguous character, which can act either as native or intruder
minimizing the number of clans or changes (in parsimony analysis) needed
to describe a tree for given data.

Set attribute labels to "old" for analysis as in Schliep et al. (2010)
or to "new" for names which are more intuitive.

`diversity` returns a data.frame with the parsimony score for each tree
and each levels of the variables in `X`. `X` has to be a `data.frame`
where each column is a factor and the rownames of `X` correspond to the
tips of the trees.

## References

Lapointe, F.-J., Lopez, P., Boucher, Y., Koenig, J., Bapteste, E. (2010)
Clanistics: a multi-level perspective for harvesting unrooted gene
trees. *Trends in Microbiology* 18: 341-347

Wilkinson, M., McInerney, J.O., Hirt, R.P., Foster, P.G., Embley, T.M.
(2007) Of clades and clans: terms for phylogenetic relationships in
unrooted trees. *Trends in Ecology and Evolution* 22: 114-115

Schliep, K., Lopez, P., Lapointe F.-J., Bapteste E. (2011) Harvesting
Evolutionary Signals in a Forest of Prokaryotic Gene Trees, *Molecular
Biology and Evolution* 28(4): 1393-1405

## See also

[`parsimony`](https://klausvigo.github.io/phangorn/reference/parsimony.md),
Consistency index
[`CI`](https://klausvigo.github.io/phangorn/reference/CI.md), Retention
index [`RI`](https://klausvigo.github.io/phangorn/reference/CI.md),
[`phyDat`](https://klausvigo.github.io/phangorn/reference/as.phyDat.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

Francois-Joseph Lapointe <francois-joseph.lapointe@umontreal.ca>

## Examples

``` r
set.seed(111)
tree <- rtree(10)
getClans(tree)
#>       t3 t8 t2 t1 t5 t6 t9 t10 t4 t7
#>  [1,]  1  0  0  0  0  0  0   0  0  0
#>  [2,]  0  1  0  0  0  0  0   0  0  0
#>  [3,]  0  0  1  0  0  0  0   0  0  0
#>  [4,]  0  0  0  1  0  0  0   0  0  0
#>  [5,]  0  0  0  0  1  0  0   0  0  0
#>  [6,]  0  0  0  0  0  1  0   0  0  0
#>  [7,]  0  0  0  0  0  0  1   0  0  0
#>  [8,]  0  0  0  0  0  0  0   1  0  0
#>  [9,]  0  0  0  0  0  0  0   0  1  0
#> [10,]  0  0  0  0  0  0  0   0  0  1
#> [11,]  0  1  1  1  1  1  1   1  1  1
#> [12,]  1  0  1  1  1  1  1   1  1  1
#> [13,]  1  1  0  1  1  1  1   1  1  1
#> [14,]  1  1  1  0  1  1  1   1  1  1
#> [15,]  1  1  1  1  0  1  1   1  1  1
#> [16,]  1  1  1  1  1  0  1   1  1  1
#> [17,]  1  1  1  1  1  1  0   1  1  1
#> [18,]  1  1  1  1  1  1  1   0  1  1
#> [19,]  1  1  1  1  1  1  1   1  0  1
#> [20,]  1  1  1  1  1  1  1   1  1  0
#> [21,]  1  1  1  0  0  0  0   0  0  0
#> [22,]  0  1  1  0  0  0  0   0  0  0
#> [23,]  0  0  0  0  1  1  1   1  1  1
#> [24,]  0  0  0  0  0  1  1   1  1  1
#> [25,]  0  0  0  0  0  1  1   1  0  0
#> [26,]  0  0  0  0  0  0  1   1  0  0
#> [27,]  0  0  0  0  0  0  0   0  1  1
#> [28,]  0  0  0  1  1  1  1   1  1  1
#> [29,]  1  0  0  1  1  1  1   1  1  1
#> [30,]  1  1  1  1  0  0  0   0  0  0
#> [31,]  1  1  1  1  1  0  0   0  0  0
#> [32,]  1  1  1  1  1  0  0   0  1  1
#> [33,]  1  1  1  1  1  1  0   0  1  1
#> [34,]  1  1  1  1  1  1  1   1  0  0
getClips(tree, all=TRUE)
#>      t3 t8 t2 t1 t5 t6 t9 t10 t4 t7
#> [1,]  1  0  1  0  0  0  0   0  0  0
#> [2,]  0  0  0  0  0  0  1   1  0  0
#> [3,]  0  0  0  0  0  0  0   0  1  1
getSlices(tree)
#>        t3 t8 t2 t1 t5 t6 t9 t10 t4 t7
#>   [1,]  0  0  1  1  1  1  1   1  1  1
#>   [2,]  0  1  0  1  1  1  1   1  1  1
#>   [3,]  0  1  1  0  1  1  1   1  1  1
#>   [4,]  0  1  1  1  0  1  1   1  1  1
#>   [5,]  0  1  1  1  1  0  1   1  1  1
#>   [6,]  0  1  1  1  1  1  0   1  1  1
#>   [7,]  0  1  1  1  1  1  1   0  1  1
#>   [8,]  0  1  1  1  1  1  1   1  0  1
#>   [9,]  0  1  1  1  1  1  1   1  1  0
#>  [10,]  0  1  1  1  0  0  0   0  0  0
#>  [11,]  0  1  1  1  1  0  0   0  0  0
#>  [12,]  0  1  1  1  1  0  0   0  1  1
#>  [13,]  0  1  1  1  1  1  0   0  1  1
#>  [14,]  0  1  1  1  1  1  1   1  0  0
#>  [15,]  1  0  1  0  1  1  1   1  1  1
#>  [16,]  1  0  1  1  0  1  1   1  1  1
#>  [17,]  1  0  1  1  1  0  1   1  1  1
#>  [18,]  1  0  1  1  1  1  0   1  1  1
#>  [19,]  1  0  1  1  1  1  1   0  1  1
#>  [20,]  1  0  1  1  1  1  1   1  0  1
#>  [21,]  1  0  1  1  1  1  1   1  1  0
#>  [22,]  1  0  1  0  0  0  0   0  0  0
#>  [23,]  1  0  1  1  0  0  0   0  0  0
#>  [24,]  1  0  1  1  1  0  0   0  0  0
#>  [25,]  1  0  1  1  1  0  0   0  1  1
#>  [26,]  1  0  1  1  1  1  0   0  1  1
#>  [27,]  1  0  1  1  1  1  1   1  0  0
#>  [28,]  1  1  0  0  1  1  1   1  1  1
#>  [29,]  1  1  0  1  0  1  1   1  1  1
#>  [30,]  1  1  0  1  1  0  1   1  1  1
#>  [31,]  1  1  0  1  1  1  0   1  1  1
#>  [32,]  1  1  0  1  1  1  1   0  1  1
#>  [33,]  1  1  0  1  1  1  1   1  0  1
#>  [34,]  1  1  0  1  1  1  1   1  1  0
#>  [35,]  1  1  0  0  0  0  0   0  0  0
#>  [36,]  1  1  0  1  0  0  0   0  0  0
#>  [37,]  1  1  0  1  1  0  0   0  0  0
#>  [38,]  1  1  0  1  1  0  0   0  1  1
#>  [39,]  1  1  0  1  1  1  0   0  1  1
#>  [40,]  1  1  0  1  1  1  1   1  0  0
#>  [41,]  1  1  1  0  0  1  1   1  1  1
#>  [42,]  1  1  1  0  1  0  1   1  1  1
#>  [43,]  1  1  1  0  1  1  0   1  1  1
#>  [44,]  1  1  1  0  1  1  1   0  1  1
#>  [45,]  1  1  1  0  1  1  1   1  0  1
#>  [46,]  1  1  1  0  1  1  1   1  1  0
#>  [47,]  1  0  0  0  1  1  1   1  1  1
#>  [48,]  1  1  1  0  1  0  0   0  0  0
#>  [49,]  1  1  1  0  1  0  0   0  1  1
#>  [50,]  1  1  1  0  1  1  0   0  1  1
#>  [51,]  1  1  1  0  1  1  1   1  0  0
#>  [52,]  1  1  1  1  0  0  1   1  1  1
#>  [53,]  1  1  1  1  0  1  0   1  1  1
#>  [54,]  1  1  1  1  0  1  1   0  1  1
#>  [55,]  1  1  1  1  0  1  1   1  0  1
#>  [56,]  1  1  1  1  0  1  1   1  1  0
#>  [57,]  0  0  0  1  0  1  1   1  1  1
#>  [58,]  1  0  0  1  0  1  1   1  1  1
#>  [59,]  1  1  1  1  0  0  0   0  1  1
#>  [60,]  1  1  1  1  0  1  0   0  1  1
#>  [61,]  1  1  1  1  0  1  1   1  0  0
#>  [62,]  1  1  1  1  1  0  0   1  1  1
#>  [63,]  1  1  1  1  1  0  1   0  1  1
#>  [64,]  1  1  1  1  1  0  1   1  0  1
#>  [65,]  1  1  1  1  1  0  1   1  1  0
#>  [66,]  0  0  0  0  1  0  1   1  1  1
#>  [67,]  0  0  0  0  0  0  1   1  1  1
#>  [68,]  0  0  0  1  1  0  1   1  1  1
#>  [69,]  1  0  0  1  1  0  1   1  1  1
#>  [70,]  1  1  1  1  1  0  1   1  0  0
#>  [71,]  1  1  1  1  1  1  0   1  0  1
#>  [72,]  1  1  1  1  1  1  0   1  1  0
#>  [73,]  0  0  0  0  1  1  0   1  1  1
#>  [74,]  0  0  0  0  0  1  0   1  1  1
#>  [75,]  0  0  0  0  0  1  0   1  0  0
#>  [76,]  0  0  0  1  1  1  0   1  1  1
#>  [77,]  1  0  0  1  1  1  0   1  1  1
#>  [78,]  1  1  1  1  1  1  0   1  0  0
#>  [79,]  1  1  1  1  1  1  1   0  0  1
#>  [80,]  1  1  1  1  1  1  1   0  1  0
#>  [81,]  0  0  0  0  1  1  1   0  1  1
#>  [82,]  0  0  0  0  0  1  1   0  1  1
#>  [83,]  0  0  0  0  0  1  1   0  0  0
#>  [84,]  0  0  0  1  1  1  1   0  1  1
#>  [85,]  1  0  0  1  1  1  1   0  1  1
#>  [86,]  1  1  1  1  1  1  1   0  0  0
#>  [87,]  0  0  0  0  1  1  1   1  0  1
#>  [88,]  0  0  0  0  0  1  1   1  0  1
#>  [89,]  0  0  0  1  1  1  1   1  0  1
#>  [90,]  1  0  0  1  1  1  1   1  0  1
#>  [91,]  1  1  1  1  1  0  0   0  0  1
#>  [92,]  1  1  1  1  1  1  0   0  0  1
#>  [93,]  0  0  0  0  1  1  1   1  1  0
#>  [94,]  0  0  0  0  0  1  1   1  1  0
#>  [95,]  0  0  0  1  1  1  1   1  1  0
#>  [96,]  1  0  0  1  1  1  1   1  1  0
#>  [97,]  1  1  1  1  1  0  0   0  1  0
#>  [98,]  1  1  1  1  1  1  0   0  1  0
#>  [99,]  0  0  0  0  1  0  0   0  1  1
#> [100,]  0  0  0  0  1  1  0   0  1  1
#> [101,]  0  0  0  0  1  1  1   1  0  0
#> [102,]  0  0  0  0  0  1  0   0  1  1
#> [103,]  0  0  0  1  1  0  0   0  0  0
#> [104,]  0  0  0  1  1  0  0   0  1  1
#> [105,]  0  0  0  1  1  1  0   0  1  1
#> [106,]  0  0  0  1  1  1  1   1  0  0
#> [107,]  1  0  0  1  0  0  0   0  0  0
#> [108,]  1  0  0  1  1  0  0   0  0  0
#> [109,]  1  0  0  1  1  0  0   0  1  1
#> [110,]  1  0  0  1  1  1  0   0  1  1
#> [111,]  1  0  0  1  1  1  1   1  0  0
#> [112,]  1  1  1  1  1  1  0   0  0  0

set.seed(123)
trees <- rmtree(10, 20)
X <- matrix(sample(c("red", "blue", "violet"), 100, TRUE, c(.5,.4, .1)),
   ncol=5, dimnames=list(paste('t',1:20, sep=""), paste('Var',1:5, sep="_")))
x <- phyDat(X, type = "USER", levels = c("red", "blue"), ambiguity="violet")
plot(trees[[1]], "u", tip.color = X[trees[[1]]$tip,1])  # intruders are blue


(divTab <- getDiversity(trees, x, var.names=colnames(X)))
#>    tree variable    E clan # natives # intruder # unknown   E slice # intruder
#> 1     1    Var_1 0.8568636         9         10         1 1.0000000          7
#> 2     1    Var_2 0.6562658        11          7         2 0.7500000          4
#> 3     1    Var_3 0.7500962         8         11         1 0.8018797          8
#> 4     1    Var_4 0.9474428         7         11         2 1.0000000          9
#> 5     1    Var_5 1.0000000         8         10         2 1.0000000          9
#> 6     2    Var_1 0.8568636         9         10         1 1.0000000          7
#> 7     2    Var_2 0.7964530        11          7         2 0.8277294          5
#> 8     2    Var_3 0.9474428         8         11         1 1.0000000          9
#> 9     2    Var_4 0.8224909         7         11         2 0.9166667          8
#> 10    2    Var_5 0.6387640         8         10         2 0.7420981          6
#> 11    3    Var_1 0.8795880         9         10         1 0.9166667          8
#> 12    3    Var_2 0.8982265        11          7         2 1.0000000          5
#> 13    3    Var_3 0.7372138         8         11         1 0.8982265          7
#> 14    3    Var_4 0.8948855         7         11         2 0.9298967          9
#> 15    3    Var_5 1.0000000         8         10         2 1.0000000          9
#> 16    4    Var_1 0.8795880         9         10         1 0.9166667          8
#> 17    4    Var_2 0.8982265        11          7         2 1.0000000          5
#> 18    4    Var_3 0.7372138         8         11         1 0.8982265          7
#> 19    4    Var_4 0.8948855         7         11         2 0.9298967          9
#> 20    4    Var_5 1.0000000         8         10         2 1.0000000          9
#> 21    5    Var_1 0.7966576         9         10         1 0.8982265          7
#> 22    5    Var_2 0.8982265        11          7         2 1.0000000          5
#> 23    5    Var_3 0.8224909         8         11         1 0.9166667          8
#> 24    5    Var_4 0.5699628         7         11         2 0.6934264          6
#> 25    5    Var_5 0.4729033         8         10         2 0.7500000          4
#> 26    6    Var_1 0.7966576         9         10         1 0.8982265          7
#> 27    6    Var_2 0.8982265        11          7         2 1.0000000          5
#> 28    6    Var_3 0.8224909         8         11         1 0.9166667          8
#> 29    6    Var_4 1.0000000         7         11         2 1.0000000         10
#> 30    6    Var_5 0.9397940         8         10         2 1.0000000          8
#> 31    7    Var_1 0.8795880         9         10         1 0.9166667          8
#> 32    7    Var_2 1.0000000        11          7         2 1.0000000          6
#> 33    7    Var_3 0.9474428         8         11         1 1.0000000          9
#> 34    7    Var_4 0.7500962         7         11         2 0.8018797          8
#> 35    7    Var_5 0.9397940         8         10         2 1.0000000          8
#> 36    8    Var_1 0.8795880         9         10         1 0.9166667          8
#> 37    8    Var_2 0.8982265        11          7         2 1.0000000          5
#> 38    8    Var_3 0.8224909         8         11         1 0.9166667          8
#> 39    8    Var_4 0.6846566         7         11         2 0.7964530          7
#> 40    8    Var_5 0.8568636         8         10         2 1.0000000          7
#> 41    9    Var_1 0.7364516         9         10         1 0.7964530          7
#> 42    9    Var_2 0.7964530        11          7         2 0.8277294          5
#> 43    9    Var_3 0.8948855         8         11         1 0.9298967          9
#> 44    9    Var_4 0.7897710         7         11         2 0.7896901          9
#> 45    9    Var_5 0.7591760         8         10         2 1.0000000          6
#> 46   10    Var_1 0.5903090         9         10         1 0.8277294          5
#> 47   10    Var_2 0.8982265        11          7         2 1.0000000          5
#> 48   10    Var_3 0.8224909         8         11         1 0.9166667          8
#> 49   10    Var_4 0.8750481         7         11         2 1.0000000          8
#> 50   10    Var_5 0.9397940         8         10         2 1.0000000          8
#>    # unknown E melange # intruder # unknown bs 1 bs 2 p-score
#> 1          1 1.0000000          6         1   NA   NA       8
#> 2          2 1.0000000          2         2   NA   NA       4
#> 3          1 1.0000000          5         1   NA   NA       6
#> 4          2 1.0000000          8         2   NA   NA       7
#> 5          2 1.0000000          8         2   NA   NA       7
#> 6          1 1.0000000          6         1   NA   NA       7
#> 7          2 1.0000000          3         1   NA   NA       5
#> 8          1 1.0000000          8         1   NA   NA       7
#> 9          2 1.0000000          6         2   NA   NA       5
#> 10         1 0.7500000          4         1   NA   NA       5
#> 11         1 1.0000000          6         1   NA   NA       7
#> 12         2 1.0000000          4         2   NA   NA       6
#> 13         1 1.0000000          5         1   NA   NA       7
#> 14         2 1.0000000          7         1   NA   NA       5
#> 15         2 1.0000000          8         2   NA   NA       8
#> 16         1 1.0000000          6         1   NA   NA       4
#> 17         2 1.0000000          4         2   NA   NA       5
#> 18         1 1.0000000          5         1   NA   NA       6
#> 19         2 1.0000000          7         2   NA   NA       6
#> 20         2 1.0000000          8         2   NA   NA       8
#> 21         1 1.0000000          5         1   NA   NA       7
#> 22         2 1.0000000          4         2   NA   NA       6
#> 23         1 1.0000000          6         1   NA   NA       7
#> 24         0 1.0000000          3         0   NA   NA       5
#> 25         2 1.0000000          2         2   NA   NA       4
#> 26         1 1.0000000          5         1   NA   NA       6
#> 27         1 1.0000000          4         1   NA   NA       6
#> 28         1 1.0000000          6         1   NA   NA       6
#> 29         2 1.0000000          9         2   NA   NA       7
#> 30         2 1.0000000          7         2   NA   NA       7
#> 31         1 1.0000000          6         0   NA   NA       6
#> 32         2 1.0000000          5         2   NA   NA       6
#> 33         1 1.0000000          8         1   NA   NA       7
#> 34         1 1.0000000          5         1   NA   NA       5
#> 35         2 1.0000000          7         2   NA   NA       7
#> 36         1 1.0000000          6         0   NA   NA       8
#> 37         2 1.0000000          4         2   NA   NA       6
#> 38         1 1.0000000          6         1   NA   NA       7
#> 39         2 0.8277294          5         2   NA   NA       6
#> 40         2 1.0000000          6         2   NA   NA       4
#> 41         1 0.8277294          5         1   NA   NA       6
#> 42         2 1.0000000          3         1   NA   NA       4
#> 43         1 1.0000000          7         1   NA   NA       7
#> 44         2 0.7964530          7         2   NA   NA       5
#> 45         2 1.0000000          5         2   NA   NA       6
#> 46         1 1.0000000          3         1   NA   NA       5
#> 47         2 1.0000000          4         2   NA   NA       6
#> 48         1 1.0000000          6         1   NA   NA       5
#> 49         2 1.0000000          7         2   NA   NA       6
#> 50         1 1.0000000          7         1   NA   NA       7
summary(divTab)
#>   Variable Natives_only Intruder_only Clan Slice Melange
#> 1    Var_1            0             0    0     0      10
#> 2    Var_2            0             0    0     0      10
#> 3    Var_3            0             0    0     0      10
#> 4    Var_4            0             0    0     0      10
#> 5    Var_5            0             0    0     0      10
```
