# Maximum agreement subtree

`mast` computes the maximum agreement subtree (MAST).

## Usage

``` r
mast(x, y, tree = TRUE, rooted = TRUE)
```

## Arguments

- x:

  a tree, i.e. an object of class `phylo`.

- y:

  a tree, i.e. an object of class `phylo`.

- tree:

  a logical, if TRUE returns a tree other wise the tip labels of the the
  maximum agreement subtree.

- rooted:

  logical if TRUE treats trees as rooted otherwise unrooted.

## Value

`mast` returns a vector of the tip labels in the MAST.

## Details

The code is derived from the code example in Valiente (2009). The
version for the unrooted trees is much slower.

## References

G. Valiente (2009). *Combinatorial Pattern Matching Algorithms in
Computational Biology using Perl and R*. Taylor & Francis/CRC Press

## See also

[`SPR.dist`](https://klausvigo.github.io/phangorn/reference/treedist.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com> based on code of Gabriel
Valiente

## Examples

``` r
tree1 <- rtree(100)
tree2 <- rNNI(tree1, 20)
tips <- mast(tree1, tree2)
```
