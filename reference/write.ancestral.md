# Export and convenience functions for ancestral reconstructions

`write.ancestral` allows to export ancestral reconstructions. It writes
out the tree, a tab delimited text file with the probabilities and the
alignment. `ancestral` generates an object of class ancestral.

## Usage

``` r
write.ancestral(x, file = "ancestral")

as.ancestral(tree, prob, align = NULL)

# S3 method for class 'ancestral'
print(x, ...)

# S3 method for class 'ancestral'
as.phyDat(x, ...)

# S3 method for class 'ancestral'
as.data.frame(x, ...)
```

## Arguments

- x:

  an object of class ancestral

- file:

  a file name. File endings are added.

- tree:

  an object of class phylo.

- prob:

  an data.frame containing a matrix of posterior probabilities for each
  state and site.

- align:

  an object of class phyDat.

- ...:

  Further arguments passed to or from other methods.

## Value

`write.ancestral` returns the input x invisibly.

## Details

This allows also to read in reconstruction made by iqtree to use the
plotting capabilities of R.

## See also

[`anc_pml`](https://klausvigo.github.io/phangorn/reference/ancestral.pml.md),
[`plotAnc`](https://klausvigo.github.io/phangorn/reference/plot.ancestral.md)

## Examples

``` r
data(Laurasiatherian)
fit <- pml_bb(Laurasiatherian[,1:100], "JC", rearrangement = "none")
#> optimize edge weights:  -2142.479 --> -2136.397 
#> optimize edge weights:  -2136.397 --> -2136.397 
#> optimize edge weights:  -2136.397 --> -2136.397 
anc_ml <- anc_pml(fit)
write.ancestral(anc_ml)
# Can be also results from iqtree
align <- read.phyDat("ancestral_align.fasta", format="fasta")
tree <- read.tree("ancestral_tree.nwk")
df <- read.table("ancestral_state.tsv", header=TRUE)
anc_ml_disc <- as.ancestral(tree, df, align)
plotAnc(anc_ml_disc, 20)

unlink(c("ancestral_align.fasta", "ancestral_tree.nwk",
         "ancestral_state.tsv"))
```
