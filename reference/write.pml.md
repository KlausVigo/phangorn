# Export pml objects

`write.pml` writes out the ML tree and the model parameters.

## Usage

``` r
write.pml(x, file = "pml", save_rds = TRUE, digits = 10, ...)
```

## Arguments

- x:

  an object of class pml.

- file:

  a file name. File endings are added.

- save_rds:

  logical, if TRUE saves the pml object as a rds file, otherwise the
  alignment is saved as a fasta file.

- digits:

  default is 10, i.e. edge length for the bootstrap trees are exported.
  For digits larger smaller than zero no edge length are exported.

- ...:

  Further arguments passed to or from other methods.

## Value

`write.pml` returns the input x invisibly.

## Details

`write.pml` creates several files. It exports the alignment as fasta
file. It writes out the ML tree in a newick file and the estimates
parameters in a txt file. It should be possible to (re-)create the pml
object up to numerical inaccuracies and this is possible with the \*.rds
file. If bootstrap trees exist these are additionally exported in a
compressed nexus file. Additionally several plots are returned. The
maximum likelihood tree, with support values, if these are available. If
an bootstrapped trees exist, a consensus tree, a consensus network (\<
200 tips) and terrace plot. And last but not least the distribution of
the rates. It might be better to adopt these on the dataset.

## See also

[`ancestral.pml`](https://klausvigo.github.io/phangorn/reference/ancestral.pml.md),
[`plotAnc`](https://klausvigo.github.io/phangorn/reference/plot.ancestral.md)

## Examples

``` r
data(woodmouse)
fit <- pml_bb(woodmouse, "JC", rearrangement = "none")
#> optimize edge weights:  -1864.042 --> -1857.165 
#> optimize edge weights:  -1857.165 --> -1857.165 
#> optimize edge weights:  -1857.165 --> -1857.165 
write.pml(fit, "woodmouse")
unlink(c("woodmouse.txt", "woodmouse_tree.nwk", "woodmouse_align.fasta",
       "woodmouse_tree.pdf", "woodmouse.rds", "woodmouse_rates.pdf"))
```
