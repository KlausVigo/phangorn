# Plot phylogeny of a pml object

`plot.pml` is a wrapper around `plot.phylo` with different default
values for unrooted, ultrametric and tip dated phylogenies.

## Usage

``` r
# S3 method for class 'pml'
plot(x, type = "phylogram", direction = "rightwards", ...,
  adj = NULL, digits = 2, method = "FBP")
```

## Arguments

- x:

  an object of class `pml`.

- type:

  a character string specifying the type of phylogeny to be drawn; it
  must be one of "phylogram" (the default), "cladogram", "fan",
  "unrooted", "radial", "tidy", or any unambiguous abbreviation of
  these.

- direction:

  a character string specifying the direction of the tree. Four values
  are possible: "rightwards" (the default), "leftwards", "upwards", and
  "downwards".

- ...:

  further parameters to be passed to `plot.phylo`.

- adj:

  one or two numeric values specifying the horizontal and vertical
  justification of the text or symbols of the support values.

- digits:

  integer indicating the number of decimal places.

- method:

  either "FBP" the classical bootstrap (default), "TBE" (transfer
  bootstrap) or "MCC" for assigning clade credibilities.

## Value

`plot.pml` returns the `pml` object x.

## See also

[`plot.phylo`](https://rdrr.io/pkg/ape/man/plot.phylo.html),
[`axisPhylo`](https://rdrr.io/pkg/ape/man/axisPhylo.html),
[`add.scale.bar`](https://rdrr.io/pkg/ape/man/add.scale.bar.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
fdir <- system.file("extdata/trees", package = "phangorn")
tmp <- read.csv(file.path(fdir,"H3N2_NA_20.csv"))
H3N2 <- read.phyDat(file.path(fdir,"H3N2_NA_20.fasta"), format="fasta")
dates <- setNames(tmp$numdate_given, tmp$name)

fit_td <- pml_bb(H3N2, model="JC", method="tipdated", tip.dates=dates,
                 rearrangement="none", control = pml.control(trace = 0))
plot(fit_td, show.tip.label = FALSE)

# Same as:
# root_time <- max(dates) - max(node.depth.edgelength(fit_td$tree))
# plot(fit_td$tree, show.tip.label = FALSE)
# axisPhylo(root.time = root_time, backward = FALSE)
plot(fit_td, show.tip.label = FALSE, direction="up")


fit_unrooted <- pml_bb(H3N2, model="JC", rearrangement="none",
                       control = pml.control(trace = 0))
plot(fit_unrooted, cex=.5)

```
