# Plot of a Sequence Alignment

This function plots an image of an alignment of sequences.

## Usage

``` r
# S3 method for class 'phyDat'
image(x, ...)

# S3 method for class 'ancestral'
image(x, ...)
```

## Arguments

- x:

  an object containing sequences, an object of class `phyDat`.

- ...:

  further arguments passed to or from other methods.

## Value

`image.phyDat` returns silently x.

## Details

A wrapper for using
[`image.DNAbin`](https://rdrr.io/pkg/ape/man/image.DNAbin.html) and
[`image.AAbin`](https://rdrr.io/pkg/ape/man/AAbin.html). Codons triplets
are transformed and handled as nucleotide sequences. So far it is not
yet possible to plot data with ` type="USER"`.

## See also

[`image.DNAbin`](https://rdrr.io/pkg/ape/man/image.DNAbin.html),
[`image.AAbin`](https://rdrr.io/pkg/ape/man/AAbin.html)

## Examples

``` r
data("chloroplast")
image(chloroplast[, 1:50], scheme="Clustal", show.aa = TRUE)
```
