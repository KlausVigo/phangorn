# Function to import partitioned data from nexus files

`read.nexus.partitions` reads in sequences in NEXUS format and splits
the data according to the charsets given in the SETS block.

## Usage

``` r
read.nexus.partitions(file, return = "list", ...)
```

## Arguments

- file:

  a file name.

- return:

  either returns a list where each element is a 'phyDat' object or an
  object of class 'multiphyDat'

- ...:

  Further arguments passed to or from other methods.

## Value

a list where each element is a 'phyDat' object or an object of class
'multiphyDat'.

## See also

[`read.nexus.data`](https://rdrr.io/pkg/ape/man/read.nexus.data.html),
[`read.phyDat`](https://klausvigo.github.io/phangorn/reference/read.phyDat.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
tree <- rtree(10)
dat <- simSeq(tree, l=24)
fcat <- function(..., file = zz) cat(..., file=file, sep="", append=TRUE)
zz <- tempfile(pattern="file", tmpdir=tempdir(), fileext=".nex")
write.phyDat(dat, file=zz, format="nexus")
fcat("BEGIN SETS;\n")
fcat("  Charset codon1 = 1-12/3;\n")
fcat("  Charset codon2 = 2-12/3;\n")
fcat("  Charset codon3 = 3-12/3;\n")
fcat("  Charset range = 16-18;\n")
fcat("  Charset range2 = 13-15 19-21;\n")
fcat("  Charset singles = 22 23 24;\n")
fcat("END;\n")

tmp <- read.nexus.partitions(zz)
tmp
#> $codon1
#> 10 sequences with 4 character and 4 different site patterns.
#> The states are a c g t 
#> 
#> $codon2
#> 10 sequences with 4 character and 4 different site patterns.
#> The states are a c g t 
#> 
#> $codon3
#> 10 sequences with 4 character and 4 different site patterns.
#> The states are a c g t 
#> 
#> $range
#> 10 sequences with 3 character and 3 different site patterns.
#> The states are a c g t 
#> 
#> $range2
#> 10 sequences with 6 character and 6 different site patterns.
#> The states are a c g t 
#> 
#> $singles
#> 10 sequences with 3 character and 3 different site patterns.
#> The states are a c g t 
#> 
unlink(zz)
```
