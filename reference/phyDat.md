# Generic functions for class phyDat

These functions help to manipulate alignments of class phyDat.

## Usage

``` r
# S3 method for class 'phyDat'
print(x, ...)

# S3 method for class 'phyDat'
cbind(..., gaps = "-", compress = TRUE)

# S3 method for class 'phyDat'
rbind(...)

# S3 method for class 'phyDat'
c(..., gaps = "-", compress = TRUE)

# S3 method for class 'phyDat'
subset(x, subset, select, site.pattern = TRUE, ...)

# S3 method for class 'phyDat'
x[i, j, ..., drop = FALSE]

# S3 method for class 'phyDat'
unique(x, incomparables = FALSE, identical = TRUE, ...)

removeUndeterminedSites(x, ...)

removeAmbiguousSites(x)

allSitePattern(n, levels = NULL, names = NULL, type = "DNA", code = 1)
```

## Arguments

- x:

  An object containing sequences.

- ...:

  further arguments passed to or from other methods.

- gaps:

  character for gaps (default is '-').

- compress:

  logical, compress data to store only site patterns.

- subset:

  a subset of taxa.

- select:

  a subset of characters.

- site.pattern:

  select site pattern or sites (see details).

- i, j:

  indices of the rows and/or columns to select or to drop. They may be
  numeric, logical, or character (in the same way than for standard R
  objects).

- drop:

  for compatibility with the generic (unused).

- incomparables:

  for compatibility with unique.

- identical:

  if TRUE (default) sequences have to be identical, if FALSE sequences
  are considered duplicates if distance between sequences is zero
  (happens frequently with ambiguous sites).

- n:

  Number of sequences.

- levels:

  Level attributes.

- names:

  Names of sequences.

- type:

  Type of sequences ("DNA", "AA" or "USER").

- code:

  The ncbi genetic code number for translation. By default the standard
  genetic code is used.

## Value

The functions return an object of class `phyDat`.

## Details

`allSitePattern` generates all possible site patterns and can be useful
in simulation studies. For further details see the vignette
AdvancedFeatures.

The generic function `c` can be used to to combine sequences and
`unique` to get all unique sequences or unique haplotypes.

`phyDat` stores identical columns of an alignment only once and keeps an
index of the original positions. This saves memory and especially
computations as these are usually need to be done only once for each
site pattern. In the example below the matrix x in the example has 8
columns, but column 1 and 2 and also 3 and 5 are identical. The `phyDat`
object y has only 6 site pattern. If argument `site.pattern=FALSE` the
indexing behaves like on the original matrix x. `site.pattern=TRUE` can
be useful inside functions.

## See also

[`DNAbin`](https://rdrr.io/pkg/ape/man/DNAbin.html),
[`as.DNAbin`](https://rdrr.io/pkg/ape/man/as.alignment.html),
[`baseFreq`](https://klausvigo.github.io/phangorn/reference/baseFreq.md),
[`glance.phyDat`](https://klausvigo.github.io/phangorn/reference/baseFreq.md),
[`dna2codon`](https://klausvigo.github.io/phangorn/reference/dna2codon.md),
[`read.dna`](https://rdrr.io/pkg/ape/man/read.dna.html),
[`read.nexus.data`](https://rdrr.io/pkg/ape/man/read.nexus.data.html)
and the chapter 1 in the
[`vignette("AdvancedFeatures", package="phangorn")`](https://klausvigo.github.io/phangorn/articles/AdvancedFeatures.md)
and the example of
[`pmlMix`](https://klausvigo.github.io/phangorn/reference/pmlMix.md) for
the use of `allSitePattern`.

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(Laurasiatherian)
class(Laurasiatherian)
#> [1] "phyDat"
Laurasiatherian
#> 47 sequences with 3179 character and 1605 different site patterns.
#> The states are a c g t 
# base frequencies
baseFreq(Laurasiatherian)
#>         a         c         g         t 
#> 0.3321866 0.1990791 0.2040652 0.2646691 
# summary statistics
summary(Laurasiatherian)
#> Alignment statistics 
#> 
#> Type:  DNA 
#> Sequences:  47 
#> Columns:  3179 
#> Site pattern:  1605 
#> Parsimony informative sites:  1400 
#> Constant sites:  1354 
#> Number of gaps:  0 
#> Number of ambiguous states:  0 
#> Duplicated sequences:  0 
#> 
#> State frequencies (empirical):
#> p(a) = 0.3321866 
#> p(c) = 0.1990791 
#> p(g) = 0.2040652 
#> p(t) = 0.2646691 
#> 
#> Composition Test (Chisq) 
#>           statistic df  p-value
#> Opposum      12.518  3 0.005805
#> Baboon        8.888  3 0.030823
#> Human         8.319  3 0.039866
#> Bandicoot     7.909  3 0.047937

# subsetting phyDat objects
# the first 5 sequences
subset(Laurasiatherian, subset=1:5)
#> 5 sequences with 3179 character and 1605 different site patterns.
#> The states are a c g t 
# the first 5 characters
subset(Laurasiatherian, select=1:5, site.pattern = FALSE)
#> 47 sequences with 5 character and 5 different site patterns.
#> The states are a c g t 
# subsetting with []
Laurasiatherian[1:5, 1:20]
#> 5 sequences with 20 character and 17 different site patterns.
#> The states are a c g t 
# short for
subset(Laurasiatherian, subset=1:5, select=1:20, site.pattern = FALSE)
#> 5 sequences with 20 character and 17 different site patterns.
#> The states are a c g t 
# the first 5 site patterns (often more than 5 characters)
subset(Laurasiatherian, select=1:5, site.pattern = TRUE)
#> 47 sequences with 454 character and 5 different site patterns.
#> The states are a c g t 

x <- matrix(c("a", "a", "c", "g", "c", "t", "a", "g",
              "a", "a", "c", "g", "c", "t", "a", "g",
              "a", "a", "c", "c", "c", "t", "t", "g"), nrow=3, byrow = TRUE,
            dimnames = list(c("t1", "t2", "t3"), 1:8))
(y <- phyDat(x))
#> 3 sequences with 8 character and 6 different site patterns.
#> The states are a c g t 

subset(y, 1:2)
#> 2 sequences with 8 character and 6 different site patterns.
#> The states are a c g t 
subset(y, 1:2, compress=TRUE)
#> 2 sequences with 8 character and 4 different site patterns.
#> The states are a c g t 

subset(y, select=1:3, site.pattern = FALSE) |> as.character()
#>    [,1] [,2] [,3]
#> t1 "a"  "a"  "c" 
#> t2 "a"  "a"  "c" 
#> t3 "a"  "a"  "c" 
subset(y, select=1:3, site.pattern = TRUE) |> as.character()
#>    [,1] [,2] [,3] [,4] [,5]
#> t1 "a"  "a"  "c"  "c"  "g" 
#> t2 "a"  "a"  "c"  "c"  "g" 
#> t3 "a"  "a"  "c"  "c"  "c" 
y[,1:3] # same as subset(y, select=1:3, site.pattern = FALSE)
#> 3 sequences with 3 character and 2 different site patterns.
#> The states are a c g t 

# Compute all possible site patterns
# for nucleotides there $4 ^ (number of tips)$ patterns
allSitePattern(5)
#> 5 sequences with 1024 character and 1024 different site patterns.
#> The states are a c g t 
```
