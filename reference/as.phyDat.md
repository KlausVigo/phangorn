# Conversion among Sequence Formats

These functions transform several DNA formats into the `phyDat` format.
`allSitePattern` generates an alignment of all possible site patterns.

## Usage

``` r
phyDat(data, type = "DNA", levels = NULL, return.index = TRUE, ...)

as.phyDat(x, ...)

# S3 method for class 'factor'
as.phyDat(x, ...)

# S3 method for class 'DNAbin'
as.phyDat(x, ...)

# S3 method for class 'AAbin'
as.phyDat(x, ...)

# S3 method for class 'alignment'
as.phyDat(x, type = "DNA", ...)

phyDat2alignment(x)

# S3 method for class 'MultipleAlignment'
as.phyDat(x, ...)

# S3 method for class 'AAStringSet'
as.phyDat(x, ...)

# S3 method for class 'DNAStringSet'
as.phyDat(x, ...)

as.StringSet(x, ...)

# S3 method for class 'phyDat'
as.StringSet(x, ...)

# S3 method for class 'phyDat'
as.MultipleAlignment(x, ...)

# S3 method for class 'phyDat'
as.character(x, allLevels = TRUE, ...)

# S3 method for class 'phyDat'
as.data.frame(x, ...)

# S3 method for class 'phyDat'
as.DNAbin(x, ...)

# S3 method for class 'phyDat'
as.AAbin(x, ...)

genlight2phyDat(x, ambiguity = NA)

acgt2ry(obj)

unalign(x)
```

## Arguments

- data:

  An object containing sequences.

- type:

  Type of sequences ("DNA", "AA", "CODON" or "USER").

- levels:

  Level attributes.

- return.index:

  If TRUE returns a index of the site patterns.

- ...:

  further arguments passed to or from other methods.

- x:

  An object containing sequences.

- allLevels:

  return original data.

- ambiguity:

  character for ambiguous character and no contrast is provided.

- obj:

  as object of class phyDat

## Value

The functions return an object of class `phyDat`.

## Details

If `type` "USER" a vector has to be give to `levels`. For example c("a",
"c", "g", "t", "-") would create a data object that can be used in
phylogenetic analysis with gaps as fifth state. There is a more detailed
example for specifying "USER" defined data formats in the vignette
"phangorn-specials".

`acgt2ry` converts a `phyDat` object of nucleotides into an binary
ry-coded dataset.

`unalign` converts a `phyDat` object of nucleotides or amino acids into
a `DNAbin` or `AAbin` object in list form removing all gaps. These
objects can be exported using
[`write.FASTA`](https://rdrr.io/pkg/ape/man/write.dna.html).

## See also

[`DNAbin`](https://rdrr.io/pkg/ape/man/DNAbin.html),
[`as.DNAbin`](https://rdrr.io/pkg/ape/man/as.alignment.html),
[`baseFreq`](https://klausvigo.github.io/phangorn/reference/baseFreq.md),
[`glance.phyDat`](https://klausvigo.github.io/phangorn/reference/baseFreq.md),
[`read.dna`](https://rdrr.io/pkg/ape/man/read.dna.html),
[`read.nexus.data`](https://rdrr.io/pkg/ape/man/read.nexus.data.html)
and the chapter 1 in the
`vignette("phangorn-specials", package="phangorn")` and the example of
[`pmlMix`](https://klausvigo.github.io/phangorn/reference/pmlMix.md) for
the use of `allSitePattern`

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
# transform as characters
LauraChar <- as.character(Laurasiatherian)
# and back
Laura <- phyDat(LauraChar)
all.equal(Laurasiatherian, Laura)
#> [1] TRUE
LauraDNAbin <- as.DNAbin(Laurasiatherian)
all.equal(Laurasiatherian, as.phyDat(LauraDNAbin))
#> [1] TRUE
```
