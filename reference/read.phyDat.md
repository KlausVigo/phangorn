# Import and export sequence alignments

These functions read and write sequence alignments.

## Usage

``` r
read.phyDat(file, format = "phylip", type = "DNA", ...)

write.phyDat(x, file, format = "phylip", colsep = "", nbcol = -1, ...)
```

## Arguments

- file:

  a file name specified by either a variable of mode character, or a
  double-quoted string.

- format:

  File format of the sequence alignment (see details). Several popular
  formats are supported: "phylip", "interleaved", "sequential",
  "clustal", "fasta" or "nexus", or any unambiguous abbreviation of
  these.

- type:

  Type of sequences ("DNA", "AA", "CODON" or "USER").

- ...:

  further arguments passed to or from other methods.

- x:

  An object of class `phyDat`.

- colsep:

  a character used to separate the columns (a single space by default).

- nbcol:

  a numeric specifying the number of columns per row (-1 by default);
  may be negative implying that the nucleotides are printed on a single
  line.

## Value

`read.phyDat` returns an object of class phyDat, `write.phyDat` write an
alignment to a file.

## Details

`write.phyDat` calls the function
[`write.dna`](https://rdrr.io/pkg/ape/man/write.dna.html) or
[`write.nexus.data`](https://rdrr.io/pkg/ape/man/write.nexus.data.html)
and `read.phyDat` calls the function
[`read.dna`](https://rdrr.io/pkg/ape/man/read.dna.html) or
`read.nexus.data`, so see for more details over there.

You may import data directly with
[`read.dna`](https://rdrr.io/pkg/ape/man/read.dna.html) or
[`read.nexus.data`](https://rdrr.io/pkg/ape/man/read.nexus.data.html)
and convert the data to class phyDat.

## References

Anonymous. FASTA format description.
<https://www.ncbi.nlm.nih.gov/blast/fasta.shtml> Felsenstein, J. (1993)
Phylip (Phylogeny Inference Package) version 3.5c. Department of
Genetics, University of Washington.
<https://phylipweb.github.io/phylip/>

## See also

[`read.dna`](https://rdrr.io/pkg/ape/man/read.dna.html),
[`read.GenBank`](https://rdrr.io/pkg/ape/man/read.GenBank.html),
[`phyDat`](https://klausvigo.github.io/phangorn/reference/as.phyDat.md),
[`read.alignment`](https://rdrr.io/pkg/seqinr/man/read.alignment.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")
```
