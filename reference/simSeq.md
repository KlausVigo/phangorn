# Simulate sequences.

Simulate sequences from a given evolutionary tree or an `pml` object.

## Usage

``` r
simSeq(x, ...)

# S3 method for class 'phylo'
simSeq(x, l = 1000, Q = NULL, bf = NULL,
  rootseq = NULL, type = "DNA", model = NULL, levels = NULL,
  rate = 1, ancestral = FALSE, code = 1, ...)

# S3 method for class 'pml'
simSeq(x, ancestral = FALSE, ...)
```

## Arguments

- x:

  a phylogenetic tree `tree` with edge lengths, i.e. an object of class
  `phylo` or and object of class `pml`.

- ...:

  Further arguments passed to or from other methods.

- l:

  The length of the sequence to simulate.

- Q:

  Either a numeric matrix of size Nstates × Nstates, giving the
  transition rates between states or a vector representing the lower
  triangular of these matrix (see details).

- bf:

  Base frequencies.

- rootseq:

  A vector of length `l` containing the root sequence. If not provided,
  the root sequence is randomly generated.

- type:

  Type of sequences ("DNA", "AA", "CODON" or "USER").

- model:

  Amino acid model of evolution to employ, for example "WAG", "JTT",
  "Dayhoff" or "LG" or structural 3Di alphabet (e.g. "Q_3Di"). For a
  full list of supported models, type `phangorn:::.aa_3Di_models`.
  Ignored if type is not equal to "AA".

- levels:

  A character vector of the different character tokens. Ignored unless
  type = "USER".

- rate:

  A numerical value greater than zero giving the mutation rate or scaler
  for edge lengths.

- ancestral:

  Logical specifying whether to return ancestral sequences.

- code:

  The ncbi genetic code number for translation (see details). By default
  the standard genetic code is used.

## Value

`simSeq` returns an object of class phyDat.

## Details

`simSeq` is a generic function to simulate sequence alignments along a
phylogeny. It is quite flexible and can generate DNA, RNA, amino acids,
codon, morphological or binary sequences. simSeq can take as input a
phylogenetic tree of class `phylo`, or a `pml` object; it will return an
object of class `phyDat`. There is also a more low level version, which
lacks rate variation, but one can combine different alignments with
their own rates (see example). The rate parameter acts like a scaler for
the edge lengths.

For codon models `type="CODON"`, two additional arguments `dnds` for the
dN/dS ratio and `tstv` for the transition transversion ratio can be
supplied.

So far `simSeq` is limited to time reversible models. `simSeq` will
normalize the rate matrix A composed from Q and bf so that every row of
A sums to zero and expected rate is one, see formulas 13.14 and 13.15 on
page 205 in Felsenstein (2004). The edge lengths are should represent
the expected number of mutations per site.

**Defaults:**

If `x` is a tree of class `phylo`, then sequences will be generated with
the default Jukes-Cantor DNA model (`"JC"`).

If `bf` is not specified, then all states will be treated as equally
probable.

If `Q` is not specified, then a uniform rate matrix will be employed.

## References

Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates,
Sunderland.

## See also

[`phyDat`](https://klausvigo.github.io/phangorn/reference/as.phyDat.md),
[`pml`](https://klausvigo.github.io/phangorn/reference/pml.md),
[`SOWH.test`](https://klausvigo.github.io/phangorn/reference/SOWH.test.md)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data(Laurasiatherian)
tree <- nj(dist.ml(Laurasiatherian))
fit <- pml(tree, Laurasiatherian, k=4)
fit <- optim.pml(fit, optNni=TRUE, model="GTR", optGamma=TRUE)
data <- simSeq(fit)
} # }


tree <- rtree(5)
plot(tree)
nodelabels()


# Example for simple DNA alignment
data <- simSeq(tree, l = 10, type="DNA", bf=c(.1,.2,.3,.4), Q=1:6,
               ancestral=TRUE)
as.character(data)
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> t1 "t"  "c"  "a"  "g"  "t"  "c"  "c"  "g"  "t"  "a"  
#> t4 "t"  "g"  "t"  "t"  "t"  "g"  "g"  "c"  "t"  "c"  
#> t5 "t"  "g"  "g"  "t"  "t"  "t"  "t"  "g"  "c"  "g"  
#> t3 "t"  "g"  "t"  "g"  "c"  "g"  "t"  "c"  "t"  "c"  
#> t2 "a"  "a"  "t"  "g"  "t"  "g"  "t"  "a"  "t"  "g"  
#> 6  "g"  "t"  "t"  "t"  "t"  "t"  "t"  "t"  "g"  "t"  
#> 7  "g"  "g"  "g"  "t"  "t"  "t"  "t"  "t"  "g"  "t"  
#> 8  "g"  "g"  "t"  "g"  "t"  "t"  "t"  "t"  "t"  "g"  
#> 9  "c"  "t"  "t"  "g"  "c"  "g"  "t"  "c"  "t"  "t"  


# Example to simulate discrete Gamma rate variation
rates <- discrete.gamma(1,4)
data1 <- simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[1])
data2 <- simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[2])
data3 <- simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[3])
data4 <- simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[4])
data <- c(data1,data2, data3, data4)

write.phyDat(data, file="temp.dat", format="sequential", nbcol = -1,
  colsep = "")
unlink("temp.dat")
```
