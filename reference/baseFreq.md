# Summaries of alignments

`baseFreq` computes the frequencies (absolute or relative) of the states
from a sample of sequences. `glance` computes some useful information
about the alignment. `composition\_test` computes a \\\chi^2\\-test
testing if the state composition for a species differs.

## Usage

``` r
baseFreq(obj, freq = FALSE, all = FALSE, drop.unused.levels = FALSE)

# S3 method for class 'phyDat'
glance(x, ...)

composition_test(obj)

# S3 method for class 'phyDat'
summary(object, ...)

# S3 method for class 'summary.phyDat'
print(x, ..., digits = max(3L, getOption("digits") -
  3L))
```

## Arguments

- freq:

  logical, if 'TRUE', frequencies or counts are returned otherwise
  proportions

- all:

  all a logical; if all = TRUE, all counts of bases, ambiguous codes,
  missing data, and alignment gaps are returned as defined in the
  contrast.

- drop.unused.levels:

  logical, drop unused levels

- ...:

  further arguments passed to or from other methods.

- object, obj, x:

  as object of class phyDat

- digits:

  minimal number of significant digits.

## Value

`baseFreq` returns a named vector and `glance` a one row `data.frame`.

## See also

[`phyDat`](https://klausvigo.github.io/phangorn/reference/as.phyDat.md)`, `[`base.freq`](https://rdrr.io/pkg/ape/man/base.freq.html)`, `[`glance`](https://generics.r-lib.org/reference/glance.html)

## Author

Klaus Schliep

## Examples

``` r
data(Laurasiatherian)
data(chloroplast)
# base frequencies
baseFreq(Laurasiatherian)
#>         a         c         g         t 
#> 0.3321866 0.1990791 0.2040652 0.2646691 
baseFreq(Laurasiatherian, all=TRUE)
#>         a         c         g         t         u         m         r         w 
#> 0.3321866 0.1990791 0.2040652 0.2646691 0.0000000 0.0000000 0.0000000 0.0000000 
#>         s         y         k         v         h         d         b         n 
#> 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
#>         ?         - 
#> 0.0000000 0.0000000 
baseFreq(Laurasiatherian, freq=TRUE)
#>     a     c     g     t 
#> 49633 29745 30490 39545 
baseFreq(chloroplast)
#>           A           R           N           D           C           Q 
#> 0.086702955 0.051669804 0.036230253 0.040435459 0.006650569 0.039678317 
#>           E           G           H           I           L           K 
#> 0.046635835 0.090775149 0.029395514 0.073585987 0.101712777 0.040148973 
#>           M           F           P           S           T           W 
#> 0.024852664 0.051260539 0.048109192 0.054176557 0.051618646 0.021905951 
#>           Y           V 
#> 0.031636245 0.072818613 

# some statistics
glance(Laurasiatherian)
#>   nseq nchar unique_site_pattern parsimony_informative_sites const_sites
#> 1   47  3179                1605                        1400        1354
#>   duplicated_seq gaps ambiguous type
#> 1              0    0         0  DNA
glance(chloroplast)
#>   nseq nchar unique_site_pattern parsimony_informative_sites const_sites
#> 1   19  5144                2775                        2032        2190
#>   duplicated_seq gaps ambiguous type
#> 1              0    0         0   AA
composition_test(Laurasiatherian)[1:10,]
#>            statistic df     p-value
#> Platypus   2.1137955  3 0.549126893
#> Wallaroo   2.4109765  3 0.491594621
#> Possum     3.5419454  3 0.315362557
#> Bandicoot  7.9087006  3 0.047936758
#> Opposum   12.5176626  3 0.005804765
#> Armadillo  1.5323619  3 0.674821614
#> Elephant   0.2412612  3 0.970668519
#> Aardvark   1.7363022  3 0.628893376
#> Tenrec     0.9995630  3 0.801357693
#> Hedghog    5.2957058  3 0.151381280
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
```
