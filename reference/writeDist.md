# Writing and reading distances in phylip and nexus format

`readDist`, `writeDist` and `write.nexus.dist` are useful to exchange
distance matrices with other phylogenetic programs.

## Usage

``` r
writeDist(x, file = "", format = "phylip", ...)

write.nexus.dist(x, file = "", append = FALSE, upper = FALSE,
  diag = TRUE, digits = getOption("digits"), taxa = !append)

readDist(file, format = "phylip")

read.nexus.dist(file)

# S3 method for class 'dist'
unique(x, incomparables, ...)
```

## Arguments

- x:

  A `dist` object.

- file:

  A file name.

- format:

  file format, default is "phylip", only other option so far is "nexus".

- ...:

  Further arguments passed to or from other methods.

- append:

  logical. If TRUE the nexus blocks will be added to a file.

- upper:

  logical value indicating whether the upper triangle of the distance
  matrix should be printed.

- diag:

  logical value indicating whether the diagonal of the distance matrix
  should be printed.

- digits:

  passed to format inside of `write.nexus.dist`.

- taxa:

  logical. If TRUE a taxa block is added.

- incomparables:

  Not used so far.

## Value

an object of class `dist`

## References

Maddison, D. R., Swofford, D. L. and Maddison, W. P. (1997) NEXUS: an
extensible file format for systematic information. *Systematic Biology*,
**46**, 590–621.

## See also

To compute distance matrices see
[`dist.ml`](https://klausvigo.github.io/phangorn/reference/dist.hamming.md)
[`dist.dna`](https://rdrr.io/pkg/ape/man/dist.dna.html) and
[`dist.p`](https://klausvigo.github.io/phangorn/reference/dist.p.md) for
pairwise polymorphism p-distances

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
data(yeast)
dm <- dist.ml(yeast)
writeDist(dm)
#> 8 
#> Scer 0 0.0867845773228219 0.137636794490986 0.161087155171287 0.177174556271658 0.345373797843901 0.381679949757669 0.537234265198945
#> Spar 0.0867845773228219 0 0.123600617669834 0.149284159950583 0.166375240744154 0.343965239042025 0.380724585630441 0.537427661847541
#> Smik 0.137636794490986 0.123600617669834 0 0.155370638334037 0.173395434179653 0.345111830071115 0.381038543686742 0.5378791145715
#> Skud 0.161087155171287 0.149284159950583 0.155370638334037 0 0.156795687803861 0.351411419602885 0.383976534826637 0.548399135029346
#> Sbay 0.177174556271658 0.166375240744154 0.173395434179653 0.156795687803861 0 0.34592156486885 0.3809600418506 0.548845409688793
#> Scas 0.345373797843901 0.343965239042025 0.345111830071115 0.351411419602885 0.34592156486885 0 0.390140463219943 0.52768391927111
#> Sklu 0.381679949757669 0.380724585630441 0.381038543686742 0.383976534826637 0.3809600418506 0.390140463219943 0 0.542132862145574
#> Calb 0.537234265198945 0.537427661847541 0.5378791145715 0.548399135029346 0.548845409688793 0.52768391927111 0.542132862145574 0
write.nexus.dist(dm)
#> #NEXUS
#> 
#> BEGIN TAXA;
#>  DIMENSIONS ntax=8;
#>  TAXLABELS Scer Spar Smik Skud Sbay Scas Sklu Calb ;
#> END;
#> 
#> BEGIN DISTANCES; 
#>  FORMAT TRIANGLE = LOWER;
#>  Matrix 
#>  Scer 0.00000000
#>  Spar 0.08678458 0.00000000
#>  Smik 0.13763679 0.12360062 0.00000000
#>  Skud 0.16108716 0.14928416 0.15537064 0.00000000
#>  Sbay 0.17717456 0.16637524 0.17339543 0.15679569 0.00000000
#>  Scas 0.34537380 0.34396524 0.34511183 0.35141142 0.34592156 0.00000000
#>  Sklu 0.38167995 0.38072459 0.38103854 0.38397653 0.38096004 0.39014046 0.00000000
#>  Calb 0.53723427 0.53742766 0.53787911 0.54839914 0.54884541 0.52768392 0.54213286 0.00000000
#>  ;
#> END; 
```
