# Pairwise Distances from a Phylogenetic Network

`cophenetic.networx` computes the pairwise distances between the pairs
of tips from a phylogenetic network using its branch lengths.

## Usage

``` r
# S3 method for class 'networx'
cophenetic(x)
```

## Arguments

- x:

  an object of class `networx`.

## Value

an object of class `dist`, names are set according to the tip labels (as
given by the element `tip.label` of the argument `x`).

## See also

[`cophenetic`](https://rdrr.io/r/stats/cophenetic.html) for the generic
function, `neighborNet` to construct a network from a distance matrix

## Author

Klaus Schliep

## Examples

``` r
example(neighborNet)
#> 
#> nghbrN> data(yeast)
#> 
#> nghbrN> dm <- dist.ml(yeast)
#> 
#> nghbrN> nnet <- neighborNet(dm)
#> 
#> nghbrN> plot(nnet)

#> 
#> nghbrN> plot(nnet, type="outline")

cophenetic(nnet)
#>            Sklu       Calb       Scas       Scer       Spar       Smik
#> Calb 0.54213286                                                       
#> Scas 0.39014281 0.52768726                                            
#> Scer 0.38167995 0.53910531 0.34350480                                 
#> Spar 0.37937970 0.53739537 0.34179485 0.08678458                      
#> Smik 0.38238343 0.54039910 0.34479859 0.13763679 0.12360062           
#> Skud 0.38397653 0.54770657 0.35210606 0.16108716 0.15082553 0.15382927
#> Sbay 0.38096004 0.54517923 0.34957872 0.17717456 0.16838347 0.17138721
#>            Skud
#> Calb           
#> Scas           
#> Scer           
#> Spar           
#> Smik           
#> Skud           
#> Sbay 0.15679569
```
