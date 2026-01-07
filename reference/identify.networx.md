# Identify splits in a network

`identify.networx` reads the position of the graphics pointer when the
mouse button is pressed. It then returns the split belonging to the edge
closest to the pointer. The network must be plotted beforehand.

## Usage

``` r
# S3 method for class 'networx'
identify(x, quiet = FALSE, ...)
```

## Arguments

- x:

  an object of class `networx`

- quiet:

  a logical controlling whether to print a message inviting the user to
  click on the tree.

- ...:

  further arguments to be passed to or from other methods.

## Value

`identify.networx` returns a splits object.

## See also

[`plot.networx`](https://klausvigo.github.io/phangorn/reference/plot.networx.md),
[`identify`](https://rdrr.io/r/graphics/identify.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data(yeast)
dm <- dist.ml(yeast)
nnet <- neighborNet(dm)
plot(nnet)
identify(nnet) # click close to an edge
} # }
```
