# Arithmetic Operators

double factorial function

## Usage

``` r
ldfactorial(x)

dfactorial(x)
```

## Arguments

- x:

  a numeric scalar or vector

## Value

`dfactorial(x)` returns the double factorial, that is \\x\\\\ = 1 \* 3
\* 5 \* \ldots \* x \\ and `ldfactorial(x)` is the natural logarithm of
it.

## See also

[`factorial`](https://rdrr.io/r/base/Special.html),
[`howmanytrees`](https://rdrr.io/pkg/ape/man/howmanytrees.html)

## Author

Klaus Schliep <klaus.schliep@gmail.com>

## Examples

``` r
dfactorial(1:10)
#>  [1]    1.000000    1.595769    3.000000    6.383076   15.000000   38.298459
#>  [7]  105.000000  306.387671  945.000000 3063.876713
```
