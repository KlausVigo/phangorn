# Parallel computing with phangorn

### Changes in phangorn 3.0

Several function in phangorn make use of parallel computing. With the
new version 3 we move to using the *future* and *future.apply* package
(Bengtsson 2021).

Before we were using the `mclapply()` function in in
[`bootstrap.pml()`](https://klausvigo.github.io/phangorn/reference/bootstrap.pml.md),
[`bootstrap.phyDat()`](https://klausvigo.github.io/phangorn/reference/bootstrap.pml.md),
[`modelTest()`](https://klausvigo.github.io/phangorn/reference/modelTest.md).
The disadvantage of `mclapply()` is that it works not for Windows and
its usage is not recommended inside a GUI like RStudio. The future
package allows the user to define its backend and can so be adjusted to
the hardware and operating system of the user. Now also on Windows one
can use parallel processing. It is very likely that several other
function will allow the future framework in future. By default processes
are running sequential. To run a function in parallel we just need to
call the function plan specifying the backend and than the function
which allows the parallel processing.

``` r
library(phangorn)
library(future)
data("Laurasiatherian")
plan(multisession, workers = 2)
mt <- modelTest(Laurasiatherian, model = c("JC", "F81", "K80", "HKY", "SYM", "GTR"))
plan(sequential) # run sequential again
```

It is good practice to clean up and set the backend to sequential after
using parallel processing.

A few things to consider. When running processes in parallel the memory
footprint might be larger compared to the sequential computation.
Furthermore some of the linear algebra routines build in R might
implicitly use parallel computing. So check how many cores are in use.
If too many cores are used parallel code might even slow down.

For more see the man page for plan
[`help(plan)`](https://future.futureverse.org/reference/plan.html), the
vignettes `vignette(package="future")` and other documentation of the
future package.

## Session info

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    ##  [5] xfun_0.56         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
    ##  [9] rmarkdown_2.30    lifecycle_1.0.5   cli_3.6.5         sass_0.4.10      
    ## [13] pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4   systemfonts_1.3.1
    ## [17] compiler_4.5.2    tools_4.5.2       ragg_1.5.0        bslib_0.10.0     
    ## [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
    ## [25] rlang_1.1.7       fs_1.6.6          htmlwidgets_1.6.4

## References

Bengtsson, Henrik. 2021. “A Unifying Framework for Parallel and
Distributed Processing in r Using Futures.” *The R Journal* 13 (2):
208–27. <https://doi.org/10.32614/RJ-2021-048>.
