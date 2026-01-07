# phangorn

phangorn is a package for phylogenetic reconstruction and analysis in
the R language. phangorn offers the possibility of reconstructing
phylogenies with distance based methods, maximum parsimony or maximum
likelihood (ML) and performing Hadamard conjugation. Extending the
general ML framework, this package provides the possibility of
estimating mixture and partition models. Furthermore, phangorn offers
several functions for exploring and comparing trees, phylogenetic models
or splits, simulating character data and performing congruence analyses.

To get an introduction into phylogenetic inference you want to look at:

``` R
vignette("Trees", package="phangorn")
```

## Installation

You can install the the latest release `phangorn` of the package from
[CRAN](https://CRAN.R-project.org/package=phangorn), or the development
version from [github](https://github.com/KlausVigo/phangorn) or
[r-universe](https://klausvigo.r-universe.dev/phangorn).

| Type        | Source     | Command                                                                                                  |
|-------------|------------|----------------------------------------------------------------------------------------------------------|
| Release     | CRAN       | `install.packages("phangorn")`                                                                           |
| Development | GitHub     | `remotes::install_github("KlausVigo/phangorn")`                                                          |
| Development | r-universe | `install.packages('phangorn', repos = c('https://klausvigo.r-universe.dev', 'https://repo.r-wasm.org'))` |

To install the development version you may need to install the
Biostrings package from bioconductor first:

``` R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
```

To use all functionality you might install all you might need to install
the `rgl` and `ggseqlogo` package.

``` R
install.packages("rgl")
install.packages("ggseqlogo")
```

The development version usually depends on the latest `ape` development
version and information to download can be found
[here](https://emmanuelparadis.github.io/ape_installation.html).

## Citation

If you use phangorn please cite:

Schliep K.P. 2011. phangorn: phylogenetic analysis in R. Bioinformatics,
27(4) 592-593

## License

phangorn is licensed under the GPLv2.
