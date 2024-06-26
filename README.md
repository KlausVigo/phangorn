[![R-CMD-check](https://github.com/KlausVigo/phangorn/workflows/R-CMD-check/badge.svg)](https://github.com/KlausVigo/phangorn/actions)
[![CRAN Status Badge](https://www.r-pkg.org/badges/version/phangorn)](https://cran.r-project.org/package=phangorn)
[![CRAN Downloads (monthly)](https://cranlogs.r-pkg.org/badges/phangorn)](https://cran.r-project.org/package=phangorn)
[![CRAN Downloads (total)](https://cranlogs.r-pkg.org/badges/grand-total/phangorn)](https://cran.r-project.org/package=phangorn)
[![Codecov test coverage](https://codecov.io/gh/KlausVigo/phangorn/branch/master/graph/badge.svg)](https://app.codecov.io/gh/KlausVigo/phangorn?branch=master)

# phangorn <img src='man/figures/logo.png' align="right" width="120" />


phangorn is a package for phylogenetic reconstruction and analysis in the R language. phangorn offers the possibility of reconstructing phylogenies with distance based methods, maximum parsimony or maximum likelihood (ML) and performing Hadamard conjugation. Extending the general ML framework, this package provides the possibility of estimating mixture and partition models. Furthermore, phangorn offers several functions for comparing trees, phylogenetic models or splits, simulating character data and performing congruence analyses. 

You can install
- the latest released version `install.packages("phangorn")`
- the latest development version `remotes::install_github("KlausVigo/phangorn")`
- [r-universe](https://r-universe.dev/) kindly provides binaries for Windows, 
Linux and OS X of the development version [here](https://klausvigo.r-universe.dev/phangorn#).

To install the development version you may need to install the Biostrings package from bioconductor first:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
```
Also the development version usually depends on the latest ape development 
version and information to download can be found 
[here](https://emmanuelparadis.github.io/ape_installation.html). 
Additionally you may need to install on windows [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and on mac [XCode](https://developer.apple.com/xcode/)
and [GFortran](https://gcc.gnu.org/wiki/GFortranBinaries).

If you use phangorn please cite:

Schliep K.P. 2011. phangorn: phylogenetic analysis in R. Bioinformatics, 27(4) 592-593 


License
-------
phangorn is licensed under the GPLv2.
