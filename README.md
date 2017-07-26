[![Build Status](https://travis-ci.org/KlausVigo/phangorn.svg?branch=master)](https://travis-ci.org/KlausVigo/phangorn)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/phangorn)](https://cran.r-project.org/package=phangorn)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/phangorn)](https://cran.r-project.org/package=phangorn)
[![Research software impact](http://depsy.org/api/package/cran/phangorn/badge.svg)](http://depsy.org/package/r/phangorn)
[![codecov.io](https://codecov.io/github/KlausVigo/phangorn/coverage.svg?branch=master)](https://codecov.io/github/KlausVigo/phangorn?branch=master)


phangorn
========================================================

phangorn is a package for phylogenetic reconstruction and analysis in the R language. phangorn offers the possibility of reconstructing phylogenies with distance based methods, maximum parsimony or maximum likelihood (ML) and performing Hadamard conjugation. Extending the general ML framework, this package provides the possibility of estimating mixture and partition models. Furthermore, phangorn offers several functions for comparing trees, phylogenetic models or splits, simulating character data and performing congruence analyses. 

You can install
- the latest released version `install.packages("phangorn")`
- the latest development version `devtools::install_github("KlausVigo/phangorn")`

To install the development version you may need to install the Biostrings and seqLogo package from bioconductor first:
```
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings", "seqLogo")
```

If you use phangorn please cite:

Schliep K.P. 2011. phangorn: phylogenetic analysis in R. Bioinformatics, 27(4) 592-593 


License
-------
phangorn is licensed under the GPLv2.
