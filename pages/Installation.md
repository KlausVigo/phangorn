---
layout: page
title: Installing phangorn
description: Instructions on how to install the phangorn package.
---

phangorn is a package for [R](http://www.r-project.org). To use it,
you first need to [download](http://cran.r-project.org/) and install
R. [RStudio](http://www.rstudio.com/) provides a nice user interface for R.

You can install the latest release version from CRAN

    install.packages("phangorn")

Then install R/qtlcharts using the `install_github` function in the
[devtools](https://github.com/hadley/devtools) package.  

    library(devtools)
    install_github("KlausVigo/phangorn")
    