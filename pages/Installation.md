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

or install the latest development version using the `install_github` function in 
the [devtools](https://github.com/hadley/devtools) package from github.
    
    library(devtools)
    install_github("KlausVigo/phangorn")

For devtools to work on windows you need addionally to have installed [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and on mac you need [Xcode](https://developer.apple.com/xcode/) to compile some C code. Once these 
programs are installed you can install phangorn with the commands above.
The development version of phangorn depends on the latest deveolpment version 
from the ape package. Instructions to for download and installation of ape can 
be found [here](http://ape-package.ird.fr/ape_installation.html).
    

    