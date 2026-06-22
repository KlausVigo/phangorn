ON_CRAN <- !identical(Sys.getenv("R_NOT_CRAN"), "true")
ON_GH <- identical(Sys.getenv("R_GH"), "true")
ON_CI <- isTRUE(ON_CRAN) || isTRUE(ON_GH)
ON_WINDOWS <- isTRUE(Sys.info()[["sysname"]] == "Windows")
ON_OSX <- isTRUE(Sys.info()[["sysname"]] == "Darwin")

library(tinytest)
Sys.setenv(LANGUAGE = "en") # Force locale
