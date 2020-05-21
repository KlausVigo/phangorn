if (require("tinytest", quietly=TRUE)){
    set.seed(42)
    home <- length(unclass(packageVersion("phangorn"))[[1]]) == 4 # 0.0.0.9000
    if(home) {
        # Only run locally, let CRAN test examples and the vignette
        test_package("phangorn", at_home = home)
    }
}
