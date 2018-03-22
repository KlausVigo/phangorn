context("densiTree")

data(Laurasiatherian)
set.seed(1)
bs_trees <- bootstrap.phyDat(Laurasiatherian, FUN = 
                   function(x) upgma(dist.hamming(x)), bs=25)

densiT <- function() densiTree(bs_trees, type="cladogram", col="blue")

# Visual tests ------------------------------------------------------------
test_that("visual appearance", {
    testthat::skip_on_cran()
    vdiffr::expect_doppelganger("densiTree plot", densiT)
})
