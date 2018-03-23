context("densiTree")

data(Laurasiatherian)
set.seed(1)
bs_trees <- bootstrap.phyDat(Laurasiatherian, FUN = 
                   function(x) upgma(dist.hamming(x)), bs=10)

densiT <- function() densiTree(bs_trees, type="phylogram", scale.bar=FALSE)

# Visual tests ------------------------------------------------------------
test_that("visual appearance", {
    skip_on_cran()
    vdiffr::expect_doppelganger("densiTree plot", densiT)
})
