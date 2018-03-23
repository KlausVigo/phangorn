context("densiTree")

data(Laurasiatherian)
set.seed(1)
bs_trees <- bootstrap.phyDat(subset(Laurasiatherian, 1:10), FUN = 
                   function(x) upgma(dist.hamming(x)), bs=5)

densiT <- function() densiTree(bs_trees, type="phylogram", scale.bar=FALSE, 
                     width=2, jitter=list(amount=.3, random=FALSE), alpha=1)


# Visual tests ------------------------------------------------------------
test_that("visual appearance", {
    skip_on_cran()
    vdiffr::expect_doppelganger("densiTree plot", densiT)
})
