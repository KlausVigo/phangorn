context("lento plot")



# Visual tests ------------------------------------------------------------
test_that("visual appearance", {
    data(yeast)
    yeast.ry <- acgt2ry(yeast)
    splits.h <- h2st(yeast.ry)
    lentoPlot <- function() lento(splits.h)
    skip_on_cran()
    vdiffr::expect_doppelganger("lento plot", lentoPlot)
})
