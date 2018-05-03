context("lento plot")



# Visual tests ------------------------------------------------------------
test_that("visual appearance", {
    skip_on_cran()
    data(yeast)
    yeast.ry <- acgt2ry(yeast)
    splits.h <- h2st(yeast.ry)
    lentoPlot <- function() lento(splits.h)
    vdiffr::expect_doppelganger("lento plot", lentoPlot)
})
