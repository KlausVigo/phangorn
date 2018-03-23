context("lento plot")

data(yeast)
yeast.ry <- acgt2ry(yeast)
splits.h <- h2st(yeast.ry)

lentoPlot <- function() lento(splits.h, trivial=TRUE)

# Visual tests ------------------------------------------------------------
test_that("visual appearance", {
    skip_on_cran()
    vdiffr::expect_doppelganger("lento plot", lentoPlot)
})
