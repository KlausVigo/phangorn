net <- as.networx(allCircularSplits(5))

test_that("plot.networx works", {
  networx_plot <- function() plot(net)
  vdiffr::expect_doppelganger("plot.networx", networx_plot)
})
