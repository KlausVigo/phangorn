net <- as.networx(allCircularSplits(5))

test_that("plot.networx works", {
  networx_plot <- function() plot(net)
  vdiffr::expect_doppelganger("plot.networx", networx_plot)
})


#' @srrstats {EA6.1} *The properties of graphical output from EDA software
#' should be explicitly tested, for example via the
#' [`vdiffr` package](https://github.com/r-lib/vdiffr) or equivalent.*
NULL
