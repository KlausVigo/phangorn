net <- as.networx(allCircularSplits(5))
outline <- as.networx(allCircularSplits(5), coord = "outline")

test_that("plot.networx works", {
  networx_plot <- function() plot(net, direction="horizontal")
  vdiffr::expect_doppelganger("plot.networx", networx_plot)
})

test_that("plot outline works", {
  outline_plot <- function() plot(outline, direction="horizontal")
  vdiffr::expect_doppelganger("outline_plot", outline_plot)
})



#' @srrstats {EA6.1} *The properties of graphical output from EDA software
#' should be explicitly tested, for example via the
#' [`vdiffr` package](https://github.com/r-lib/vdiffr) or equivalent.*
NULL
