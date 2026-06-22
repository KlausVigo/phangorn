source("helpers.R")
if (ON_WINDOWS || ON_OSX) exit_file("plot on Linux only")
using("tinysnapshot")

net <- as.networx(allCircularSplits(5))
outline <- as.networx(allCircularSplits(5), coord = "outline")

networx_plot <- function() plot(net, direction="horizontal")
expect_snapshot_plot(networx_plot, "networx_plot")

outline_plot <- function() plot(outline, direction="horizontal")
expect_snapshot_plot(outline_plot, "outline_plot")

#' @srrstats {EA6.1} *The properties of graphical output from EDA software
#' should be explicitly tested, for example via the
#' [`vdiffr` package](https://github.com/r-lib/vdiffr) or equivalent.*
NULL
