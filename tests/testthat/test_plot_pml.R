tree <- read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
dat <- matrix(c("a", "a",
                "a", "t",
                "t", "a",
                "t", "t"), byrow = TRUE, nrow = 4L,
              dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
dna <- phyDat(dat)
fit <- pml(tree, dna)

test_that("plot.pml works", {
  pml_plot <- function() plot(fit)
  vdiffr::expect_doppelganger("plot.pml", pml_plot)
})
