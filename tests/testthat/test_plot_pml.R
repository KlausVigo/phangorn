tree <- read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
tree2 <- read.tree(text = "((t1:1.1,t2:1.1):1,(t3:1,t4:1):1.1);")
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


trees <- c(tree, tree2) |> .compressTipLabel()

test_that("densiTree works", {
  densi_plot <- function() densiTree(trees, type="phylogram", width=2,
          jitter=list(amount=.1, random=FALSE), alpha=1)
  vdiffr::expect_doppelganger("densiTree", densi_plot)
})
