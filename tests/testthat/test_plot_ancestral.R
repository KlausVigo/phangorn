tree <- read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
tree2 <- read.tree(text = "(((t1:1,t2:1):1,t3:2):1,t4:3);")
dat <- matrix(c("a", "a",
                "a", "t",
                "t", "a",
                "t", "t"), byrow = TRUE, nrow = 4L,
              dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
dna <- phyDat(dat)
fit <- pml(tree, dna)
test_ml <- anc_pml(fit)


test_that("Pie_plots works", {
  pie_ancestral <- function() plotAnc(test_ml)
  vdiffr::expect_doppelganger("Pie plots", pie_ancestral)
})


test_that("plotSeqLogo works", {
  seq_logo <- function() plotSeqLogo(test_ml)
  vdiffr::expect_doppelganger("SeqLogo", seq_logo)
})
