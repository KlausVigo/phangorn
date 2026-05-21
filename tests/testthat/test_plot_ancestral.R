#tree <- read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
tree <- read.tree(text = "(((t1:1,t2:1):1,t3:2):1,(t4:1.5, t5:1.5):1.5);")
dat <- matrix(c("a", "c",
                "a", "g",
                "t", "c",
                "t", "g",
                "t", "g"), byrow = TRUE, nrow = 5L,
              dimnames = list(c("t1", "t2", "t3", "t4", "t5"), NULL))
dna <- phyDat(dat)
fit <- pml(tree, dna)
test_ml <- anc_pml(fit)



test_that("Pie_plots works", {
  pie_ancestral <- function() {
    plotAnc(test_ml, pos="bottomleft")
    add_mutations(test_ml, adj = c(0.5, -0.3))
  }
  vdiffr::expect_doppelganger("Pie plots", pie_ancestral)
})

#old_par <- par(mar=c(5,4,4,4))
#test_that("anc_heatmap works", {
#  heat <- function() anc_heatmap(test_ml)
#  vdiffr::expect_doppelganger("anc_heatmap", heat)
#})
#par(old_par)

#test_that("plotSeqLogo works", {
#  p <- plotSeqLogo(test_ml)
#  vdiffr::expect_doppelganger("SeqLogo", p)
#})
