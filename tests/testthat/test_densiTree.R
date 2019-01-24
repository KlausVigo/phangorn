context("densiTree")

test_that("minimal use", {
  phylogeny_1 <- ape::read.tree(text = "((A:2, B:2):1, C:3);")
  phylogeny_2 <- ape::read.tree(text = "((A:1, B:1):2, C:3);")
  trees <- c(phylogeny_1, phylogeny_2)

  # Fails, with error 'Error in temp[[j + 1]] : subscript out of bounds'
  testthat::expect_silent(densiTree(trees))
})

# Visual tests ------------------------------------------------------------
test_that("visual appearance", {
    skip_on_cran()
    data(Laurasiatherian)
    set.seed(1)
    bs_trees <- bootstrap.phyDat(subset(Laurasiatherian, 1:10), FUN = 
                                     function(x) upgma(dist.hamming(x)), bs=5)
    densiT <- function() densiTree(bs_trees, type="phylogram", scale.bar=FALSE, 
                        width=2, jitter=list(amount=.3, random=FALSE), alpha=1)
    
    vdiffr::expect_doppelganger("densiTree plot", densiT)
})
