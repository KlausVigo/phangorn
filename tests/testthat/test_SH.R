context("SH-test")

tree <- read.tree(text="(t1:0.1,t2:0.1,(t3:0.1,(t4:0.1,t5:0.1):0.1):0.1);")
set.seed(42)
dat <- simSeq(tree, l=500)

trees <- c(tree, nni(tree))
fits <- lapply(trees, pml, data=dat)

test_that("SH-test works properly", {
    skip_on_cran()
    tmp <- SH.test(fits[[1]], fits[[2]])
    expect_gt(tmp[1,"p-value"], 0)
    expect_lt(tmp[2,"p-value"], 0.05)
    tmp <- SH.test(fits)
    expect_gt(tmp[1,"p-value"], 0)
    expect_lt(tmp[2,"p-value"], 0.05)
    expect_lt(tmp[3,"p-value"], 0.05)
})