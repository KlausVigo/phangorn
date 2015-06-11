context("parsimony")

tree1 = read.tree(text = "((t1,t2),t3,t4);")
tree2 = read.tree(text = "((t1,t3),t2,t4);")
dat <- phyDat(c(t1="a", t2="a",t3="t",t4="t"), type="USER", levels=c("a","c","g","t"))
#tr_acctran = acctran(tree1, dat)
#tr_ratchet = pratchet(dat, trace=0)
#bab(dat)
test_that("parsimony works properly", {
    skip_on_cran()
    expect_that(fitch(tree1, dat), equals(1))
    expect_that(fitch(tree2, dat), equals(2))
    expect_that(sankoff(tree1, dat), equals(1))
    expect_that(sankoff(tree2, dat), equals(2))
    expect_that(parsimony(tree1, dat), equals(1))
})





