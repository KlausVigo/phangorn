context("distances")

X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
fit <- pml(tree,X, k=4, shape=0.5)

weights <- as.vector(1000*exp(fit$site)) 
attr(X, "weight") <- weights
dm <- cophenetic(tree)


test_that("dist.ml works properly", {
    skip_on_cran()
    expect_that(dist.logDet(X), is_a("dist"))
    expect_that(dist.hamming(X), is_a("dist"))
    expect_that(dist.ml(X), is_a("dist"))
    all.equal(as.matrix(dist.ml(X, k=4, shape=.5)), dm)
})