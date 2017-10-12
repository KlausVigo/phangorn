context("modelTest")

X <- allSitePattern(4)
attr(X, "type") <- "DNA"
tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,t3:0.2,t4:0.3);")
fit <- pml(tree, X, k=4, bf = c(1:4)/10)  # F81 + Gamma
weights <- 1000*exp(fit$site) 
attr(X, "weight") <- weights


test_that("modelTest works properly", {
    skip_on_cran()
    MT <- modelTest(X, tree = tree, 
                control = pml.control(epsilon = 1e-08, maxit = 10, trace = 0))
    expect_equal(MT$Model[which.min(MT$BIC)], "F81+G")
})