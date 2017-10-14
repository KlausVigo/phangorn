context("modelTest")

X <- allSitePattern(4)
attr(X, "type") <- "DNA"
tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,t3:0.2,t4:0.3);")
fit <- pml(tree, X, k=4, bf = c(1:4)/10)  # F81 + Gamma
weights <- 1000*exp(fit$site) 
attr(X, "weight") <- weights

set.seed(42)
# tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,t3:0.2);")
Y <- simSeq(tree, l=500, type = "AA", model="WAG")

test_that("modelTest works properly", {
    skip_on_cran()
    MT <- modelTest(X, tree = tree, 
                control = pml.control(epsilon = 1e-08, maxit = 10, trace = 0))
    expect_equal(MT$Model[which.min(MT$BIC)], "F81+G")
    # amino acid models
    MT_AA <- modelTest(Y, tree = tree, model=c("JTT", "WAG"), FREQ = TRUE,
            control = pml.control(epsilon = 1e-08, maxit = 10, trace = 0))
    expect_equal(MT_AA$Model[which.min(MT_AA$BIC)], "WAG")
})
