context("Partitioned likelihood models")

X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
fit0 <- pml(tree, X, k=4)

fit1 <- update(fit0, rate=.5)
fit2 <- update(fit0, rate=2)

weights0 <- 1000*exp(fit0$site) 
weights1 <- 1000*exp(fit1$site) 
weights2 <- 1000*exp(fit2$site) 


W = cbind(weights0, weights1, weights2) 
colnames(W) = c("g1", "g2", "g3")


# rate
test_that("rate optimisation works properly", {
    skip_on_cran()
    sp <- pmlPart(edge ~ rate, fit0, weight=W)
    expect_equal( sp$fits[[1]]$rate  / sp$fits[[2]]$rate , 2, tolerance = 1e-5)
    expect_equal( sp$fits[[1]]$rate  / sp$fits[[3]]$rate , 0.5, tolerance = 1e-5)
})    

# nni

# bf + Q

# model

# rate
