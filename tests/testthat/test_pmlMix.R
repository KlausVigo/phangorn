context("Mixture models")

X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
fit <- pml(tree,X, k=4)
weights <- 1000*exp(fit$site) 

ll0 <- sum(weights*log(weights/sum(weights)))

attr(X, "weight") <- weights
fit1 <- update(fit, data=X, k=1)

#fit2 <- update(fit, data=X)
#(fit2 <- optim.pml(fit2, optGamma=TRUE))

test_that("rate optimisation works properly", {
    skip_on_cran()
    fitMixture <- pmlMix(~rate, fit1 , m=4, control=pml.control(trace=0))
    expect_equal(fitMixture$logLik, ll0)
})    
