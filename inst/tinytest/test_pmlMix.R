X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
tree <- reorder(tree, "postorder")
fit <- pml(tree,X, k=4)
weights <- 1000*exp(fit$siteLik)

ll0 <- sum(weights*log(weights/sum(weights)))

attr(X, "weight") <- weights
fit1 <- update(fit, data=X, k=1)

#fit2 <- update(fit, data=X)
#(fit2 <- optim.pml(fit2, optGamma=TRUE))

# test rate optimisation works properly
fitMixture <- pmlMix(~rate, fit1 , m=4, control=pml.control(trace=0))
expect_equal(fitMixture$logLik, ll0, tolerance = 1e-4)
rates <- sapply(fitMixture$fits, \(x)x$rate )
expect_true( cor(rates, discrete.gamma(1,4)) > .95)

# test edge optimisation works properly
tree_tmp <- tree
tree_tmp$edge.length[] <- .3
rates <- discrete.gamma(1,4)
fits <- list()
for(i in 1:4) fits[[i]] <- update(fit1, tree=tree_tmp, rate=rates[i])
fitMixture <- pmlMix(edge~., fits, m=4, control=pml.control(trace=0))
expect_true( cor(fitMixture$fits[[1]]$tree$edge.length, tree$edge.length) > .99)

