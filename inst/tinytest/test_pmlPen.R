X <- allSitePattern(5)
attr(X, "index") <- NULL
tree1 <- read.tree(text = "((t1:0.1,t2:0.5):0.1,(t3:0.1,t4:0.5):0.1,t5:0.5);")
tree2 <- read.tree(text = "((t1:0.5,t2:0.1):0.1,(t3:0.5,t4:0.1):0.1,t5:0.5);")

fit1 <- pml(tree1,X)
fit2 <- pml(tree2,X)

attr(X, "weight") <- 1000*exp(fit1$siteLik)
Y <- X
attr(Y, "weight") <- 1000*exp(fit2$siteLik)

fit1 <- update(fit1, data=X)
fit2 <- update(fit2, data=Y)

sp <- pmlPart(~ edge, list(fit1, fit2), pml.control(trace = 0))
pp0 <- pmlPen(sp, lambda = 0, pml.control(trace = 0))
ppInf <- pmlPen(sp, lambda = 1e6, pml.control(trace = 0))

Z <- cbind(X,Y)
fit3 <- update(fit1, data=Z)
fit3 <- optim.pml(fit3, control = pml.control(trace = 0))

# test penalized partition model
expect_equal(pp0$logLik[1], fit1$logLik + fit2$logLik) #, tolerance = 0.002)
expect_equal(ppInf$logLik[1], fit3$logLik, tolerance = 1e-5)


#

#
# simulation of mixture models
#
X <- allSitePattern(5)
tree1 <- read.tree(text = "((t1:0.1,t2:0.5):0.1,(t3:0.1,t4:0.5):0.1,t5:0.5);")
tree2 <- read.tree(text = "((t1:0.5,t2:0.1):0.1,(t3:0.5,t4:0.1):0.1,t5:0.5);")
tree1 <- unroot(tree1)
tree2 <- unroot(tree2)
fit1 <- pml(tree1,X)
fit2 <- pml(tree2,X)

weights <- 2000*exp(fit1$siteLik) + 1000*exp(fit2$siteLik)
attr(X, "weight") <- weights

ll0 <- sum(weights * log(weights / sum(weights)))


fit1 <- pml(tree1, X)
fit2 <- optim.pml(fit1)
logLik(fit2)
AIC(fit2, k=log(3000))

fitMixEdge <- pmlMix( ~ edge, fit1, m=2)
logLik(fitMixEdge)
AIC(fitMixEdge, k=log(3000))

fit.p <- pmlPen(fitMixEdge, .01, control=pml.control(trace = 0))
expect_equal(fit.p$logLik, ll0, tolerance = 0.4)

