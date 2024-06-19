X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
tree <- reorder(tree, "postorder")
tree2 <- read.tree(text = "((t1:0.3,t3:0.3):0.1,(t2:0.3,t4:0.3):0.1,t5:0.5);")
tree2 <- reorder(tree2, "postorder")
fit <- pml(tree,X, k=4)
weights <- 1000*exp(fit$siteLik)

ll0 <- sum(weights*log(weights/sum(weights)))

attr(X, "weight") <- weights
fit1 <- update(fit, data=X, k=1)


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

for(i in 1:4) fits[[i]] <- update(fit1, tree=tree_tmp)
fitMixture <- pmlMix(edge~rate, fits, m=4, control=pml.control(trace=0))
expect_true( cor(fitMixture$fits[[1]]$tree$edge.length, tree$edge.length) > .99)
rates <- sapply(fitMixture$fits, \(x)x$rate )
expect_true( cor(rates, discrete.gamma(1,4)) > .95)


fit1 <- pml(tree, X)
fit2 <- pml(tree2, X)
weights <- 600*exp(fit1$siteLik) + 400*exp(fit2$siteLik)
attr(X, "weight") <- weights
tmp_tree1 <- tree
tmp_tree2 <- tree2
tmp_tree1$edge.length[] <- tmp_tree1$edge.length[] <- 1
fits <- list()
fits[[1]] <- pml(tmp_tree1, X)
fits[[2]] <- pml(tmp_tree2, X)
fitMixture <- pmlMix( ~ edge, fits, m=2, control=pml.control(trace=0))



# test tree rearrangements work properly, needs bug fixes
#fit1 <- pml(tree, X)
#fit2 <- pml(tree2, X)
#weights <- 500*exp(fit1$siteLik) + 500*exp(fit2$siteLik)
#attr(X, "weight") <- weights
#fits <- list()
#for(i in 1:2) fits[[i]] <- fit1
#fitMixture <- pmlMix( ~ nni, fits, m=2, control=pml.control(trace=0))


# test rate matrix optimization works properly
fit <- pml(tree, X, k=4, Q=6:1)
weights <- 1000*exp(fit$siteLik)
attr(X, "weight") <- weights
fit <- pml(tree, X, k=4, Q=6:1)
fits <- list()
for(i in 1:4) fits[[i]] <-  pml(tree, X, rate=rates[i])
fitMixture <- pmlMix(Q ~ ., fits, m=4, control=pml.control(trace=0))
Q_est <- fitMixture$fits[[1]]$Q
expect_true( cor(Q_est, 6:1) > .999)


# test base frequency optimization works properly
fit <- pml(tree, X, k=4, bf=(1:4)/10)
weights <- 1000*exp(fit$siteLik)
attr(X, "weight") <- weights
fits <- list()
for(i in 1:4) fits[[i]] <-  pml(tree, X, rate=rates[i])
fitMixture <- pmlMix(bf ~ ., fits, m=4, control=pml.control(trace=0))
bf_est <- fitMixture$fits[[1]]$bf
expect_true( cor(bf_est, (1:4)/10) > .999)


# test base frequency optimization works properly
fit1 <- pml(tree, X, bf=(1:4)/10)
fit2 <- pml(tree, X, bf=(4:1)/10)
weights <- 600*exp(fit1$siteLik) + 400*exp(fit2$siteLik)
attr(X, "weight") <- weights
fits <- list()
for(i in 1:2) fits[[i]] <- pml(tree, X)
fitMixture <- pmlMix( ~ bf, fits, m=2, control=pml.control(trace=0))
bf1_est <- fitMixture$fits[[1]]$bf
bf2_est <- fitMixture$fits[[2]]$bf
if(bf1_est[1] < bf1_est[2]) {
  expect_true( cor(bf1_est, (1:4)/10) > .999)
} else expect_true( cor(bf1_est, (4:1)/10) > .999)



# test invariant site optimization works properly
#rates <- discrete.gamma(1,2)
fit1 <- pml(tree, X, inv=.3)
weights <- 1000*exp(fit1$siteLik)
attr(X, "weight") <- weights
fits <- list()
for(i in 1:2) fits[[i]] <- pml(tree, X)
fitMixture <- pmlMix(inv ~ ., fits, m=2, control=pml.control(trace=0))
inv1_est <- fitMixture$fits[[1]]$inv
inv2_est <- fitMixture$fits[[2]]$inv
expect_equal( inv1_est, .3, 1e-3)
expect_equal( inv2_est, .3, 1e-3)


# Not identifiable, so only check logLik
fit1 <- pml(tree, X, inv=0.3)
fit2 <- pml(tree, X, inv=.0)
weights <- 500*exp(fit1$siteLik) + 500*exp(fit2$siteLik)
attr(X, "weight") <- weights
ll <- sum(weights * log(fit1$lv *.5 + fit2$lv*.5))
fits <- list()
for(i in 1:2) fits[[i]] <- pml(tree, X, wMix=.5)
fitMixture <- pmlMix( ~ inv, fits, m=2,
                      control=pml.control(maxit = 25, trace=0))
expect_equal(logLik(fitMixture)[1], ll, 1e-3)
