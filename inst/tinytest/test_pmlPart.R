X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
tree2 <- read.tree(text = "((t1:0.3,t3:0.3):0.1,(t2:0.3,t4:0.3):0.1,t5:0.5);")

fit0 <- pml(tree, X)
fit1 <- update(fit0, rate=0.5)
fit2 <- update(fit0, rate=2)

weights0 <- 1000*exp(fit0$siteLik)
weights1 <- 1000*exp(fit1$siteLik)
weights2 <- 1000*exp(fit2$siteLik)

W <- cbind(weights0, weights1, weights2)
colnames(W) <- c("g1", "g2", "g3")


# test_that rate optimisation works properly
sp <- pmlPart(edge ~ rate, fit0, weight=W, control = pml.control(trace=0))
expect_equal( sp$fits[[1]]$rate / sp$fits[[2]]$rate , 2, tolerance = 1e-5)
expect_equal( sp$fits[[1]]$rate / sp$fits[[3]]$rate , 0.5, tolerance = 1e-5)


# extract trees
trees <- pmlPart2multiPhylo(sp)
expect_true(inherits(trees, "multiPhylo"))

# test multiphyDat objects
XX <- as.character(X)
p0 <- sample(1024, 1000, replace = TRUE, prob = exp(fit0$siteLik))
X0 <- phyDat(XX[,p0])
p1 <- sample(1024, 1000, replace = TRUE, prob = exp(fit1$siteLik))
X1 <- phyDat(XX[,p1])
p2 <- sample(1024, 1000, replace = TRUE, prob = exp(fit2$siteLik))
X2 <- phyDat(XX[,p2])

if(suppressPackageStartupMessages(requireNamespace('apex'))){
  mp <- new("multiphyDat", list(X0, X1, X2))
  sp_mp <- pmlPart(~ rate, mp, control = pml.control(trace = 0))
  expect_true(inherits(sp_mp, "pmlPart"))
  expect_true(sp_mp$fits[[1]]$rate > sp_mp$fits[[2]]$rate)
  expect_true(sp_mp$fits[[1]]$rate < sp_mp$fits[[3]]$rate)
}


# nni


# test_that transition rate optimisation works properly
Q <- c(6:1)

fit0 <- pml(tree, X, k=4)
fit1 <- pml(tree, X, k=4, Q=Q)
weights1 <- 1000*exp(fit1$siteLik)
Y <- X
attr(Y, "weight") <- weights1

fit1 <- pml(tree, Y, k=4, Q=Q)
weights0 <- weights1
weights2 <- weights1

W <- cbind(weights0, weights1, weights2)
colnames(W) <- c("g1", "g2", "g3")

# linked parameter

sp <- pmlPart(edge + Q ~ ., fit0, weight=W, control = pml.control(trace=0))
expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
expect_equal(Q, sp$fits[[1]]$Q, tolerance=5e-4)

# unlinked parameter
    sp <- pmlPart( ~ Q, fit0, weight=W, control = pml.control(trace=0))
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(Q, sp$fits[[1]]$Q, tolerance=5e-4)


# test_that base frequency optimisation works properly
    bf <- (1:4)/10

    fit0 <- pml(tree, X, k=4)
    fit1 <- pml(tree, X, k=4, bf=bf)
    weights1 <- 1000*exp(fit1$siteLik)
    Y <- X
    attr(Y, "weight") <- weights1
    fit1 <- pml(tree, Y, k=4, bf=bf)

    weights0 <- weights1
    weights2 <- weights1

    W <- cbind(weights0, weights1, weights2)
    colnames(W) <- c("g1", "g2", "g3")

    # linked parameter

    sp <- pmlPart(edge + bf ~ ., fit0, weight=W, control = pml.control(trace=0))
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(bf, sp$fits[[1]]$bf, tolerance=5e-4)

    # unlinked parameter
    sp <- pmlPart( ~ bf, fit0, weight=W, control = pml.control(trace=0))
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(bf, sp$fits[[1]]$bf, tolerance=5e-4)

# Gamma
# test_that shape parameter optimisation works properly
    shape <- 2

    fit0 <- pml(tree, X, k=4)
    fit1 <- pml(tree, X, k=4, shape=shape)
    weights1 <- 1000*exp(fit1$siteLik)
    Y <- X

    attr(Y, "weight") <- weights1
    fit1 <- pml(tree, Y, k=4, shape=shape)

    weights0 <- weights1
    weights2 <- weights1

    W <- cbind(weights0, weights1, weights2)
    colnames(W) <- c("g1", "g2", "g3")

    # linked parameter

    sp <- pmlPart(edge + shape ~ ., fit0, weight=W,
                  control=pml.control(trace=0))
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(shape, sp$fits[[1]]$shape, tolerance=5e-3)

    # unlinked parameter
    sp <- pmlPart( ~ shape, fit0, weight=W, control=pml.control(trace=0))
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(shape, sp$fits[[1]]$shape, tolerance=5e-4)

# Invariant sites
# test_that Invariant sites optimisation works properly
    inv <- 0.2

    fit0 <- pml(tree, X, k=4)
    fit1 <- pml(tree, X, k=4, inv=inv)
    weights1 <- 1000*exp(fit1$siteLik)
    Y <- X
    attr(Y, "weight") <- weights1
    fit1 <- pml(tree, Y, k=4, inv=inv)

    weights0 <- weights1
    weights2 <- weights1

    W <- cbind(weights0, weights1, weights2)
    colnames(W) <- c("g1", "g2", "g3")

    # linked parameter

    sp <- pmlPart(edge + inv ~ ., fit0, weight=W, control=pml.control(trace=0))
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(inv, sp$fits[[1]]$inv, tolerance=5e-5)

    # unlinked parameter
    sp <- pmlPart( ~ inv, fit0, weight=W, control = pml.control(trace=0))
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(inv, sp$fits[[1]]$inv, tolerance=5e-4)



# linked parameters
# test_that Linked parameters optimisation works properly
    Z <- X

    fit0 <- pml(tree, X, k=4)
    weights0 <- 1000*exp(fit0$siteLik)
    weights1 <- 1000*exp(update(fit0, rate=0.5)$siteLik)
    weights2 <- 1000*exp(update(fit0, tree=tree2)$siteLik)

    attr(Z, "weight") <- weights0 + weights1 + weights2
    W <- cbind(weights0, weights1, weights2)

    fit_Z <- update(fit0, data=Z)
    fit_Z <- optim.pml(fit_Z, model="GTR", rearrangement="NNI", optGamma=TRUE,
                   optInv=TRUE, control=pml.control(trace=0))
    sp <- pmlPart(edge + bf + Q + shape + inv + nni ~ ., fit_Z, weight=W,
                  control = pml.control(trace=0))

    expect_equal(sp$logLik[[1]], fit_Z$logLik, tolerance = 1e-5)

