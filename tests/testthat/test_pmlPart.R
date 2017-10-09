context("Partitioned likelihood models")

X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
tree2 <- read.tree(text = "((t1:0.3,t3:0.3):0.1,(t2:0.3,t4:0.3):0.1,t5:0.5);")

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
    sp <- pmlPart(edge ~ rate, fit0, weight=W, control = pml.control(trace=0))
    expect_equal( sp$fits[[1]]$rate  / sp$fits[[2]]$rate , 2, tolerance = 1e-5)
    expect_equal( sp$fits[[1]]$rate  / sp$fits[[3]]$rate , 0.5, tolerance = 1e-5)
})    

# nni



# Q
test_that("transition rate optimisation works properly", {
    skip_on_cran()

    Q <- c(6:1)
    
    fit0 <- pml(tree, X, k=4)
    fit1 <- pml(tree, X, k=4, Q=Q)
    weights1 <- 1000*exp(fit1$site) 
    Y <- X 
    attr(Y, "weight") <- weights1
    fit1 <- pml(tree, Y, k=4, Q=Q)
    
    weights0 <- weights1
    weights2 <- weights1

    W = cbind(weights0, weights1, weights2) 
    colnames(W) = c("g1", "g2", "g3")
    
    # linked parameter
    
    sp <- pmlPart(edge + Q ~ ., fit0, weight=W)    
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(Q, sp$fits[[1]]$Q, tolerance=5e-4)
    
    # unlinked parameter
    # TODO more complicated models
    # weights0 <- 1000*exp(fit0$site)
    sp <- pmlPart( ~ Q, fit0, weight=W)   
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(Q, sp$fits[[1]]$Q, tolerance=5e-4)
})


# bf
test_that("base frequency optimisation works properly", {
    skip_on_cran()
    
    bf <- (1:4)/10 
    
    fit0 <- pml(tree, X, k=4)
    fit1 <- pml(tree, X, k=4, bf=bf)
    weights1 <- 1000*exp(fit1$site) 
    Y <- X 
    attr(Y, "weight") <- weights1
    fit1 <- pml(tree, Y, k=4, bf=bf)
    
    weights0 <- weights1
    weights2 <- weights1
    
    W = cbind(weights0, weights1, weights2) 
    colnames(W) = c("g1", "g2", "g3")
    
    # linked parameter
    
    sp <- pmlPart(edge + bf ~ ., fit0, weight=W)    
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(bf, sp$fits[[1]]$bf, tolerance=5e-4)
    
    # unlinked parameter
    # TODO more complicated models
    # weights0 <- 1000*exp(fit0$site)
    sp <- pmlPart( ~ bf, fit0, weight=W)   
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(bf, sp$fits[[1]]$bf, tolerance=5e-4)
})

# Gamma
test_that("shape parameter optimisation works properly", {
    skip_on_cran()
    
    shape=2
    
    fit0 <- pml(tree, X, k=4)
    fit1 <- pml(tree, X, k=4, shape=shape)
    weights1 <- 1000*exp(fit1$site) 
    Y <- X 
    attr(Y, "weight") <- weights1
    fit1 <- pml(tree, Y, k=4, shape=shape)
    
    weights0 <- weights1
    weights2 <- weights1
    
    W = cbind(weights0, weights1, weights2) 
    colnames(W) = c("g1", "g2", "g3")
    
    # linked parameter
  
    sp <- pmlPart(edge + shape ~ ., fit0, weight=W)    
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(shape, sp$fits[[1]]$shape, tolerance=5e-3)
    
    # unlinked parameter
    # TODO more complicated models
    # weights0 <- 1000*exp(fit0$site)
    sp <- pmlPart( ~ shape, fit0, weight=W)   
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(shape, sp$fits[[1]]$shape, tolerance=5e-4)
})

# Invariant sites
test_that("Invariant sites optimisation works properly", {
    skip_on_cran()
    
    inv=.2
    
    fit0 <- pml(tree, X, k=4)
    fit1 <- pml(tree, X, k=4, inv=inv)
    weights1 <- 1000*exp(fit1$site) 
    Y <- X 
    attr(Y, "weight") <- weights1
    fit1 <- pml(tree, Y, k=4, inv=inv)
    
    weights0 <- weights1
    weights2 <- weights1
    
    W = cbind(weights0, weights1, weights2) 
    colnames(W) = c("g1", "g2", "g3")
    
    # linked parameter
    
    sp <- pmlPart(edge + inv ~ ., fit0, weight=W)    
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(inv, sp$fits[[1]]$inv, tolerance=5e-5)
    
    # unlinked parameter
    # TODO more complicated models
    # weights0 <- 1000*exp(fit0$site)
    sp <- pmlPart( ~ inv, fit0, weight=W)   
    expect_equal(logLik(sp)[1], logLik(fit1)[1]*3, tolerance=5e-4 )
    expect_equal(inv, sp$fits[[1]]$inv, tolerance=5e-5)
})


# model




