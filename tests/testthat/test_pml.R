context("likelihood")

treeU1 <- read.tree(text = "((t1:.1,t2:.1):.1,t3:.2,(t4:.1,t5:.1):.3);")
# edge length changed
treeU2 <- read.tree(text = "((t1:.15,t2:.15):.05,t3:.2,(t4:.1,t5:.1):.3);") 
# 1 NNI move
treeU3 <- read.tree(text = "((t1:.1,t3:.1):.1,t2:.2,(t4:.1,t5:.1):.3);") 

treeR1 <- read.tree(text = "(((t1:.1,t2:.1):.1,t3:.2):.1,(t4:.1,t5:.1):.2);")
# edge length changed
treeR2 <- read.tree(text = "(((t1:.15,t2:.15):.05,t3:.2):.1,(t4:.1,t5:.1):.2);") 
# 1 NNI move
treeR3 <- read.tree(text = "(((t1:.1,t3:.1):.1,t2:.2):.1,(t4:.1,t5:.1):.2);")


# dat <- phyDat(c(t1="a", t2="a",t3="t",t4="t"), type="USER", 
#                                  levels=c("a","c","g","t"))
#tree2 <- read.tree(text = "((t1,t3),t2,t4);")

dat <- allSitePattern(5)
weights <- as.vector(1000 * exp(pml(treeR1, dat)$siteLik))
attr(dat, "weight") <- weights

dat_Mk <- subset(dat, select = -c(1,342,683, 1024))


Q <- c(6:1)


pmlU1 <- pml(treeU1, dat)
pmlU2 <- pml(treeU2, dat)
pmlU2.fitted <- optim.pml(pmlU2, control = pml.control(trace=0))
pmlU3 <- pml(treeU3, dat)
pmlU3.fitted <- optim.pml(pmlU3, TRUE, control = pml.control(trace=0))

pmlR1 <- pml(treeR1, dat)
pmlR2 <- pml(treeR2, dat)
pmlR2.fitted <- optim.pml(pmlR2, optRooted = TRUE, control = 
                              pml.control(epsilon=1e-10, trace=0))
pmlR3 <- pml(treeR3, dat)
pmlR3.fitted <- optim.pml(pmlR3, TRUE, optRooted = TRUE,  control = 
                              pml.control(epsilon=1e-10, trace=0))


test_that("edge length optimisation works properly", {
    skip_on_cran()
    expect_equal(logLik(pmlU2.fitted), logLik(pmlU1))
    expect_equal(logLik(pmlR2.fitted), logLik(pmlR1))  
    expect_equal(pmlU2.fitted$tree, pmlU1$tree, tolerance=1e-6)
#    expect_equal(pmlR2.fitted$tree, pmlR1$tree, tolerance=5e-5)
})

test_that("NNI optimisation works properly", {
    skip_on_cran()
    expect_equal(logLik(pmlU3.fitted), logLik(pmlU1))
    expect_equal(logLik(pmlR3.fitted), logLik(pmlR1))
    expect_equal(pmlU3.fitted$tree, pmlU1$tree, tolerance=1e-6)
    expect_equal(storage.mode(pmlU3.fitted$tree$edge), "integer")
#    expect_equal(pmlR3.fitted$tree, pmlR1$tree, tolerance=5e-6)
})



test_that("bf optimisation works properly", {
    skip_on_cran()
    bf <- c(.1,.2,.3,.4)
    fit_T <- pml(treeU1, dat, bf=bf)
    weights <- as.vector(1000 * exp(fit_T$siteLik))
    dat_tmp <- dat
    attr(dat_tmp, "weight") <- weights
    fit0 <- pml(treeU1, dat_tmp)
    fit.bf <- optim.pml(fit0, optEdge=FALSE, optBf = TRUE, 
                        control = pml.control(epsilon=1e-10, trace=0))
    expect_equal(logLik(fit.bf), logLik(pml(treeU1, dat_tmp, bf=bf)))
    expect_equal(bf, fit.bf$bf, tolerance=5e-4)
})


test_that("Q optimisation works properly", {
    skip_on_cran()
    Q <- c(6:1)
    fit_T <- pml(treeU1, dat, Q=Q)
    weights <- as.vector(1000 * exp(fit_T$siteLik))
    dat_tmp <- dat
    attr(dat_tmp, "weight") <- weights
    fit0 <- pml(treeU1, dat_tmp)
    fit.Q <- optim.pml(fit0, optEdge=FALSE, optQ = TRUE, 
                       control = pml.control(epsilon=1e-10, trace=0))
    expect_equal(logLik(fit.Q), logLik(pml(treeU1, dat_tmp, Q=Q)))
    expect_equal(Q, fit.Q$Q, tolerance=5e-4)
})


test_that("Inv optimisation works properly", {
    skip_on_cran()
    inv <- 0.25
    fit_T <- pml(treeU1, dat, inv=inv)
    weights <- as.vector(1000 * exp(fit_T$siteLik))
    dat_tmp <- dat
    attr(dat_tmp, "weight") <- weights
    fit0 <- pml(treeU1, dat_tmp)
    fit.Inv <- optim.pml(fit0, optEdge=FALSE, optInv = TRUE, 
                         control = pml.control(epsilon=1e-10, trace=0))
    expect_equal(logLik(fit.Inv), logLik(pml(treeU1, dat_tmp, inv=inv)))
    expect_equal(inv, fit.Inv$inv, tolerance=5e-4)
})


test_that("Gamma optimisation works properly", {
    skip_on_cran()
    shape <- 2
    fit_T <- pml(treeU1, dat, shape=shape, k=4)
    weights <- as.vector(1000 * exp(fit_T$siteLik))
    dat_tmp <- dat
    attr(dat_tmp, "weight") <- weights
    fit0 <- pml(treeU1, dat_tmp, k=4)
    fit.Gamma <- optim.pml(fit0, optEdge=FALSE, optGamma = TRUE, 
                           control = pml.control(epsilon=1e-10, trace=0))
    expect_equal(logLik(fit.Gamma), logLik(pml(treeU1, dat_tmp, shape=shape, k=4)))
    expect_equal(shape, fit.Gamma$shape, tolerance=5e-4)
})


test_that("rate optimisation works properly", {
    skip_on_cran()
    rate <- 2
    fit_T <- pml(treeU1, dat, rate=rate)
    weights <- as.vector(1000 * exp(fit_T$siteLik))
    dat_tmp <- dat
    attr(dat_tmp, "weight") <- weights
    fit0 <- pml(treeU1, dat_tmp)
    fit.rate <- optim.pml(fit0, optEdge=FALSE, optRate = TRUE, 
                          control = pml.control(epsilon=1e-10, trace=0))
    expect_equal(logLik(fit.rate), logLik(pml(treeU1, dat_tmp, rate=rate)))
    expect_equal(rate, fit.rate$rate, tolerance=5e-4)
})

test_that("Mkv model works properly", {
    skip_on_cran()
    expect_equal(logLik(pmlU2.fitted), logLik(pmlU1))
})

