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

dat_Mkv <- subset(dat, select = -c(1,342,683, 1024), site.pattern = TRUE)


Q <- c(6:1)


pmlU1 <- pml(treeU1, dat)
pmlU2 <- pml(treeU2, dat)
pmlU2.fitted <- optim.pml(pmlU2, control = pml.control(trace=0))
pmlU3 <- pml(treeU3, dat)
pmlU3.fitted <- optim.pml(pmlU3, TRUE, control =
                              pml.control(epsilon=1e-10, trace=0))

pmlR1 <- pml(treeR1, dat)
pmlR2 <- pml(treeR2, dat)
pmlR2.fitted <- optim.pml(pmlR2, optRooted = TRUE, control =
                              pml.control(epsilon=1e-10, trace=0))
pmlR3 <- pml(treeR3, dat)
pmlR3.fitted <- optim.pml(pmlR3, TRUE, optRooted = TRUE,  control =
                              pml.control(epsilon=1e-10, trace=0))


# test input parameters
# missing model
expect_error(pml_bb(dat))
# missing tip,dates
expect_error(pml_bb(dat, model="GTR", method="tipdated"))

## Parameter Optimisation
# test edge length optimisation
    expect_equal(logLik(pmlU2.fitted), logLik(pmlU1))
    expect_equal(logLik(pmlR2.fitted), logLik(pmlR1))
    expect_equal(pmlU2.fitted$tree, pmlU1$tree, tolerance=1e-6)
#    expect_equal(pmlR2.fitted$tree, pmlR1$tree, tolerance=5e-5)

# test NNI optimisation
    expect_equal(logLik(pmlU3.fitted), logLik(pmlU1))
    expect_equal(logLik(pmlR3.fitted), logLik(pmlR1), tolerance = 1e-06)
    expect_equal(pmlU3.fitted$tree, pmlU1$tree, tolerance=1e-6)
    expect_equal(storage.mode(pmlU3.fitted$tree$edge), "integer")
#    expect_equal(pmlR3.fitted$tree, pmlR1$tree, tolerance=5e-6)

# test bf optimisation
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
    fit.F81 <- pml_bb(dat_tmp, model="F81", control = pml.control(trace=0))
    expect_equal(logLik(fit.F81), logLik(pml(treeU1, dat_tmp, bf=bf)))
    expect_equal(bf, fit.F81$bf, tolerance=5e-4)
# test Q optimisation
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


# test Inv optimisation
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
    fit.JC.I <- pml_bb(dat_tmp, model="JC+I", control = pml.control(trace=0))
    expect_equal(logLik( fit.JC.I), logLik(pml(treeU1, dat_tmp, inv=inv)))
    expect_equal(inv, fit.JC.I$inv, tolerance=5e-4)

# test Gamma optimisation
    shape <- 2
    fit_T <- pml(treeU1, dat, shape=shape, k=4)
    weights <- as.vector(1000 * exp(fit_T$siteLik))
    dat_tmp <- dat
    attr(dat_tmp, "weight") <- weights
    fit0 <- pml(treeU1, dat_tmp, k=4)
    fit.Gamma <- optim.pml(fit0, optEdge=FALSE, optGamma = TRUE,
                           control = pml.control(epsilon=1e-10, trace=0))
    expect_equal(logLik(fit.Gamma),
                 logLik(pml(treeU1, dat_tmp, shape=shape, k=4)))
    expect_equal(shape, fit.Gamma$shape, tolerance=5e-4)
    fit.JC.G4 <- pml_bb(dat_tmp, model="JC+G(4)", control=pml.control(trace=0))
    expect_equal(logLik( fit.JC.G4),
                 logLik(pml(treeU1, dat_tmp, shape=shape, k=4)))
    expect_equal(shape, fit.JC.G4$shape, tolerance=5e-4)


# test free_rate
    fit0 <- pml(treeU1, dat_tmp, k=4, site.rate = "free_rate")
    fit.freerate <- optim.pml(fit0, optEdge=FALSE, optGamma = TRUE,
                              control = pml.control(epsilon=1e-10, trace=0))
#    expect_equal(discrete.gamma(2,4), fit.freerate$g, tolerance=1e-4)
    expect_equal(logLik(fit.freerate),
                 logLik(pml(treeU1, dat_tmp, shape=shape, k=4)))
    expect_equal(phangorn:::guess_model(fit.freerate), "JC+R(4)")

# test Laguerre quadrature
    fit0 <- pml(treeU1, dat_tmp, k=4, site.rate = "gamma_quadrature")
    fit.quadrature <- optim.pml(fit0, optEdge=FALSE, optGamma = TRUE,
                              control = pml.control(epsilon=1e-10, trace=0))
    expect_equal(fit.quadrature$shape, shape, tolerance=1e-4)

# test gamma rw
    fit0 <- pml(treeU1, dat_tmp, k=4, site.rate = "gamma_unbiased")
    fit.rw <- optim.pml(fit0, optEdge=FALSE, optGamma = TRUE,
                                control = pml.control(epsilon=1e-10, trace=0))
    expect_equal(fit.rw$shape, shape, tolerance=1e-4)



# test rate optimisation
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


# test Mkv model
    fit_Mk <- pml(treeR2, dat_Mkv)
    fit_Mk <- optim.pml(fit_Mk, optRooted = TRUE, control=pml.control(trace=0))
    fit_Mkv_1 <- pml(treeR2, dat_Mkv, ASC=TRUE)
    fit_Mkv_1 <- optim.pml(fit_Mkv_1, optRooted = TRUE,
                           control=pml.control(trace=0))
    fit_Mkv_2 <- pml(treeR3, dat_Mkv, ASC=TRUE)
    fit_Mkv_2 <- optim.pml(fit_Mkv_2, optRooted = TRUE, rearrangement = "NNI",
                           control=pml.control(trace=0))
    fit_Mkv_3 <- pml_bb(dat_Mkv, "JC+ASC", method="ultrametric",
                        rearrangement = "NNI", control=pml.control(trace=0))
    expect_equal(fit_Mkv_1$tree, treeR1, tolerance=1e-3)
    expect_equal(fit_Mkv_2$tree, treeR1, tolerance=1e-3)
    expect_equal(fit_Mkv_3$tree, treeR1, tolerance=1e-3)
    expect_true(sum(fit_Mk$tree$edge.length) > sum(fit_Mkv_1$tree$edge.length))
    expect_true(sum(fit_Mk$tree$edge.length) > sum(fit_Mkv_2$tree$edge.length))
    expect_true(sum(fit_Mk$tree$edge.length) > sum(fit_Mkv_3$tree$edge.length))
