context("likelihood")

treeU1 = read.tree(text = "((t1:.1,t2:.1):.1,t3:.2,(t4:.1,t5:.1):.3);")
treeU2 = read.tree(text = "((t1:.15,t2:.15):.05,t3:.2,(t4:.1,t5:.1):.3);") # edge length changed
treeU3 = read.tree(text = "((t1:.1,t3:.1):.1,t2:.2,(t4:.1,t5:.1):.3);") # 1 NNI move

treeR1 = read.tree(text = "(((t1:.1,t2:.1):.1,t3:.2):.1,(t4:.1,t5:.1):.2);")
treeR2 = read.tree(text = "(((t1:.15,t2:.15):.05,t3:.2):.1,(t4:.1,t5:.1):.2);") # edge length changed
treeR3 = read.tree(text = "(((t1:.1,t3:.1):.1,t2:.2):.1,(t4:.1,t5:.1):.2);") # 1 NNI move


# dat <- phyDat(c(t1="a", t2="a",t3="t",t4="t"), type="USER", levels=c("a","c","g","t"))
#tree2 = read.tree(text = "((t1,t3),t2,t4);")

dat = allSitePattern(5)
weights = as.vector(1000 * exp(pml(treeR1, dat)$siteLik))
attr(dat, "weight") = weights

bf = c(.1,.2,.3,.4)
Q = c(6:1)


pmlU1 = pml(treeU1, dat)
pmlU2 = pml(treeU2, dat)
pmlU2.fitted = optim.pml(pmlU2, control = pml.control(trace=0))
pmlU3 = pml(treeU3, dat)
pmlU3.fitted = optim.pml(pmlU3, TRUE, control = pml.control(trace=0))

pmlR1 = pml(treeR1, dat)
pmlR2 = pml(treeR2, dat)
pmlR2.fitted = optim.pml(pmlR2, optRooted = TRUE, control = pml.control(epsilon=1e-10, trace=0))
pmlR3 = pml(treeR3, dat)
pmlR3.fitted = optim.pml(pmlR3, TRUE, optRooted = TRUE,  control = pml.control(epsilon=1e-10, trace=0))



# optim.pml:
# bf F81
# Q  Sym
# Gamma
# Inv


#tr_acctran = acctran(tree1, dat)
#tr_ratchet = pratchet(dat, trace=0)
#bab(dat)
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
#    expect_equal(pmlR3.fitted$tree, pmlR1$tree, tolerance=5e-6)
})
