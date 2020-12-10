# use lysin data & tree ?
tree <- read.tree(text = "(((t1:.1,t2:.1):.1,t3:.2):.1,(t4:.1,t5:.1):.2);")
set.seed(1)

dat_1 <- simSeq(tree, l=1000, type = "CODON", dnds=0.5, tstv=2)
dat_2 <- simSeq(tree, l=1000, type = "CODON", dnds=1, tstv=2)
dat_3 <- simSeq(tree, l=500, type = "CODON", dnds=3, tstv=2)

dat_4 <- c(dat_1, dat_2, dat_3)

fit_F1x4 <- pml(tree, dat_1, bf="F1x4")
fit_F3x4 <- pml(tree, dat_1, bf="F3x4")
fit_GY <- pml(tree, dat_1, bf="empirical")

# test dn/ds optimisation works properly
fit_GY_opt <- optim.pml(fit_GY, model="codon1",
                        control=pml.control(trace=0))
expect_true(fit_GY_opt$dnds < 1)
expect_true(fit_GY_opt$tstv > 1)
# fit_selection <- codonTest(tree, dat_4, control=pml.control(trace=0))


library(ape)
data(woodmouse)
dat_codon <- dna2codon(as.phyDat(woodmouse))
tree <- NJ(dist.ml(dat_codon))

# test M0, M1a optimisation works properly
fit_codon <- codonTest(tree, dat_codon, model = c("M0", "M1a"),
                       control = pml.control(maxit = 20))
expect_true(inherits(fit_codon, "codonTest"))

