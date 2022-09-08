X <- allSitePattern(4)
attr(X, "type") <- "DNA"
tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,t3:0.2,t4:0.3);")
fit <- pml(tree, X, k=4, bf = c(1:4)/10)  # F81 + Gamma
weights <- 1000*exp(fit$siteLik)
attr(X, "weight") <- weights

set.seed(42)
# tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,t3:0.2);")
Y <- simSeq(tree, l=500, type = "AA", model="WAG")

# test modelTest
MT <- modelTest(X, tree = tree, I=FALSE,
                control = pml.control(epsilon = 1e-08, maxit = 10, trace = 0),
                multicore = TRUE, mc.cores = 2L)
expect_equal(MT$Model[which.min(MT$BIC)], "F81+G(4)")

fitMT <- as.pml(MT)
expect_true(inherits(fitMT, "pml"))
expect_equal(fitMT$model, "F81")

MT2 <- modelTest(X, tree = tree, model = c("JC", "F81", "K80", "HKY", "SYM",
          "GTR"), control = pml.control(epsilon = 1e-08, maxit = 10, trace = 0),
          multicore = TRUE, mc.cores = 2L)
expect_equal(MT2$Model[which.min(MT2$BIC)], "F81+G(4)")

# amino acid models with FREQ
MT_AA <- modelTest(Y, model=c("JTT", "WAG"), FREQ = TRUE,
            control = pml.control(epsilon = 1e-08, maxit = 5, trace = 0),
            multicore = TRUE, mc.cores = 2L)
expect_equal(MT_AA$Model[which.min(MT_AA$BIC)], "WAG")
# test all
MT_AA_all <- modelTest(Y, I=FALSE, G=FALSE,
             control = pml.control(epsilon = 1e-08, maxit = 5, trace = 0),
             multicore = TRUE, mc.cores = 2L)
expect_equal(MT_AA_all$Model[which.min(MT_AA_all$BIC)], "WAG")


# test user defined states
tree <- rcoal(10)
Z <- simSeq(tree, Q = c(1,0,0,1,0,1), type = "USER",
            levels=c("A", "B", "C", "D"))
MT_USER <- modelTest(Z, I=TRUE, G=TRUE,
                  control = pml.control(epsilon = 1e-08, maxit = 5, trace = 0),
                  multicore = TRUE, mc.cores = 2L)
expect_equal(MT_USER$Model[which.min(MT_USER$BIC)], "ORDERED")

