data(woodmouse)
set.seed(42)
fit <- pml_bb(woodmouse, "HKY+I", rearrangement = "NNI",
              control=pml.control(trace=0))

expect_true(AICc(fit) > AIC(fit))
expect_true(BIC(fit) > AIC(fit))

c_mat <- vcov(fit)
expect_true( isSymmetric(c_mat, tol=1e-8) )

# test write.pml
write.pml(fit, save_rds = TRUE)
fit2 <- readRDS("pml.rds")
expect_equal(fit, fit2)

# Add more checks (best expect_equal(fit, fit3))
# order of tip label changes (and of data object)
# attribute method might missing (when optim.pml was not called)
tree <- read.tree("pml_tree.nwk")
tree <- relabel(tree, fit$tree$tip.label)
align <- read.phyDat("pml_align.fasta", format="fasta")
fit3 <- pml(tree = tree, data = align, k = 1L, site.rate = "gamma", ASC = FALSE,
           bf = c(0.306541405706333, 0.261308281141267, 0.126026443980515,
                  0.306123869171886), Q = c(1, 23.0401029286372, 1, 1,
                  23.0401029286372, 1), inv = 0.850681184138605, model = "HKY")
expect_equal(logLik(fit), logLik(fit3))

unlink(c("pml_align.fasta", "pml_tree.nwk", "pml.rds", "pml.txt",
         "pml_rates.pdf", "pml_tree.pdf"))
