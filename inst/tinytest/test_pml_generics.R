data(woodmouse)

fit <- pml_bb(woodmouse, "HKY+I", rearrangement = "NNI",
              control=pml.control(trace=0))

expect_true(AICc(fit) > AIC(fit))
expect_true(BIC(fit) > AIC(fit))

c_mat <- vcov(fit)
expect_true( isSymmetric(c_mat, tol=1e-12) )

# test write.pml
write.pml(fit, save_rds = TRUE)
fit2 <- readRDS("pml.rds")
expect_equal(fit, fit2)



#write.pml(fit, save_rds = FALSE)

#tree <- read.tree("pml_tree.nwk")
#align <- read.phyDat("pml_align.fasta", format="fasta")
#fit3 <- pml(tree = tree, data = align, k = 1L, site.rate = "gamma", ASC = FALSE,
#           bf = c(0.306541405706333, 0.261308281141267, 0.126026443980515,
#                  0.306123869171886),
#           Q = c(1, 23.0401288315219, 1, 1, 23.0401288315219, 1),
#           inv = 0.850681761772363, model = "HKY")

#fit$method <- NULL
#expect_equal(fit, fit3)
unlink(c("pml_align.fasta", "pml_tree.nwk", "pml.rds", "pml.txt"))
