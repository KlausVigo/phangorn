data(woodmouse)

fit <- pml_bb(woodmouse, "HKY+I", rearrangement = "NNI",
              control=pml.control(trace=0))

expect_true(AICc(fit) > AIC(fit))
expect_true(BIC(fit) > AIC(fit))

c_mat <- vcov(fit)
expect_true( isSymmetric(c_mat, tol=1e-12) )
