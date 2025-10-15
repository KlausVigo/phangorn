data(woodmouse)
library(rgl)
fit <- pml_bb(woodmouse, model="JC")
trees <- unique(c(fit$tree, fit$bs))

n_tree <- length(trees)
ll <- logLik(fit)
terrace <- terraces(fit, plot=FALSE)

expect_equal(dim(terrace), c(n_tree, 3L))
expect_true(all(terrace[,3] <= ll + 1e-6))
