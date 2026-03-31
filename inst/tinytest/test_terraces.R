data(woodmouse)
fit <- pml_bb(woodmouse, model="JC", control=pml.control(trace=0))
trees <- unique(c(fit$tree, fit$bs))

dm <- dist.dna(woodmouse)


n_tree <- length(trees)
ll <- logLik(fit)
terrace_pml <- terraces(fit, plot=FALSE)

expect_equal(dim(terrace_pml), c(n_tree, 3L))
expect_true(all(terrace_pml[,3] <= ll + 1e-6))

terrace_dm <- terraces(dm, trees, plot=FALSE)
expect_equal(dim(terrace_dm), c(n_tree, 3L))

terrace_pars <- terraces(as.phyDat(woodmouse), trees, plot=FALSE)
expect_equal(dim(terrace_pars), c(n_tree, 3L))

