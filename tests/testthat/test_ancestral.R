context("ancestral")

# 
tree <- read.tree(text="((t1:1,t2:1):1,(t3:1,t4:1):1);")
dat <- matrix( c("a", "a",
                 "a", "t",
                 "t", "a",
                 "t", "t"), byrow = TRUE, nrow = 4L, 
               dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
dna <- phyDat(dat)
fit <- pml(tree, dna)



# dna tests differs from other data types as it may returns ambiguous data
test_that("test ancestral dna", {
    test.ml1 <- ancestral.pml(fit, type = "ml")
    test.ml2 <- ancestral.pml(fit, type = "ml", return = "phyDat")
    test1 <- ancestral.pars(tree, dna, "MPR", return = "prob")
    test2 <- ancestral.pars(tree, dna, "MPR", return = "phyDat")
    test3 <- ancestral.pars(tree, dna, "ACCTRAN", return = "prob")
    test4 <- ancestral.pars(tree, dna, "ACCTRAN", return = "phyDat")
    expect_equal(as.character(test2), as.character(test4))
    expect_equal(as.character(test.ml2), as.character(test2))
})

