tree <- read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
tree2 <- read.tree(text = "(((t1:1,t2:1):1,t3:2):1,t4:3);")


dat <- matrix(c("a", "a",
                "a", "t",
                "t", "a",
                "t", "t"), byrow = TRUE, nrow = 4L,
               dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
dna <- phyDat(dat)
fit <- pml(tree, dna)



# dna tests differs from other data types as it may returns ambiguous data
# test ancestral generics
#    test.ml1 <- ancestral.pml(fit, type = "ml")
test_ml <- ancestral.pml(fit, type = "ml", return = "phyDat")
test_mpr <- ancestral.pars(tree, dna, "MPR", return = "phyDat")
test_acctran <- ancestral.pars(tree, dna, "ACCTRAN", return = "phyDat")

expect_equal(as.character(test_ml), as.character(test_acctran))
expect_equal(as.character(test_ml), as.character(test_mpr))

test_mpr_2 <- ancestral.pars(tree2, dna, "MPR", return = "phyDat")
test_acctran_2 <- ancestral.pars(tree2, dna, "ACCTRAN", return = "phyDat")

expect_equal(test_mpr_2[,1], test_acctran_2[,1], check.attributes = FALSE)

