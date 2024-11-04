tree <- read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
tree2 <- read.tree(text = "(((t1:1,t2:1):1,t3:2):1,t4:3);")


dat <- matrix(c("a", "a",
                "a", "t",
                "t", "a",
                "t", "t"), byrow = TRUE, nrow = 4L,
               dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
dna <- phyDat(dat)
fit <- pml(tree, dna)

anc_p <- anc_pars(tree, dna)
expect_true( inherits(anc_p, "ancestral"))
dna_df <- as.data.frame(anc_p)
expect_equal(nrow(dna_df), 2 * (Ntip(tree) + Nnode(tree)))

# dna tests differs from other data types as it may returns ambiguous data
# test ancestral generics
#    test.ml1 <- ancestral.pml(fit, type = "ml")
#test_ml <- ancestral.pml(fit, type = "ml")
#test_mpr <- ancestral.pars(tree, dna, "MPR")
#test_acctran <- ancestral.pars(tree, dna, "ACCTRAN")

#expect_equal(as.character(test_ml), as.character(test_acctran))
#expect_equal(as.character(test_ml), as.character(test_mpr))

#test_mpr_2 <- ancestral.pars(tree2, dna, "MPR")
#test_acctran_2 <- ancestral.pars(tree2, dna, "ACCTRAN")

#expect_equal(test_mpr_2[,1], test_acctran_2[,1], check.attributes = FALSE)


data(Laurasiatherian)
fit <- pml_bb(Laurasiatherian[,1:100], "JC", rearrangement = "none",
              control=pml.control(trace=0))
anc_ml <- anc_pml(fit)
write.ancestral(anc_ml)
align <- read.phyDat("ancestral_align.fasta", format = "fasta")
tree <- read.tree("ancestral_tree.nwk")
df <- read.table("ancestral.state", header=TRUE)
anc_ml_disc <- as.ancestral(tree, df, align)
expect_equal(anc_ml$tree, anc_ml_disc$tree)
expect_equal(anc_ml$data[tree$tip.label], anc_ml_disc$data)
#expect_equal(anc_ml[[3]], anc_ml_disc[[3]])
expect_equal(anc_ml$state, anc_ml_disc$state)
unlink(c("ancestral_align.fasta", "ancestral_tree.nwk", "ancestral.state"))


## Felsenstein example
